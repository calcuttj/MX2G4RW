#include "TFile.h"
#include "TTree.h"
#include "EDepSim/TG4Event.h"
#include "EDepSim/TG4Trajectory.h"
//#include "TG4Event.h"
//#include "TG4Trajectory.h"


#include "boost/program_options/options_description.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/parsers.hpp"

#include "cetlib/filepath_maker.h"
#include "cetlib/filesystem.h"
#include "cetlib/search_path.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"

#include "geant4reweight/src/ReweightBase/G4ReweightManager.hh"
#include "geant4reweight/src/ReweightBase/G4MultiReweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightStep.hh"

#include "mx2g4rw/Base/ReweightSuite.h"
#include "mx2g4rw/Base/Utilities.h"


using PDGMatParam_t = std::tuple<int, std::string, std::string>;
using WeightMap_t = std::map<PDGMatParam_t, std::vector<std::vector<double>>>;
using IDMap_t = std::map<PDGMatParam_t, std::vector<int>>;

class EDepTraj {
 public:
  EDepTraj(const TG4Trajectory & traj) {
    fParentID = traj.GetParentId();
    fPDG = traj.GetPDGCode();
    fID = traj.GetTrackId();
    fTraj = &traj;
  }
  int fParentID, fPDG, fID;
  const TG4Trajectory * fTraj;
};

std::vector<int> ValidPDGs {
  211, -211, 111,
  2212, 2112,
  321, -321, 311
};

std::map<int, std::string> InteractionNames {
  {211, "pi+Inelastic"},
  {-211, "pi-Inelastic"},
  {2212, "protonInelastic"},
  {2112, "neutronInelastic"},
  {321, "kaon+Inelastic"},
  {-321, "kaon-Inelastic"}
};
std::map<int, std::string> PDGStrings {
  {211, "piplus"},
  {-211, "piminus"},
  {2212, "proton"},
  {2112, "neutron"},
  {321, "kplus"},
  {-321, "kminus"}
};
//TODO
//Might have to check the processes in edep docs
bool IsInteraction(const TG4Trajectory & traj) {
  return (traj.Points.back().GetProcess() == 4 &&
          traj.Points.back().GetSubprocess() == 121);
};
class EDepSimEvent {
 public:
  EDepSimEvent(TG4Event * event) {
    std::map<int, std::vector<EDepTraj>> children;
    for (size_t i = 0; i < event->Trajectories.size(); ++i) {

      const auto & traj = event->Trajectories[i];
      EDepTraj edep_traj(traj);
      if (edep_traj.fParentID == -1) {
        fPrimaryPDG = edep_traj.fPDG;
        fPrimaryID = edep_traj.fID;
        fPrimaryTraj = &traj;
      }
      //std::cout << "Checking " << edep_traj.fPDG << std::endl;
      if (std::find(ValidPDGs.begin(), ValidPDGs.end(), edep_traj.fPDG)
          != ValidPDGs.end()) {
        children[edep_traj.fParentID].push_back(edep_traj);
      }
    }
    fPrimaryChildren = children[fPrimaryID];
  };

  bool IsInteraction() {
    return (fPrimaryTraj->Points.back().GetProcess() == 4 &&
            fPrimaryTraj->Points.back().GetSubprocess() == 121);
  };

  int fPrimaryPDG, fPrimaryID = -999;
  const TG4Trajectory * fPrimaryTraj = 0x0;
  std::vector<EDepTraj> fPrimaryChildren;
};

bool IsValidChild(int pdg) {
  return (std::find(ValidPDGs.begin(), ValidPDGs.end(), pdg)
          != ValidPDGs.end());
};

//Have to manually take care of the memory here
//I'm sorry this sucks so bad
G4ReweightTraj * MakeChild(const TG4Trajectory * edep_child) {
  auto * result = new G4ReweightTraj(
      edep_child->GetTrackId(),
      edep_child->GetPDGCode(),
      edep_child->GetParentId(),
      0, {0, 0});

  //Make a single step in case the reweighter depends on the specific
  //momentum 
  const auto & pt0 = edep_child->Points[0];
  const auto & pt1 = edep_child->Points[1];
  double pt0_p[3] = {
    pt0.GetMomentum()[0],
    pt0.GetMomentum()[1],
    pt0.GetMomentum()[2]
  };
  double pt1_p[3] = {
    pt1.GetMomentum()[0],
    pt1.GetMomentum()[1],
    pt1.GetMomentum()[2]
  };
  //Get first and second point and add pre/post step momentum
  G4ReweightStep * step = new G4ReweightStep(
      edep_child->GetTrackId(),
      edep_child->GetPDGCode(),
      edep_child->GetParentId(),
      0, pt0_p, pt1_p,
      0., "default");
  result->AddStep(step);
  return result;
}
void ClearChildren(G4ReweightTraj & traj) {
  for (auto * child : traj.GetChildren())
    delete child;
}

G4ReweightTraj MakeG4RWTraj(
    const TG4Trajectory & traj,
    const std::vector<const TG4TrajectoryPoint *> & segment,
    const std::string & last_proc) {
  //Sorry for the extra bullshit here
  G4ReweightTraj result(
      traj.GetTrackId(), traj.GetPDGCode(),
      traj.GetParentId(), 0, {0,0});

  for (size_t i = 1; i < segment.size(); ++i) {
    //If we're at the last step, set the proc
    //This will usually also be default, but the caller will 
    //set this to the actual interaction when needed
    std::string proc = (i == segment.size()-1 ? last_proc : "default");
    const auto * pt = segment[i];
    //std::string material = pt.GetMaterial();

    const auto * prev_pt = segment[i-1];
    auto dist = (pt->GetPosition() - prev_pt->GetPosition());
    double len = sqrt(dist.X()*dist.X() + dist.Y()*dist.Y() +
                      dist.Z()*dist.Z())/10.; //TODO -- Double check units
    double preStepP[3] = {
      prev_pt->GetMomentum()[0],
      prev_pt->GetMomentum()[1],
      prev_pt->GetMomentum()[2]
    };
    double postStepP[3] = {
      pt->GetMomentum()[0],
      pt->GetMomentum()[1],
      pt->GetMomentum()[2]
    };

    double total_pre_p = sqrt(preStepP[0]*preStepP[0] +
                              preStepP[1]*preStepP[1] +
                              preStepP[2]*preStepP[2]);

    double total_post_p = sqrt(postStepP[0]*postStepP[0] +
                              postStepP[1]*postStepP[1] +
                              postStepP[2]*postStepP[2]);

    std::cout << i << " " << pt->GetPosition().X()/10. << " " <<
                 pt->GetPosition().Y()/10. << " " <<
                 pt->GetPosition().Z()/10. << " " << pt->GetProcess() <<
                 " " << pt->GetSubprocess() << " " <<
                 len <<
                 " " << total_pre_p << " " << total_post_p << " " <<
                 pt->GetMaterial() << std::endl;

    //Sorry for my poor memory management back in the day
    G4ReweightStep * step = new G4ReweightStep(
        traj.GetTrackId(), traj.GetPDGCode(), traj.GetParentId(),
        0, preStepP, postStepP,
        len, proc);
    result.AddStep(step);
  }
  return result;
};

std::string MakeBranchName(PDGMatParam_t pdg_mat_param) {
  return (PDGStrings[std::get<0>(pdg_mat_param)] + "_" +
          std::get<1>(pdg_mat_param) + "_" +
          std::get<2>(pdg_mat_param));
};

int main(int argc, char ** argv) {
  namespace po = boost::program_options;

  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("input,i", po::value<std::string>(), "Input file")
      ("config,c", po::value<std::string>(), "Fcl file")
      ("output,o", po::value<std::string>(), "Output file");
  
  po::variables_map vm;
  po::store(
      po::command_line_parser(argc, argv).options(desc).style(
          po::command_line_style::unix_style |
          po::command_line_style::allow_dash_for_short).run(),
      vm);
  po::notify(vm);    
  
  if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
  }

  std::string input_file = "";
  if (vm.count("input")) {
    input_file = vm["input"].as<std::string>();
  }
  else {
    std::cout << "Need file" << std::endl;
    return 1;
  }

  std::string output_file = "";
  if (vm.count("output")) {
    output_file = vm["output"].as<std::string>();
  }
  else {
    std::cout << "Need output file" << std::endl;
    return 1;
  }

  std::string fcl_file = "";
  if (vm.count("config")) {
    fcl_file = vm["config"].as<std::string>();
  }
  else {
    std::cout << "Need fcl file" << std::endl;
    return 1;
  }

  //
  //Get Parameters
  auto fcl_pars = mx2g4rw::MakeFCLPars(fcl_file);

  mx2g4rw::ReweightSuite reweight_suite(fcl_pars);

  //Inputfile 
  TFile f(input_file.c_str(), "open");
  //TODO -- ADD CHECK HERE FOR GOOD FILE TYPE/HAS TREE ETC.
  TTree * tree = (TTree*)f.Get("EDepSimEvents");

  TG4Event * event = 0x0;
  tree->SetBranchAddress("Event", &event);
  std::cout << "Opened file" << std::endl;

  //Initialize Weight Map
  
  ////Output File
  TFile fOut(output_file.c_str(), "recreate");
  TTree outtree("tree", "");

  WeightMap_t weight_map;
  IDMap_t id_map;
  for (const auto & [part_mat, par_names] :
       reweight_suite.GetParameterNames()) {
    int pdg = part_mat.first;
    auto material = part_mat.second;

    for (const auto & name : par_names) {
      auto part_mat_param = std::make_tuple(pdg, material, name);
      weight_map[part_mat_param] = {};
      id_map[part_mat_param] = {};
      //coeff_map[part_mat_param] = {};
      std::string weight_br_name = MakeBranchName(part_mat_param) + "_weights";
      std::string id_br_name = MakeBranchName(part_mat_param) + "_ID";
      outtree.Branch(weight_br_name.c_str(),
                     &weight_map[part_mat_param]);
      outtree.Branch(id_br_name.c_str(),
                     &id_map[part_mat_param]);
    }
  }

  //double weight = 1., old_weight = 1.;
  //double out_len = 0.;
  //int npiplus = 0, npiminus = 0, npi0 = 0, nproton = 0, nneutron = 0, nkaon = 0;
  //int nchildren = 0;
  //bool interacted = false, interacted_in_lar = false, interacted_in_fr4 = false;
  //double max_p = -999, max_p_costheta = -999, final_p = -999;
  //double end_x = -999.;
  //std::vector<double> piplus_p, product_p;
  //std::vector<double> lar_incidents, fr4_incidents;
  //std::map<int, int> n_children;
  //outtree.Branch("weight", &weight);
  //outtree.Branch("old_weight", &old_weight);
  //outtree.Branch("len", &out_len);
  //outtree.Branch("max_p", &max_p);
  //outtree.Branch("end_x", &end_x);
  //outtree.Branch("max_p_costheta", &max_p_costheta);
  //outtree.Branch("interacted", &interacted);
  //outtree.Branch("interacted_in_lar", &interacted_in_lar);
  //outtree.Branch("interacted_in_fr4", &interacted_in_fr4);
  //outtree.Branch("nchildren", &nchildren);

  //outtree.Branch("npiplus", &npiplus);
  //outtree.Branch("piplus_p", &piplus_p);
  //outtree.Branch("product_p", &product_p);
  //outtree.Branch("npiminus", &npiminus);
  //outtree.Branch("npi0", &npi0);
  //outtree.Branch("nproton", &nproton);
  //outtree.Branch("nneutron", &nneutron);
  //outtree.Branch("nkaon", &nkaon);
  //outtree.Branch("final_p", &final_p);

  //outtree.Branch("lar_incidents", &lar_incidents);
  //outtree.Branch("fr4_incidents", &fr4_incidents);



  //Loop over edep sim events
  int nentries = tree->GetEntries();
  for (int i = 0; i < nentries; ++i) {
    if (!(i%1000)) std::cout << i << "/" << nentries << "\r";
    std::cout << "############" << std::endl;
    tree->GetEntry(i);

    //Build-up map to children, we'll need to look up the 
    //children of each particle in question later on
    std::map<int, std::vector<const TG4Trajectory *>> ChildrenMap;
    for (size_t i = 0; i < event->Trajectories.size(); ++i) {
      const auto & traj = event->Trajectories[i];

      //Skip the ones we don't care about and primaries (id == -1)
      if (!IsValidChild(traj.GetPDGCode()) || traj.GetParentId() == -1)
        continue;
      ChildrenMap[traj.GetParentId()].push_back(&traj);
    }

    //Start iterating again to actually reweight
    for (size_t i = 0; i < event->Trajectories.size(); ++i) {
      const auto & traj = event->Trajectories[i];
      int fParentID = traj.GetParentId();
      int fPDG = traj.GetPDGCode();
      int fID = traj.GetTrackId();
      
      //Check if we want to reweight this particle at all
      bool need_to_weight = reweight_suite.CheckPDG(fPDG);
      std::cout << fID << " " << fPDG << " " << fParentID << " " <<
                   need_to_weight << std::endl;
      if (!need_to_weight) continue;

      ////TODO fix this comment
      //Loop over the children
      std::cout << "Reweightable " << fPDG << " Has children: " << std::endl;
      if (ChildrenMap.find(fID) != ChildrenMap.end()) {
        for (const auto * child : ChildrenMap.at(traj.GetTrackId())) {
          std::cout << "\t" << child->GetPDGCode() << std::endl;
        }
        if (!IsInteraction(traj) && ChildrenMap.at(traj.GetTrackId()).size())
          std::cout << "Warning. No interaction but " <<
                        ChildrenMap.at(traj.GetTrackId()).size() << std::endl;
      }

      //Build up the reweightable trajectories.
      //We'll look for contiguous steps in a given (reweightable) material and
      //consider these as a single 'track segment' in that material'.
      //If a particle passes through mutliple materials and within that,
      //a given reweightable material multiple times, we need to multiply those
      //'sub weights' together to give a single weight for that material

      size_t npts = traj.Points.size();
      std::cout << "Traj has npts: " << npts << std::endl;

      //Maybe TODO? replace vector with a deque and pop front to save memory?
      //If that's necessary I guess...
      std::map<std::string,
               std::vector<std::vector<const TG4TrajectoryPoint *>>>
          reweightable_segments;

      bool new_segment = true; //Use this to bookkeep whether a point is in a
                               //contiguous segment

      //Loop over trajectory points for this trajectory
      for (const auto & pt : traj.Points) {
        const std::string & material = pt.GetMaterial();
        std::cout << "\t"<< material << std::endl;


        //Check if we're interested in this material 
        bool is_reweightable_material = reweight_suite.CheckMaterial(material);
        std::cout << "\t\tReweightable? " <<
                     is_reweightable_material << std::endl;
        if (!is_reweightable_material) {
          //Set to true, next time we're in a reweightable space,
          //we'll have to add a segment
          new_segment = true;
          continue;
        }

        //Make a new sub-segment in the material
        if (new_segment)
          reweightable_segments[material].push_back(
            std::vector<const TG4TrajectoryPoint *>()
          );

        //Set this to false now,
        //we'll stay in this segment until we hit a
        //non-reweightable material
        new_segment = false;

        //Add to the vector at the back (the newest segment)
        reweightable_segments[material].back().push_back(&pt);
        std::cout << "Segment now has " <<
                     reweightable_segments[material].back().size() << " pts" <<
                     std::endl;
      }

      //No interesting segments?
      //Skip this trajectory
      //TODO -- NEED TO SET WEIGHT TO 1. ONCE WE ACTUALLY START
      //IMPLEMENTING THAT PART
      if (reweightable_segments.size() == 0)
        continue;

      //Determine if 1) the particle interacted and
      //2) the interaction is in an interesting material
      const std::string last_material = traj.Points.back().GetMaterial();
      std::cout << "Interaction? " << IsInteraction(traj) << std::endl;
      std::cout << "Last Material: " << last_material << std::endl;
      std::string interaction_material = "";
      if (IsInteraction(traj) && reweight_suite.CheckMaterial(last_material)) {
        std::cout << "\tIs Reweightable" << std::endl;
        interaction_material = last_material; 
      }


      //Iterate through the materials and
      //make reweightable trajectories for each (sub)segment
      //Along the weight, check that the material in question
      //Is one in which the particle interacted
      for (const auto & [material, segments] : reweightable_segments) {

        //Loop over the sub segments and get the full weights
        //We have a 2D vector of weights
        //(N parameters) x (N Scan Steps)
        //And we'll multiply segment by segment
        size_t n_pars = reweight_suite.GetNPars({traj.GetPDGCode(), material});
        std::vector<std::vector<double>> weight_scans;
        for (size_t j = 0; j < segments.size(); ++j) {
          const auto & segment = segments[j];
          bool set_interaction = (
              (material == interaction_material) && (j == segments.size()-1)
          );
          std::string proc = (set_interaction ?
                              InteractionNames[traj.GetPDGCode()] :
                              "default"
                             );
          auto reweightable_traj = MakeG4RWTraj(traj, segment, proc);
          std::cout << "Made rw traj: " << reweightable_traj.GetTrackID() <<
                       " " << reweightable_traj.GetPDG() << " " <<
                       reweightable_traj.GetParID() << std::endl;
          //If this is an interaction, add the children 
          if (set_interaction) {
            const auto & children = ChildrenMap.at(traj.GetTrackId());
            for (const auto & child : children) {
              reweightable_traj.AddChild(MakeChild(child));
            }
          }


          for (size_t k = 0; k < n_pars; ++k) {
            //1 weight per scan step
            auto weights_k = reweight_suite.Scan(
                reweightable_traj,
                traj.GetPDGCode(), material, k
                //Scan steps, start, end, configurable here
            );
            if (j == 0) {
              weight_scans.push_back(weights_k);
            }
            else {
              for (size_t l = 0; l < weights_k.size(); ++l) {
                //Multiply weight from lth scan step for parmaeter k
                //by new subweight
                weight_scans[k][l] *= weights_k[l];
              }
            }
          }

          //Clear out any children when we're done
          ClearChildren(reweightable_traj);
        }
        std::cout << "Finished getting weights for trajectory " <<
                     traj.GetTrackId() << std::endl;
        for (size_t j = 0; j < weight_scans.size(); ++j) {
          std::string par_name
              = reweight_suite.GetParameterNames().at({fPDG, material})[j];
          std::cout << "\tPar " << j << " " << par_name << " has weights" <<
                       std::endl;


          auto pdg_mat_param = std::make_tuple(
              fPDG, material, par_name);
          weight_map[pdg_mat_param].push_back(weight_scans[j]);
          id_map[pdg_mat_param].push_back(fID);
          /*for (size_t k = 0; k < weight_scans[j].size(); ++k) {
            std::cout << "\t\t" << k << " " << weight_scans[j][k] << std::endl;
          }*/
        }

      }
    }
  //  EDepSimEvent edep_event(event);

  //  lar_incidents.clear();
  //  fr4_incidents.clear();
  //  interacted_in_lar = false;
  //  interacted_in_fr4 = false;
  //  piplus_p.clear();
  //  product_p.clear();
  //  npiplus = 0;
  //  nkaon = 0;
  //  npiminus = 0;
  //  npi0 = 0;
  //  nproton = 0;
  //  nneutron = 0;
  //  nchildren = 0;

  //  max_p = -999;
  //  end_x = -999.;
  //  max_p_costheta = -999;


  //  auto last_pt = edep_event.fPrimaryTraj->Points[edep_event.fPrimaryTraj->Points.size()-2];
  //  double traj_p = sqrt(last_pt.GetMomentum()[0]*last_pt.GetMomentum()[0] +
  //                         last_pt.GetMomentum()[1]*last_pt.GetMomentum()[1] +
  //                         last_pt.GetMomentum()[2]*last_pt.GetMomentum()[2]);

  //  final_p = traj_p;
  //  interacted = edep_event.IsInteraction();
  //  nchildren = edep_event.fPrimaryChildren.size();
  //  for (const auto & ct : edep_event.fPrimaryChildren) {
  //    //std::cout << "\tChild PDG: " << ct.fPDG << std::endl;
  //    int pdg = ct.fPDG;

  //    auto pt0 = ct.fTraj->Points[0];
  //    std::cout << pdg << " Prim child npts: " << ct.fTraj->Points.size() << std::endl;
  //    size_t npts = ct.fTraj->Points.size();
  //    for (size_t j = 0; j < npts; ++j) {
  //      const auto & pt = ct.fTraj->Points[j];
  //      double p = sqrt(pt.GetMomentum()[0]*pt.GetMomentum()[0] +
  //                      pt.GetMomentum()[1]*pt.GetMomentum()[1] +
  //                      pt.GetMomentum()[2]*pt.GetMomentum()[2]);
  //      std::cout << "\t" << p << " " << pt.GetMaterial() << std::endl;
  //    }

  //    double p = sqrt(pt0.GetMomentum()[0]*pt0.GetMomentum()[0] +
  //                           pt0.GetMomentum()[1]*pt0.GetMomentum()[1] +
  //                           pt0.GetMomentum()[2]*pt0.GetMomentum()[2]);
  //    product_p.push_back(p);
  //    if (p > max_p) {
  //      max_p = p;
  //      max_p_costheta = (
  //        pt0.GetMomentum()[0]*last_pt.GetMomentum()[0] +
  //        pt0.GetMomentum()[1]*last_pt.GetMomentum()[1] +
  //        pt0.GetMomentum()[2]*last_pt.GetMomentum()[2]
  //      )/(p*traj_p);
  //    }


  //    if (pdg == 211) {
  //      ++npiplus;
  //      piplus_p.push_back(p);
  //    }
  //    else if (pdg == -211)
  //      ++npiminus;
  //    else if (pdg == 111)
  //      ++npi0;
  //    else if (pdg == 2212)
  //      ++nproton;
  //    else if (pdg == 2112)
  //      ++nneutron;
  //    else if (abs(pdg) == 321)
  //      ++nkaon;
  //    else
  //      continue;
  //  }

  //  if (edep_event.IsInteraction()) ++ninteractions;

  //  //Start building reweightable trajectory
  //  G4ReweightTraj primary_traj(
  //      edep_event.fPrimaryID,
  //      edep_event.fPrimaryPDG,
  //      -1, i, {0,0});
  //  for (auto & child_traj : edep_event.fPrimaryChildren) {
  //    primary_traj.AddChild(
  //      new G4ReweightTraj(child_traj.fID, child_traj.fPDG,
  //                         edep_event.fPrimaryID, i, {0, 0})
  //    );
  //  }

  //  //Build the path it took
  //  const auto & pt0 = edep_event.fPrimaryTraj->Points[0];
  //  //std::cout << 0 << " " << pt0.GetPosition().X()/10. << " " <<
  //  //             pt0.GetPosition().Y()/10. << " " << pt0.GetPosition().Z()/10. <<
  //  //             " " << pt0.GetProcess() << " " << pt0.GetSubprocess() <<
  //  //             " " << sqrt(pt0.GetMomentum()[0]*pt0.GetMomentum()[0] +
  //  //                         pt0.GetMomentum()[1]*pt0.GetMomentum()[1] +
  //  //                         pt0.GetMomentum()[2]*pt0.GetMomentum()[2]) <<
  //  //             std::endl;

  //  //bool prev_is_LAr = pt0.GetMaterial() == "LAr";

  //  out_len = 0.;
  //  for (size_t j = 1; j < edep_event.fPrimaryTraj->Points.size(); ++j) {
  //    std::string proc = "default";
  //    if ((j == edep_event.fPrimaryTraj->Points.size()-1) &&
  //        edep_event.IsInteraction()) {//Last step
  //      proc = "pi+Inelastic";
  //    }
  //    const auto & pt = edep_event.fPrimaryTraj->Points[j];

  //    std::string material = pt.GetMaterial();
  //    bool is_LAr = pt.GetMaterial() == "LAr";
  //    bool is_FR4 = pt.GetMaterial() == "FR4";
  //    if (proc == "pi+Inelastic" && is_LAr) interacted_in_lar = true;
  //    if (proc == "pi+Inelastic" && is_FR4) interacted_in_fr4 = true;

  //    const auto & prev_pt = edep_event.fPrimaryTraj->Points[j-1];
  //    auto dist = (pt.GetPosition() - prev_pt.GetPosition());
  //    double len = sqrt(dist.X()*dist.X() + dist.Y()*dist.Y() +
  //                      dist.Z()*dist.Z())/10.;
  //    out_len += len;
  //    double preStepP[3] = {
  //      prev_pt.GetMomentum()[0],
  //      prev_pt.GetMomentum()[1],
  //      prev_pt.GetMomentum()[2]
  //    };
  //    double postStepP[3] = {
  //      pt.GetMomentum()[0],
  //      pt.GetMomentum()[1],
  //      pt.GetMomentum()[2]
  //    };

  //    double total_pre_p = sqrt(preStepP[0]*preStepP[0] +
  //                              preStepP[1]*preStepP[1] +
  //                              preStepP[2]*preStepP[2]);

  //    double total_post_p = sqrt(postStepP[0]*postStepP[0] +
  //                              postStepP[1]*postStepP[1] +
  //                              postStepP[2]*postStepP[2]);

  //    std::cout << j << " " << pt.GetPosition().X()/10. << " " <<
  //                 pt.GetPosition().Y()/10. << " " <<
  //                 pt.GetPosition().Z()/10. << " " << pt.GetProcess() <<
  //                 " " << pt.GetSubprocess() << " " <<
  //                 len <<
  //                 " " << total_pre_p << " " << total_post_p << " " <<
  //                 pt.GetMaterial() << std::endl;

  //    G4ReweightStep * step = new G4ReweightStep(edep_event.fPrimaryID,
  //                                               edep_event.fPrimaryPDG,
  //                                               0, i, preStepP, postStepP,
  //                                               len, proc);
  //    primary_traj.AddStep(step);
  //  }

  //  auto & end_pt
  //      = edep_event.fPrimaryTraj->Points[
  //          edep_event.fPrimaryTraj->Points.size()-1];
  //  end_x = end_pt.GetPosition().X()/10.;

  //  
  // // std::vector<std::vector<size_t>> segments {};
  // // bool prev_is_LAr = (pt0.GetMaterial() == "LAr");
  // // std::cout << "0 " << pt0.GetMaterial() << std::endl;
  // // for (size_t j = 1; j < edep_event.fPrimaryTraj->Points.size(); ++j) {
  // //   const auto & pt = edep_event.fPrimaryTraj->Points[j];
  // //   std::cout << j << " " << pt.GetMaterial() << std::endl;
  // //   bool is_LAr = pt.GetMaterial() == "LAr";
  // //   if (is_LAr) {
  // //     if (!prev_is_LAr || (segments.size() == 0)) {
  // //       segments.push_back({j-1});
  // //     }
  // //     else {
  // //       segments.back().push_back(j-1);
  // //     }
  // //   }
  // //   else {
  // //     if (prev_is_LAr) {
  // //       segments.back().push_back(j-1);
  // //     }
  // //   }

  // //   if (j == edep_event.fPrimaryTraj->Points.size()-1 && prev_is_LAr) {
  // //     segments.back().push_back(j);
  // //   }
  // //   prev_is_LAr = is_LAr;
  // // }

  // // for (auto & seg : segments) {
  // //   if (seg.size() < 2) {
  // //     throw cet::exception("do_weighting") << "Found segment < 2 size " << seg.size();
  // //   }
  // // }

  //  //Split up the traj by contiguous weightable material (LAr)
  //  G4ReweightTraj traj(
  //      edep_event.fPrimaryID,
  //      edep_event.fPrimaryPDG,
  //      -1, i, {0,0});
  //  for (auto & child_traj : edep_event.fPrimaryChildren) {
  //    traj.AddChild(
  //      new G4ReweightTraj(child_traj.fID, child_traj.fPDG,
  //                         edep_event.fPrimaryID, i, {0, 0})
  //    );
  //  }

  //  size_t npts = edep_event.fPrimaryTraj->Points.size();
  //  for (size_t j = 1; j < npts; ++j) {
  //    const auto & pt = edep_event.fPrimaryTraj->Points[j];
  //    if (pt.GetMaterial() != "LAr") continue;
  //    const auto & prev_pt = edep_event.fPrimaryTraj->Points[j-1];
  //    auto dist = (pt.GetPosition() - prev_pt.GetPosition());
  //    double len = sqrt(dist.X()*dist.X() + dist.Y()*dist.Y() +
  //                      dist.Z()*dist.Z())/10.;
  //    double preStepP[3] = {
  //      prev_pt.GetMomentum()[0],
  //      prev_pt.GetMomentum()[1],
  //      prev_pt.GetMomentum()[2]
  //    };
  //    double postStepP[3] = {
  //      pt.GetMomentum()[0],
  //      pt.GetMomentum()[1],
  //      pt.GetMomentum()[2]
  //    };

  //    double total_pre_p = sqrt(preStepP[0]*preStepP[0] +
  //                              preStepP[1]*preStepP[1] +
  //                              preStepP[2]*preStepP[2]);

  //    double total_post_p = sqrt(postStepP[0]*postStepP[0] +
  //                              postStepP[1]*postStepP[1] +
  //                              postStepP[2]*postStepP[2]);

  //    std::string proc = "default";
  //    if ((j == npts-1) && edep_event.IsInteraction()) {//Last step
  //      proc = "pi+Inelastic";
  //    }
  //    G4ReweightStep * step = new G4ReweightStep(edep_event.fPrimaryID,
  //                                               edep_event.fPrimaryPDG,
  //                                               0, i, preStepP, postStepP,
  //                                               len, proc);
  //    traj.AddStep(step);
  //    std::cout << "Added step " << j-1 << " " << j << " " << proc << std::endl;
  //  }
  //  std::cout << std::endl;

  //  //Get the weight for this trajectory/event
  //  weight = reweighter->GetWeightFromSetParameters(traj);
  //  //for (auto & traj : trajs) {
  //  //  weight *= reweighter->GetWeightFromSetParameters(*traj);
  //  //}
  //  old_weight = reweighter->GetWeightFromSetParameters(primary_traj);
  //  std::cout << "Weight: " << weight << " " << old_weight << std::endl;;
  //  if (edep_event.IsInteraction()) nweighted += weight;
  //  ntotal += weight;


    outtree.Fill();
  }
  outtree.Write();
  fOut.Close();

  f.Close();

  return 0;
}
