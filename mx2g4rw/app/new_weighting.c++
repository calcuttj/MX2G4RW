#include "TFile.h"
#include "TTree.h"
//#include "EDepSim/TG4Event.h"
//#include "EDepSim/TG4Trajectory.h"
#include "TG4Event.h"
#include "TG4Trajectory.h"


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
  321, -321
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

int main(int argc, char ** argv) {
  namespace po = boost::program_options;

  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produce help message")
      ("i", po::value<std::string>(), "Input file")
      ("c", po::value<std::string>(), "Fcl file")
      ("o", po::value<std::string>(), "Output file");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    
  
  if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
  }

  std::string input_file = "";
  if (vm.count("i")) {
    input_file = vm["i"].as<std::string>();
  }
  else {
    std::cout << "Need file" << std::endl;
    return 1;
  }

  std::string output_file = "";
  if (vm.count("o")) {
    output_file = vm["o"].as<std::string>();
  }
  else {
    std::cout << "Need output file" << std::endl;
    return 1;
  }

  std::string fcl_file = "";
  if (vm.count("c")) {
    fcl_file = vm["c"].as<std::string>();
  }
  else {
    std::cout << "Need fcl file" << std::endl;
    return 1;
  }

  //
  //Get Parameters
  auto fcl_pars = mx2g4rw::MakeFCLPars(fcl_file);

  auto reweight_suite = mx2g4rw::ReweightSuite(fcl_pars);
  return 0;

  //auto material = fcl_pars.get<fhicl::ParameterSet>("Material");
  //auto fracs_file_name = fcl_pars.get<std::string>("FracsFile");
  //auto g4rw_pars = fcl_pars.get<std::vector<fhicl::ParameterSet>>("ParameterSet");
  //auto set_weight = fcl_pars.get<double>("Weight", 2.);

  ////Open the file holding the final state fractions
  //TFile * fracs_file = OpenFile(fracs_file_name);
  //std::cout << "Fracs file: " << fracs_file << std::endl;

  ////Set up reweighting objects
  //G4ReweightManager * RWManager = new G4ReweightManager({material});
  //G4MultiReweighter * reweighter = new G4MultiReweighter(
  //  211, *fracs_file, g4rw_pars, material, RWManager);
  //reweighter->SetAllParameterValues({set_weight}); //Set the weight

  //Inputfile 
  TFile f(input_file.c_str(), "open");
  TTree * tree = (TTree*)f.Get("EDepSimEvents");
  TG4Event * event = 0x0;
  tree->SetBranchAddress("Event", &event);

  int ninteractions = 0;
  double nweighted = 0., ntotal = 0.;

  //Output File
  TFile fOut(output_file.c_str(), "recreate");
  TTree outtree("tree", "");
  double weight = 1., old_weight = 1.;
  double out_len = 0.;
  int npiplus = 0, npiminus = 0, npi0 = 0, nproton = 0, nneutron = 0, nkaon = 0;
  int nchildren = 0;
  bool interacted = false, interacted_in_lar = false, interacted_in_fr4 = false;
  double max_p = -999, max_p_costheta = -999, final_p = -999;
  double end_x = -999.;
  std::vector<double> piplus_p, product_p;
  std::vector<double> lar_incidents, fr4_incidents;
  std::map<int, int> n_children;
  outtree.Branch("weight", &weight);
  outtree.Branch("old_weight", &old_weight);
  outtree.Branch("len", &out_len);
  outtree.Branch("max_p", &max_p);
  outtree.Branch("end_x", &end_x);
  outtree.Branch("max_p_costheta", &max_p_costheta);
  outtree.Branch("interacted", &interacted);
  outtree.Branch("interacted_in_lar", &interacted_in_lar);
  outtree.Branch("interacted_in_fr4", &interacted_in_fr4);
  outtree.Branch("nchildren", &nchildren);

  outtree.Branch("npiplus", &npiplus);
  outtree.Branch("piplus_p", &piplus_p);
  outtree.Branch("product_p", &product_p);
  outtree.Branch("npiminus", &npiminus);
  outtree.Branch("npi0", &npi0);
  outtree.Branch("nproton", &nproton);
  outtree.Branch("nneutron", &nneutron);
  outtree.Branch("nkaon", &nkaon);
  outtree.Branch("final_p", &final_p);

  outtree.Branch("lar_incidents", &lar_incidents);
  outtree.Branch("fr4_incidents", &fr4_incidents);

  //Loop over edep sim events
  int nentries = tree->GetEntries();
  for (int i = 0; i < nentries; ++i) {
    if (!(i%1000)) std::cout << i << "/" << nentries << "\r";
    //std::cout << "############" << std::endl;
    tree->GetEntry(i);
    EDepSimEvent edep_event(event);

    lar_incidents.clear();
    fr4_incidents.clear();
    interacted_in_lar = false;
    interacted_in_fr4 = false;
    piplus_p.clear();
    product_p.clear();
    npiplus = 0;
    nkaon = 0;
    npiminus = 0;
    npi0 = 0;
    nproton = 0;
    nneutron = 0;
    nchildren = 0;

    max_p = -999;
    end_x = -999.;
    max_p_costheta = -999;


    auto last_pt = edep_event.fPrimaryTraj->Points[edep_event.fPrimaryTraj->Points.size()-2];
    double traj_p = sqrt(last_pt.GetMomentum()[0]*last_pt.GetMomentum()[0] +
                           last_pt.GetMomentum()[1]*last_pt.GetMomentum()[1] +
                           last_pt.GetMomentum()[2]*last_pt.GetMomentum()[2]);

    final_p = traj_p;
    interacted = edep_event.IsInteraction();
    nchildren = edep_event.fPrimaryChildren.size();
    for (const auto & ct : edep_event.fPrimaryChildren) {
      //std::cout << "\tChild PDG: " << ct.fPDG << std::endl;
      int pdg = ct.fPDG;

      auto pt0 = ct.fTraj->Points[0];
      std::cout << pdg << " Prim child npts: " << ct.fTraj->Points.size() << std::endl;
      size_t npts = ct.fTraj->Points.size();
      for (size_t j = 0; j < npts; ++j) {
        const auto & pt = ct.fTraj->Points[j];
        double p = sqrt(pt.GetMomentum()[0]*pt.GetMomentum()[0] +
                        pt.GetMomentum()[1]*pt.GetMomentum()[1] +
                        pt.GetMomentum()[2]*pt.GetMomentum()[2]);
        std::cout << "\t" << p << " " << pt.GetMaterial() << std::endl;
      }

      double p = sqrt(pt0.GetMomentum()[0]*pt0.GetMomentum()[0] +
                             pt0.GetMomentum()[1]*pt0.GetMomentum()[1] +
                             pt0.GetMomentum()[2]*pt0.GetMomentum()[2]);
      product_p.push_back(p);
      if (p > max_p) {
        max_p = p;
        max_p_costheta = (
          pt0.GetMomentum()[0]*last_pt.GetMomentum()[0] +
          pt0.GetMomentum()[1]*last_pt.GetMomentum()[1] +
          pt0.GetMomentum()[2]*last_pt.GetMomentum()[2]
        )/(p*traj_p);
      }


      if (pdg == 211) {
        ++npiplus;
        piplus_p.push_back(p);
      }
      else if (pdg == -211)
        ++npiminus;
      else if (pdg == 111)
        ++npi0;
      else if (pdg == 2212)
        ++nproton;
      else if (pdg == 2112)
        ++nneutron;
      else if (abs(pdg) == 321)
        ++nkaon;
      else
        continue;
    }

    if (edep_event.IsInteraction()) ++ninteractions;

    //Start building reweightable trajectory
    G4ReweightTraj primary_traj(
        edep_event.fPrimaryID,
        edep_event.fPrimaryPDG,
        -1, i, {0,0});
    for (auto & child_traj : edep_event.fPrimaryChildren) {
      primary_traj.AddChild(
        new G4ReweightTraj(child_traj.fID, child_traj.fPDG,
                           edep_event.fPrimaryID, i, {0, 0})
      );
    }

    //Build the path it took
    const auto & pt0 = edep_event.fPrimaryTraj->Points[0];
    //std::cout << 0 << " " << pt0.GetPosition().X()/10. << " " <<
    //             pt0.GetPosition().Y()/10. << " " << pt0.GetPosition().Z()/10. <<
    //             " " << pt0.GetProcess() << " " << pt0.GetSubprocess() <<
    //             " " << sqrt(pt0.GetMomentum()[0]*pt0.GetMomentum()[0] +
    //                         pt0.GetMomentum()[1]*pt0.GetMomentum()[1] +
    //                         pt0.GetMomentum()[2]*pt0.GetMomentum()[2]) <<
    //             std::endl;

    //bool prev_is_LAr = pt0.GetMaterial() == "LAr";

    out_len = 0.;
    for (size_t j = 1; j < edep_event.fPrimaryTraj->Points.size(); ++j) {
      std::string proc = "default";
      if ((j == edep_event.fPrimaryTraj->Points.size()-1) &&
          edep_event.IsInteraction()) {//Last step
        proc = "pi+Inelastic";
      }
      const auto & pt = edep_event.fPrimaryTraj->Points[j];

      std::string material = pt.GetMaterial();
      bool is_LAr = pt.GetMaterial() == "LAr";
      bool is_FR4 = pt.GetMaterial() == "FR4";
      if (proc == "pi+Inelastic" && is_LAr) interacted_in_lar = true;
      if (proc == "pi+Inelastic" && is_FR4) interacted_in_fr4 = true;

      const auto & prev_pt = edep_event.fPrimaryTraj->Points[j-1];
      auto dist = (pt.GetPosition() - prev_pt.GetPosition());
      double len = sqrt(dist.X()*dist.X() + dist.Y()*dist.Y() +
                        dist.Z()*dist.Z())/10.;
      out_len += len;
      double preStepP[3] = {
        prev_pt.GetMomentum()[0],
        prev_pt.GetMomentum()[1],
        prev_pt.GetMomentum()[2]
      };
      double postStepP[3] = {
        pt.GetMomentum()[0],
        pt.GetMomentum()[1],
        pt.GetMomentum()[2]
      };

      double total_pre_p = sqrt(preStepP[0]*preStepP[0] +
                                preStepP[1]*preStepP[1] +
                                preStepP[2]*preStepP[2]);

      double total_post_p = sqrt(postStepP[0]*postStepP[0] +
                                postStepP[1]*postStepP[1] +
                                postStepP[2]*postStepP[2]);

      std::cout << j << " " << pt.GetPosition().X()/10. << " " <<
                   pt.GetPosition().Y()/10. << " " <<
                   pt.GetPosition().Z()/10. << " " << pt.GetProcess() <<
                   " " << pt.GetSubprocess() << " " <<
                   len <<
                   " " << total_pre_p << " " << total_post_p << " " <<
                   pt.GetMaterial() << std::endl;

      G4ReweightStep * step = new G4ReweightStep(edep_event.fPrimaryID,
                                                 edep_event.fPrimaryPDG,
                                                 0, i, preStepP, postStepP,
                                                 len, proc);
      primary_traj.AddStep(step);
    }

    auto & end_pt
        = edep_event.fPrimaryTraj->Points[
            edep_event.fPrimaryTraj->Points.size()-1];
    end_x = end_pt.GetPosition().X()/10.;

    
   // std::vector<std::vector<size_t>> segments {};
   // bool prev_is_LAr = (pt0.GetMaterial() == "LAr");
   // std::cout << "0 " << pt0.GetMaterial() << std::endl;
   // for (size_t j = 1; j < edep_event.fPrimaryTraj->Points.size(); ++j) {
   //   const auto & pt = edep_event.fPrimaryTraj->Points[j];
   //   std::cout << j << " " << pt.GetMaterial() << std::endl;
   //   bool is_LAr = pt.GetMaterial() == "LAr";
   //   if (is_LAr) {
   //     if (!prev_is_LAr || (segments.size() == 0)) {
   //       segments.push_back({j-1});
   //     }
   //     else {
   //       segments.back().push_back(j-1);
   //     }
   //   }
   //   else {
   //     if (prev_is_LAr) {
   //       segments.back().push_back(j-1);
   //     }
   //   }

   //   if (j == edep_event.fPrimaryTraj->Points.size()-1 && prev_is_LAr) {
   //     segments.back().push_back(j);
   //   }
   //   prev_is_LAr = is_LAr;
   // }

   // for (auto & seg : segments) {
   //   if (seg.size() < 2) {
   //     throw cet::exception("do_weighting") << "Found segment < 2 size " << seg.size();
   //   }
   // }

    //Split up the traj by contiguous weightable material (LAr)
    G4ReweightTraj traj(
        edep_event.fPrimaryID,
        edep_event.fPrimaryPDG,
        -1, i, {0,0});
    for (auto & child_traj : edep_event.fPrimaryChildren) {
      traj.AddChild(
        new G4ReweightTraj(child_traj.fID, child_traj.fPDG,
                           edep_event.fPrimaryID, i, {0, 0})
      );
    }

    size_t npts = edep_event.fPrimaryTraj->Points.size();
    for (size_t j = 1; j < npts; ++j) {
      const auto & pt = edep_event.fPrimaryTraj->Points[j];
      if (pt.GetMaterial() != "LAr") continue;
      const auto & prev_pt = edep_event.fPrimaryTraj->Points[j-1];
      auto dist = (pt.GetPosition() - prev_pt.GetPosition());
      double len = sqrt(dist.X()*dist.X() + dist.Y()*dist.Y() +
                        dist.Z()*dist.Z())/10.;
      double preStepP[3] = {
        prev_pt.GetMomentum()[0],
        prev_pt.GetMomentum()[1],
        prev_pt.GetMomentum()[2]
      };
      double postStepP[3] = {
        pt.GetMomentum()[0],
        pt.GetMomentum()[1],
        pt.GetMomentum()[2]
      };

      double total_pre_p = sqrt(preStepP[0]*preStepP[0] +
                                preStepP[1]*preStepP[1] +
                                preStepP[2]*preStepP[2]);

      double total_post_p = sqrt(postStepP[0]*postStepP[0] +
                                postStepP[1]*postStepP[1] +
                                postStepP[2]*postStepP[2]);

      std::string proc = "default";
      if ((j == npts-1) && edep_event.IsInteraction()) {//Last step
        proc = "pi+Inelastic";
      }
      G4ReweightStep * step = new G4ReweightStep(edep_event.fPrimaryID,
                                                 edep_event.fPrimaryPDG,
                                                 0, i, preStepP, postStepP,
                                                 len, proc);
      traj.AddStep(step);
      std::cout << "Added step " << j-1 << " " << j << " " << proc << std::endl;
    }
    std::cout << std::endl;

    //Get the weight for this trajectory/event
    weight = reweighter->GetWeightFromSetParameters(traj);
    //for (auto & traj : trajs) {
    //  weight *= reweighter->GetWeightFromSetParameters(*traj);
    //}
    old_weight = reweighter->GetWeightFromSetParameters(primary_traj);
    std::cout << "Weight: " << weight << " " << old_weight << std::endl;;
    if (edep_event.IsInteraction()) nweighted += weight;
    ntotal += weight;


    outtree.Fill();
  }

  std::cout << std::endl;
  std::cout << "Interaction rate: " << 100.*ninteractions/tree->GetEntries() << std::endl;
  std::cout << "Weighted rate: " << 100.*nweighted/tree->GetEntries() << std::endl;
  std::cout << "Total weights: " << ntotal << std::endl;

  outtree.Write();
  fOut.Close();

  f.Close();

  return 0;
}
