#include "ReweightSuite.h"
#include "Utilities.h"
#include "cetlib/filepath_maker.h"
#include "cetlib/filesystem.h"
#include "cetlib/search_path.h"

namespace mx2g4rw{

using ParSetVector = std::vector<fhicl::ParameterSet>;

ReweightSuite::~ReweightSuite() {
  for (auto & rw : fReweighters) {
    delete rw.second;
  }
  for (auto & f : fFracsFiles) {
    f.second->Close();
    delete f.second;
  }
  delete fManager;
}

ReweightSuite::ReweightSuite(fhicl::ParameterSet & pset) {
  //Get the material list and make the manager
  auto materials = pset.get<ParSetVector>("Materials");
  fManager = new G4ReweightManager(materials);

  //fParameters = pset.get<std::vector<fhicl::ParameterSet>>("ParameterSet");
  auto all_reweights = pset.get<ParSetVector>("Reweights");

  for (auto & rw : all_reweights) {
    std::string name = rw.get<std::string>("Name");
    std::cout << "Making " << name << std::endl;
    
    int pdg = rw.get<int>("PDG");
    auto this_material = rw.get<fhicl::ParameterSet>("Material");
    auto mat_name = this_material.get<std::string>("Name");
    auto this_part_mat = std::make_pair(pdg, mat_name);

    auto parameter_set = rw.get<ParSetVector>("ParameterSet");

    fFracsFiles[this_part_mat] = OpenFile(mat_name);
    fReweighters[this_part_mat] = new G4MultiReweighter(
        pdg, *fFracsFiles[this_part_mat], parameter_set, this_material,
        fManager);
  }
}


}
