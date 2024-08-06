#ifndef ReweightSuite_h
#define ReweightSuite_h

#include "TFile.h"

#include "geant4reweight/src/ReweightBase/G4ReweightManager.hh"
#include "geant4reweight/src/ReweightBase/G4MultiReweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightStep.hh"

#include "fhiclcpp/ParameterSet.h"

namespace mx2g4rw {

using PartMat_t = std::pair<int, std::string>;
using ParameterMap_t = std::map<PartMat_t, std::vector<fhicl::ParameterSet>>;
using ReweighterMap_t = std::map<PartMat_t, G4MultiReweighter*>;
using FracsFileMap_t = std::map<PartMat_t, TFile*>;
using NParMap_t = std::map<PartMat_t, size_t>;

using ParNameMap_t = std::map<PartMat_t, std::vector<std::string>>;

class ReweightSuite {
 public:
  ReweightSuite(fhicl::ParameterSet & pset);
  ~ReweightSuite();

  bool CheckPDG(int pdg) const;
  bool CheckMaterial(const std::string & material) const;

  const std::vector<PartMat_t> & GetAllPartMats() const {
    return fPartMats;
  };

  const ParNameMap_t & GetParameterNames() const {
    return fParameterNames;
  };

  std::vector<double> Scan(
      G4ReweightTraj & traj,
      int pdg, std::string material, size_t param_number,
      size_t nsteps=20, double start=.1, double end=2.);
  std::vector<double> Scan(
      G4ReweightTraj & traj,
      PartMat_t part_mat, size_t param_number,
      size_t nsteps=20, double start=.1, double end=2.);
  size_t GetNPars(PartMat_t part_mat) const {
    //TODO -- throw exception
    return fNParameters.at(part_mat);
  };
 private:
  G4ReweightManager * fManager;
  ReweighterMap_t fReweighters;
  FracsFileMap_t fFracsFiles;
  ParameterMap_t fParameters;
  NParMap_t fNParameters;
  ParNameMap_t fParameterNames;
  std::vector<PartMat_t> fPartMats;
};
}

#endif
