#ifndef ReweightSuite_h
#define ReweightSuite_h

#include "TFile.h"

#include "geant4reweight/src/ReweightBase/G4ReweightManager.hh"
#include "geant4reweight/src/ReweightBase/G4MultiReweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightStep.hh"

#include "fhiclcpp/ParameterSet.h"

namespace mx2g4rw {

using PartMat = std::pair<int, std::string>;
using ParameterMap = std::map<PartMat, std::vector<fhicl::ParameterSet>>;
using ReweighterMap = std::map<PartMat, G4MultiReweighter*>;
using FracsFileMap = std::map<PartMat, TFile*>;

class ReweightSuite {
 public:
  ReweightSuite(fhicl::ParameterSet & pset);
  ~ReweightSuite();
 private:
  G4ReweightManager * fManager;
  ReweighterMap fReweighters;
  FracsFileMap fFracsFiles;
  ParameterMap fParameters;
};
}

#endif
