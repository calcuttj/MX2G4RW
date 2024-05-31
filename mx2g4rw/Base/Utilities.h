#ifndef Utilities_h
#define Utilities_h

#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"

#include "TFile.h"

namespace mx2g4rw {
fhicl::ParameterSet MakeFCLPars(std::string & fcl_file);
TFile * OpenFile(const std::string & fracs_file);
}
#endif
