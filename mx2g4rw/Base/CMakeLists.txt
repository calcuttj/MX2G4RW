#include_directories()

add_library(
  mx2g4rw_base_lib
  SHARED
  Utilities.cxx
  ReweightSuite.cxx
)

target_include_directories(
  mx2g4rw_base_lib
  PUBLIC
  $ENV{FHICLCPP_INC}
  $ENV{CETLIB_INC}
  $ENV{GEANT4REWEIGHT_INC}
  $ENV{CETLIB_EXCEPT_INC}
  $ENV{BOOST_INC}
  #INTERFACE

)

find_library(REWEIGHTLIB ReweightBaseLib REQUIRED)
find_library(FHICLLIB fhiclcpp REQUIRED)
find_library(CETLIB cetlib REQUIRED)
find_library(CETLIB_EXCEPT cetlib_except REQUIRED)
find_library(ROOTCORE Core PATHS ${ROOTSYS}/lib/ REQUIRED)
#message(STATUS "edepsim dir: ${EDEPSIM}")
set(EDEPSIMLIBS $ENV{EDepSim_LIB})
#find_library(EDEPSIM E
find_library(
  EDEPSIMIO
  edepsim_io
  PATHS
  ${EDEPSIMLIBS}
  REQUIRED
)
find_library(BOOSTOPTS boost_program_options REQUIRED)

#message($

#target_include_directories()

target_link_libraries(
  mx2g4rw_base_lib
  PUBLIC
  ${REWEIGHTLIB}
  ${FHICLLIB}
  ${CETLIB}
  ${CETLIB_EXCEPT}
  ${ROOTCORE}
  ${EDEPSIMIO}
  ${BOOSTOPTS}
)

install(TARGETS mx2g4rw_base_lib) #"Utilities.h" "ReweightSuite.h")
