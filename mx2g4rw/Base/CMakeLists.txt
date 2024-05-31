#include_directories()

add_library(
  mx2g4rw_base_lib
  SHARED
  Utilities.cxx
  ReweightSuite.cxx
)

find_library(REWEIGHTLIB ReweightBaseLib REQUIRED)
find_library(FHICLLIB fhiclcpp REQUIRED)
find_library(CETLIB cetlib REQUIRED)
find_library(CETLIB_EXCEPT cetlib_except REQUIRED)
find_library(ROOTCORE Core PATHS ${ROOTSYS}/lib/ REQUIRED)
#message(STATUS "edepsim dir: ${EDEPSIM}")
set(EDEPSIMLIBS $ENV{EDepSim_LIB})
find_library(
  EDEPSIMIO
  edepsim_io
  PATHS
  ${EDEPSIMLIBS}
  REQUIRED
)
find_library(BOOSTOPTS boost_program_options REQUIRED)

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
