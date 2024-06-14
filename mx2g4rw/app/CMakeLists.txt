
find_library(ROOTCORE Core PATHS ${ROOTSYS}/lib/ REQUIRED)
find_library(ROOTIO RIO PATHS ${ROOTSYS}/lib/ REQUIRED)
find_library(ROOTMATH MathCore PATHS ${ROOTSYS}/lib/ REQUIRED)
add_executable(
  new_weighting
  new_weighting.c++
)

message("source dir: ${CMAKE_SOURCE_DIR}")

target_include_directories(
  new_weighting
  PUBLIC
  $ENV{EDepSim_INC}
  ${CMAKE_SOURCE_DIR}
)
target_link_libraries(
  new_weighting
  mx2g4rw_base_lib
  ${ROOTCORE}
  RIO
  MathCore 
  Physics
  Tree
)

install(TARGETS new_weighting)

# TODO -- install FCL
