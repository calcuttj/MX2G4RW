# ======================================================================
#  mx2g4rw main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/mx2g4rw/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================

cmake_minimum_required (VERSION 3.20 FATAL_ERROR)

project(mx2g4rw VERSION 00.00.00)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(geant4reweight REQUIRED )
find_package(cetlib REQUIRED )
find_package(cetlib_except REQUIRED )
find_package(fhiclcpp REQUIRED )

find_package(EDepSim REQUIRED)
# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(mx2g4rw)
