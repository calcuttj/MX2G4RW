#setup edepsim v3_2_0b -q e20:prof
#setup geant4reweight v01_18_01 -q e20:prof:s120a
#setup dune_pardata v01_84_00

setup geant4reweight v01_20_08 -q e26:prof:s131
setup cmake v3_27_4
setup dune_pardata  v01_84_00


##For building and link against edepsim
export EDepSim_INC=/exp/dune/app/users/calcuttj/edepsim/edep-sim/install/include
export EDepSim_LIB=/exp/dune/app/users/calcuttj/edepsim/edep-sim/install/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${EDepSim_LIB} ### Why did I do this backward?
export PATH=/exp/dune/app/users/calcuttj/edepsim/edep-sim/install/bin/:$PATH
export CMAKE_PREFIX_PATH=$EDepSim_LIB/cmake/:${CMAKE_PREFIX_PATH}

##Setting up Mx2 software
export PATH=/exp/dune/app/users/calcuttj/mx2_g4rw_dev/cmake_build/install/bin/:$PATH
export LD_LIBRARY_PATH=/exp/dune/app/users/calcuttj/mx2_g4rw_dev/cmake_build/install/lib:$LD_LIBRARY_PATH
#export FHICL_FILE_PATH=
