# MX2G4RW

Repository to start developing utils for using Geant4Reweight for studies at the MINERvA + 2x2 ND prototype.

# Setting up
This code can be ran easily on a dunegpvm. The required packages are 
<pre>
  dune_pardata
  geant4reweight
  edepsim
  boost
</pre>
which can be set up by sourcing the provided `setup.sh` script.
Note: this includes a modified edepsim which adds some necessary information allowing for use with geant4reweight. It currently is hardcoded to a path under the /exp/dune/app area

TODO -- update how to install!

# Configuration
The `reweight.fcl` file contains the configuration used to run the weighting code in `new_weighting.c++`.

# Running weighting/processing code
The processing can be run on a file produced by edep-sim.
<pre>
  new_weighting -c reweight.fcl -i [input].root -o [outputname].root
</pre>


# Output Format
The output of the weighting routine is a file with a single weight table stored as a TTree. Each entry in the TTree/weighting table corresponds to an event simulated within edep-sim. The weights are separated in several branches (determined by configuration) and are separated according to

1) The particle species in question (i.e. pi+, proton)
2) The target material in which we are weighting (i.e. LAr, FR4)
3) The 'parameter' of interest (a component of the total cross section within some defined momentum range).

The branch names are of the form "{particle}\_{material}\_{parameter}\_weights", and is a 2D vector of doubles in which the first index corresponds to a trajectory within the edep-sim event, and the second corresponds to several parameter values over which we have "scanned" the weighting routines. Note: only the trajectories which match the particle species and which travelled through the target material during edep-sim are saved in each vector. That way, weights for protons will not be saved in a branch representing a change to some pion cross section, nor will weights be saved for a pi+ trajectory which travelled only through FR4 when the target material was LAr. 

In addition to the 2D weight branch is a branch holding a vector of the relevant trajectory IDs corresponding to a given weight branch. It will have the name "{particle}\_{material}\_{parameter}\_IDs" and the particle, material, and parameter will match some weight branch. The first index of the ID vector will correspond to the first index of the weight vector event-by-event.

# Demonstration of Applying Weights
