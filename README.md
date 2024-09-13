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
