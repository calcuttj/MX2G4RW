# MX2G4RW

Repository to start developing utils for using Geant4Reweight for studies at the MINERvA + 2x2 ND prototype.I'm ready whenever you are to 

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
For this demonstration, several 'observables' branches have been added to the output tree. It uses an edep-sim simulation of a pion gun generated in a single direction within a homogenous block of LAr.

## Primary Pion Length
Below is a snippet to draw distributions of the primary pion length using the branches stored in the weight table "tree". Note that the specific access of the primary pions using the 0th index only works because of the specific simulation produced by a pion gun (thus the primary pion is defined to be the 0th trajectory in the tree). The weighting routines were configured to vary the total piplus cross section in LAr by a factor of 0.1 -- 2.0 (separated by 0.1). Thus, the 0th element in the weight vector's 2nd dimension coincides to a factor of 0.1, element 4 is 0.5, element 9 is 1.0, element 14 is 1.5, etc.

<pre>
  //Draw the length (divide by 10 to fix units mm --> cm
  //Use the 1st entry in each vector to correspond to the primary traj
  //  
  //NOTE: THIS IS ONLY FOR THIS CASE AND IS ONLY BECAUSE I GENERATED PIONS
  //SO BE CAREFUL WHEN YOU MAKE YOUR OWN APPLICATION

  //Unweighted
  tree->Draw("len[0]/10.>>hNom", "");

  //Weighted up by 50% -- 2nd index corresponds to parameter value
  tree->Draw("len[0]/10.>>hUp", "piplus_LAr_fReac_weights[0][14]");

  //Weighted down by 50%
  tree->Draw("len[0]/10.>>hDown", "piplus_LAr_fReac_weights[0][4]");
</pre>

This produces the following image. Several features are evident: 1) by increasing the cross section (red), the pions tend to interact earlier in their traversal through the LAr. The opposite is true when decreasing the cross section. 2) The normalization remains the same.* This reflects the fact that we are not changing the number of primary particles in our simulation, we are simply changing how they behave. This is an important concept to start thinking about, especially when considering more complicated demonstrations. When we start to weight over multiple trajectories in the event, a given distribution's normalization could in fact change when weighting because we might increase or decrease the chance for particles to be created by the model we are varying (as opposed to the particle gun). 

\* This is true to within a few percent. Any inaccuracies can be caused by i.e. the post-simulation merging of steps in the edep-sim output. 

![example_len](docs/images/example_len.png)



## Number of Pions
This time, I want to draw the number of pions produced within the event. In order to encompass the effect of changing the cross section on all of the pions in the event, one must multiply the weights from every trajectory together. 

<pre >
  f = RT.TFile.Open(args.i);
  t = f.Get("tree")
  
  #Make some hists for storage
  hNom = RT.TH1D("hNom", "", 15, 0, 15) 
  hUp = RT.TH1D("hUp", "", 15, 0, 15) 
  hDown = RT.TH1D("hDown", "", 15, 0, 15);

  for e in t:
    #For each entry, we have to loop over all of the pions in the event 
    #and combine all of their weights
    weight_up = 1.
    weight_down = 1.
    for wv in e.piplus_LAr_fReac_weights:
      weight_up *= wv[14]
      weight_down *= wv[4]

    #Count the number of piplus 
    npions = [i for i in e.pdg].count(211)

    #Weight accordingly
    hUp.Fill(npions, weight_up);
    hNom.Fill(npions);
    hDown.Fill(npions, weight_down);

</pre>

![example_pdg](docs/images/npions.png) ![example_pdg_log](docs/images/npions_log.png)
