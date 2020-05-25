# Detection Process {#detection-process}
# The GCG Detection Process
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

Generally, there are two options of how the detection can be conducted.
1. GCG detects structures in your program **without you** giving additional decomposition information.
2. GCG uses a decomposition file **given by you** as described in the @ref u5 "corresponding Use Case", which can be @ref create-own-decfile "created by yourself" or [exported with GCG](FAQ.html#createsettingsfile). Then, this decomposition can be
  - **complete**, i.e. assigning every variable and constraint to a block (or the master) or
  - **partial**, i.e. assigns only a subset of variables and constraints.

If your decomposition is complete, no detector can contribute anything anymore.
Thus, GCG will not start @ref loop, but instead use the decomposition you gave
and immediately start the solving.  
If your decomposition is just partial, it is added to the (previously empty) pool of 
partial decompositions to use and GCG starts @ref loop. The procedure is clarified
in the following section.  

## The Detection Loop {#loop}
### Brief Description
All detectors implemented in GCG assign variables/constraints to a block or the master.
When the detection loop is started, a pool with partial decompositions is initialized,
first with only one empty decomposition. 
Then, we start the loop. During this detection loop, we iterate over pairs of 
partial decompositions and detectors. This means that for each partial decomposition
a detector yields (and the initial empty one), all other detectors will 
(potentially) work on this decomposition, completing it more and more with each
iteration of the loop. It can happen that only one detector already completes a decomposition.
After a number of iterations, the decomposition is complete and can be rated according
to a score (see @ref scores for a list of scores that can be set to be used). Then,
the decomposition with the best score is used to proceed to the solving.

@todo add figure

### Long Description

Todo:
- write about "standalone" capability of the detection loop
