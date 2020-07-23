# Detection Process {#detection-process}
# The GCG Detection Process
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

Generally, there are two options of how the detection can be conducted.
1. GCG detects structures in your program **without you** giving additional decomposition information.
2. GCG uses a decomposition file **given by you** as described in the @ref u5 "corresponding Use Case", which can be @ref own-dec-file "created by yourself" or [exported with GCG](FAQ.html#createsettingsfile). Then, this decomposition can be
  - **complete**, i.e. assigning every variable and constraint to a block (or the master) or
  - **partial**, i.e. assigns only a subset of variables and constraints.

If your decomposition is complete, no detector can contribute anything anymore.
Thus, GCG will not start @ref loop, but instead use the decomposition you gave
and immediately start the solving.  
If your decomposition is just partial, it is added to the (previously empty) pool of 
partial decompositions to use and GCG starts @ref loop. The procedure is clarified
in the following section.  

## The Detection Loop {#loop}
All detectors implemented in GCG assign variables/constraints to a block or the master.
When the detection loop is started, a pool with partial decompositions is initialized,
first with only one empty decomposition. 
Then, we start the loop. During this detection loop, we iterate over pairs of 
partial decompositions and detectors. This means that for each partial decomposition
a detector yields (and the initial empty one), all other detectors will 
(potentially) work on this decomposition, completing it more and more with each
iteration of the loop. It can happen that only one detector already completes a decomposition.
After a number of iterations, the decomposition is complete and can be rated according
to a score (see @subpage detection-scores for a list of scores that can be set to be used). Then,
the decomposition with the best score is used to proceed to the solving.

### Before Presolving
The instance (the problem) you read into GCG is usually not presolved yet (see @ref presolving).
Then there are three different cases that can happen:
- You don't give any decomposition.
- You have a rough idea of how your problem should be decomposed and give a partial decomposition.
- You know exactly how to decompose and give a complete decomposition.

Depending on those three cases, the detection runs differently. The exact process can be found
in Figure 1.

<img src="non-presolved.jpg" alt="Detection Process: non-presolved" style="display: block;  margin-left: auto;  margin-right: auto;  width: 90%;">
<div style="text-align: center;"><b>Fig. 1:</b> First stage of detection (problem is not yet presolved)</div>

### During/After Presolving
After the process shown in Figure 1 was executed, the problem you gave will be presolved.
After presolving, some variables might have been aggregated or even removed and thus, the
decompositions we found in the first run will be translated to match with the new, semantically equal,
yet different problem. If you already knew how your problem would look like after presolving,
you can also give a decomposition file for your problem - but presolved.\n
The detection process will then be executed as shown in Figure 2.

<img src="presolved.jpg" alt="Detection Process: presolved" style="display: block;  margin-left: auto;  margin-right: auto;  width: 90%;">
<div style="text-align: center;"><b>Fig. 2:</b> Nonzero entries after decomposition</div>


## Standalone Capability of the Loop
Todo:
- write about "standalone" capability of the detection loop
