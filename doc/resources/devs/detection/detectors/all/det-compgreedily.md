# Complete Greedily Detector {#det-compgreedily}

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| `g` | compgreedily                | ✓ | ✓ |   |

This detector assigns the open cons and open vars of the partialdec greedily.

### Algorithmic Details
There are two different assignments performed. For both, i.e. constraints and variables, they are tried to be added to a fitting block first and if none could have been found, added to the master constraints.

* If there are no blocks (just master constraints) yet, create a new block and add all open constraints to this block.
Behavior for Variables:
* Iterate over all open variables. For each, reset the `varinblocks` vector and remove it from the open vars once it was assigned.
  * Iterate over all existing blocks and get their constraints.
    * Check if in any of these constraints, the variable occurs. If it does, add it to `varinblocks`. 
  * If the variable only occured in constraints of a single block, add it to this block
  * Else if the variable occured twice, i.e. in the constraints of two different blocks:
    * if those two blocks are subsequent: add it as a staircase linking variable
    * else add it as a usual linking variable
  * Else if the variable occured more than twice, add it as a usual linking variable
  * If the variable does not occur in an open constraint:
    * Check if it occurs in a master constraint. If it does, add it to the master problem.
  * If the variable couldn't be added to a block, assign it to the master, if it occurs in a master constraint. (It must after the loop for the constraints (see below) ran)

Behavior for constraints:
* Iterate over all open constraints.
  * Iterate over all existing blocks and check if all variables of the current constraint are in the current block, or linking variables, or master variables, or still open.
  * If so, add the constraint to the current block and add all variables that were open to this block as well.
  * If the constraint couldn't be added to any block, add it to the master problem.

### Theoretical Details

### Parameters

    enabled               flag to indicate whether detector <compgreedily> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <compgreedily> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <compgreedily> [1]
    maxcallround          maximum round the detector gets called in detection loop <compgreedily> [2147483647]
    mincallround          minimum round the detector gets called in detection loop <compgreedily> [0]
    origenabled           flag to indicate whether detector <compgreedily> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <compgreedily> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <compgreedily> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <compgreedily> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <compgreedily> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <compgreedily> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <compgreedily> [0]
    skip                  flag to indicate whether detector <compgreedily> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <compgreedily> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_compgreedily.cpp
