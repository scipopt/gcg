# Connected Detector (nonewlinkingvars) {#det-connected_nonewlinkingvars}

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| `?` | connected_noNewLinkingVars  | ✓ | ✓ |   |

This detector assigns all dependent open constraints and variables and completes the partial decomposition by breadth-first search.

### Algorithmic Details
* If the number of blocks in the given partial decomposition is below 0, set it to 0.
* Create two boolean vectors with number of cons/vars as length. Vector value is true iff the con/var is open.
* While there are open constraints remaining (initialize queue with the first open constraint)
  * add the current constraint to the neighborConss vector
  * Iterate over open variables in this constraint (all other variables are linking variables and are skipped)
    * Iterate over open constraints using this variable
      * add this constraint to the neighborConss vector
      * add the current constraint to the queue
    * add this variable to the neighborVars vector
  * define a new block with the neighborConss and neighborVars
  * remove all neighborConss and neighborVars from the open variables
* assign all remaining variables to block 0, if it exists, and to master otherwise

### Theoretical Details

### Parameters

    enabled               flag to indicate whether detector <connected_nonewlinkingvars> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <connected_nonewlinkingvars> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <connected_nonewlinkingvars> [1]
    maxcallround          maximum round the detector gets called in detection loop <connected_nonewlinkingvars> [2147483647]
    mincallround          minimum round the detector gets called in detection loop <connected_nonewlinkingvars> [0]
    origenabled           flag to indicate whether detector <connected_nonewlinkingvars> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <connected_nonewlinkingvars> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <connected_nonewlinkingvars> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <connected_nonewlinkingvars> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <connected_nonewlinkingvars> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <connected_nonewlinkingvars> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <connected_nonewlinkingvars> [0]
    skip                  flag to indicate whether detector <connected_nonewlinkingvars> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <connected_nonewlinkingvars> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_connected_noNewLinkingVars.cpp
