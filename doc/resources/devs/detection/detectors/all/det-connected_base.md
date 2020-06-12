# Connected Detector (base) {#det-connected_base}

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| `C` | connectedbase               |   | ✓ |   |

This detector assigns all constraints and variables to the same block if they are connected,
A constraint and a variable are adjacent if the variable appears in the constraint.

### Algorithmic Details
This detector is only different from "connected_nonewlinkingvars", if the parameter `useconssadj` is true and the function `GCGconshdlrDecompGetConssAdjCalculated(scip)` returns true and the number of linking variables is 0. Else, the `completeByConnected()` function is called, invoking the same procedure as the "connected_nonewlinkingvars" detector.

* If the number of blocks in the given partial decomposition is below 0, set it to 0.
* Create a boolean vector with number of constraints as length. Vector value is true iff the constraint is open.
* While there are open constraints remaining (initialize queue with the first open constraint)
  * add the current constraint to the neighborConss vector
  * Iterate over open constraints that are adjacent to the current one according to the Seeedpool::createConssAdjacency() function
    * if this constraint is not open or in the master constraints or was already visited, continue
    * else add this constraint to the neighborConss vector
  * get neighborVars according to the calculated neighborConss
  * define a new block with the neighborConss and neighborVars
  * remove all neighborConss and neighborVars from the open variables
* assign all remaining variables to block 0, if it exists, and to master otherwise

### Theoretical Details

### Parameters

    enabled               flag to indicate whether detector <connectedbase> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <connectedbase> is enabled for finishing of incomplete decompositions [TRUE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <connectedbase> [1]
    maxcallround          maximum round the detector gets called in detection loop <connectedbase> [2147483647]
    mincallround          minimum round the detector gets called in detection loop <connectedbase> [0]
    origenabled           flag to indicate whether detector <connectedbase> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <connectedbase> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <connectedbase> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <connectedbase> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <connectedbase> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <connectedbase> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <connectedbase> [0]
    skip                  flag to indicate whether detector <connectedbase> should be skipped if others found decompositions [FALSE]
    useconssadj           should the constraint adjacency be used [TRUE]
    usefullrecall         flag to indicate whether detector <connectedbase> should be called on descendants of the current partialdec [FALSE]



### Links
 * Documentation: dec_connectedbase.cpp
