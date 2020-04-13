# Connected Detector (base) {#det-connected_base}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | connectedbase               |   | ✓ |   |

This detector assigns all constraints and variables to the same block if they are connected,
A constraint and a variable are adjacent if the variable appears in the constraint.

### Algorithmic Details

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
