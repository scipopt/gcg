# General Set Covering Detector {#det-generalmastersetcover}
### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| `d` | generalmastersetcover       | ✓ |   |   |

This detector sets the following constraint types as master constraints:
- set covering constraints
- logical OR constraints
- constraints with \f$\text{rhs}=\infty\f$ and \f$\text{lhs}\geq 0\f$


### Algorithmic Details
This detector adds the aforementioned constraints to the master. This is done as follows:
* Iterate over all open constraints
 * If the constraint's type (determined in [preprocessing](#preprocessing)) is `setcovering` or `logicor`, fix it to the master
 * If the constraint's type is not `logicor` and not `setpacking` and not `setpartitioning`, but its right hand side is infinity and its left hand side is non-negative, fix it to the master

### Theoretical Details

### Parameters

    enabled               flag to indicate whether detector <generalmastersetcover> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <generalmastersetcover> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <generalmastersetcover> [1]
    maxcallround          maximum round the detector gets called in detection loop <generalmastersetcover> [0]
    mincallround          minimum round the detector gets called in detection loop <generalmastersetcover> [0]
    origenabled           flag to indicate whether detector <generalmastersetcover> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <generalmastersetcover> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <generalmastersetcover> [0]
    origmincallround      minimum round the detector gets called in detection loop <generalmastersetcover> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <generalmastersetcover> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <generalmastersetcover> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <generalmastersetcover> [0]
    skip                  flag to indicate whether detector <generalmastersetcover> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <generalmastersetcover> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_generalmastersetcover.cpp
