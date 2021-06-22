# General Set Partitioning Detector {#det-generalmastersetpart}
### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| `?` | generalmastersetpart        | ✓ |   |   |

This detector sets the following constraint types as master constraints:
- set partitioning constraints
- constraints with \f$\text{rhs}\geq 0\f$ and \f$\text{rhs}\geq 0\f$ and \f$\text{rhs} = \text{lhs}\f$


### Algorithmic Details
This detector adds the aforementioned constraints to the master. This is done as follows:
* Iterate over all open constraints
 * If the constraint's type (determined in [preprocessing](#preprocessing)) is `setpartitioning`, fix it to the master
 * If the constraint's type is not `logicor` and not `setcovering` and not `setpacking`, but its right hand side is non-negative and its left hand side is non-negative, fix it to the master

### Theoretical Details

### Parameters

    enabled               flag to indicate whether detector <generalmastersetpart> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <generalmastersetpart> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <generalmastersetpart> [1]
    maxcallround          maximum round the detector gets called in detection loop <generalmastersetpart> [0]
    mincallround          minimum round the detector gets called in detection loop <generalmastersetpart> [0]
    origenabled           flag to indicate whether detector <generalmastersetpart> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <generalmastersetpart> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <generalmastersetpart> [0]
    origmincallround      minimum round the detector gets called in detection loop <generalmastersetpart> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <generalmastersetpart> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <generalmastersetpart> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <generalmastersetpart> [0]
    skip                  flag to indicate whether detector <generalmastersetpart> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <generalmastersetpart> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_generalmastersetpart.cpp
