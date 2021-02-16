# Set Partitioning Detector {#det-mastersetpart}

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| `?` | mastersetpart               | ✓ |   |   |

This detector sets the following constraint types as master constraints:
- set partitioning constraints

### Algorithmic Details
This detector adds the aforementioned constraints to the master. This is done as follows:
* Iterate over all open constraints
 * If the constraint's type (determined in [preprocessing](#preprocessing)) is `setpartitioning`, fix it to the master

### Theoretical Details

### Parameters

    enabled               flag to indicate whether detector <mastersetpart> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <mastersetpart> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <mastersetpart> [1]
    maxcallround          maximum round the detector gets called in detection loop <mastersetpart> [2147483647]
    mincallround          minimum round the detector gets called in detection loop <mastersetpart> [0]
    origenabled           flag to indicate whether detector <mastersetpart> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <mastersetpart> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <mastersetpart> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <mastersetpart> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <mastersetpart> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <mastersetpart> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <mastersetpart> [0]
    skip                  flag to indicate whether detector <mastersetpart> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <mastersetpart> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_mastersetpart.cpp
