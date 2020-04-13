# Constraint Class Detector {#det-consclass}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | consclass                   | ✓ |   |   |

This detector uses the classes identified by the **constraint classifiers**.

### Algorithmic Details
The Consclass Detector creates partial decompositions by fixing some constraints to the master problem and leaving the remaining constraints open (for finishing). This is done as follows:
 * Iterate over all used classifiers
 * For each classifier, iterate over all possible subsets of classes (with `ConsClassDecompInfo == BOTH`). For each subset the following partial decomposition is created:
   * The constraints beloning to classes in the current subset are fixed to the master problem
   * All constraints from classes with `ConsClassDecompInfo == ONLY_MASTER` are fixed to the master problem
   * All other constraints are left open

This resuls in all possible partial decompositions where the constraints of one class are either all in the master problem or all in pricing problems.

NOTE: In order to get sensible results with this detector a finishing detector needs to be used. Otherwise it results in only one decomposition (where all constraints are in the master problem). This is because all other partial decompositions are not finished, because the constraints which are not fixed to the master, are left open and are not assigned to blocks.

### Theoretical Details

### Parameters

    enabled               flag to indicate whether detector <consclass> is enabled [TRUE]
    finishingenabled      flag to indicate whether detector <consclass> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <consclass> [1]
    maxcallround          maximum round the detector gets called in detection loop <consclass> [0]
    maxnclasses           maximum number of classes  [5]
    mincallround          minimum round the detector gets called in detection loop <consclass> [0]
    origenabled           flag to indicate whether detector <consclass> is enabled for detecting in the original problem [TRUE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <consclass> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <consclass> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <consclass> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <consclass> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <consclass> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <consclass> [0]
    skip                  flag to indicate whether detector <consclass> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <consclass> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_consclass.cpp
