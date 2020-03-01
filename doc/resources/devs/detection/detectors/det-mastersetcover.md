# Mastersetcover Detector {#det-mastersetcover}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

Brief description missing.

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | mastersetcover              | ✓ |   |   |

### Details

### Parameters

    enabled               flag to indicate whether detector <mastersetcover> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <mastersetcover> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <mastersetcover> [1]
    maxcallround          maximum round the detector gets called in detection loop <mastersetcover> [2147483647]
    mincallround          minimum round the detector gets called in detection loop <mastersetcover> [0]
    origenabled           flag to indicate whether detector <mastersetcover> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <mastersetcover> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <mastersetcover> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <mastersetcover> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <mastersetcover> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <mastersetcover> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <mastersetcover> [0]
    skip                  flag to indicate whether detector <mastersetcover> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <mastersetcover> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_mastersetcover.cpp
