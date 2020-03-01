# Mastersetpart Detector {#det-mastersetpart}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

Brief description missing.

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | mastersetpart               | ✓ |   |   |


### Details

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
