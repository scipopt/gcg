# Densemasterconss Detector {#det-densemasterconss}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | densemasterconss            | ✓ |   |   |

Brief description missing.

### Details

### Parameters

    enabled               flag to indicate whether detector <densemasterconss> is enabled [TRUE]
    finishingenabled      flag to indicate whether detector <densemasterconss> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <densemasterconss> [1]
    maxcallround          maximum round the detector gets called in detection loop <densemasterconss> [0]
    mincallround          minimum round the detector gets called in detection loop <densemasterconss> [0]
    origenabled           flag to indicate whether detector <densemasterconss> is enabled for detecting in the original problem [TRUE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <densemasterconss> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <densemasterconss> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <densemasterconss> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <densemasterconss> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <densemasterconss> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <densemasterconss> [0]
    skip                  flag to indicate whether detector <densemasterconss> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <densemasterconss> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_densemasterconss.cpp
