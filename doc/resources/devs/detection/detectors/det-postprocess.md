# Postprocess Detector {#det-postprocess}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

This detector is our default postprocessing detector. It is called to add all unassigned constraints/variables to the master.

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | postprocess                 |   |   | ✓ |

### Details

### Parameters

    enabled               flag to indicate whether detector <postprocess> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <postprocess> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <postprocess> [1]
    maxcallround          maximum round the detector gets called in detection loop <postprocess> [2147483647]
    mincallround          minimum round the detector gets called in detection loop <postprocess> [0]
    origenabled           flag to indicate whether detector <postprocess> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <postprocess> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <postprocess> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <postprocess> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <postprocess> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <postprocess> is enabled for postprocessing of finished decompositions [TRUE]
    priority              priority of detector <postprocess> [1000000]
    skip                  flag to indicate whether detector <postprocess> should be skipped if others found decompositions [FALSE]
    useconssadj           should the constraint adjacency be used [TRUE]
    usefullrecall         flag to indicate whether detector <postprocess> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_postprocess.cpp
