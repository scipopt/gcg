# Generalmastersetpack Detector {#det-generalmastersetpack}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

Brief description missing.

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | generalmastersetpack        | ✓ |   |   |


### Details

### Parameters

    enabled               flag to indicate whether detector <generalmastersetpack> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <generalmastersetpack> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <generalmastersetpack> [1]
    maxcallround          maximum round the detector gets called in detection loop <generalmastersetpack> [0]
    mincallround          minimum round the detector gets called in detection loop <generalmastersetpack> [0]
    origenabled           flag to indicate whether detector <generalmastersetpack> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <generalmastersetpack> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <generalmastersetpack> [0]
    origmincallround      minimum round the detector gets called in detection loop <generalmastersetpack> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <generalmastersetpack> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <generalmastersetpack> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <generalmastersetpack> [0]
    skip                  flag to indicate whether detector <generalmastersetpack> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <generalmastersetpack> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_generalmastersetpack.cpp
