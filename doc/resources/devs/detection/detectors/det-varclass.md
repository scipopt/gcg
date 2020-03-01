# Varclass Detector {#det-varclass}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

This detector uses the classes identified by the **variable classifiers**.

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | varclass                    | ✓ |   |   |

### Details

### Parameters

    enabled               flag to indicate whether detector <varclass> is enabled [TRUE]
    finishingenabled      flag to indicate whether detector <varclass> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <varclass> [1]
    maxcallround          maximum round the detector gets called in detection loop <varclass> [0]
    maxnclasses           maximum number of classes  [8]
    mincallround          minimum round the detector gets called in detection loop <varclass> [0]
    origenabled           flag to indicate whether detector <varclass> is enabled for detecting in the original problem [TRUE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <varclass> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <varclass> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <varclass> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <varclass> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <varclass> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <varclass> [0]
    skip                  flag to indicate whether detector <varclass> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <varclass> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_varclass.cpp
