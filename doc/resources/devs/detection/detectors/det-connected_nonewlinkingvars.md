# Connected Detector (nonewlinkingvars) {#det-connected_nonewlinkingvars}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

Brief description missing.

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | connected_noNewLinkingVars  | ✓ | ✓ |   |

### Details

### Parameters

    enabled               flag to indicate whether detector <connected_nonewlinkingvars> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <connected_nonewlinkingvars> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <connected_nonewlinkingvars> [1]
    maxcallround          maximum round the detector gets called in detection loop <connected_nonewlinkingvars> [2147483647]
    mincallround          minimum round the detector gets called in detection loop <connected_nonewlinkingvars> [0]
    origenabled           flag to indicate whether detector <connected_nonewlinkingvars> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <connected_nonewlinkingvars> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <connected_nonewlinkingvars> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <connected_nonewlinkingvars> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <connected_nonewlinkingvars> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <connected_nonewlinkingvars> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <connected_nonewlinkingvars> [0]
    skip                  flag to indicate whether detector <connected_nonewlinkingvars> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <connected_nonewlinkingvars> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_connected_noNewLinkingVars.cpp
