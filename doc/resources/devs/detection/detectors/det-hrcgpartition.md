# Hypergraph Partitioning Detector (Columns & Rows) {#det-hrcgpartition}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | hrcgpartition               | ✓ | ✓ |   |


This detector detects [arrowhead](#arrowhead) structures through graph partitioning on **column-row hypergraphs**.

**Note:** This detector needs hmetis and works only under Linux/MacOS, it further needs the Z-shell (zsh)
to enforce memory and time limits on hmetis as this is the only shell reliably doing that.

### Details

### Parameters

    alpha                 Factor on how heavy the standard deviation of the coefficients is measured [0.0]
    beta                  Factor on how heavy equality (beta) and inequality constraints are measured [0.5]
    consWeight            Weight of a constraint hyperedge [1]
    consWeightSetppc      Weight for constraint hyperedges that are setpartitioning or covering constraints [5]
    dummynodes            Percentage of dummy nodes for metis [0.2]
    enabled               flag to indicate whether detector <hrcgpartition> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <hrcgpartition> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <hrcgpartition> [1]
    maxblocks             The maximal number of blocks (detector is called for all block numbers in [minblocks,maxblocks]) [20]
    maxcallround          maximum round the detector gets called in detection loop <hrcgpartition> [1]
    maxnblockcandidates   The maximal number of block number candidates [3]
    metisuseptyperb       Should the rb or kway method be used for partitioning by metis [TRUE]
    metisverbose          Should the metis output be displayed [FALSE]
    minblocks             The minimal number of blocks (detector is called for all block numbers in [minblocks,maxblocks]) [2]
    mincallround          minimum round the detector gets called in detection loop <hrcgpartition> [0]
    origenabled           flag to indicate whether detector <hrcgpartition> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <hrcgpartition> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <hrcgpartition> [1]
    origmincallround      minimum round the detector gets called in detection loop <hrcgpartition> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <hrcgpartition> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <hrcgpartition> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <hrcgpartition> [1000]
    randomseed            Random seed for hmetis [1]
    realname              Should the problem be used for metis files or a temporary name [FALSE]
    skip                  flag to indicate whether detector <hrcgpartition> should be skipped if others found decompositions [FALSE]
    tidy                  Whether to clean up temporary files [TRUE]
    ubfactor              Unbalance factor for metis [5.0]
    usefullrecall         flag to indicate whether detector <hrcgpartition> should be called on descendants of the current partialdec [TRUE]
    varWeight             Weight of a variable hyperedge [2]
    varWeightBinary       Weight of a binary variable hyperedge [3]
    varWeightContinous    Weight of a continuos variable hyperedge [2]
    varWeightImplint      Weight of a implicit integer variable hyperedge [3]
    varWeightInteger      Weight of a integer variable hyperedge [3]


### Links
 * Documentation: dec_hrcgpartition.cpp
