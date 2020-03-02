# Hrgpartition Detector {#det-hrgpartition}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | hrgpartition                | ✓ |   |   |


This detector detects [arrowhead](#arrowhead) structures as well as structures
with only linking variables or linking constraints.

**Note:** This detector needs hmetis and works only under Linux/MacOS, it further needs the Z-shell (zsh)
to enforce memory and time limits on hmetis as this is the only shell reliably doing that.

### Details

### Parameters

    advanced/enabledforlargeproblems
                     -->  Should this detector be enabled even the limit nconssnvars is exceeded [FALSE]
    advanced/limitnconssnvarsdefault
                     -->  Limit for sum of nvars and nconss for enabling this detector in default [10000]
    alpha                 Factor on how heavy the standard deviation of the coefficients is measured [0.0]
    beta                  Factor on how heavy equality (beta) and inequality constraints are measured [0.5]
    consWeight            Weight of a constraint hyperedge [5]
    consWeightSetppc      Weight for constraint hyperedges that are setpartitioning or covering constraints [5]
    dummynodes            Percentage of dummy nodes for metis [0.2]
    enabled               flag to indicate whether detector <hrgpartition> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <hrgpartition> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <hrgpartition> [1]
    maxblocks             The maximal number of blocks (detector is called for all block numbers in [minblocks,maxblocks]) [20]
    maxcallround          maximum round the detector gets called in detection loop <hrgpartition> [0]
    maxnblockcandidates   The maximal number of block number candidates [3]
    metisuseptyperb       Should the rb or kway method be used for partitioning by metis [TRUE]
    metisverbose          Should the metis output be displayed [FALSE]
    minblocks             The minimal number of blocks (detector is called for all block numbers in [minblocks,maxblocks]) [2]
    mincallround          minimum round the detector gets called in detection loop <hrgpartition> [0]
    origenabled           flag to indicate whether detector <hrgpartition> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <hrgpartition> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <hrgpartition> [0]
    origmincallround      minimum round the detector gets called in detection loop <hrgpartition> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <hrgpartition> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <hrgpartition> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <hrgpartition> [1000]
    randomseed            Random seed for hmetis [1]
    realname              Should the problem be used for metis files or a temporary name [FALSE]
    skip                  flag to indicate whether detector <hrgpartition> should be skipped if others found decompositions [FALSE]
    tidy                  Whether to clean up temporary files [TRUE]
    ubfactor              Unbalance factor for metis [5.0]
    usefullrecall         flag to indicate whether detector <hrgpartition> should be called on descendants of the current partialdec [TRUE]
    varWeight             Weight of a variable hyperedge [1]
    varWeightBinary       Weight of a binary variable hyperedge [2]
    varWeightContinous    Weight of a continuos variable hyperedge [1]
    varWeightImplint      Weight of a implicit integer variable hyperedge [2]
    varWeightInteger      Weight of a integer variable hyperedge [2]


### Links
 * Documentation: dec_hrgpartition.cpp
