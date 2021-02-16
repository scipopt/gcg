# Staircase Detector (Matrix Reordering) {#det-stairheur}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| `s` | stairheur                   | ✓ |   |   |

This detector finds staircase structures via a modified version of a
 **matrix reordering** algorithm, Rank-Order-Clustering (ROC)[2].

### Algorithmic Details
The algorithmic details can be found in [2].

### Theoretical Details
The theoretical details can be found in [1].

### Parameters

    blockingassoonaspossible
                     -->  Enable blocking type 'as soon as possible [FALSE]
    desiredblocks         The desired number of blocks. 0 means automatic determination of the number of blocks. [0]
    dynamicblocking       Enable blocking type 'dynamic' [FALSE]
    enabled               flag to indicate whether detector <stairheur> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <stairheur> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <stairheur> [1]
    maxblocks             The maximal number of blocks [20]
    maxcallround          maximum round the detector gets called in detection loop <stairheur> [2147483647]
    maxiterationsROC      The maximum number of iterations of the ROC-algorithm. -1 for no limit [1000000]
    minblocks             The minimal number of blocks [2]
    mincallround          minimum round the detector gets called in detection loop <stairheur> [0]
    multipledecomps       Enables multiple decompositions for all enabled blocking types. Ranging from minblocks to maxblocks [TRUE]
    nconssperblock        The number of constraints per block (static blocking only) [32]
    origenabled           flag to indicate whether detector <stairheur> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <stairheur> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <stairheur> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <stairheur> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <stairheur> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <stairheur> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <stairheur> [1200]
    skip                  flag to indicate whether detector <stairheur> should be skipped if others found decompositions [FALSE]
    staticblocking        Enable blocking type 'static' [TRUE]
    usefullrecall         flag to indicate whether detector <stairheur> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_stairheur.cpp
 * [1] [Luers, M., Lübbecke, M. "Eine Heuristik zum Erkennen von
Staircase-Strukturen in Matrizen" RWTH Aachen University. Aachen, 2012.](https://www.or.rwth-aachen.de/de/abschlussarbeiten.html?file=files/research/theses/2012/MA_Luers.pdf)
 * [2] [Jayakumar, Maliyakal D., and Ranga V. Ramasesh. "A clustering heuristic to detect staircase structures in large scale linear programming models." European journal of operational research 76.1 (1994): 229-239.](https://www.sciencedirect.com/science/article/abs/pii/0377221794900191)
