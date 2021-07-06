# Staircase Detector (Longest-Shortest Path) {#det-staircase_lsp}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| `S` | staircase_lsp               | ✓ | ✓ |   |

This detector detects staircase structures in the constraint matrix by searching for the **longest shortest path in the row graph** of the matrix.

### Algorithmic Details

### Theoretical Details

### Parameters

    enabled               flag to indicate whether detector <staircase_lsp> is enabled [FALSE]
    finishingenabled      flag to indicate whether detector <staircase_lsp> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <staircase_lsp> [1]
    maxcallround          maximum round the detector gets called in detection loop <staircase_lsp> [2147483647]
    mincallround          minimum round the detector gets called in detection loop <staircase_lsp> [0]
    origenabled           flag to indicate whether detector <staircase_lsp> is enabled for detecting in the original problem [FALSE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <staircase_lsp> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <staircase_lsp> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <staircase_lsp> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <staircase_lsp> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <staircase_lsp> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <staircase_lsp> [200]
    skip                  flag to indicate whether detector <staircase_lsp> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <staircase_lsp> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_staircase_lsp.cpp
