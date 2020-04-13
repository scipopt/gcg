# Dense Constraints Detector {#det-densemasterconss}

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | densemasterconss            | ✓ |   |   |

This detector adds all constraints that have at least two variable coefficients that are closer than \f$1\f$ apart from another to the master constraints.

### Algorithmic Details

* Iterate over all open constraints
 * Get variables of the constraint
 * Sort variables in ascending order.
  * Iterate over all variables
   * If the difference between two successive variable coefficients is greater than or equal to \f$\Delta_{\text{max}}\f$ (per default \f$1\f$), stop counting
   * Else set \f$\Delta_{\text{max}}\f$ to the current difference and continue counting. (this leads to a stop in counting once the difference increases for the first time)
   * If a maximal ratio (\f$0.2\f$) of variables being very similar to each other is reached, stop counting

### Theoretical Details

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
