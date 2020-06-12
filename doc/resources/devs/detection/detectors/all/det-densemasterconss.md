# Constraint Density Detector {#det-densemasterconss}

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| `D` | densemasterconss            | ✓ |   |   |

This detector adds all constraints that have at least two variable coefficients that are closer than \f$1\f$ apart from another to the master constraints.

### Algorithmic Details

* Iterate over all open constraints (\f$N\f$ many)
 * For each constraint, get the number of nonzeros (\f$n_{\text{nonzeros}}=\f$`getNVarsForCons()`) and make a tuple of the two: \f$(n,c)\f$, with \f$c\f$ the constraint and add it to the list `l`
* Sort the tuples in `l` in ascending order by the first element (the neighborhood, higher means more dense).
* Set the `last index` to be \f$r_{max}*N\f$, with \f$r_{max}\f$ being the parameter "maximal ratio" (by default \f$0.2\f$)
* Iterate over all constraints
  * If the `last index` (i.e. by default the first \f$20\%\f$ of the constraints, starting with the greatest number of nonzeros) is reached, stop.
  * If the difference in nonzeros between the current item `l[i]` and the next item `l[i+1]` is smaller than the default `maximum difference index`=\f$1\f$, set `i` to be the `maximum difference index`, i.e. if there is a jump between two constraints in their number of nonzero entries, stop.
* Set all constraints up until the `maximum difference index` to be master constraints, i.e. at max the first \f$20\%\f$ of the constraints with the highest and similar neighborhood size.

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
