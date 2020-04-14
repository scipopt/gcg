# Constraint Neighborhood Detector {#det-neighborhoodmaster}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | neighborhoodmaster          | ✓ |   |   |

This detector calculates cons-cons adjacency (if not already done) and sorts constraints according to the size of their neighborhood. It then searches for two consecutive constraints with largest size difference (according neighborhood size) in the sorted constraints. All constraints having a larger neighborhood than the second one are then assigned to the master.

### Algorithmic Details
* Iterate over all open constraints (\f$N\f$ many)
 * For each constraint, get the neighborhood size (\f$n_{\text{neighborhood}}=\f$getNConssForCons()) and make a tuple of the two: \f$(n,c)\f$, with \f$c\f$ the constraint and add it to the list `l`
* Sort the tuples in `l` in ascending order by the first element (the neighborhood, higher means more dense).
* Set the `last index` to be \f$r_{max}*N\f$, with \f$r_{max}\f$ being the parameter "maximal ratio" (by default \f$0.2\f$)
* Iterate over all constraints
  * If the `last index` (i.e. by default the first \f$20\%\f$ of the constraints, starting with the greatest neighborhood size) is reached, stop.
  * If the difference in neighborhood size between the current item `l[i]` and the next item `l[i+1]` is smaller than the default `maximum difference index`=\f$1\f$, set `i` to be the `maximum difference index`, i.e. if there is a jump between two constraints in their neighborhood size, stop.
* Set all constraints up until the `maximum difference index` to be master constraints, i.e. at max the first \f$20\%\f$ of the constraints with the highest and similar number of nonzeros.

### Theoretical Details

### Parameters

    enabled               flag to indicate whether detector <neighborhoodmaster> is enabled [TRUE]
    finishingenabled      flag to indicate whether detector <neighborhoodmaster> is enabled for finishing of incomplete decompositions [FALSE]
    freqcallround         frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <neighborhoodmaster> [1]
    maxcallround          maximum round the detector gets called in detection loop <neighborhoodmaster> [0]
    maxratio              the maximal ratio of open constraints that are assigned to the master problem [0.2]
    mincallround          minimum round the detector gets called in detection loop <neighborhoodmaster> [0]
    origenabled           flag to indicate whether detector <neighborhoodmaster> is enabled for detecting in the original problem [TRUE]
    origfreqcallround     frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <neighborhoodmaster> [1]
    origmaxcallround      maximum round the detector gets called in detection loop <neighborhoodmaster> [2147483647]
    origmincallround      minimum round the detector gets called in detection loop <neighborhoodmaster> [0]
    overruleemphasis      flag to indicate whether emphasis settings for detector <neighborhoodmaster> should be overruled by normal settings [FALSE]
    postprocessingenabled flag to indicate whether detector <neighborhoodmaster> is enabled for postprocessing of finished decompositions [FALSE]
    priority              priority of detector <neighborhoodmaster> [0]
    skip                  flag to indicate whether detector <neighborhoodmaster> should be skipped if others found decompositions [FALSE]
    usefullrecall         flag to indicate whether detector <neighborhoodmaster> should be called on descendants of the current partialdec [FALSE]


### Links
 * Documentation: dec_neighborhoodmaster.cpp
