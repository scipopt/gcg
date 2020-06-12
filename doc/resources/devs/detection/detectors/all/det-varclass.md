# Variable Class Detector {#det-varclass}

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| `v` | varclass                    | ✓ |   |   |

This detector uses the classes identified by the **variable classifiers**.

### Algorithmic Details

The Variable Class Detector creates partial decompositions by fixing some variables to the master problem and some to the linking variables, leaving the remaining variables open (for finishing). This is done as follows:
 * Iterate over all used variable classifiers
   * Iterate over all classes this classifier found and sort the respective index into the `linking` and `master` index vectors if the `DecompInfo` was set to `linking`/`master`.
   * Iterate over all possible subsets of classes (`MASTER`, `BLOCK` and `LINKING`) for this classifier.
    * Iterate over all open variables
     * If the variable class is inside the current subset of variable classes, fix it to the linking variables
     * If the decomposition info says it belongs into the linking variables, fix it to the linking variables
     * If the decomposition info says it belongs into the master variables, fix it to the master variables

This results in all possible partial decompositions where the constraints of one class are either all in the master problem or all in pricing problems.

### Theoretical Details

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
