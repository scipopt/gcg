# Constraint Class Detector {#det-consclass}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

This detector uses the classes identified by the **constraint classifiers**.

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| ?  | consclass                   | ✓ |   |   |


### Details
The Consclass Detector creates partial decompositions by fixing some constraints to the master problem and leaving the remaining constraints open (for finishing). This is done as follows:
 * Iterate over all used classifiers
 * For each classifier, iterate over all possible subsets of classes (with `ConsClassDecompInfo == BOTH`). For each subset the following partial decomposition is created:
   * The constraints beloning to classes in the current subset are fixed to the master problem
   * All constraints from classes with `ConsClassDecompInfo == ONLY_MASTER` are fixed to the master problem
   * All other constraints are left open

This resuls in all possible partial decompositions where the constraints of one class are either all in the master problem or all in pricing problems.

NOTE: In order to get sensible results with this detector a finishing detector needs to be used. Otherwise it results in only one decomposition (where all constraints are in the master problem). This is because all other partial decompositions are not finished, because the constraints which are not fixed to the master, are left open and are not assigned to blocks.

### Parameters

### Links
 * Documentation: dec_consclass.cpp
