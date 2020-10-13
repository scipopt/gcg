# Presolving Process {#presolving-process}

# Presolving in GCG
## Overview
Apart from the usual presolving process that is very similar to the one that SCIP is using (see "Differences
to SCIP" for differences), it is worth noting that **variable types and constraint types are set during
presolve.**

## Problems with Presolving
When using usual presolving, there arises the problem that **blocks can be damaged**, leading to a problem
that was decomposable, but is not anymore due to presolving. This will result in much higher solving durations
if a decomposition that was made on the presolved problem is finally chosen by GCG.

## Differences to SCIP
### Dual Presolving
There are different types of presolving and one of them is dual presolving. Using the reduced cost of variables,
we can preprocess models in such a way that they are strengthened. In detail, if a variable yields higher reduced cost
than the difference to the lower bound currently is, it cannot be in the basis of a feasible solution. Thus, we do not take it.
For different reasons, in GCG, dual presolving is off.