# Type Classifier (SCIP Constraint Type) {#clscons_scipconstypes}

This classifier adds constraints to classes according to their **type given by SCIP**.
There are different constraint types supported:

 * Linear Constraints (`linear`)
 * Knapsack Constraints (`knapsack`)
 * Variable Bounds (`varbound`)
 * Set Packing Constraints (`setpacking`)
 * Set Covering Constraints (`setcovering`)
 * Set Partitioning Constraints (`setpartitioning`)
 * Logical OR Constraints (`logicor`)
 * Special Ordered Set (at most one variable non-zero) (`sos1`)
 * Special Ordered Set (at most two variables non-zero (must be adjacent)) (`sos2`)
 * Unknown Type Constraints (`unknown`)

### Classification
#### Dantzig-Wolfe and Benders Decomposition
In both modes, the handling is as follows:
  * A new class is generated for each of the aforementioned types.
  * By default, the class with name `newConstype` is generated

### Example
