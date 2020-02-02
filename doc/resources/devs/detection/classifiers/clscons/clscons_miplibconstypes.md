# Type Classifier (MIPLIB Constraint Type) {#clscons_miplibconstypes}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

This classifier adds constraints to classes according to their **type as defined by the MIPLIB and given by SCIP**.
There are different constraint types supported:

 * `SCIP_CONSTYPE_EMPTY`
 * `SCIP_CONSTYPE_FREE`
 * `SCIP_CONSTYPE_SINGLETON`
 * `SCIP_CONSTYPE_AGGREGATION`
 * `SCIP_CONSTYPE_{VARBOUND}`
 * `SCIP_CONSTYPE_SETPARTITION`
 * `SCIP_CONSTYPE_SETPACKING`
 * `SCIP_CONSTYPE_SETCOVERING`
 * `SCIP_CONSTYPE_CARDINALITY`
 * `SCIP_CONSTYPE_INVKNAPSACK`
 * `SCIP_CONSTYPE_EQKNAPSACK`
 * `SCIP_CONSTYPE_BINPACKING`
 * `SCIP_CONSTYPE_KNAPSACK`
 * `SCIP_CONSTYPE_INTKNAPSACK`
 * `SCIP_CONSTYPE_MIXEDBINARY`

Descriptions of the constraint types are given [here](https://miplib.zib.de/statistics.html).

### Classification
#### Dantzig-Wolfe and Benders Decomposition
In both modes, the handling is as follows:
  * A new class is generated for each of the aforementioned types.
  * By default, the class with name `newConstype` is generated

### Example
