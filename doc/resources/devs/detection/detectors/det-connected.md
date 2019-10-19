# Connected Detector {#det-connected}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

The Connected Detector finds variables shared by constraints and places those constraints into the same block.

### Details

To detect blocks the detector uses the following steps:

1. The constraints of the problem are partitioned into a master and a rest problem. The master problem contains all set covering and set partitioning constraints.
2. All non-master constraints are checked for connectivity as follows:
  * Initially, there are no blocks.
  * The variables of a constraint are checked for assignment to a block. If so, the constraint is assigned to the same block. Else, the constraint and its variables are assigned to a new block.
  * The step beforehand is applied to all constraints.

### Parameters
```
  enabled               flag to indicate whether detector <connected> is enabled [FALSE]
  priority              priority of detector <connected> [0]
  setppcinmaster        Controls whether SETPPC constraints chould be ignored while detecting and be directly placed in the master [TRUE]
  skip                  flag to indicate whether detector <connected> should be skipped if others found decompositions [FALSE]
```
### Future work

There is current work on
 * finding bordered structures with the Connected Detector as well as any other connector that defines a master problem. This will be done by finding other structures within the non-master problem by additionally to finding connected components applying other detectors to a subscip.
Additional parameters will be
  ```
  maxnvarsstaircase    variable limit for finding bordered staircase structures
  maxnconsstaircase    constraint limit for finding bordered staircase structures
  ```
Currently this is being implemented in the method DECcreateDecompFromMasterconss in decomp.c. As this method originally only finds connected components and is used by other detectors, too, there is discussion on moving this feature to a new function instead.
 * allowing decompositions with only one pricing problem by just removing generalized covering and partitioning constraints.

### Links
Documentation: http://www.or.rwth-aachen.de/gcg/doc/dec__connected_8c.html
