# How to obtain root node dual bounds {#dual-bounds}

[TOC]

# How to obtain root node dual bounds
In its default mode, GCG performs an automatic Dantzig-Wolfe reformulation to solve a problem. A key reason to perform a DW reformulation is a stronger dual bound obtained from the LP relaxation of the DW master problem.

During the solving process, GCG will output the latest dual bound. However, the reported dual bound is *not* necessarily identical to the bound obtained from the LP relaxation. This is due to internal optimizations in SCIP and GCG which have an impact on the dual bound.

## Obtaining LP relaxation of root node
The following settings can be used to obtain the LP relaxation of the root node:

```
set limits nodes 1
set sepa emphasis off
set pricing masterpricer advanced abortpricingint FALSE
```

These parameters result in GCG (1) solving only the root node, (2) not including cuts in the master and pricing problems, and (3) not aborting column generation *before* the restricted master problem is optimal.

After solving the problem with `optimize`, the reported final dual bound is the objective function value of the LP relaxation of the restricted master problem.
