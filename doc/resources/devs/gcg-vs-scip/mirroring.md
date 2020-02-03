# Mirroring of Branching Decisions to SCIP {#mirroring}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

GCG mirrors its branching decisions over to SCIP such that SCIP can execute them.

## Technical Details
The code for the communication to SCIP during branching on original variables is inside the
`cons_masterbranch.c` and `cons_origbranch.c` source files. The process is as follows:
We (since we can make better branching decisions in most cases) branch ourselves
(`cons_masterbranch`) and then mirror those decisions to SCIP (`cons_origbranch`)
where they are reconstructed.

In the case that an aggregation took place, we do not do branching on original variables.

@todo add stuff from presentations?
