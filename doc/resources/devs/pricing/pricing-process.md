# Pricing Process {#pricing-process}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

In every iteration of the branch-and-bound tree, instead of an LP solver, the GCG relaxator
(see `relax_gcg.c`) and perform a Dantzig-Wolfe reformulation and then do Column Generation.

# Dantzig-Wolfe Reformulation {#dantzig}

# The Restricted Master Problem (RMP) {#rmp}
`SCIPsolve` is called again for the RMP.

# The Pricing Problem (PP) {#pp}

# Column Generation {#cg}

## Pricing Solver
The pricing problem is then solved as follows:
- Loop through pricing subproblems in rounds
- Combine next pricing subproblem in queue with pricing solver
- Queue can be processed in parallel

If the chosen pricing solver is MIP, `SCIPsolve` is called again for the PP.
