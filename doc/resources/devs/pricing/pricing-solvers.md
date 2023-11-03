# Pricing Solvers {#pricing-solvers}

[TOC]

@todo add information on pricing solver priorities

# The GCG Pricing Solvers
During the detection, we have already identified the master and subproblem (i.e. pricing problem). In our GCG-typical
illustration, the master problem is the darker blue one which is on the top and the blocks in brighter blue are the
pricing problems. In each node of the Branch-and-Bound tree, we execute Column Generation, which means solving a
pricing problem and using this solution to 1) increase the lower bound (dual bound) as given by the pricing problem and
2) the upper bound (primal bound) as given by the master problem. Note that those two are only bounds on the linear programming relaxation of the instance, not on a possible (mixed-)integer solution.

## Solving the Pricing Problem
In order to solve these pricing problem(s), we want to use an algorithm that can solve the problem more efficiently than
usual MIP solving methods. This is why using the MIP pricing solver is in many cases an indicator that either the decomposition
of the problem was imperfect (i.e. there are "complicating constraints" in the pricing problem) or there simply is no
way known to GCG yet to efficiently solve the pricing problem at hand. However, a highly performant pricing solver is
crucial for the performance of every Branch-and-Price solver, since we have to solve a high number of pricing problems,
and that in each node of the branch-and-bound tree.

## List of Pricing Solvers
We currently have three different pricing solvers:

| Name | Description | Priority |
| -- | -- | -- |
| Knapsack ([Documentation](solver__knapsack_8c_source.html)) | Knapsack program solving using a dynamic programming approach | 200 |
| Cliquer ([Documentation](solver__knapsack_8c_source.html))  | Independent set problem solving using a heuristic (external) solver | 150 |
| MIP ([Documentation](solver__knapsack_8c_source.html))      | Simply call MIP solver again (expensive) | 0 |

All three of those can solve pricing problems heuristically or exact.

### Knapsack Solver
In very many cases, we can see variants of the so-called knapsack constraint \f$\sum_{i=1}^n x_i w_i \leq b\f$
(the knapsack's capacity is respected for the set of items we are packing into it)
and use this constraint as the only constraint of the pricing problem. By doing that, we can use
a commonly known dynamic programming approach to solve a knapsack problem optimally and very efficiently.

### Cliquer Solver
This solver can be used for independent set problems. It makes use of the external software "Cliquer".  
For installation instructions, please consult the @ref install-manually "installation guide". Note that
for the Cliquer solver to work, GCG must have been compiled with `CLIQUER=TRUE`

### Mixed Integer Program Solver
This pricing solver is simply a call of SCIP's MIP solving functionality, meaning it is no faster than
the usual solving. This consequentially leads to poor performance on many instances.

<hr>

## Adding your own Pricing Solver
If you want to **write your own pricing solver**, i.e. implement algorithmics to solve your
subproblem efficiently, please consult the "How to" for that.

@ref own-pricing-solver
