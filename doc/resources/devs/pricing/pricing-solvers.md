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
| Knapsack ([Documentation](solver__knapsack_8c.html)) | Knapsack program solving using a dynamic programming approach | 200 |
| Cliquer ([Documentation](solver__knapsack_8c.html))  | Independent set problem solving using a heuristic (external) solver | 150 |
| GCG ([Documentation](solver__gcg_8c.html))      | Call GCG solver | 110 |
| HiGHS ([Documentation](solver__highs_8c.html))      | Call (external) HiGHS solver | 100 |
| CPLEX ([Documentation](solver__cplex_8c.html))      | Call (external) CPLEX solver | 100 |
| MIP ([Documentation](solver__knapsack_8c.html))      | Simply call MIP solver again (expensive) | 0 |

All three of those can solve pricing problems heuristically or exact.

### Knapsack Solver
In very many cases, we can see variants of the so-called knapsack constraint \f$\sum_{i=1}^n x_i w_i \leq b\f$
(the knapsack's capacity is respected for the set of items we are packing into it)
and use this constraint as the only constraint of the pricing problem. By doing that, we can use
a commonly known dynamic programming approach to solve a knapsack problem optimally and very efficiently.

### Cliquer Solver
This solver can be used for independent set problems. It makes use of the external software "Cliquer".  
For installation instructions, please consult the @ref install-manually "installation guide". Note that
for the Cliquer solver to work, GCG must have been compiled with the enabled `CLIQUER` flag.

### Mixed Integer Program Solver
This pricing solver is simply a call of SCIP's MIP solving functionality, meaning it is no faster than
the usual solving. This consequentially leads to poor performance on many instances.

### GCG Solver
This pricing solver uses GCG recursively to solve the pricing problems.
The solver is disabled by default. It can be enabled by setting the maximum (recursion) depth parameter `pricingsolver/gcg/maxdepth` to a value greater than 0.
I.e., setting this parameter to 1 means that the GCG solver is used to solve the root pricing problems, but is disabled for the nested pricing problems.
The solver either uses a nested decomposition provided by the user, or tries to detect a structure if no nested structure information is available.

### HiGHS Solver
This pricing solver uses the HiGHS shared library to solve the pricing problems.
See the [documentation](https://ergo-code.github.io/HiGHS/dev/) and the article (DOI: [10.1007/s12532-017-0130-5](https://link.springer.com/article/10.1007/s12532-017-0130-5)):
```
Parallelizing the dual revised simplex method, Q. Huangfu and J. A. J. Hall, Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: 10.1007/s12532-017-0130-5
```
Note that for the HiGHS solver to work, GCG must have been compiled with the enabled `HIGHS` flag.

### CPLEX Solver
This pricing solver uses the CPLEX shared library to solve the pricing problems.
See the CPLEX [website](https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer).
Note that for the CPLEX solver to work, GCG must have been compiled with the enabled `CPLEX` flag.

<hr>

## Adding your own Pricing Solver
If you want to **write your own pricing solver**, i.e. implement algorithmics to solve your
subproblem efficiently, please consult the "How to" for that.

@ref own-pricing-solver
