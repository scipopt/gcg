# How to add pricing solvers {#own-pricing-solver}
> **This page is currently being refactored. Some things might still be outdated.**

[TOC]

Pricing problem solvers are **called by the pricer to solve a single pricing
problem** either heuristically or to optimality and return a set of solutions.
\n
A complete list of all pricing problem solvers contained in this release can be found [here](#pricing-solvers).

With the following steps, we explain how you can **add your own pricing problem solver plug-in**:
1. **Preparations**
  1. Choose a name `mysolver` for your solver.
  2. Copy the template files `src/solver_xyz.c` and `src/solver_xyz.h`
   while renaming `xyz` to `mysolver`.
  3. Open the new files with a text editor and replace all occurrences of `Xyz` by `Mysolver` and `xyz` by `mysolver`.
2. **Creating your Pricing Solver**
  1. Adjust the properties of the pricing solver (see @ref SOLVER_PROPERTIES).
  2. [optional] Define the pricing solver data (see @ref SOLVER_DATA).
  3. Implement the interface methods (see @ref SOLVER_INTERFACE).
  4. Implement the fundamental callback methods (see @ref SOLVER_FUNDAMENTALCALLBACKS).
  5. [optional] Implement the additional callback methods (see @ref SOLVER_ADDITIONALCALLBACKS).
3. **Make GCG use it**
  1. Add it to masterplugins.c by adding
    1. the line <tt>\#include solver_mysolver.h</tt>.
    2. the line `SCIP_CALL( GCGincludeSolverMysolver(scip) );` in the method SCIPincludeGcgPlugins().
  2. Add it to your build system:
    1. _Using Makefile:_ Add your solver `.o` (`solver_mysolver.o`) to the list below `LIBOBJ =` in the file `Makefile` in the root folder.
    2. _Using CMake:_ In `src/CMakeLists.txt`, add your `solver_mysolver.c` below `set(gcgsources` and your
   `solver_mysolver.h` below `set(gcgheaders`.


# Properties of a pricing problem solver {#SOLVER_PROPERTIES}

At the top of the new file "solver_mysolver.c", you can find the solver properties.
These are given as compiler defines.
The properties you have to set have the following meaning:

\par SOLVER_NAME: the name of the solver.
This name is used in the interactive shell to address the solver.
Names have to be unique: no two solvers may have the same name.

\par SOLVER_DESC: the description of the solver.
This string is printed as description of the solver in the interactive shell.

\par SOLVER_PRIORITY: the priority of the solver.
Whenever a pricing problem has to be solved, the solvers are called in a predefined order, until one of the solvers solved the
pricing problem. This order is given by the priorities of the solvers; they are called in the order of decreasing priority.
<br>
The priority of the solver should be set according to the range of problems it can solve. The more specialized a solver is,
the higher its priority should be in order to be called before more general solvers that also cover this type of problem.
The default pricing solver that solves the pricing problem as a MIP using SCIP has priority 0, so all other pricing solvers
should have positive priority.
An easy way to list the priorities and descriptions of all solvers is "display solvers" in the interactive shell of GCG.

# Solver Data {#SOLVER_DATA}

Below the header "Data structures" you can find the struct "struct GCG_SolverData".
In this data structure, you can store the data of your solver. For example, you should store the adjustable parameters
of the solver in this data structure and also arrays to store the solutions that you return after solving a pricing problem
(see solver_mip.c).
Defining solver data is optional. You can leave this struct empty.


# Interface Methods {#SOLVER_INTERFACE}

At the bottom of "solver_mysolver.c", you can find the interface method GCGincludeSolverMysolver(),
which also appears in "solver_mysolver.h".
\n
This method has to be adjusted only slightly.
It is responsible for notifying GCG (and especially the pricer in the master SCIP instance) of the presence of the solver by
calling the method GCGpricerIncludeSolver().
SCIPincludeSolverMysolver() is called by the user to include the solver,
i.e., to use the solver in the application (see 3.1.1. at the top of the page).

If you are using solver data, you have to allocate the memory for the data at this point.
You can do this by calling
```C
SCIP_CALL( SCIPallocMemory(scip, &solverdata) );
```
You can also initialize the fields in struct SCIP_SolverData afterwards. For freeing the
solver data, see [SOLVEREXIT](#SOLVER_EXIT). Alternatively, you can also initialize and free the fields in the solver data
in the SOLVERINITSOL and SOLVEREXITSOL callback methods, respectively. For an example, see solver_mip.c.

You may also add user parameters for your solver, see the SCIP documentation for how to add user parameters and
the method SCIPincludeSolverMip() in solver_mip.c for an example.


# Fundamental Callback Methods of a Solver {#SOLVER_FUNDAMENTALCALLBACKS}

The fundamental callback methods of the plug-ins are the ones that have to be implemented in order to obtain
an operational algorithm. Pricing problem solvers do not have fundamental callback methods, but they should
implement at least of the [SOLVERSOLVE](#SOLVER_SOLVE) and [SOVLERSOLVEHEUR](#SOLVER_SOLVEHEUR) methods.

Additional documentation to the callback methods, in particular to their input parameters,
can be found in type_solver.h.

# Additional Callback Methods of a Solver {#SOLVER_ADDITIONALCALLBACKS}

## SOLVERSOLVE {#SOLVER_SOLVE}

The SOLVERSOLVE callback is called when the variable pricer in the master SCIP instance wants to solve a pricing problem
to optimality.
It is given a SCIP instance representing the pricing problem that should be solved and should check whether the pricing
problem has a structure that can be handled by the solver. If so, it should solve the pricing problem to optimality
and return the optimal objective value of the problem (to be stored in the lowerbound pointer) along with a set of
primal feasible solutions (including the optimal one).

The solutions should be returned by setting some given pointers:
solvars should store the number of solutions returned, nsolvars should point to an array storing for each solution
the number of variables with non-zero value, and solisray should point to an array storing for each solution whether it
represents an extreme point of the pricing problem or an extreme ray. Furthermore, solvars and solvals should point to
arrays storing for each solution an array of variables and corresponding solution values, respectively (leaving out
variables with value 0). Therefore, a pricing problem solver should manage these arrays in its own data structure,
fill them after solving a pricing problem and set the return pointers to these arrays. Have a look at solver_mip.c for
an example of how this can be done.

Last, the callback should adjust the given result pointer to SCIP_STATUS_OPTIMAL if the problem was solved to optimality,
to SCIP_STATUS_UNBOUNDED if the problem was solved and is unbounded, or SCIP_STATUS_UNKNOWN if the solver was not applicable
to the pricing problem or if the solving was stopped.

The given SCIP instance representing the pricing problem can be seen as a container to store all information about the pricing
problem. During the solving process, especially the objective function coefficients change according to the current dual
solution and branching might change bounds in the pricing problem or add constraints.
If the structure of the problem is independent of the changes that can occur with the selected branching rule, the structure
detection for the pricing problems can also be done in the [SOLVERINITSOL](#SOLVER_INITSOL) callback and stored internally to avoid doing it every
time a pricing problem is solved.

## SOLVERSOLVEHEUR {#SOLVER_SOLVEHEUR}

The SOLVERSOLVEHEUR callback is called during heuristic pricing when the variable pricer in the master SCIP instance
wants to solve a pricing problem heuristically.
It has the same input and return parameters as the [SOLVERSOLVE](#SOLVER_SOLVE) callback. It does not need to solve the pricing problem to
optimality, but should try to construct good feasible solutions using fast methods. Nevertheless, it can return a lower bound
for the optimal solution value of the problem, if it computes one.

## SOLVERFREE {#SOLVER_FREE}

If you are using solver data, you have to implement this method in order to free the solver data.
This can be done by the following procedure:
```C
static
GCG_DECL_SOLVERFREE(solverFreeMysolver)
{
  GCG_SOLVERDATA* solverdata;

  assert(scip != NULL);
  assert(solver != NULL);

  solverdata = GCGpricerGetSolverdata(scip, solver);
  assert(solverdata != NULL);

  SCIPfreeMemory(scip, &solverdata);

  GCGpricerSetSolverdata(scip, solver, NULL);

  return SCIP_OKAY;
}
```
If you have allocated memory for fields in your solver data, remember to free this memory
before freeing the pricing solver data itself.

## SOLVERINIT {#SOLVER_INIT}

The SOLVERINIT callback is executed after the problem is transformed.
The pricing problem solver may, e.g., use this call to replace the original constraints stored in its solver data by transformed
constraints, or to initialize other elements of his solver data.

## SOLVEREXIT {#SOLVER_EXIT}

The SOLVEREXIT callback is executed before the transformed problem is freed.
In this method the pricing problem solver should free all resources that have been allocated for the solving process in [SOLVERINIT](#SOLVER_INIT).

## SOLVERINITSOL {#SOLVER_INITSOL}

The SOLVERINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to begin.
The pricing problem solver may use this call to initialize its branch-and-bound specific data.

## SOLVEREXITSOL {#SOLVER_EXITSOL}

The SOLVEREXITSOL callback is executed before the branch-and-bound process is freed.
The pricing problem solver should use this call to clean up its branch-and-bound data, which was allocated in SOLVERINITSOL.
