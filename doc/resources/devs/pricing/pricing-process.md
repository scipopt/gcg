# Pricing Process {#pricing-process}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

Our detection will yield a decomposition of the problem into master constraints, linking variables and 
blocks (see @ref detection). Using this structure, we can apply a Dantzig-Wolfe Reformulation using 
the theorem of Minkowski-Weyl to reformulate the constraints that we "can handle well". 
This will create a stronger model already.

## Branch-and-Price
In every iteration of the branch-and-bound tree, instead of an LP solver, the GCG relaxator
(see `relax_gcg.c`) is called to perform a Dantzig-Wolfe reformulation on the current problem,
with the current bounds of the branch-and-bound tree that are on the LP.

### Enforcing Bounds during Branching
#### Theoretical Point of View
In GCG, we branch on the original variables. Now, if we want to branch down (or "left"), we want 
to impose \f$x_i\leq\lfloor x_i^*\rfloor\f$ in the master problem. But, as the name "original variables"
suggests: we do not have \f$x_i\f$s anymore, instead we have convex and conic combinations of 
points \f$q\in Q\f$ and rays \f$r\in R\f$, thus we have to enforce 
\f{align}{
  \sum_{q\in Q} x_{qi}\lambda_q+\sum_{r\in R} x_{ri}\lambda_r \leq \lfloor x_i^*\rfloor
\f}

#### Practical Point of View
It would be **infeasible to store an LP for every node** of the Branch-and-Bound tree, since the number 
of nodes in the tree is growing exponentially and the LP can already be quite big in itself. This
is why we **only save which bounds we enforce in each node**. Then, if we enter a node that is 
succeeding the current node, we just have to add a new constraint, reformulate and solve each
pricing problem (see next section). But if, possibly due to a node selection technique 
that we are using, we want to go to a "whole different part" in the Branch-and-Bound tree, 
we have to **iteratively undo each bound we imposed** while going there. 


### Pricing Loop
Each block represents one pricing problem that can be solved independently from all the others and 
at first even independent of the master constraints. Since we want to generate columns (variables)
for the Restricted Master Problem that is yet lacking an optimal integer solution, **we solve the 
pricing problems and add the best variables** ("most negative" reduced costs) to the RMP. For more
information on how we do that, please read the section on "Storing and Using Columns".
The SCIP-Instance corresponding to a pricing problem is created inside `relax_gcg.c` and freed in exitSol().

The **pricing problem is then solved** as follows:
- Loop through pricing subproblems in rounds
- Combine next pricing subproblem in queue with pricing solver
- Queue can be processed in parallel

If the chosen pricing solver is MIP, `SCIPsolve` is called again for the pricing problem.

#### Storing and Using Columns
When performing Column Generation, it often happens that it is **very expensive to generate variables**.
Thus, GCG will take _all_ columns that are generated that have negative reduced cost, because we might
be able to use them at some point. That point might be at different times. For that purpose, we have 
two different structures storing columns (variables) that are generated during the pricing.
- `Pricestore`: contained columns survive only in the _current round_.
- `Column Pool`: contained columns survive into _future rounds_.