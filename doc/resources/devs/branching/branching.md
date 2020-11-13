# The GCG Branching {#branching}
> The **branching on fractional variables to achieve integrality** is what guides the Branch-and-Price-and-Cut process, finally
> making sure that the solution is not just a valid solution to the LP relaxation, but also to the original MIP formulation.

During the **execution of Branch-and-Price-and-Cut**, GCG will **select nodes** to open and **select variables** to branch
on. The whole process can be lead by **tree heuristics (diving heuristics)** that heuristically determine a good path to
an integer solution. Most importantly, there are some things one has to pay attention to (GCG uses two
copies of the tree), which will be explained in the page on the branching process. Then, we also explain
what branching rules and tree heuristics are implemented in GCG.

@subpage branching-process \n
@subpage branching-rules \n
@subpage diving-heuristics \n
\n

As always with GCG, this algorithmic component is highly configurable. Thus, we want to present 
**branching parameters** that exist to tweak the process if you have knowledge of how you want to
branch on the variables of your problem.

@subpage branching-params \n
\n

If you want to **write your own plug-ins** for the branching, you might want to consider our and SCIP's how-to
pages on that, since some algorithmic components are implemented in and executed by SCIP only.

@ref own-branching-rule \n
@ref own-primal-heuristic \n
[**How to add node selectors (SCIP guide)**](https://scipopt.org/doc/html/NODESEL.php) \n
