# The GCG Pricing {#pricing}

> The **pricing of variables (columns)** in a given model is one of the main features that
> distinguishes GCG from SCIP. Just like the detection and decomposition,
> it accelerates the solving speed of many problems.

Since GCG uses a Branch-and-Price approach to solve Mixed-Integer Programs, we have an own
**Pricing algorithm** at hand. In simple terms, GCG creates a **Restricted Master Problem**
with only very few variables and a Sub-Problem (**Pricing Problem**).
The latter is than undergoing a pricing algorithm in each node of the Branch-and-Bound tree.
There, we implicitly (heuristically) go through "all" possible variables (*columns* in the coefficient matrix) and
*price* them, i.e. we calculate their reduced cost (*reduced cost pricing*).
If a column yields a negative reduced cost, it will improve our currently best
solution and thus we take it into our Restricted Master Problem.

If you want to know more about this process, please check out the respective explanation pages
for the components used during the pricing.

@subpage pricing-process \n
@subpage pricing-solvers \n
@subpage recursive-dw \n


There are many different parameters that can influence the pricing of your problem.
We collected some information about those and also some hints about possibly improving
**parameter configurations** for specific problems.
It is assumed that you are rather familiar with the technical details of Branch-and-Price here.

@subpage pricing-params \n

<hr>

If you want to **write your own** pricing solver, please consult the "How to"
for that.

@ref own-pricing-solver \n
