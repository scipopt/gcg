# Why use GCG? {#why-gcg}
> On this page, we give a brief overview of **what GCG is capable of** (and what it is not). 

# Why should you use GCG?
If you have solved mixed integer linear programs already, you will know that it can take _forever_ to solve
them, even if they are well thought out. Whether you did think your model through or you didn't, it will always
be worth a try to solve them with GCG, because we deploy a variety of tools to possibly increase solving speed
in orders of magnitudes.

@todo add performance profile between SCIP and GCG

## Introduction to GCG
GCG will decompose your problem into problems that can be solved far more easily (decomposition). Then, we 
apply a Dantzig-Wolfe reformulation to strengthen your formulation. This reformulation will convexify the 
"easy" constraints, meaning that the polyhedron (i.e. the space where your solution can lie in) is 
shrinked down as much as possible. Finally, we solve your problem using Branch-and-Price. In each node
of the Branch-and-Bound tree we solve the pricing problems, getting variables that are promising.

### Features
**As an extendible and easy-to-use toolkit, GCG...**
- solves your problem using complex solving algorithmics **without requiring any input or knowledge from your side**
- has an **explore** function to investigate identified decompositions
- can read your **custom decompositions** using the DEC file standard
- is a framework for **implementing** column generation algorithms
For further use cases of GCG, please consult the @ref users.

**As a solver, GCG...**
- **detects hidden or apparent structures** in the constraint matrix in order to apply DWR. Among others, it detects
  -  single-bordered structures
  -  arrowhead structures (using [hmetis](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview))
  -  staircase structures
  -  set partitioning master structures
  -  clustered and graph connectivity based structures 
- uses a **Dantzig-Wolfe reformulation** (DWR) to solve arbitrary MIPs.
- automatically aggregates subproblems if possible
- offers an automated **Benders'** decomposition algorithm
- has **branching rules** for automatic branching on any problem, e.g. branching on original, Ryan-Foster branching or generic branching
- applies a wide variety of **cuts** to the problem, e.g. combinatorially or from basis)
- rates columns in its (parallel) exact or heuristic **pricing**
- has a large number of **primal heuristics** both in the original and the reformulated space
- applies **dual variable stabilization** by dual value smoothing

For more details on what and how GCG does things to solve your problem quickly, please consult the @ref devs.

## When is it sensible to use GCG?
In order for a decomposition to make sense, **your problem has to exhibit some kind of structure**. In many
cases, it will have a certain structure (maybe even one that is not known to you), but if it does not have 
one at all, it can be detrimental to performance to try to use a decomposition anyways. 
However, if it has a structure there is a **good chance that GCG will be able to solve it much faster** than 
other state-of-the-art open-source solvers, 
such as SCIP, and often at least as fast as commercial solvers such as Gurobi or CPLEX.  

<table>
  <tr>
    <th>GCG</th>
    <th>SCIP</th>
  </tr>
  <tr>
    <td>
    <div class="fragment">
      <div class="line">SCIP Status        : problem is solved [optimal solution found]</div>
      <div class="line">Solving Time (sec) : 28.48</div>
      <div class="line">Solving Nodes      : 1</div>
      <div class="line">Primal Bound       : -4.10000000000000e+01 (5 solutions)</div>
      <div class="line">Dual Bound         : -4.10000000000000e+01</div>
      <div class="line">Gap                : 0.00 %</div>
    </div>
    </td>
    <td>
    <div class="fragment">
      <div class="line">SCIP Status        : problem is solved [optimal solution found]</div>
      <div class="line">Solving Time (sec) : 337.72</div>
      <div class="line">Solving Nodes      : 1088063 (total of 1088983 nodes in 2 runs)</div>
      <div class="line">Primal Bound       : -4.10000000000000e+01 (389 solutions)</div>
      <div class="line">Dual Bound         : -4.10000000000000e+01</div>
      <div class="line">Gap                : 0.00 %</div>
    </div>
    </td>
  </tr>
</table>
