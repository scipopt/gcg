# Diving Heuristics {#diving-heuristics}

# The GCG Diving Heuristics
GCG brings its own diving heuristics to the table. They define how one can proceed
in the tree (dive) to explore promising nodes as fast as possible by selecting respective variables.
Diving heuristics **belong to the group of primal heuristics**, which always aim at finding a 
solution quickly during the solving process, which here is specifically during the branching. 

## List of Diving Heuristics
We currently offer 5 different (GCG-specific) diving heuristics (7 more are included from SCIP):
- **Coefficient Diving** \n
LP diving heuristic that chooses fixings w.r.t. the matrix coefficients.
- **Fractional Diving** \n
LP diving heuristic that chooses fixings w.r.t. the fractionalities.
- **Guided Diving** \n
LP diving heuristic that chooses fixings in direction of incumbent solutions.
- **Pseudo Cost Diving** \n
LP diving heuristic that chooses fixings w.r.t. the pseudo cost values.
- **Vector Length Diving** \n
LP diving heuristic that rounds variables with long column vectors.

Descriptions of those and their code documentation can be
found @ref DIVINGHEURISTICS "here".

## Adding own Diving Heuristics
If you want to add your own diving heuristics (or other primal heuristics), 
i.e. define exactly **what variables to branch on**, please consider our "How to add" for that.

 ⇨ @ref own-primal-heuristic \n