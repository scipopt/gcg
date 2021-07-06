# Isomorph Detector {#det-isomorph}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

### Overview

| ID |          Full Name          | Propagate | Finish | Postprocess |
|----|-----------------------------|:---------:|:------:|:-----------:|
| `I` | isomorph                    | ✓ |   |   |

This detector finds subproblems that can be aggregated thus reducing the symmetry of the problem using color preserving automorphisms and bliss.

Note: This detector requires BLISS, i.e. a `make` with `BLISS=true`.

### Details

To find subproblems the detector proceeds as follows:

1. Build a tripartite graph:
 * Add nodes for all constraints with the color of its equivalence class
 * Add nodes for all variables with the color of its equivalence class
 * Add nodes for nonzero entries the color of its equivalence class
 * connect nonzero entries with their corresponding constraint and variable nodes
2. Search for a color-preserving automorphism of the graph using BLISS
3. If a color-preserving automorphism is found, it gives rise to permutations of the constraints and variables. If the permutation of the constraints is non-trivial, subsets of nodes with a bijection of one to another can be identified.
4. Put the constraints of the different subsets in different pricing problems. These problems can be aggregated.


### Parameters
```

```

### Future work

If the graph is not coherent, BLISS returns more than one generator. To get all possible automorphisms, all generators have to be "connected" in all possible ways.


### Links
 * Documentation: dec_isomorph.cpp
