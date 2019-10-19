# Gomory-Hu Detector {#det-gomory}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

This detector constructs a Gomory-Hu tree of a hypergraph where the vertices correspond to rows of the LP and the hyperedges correspond to columns of the LP. Each edge of the tree has a weight which is the size of a minimum cut between the two nodes incident to the edge. Then, a decomposition is constructed by merging nodes of the tree and assigning each node (which represents several vertices) to a block.

### Details
1. Construct a hypergraph H = (V, E) as follows: add a vertex for every row of the constraint matrix, and add an edge for every column of the constraint matrix. A vertex is incident to an edge if and only if the corresponding matrix coefficient is nonzero.
2. Using Gusfield's algorithm [1], construct a Gomory-Hu tree of the hypergraph by computing a sequence of min-cuts between pairs of vertices.
3. Construct a decomposition as follows:
  1. Initially, assign each node of the tree to its own block.
  2. As long as there are "too many" linking edges (in the hypergraph) that are incident to nodes from different blocks (in the tree), find a tree-edge of maximum weight, contract it, and merge the blocks to which the endpoints of the edge belong.

Currently, there are "too many" linking variables if there are more linking variables than a constant times the total number of variables.

### Parameters
```
  deckpath              whether to construct a decomposition from BFS layers of a traversal starting at an endpoint of a shortest k-path [FALSE]
  declongestpath        whether to construct a decomposition from BFS layers of a traversal starting at a node of maximum eccentricity [FALSE]
  decmerge              whether to construct a decomposition by merging tightly connected components [TRUE]
  enabled               flag to indicate whether detector <gomoryhu> is enabled [TRUE]
  kpath                 length of a shortest path of the gomory-hu tree [10]
  percentage            Fraction of variables that is allowed as linking variables before merging is stopped [0.025]
  priority              priority of detector <gomoryhu> [0]
  skip                  flag to indicate whether detector <gomoryhu> should be skipped if others found decompositions [FALSE]
```
### Finding Min-Cuts in Hypergraphs
There exists an algorithm by Klimmek and Wagner [2] that directly works on hypergraphs, and there exists a reduction (see Lawler [3]) from hypergraphs to simple graphs to which the usual max-flow algorithms like Push-Relabel can be applied.

In the current implementation, a Push-Relabel algorithm with the Gap Relabeling heuristic (see Cherkassky and Goldberg [4]) is implemented and the reduction to a simple graph is only performed implicitly.

### Future work
Different ways to construct decompositions from Gomory-Hu trees need to be developed as contracting tree-edges of maximum weight does not always lead to good decompositions. For example, this rule does not take into account the sizes of the blocks that are merged together. This often results in decompositions with few big blocks and some blocks that only contain single constraints.

### Links
* [[1] Dan Gusfield, Very Simple Methods for All Pairs Network Flow Analysis](http://epubs.siam.org/doi/abs/10.1137/0219009)
* [[2] Regina Klimmek, Frank Wagner, A Simple Hypergraph Min Cut Algorithm](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.31.4535)
* [[3] E.L. Lawler, Cutsets and partitions of hypergraphs](http://onlinelibrary.wiley.com/doi/10.1002/net.3230030306/abstract)
* [[4] Boris V. Cherkassky, Andrew V. Goldberg, On implementing push-relabel method for the maximum flow problem](http://link.springer.com/chapter/10.1007%2F3-540-59408-6_49)
