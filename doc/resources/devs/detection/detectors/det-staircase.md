# Staircase Detector {#det-staircase}

This detector detects staircase structures in the constraint matrix by searching for the longest shortest path in the row graph of the matrix.

### Details
1. Given a constraint matrix, construct an undirected graph as follows: a vertex is added to the graph for every row of the matrix, and two vertices are connected by an edge if and only if there exists a column in the matrix that has a nonzero coefficient in both rows that are represented by the vertices.
2. Find the connected components of the graph.
3. For each connected component, compute the diameter and find an endpoint of a longest shortest path.
4. For each component, start BFS at one of the endpoints and assign all vertices in each BFS layer to a block of the decomposition.

### Parameters
```
  enabled               flag to indicate whether detector <staircase> is enabled [TRUE]
  priority              priority of detector <staircase> [200]
  skip                  flag to indicate whether detector <staircase> should be skipped if others found decompositions [FALSE]
```
### Future work
Each BFS layer does not necessarily form a connected connected component of the graph and therefore it could be split up into different blocks if it is not connected.

### Links
 * Documentation:http://www.or.rwth-aachen.de/gcg/doc/dec__staircase_8c.html
