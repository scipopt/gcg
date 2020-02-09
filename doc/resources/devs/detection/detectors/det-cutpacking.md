# Cutpacking Detector (not in refactoring) {#det-cutpacking}
> This detector will be removed in a future version.\n
> It was replaced by the @ref det-staircase-lsp and @ref det-staircaseheur.

This detector tries to find a staircase structure by recursively partitioning the row graph of the matrix.

### Details

The detector proceeds as follows:

 1. Define the row graph as follows:
    * Each vertex corresponds to a row of the constraint matrix
    * There is an edge between to vertices iff the corresponding constraints share at least one variable with nonzero coefficient
    * The weight of an edge is the number of these variables
 2. Compute a minimum cut by either
    * using the graph partitioner `hmetis` [1] as an external tool
    * or applying the Stoer-Wagner algorithm [2]
 3. Partition the graph according to the cut; each side of the partition will define a block in the matrix, with the linking constraints being those whose vertices are adjacent to a vertex from the other side if the partition
 4. For each side of the partition, if the corresponding block contains at least one non-linking constraint:
    * Define a new graph by merging the vertices representing the linking constraints to a single one
    * Try to further partition this graph by applying steps 2-4.

### Parameters
```
  algorithm             should the Stoer-Wagner algorithm or metis be used for finding a minimal cut [TRUE]
  blocksize             number of constraints per block [200]
  enabled               flag to indicate whether detector <cutpacking> is enabled [FALSE]
  fixedblocks           Should the blocks consist of a certain number of constraints [FALSE]
  metisuseptyperb       Should the rb or kway method be used for partitioning by metis [TRUE]
  metisverbose          Should the metis output be displayed [FALSE]
  priority              priority of detector <cutpacking> [1100]
  randomseed            random seed for hmetis [1]
  skip                  flag to indicate whether detector <cutpacking> should be skipped if others found decompositions [FALSE]
  tidy                  Whether to clean up temporary files [TRUE]
  ubfactor              Unbalance factor for metis [5.0]
```
### Future Work

### Known Issues

The detector is still quite slow :snail: on larger instances, see issue #78, but there is no big hope to fix this.

### Links
 * Documentation: dec_cutpacking.cpp
 * [[1] hMETIS - Hypergraph & Circuit Partitioning](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview)
 * [[2] Mechthild Stoer, Frank Wagner: A simple min-cut algorithm](http://dl.acm.org/ft_gateway.cfm?id=263872&ftid=11827&dwn=1&CFID=576236458&CFTOKEN=91189897)
