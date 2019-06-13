# the model is taken from
# R. Li, S. Hu, P. Zhao, Y. Zhou and M. Yin
# A novel local search algorithm for the minimum capacitated dominating set
# Journal of the Operational Research Society, 2017
# DOI: 10.1057/s41274-017-0268-6

# the example graph originally appears in
# A. Potluri and A. Singh
# A greedy heuristic and its variants for minimum capacitated dominating set
# In M. Parashar et al. (Eds.): IC3 2012, CCIS 306, pp. 28–39, Springer, 2012.
# DOI: 10.1007/978-3-642-32129-0_9


# vertices 
set V := {1 .. 14};
# (undirected) edges
set E := {<1,2>, <2,3>, <2,4>, <4,6>, <3,8>, <7,8>, <5,6>, <6,9>, <6,10>, <6,11>, <7,11>, <7,12>, <7,13>, <7,14>, <10,12>, <11,12>};
# (directed) arcs, produced from (undirected) edges
set A := E union {<j,i> in V*V with <i,j> in E};

# uniform capacity of 2 for each vertex
param cap[<v> in V] := 2;

# neighbors of v
defset N(v) := {<i> in V with <v,i> in A};


# is v selected in the dominating set?
var x[V] binary;

# for <i,j> in A: is j dominated by i?
var y[A] binary;


minimize cost: sum<v> in V: x[v];

# every vertex v must be dominated by some vertex i (unless v is selected)
subto dominate:
	forall <v> in V: sum<i> in N(v): y[i,v] >= 1 - x[v];

# the domination capacity of each vertex must be respected
subto capacity:
	forall <v> in V: sum<i> in N(v): y[v,i] <= x[v]*cap[v];

# link the y and the x variables
subto link:
	forall <v> in V: forall <i> in N(v): y[v,i] <= x[v];
