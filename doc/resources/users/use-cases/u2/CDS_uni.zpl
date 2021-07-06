# the model is a corrected version from
# R. Li, S. Hu, P. Zhao, Y. Zhou and M. Yin
# A novel local search algorithm for the minimum capacitated dominating set
# Journal of the Operational Research Society, 2017
# DOI: 10.1057/s41274-017-0268-6


###### this file only works with the input format of the *uniform* capacity instances ######


# this can be overwritten by a -D file=<filename> option
# you can download all instances from https://scis.uohyd.ac.in/~apcs/capmds/
param file := "CAPMDS-data/uni/UDG/2/50/150/graphS0N50"; 

# how many vertices?
param n := read file as "1n" use 1;

# uniform capacity for each vertex
param cap := read file as "2n" use 1;

# vertices
set V := {1 .. n};
# edges
set E := {read file as "<1n,2n>" skip 1};

# set of directed arcs, produced from undirected edges
set A := E union {<j,i> in V*V with <i,j> in E};

# neighbors of v
defset N(v) := {<i> in V with <v,i> in A};


# is v selected in the dominating set?
var x[V] binary;

# for <i,j> in A: is j dominated by i?
var y[A] binary;


minimize cost: sum<v> in V: x[v];

# every vertex i must be dominated by some vertex j
subto dominate:
	forall <v> in V: sum<i> in N(v): y[i,v] >= 1 - x[v];

# the domination capacity of each vertex must be respected
subto capacity:
	forall <v> in V: sum<i> in N(v): y[v,i] <= x[v]*cap;
