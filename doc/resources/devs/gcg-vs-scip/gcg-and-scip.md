# Interplay: GCG and SCIP {#gcg-and-scip}

GCG is built on top of SCIP and uses SCIP's methods, data structures and interfaces.
On this page, you can get an overview of the advantages of GCG over SCIP,
as well as differences and the interaction between the two.

First, we give a quick overview **for whom GCG is perfectly suited** and why, also
giving references to our Use Cases.

@subpage why-gcg "Why should you use GCG?"  

\n
There are are some **methods that GCG applies differently** than SCIP that
might increase solving speed for _structured_ programs.

@subpage presolving  

\n
Other times, **GCG needs to interact with SCIP** directly. Note that this can only
be done within the limits of the current SCIP stage.

@subpage stages "SCIP Solving Stages"  \n
@subpage original-vs-transformed \n
@subpage mirroring "Interaction during Branching"  
