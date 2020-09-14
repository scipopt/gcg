# Getting Started: Developers Edition {#dev-getting-started}
> <div style="width:80%">
> <img src="expert.png" style="vertical-align:middle; height:6%; position:absolute; right:40px;">
> With this getting started, we do not want to explain how to use GCG and its features (this was
> done in the @ref users), but rather how GCG is implemented. \n
> GCG is built on top of SCIP and uses SCIP's methods, data structures and interfaces and thus 
> it is necessary to understand **how SCIP and GCG generally do things**
> and secondly, **how GCG communicates with SCIP**.
> </div>

@todo we should describe the general *architecture* of GCG with the *original* instance and an *extended* instance which "mirrors" the branching decisions done in the original instance, etc.
@todo we should describe how the original problem is reformulated

# Getting Started as a Developer
There are special facets to the **way that things are generall implemented**
in SCIP and GCG and also some **methods that GCG applies differently** 
than SCIP to try to increase solving speed for _structured_ programs.
All of these are necessary to know for experts and developers and are 
introduced and explained briefly on the following page.


## Implementation Style Overview
In this section, we explain the characteristics that SCIP has that are most important to know about.

### Coding Style Guidelines
Both SCIP and GCG (aim to) comply with a **common set of coding style guidelines**. Those are given by
the [SCIP documentation](https://www.scipopt.org/doc/html/CODE.php).

### SCIP Stages
At times, **GCG needs to interact with SCIP** directly. This can only be done within the limits of the current SCIP stage, because
the solving process within SCIP is executed in stages (see Figure 1). For more information, please check the
SCIP documentation or the [SCIP intro presentation](https://www.scipopt.org/workshop2018/SCIP-Intro.pdf).

\image html SCIP-stages.png "Figure 1: A diagram showing the stages that SCIP works in." width=50% 

@todo add GCG/SCIP stages/interaction from GCG presentation slides

### Original and Transformed Problems {#original-vs-transformed}
As you read in your instance, it **will be kept in SCIP and GCG as the "original" problem**. 
Everything you do to it after reading in is performed on the "transformed" problem 
(presolving is applied on the "transformed" one). 
The original problem is **used as a safe copy** to check the feasibility of solutions. In particular,
it cannot be manipulated.
GCG is detecting on the transformed (i.e. also presolved) problem (`opt`), but can also detect on the original 
(`detect` without `presolve` before it).

## Interplay between GCG and SCIP {#gcg-and-scip}
### Mirroring of Branching Decisions to SCIP {#mirroring}
One of the core features of GCG, the generic column generation, leads to the
fact that GCG sometimes wants to branch differently than SCIP wants to.
This is why we synchronize the branch-and-bound tree between the underlying
SCIP instance and GCG, such that SCIP can execute them.

The code for the communication to SCIP during branching on original variables is inside the
cons_masterbranch.c and cons_origbranch.c source files. The process is as follows:
We (since we can make better branching decisions in most cases) branch ourselves
(`cons_masterbranch`) and then mirror those decisions to SCIP (`cons_origbranch`)
where they are reconstructed.

In the case that an aggregation took place, we do not do branching on original variables.

@todo add stuff from presentations?
