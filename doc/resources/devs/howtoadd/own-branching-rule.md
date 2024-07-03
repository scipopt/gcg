# How to add branching rules & selection heuristics{#own-branching-rule}

[TOC]

Branching rules are used to **intelligently choose variables** to speed up the overall branching process of a problem. They add custom constraints, i.e. introduce
new subproblems that contain more constraints instead of e.g. only dichotomous branching. Branching candidate selection heuristics are used to select one of the possible candidates for a given branching rule.
\n
A complete list of all branching rules contained in this release can be found @ref branching-rules "here". A complete list of all branching candidate selection heuristics contained in this release can be found @ref branching-selection-heuristics "here".

# Adding your own Branching Rule
Before starting to add your branching rule, we strongly encourage you to understand how branching works inside a branch-and-cut context and
in particular how GCG handles two constraint handlers, one for the original problem and one for the master problem (see @ref GCG_BRANCHING).

With the following steps, we explain how you can **add your own branching rule plugin**:
1. **Preparations**
  1. Choose a name `mybranchrule` for your branching rule.
  2. Copy the template files `src/branch_xyz.cpp` and `src/branch_xyz.h`
   while renaming `xyz` to `mybranchrule`.
  3. Open the new files with a text editor and replace all occurrences of `Xyz` by `Mybranchrule`.
2. **Creating your Branching Rule**
  1. Adjust the properties of the branching rule (see @ref BRANCHRULE_PROPERTIES).
  2. [optional] Define the branching rule data (see @ref BRANCHRULE_DATA).
  3. Implement the interface methods (see @ref BRANCHRULE_INTERFACE).
  4. Implement the fundamental callback methods (see @ref BRANCHRULE_FUNDAMENTALCALLBACKS).
  5. [optional] Implement the additional callback methods (see [additional callbacks (SCIP docu)](https://scipopt.org/doc/html/BRANCH.php#BRANCHRULE_ADDITIONALCALLBACKS)).
3. **Make GCG use it**
  1. Add it to gcgplugins.c by adding
    1. the line <tt>\#include branch_mybranchrule.h</tt> in the `/* GCG specific stuff */` section.
    2. the line `SCIP_CALL( SCIPincludeBranchMybranchrule(scip) );`.
  2. the branching rule and the additional callbacks have to be registered by a call of GCGrelaxIncludeBranchrule(), as it is done in the BRANCHINIT callback of branch_mybranchrule.
  3. Add it to your build system:
    1. _Using Makefile:_ Adjust your Makefile such that these files are compiled and linked to your project by adding your branching rule with ending `.o` (`branch_mybranchrule.o`) to the list under `LIBOBJ =` in the file `Makefile` in the root folder.
    2. _Using CMake:_ In `src/CMakeLists.txt`, add your `branch_mybranchrule.cpp` below `set(gcgsources` and your `branch_mybranchrule.h` below the line `set(gcgheaders`.


Additional documentation for the GCG-specific callback methods of branching rules, in particular for their input parameters,
can be found in the file type_branchgcg.h.

# Adding your own Selection Heuristic to an Existing Branching Rule
Unlike branching rules, there is no specific template to add selection heuristics, as they are usually implemented together with the corresponding branching rule. Hence, how to add a new selection heuristics depends on the branching rule to which it is added. However, the selection process is usually implemented in ``branchExeclpbranchrule``, ``branchExecextbranchrule`` or ``branchExecpsbranchrule``.

## Strong Branching
Several strong branching-based heuristics are already implemented inside the file ``branch_bpstrong``. If your new selection heuristic uses strong branching as well, it might be worthwhile to extend the existing implementation instaed of starting from scratch. Naturally, the specific steps again depend on your selection heuristic, however the following might be useful:

**New heuristics/additional features**

* If you want to add a new "quick" heuristic to the hybrid strong branching based heuristics, it likely suffices to add a corresponding case to the ``score_function`` and a corresponding flag whenever the score function is called in main loop of the function ``selectCandidate``.
* If you want to add a heuristic that changes which heuristic should evaluate a candidate at a given point in time, you probably need to alter something inside the main loop of the function ``selectCandidate``.
* Features that change the way that strong branching works need to be added to ``executeStrongBranching``.

**Support for a new branching rule**

In order to add connect the strong branching implementation with a new branching rule, you need to add an interface method, a method that constructs probing child nodes, _potentially_ a new way to store information about the candidates, and _potentially_ additional features specific to the branching rule:
* The interface method is the method that the branching rule calls when it wants to use one of the heuristics. It mainly needs to intialize some of the parameters and provide the set of candidates to the method `selectCandidate`. See e.g. ``GCGbranchSelectCandidateStrongBranchingRyanfoster`` or ``GCGbranchSelectCandidateStrongBranchingOrig``.
* To perform strong branching, we need to know how the potential child nodes are created to compute their scores. In particular, ``executeStrongBranching`` needs to create a new probing node for each potential child node. See ``executeStrongBranching`` and e.g. ``newProbingNodeRyanfosterMaster``.
* Certain information about (the variables corresponding to) candidates needs to be stored over multiple selection rounds, e.g. the reevaluation age. For this purpose, (the variables corresponding to) candidates are stored inside a hash map. How their identifier inside this hash map looks like depends again on the branching rule. We can already handle branching rules where candidates are identified by 1 or 2 candidates, up to a certain point. If this is not sufficient, you need to extend the method ``buildIdentifier``, or possibly change the way that this information is stored.
* Certain behavior is specific to the branching rule that is used. To identify the current branching rule, the ``branchruleData`` in this implementation stores the "initiator", identified by a fixed number. For an impression of how this is done for Ryan-Foster branching, search for "RYANFOSTER" inside the document.

# Introduction to Branching in GCG {#GCG_BRANCHING}
Branching in the branch-cut-and-price context is a bit more complicated than in a branch-and-cut algorithm.
GCG manages two SCIP instances, one for the original instance and one for the reformulated instance, that are solved jointly.
The original SCIP instance coordinates the solving process, while the reformulated instance builds the tree in the same
way, with each node representing the Dantzig-Wolfe reformulated subproblem of the related node in the original instance.

Therefore, if you want to implement a new branching rule, it has to be added to the original SCIP instance, but it gets some
more callback methods, e.g., to transfer the branching decision to the corresponding nodes of the master problem and
enforce these changes during propagation.
\n

The linking of the nodes of the two trees is done via two constraint handler, cons_origbranch.h and cons_masterbranch.h,
that add local constraints to the nodes which know about the corresponding constraint in the other SCIP instance and by
this also the corresponding node. Therefore, each branching rule in GCG has to add one of the origbranch constraints to each
child node it creates. This origbranch constraint also stores a branching data, as described @ref BRANCHRULE_DATA "below", that
can be used to store information about the branching decisions.

```C
SCIP_NODE* childleft;
SCIP_NODE* childright;
GCG_BRANCHDATA* branchleftdata;
GCG_BRANCHDATA* branchrightdata;
SCIP_CONS* origbranchleft;
SCIP_CONS* origbranchright;

/* create the branching data for the children */
SCIP_CALL( SCIPallocMemory(scip, &(branchleftdata)) );
SCIP_CALL( SCIPallocMemory(scip, &(branchrightdata)) );

/* add branching decision to the subproblems and store them in the branching data */
(...)

/* create the origbranch constraints */
SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchleft, "left", childleft,
      GCGconsOrigbranchGetActiveCons(scip), branchrule, branchleftdata) );
SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchright, "right", childright,
      GCGconsOrigbranchGetActiveCons(scip), branchrule, branchrightdata) );

/* add constraints to nodes */
SCIP_CALL( SCIPaddConsNode(scip, childleft, origbranchleft, NULL) );
SCIP_CALL( SCIPaddConsNode(scip, childright, origbranchright, NULL) );

/* release constraints */
SCIP_CALL( SCIPreleaseCons(scip, &origbranchleft) );
SCIP_CALL( SCIPreleaseCons(scip, &origbranchright) );
```

## Properties of a Branching Rule {#BRANCHRULE_PROPERTIES}
At the top of the new file `branch_mybranchrule.cpp`, you can find the branching rule properties.
These are given as compiler defines.
The properties you have to set have the following meaning:

At the top of the new file "branch_mybranchingrule.c" you can find the branching rule properties.
These are given as compiler defines.
In the C++ wrapper class, you have to provide the branching rule properties by calling the constructor
of the abstract base class scip::ObjBranchrule from within your constructor.
The properties you have to set have the following meaning:

\par BRANCHRULE_NAME: the name of the branching rule.
This name is used in the interactive shell to address the branching rule.
Additionally, if you are searching for a branching rule with SCIPfindBranchrule(), this name is looked up.
Names have to be unique: no two branching rules may have the same name.

\par BRANCHRULE_DESC: the description of the branching rule.
This string is printed as a description of the branching rule in the interactive shell.

\par BRANCHRULE_PRIORITY: the default value for the priority of the branching rule.
In the subproblem processing, the branching rules are called in decreasing order of their priority until
one succeeded to branch. Since most branching rules are able to generate a branching in all situations,
only the rule of highest priority is used. In combination with the BRANCHRULE_MAXDEPTH and
BRANCHRULE_MAXBOUNDDIST settings, however, interesting strategies can be easily employed. For example,
the user can set the priority of the "full strong branching" strategy to the highest value and assign the
second highest value to the "reliable pseudo cost" rule. If one also sets the maximal depth for the
"full strong branching" to 5, in the top 5 depth levels of the search tree the "full strong branching" is
applied, while in the deeper levels "reliable pseudo cost branching" is used.
\n
Note that the BRANCHRULE_PRIORITY property only specifies the default value of the priority. The user can
change this value arbitrarily.

\par BRANCHRULE_MAXDEPTH: the default value for the maximal depth level of the branching rule.
This parameter denotes the maximal depth level in the branch-and-bound tree up to which the branching method of the
branching rule will be applied. Use -1 for no limit.
\n
Note that this property only specifies the default value. The user can change this value arbitrarily.

\par BRANCHRULE_MAXBOUNDDIST: the default value for the maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for applying branching.
At the current branch-and-bound node, the relative distance from its dual bound (local dual bound)
to the primal bound compared to the best node's dual bound (global dual bound) is considered. The branching method of
the branching rule will only be applied at the node if this relative distance does not exceed BRANCHRULE_MAXBOUNDDIST.
\n
For example, if the global dual bound is 50 and the primal bound is 60, BRANCHRULE_MAXBOUNDDIST = 0.25 means that
branching is only applied if the current node's dual bound is in the first quarter of the interval [50,60], i.e., if it
is less than or equal to 52.5. In particular, the values 0.0 and 1.0 mean that the branching rule is applied at the
current best node only or at all nodes, respectively.
\n
Note that the BRANCHRULE_MAXBOUNDDIST property only specifies the default value of the maximal bound distance.
The user can change this value arbitrarily.


## GCG specific branching data {#BRANCHRULE_DATA}
Below the header "Data structures" you can find the new struct "struct GCG_BranchData".
In this data structure, each branching rule can store information about the branching decisions applied to a newly created node.
Later, this information is used to transfer the branching decisions to the corresponding nodes in the reformulated SCIP instance.
Defining branching data is optional. You can leave this struct empty.


## Interface Methods {#BRANCHRULE_INTERFACE}
At the bottom of "branch_mybranchingrule.c", you can find the interface method SCIPincludeBranchruleMybranchingrule(),
which also appears in "branch_mybranchingrule.h"
SCIPincludeBranchruleMybranchingrule() is called by the user, if one wants to include the branching rule,
i.e., if one wants to use the branching rule in their application.

This method only has to be adjusted slightly.
It is responsible for notifying SCIP of the presence of the branching rule. For this, you can either call
SCIPincludeBranchrule(),
or SCIPincludeBranchruleBasic() since SCIP version 3.0. In the latter variant, [additional callbacks (SCIP docu)](https://scipopt.org/doc/html/BRANCH.php#BRANCHRULE_ADDITIONALCALLBACKS)
must be added via setter functions as, e.g., SCIPsetBranchruleCopy(). We recommend this latter variant because
it is more stable towards future SCIP versions which might have more callbacks, whereas source code using the first
variant must be manually adjusted with every SCIP release containing new callbacks for branchrule in order to compile.


If you are using branching rule data, you have to allocate the memory for the data at this point.
You can do this by calling:
```C
SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );
```
You also have to initialize the fields in struct SCIP_BranchruleData afterwards.

You may also add user parameters for your branching rule, see the method SCIPincludeBranchruleRelpscost() in
src/scip/branch_relpscost.c for an example.


## GCG-specific callbacks of branching rules {#BRANCHRULE_FUNDAMENTALCALLBACKS}
### BRANCHACTIVEMASTER
The BRANCHACTIVEMASTER callback is called whenever a node in the master problem is activated.
It should transfer the branching decisions that were made at the corresponding node in the original SCIP instance to the
current node of the reformulated SCIP instance. Therefore, it gets branching data stored at the corresponding node in the
original SCIP instance as a parameter. It should then for example reformulate a banching constraint to the
reformulated problem space and add it locally to the current node or update the pricing problem according to the
branching decision.

### BRANCHDEACTIVEMASTER
The BRANCHDEACTIVEMASTER callback is called whenever a node in the master problem is deactivated.
It should undo all changes that were performed in the BRANCHACTIVEMASTER callback, e.g., remove branching restrictions
from the pricing problems.

### BRANCHPROPMASTER
The BRANCHPROPMASTER callback is called whenever a node in the master problem is propagated by the masterbranch
constraint handler.
It should perform propagation according to the branching changes done at the related original node, which are stored in
the given banching data. In particular, it should fix all variable locally to 0, that contradict the branching decision
and represent an extreme point or extreme ray which is not feasible for the current pricing problem.

### BRANCHMASTERSOLVED
The BRANCHMASTERSOLVED callback is called whenever the relaxation of a node is solved for the first time. It can be used
to collect statistical information about branching decision, e.g., update pseudocost values.

### BRANCHDATADELETE
The BRANCHDATADELETE callback is called when a branch-and-bound node is freed and its branching data should be freed, too.

### BRANCHNEWCOL
The BRANCHNEWCOL callback is called whenever a new column is added to the master problem. It can be used to update any branching constraints that were added to the master problem, e.g. to determine and add the coefficient of the new master variable to the corresponding branching constraint.

### BRANCHGETMASTERCUT
If the branching rule branches in the master problem, generating a constraint that can be enforced only if the pricing problem is modified, then the BRANCHGETMASTERCUT callback will be used by GCG to retrieve the master branching constraint as well as the corresponding modifications to be applied in the pricing problems.