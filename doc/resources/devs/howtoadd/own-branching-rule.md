# Your own Branching Rule (deprecated) {#own-branching-rule}

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
 child node it creates. This origbranch constraint also stores a branching data, as described \ref BRANCHDATA "below", that
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

 A complete list of all branching rules contained in this release can be found "here" (please add).

 In the following, we explain how the user can add his own branching rule.
 Start by copying the template files src/branch_xyz.c and src/branch_xyz.h into files named "branch_mybranchrule.c"
    and "branch_mybranchrule.h".

 Take the branching rule that branches on original variables (src/branch_orig.c) as an example.
 As all other default plug-ins, it is written in C. There is currently no C++ wrapper available.

 In the following, we will only explain the additional callback methods which a branching rule in GCG can implement.
 For a general introduction about how to write a branching rule in SCIP and a description of the default callback methods
 of a branching rule in SCIP, we refer to the SCIP documentation and the "How to add constraint handlers" there.

 Additional documentation for the GCG-specific callback methods of branching rules, in particular for their input parameters,
 can be found in the file type_branchgcg.h.

 Here is what you have to do to implement a branching rule in GCG:
 - Copy the template files src/branch_xyz.c and src/branch_xyz.h into files named "branch_mybranchrule.c"
    and "branch_mybranchrule.h".
    \n
    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 - Open the new files with a text editor and replace all occurrences of "xyz" by "mybranchrule".
 - Adjust the properties of the branching rule (see SCIP manual).
 - Define the branching rule data (see SCIP manual) and the branching data (see \ref BRANCHDATA). This is optional.
 - Implement the interface methods and the callback methods as described in the SCIP manual.
 Besides including the branching rule in SCIP via SCIPincludeBranchrule(), the branching rule and the additional callbacks
 have to be registered by a call of GCGrelaxIncludeBranchrule(), as it is done in the BANCHINIT callback of branch_mybranchrule.
 - Implement the GCG-specific callback methods (see \ref BRANCH_CALLBACKS).

 # GCG specific branching data {#BRANCHDATA}

 Below the header "Data structures" you can find the new struct "struct GCG_BranchData".
 In this data structure, each branching rule can store information about the branching decisions applied to a newly created node.
 Later, this information is used to transfer the branching decisions to the corresponding nodes in the reformulated SCIP instance.
 Defining branching data is optional. You can leave this struct empty.

 # GCG-specific callbacks of branching rules {#BRANCH_CALLBACKS}

 ## BRANCHACTIVEMASTER

 The BRANCHACTIVEMASTER callback is called whenever a node in the master problem is activated.
 It should transfer the branching decisions that were made at the corresponding node in the original SCIP instance to the
 current node of the reformulated SCIP instance. Therefore, it gets branching data stored at the corresponding node in the
 original SCIP instance as a parameter. It should then for example reformulate a banching constraint to the
 reformulated problem space and add it locally to the current node or update the pricing problem according to the
 branching decision.

 ## BRANCHDEACTIVEMASTER

 The BRANCHDEACTIVEMASTER callback is called whenever a node in the master problem is deactivated.
 It should undo all changes that were performed in the BRANCHACTIVEMASTER callback, e.g., remove branching restrictions
 from the pricing problems.

 ## BRANCHPROPMASTER

 The BRANCHPROPMASTER callback is called whenever a node in the master problem is propagated by the masterbranch
 constraint handler.
 It should perform propagation according to the branching changes done at the related original node, which are stored in
 the given banching data. In particular, it should fix all variable locally to 0, that contradict the branching decision
 and represent an extreme point or extreme ray which is not feasible for the current pricing problem.

 ## BRANCHMASTERSOLVED

 The BRANCHMASTERSOLVED callback is called whenever the relaxation of a node is solved for the first time. It can be used
 to collect statistical information about branching decision, e.g., update pseudocost values.

 ## BRANCHDATADELETE

 The BRANCHDATADELETE callback is called when a branch-and-bound node is freed and its branching data should be freed, too.
