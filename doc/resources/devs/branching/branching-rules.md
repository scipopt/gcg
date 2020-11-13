# Branching Rules {#branching-rules}

# The GCG Branching Rules
Branching Rules describe **how one can branch in a given node** in the Branch-and-Bound tree. For example,
one can branch on constraints or on a set of variables all at once. Adjusting the branching rule according
to your problem will be done automatically by GCG and can significantly influence the runtime required
to solve it.

## List of Branching Rules
We currently offer 4 different branching rules, as well as one empty branching rule to coordinate the
branching process between SCIP and GCG.

- **Vanderbeck Generic Branching** \n
branching rule for original problem in GCG while real branching is in the master 
- **Original Variable Branching** \n
branching rule for original problem in GCG 
- **Pseudo Cost Branching** \n
generalized reliable pseudo costs branching rule 
- **Ryan Foster Branching** \n
branching rule for original problem in GCG implementing the Ryan and Foster branching scheme 

Descriptions of those and their code documentation can be
found @ref BRANCHINGRULES "here".

## Adding own Branching Rules
If you want to add your own branching rules, i.e. define exactly **how to branch on variables**, 
please consider our "How to add" for that.

 ⇨ @ref own-branching-rule \n