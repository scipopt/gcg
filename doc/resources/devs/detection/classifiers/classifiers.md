# Classifiers {#classifiers}

# The GCG Classification
To be able to decompose a mixed-integer program into blocks, we have to understand
what the underlying structure of it is. This is done using two different types of
classifiers, namely **variable classifiers** and **constraint classifiers**.
As the name suggests, a variable classifier groups the variables into
different categories, e.g. according to its type (binary/integer/continuous),
while a constraint classifier does the same for constraints.

In addition to grouping the variables/constraints into different (disjoint) classes,
it is also possible to precompute how the variables/constraints of a class should be 
handled in the later detection. This information is stored for each variable/constraint class 
in `VarClassDecompInfo` and `ConsClassDecompInfo` respectively.

### Variable Classifiers
Variable classifiers can store in `VarClassDecompInfo` how the variables of a class should be handled by a detector:
 * `ALL` = 0 -> not specified
 * `LINKING` = 1 -> as linking variables
 * `MASTER` = 2 -> as master variables
 * `BLOCK` = 3  -> as block variables

 ### Constraint Classifiers 
Constraint classifiers can store in `ConsClassDecompInfo` how the constraints of a class should be handled by a detector:
 * `BOTH` -> not specified
 * `ONLY_MASTER` -> assigned to master problem
 * `ONLY_PRICING` -> assigned to pricing problem

## List of Classifiers
Here you find a list of all classifiers available in GCG and
a short description of their functionality. They are divided
into variable and constraint classifiers.

The following **variable classifiers** are available:
- @subpage clsvar_scipvartypes
- @subpage clsvar_objvalues
- @subpage clsvar_objvaluesigns

The following **constraint classifiers** are available:
- @subpage clscons_scipconstypes
- @subpage clscons_miplibconstypes
- @subpage clscons_nnonzeros
- @subpage clscons_consnamenonumbers
- @subpage clscons_consnamelevenshtein

<hr>

## Adding own Classifiers
If you want to **write your own classifier**, i.e. define by what criteria you want to classify
variables or constraints, please consult the "How to" for that.

@ref own-classifier
