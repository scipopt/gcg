# Classifiers {#classifiers}
To be able to decompose a mixed-integer program into blocks, we have to understand
what the underlying structure of it is. This is done using two different types of
classifiers, namely **variable classifiers** and **constraint classifiers**.
As the name suggests, a variable classifier groups the variables into
different categories, e.g. according to its type (binary/integer/continuous),
while a constraint classifier does the same for constraints.

There are generally three kinds of variables and constraints that we want the classification to end up with:
* Linking variables/constraints,
* Block variables/constraints and
* Master variables/constraints

Here you find a list of all classifiers available in GCG and
a short description of their functionality.

- @subpage clsvar "Overview: Variable Classifiers"
- @subpage clscons "Overview: Constraint Classifiers"
- @subpage clsindex "Overview: Index Classifiers"
