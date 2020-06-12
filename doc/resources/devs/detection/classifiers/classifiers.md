# Classifiers {#classifiers}
To be able to decompose a mixed-integer program into blocks, we have to understand
what the underlying structure of it is. This is done using two different types of
classifiers, namely **variable classifiers** and **constraint classifiers**.
As the name suggests, a variable classifier groups the variables into
different categories, e.g. according to its type (binary/integer/continuous),
while a constraint classifier does the same for constraints.

In addition to grouping the variables/constraints into different (disjoint) classes,
it is also possible to precompute how the variables/constraints of a class should be handled in the later detection. This information is stored for each variable/constraint class in `VarClassDecompInfo` and `ConsClassDecompInfo` respectively.

Here you find a list of all classifiers available in GCG and
a short description of their functionality.

- @subpage clsvar "Overview: Variable Classifiers"
- @subpage clscons "Overview: Constraint Classifiers"

\n
If you want to **write your own** classifier, please consult the "How to"
for that.

@ref own-classifier
\n
