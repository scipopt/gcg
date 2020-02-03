# The GCG Detection {#detection}

> The **detection of structures** in a given model is one of the main features that
> distinguishes GCG from SCIP. It accelerates the solving speed of many problems.
> For further reading, please consult the @ref use-cases "Use-Cases", in particular @ref u2 "the second one".

There are many different kinds of constraints and variables inside a model and each has a certain
form. Since we want to use those, we let our **classifiers** run on your model to check
for specific classes of constraints and variables. Then, there are many structures
that can be found in models and for each type, we have a **detector** to find (detect) it.
In the following, you can find the pages that lead to the description of all classifiers and
detectors that GCG offers. For each, we will present what exact structure it finds
(or what classes, respectively), how it does that and for what instances this might be suited.

@subpage classifiers "Classifiers: Overview Page"\n
@subpage detectors "Detectors: Overview Page"\n

\n
Our **detection** is executed in a **loop**, where all detectors are called one by one.
Furthermore, one can also activate an alternative different detection mode, namely the
**benders detection**. A description of both is given on the following subpages.

@subpage detection-process\n
@subpage benders\n

\n
If you want to **write your own** classifiers and detectors, please consult the "How to"
for that.

@ref own-classifier\n
@ref own-detector\n
