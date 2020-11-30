# The GCG Detection {#detection}
> The **detection of structures** in a given model is one of the main features that
> distinguishes GCG from SCIP. It accelerates the solving speed of many problems.
> For further reading, please consult the @ref use-cases "Use-Cases", in particular @ref u2 "the second one".

To apply Branch-and-Price efficiently, it is **required that your problem is structured**. 
This means that it can be decomposed by our detection mechanism. 

@subpage structure-types "List of different structures GCG can find" \n
\n

Since we want to find these structures, we let our classifiers run on your model to 
check for specific classes of constraints and variables. This information is further used
by different detectors to finally complete decompositions.
In the following, you can find pages describing **how exactly this process is carried out**
as well as to all **classifiers and detectors** that GCG offers. 
For each, we will present what exact structure it finds (or what classes, respectively), 
how it does that and for what instances this might be suited.

@subpage detection-process \n
@subpage classifiers \n
@subpage detectors \n
\n

When trying to find the best possible decomposition, you will stumble across different 
**detection scores** and a different mode for detecting structures, the **Benders Mode**. 

@subpage detection-scores \n
@subpage benders \n
\n

If you want to tweak the detection even more (see also @ref u5 for a use case on that subject),
you can also change **detection parameters** manually. 

@subpage detection-params \n

<hr>

If you want to **write your own classifiers and detectors**, i.e. define how you want
to decompose your problem, please consult the "How to" for that.

@ref own-classifier \n
@ref own-detector \n
