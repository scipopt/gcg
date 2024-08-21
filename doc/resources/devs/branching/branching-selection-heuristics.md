# Branching Candidate Selection Heuristics {#branching-selection-heuristics}

# The GCG Branching Candidate Selection Heuristics
The branching candidate selection heuristic decides **which of the candidates** provided by the branching rule is
selected. As this influences the size of the branch-and-bound tree, it can be highly influential on the runtime of the
overall algorithm.

## List of Selection Heuristics
 As the the type of candidates can differ between the branching rules, the concrete implementation of the heuristics
 changes as well. Which heuristics are implemented for which branching rules in GCG is listed in the following table,
 more detailed explanations can be found below:

<center>
|                                                             | Original Variable | Ryan-Foster     | Vanderbeck Generic | Component Bound |
|:------------------------------------------------------------|:-----------------:|:---------------:|:------------------:|:---------------:|
| Random Branching                                            | X                 | X               | X                  | X               |
| Most Fractional/Infeasible <br> Branching                   | X                 | (X)<sup>2</sup> |                    |                 |
| Pseudocost Branching                                        | X                 | (X)<sup>2</sup> |                    |                 |
| Strong Branching <br> with Column Generation<sup>1</sup>    | X                 | X               |                    |                 |
| Strong Branching <br> without Column Generation<sup>1</sup> | X                 | X               |                    |                 |
| Hybrid Branching<sup>1</sup>                                | X                 | (X)<sup>3</sup> |                    |                 |
| Reliability Branching<sup>1</sup>                           | X                 | (X)<sup>3</sup> |                    |                 |
| Hierarchical Branching<sup>1</sup>                          | X                 | (X)<sup>3</sup> |                    |                 |
</center>
<sup>1</sup> : <small>GCG offers a highly customizable implementation of strong branching-based heuristics (currently
in an experimental state), which, among other things, allows to combine them. Refer to this page for more information.</small>
<sup>2</sup> : <small>GCG can only aggregate the respective scores of the individual variables.</small>
<sup>3</sup> : <small>These heuristics where originally meant to use pseudocost branching, see <sup>2</sup>. However,
pseudocost branching can also be substituted by any other heuristic, with varying performance.</small>

## Discussion of Available Selection Heuristics
The goal of any selection heuristic is to select candidates in a way that will result in the **smallest branch-and-bound
tree**. As determining these candidates is (to our current knowledge) computationally infeasible, the heuristics instead
pick candidates that optimize some other criterion that is potentially thought to be a good indicator of the resulting
tree size.

### Random Branching
Randomly selects a candidate.

### Most Fractional/Infeasible Branching
Selects one of the candidates where the distance of the current solution value of the corresponding variable to the
next candidate is maximal.

<hr>

> The following heuristics all try to maximize the aggregated _gain in dual bound_ in the child nodes that would be
> created if a given candidate was branched on. Letting \f$ z \f$ and \f$ z_i \f$ be the dual bound of the current node
> and child node \f$ i \f$, respectively, the gain of child node \f$ z_i \f$ is equal to \f$ z - z_i \f$.

### Pseudocost Branching
> Pseudocost branching was designed to be used in conjunction with branching on (original/master) variables, and will
> likely not be as effective if used with other branching rules.

Pseudocost branching estimates the gain by **averaging previous** nodes where a given candidate was branched on: given some
candidate that introduces bounds on a single variable \f$ x_i \f$, the projected up gain ("upwards pseudocost") is
\f{align}{
    \Phi^+(x_i) := \frac{1}{|I^+_i|}\sum_{i \in I^+_i}\frac{z'_{u_i(y)} - z'_y}{\bar{x}_i - \lfloor\bar{x}_i\rfloor},
\f}
where \f$ I^+_i \f$ is the set of nodes where \f$ x_i \f$ already was branched on and the up branch has already been
solved, \f$ u_i(y) \f$ is the up branch when branching on \f$ x_i \f$ in node \f$ y \f$, \f$ z'_y \f$ is the local
lower bound of node \f$ y \f$, and \f$ \bar{x}_i \f$ is the solution value of \f$ x_i \f$ in the current node. The
projected down gain can be calculated analogously.

### Strong Branching (with or without Column Generation)
Strong branching estimates the gain by actually **solving the respective child nodes**. _Fully_ evaluating _all_ child
nodes is called (all) full strong branching. We additionally differentiate between strong branching _with_ and strong
branching _without_ column generation, i.e. strong branching can be relaxed by not performing column generation to
evaluate the child nodes. Using strong branching generally results in **smaller trees**, however its significant
computational effort means that it is usually not worth it to use just strong branching.

### Hybrid Pseudocost/Strong Branching
Hybrid pseudocost/strong branching uses strong branching for **nodes up to a given depth**, and pseudocost branching for
all nodes after that. GCG can substitute any of the aformentioned heuristics for the latter, and both strong branching
with and without column generation for the former. Furthermore, GCG can dynamically adjust the height parameter based on the total number of binary and integer variables.

### Reliability Branching
Reliability branching uses strong branching for some candidate \f$ if \f$ if its pseudocost values are not _reliable_
enough yet, i.e. if the size of either I^+_i or I^-_i (see the section on pseudocost branching) is below a given
threshold. Pseudocost branching is applied to all other candidates. GCG can substitute any of the aformentioned
heuristics (except hybrid branching) for the latter, and both strong branching with and without column generation for
the former.

### Hierarchical Branching
In hierarchical branching, the selection process is divided into **three phases**, where
- in phase 0, the candidates are filtered based on some heuristic that is quick to compute (e.g. pseudocost branching),
then
- in phase 1, the remaining candidates are filtered based their strong branching without column generation scores, and
finally
- in phase 2, a candidate is selected out of the remaining candidates based on the score the candidates receive from
strong branching with column generation.

## Parameters for Branching Candidate Selection Heuristics
Some of the selection heuristics, in GCG primarily those that are based on strong branching, can further be adjusted
by **setting values for certain parameters**. These parameters include the depth until which strong branching is performed
for hybrid branching, the reliability parameter for reliability branching, or the number of candidates that pass
through each phase for hierarchical branching. See @ref selection-heuristics-params "here" for a full list of the available parameters, together with
descriptions, value ranges and default values.
While in GCG, the parameters can individually be changed before optimizing a given problem, the best way to do so is by
loading a settings file. This can either be done from inside GCG by entering
 ```
set load [path/to/settings/file].set
 ```
 or when starting GCG from the command line by adding the option
```
-s "[path/to/settings/file].set"
```
You can find further information on settings files @ref getting-started "here".
GCG already has several **preconfigured settings files** for the heuristics introduced on this page, which can be found in
the settings folder, where ``orig`` refers to original variable branching and ``rf`` to Ryan-Foster branching:
```
[orig/rf]Random              : random branching
origMostFrac                 : most fractional/infeasible branching
origPseudo                   :
[orig/rf]FullSBCG            : full strong branching with column generation
[orig/rf]FullSBnoCG          : full strong branching without column generation
[orig/rf]HybridCGFixedParams : hybrid pseudocost/strong branching with column generation with fixed depth parameters
[orig/rf]HybridCGDynParams   : hybrid pseudocost/strong branching with column generation with dynamic depth parameters
[orig/rf]ReliabilityCG       : reliability branching with column generation
[orig/rf]Hierarchical        : hierarchical strong branching
[orig/rf]HybridHierarchical  : hierarchical strong branching combined with hybrid branching
rfHierarchicalNoP0           : hierarchical strong branching without the first phase
```
These files might provide a good starting point for the configuration you wish to achieve. Additionally, the file
``strongBranchingParams.set`` again contains all parameters together with descriptions, value ranges and default
values.
As you might have noticed, the above list also contains settings files that combine multiple of the aforementioned heuristics (``[orig/rf]HybridHierarchical.set``). This is possible because all of the strong branching-based heuristics are implemented inside the same method. Combining heuristics this way can be beneficial for some instances.

<hr>

## Adding own Branching Candidate Selection Heuristics
If you want to **add your own branching candidate selection heuristics**, i.e. define exactly which candidates are selected
please consider our "How to add" for that.

@ref own-branching-rule