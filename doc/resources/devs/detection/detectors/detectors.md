# Detectors {#detectors}

# The GCG Detectors
During the detection, the included detectors will try to **find a suitable decomposition**
for your problem. For more information on how they are called, please visit the page on
@ref detection-process.

## List of Detectors
Here you find a list of all detectors available in GCG. Each of those fulfills at least one
function out of
* **propagating** (**assigning** variables/constraints to a block or the master)
* **finishing** (**assigning all remaining** variables/constraints to a block or the master)
* **postprocessing** (**reassigning all remaining** variables/constraints to a block or the master to get a better finished decomposition).

For an overview of which detector does what of those three core functionalities, please check the @subpage detectors-overview.

Detectors based on **constraint or variable classes**:
- @subpage det-consclass
- @subpage det-varclass

Detectors finding **staircase structures**:
- @subpage det-staircase_lsp
- @subpage det-stairheur

Detectors finding **inter- and intra-constraint connections**:
- @subpage det-connected_base
- @subpage det-connected_nonewlinkingvars
- @subpage det-neighborhoodmaster
- @subpage det-densemasterconss

Detectors finding **set structures** (covering, partitioning and packing):
- @subpage det-mastersetcover
- @subpage det-generalmastersetcover
- @subpage det-mastersetpack
- @subpage det-generalmastersetpack
- @subpage det-mastersetpart
- @subpage det-generalmastersetpart

Detectors finding **arrowhead** and **single-bordered structures** (require [hMETIS](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview)):
- @subpage det-hcgpartition
- @subpage det-hrgpartition
- @subpage det-hrcgpartition

Detectors performing **clustering** (require [GSL](https://gist.github.com/TysonRayJones/af7bedcdb8dc59868c7966232b4da903)):
- @subpage det-mst
- @subpage det-dbscan

**Symmetry** detectors (require [bliss](http://www.tcs.hut.fi/Software/bliss/)):
- @subpage det-isomorph 

Detectors used for detection algorithmics:
- @subpage det-compgreedily
- @subpage det-postprocess

<hr>

## Adding own Detectors
If you want to **write your own detector**, i.e. define how you want to use the information 
given by the classifiers to decompose your problem, please consult the "How to" for that.

@ref own-detector
