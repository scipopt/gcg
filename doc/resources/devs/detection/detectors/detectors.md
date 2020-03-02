# Detectors {#detectors}

> **Important: The documentation here is fully based on the refactoring branch as of commit 6237731b.**  \n
> Thus, links to source files might not work, since the documentation branch branches from the master.

Here you find a list of all detectors available in GCG. Each of those fulfills at least one
function out of
* **propagating** (**assigning** variables/constraints to a block or the master)
* **finishing** (**assigning all remaining** variables/constraints to a block or the master)
* **postprocessing** (**reassigning all remaining** variables/constraints to a block or the master to get a better finished decomposition).

For an overview of which detector does what of those three core functionalities, please check the @subpage detectors-overview.

### Detectors
#### Working detectors
Detectors based on **constraint or variable classes**:
- @subpage det-consclass
- @subpage det-varclass


Detectors finding **staircase structures**:
- @subpage det-staircase_lsp
- @subpage det-stairheur

Detectors finding **connections** between constraints via variables:
- @subpage det-connected_base
- @subpage det-connected_nonewlinkingvars

Detectors finding **set structures** (covering, partitioning and packing):
- @subpage det-mastersetcover
- @subpage det-generalmastersetcover
- @subpage det-mastersetpack
- @subpage det-generalmastersetpack
- @subpage det-mastersetpart
- @subpage det-generalmastersetpart

Detectors finding **arrowhead structures**:
- @subpage det-hcgpartition
- @subpage det-hrcgpartition
- @subpage det-hrgpartition

Currently uncategorized detectors:
- @subpage det-densemasterconss
- @subpage det-neighborhoodmaster

Detectors used for detection algorithmics:
- @subpage det-compgreedily
- @subpage det-postprocess

<hr>
#### Existing in the refactoring, but currently not in use

Detectors performing clustering (currently not in refactoring):
- @subpage det-mst
- @subpage det-dbscan

Other detectors:
- @subpage det-isomorph

<hr>
#### Will be removed after the refactoring or are already removed

Legacy Mode detectors:
- @subpage det-staircase
- @subpage det-colors
- @subpage det-cutpacking
- @subpage det-connected
- @subpage det-random

Duplicates:
- @subpage det-constype
- @subpage det-consname
