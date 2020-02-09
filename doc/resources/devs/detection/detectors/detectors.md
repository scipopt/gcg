# Detectors {#detectors}

> **Important: The documentation here is fully based on the refactoring branch as of commit 6237731b.**  \n
> Thus, links to source files might not work, since the documentation branch branches from the master.

Here you find a list of all detectors available in GCG. Each of those fulfills at least one
function out of
* **propagating** (**assigning** variables/constraints to a block or the master)
* **finishing** (**assigning all remaining** variables/constraints to a block or the master)
* **postprocessing** (**reassigning all remaining** variables/constraints to a block or the master to get a better finished decomposition).

For an overview of which detector does what of those three core functionalities, please check the @subpage detectors-overview.

## Detectors
Detectors based on constraint or variable classes:
- @subpage det-consclass
- @subpage det-varclass


Detectors for staircase structures:
- @subpage det-staircase_lsp
- @subpage det-stairheur


Detectors performing clustering:
- @subpage det-mst
- @subpage det-dbscan


Detectors finding **connections** between constraints via variables:
- @subpage det-connected_base
- @subpage det-connected_nonewlinkingvars


Currently uncategorized detectors:
- @subpage det-densemasterconss
- @subpage det-generalmastersetcover
- @subpage det-generalmastersetpack
- @subpage det-generalmastersetpart
- @subpage det-hcgpartition
- @subpage det-hrcgpartition
- @subpage det-hrgpartition
- @subpage det-isomorph
- @subpage det-mastersetcover
- @subpage det-mastersetpack
- @subpage det-mastersetpart
- @subpage det-neighborhoodmaster


Detectors used for detection algorithmics:
- @subpage det-compgreedily
- @subpage det-postprocess


Detectors that will be removed after the refactoring or are already removed:
- @subpage det-staircase
- @subpage det-colors
- @subpage det-cutpacking
- @subpage det-connected
- @subpage det-random

Detectors that will be removed:
- @subpage det-constype
- @subpage det-consname
