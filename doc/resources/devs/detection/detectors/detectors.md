# Detectors {#detectors}

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

Detectors finding **arrowhead** and **single-bordered structures**:
- @subpage det-hcgpartition
- @subpage det-hrgpartition
- @subpage det-hrcgpartition

Detectors used for detection algorithmics:
- @subpage det-compgreedily
- @subpage det-postprocess

Detectors performing clustering (require GSL):
- @subpage det-mst
- @subpage det-dbscan

Other detectors (require bliss):
- @subpage det-isomorph 

<hr>
#### Will be removed in a future version:

Duplicates of the consclass detector:
- @subpage det-constype
- @subpage det-consname

\n
If you want to **write your own** detector, please consult the "How to"
for that.

@ref own-detector\n
