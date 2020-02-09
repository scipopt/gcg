# Detector Functionality Overview {#det-function-overview}

The following list provides an overview over all detectors with their respective
functionality, i.e. whether it can propagate/finish/postprocess. For
more information, see @ref detection.

|         Detector (ID)       | Propagate | Finish | Postprocess |
|-----------------------------|:---------:|:------:|:-----------:|
| compgreedily                | ✓ | ✓ |   |
| connectedbase               |   | ✓ |   |
| connected_noNewLinkingVars  | ✓ | ✓ |   |
| consclass                   | ✓ |   |   |
| dbscan                      | ✓ |   |   |
| densemasterconss            | ✓ |   |   |
| generalmastersetcover       | ✓ |   |   |
| generalmastersetpack        | ✓ |   |   |
| generalmastersetpart        | ✓ |   |   |
| hcgpartition                | ✓ | ✓ |   |
| hrcgpartition               | ✓ | ✓ |   |
| hrgpartition                | ✓ |   |   |
| isomorph                    | ✓ |   |   |
| mastersetcover              | ✓ |   |   |
| mastersetpack               | ✓ |   |   |
| mastersetpart               | ✓ |   |   |
| mst                         | ✓ |   |   |
| neighborhoodmaster          | ✓ |   |   |
| postprocess                 |   |   | ✓ |
| staircase_lsp               | ✓ | ✓ |   |
| stairheur                   | ✓ |   |   |
| varclass                    | ✓ |   |   |
