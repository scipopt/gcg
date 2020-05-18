# Detector Functionality Overview {#detectors-overview}

The following list provides an overview over all detectors with their respective
functionality, i.e. whether it can propagate/finish/postprocess. For
more information, see @ref detection.

| ID  |          Full Name          | Propagate | Finish | Postprocess |
|-----|-----------------------------|:---------:|:------:|:-----------:|
| `U` | (***Given by User***)       |   |   |   |
| `g` | compgreedily                | ✓ | ✓ |   |
| `C` | connectedbase               |   | ✓ |   |
| `?` | connected_noNewLinkingVars  | ✓ | ✓ |   |
| `c` | consclass                   | ✓ |   |   |
| `t` | dbscan                      | ✓ |   |   |
| `D` | densemasterconss            | ✓ |   |   |
| `d` | generalmastersetcover       | ✓ |   |   |
| `?` | generalmastersetpack        | ✓ |   |   |
| `?` | generalmastersetpart        | ✓ |   |   |
| `G` | hcgpartition                | ✓ | ✓ |   |
| `a` | hrcgpartition               | ✓ | ✓ |   |
| `r` | hrgpartition                | ✓ |   |   |
| `I` | isomorph                    | ✓ |   |   |
| `?` | mastersetcover              | ✓ |   |   |
| `?` | mastersetpack               | ✓ |   |   |
| `?` | mastersetpart               | ✓ |   |   |
| `M` | mst                         | ✓ |   |   |
| `n` | neighborhoodmaster          | ✓ |   |   |
| `p` | postprocess                 |   |   | ✓ |
| `S` | staircase_lsp               | ✓ | ✓ |   |
| `s` | stairheur                   | ✓ |   |   |
| `v` | varclass                    | ✓ |   |   |
