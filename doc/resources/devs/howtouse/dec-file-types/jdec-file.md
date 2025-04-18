# .jdec file format {#jdec-file}

> This file format was introduced with GCG 4.0 and is still in an experimental state.
> It may be extended in the future.

This page will give a short introduction to the `.jdec` file format by means of an example
and an overview of the structure of the file format.

## Features
This file format was introduced with GCG 4.0. Currently it supports the following features:
* JSON compatible file format (i.e., all general JSON file format rules apply)
* complete and partial decompositions
* decompositions of presolved problems
* block constraints and master constraints
* symmetry information used to aggregate pricing problems
* nested decompositions (decompositions for blocks can be specified, utilized by the @ref pricing-solvers "GCG pricing solver")

The reading and writing process works as for the `.dec` @ref dec-file "file format".
As it is based on JSON, you can read and write it with JSON libraries.

## Example
Do define a decomposition, you need a given problem instance, whose constraints and variables are referred to by name:
```
Minimize
 cost:  + y#0 + y#1 + y#2 + y#3
Subject to
 assign_0:
  + x#0#0 + x#1#0 + x#2#0 + x#3#0 >= 1
 assign_1:
  + x#0#1 + x#1#1 + x#2#1 + x#3#1 >= 1
 link_0:
  - y#0 + y#1 >= 0
 link_1:
  - y#2 + y#3 >= 0
 cap_0:
  -10 y#0 + 2 x#0#0 + 3 x#0#1 <= 0
 cap_1:
  -5 y#1 + 2 x#1#0 + 3 x#1#1 <= 0
 cap_2:
  -10 y#2 + 2 x#2#0 + 3 x#2#1 <= 0
 cap_3:
  -5 y#3 + 2 x#3#0 + 3 x#3#1 <= 0
```
Using the names defined in the model, you can specify a decomposition like that:
```json
{
  "version": 1,
  "name": "example_nested_dec",
  "problem_name": "example.lp",
  "description": "nested decomposition with aggregated blocks",
  "decomposition": {
    "presolved": false,
    "n_blocks": 2,
    "master_constraints": [
      "assign_0"
    ],
    "blocks": [
      {
        "index": 0,
        "constraints": [
          "cap_0",
          "cap_1",
          "link_0"
        ],
        "decomposition": {
          "presolved": false,
          "n_blocks": 2,
          "master_constraints": [
            "link_0"
          ],
          "blocks": [
            {
              "index": 0,
              "constraints": [
                "cap_0"
              ]
            },
            {
              "index": 1,
              "constraints": [
                "cap_1"
              ]
            }
          ]
        },
        "symmetry_representative_block": 0
      },
      {
        "index": 1,
        "constraints": [
          "cap_2",
          "cap_3",
          "link_1"
        ],
        "decomposition": {
          "presolved": false,
          "n_blocks": 2,
          "master_constraints": [
            "link_1"
          ],
          "blocks": [
            {
              "index": 0,
              "constraints": [
                "cap_2"
              ]
            },
            {
              "index": 1,
              "constraints": [
                "cap_3"
              ]
            }
          ]
        },
        "symmetry_representative_block": 0,
      }
    ],
    "symmetry_var_mapping": {
      "y_0": "y_0",
      "x_0_0": "x_0_0",
      "x_0_1": "x_0_1",
      "y_1": "y_0",
      "x_1_0": "x_0_0",
      "x_1_1": "x_0_1",
      "y_2": "y_0",
      "x_2_0": "x_0_0",
      "x_2_1": "x_0_1",
      "y_3": "y_0",
      "x_3_0": "x_0_0",
      "x_3_1": "x_0_1"
    }
  }
}
```

## .jdec file structure
The root object of a `.jdec` file contains meta information about the decomposition and the decomposition itself.
The following data fields are supported:
| Name                      | Mandatory | Type      | Default | Description |
|:--------------------------|:---------:|:----------|:--------|:----|
| `version`                 | X         | integer   |         | version of the file format |
| `problem_name`            |           | string    |         | name of the problem the decomposition belongs to |
| `description`             |           | string    |         | description of the decomposition |
| `decomposition`           | X         | object    |         | the decomposition (see below) |
| `decomposition_id`        |           | integer   |         | internal ID of the decomposition, only written by GCG |

A decomposition object supports the following data fields:
| Name                      | Mandatory | Type      | Default | Description |
|:--------------------------|:---------:|:----------|:--------|:----|
| `n_blocks`                |           | integer   |         | number of blocks of the decomposition, ignored by GCG |
| `presolved`               |           | boolean   | `false` | indicates whether the decomposition refers to a presolved problem (must be `false` for a decomposition of a block) |
| `master_constraints`      | X         | list      |         | list of master constraints |
| `blocks`                  | X         | list      |         | list of block objects (see below) |
| `symmetry_var_mapping`    |           | object    |         | a mapping that maps all variables to their representatives |

A block object supports the following data fields:
| Name                      | Mandatory | Type      | Default | Description |
|:--------------------------|:---------:|:----------|:--------|:----|
| `index`                   |           | integer   |         | index (starting at 0) of the block (note: if set, it has to be set for all blocks of a decomposition) |
| `constraints`             | X         | list      |         | list of constraints assigned to the block |
| `decomposition`           |           | object    |         | decomposition (object) of the block |
| `symmetry_representative_block`|      | integer   |         | index of the representative block this block should be mapped to |

> Note: The provided symmetry information has to be consistent. GCG does not validate the mappings in detail.
> Representative blocks must be the blocks with the lowest indices of the blocks that should be aggregated (and are symmetrical).
> Variable mappings have to map the block variables of the blocks to the corresponding representative variable of the representative block.
> The mappings of the representative blocks/variables to themself can be omitted.
