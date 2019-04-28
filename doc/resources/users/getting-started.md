# How to get started with GCG {#getting-started}
## Start GCG
To start GCG, you just type `gcg` into your console, preferably inside the GCG
root folder. If the console prompts you to install packages or throws an error,
that such command does not exist, you have not installed everything correctly.
Otherwise, you're in and you should see something like:
```
GCG version 3.0.1 [GitHash: dbc6c61]
Copyright (C) 2010-2019 Operations Research, RWTH Aachen University
                        Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)

...

GCG>

```

You can now enter commands into the interactive console.

## Basic commands
### `read`
GCG reads your LP file. If you go to `/check/instances/bpp` from your GCG root
directory, you can read the problem "N1C3W1_A". Output:
```
GCG> read
enter filename: N1C3W1_A.lp

read problem <N1C3W1_A.lp>
============

original problem has 1224 variables (0 bin, 1224 int, 0 impl, 0 cont) and 74 constraints

```

### `optimize`
GCG optimizes what you previously read.
```
GCG> optimize

...

A Dantzig-Wolfe reformulation is applied to solve the original problem.

  time | node  | left  |SLP iter|MLP iter|LP it/n| mem |mdpt |mvars|ocons|mcons|mcuts|sepa|  dualbound   | primalbound  |  deg   |  gap   
T  0.1s|     1 |     0 |      0 |      0 |     - |6428k|   0 |   0 |  74 |   0 |   0 |  0 | 0.000000e+00 | 2.400000e+01 |   --   |    Inf
b  0.1s|     1 |     0 |      0 |      0 |     - |6279k|   0 |   0 |  74 |   0 |   0 |  0 | 0.000000e+00 | 2.200000e+01 |   --   |    Inf
   0.1s|     1 |     0 |      0 |      0 |     - |6262k|   0 |   0 |  74 |   0 |   0 |  0 | 1.587333e+01 | 2.200000e+01 |   --   |  38.60%
Chosen structure has 24 blocks and 50 linking constraints.
This decomposition has a maxwhite score of 0.310811.
Master problem is a set covering problem.
Matrix has 24 blocks, using 1 aggregated pricing problem.


   0.2s|     1 |     0 |      0 |      0 |     - |  14M|   0 |  24 |  75 |  52 |   0 |  0 | 1.587333e+01 | 2.200000e+01 |   0.00%|  38.60%
   0.2s|     1 |     0 |      0 |      0 |     - |  14M|   0 |  48 |  75 |  52 |   0 |  0 | 1.587333e+01 | 2.200000e+01 |   0.00%|  38.60%
Starting reduced cost pricing...
*r 0.2s|     1 |     0 |    515 |    515 |     - |  15M|   0 | 122 |  75 |  52 |   0 |  0 | 1.587333e+01 | 1.600000e+01 |  12.13%|   0.80%
   0.2s|     1 |     0 |    515 |    515 |     - |  15M|   0 | 122 |  75 |  52 |   0 |  0 | 1.587333e+01 | 1.600000e+01 |  12.13%|   0.80%
   0.2s|     1 |     0 |    516 |    516 |     - |  15M|   0 | 122 |  75 |  52 |   0 |  0 | 1.600000e+01 | 1.600000e+01 |  12.13%|   0.00%
   0.2s|     1 |     0 |    516 |    516 |     - |  15M|   0 | 122 |  75 |  52 |   0 |  0 | 1.600000e+01 | 1.600000e+01 |   --   |   0.00%

SCIP Status        : problem is solved [optimal solution found]
Solving Time (sec) : 0.19
Solving Nodes      : 1
Primal Bound       : +1.60000000000000e+01 (3 solutions)
Dual Bound         : +1.60000000000000e+01
Gap                : 0.00 %


```

#### `detect` (included in `optimize`)
Detects structure and applies Dantzig-Wolfe Reformulation.

#### `presolve` (included in `optimize`)
Presolves your LP.

### `free`
Removes current LP from memory. You can then read another.

### `quit`
Ends GCG.

## Advanced Commands
### `validatesolve`
Checks your solution against a `.solu` file.

### `explore`
The explore menu offers you insight into the different decompositions created
in the `detect`-step.

### `write`
Can write out multiple things, e.g. decompositions.
