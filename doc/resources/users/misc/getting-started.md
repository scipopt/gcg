# Getting started {#getting-started}
## Start GCG
If you did a "global" install with `make install` of the SCIP Optimization Suite,
just type `gcg` into your console. \n
Otherwise, execute `./bin/gcg` from the GCG root folder.
You should then see something like:
```
GCG version 3.0.1 [GitHash: dbc6c61]
Copyright (C) 2010-2019 Operations Research, RWTH Aachen University
                        Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)

...

GCG>

```

You can now enter commands into the interactive console.

## Basic commands
> In the interactive console, you can always leave a submenu with `..`. \n
> If you want to quit GCG, just type `quit`.

Now let's go ahead and read some problem:

### Make GCG know your problem with `read`
With this command, GCG reads your LP file. If you, for example, go to `/check/instances/bpp` from your GCG root
directory, you can read the problem "N1C3W1_A.lp". The output looks as following:
```
GCG> read
enter filename: N1C3W1_A.lp

read problem <N1C3W1_A.lp>
============

original problem has 1224 variables (0 bin, 1224 int, 0 impl, 0 cont) and 74 constraints

```

You can obviously also read your own problem, as long as it is in a format as described in @subpage input-formats.

### Solve your problem with `optimize`
Here the magic happens: GCG optimizes what you previously read.
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

### Look at your decompositions with `explore`
This menu is probably the most interesting for you, thus it earned itself its own page: @ref explore-menu.

### Start from new with `free`
If you now want to start all over and read a different problem, that's no problem. `free` will happily
remove the current LP from memory, such that you can then read another.


## Advanced Commands
After these basic commands, you might want to change some more settings or give GCG information about your problem.

### GCG decompositions
**Detecting without solving:** This is possible with a simple `detect`. The instance already has to be read.\n
**Marking your own .dec file:** \n
**Reading your decomposition:** You can read decomposition files (with formats as given @ref input-formats "here") with a simple `read`.
The behavior when giving own decomposition files is described @ref detectionloop "here".\n
**Saving your decomposition:** You can save all decompositions found with `write alldecompositions`. To just save single ones, you have
to select a decomposition in the @ref explore-menu "Explore Menu" with `explore select <nr>` and then `write selecteddecompositions`.

### GCG settings
Everything related to settings can be found in the `set` submenu. All possible settings that can be done there, can be found @ref params "here".\n
**Making your own settings file:** If you want to create a whole settings file from the scratch, you can find information
about that <a href="FAQ.html#createsettingsfile">here</a>.\n
**Saving your settings:** Otherwise, you can just set everything you want (e.g. `set display 1` to make GCG display solving information
  every second) and then save those settings with `save`. You will
be prompted to give a settings name (which can also be a path), that should end
with `.set`.\n
**Loading your settings:** You can also load your own settings file with `load`. You have to give the path to the settings file (relative to the folder that
you are in), with the `.set` extension of the file. If some settings could not
be applied, you will be reminded of that.
