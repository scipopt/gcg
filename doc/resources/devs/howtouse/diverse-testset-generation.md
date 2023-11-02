# How to generate diverse test sets {#diverse-test-sets}
> With this guide we describe how one can create (diverse) test sets from
> a large number of instances and decompositions to carry out experiments on.

# Curate Diverse Test Sets
As presented by Donkiewicz (2021) with the help of the strIPlib,
it is possible to generate diverse test sets. For theoretical details,
we refer the reader to the thesis.

## Preparations
### Basis Test Set
First, you will need a test set that the script will use as ground set to filter from.
This can be done e.g., by downloading a filtered test set from the [strIPlib page](https://striplib.or.rwth-aachen.de).
Due to the script's runtime, we recommend you not to exceed 1000 instances with the ground set.
We will consider the filtered ground test set as `filtered.test` in the following.

### Required Python Packages
We require python3.6 upwards and the following packages to be installed:
```
pandas scipy sklearn pyscipopt
```
For a tutorial on how to install PySCIPOpt, we refer the reader to the [respective installation page](https://pypi.org/project/PySCIPOpt/).

## Experiment Test Sets Generation
If you are a GCG developer and call the script inside `/shared/gcg_runtime_data/diverse_testset_generation`,
the following command suffices:

    ./striplib.py --hours 30 -t filtered.test

Otherwise, you simply call the `striplib.py` script with your static and runtime data files
(either `.pkl` or `.csv`). A command would then look like that:

    ./striplib.py --experiment --data results/static.csv --runtimes results/runtime_data.pkl --hours 30 -t filtered.test

where the runtimes pickle (here, `striplib.pkl`) must be in the format as exported
by the general plotter (see @ref visu-suite "visualization suite documentation" for more information).
The `hours`-flag defines how long the experiment should take (approximated, based
on the runtime data you provided with the runtimes pickle).
The given test set file defines the ground set of instances to select a diverse
set from. If no decompositions are given in it (`instance.lp` instead of `instance.lp;instance.dec`),
it is assumed that GCG should use these instances, detecting without any given decomposition.

Additional parameters are available to tune the selection process:
```
--kappa_min KAPPA_MIN
                      [Experiment] Minimum number of instances chosen from each suitable model group.
--kappa_max KAPPA_MAX
                      [Experiment] Maximum number of instances chosen from each suitable model group.
--mininstances MININSTANCES
                      [Experiment] Minimum number of instances chosen.
--maxinstances MAXINSTANCES
                      [Experiment] Maximum number of instances chosen.
-cs CLUSTER_STRICTNESS, --cluster_strictness CLUSTER_STRICTNESS
                      [Experiment] Diversity strictness ([0,1]).
```

The script will execute a modified variant of the strIPlib benchmark selection
and finally export a test set (including decompositions).
