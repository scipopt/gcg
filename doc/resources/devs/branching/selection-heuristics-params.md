# Parameters for Strong Branching-Based Branching Candidate Selection Heuristics {#selection-heuristics-params}

See @ref branching-selection-heuristics "here" for explanations regarding branching candidate selection heuristics.
## Activating selection heuristics
Per default, **strong branching is disabled**. Furthermore, GCG will use the first active non-strong branching heuristics in place of 
pseudocost branching for hybrid, reliability and hierarchical branching.

### Original Variable Branching
```
branching/orig/usepseudocosts =  TRUE
branching/orig/mostfrac       = FALSE
branching/orig/random         = FALSE
branching/orig/usestrong      = FALSE
```

### Ryan-Foster Branching
```
branching/bp_strong/ryanfoster/usepseudocosts = FALSE
branching/bp_strong/ryanfoster/mostfrac       = FALSE
# (else random branching)
branching/ryanfoster/usestrong                =  TRUE
```


## Branching rule-dependent settings for hierarchical strong branching

### Number of candidates in phase 0
```
# minimum number of output candidates from phase 0 [1..MAX]
branching/orig/minphase0outcands           =  10
branching/ryanfoster/minphase0outcands     =  10

# maximum number of output candidates from phase 0 [1..MAX]
branching/orig/maxphase0outcands           =  50
branching/ryanfoster/maxphase0outcands     =  50

# maximum number of output candidates from phase 0 as fraction of total candidates [0,1]
branching/orig/maxphase0outcandsfrac       = 0.7
branching/ryanfoster/maxphase0outcandsfrac = 0.7

# how much impact should the node gap have on the number of precisely evaluated
# candidates in phase 1 during strong branching? [0,1]
branching/orig/phase1gapweight             = 0.25
branching/ryanfoster/phase1gapweight       = 0.25
```

### Number of candidates in phase 1
```
# minimum number of output candidates from phase 1 [1..MAX]
branching/orig/minphase1outcands           =   3
branching/ryanfoster/minphase1outcands     =   3

# maximum number of output candidates from phase 1 [1..MAX]
branching/orig/maxphase1outcands           =  20
branching/ryanfoster/maxphase1outcands     =  20

# maximum number of output candidates from phase 1 as fraction of total
# candidates [0,1]
branching/orig/maxphase1outcandsfrac       = 0.7
branching/ryanfoster/maxphase1outcandsfrac = 0.7

# how much impact should the node gap have on the number of precisely evaluated
# candidates in phase 2 during strong branching? [0,1]
branching/orig/phase2gapweight             =   1
branching/ryanfoster/phase2gapweight       =   1
```

## Additional hierarchical strong branching parameters
```
# minimum number of variables for phase 2 to be executed, otherwise the best
# candidate from phase 1 will be chosen [0..MAX]
branching/bp_strong/mincolgencands =   4

# how many candidates should be chosen based on historical strong branching
# scores as opposed to current heuristic scores in phase 0 (e.g. 0.5 = 50%)? [0,1]
branching/bp_strong/histweight     = 0.5
```

## Hybrid Branching
```
# minimum tree depth from which on phase 0 is performed [0..MAX]
branching/bp_strong/minphase0depth     =    0

# maximum tree depth up to which phase 1 is performed [0..MAX]
branching/bp_strong/maxphase1depth     =    4

# maximum tree depth up to which phase 2 is performed [0..MAX]
branching/bp_strong/maxphase2depth     =    3

# how much should the logarithm of the number of variables influence the depth
# for hybrid branching? (0 = not at all, 1 = fully) [0,1]
branching/bp_strong/depthlogweight     = 0.5

# what should be the base of the logarithm that is used to compute the depth of
# hybrid branching? [0,MAX]
branching/bp_strong/depthlogbase       = 3.5

# if using a logarithm to compute the depth of hybrid branching, what should be
# the fraction of the depth assigned to phase 1 that is assigned to phase 0? [0,1]
branching/bp_strong/depthlogphase0frac =    0

# if using a logarithm to compute the depth of hybrid branching, what should be
# the fraction of the depth assigned to phase 1 that is assigned to phase 2? [0,1]
branching/bp_strong/depthlogphase2frac = 0.75
```

## Reliability Branching
```
# min count of pseudocost scores for a variable to be considered reliable in
# phase 1 (-1 = never) [-1..MAX]
branching/bp_strong/phase1reliable = 2147483000

# min count of pseudocost scores for a variable to be considered reliable in
# phase 2 (-1 = never) [-1..MAX]
branching/bp_strong/phase2reliable = 2147483000
```

## Additional strong branching parameters
```
# should infeasibility detected during strong branching be handled immediately,
# or only if the candidate is selected? [TRUE/FALSE]
branching/bp_strong/immediateinf       =  TRUE

# how many times can bounds be changed due to infeasibility during strong
# branching until an already evaluated variable needs to be reevaluated? [0..MAX]
branching/bp_strong/reevalage          =     1

# maximum number of strong branching lp iterations,
# set to 2*avg lp iterations if <= 0 [0..MAX]
branching/bp_strong/maxsblpiters       = 2147483000

# maximum number of strong branching price rounds,
# set to 2*avg lp iterations if <= 0 [0..MAX]
branching/bp_strong/maxsbpricerounds   = 2147483000

# maximum number of non-improving candidates until phase 2 is stopped [0..MAX]
branching/bp_strong/maxlookahead       =     8

# how much should the lookahead scale with the overall evaluation effort?
# (0 = not at all, 1 = fully) [0,1]
branching/bp_strong/lookaheadscales    =   0.5

# what percentage of the strong branching score of the candidate that was selected
# does the heuristic's incumbent need to be considered close (e.g. 0.5 = 50%)? [0,1]
branching/bp_strong/closepercentage    =   0.9

# how many times in a row can the heuristic be close before strong branching is
# stopped? (-1 = forever) [-1..MAX]
branching/bp_strong/maxconsecheurclose =     4

# with how much weight should strong branching scores be considered for
# pseudocost scores? [0,1]
branching/bp_strong/sbpseudocostweight =     1

# should phase 0 be performed even if the number of input candidates is already
# lower or equal to the number of output candidates? [TRUE/FALSE]
branching/bp_strong/forcep0            = FALSE
```