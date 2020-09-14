# Pricing Parameters {#pricing-params}
# Modifying Parameters of the GCG Pricing
We have around 2000 parameters for the pricing available, so it might be difficult
to choose the right ones for your class of problems. In the following, some hints
about good settings will be given.

### Parameters often having an impact
@todo @schilling

### Change other Parameters of the Pricing
You can also modify other parameters of the pricing. 
You can set them by entering `set` and then `pricing`. 
In the following, we give a list of all parameters that can then be changed.

# List of Pricing Parameters
A list of detection parameters can also be seen (and searched) when inserting `set pricing` or `set pricingsolver` into the search box on the page of @ref interactive-menu.

### Pricer
```
  abortfac              pricing is aborted, if fac * pricing/maxvars pricing candidates were found [2.0]
  delvars               should variables created at the current node be deleted when the node is solved in case they are not present in the LP anymore? [FALSE]
  delvarsroot           should variables created at the root node be deleted when the root is solved in case they are not present in the LP anymore? [FALSE]
  maxvars               maximal number of variables priced in per pricing round [100]
  maxvarsroot           maximal number of priced variables at the root node [2000]
```

### Masterpricer
```
  <advanced>            advanced parameters
    abortpricinggap       gap between dual bound and RMP objective at which pricing is aborted [0.0]
    abortpricingint       should pricing be aborted due to integral objective function? [TRUE]
    chunksize             maximal number of pricing problems to be solved during one pricing loop [2147483647]
    nroundscol            number of previous pricing rounds for which the number of improving columns should be counted [15]

  <colpool>             parameters for <colpool>
    agelimit              age limit for columns in column pool? (-1 for no limit) [100]

  <pricestore>          parameters for <pricestore>
    efficiacychoice       choice to base efficiacy on [0]
    mincolorth            minimal orthogonality of columns to add [0.0]
    objparalfac           factor of objective parallelism in score function [0.0]
    orthofac              factor of orthogonalities in score function [0.0]
    redcostfac            factor of -redcost/norm in score function [1.0]

  bigmartificial        value for for big M objective of artificial variables (negative if max obj should be used) [1000.0]
  disablecutoff         should the cutoffbound be applied in master LP solving (0: on, 1:off, 2:auto)? [2]
  dispinfos             should additional informations concerning the pricing process be displayed? [FALSE]
  eagerfreq             frequency at which all pricingproblems should be solved (0 to disable) [10]
  factorunreliable      factor to use for objective of unbounded variables [1000.0]
  heurpricingiters      maximum number of heuristic pricing iterations per pricing call and problem [1]
  maxcolsprobfarkas     maximal number of columns per problem to be generated during Farkas pricing [10]
  maxcolsprobredcost    maximal number of columns per problem to be generated during red. cost pricing [10]
  maxcolsprobredcostroot 
                   -->  maximal number of columns per problem to be generated during red. cost pricing at root node [10]
  maxcolsroundfarkas    maximal number of columns per Farkas pricing round [10]
  maxcolsroundredcost   maximal number of columns per reduced cost pricing round [100]
  maxcolsroundredcostroot 
                   -->  maximal number of columns per reduced cost pricing round at root node [100]
  maxheurdepth          maximum depth at which heuristic pricing should be performed (-1 for infinity) [-1]
  maxroundsredcost      maximal number of pricing rounds per node after the root node [2147483647]
  maxsuccessfulprobsredcost 
                   -->  maximal number of successfully solved red. cost pricing problems until pricing loop is aborted [2147483647]
  onlyreliablebigm      only use maxobj for big M objective of artificial variables if it is reliable [TRUE]
  relmaxprobsfarkas     maximal percentage of Farkas pricing problems that are solved if variables have already been found [1.0]
  relmaxprobsredcost    maximal percentage of red. cost pricing problems that are solved if variables have already been found [1.0]
  relmaxprobsredcostroot 
                   -->  maximal percentage of red. cost pricing problems that are solved at root node if variables have already been found [1.0]
  relmaxsuccessfulprobsredcost 
                   -->  maximal percentage of successfully solved red. cost pricing problems until pricing loop is aborted [1.0]
  sorting               order by which the pricing problems should be sorted ('i'ndices, 'd'ual solutions of convexity constraints, 'r'eliability from previous rounds, reliability from the 'l'ast nroundscol rounds) [r]
  stabilization         should stabilization be performed? [TRUE]
  stabilizationtree     should stabilization be performed in the tree (in nodes other than the root node)? [FALSE]
  threads               how many threads should be used to concurrently solve the pricing problem (0 to guess threads by OpenMP) [0]
  useartificialvars     should artificial variables be used to make the RMP feasible (instead of applying Farkas pricing)? [FALSE]
  usecolpool            should the colpool be checked for negative redcost cols before solving the pricing problems? [TRUE]
  usemaxobj             use maxobj for big M objective of artificial variables [TRUE]
```

### Pricing Solvers
```
  <knapsack>            parameters for <knapsack>
    exactenabled          flag to indicate whether exact solving method of solver <knapsack> is enabled [TRUE]
    heurenabled           flag to indicate whether heuristic solving method of solver <knapsack> is enabled [FALSE]
    priority              priority of solver <knapsack> [200]

  <mip>                 parameters for <mip>
    <advanced>            advanced parameters
        checksols             should solutions of the pricing MIPs be checked for duplicity? [TRUE]
        gaplimitfac           factor by which to decrease gap limit for heuristic pricing (1.0: subtract start limit) [0.8]
        nodelimitfac          factor by which to increase node limit for heuristic pricing (1.0: add start limit) [1.0]
        settingsfile          settings file for pricing problems [-]
        sollimitfac           factor by which to increase solution limit for heuristic pricing (1.0: add start limit) [1.0]
        stallnodelimitfac     factor by which to increase stalling node limit for heuristic pricing (1.0: add start limit) [1.0]
        startgaplimit         start gap limit for heuristic pricing [0.2]
        startnodelimit        start node limit for heuristic pricing [1000]
        startsollimit         start solution limit for heuristic pricing [10]
        startstallnodelimit   start stalling node limit for heuristic pricing [100]
    exactenabled          flag to indicate whether exact solving method of solver <mip> is enabled [TRUE]
    heurenabled           flag to indicate whether heuristic solving method of solver <mip> is enabled [TRUE]
    priority              priority of solver <mip> [0]
```