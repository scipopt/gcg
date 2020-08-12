# Modifying Parameters of the GCG Presolving {#presolving-params}
There are different options to influence the presolving done to your problem. You can either
   * **Deactivate** presolving for your instance completely
   * **Change parameters** for the presolving

### Deactivate Presolving
The easiest way to not use presolving is to simply read in your instance and then, instead of doing an
`optimize`, perform `detect` and `optimize` instead. We will work on the unpresolved instance with its
respective decomposition. You will be able to read the following message in the log: 

```
GCG> optimize
there is an original decomposition and problem is not presolved yet -> disable presolving and start optimizing (rerun with presolve command before detect command for detecting in presolved problem)
```

### Change Parameters of the Presolving
You can modify different parameters of the presolving. 
You can set them by entering `set` and then `presolving`. 
In the following, we give a list of all parameters that can then be changed.

```
  <advanced>            advanced parameters
  <boundshift>          converts variables with domain [a,b] to variables with domain [0,b-a]
  <convertinttobin>     converts integer variables to binaries
  <domcol>              dominated column presolver
  <emphasis>            predefined parameter settings
  <implics>             implication graph aggregator
  <inttobinary>         converts integer variables with domain [a,a+1] to binaries
  <trivial>             round fractional bounds on integers, fix variables with equal bounds
  maxrestarts           maximal number of restarts (-1: unlimited) [0]
  maxrounds             maximal number of presolving rounds (-1: unlimited, 0: off) [0]
```