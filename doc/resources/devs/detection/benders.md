# Benders Detection Mode {#benders}

# Using GCG's Benders Mode

## Introduction into Benders Decomposition
We are offering a Benders decomposition mode (also called variable/Benders partitioning or "L-shaped method" 
in two-stage stochastic programming), which is said to be the "dual version" of Dantzig-Wolfe and Lagrange.
With a Benders decomposition, we are generating **cuts instead of variables** (also called row generation)
and thus we also have a (mixed-) integer master problem, as opposed to a relaxed (continuous) one, as well
as a linear subproblem. You can find more information on the theory in the paper by Geoffrion (1972) \cite geoffrion1972.\n
It is particularly sensible to apply Benders decomposition when there are **very small links between subproblems**,
for example if you have a production planning program that should be linked to a vehicle routing problem for delivering
the produced goods.

## Activating GCG's Benders functionality automatically
In order to activate the full Benders capabilities, you first have to switch the mode. After starting GCG,
perform the following steps.

#### 1. Activating Benders Relaxation
Switch to the relaxing mode "Benders", which is number 1.
```
GCG> set relaxing gcg mode
current value: 0, new value [0,2]: 1
relaxing/gcg/mode = 1

```

#### 2. Loading a settings file
Then, a `.set` file should be loaded by executing
```
set load settings/SETNAME
```

Depending on what you know about your problem, there are different predefined options for the Benders Mode:

- Standard Mode
  - Default Mode `detect-benders.set`
  - Fast Mode `detect-benders-fast.set`
- Binary Master
  - Default Mode `detect-benders-bin_master.set`
  - Fast Mode `detect-benders-bin_master-fast.set`
- Continuous Subproblems
  - Default Mode `detect-benders-cont_subpr.set`
  - Fast Mode `detect-benders-cont_subpr-fast.set`

#### 3. Solve your problem!
If you want to see some more statistics in the statistics table, you can activate additional Benders
information. 

```
GCG> set table benders active
current value: TRUE, new value (TRUE/FALSE): true
table/benders/active = TRUE
```

## Activating GCG's Benders functionality manually
#### 1. Activate Benders Relaxation
```
GCG> set relaxing gcg mode
current value: 0, new value [0,2]: 1
relaxing/gcg/mode = 1
```

#### 2. Activate Benders Detection
```
GCG> set detection benders enabled
current value: FALSE, new value (TRUE/FALSE): true
detection/benders/enabled = TRUE
```

#### 3. Set a different, more suitable decomposition score (e.g. Benders experimental score)
In order to find better suitable decompositions, you should use a different score than the default one.
```
GCG> set detection score scoretype 7
detection/score/scoretype = 7
```

#### 4. Modify the Detection
You should take care of well-chosen detection settings. The settings files above already feature them.
They e.g. include setting the `maxrounds` parameter.

#### 5. Deactivate/Activate Detectors
Along with general detection parameters, you should also disable and enable detectors depending on
your problem. Again, the settings files above already feature these settings.

#### 6. [optional] Show additional statistics for Benders in the log
```
GCG> set table benders active
current value: TRUE, new value (TRUE/FALSE): true
table/benders/active = TRUE
```