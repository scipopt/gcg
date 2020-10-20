# How to conduct experiments with GCG {#conduct-experiments}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

[TOC]

> When trying to find out if own modifications make a difference on solver behavior,
> especially on runtimes, one has to do let GCG run for both and compare the two (or more). 
> On this page, we explain **best practices of and give a guide for conducting experiments** 
> with GCG.   

# Generating Runtime Data {#generate-data}
The first step to start experimenting is to collect runtime data that the solver outputs.
This is in particular the log (explained in @ref u1), as well as `vbc` files 
("Visualization of Branch Cut algorithms").

## 1. Checking your installation
Before starting with your (potentially quite big) testset, we recommend 
**checking that your installation and configuration of GCG is correct** by performing a
test on the default (called `short`) testset. This can be done by executing

```
make test
```

from inside the GCG root directory.

If it went fine, the last lines of your output should look somewhat like this:

```
------------------------------[Nodes]---------------[Time]----------[Pricing-Time]--------[LP-Time]-------[Pricing-Calls]--
Cnt  Pass  Time  Fail  total(k)     geom.     total     geom.     total     geom.     total     geom.    total     geom. 
---------------------------------------------------------------------------------------------------------------------------
14    14     0     0         0       5.6      12.9       1.3       8.4       1.2       0.7       1.0      1229      26.8 
shifted geom. [ 100/10.0/5.0/2.0]     18.3                 0.9                 5.5                 2.0                 0.0 
---------------------------------------------------------------------------------------------------------------------------
@02 timelimit: 3600
@01 GCG(3.1.0)SCIP(7.0.0)spx(5.0.0):default
```

If it does, go ahead and continue with the next step.

## 2. Defining a test set
If you have a problem file, e.g. an `lp` file, this is considered to be an "instance"
of your problem. Multiple instances assemble a "test set". In a test set file, ending
with `.test`, each line gives a path to an instance file. For checking plausibility 
and feasibility of solutions, there exist `solu` files that specify each instance's solution.\n

It is crucial that you **determine with which instances your modifications of settings
or code can be seen**, e.g. in the total runtime. This could be, for example, instances that are
solved in the root node, or those that exhibit a staircase structure etc. Depending on
what you are filtering for, there are different options one usually takes:

-   In some cases, you have the objective to solve a given set of instances quicker. In
    this case, you can just **assemble the set using these exact instances**, e.g. because 
    you already have them lieing around. However, to show that your changes do not influence 
    other instances negatively, you should also test using a 
    @ref general-purpose-testsets "general-purpose testset".
-   In most cases, you will know some properties (the most common property being a 
    sufficiently low, but not too low runtime) that the instances have to fulfill, 
    but don't have instances at hand. For this case, using our 
    **[strIPlib](https://striplib.or.rwth-aachen.de)** (structured integer programming library), 
    you can **filter for a range of different characteristics** to get a set of instances.
    If you don't want to use the web frontend, please refer to our guide on @ref testset-selection.
-   In some cases, you do not know anything about how the changes might affect GCG's
    behavior. For these, we have some (partly very specific) **predefined testset files** 
    inside the folder `/check/testset` (in particular, the @ref general-purpose-testsets "general-purpose testsets"
    might be interesting for you). They mostly require a link to the strIPlib, 
    so you will have to download the respective instances or, as a GCG developer, 
    link the strIPlib directory (see [Git Wiki](https://git.or.rwth-aachen.de/gcg/gcg/-/wikis/GCG-Wiki/Testing#setting-links-if-required) for more information).

### General Purpose Testsets {#general-purpose-testsets}
In every case, test sets should **cover a variety of characteristics** in order
to correctly represent the changes done. To have a diverse set that gives
an **impression on improvements and deteriorations**, we have generated some sets
that might be suitable.

> Currently, the general purpose test sets are not very "reliable", i.e. diverse
> and still running short enough. We are working on a solution.

In the following, you will find a list of testsets with descriptions of their
intended application:
* `goodinstances`: instances with at most 1 hour of runtime that cover most of
the algorithmic machinery

## 3. Perform automated tests {#testing}
> Just like SCIP, GCG supports automatic testing. It uses the same command
> structure as SCIP, so if any information is missing, you can also consult 
> the <a href="https://scipopt.org/doc/html/TEST.php">SCIP documentation</a>.

Finally, we want to **perform the test** with your testset `mytest.test`, located in the 
folder `check/testset` and your customized settings `modified.set`, located in the
folder `settings`. This can be done easily with 

```
make TEST=mytest SETTINGS=modified test
```

As you can see, you must not give the file endings for the testset and settings file.
If you just want to read into what exactly GCG is doing, this should already be the
end of your journey. Otherwise, you might need to set some flags to work with the data
some more.

### Compilation and Runtime Flags
There exist more **arguments and flags** apart from the testset, statistics and 
setting **that can be configured and set**.
All possible arguments for automatic testing in GCG can be found 
@ref makefiles-args "here (makefiles)" or @ref cmake-args "here (cmake)". \n
Some information about important flags when testing: 

-   If you plan on **visualizing** the results, you should _set_ the
    variable `STATISTICS=true`. This will lead to, for each instance, GCG printing
    out the additional statistics that you can also get when doing `display additionalstatistics`.
-   If you want to **analyze algorithmic behavior** very deeply, you
    should _compile_ GCG and SCIP with `STATISTICS=true`. This will
    print more output that is flagged accordingly inside the code.
    For some visualization scripts, this is required (e.g. pricing plotter).

Note that the more lines are printed, the bigger the logs get and, for 
a significant subset of instances, this will definitely slow down the
whole solving process, since printing, if it is very often, will use
computation time.


## 4. Checking generated output files {#what-files}
After executing the test, GCG will **automatically have exported
some files**. After navigating to the folder `check/results`,
you will see multiple files per test. The output files for a test 
look somewhat similar to this format, giving you meta information about
the test run:
```
check.TESTSETNAME.gcg-GCGVERSION.linux.x86_64.gnu.opt.spx2.COMPUTERNAME.SETTINGSNAME.default{.err,.out,.set}
```

The three files `{.err,.out,.set}` are always generated during runtime. 
Additional files with endings `{.res,.pav,.tex}` are generated after the test was 
completed. However, if you aborted the test before it finished, you can still 
generate those files (possibly incompletely) using the command `./evalcheck.sh` 
in the `/check` folder. Note: changing the files' names is not recommended, since 
some visualization scripts might try to get their information from the file name.

#### `out`-files
The .out files contain the `stdout`-output and therefore all data that was
printed during the execution of GCG. From them, data can be parsed more or less
easily. In them, you can read details about the run (what details exactly depend
on the flags you set (see "Compilation and Runtime Flags").<br>
#### `res`-and `tex`-files
The .res files only contain an overview table with primal bound, dual bound, gap,
used heuristics more information. They do not show any temporal runtime data. Most
commonly, they can be used for checking the status (whether a solution was proven to
be correct using a `solu` file) and the time GCG took to solve this instance.<br>
The `res` file is also given as embeddable LaTeX code in the `tex` file.<br>
#### `vbc`-files
These files show how the branch-and-bound tree developed. They can be played
back using [vbctool](https://cs.uni-koeln.de/ls-juenger/software/vbctool),
allowing us to see how the tree developed and also how e.g. node selection
heuristics worked. `vbc` files are only generated when the test is conducted with
`STATISTICS=true`.

> For comparing different runs, we will only need `out` and `vbc` files. 
> The other files are handy for "manual" comparisons, but will not be required further
> in this guide.

# Evaluating a Single Test Set
> In this part, we explain how to generate a report for a test set, not comparing
> anything, but just **exploring your instances**.

After collecting the data, you can start analyzing and in particular **visualizing it**.
We assume that you have collected your runtime data using e.g. `make test`. To then
generate a **report for the whole testset** including basic statistics and all
possible visualizations, you can call `make visu`. More details on the test set
report feature are explained @ref testset-report "here".
The test set report can be found in `check/reports/` in a uniquely named (using a time stamp) report folder.

# Comparing Runtime Data
If you do not only have one test run, but multiple ones, you will want to compare it.
In the following, we will give some hints on how to start, depending on what you want to compare.
The most common cases are comparing default GCG...

A. against a GCG with different **settings** (continue @ref compare-settings "here") \n
B. against an **algorithmically** modified GCG version (continue @ref compare-code "here") \n
C. against a **different version** of default GCG (continue @ref compare-versions "here") \n
D. against a **different solver** (continue @ref compare-solvers "here") \n


Depending on what you want to test, skip to the respective section.

In general, please again make sure that you compiled your GCG correctly and executed
the test with `STATISTICS=true`. Furthermore, have your runtime data
(`.out`, `.res` and `vbc/` files) ready in a known directory.

## A. Comparing Settings {#compare-settings}
> In this part, we explain how to **compare different settings** (with unmodified code)
> using our automatic comparison report script.

If you want to **compare settings** with otherwise unmodified (and similarly versioned)
installations of GCG, then you have already collected runtime data, e.g. by two calls,
`make test SET=settings1` and `make test SET=settings2` and have your runtime data
(2 or more sets of `.out`, `.res` and `vbc/` files) ready.\n
To then generate a **comparison report for the two runs** including basic statistics
and all possible _comparing_ statistics, you can call

    make compare

More details on the comparison report feature are explained @ref comparison-report "here".
The version comparison report can be found in `check/reports/` in a uniquely named (using a time stamp) report folder.

## B. Comparing Code {#compare-code}
> In this part, we explain how to parse own expressions that your custom plug-ins and other
> code might have generated that you want to visualize.

If you want to **compare different code** with similarly versioned installations of GCG,
then you have already collected runtime data, e.g. by **a `make test` for the default GCG**
and **a `make test` for your custom version of GCG** and have your runtime data
(2 or more sets of `.out`, `.res` and `vbc/` files) ready.\n
For our visualization scripts to keep on working, you may not change any existing
GCG logging formats to not get differently formatted `.out`-files that will not be
parseable anymore by the parsers. Instead, you should print output that you later
want to visualize in the ways presented in the following:

### Logging Custom Data: Single Data Points
> **The custom logging feature (single data points) is not yet implemented. Please stay tuned.**
> This is a sketch of what it will be capable of.

If you want to log just single events, like the overall runtime inside your code,
which you just have to print once per instance, we recommend that you print lines
in the following format:

    codeIdentifier id1:val1 id2:val2 id3:val3

Here, `codeIdentifier` must be a unique line start, identifying your code. Commonly,
this is some short form of your project's name, e.g. `subGCG` or `[subGCG]`.
The general parser will search for this identifier and **save all data points that come
after it**. To save them correctly, you have to give an (again unique) data point identifier,
e.g. `totaltime`, at the place of `id1` and then print the value `val1`, e.g. `42`.
The same holds for the rest of the ID-value pairs. In the dataframe (where for each
instance in the test run there is a line), a column with the name `id1` will then
be introduced for all instances.\n
In order to make the identifiers known to the parser, you have to give GCG the
`codeIdentifier` (the other identifiers are always in the line parsed, no need to give them)
using the flag `-ci [CID]` in one of the four general plotters. A sample call could then be

    python3 general/plot.py check.gcg.out -ci "subGCG" -t "PRICING TIME" 10 20 "totaltime"

As you can see, you can use the data point identifiers directly (caps sensitive) to plot.
In order to automatically generate comparison reports with your own printed data,
you can give this setting as
    
    GENERALARGS=-ci "subGCG" -t "PRICING TIME" 10 20 "totaltime"

in the @ref comparison-report-settings "script settings" of the @ref comparison-report "comparison report".
The comparison report can be found in `check/reports/` in a uniquely named (using a time stamp) report folder.

### Logging Custom Data: Streaming Data
> **The custom logging feature (streaming data) is not yet implemented. Please stay tuned.**

## C. Comparing (unmodified) Versions {#compare-versions}
> In this part, we explain how one can **compare different (unmodified) versions of GCG** and its
> submodules with each other. This requires access to our git and thus is targeted at developers.

### Automated Version Comparison Script
The version comparison should be started from the chair's computers. The script is localized
in the check directory and can be called as follows:

    ./compareversions.sh "<GlobalFlags>" "<Version1>" "<Flags1>" "<Version2>" "<Flags2>" ...

The global flags are applied to all compilations (could include e.g. `STATISTICS=true`). The numbered
flags are only applied to the respective versions of GCG. Each string of flags must be in quotation marks.
The versions `Version1`, ... correspond to the version numbers that are tagged in git. If you want to
compare a version vX of GCG with a version vY of GCG (this time with CPLEX as solver), then you can execute   

    ./compareversions.sh "TEST=mytestset SETTINGS=mysettings -j" "vX" "" "vY" "LPS=cpx CPLEXSOLVER=true"

The version comparison report can be found in `check/results/` in a uniquely named (using a time stamp) report folder.

### Manual Version Comparison
In the results folder generated by the `compareversions.sh` script you can find files with ending `.pkl`.
In order to create a new report with existing runtime data, move those `.pkl` files into an own folder and call
the script `plotcomparedres.py` that is located inside the `check` folder:

    ./plotcomparedres.py <folder with pkl files> <output folder>

All files with `.pkl` ending in this folder will be concerned.
In this manual mode, you yourself have to pay attention to the general comparability of your runtime data
(same instances, settings, time limits, computer used, ...).

## D. Comparing different Solvers {#compare-solvers}
> In this part, we explain how one can **compare completely different solvers**, with the
> most likely example being GCG vs. SCIP.

For SCIP, the format of the output files should be similar, so the above instructions apply,
since most (if not all) visualization scripts will still function.
