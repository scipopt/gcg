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

# Comparing Runtime Data
After collecting the data, you can start comparing it. In the following, we will give
some hints on how to start, depending on what you want to compare. 
The most common cases are comparing default GCG...

A. against a GCG with different **settings** (continue @ref compare-settings "here") \n
B. against an **algorithmically** modified GCG version (continue @ref compare-code "here") \n
C. against a **different version** of default GCG (continue @ref compare-versions "here") \n
D. against a **different solver** (continue @ref compare-solvers "here") \n


Depending on what you want to test, skip to the respective section.

## A. Comparing Settings {#compare-settings}
## B. Comparing Code {#compare-code}
## C. Comparing (unmodified) Versions {#compare-versions}
## D. Comparing different Solvers {#compare-solvers}