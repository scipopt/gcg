# do automatic testing # {#dotesting}
Just like SCIP, GCG supports automatic testing. It uses the same command
structure as SCIP, so if any information is missing,
you can also consult the <a href="https://scip.zib.de/doc-6.0.1/html/">SCIP documentation</a>.

# Testing Basics #

### How do I test?
#### Simple Testing
If you want to check whether your own modules generate the right results,
you can work with the testset files `*.test` in the folder `/check/testset`
to specify which instances you want to test. The `*.solu` files specify their
respective solution. <i>For further syntax and creation information on the
`*.solu`-files, please check the
<a href="https://scip.zib.de/doc-6.0.1/html/TEST.php">SCIP documentation</a>.</i><br>
A simple test on `testrun.test` is executed with the following command
in the GCG main folder:

    make TEST=testrun test

Note that the default test is `short.test`, which will be executed if no
`TEST`-argument is given.

#### Testing with statistics
In order to use the test data to create visualizations, the testing has to be
done with statistics enabled. You can do that with:

    make TEST=testrun STATISTICS=true test


#### Testing with defined settings
To test with your own settings, e.g. `my.set` (located in `/settings`)
you can use the following command:

    make TEST=testrun SETTINGS=my test


### What files do I get?
In the folder `/check/results` you will see multiple files per test.
The output files for a test executed by
`make TEST=testrun SETTINGS=my STATISTICS=true test` are as follows:

    check.testrun.gcg-GCGVERSION.linux.x86_64.gnu.opt.spx2.COMPUTERNAME.my.default{.err,.out,.set}

and the files `{.res,.pav,.tex}` are added after the test was completed.
If you aborted the test before it finished, you can still generate those files
(incompletely) using the command `./evalcheck.sh` in the `/check` folder.

### What do the files contain?
The .out files contain the `stdout`-output and therefore all data that was
printed during the execution of GCG.<br>
The .res files only contain the execution table with
primal bound, dual bound, gap and so on.<br>
The latter file is also given as (incomplete) LaTeX code in the .tex file.<br>

# Using the results
The previous chapter only covered the absolute basics of testing, mainly used
for verification purposes. But heads up, the output files still have a purpose,
even if you verified that your GCG works (and how fast it works). At the moment,
you can do the following things:
- [Make overviews with the raw data](#raw)
- [Create Performance Plots to compare different runs](#plot)

Those things will be explained in the following subsections.<br>

#### Prerequisites
**Python:** <br>
Some scripts are just raw shell scripts, but some are also python
ones. For those, you need python to be installed on your system:

    sudo apt-get install python3

**Python-Packages:** <br>
For some python-scripts, you will need packages, that are probably
not yet installed on your computer. You can do that by first installing `pip`:

    sudo apt-get install python3-pip

And then the needed package via

    pip install <packackename>

### <a name="raw">Create Tables with raw data</a>
A (quite raw) comparison of testruns can be done using the `allcmpres.sh` script
in the `check`-folder. This script just puts the statistics of all runs that are
given as arguments into a `.tex`-file and prints it as ASCII on the console.
![Comparison Table (run with settings earlybranching was stopped)](cmpres.png)
#### Usage

    ./allcmpres.sh run1.res run2.res run3.res ...

with run1, run2, ... being a `.res` file in the format as shown in "What files do I get".

### <a name="plot">Plot performance profiles</a>
Using the tool `plotperprof.py`, you can plot performance profiles
(which will be saved as `perprof_plot.pdf` in the `check`-folder) like this one:
![Performance Profile](perfprofile.png)
#### Interpretation
The x-axis represents the factor by which the respective run is worse than the
optimal run, while the y axis is the corresponding probability. so for the green
line in `x=1.05`, we see `y≈0.7`, therefore, with a probability of 70%, setting 3
will be worse than setting 2 by 5%.

#### Usage

    ./plotperprof.py run1.res run2.res run3.res ...

with run1, run2, ... being a `.res` file in the format as shown in "What files do I get".

#### Issues
As described in the SCIP documentation, one should be able to run a test with
`SETTINGS="set1,set2"` and SCIP should run the test on both settings.
**This apparently does not work.**


### Experimental Scripts
**Script:** bounds_script.py<br>
**Error:** Lines with rootbounds in .out files don't start with `iter    pd   db` anymore<br>
<br>
**Script:** plotclassifier.py<br>
**Error:** Does not plot anything, same error as with bounds_script<br>
<br>
**Script:** pricing_statistics.py<br>
**Error:** No pricing information found (same as with bounds_script, just with pricing)<br>
<br>
**Script:** createInstanceFileMittelmann.py<br>
**Error:** Well, it probably works as intended, but the intention seems rather
obscure. It copies the first word of every line from instancesMittelmann.txt to instancesMittelmannClean.txt.
Whatever this is used for...
