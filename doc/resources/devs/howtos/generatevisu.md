# generate visualizations {#generatevisu}

### Prerequisites
#### Run tests
First, you have to do a testrun to gather the data for the visualizations. A guide on how to
do that can be found under \ref dotesting.
**Note:**
- for scripts in "Plotting (1)": compile GCG with `STATISTICS=false`, run tests with `STATISTICS=true`,
- for scripts in "Plotting (2)", compile GCG with `STATISTICS=true`, run tests with `STATISTICS=true`.

#### Copy vbc files
For the `bounds_script.py` to work, you have to move your vbc folder from check/results/
to stats/visuscripts/. You can do that with `cp -r check/results/vbc stats/visuscripts`.

#### Software requirements for the scripts
**Python:** <br>
Some scripts are just raw shell scripts, but some are also python
ones. For those, you need python to be installed on your system:

    sudo apt-get install python3

**Python-Packages:** <br>
For some python-scripts, you will need packages, that are probably
not yet installed on your computer. You can do that by first installing `pip`:

    sudo apt-get install python3-pip

And then the needed package via

    pip install <packagename>


### Plotting (1)
You can create different plots with the included scripts in the `stats/` folder. Every script
needs to be executed with at least one pickle or outfile as argument, just as following:

  python3 visuscripts/{bubble,plot,twin,time}.py [args] FILES

with FILES being either some pickles, or some outfiles. \n
**Common arguments** for all scripts are the following:

    -h, --help            show this help message and exit
    -o OUTDIR, --outdir OUTDIR
                          output directory (default: "plots")

All other arguments, along with a short description, are listed in the following subsections.

#### bubble.py
**Arguments to choose data to create plots from:**

    -t [TIMES [TIMES ...]], --times [TIMES [TIMES ...]]
                          times to be plotted. Use exactly two times, all others
                          will be ignored.

#### plot.py
**Arguments to choose data to create plots from:**

    -t [TIMES [TIMES ...]], --times [TIMES [TIMES ...]]
                          times to be plotted. Use arbitrary amount.

#### twin.py
**Arguments to choose data to create plots from:**

    -t [TIMES [TIMES ...]], --times [TIMES [TIMES ...]]
                          times to be plotted. Use exactly two times, all others
                          will be ignored.


#### time.py
Creates a time distribution plot for a testrun, showing where the instances
needed the most time.
**Arguments to choose data to create plots from:**

    -t [TIMES [TIMES ...]], --times [TIMES [TIMES ...]]
                          times to be plotted. "OTHERS" time will be added
                          automatically. The first argument is the one to be
                          sorted by.

**Arguments to choose plots to be generated:**

    -A, --all             create all visualizations
    -B, --bar             create barchart
    -G, --grouped-bar     create barchart, grouped by CPU times
    -P, --plot            create simple plot

**Argument to change appearance of the plots:**

    --buckets BUCKETS     amount of buckets (resolution of grouped_bar plot)


### Plotting (2)
The following scripts have to be executed on an outfile, since they use a different
pickle format.

#### pricing_statistics.py
Plots different statistics about pricing.
execute with argument `--complete-only`

#### bounds_script.py
Plots development of primal and dual bound, as well as added LP and IP vars.

### Plotting (3)

### <a name="raw">Create Tables with raw data</a> (allcmpres.sh)
A (quite raw) comparison of testruns can be done using the `allcmpres.sh` script
in the `check`-folder. This script just puts the statistics of all runs that are
given as arguments into a `.tex`-file and prints it as ASCII on the console.
![Comparison Table (run with settings earlybranching was stopped)](cmpres.png)
#### Usage

    ./allcmpres.sh run1.res run2.res run3.res ...

with run1, run2, ... being a `.res` file in the format as shown in "What files do I get".

### <a name="plot">Plot performance profiles</a> (plotperprof.py)
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

### Experimental scripts
#### plotclassifier.py
This script does not plot anything at the moment.

#### createInstanceFileMittelmann.py
Copies the first word of every line from instancesMittelmann.txt to instancesMittelmannClean.txt.

### Other tools
#### Parse outfiles without plotting
For the scripts `bubble.py`, `plot.py`, `twin.py` and `time.py`, the outfile(s)
are first parsed, if you don't already execute them with a `.pkl` file.
If you want to parse your outifles by hand, you can do that by executing

  python3 visuscripts/parseout.py myoutfile.out

in the `stats/` folder. A `parseout.pkl` file will be saved into this same folder.

#### Testset selection
See pdf

### Troubleshooting
**Q: Why don't I get any detection times?**
A: You probably did not run the test with `MODE=0` (use .dec files instead of detecting).
