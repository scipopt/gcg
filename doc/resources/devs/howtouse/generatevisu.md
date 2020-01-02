# Visualization Generation {#generatevisu}
[TOC]

# Prerequisites
### Software requirements
**Python:** \n
Except for the comparison_table script, all visualizations are generated using python.
For those, you need python3 to be installed correctly on your system:

    sudo apt-get install python3

\n
**Python-Packages:** \n
For each script, you will need packages and some are probably
not yet installed on your computer. You can do that by first installing `pip`:

    sudo apt-get install python3-pip

And then the needed package via

    pip install <packagename>

The packages needed for the scripts are (space-separated):

    numpy matplotlib pandas sklearn


### Run tests
First, you have to do a testrun to gather the data for the visualizations. A guide on how to
do that can be found @ref dotesting "here".\n
The scripts sometimes expect a certain flag to be set. The following table gives an overview over
the requirements for the runtime data for each visualization script.

|                       | Expected Input | Required  compile flags | Required test flags |
|-----------------------|----------------|-------------------------|---------------------|
| Performance Profile   | .out/.res      | -                       | -                   |
| General               | .out           |                         | STATISTICS=true     |
| Bounds                | .out           | -                       | STATISTICS=true     |
| Classifier/ Detection | .out           | -                       | MODE=detectionstatistics         |
| Pricing               | .out, .vbc     | STATISTICS=true         | -                   |
| Tree                  | .vbc           |                         | -                   |
| Comparison Table      | .res           | -                       | -                   |

# Plotting
> You can create different plots with the included scripts in the `stats` folder. All scripts are briefly explained in the following.
> When describing how to execute them, it is always assumed that you are inside the `stats` folder.

### Common arguments
The following arguments are common across all following scripts (except for the performance profile plotter and comparison table).
\n
#### Defining the Output Directory

    -o OUTDIR, --outdir OUTDIR
                          output directory (default: "plots")

\n
#### Parsing without Plotting, Plotting without Parsing

    -save, --savepickle      parses the given .out-file without plotting
    -load, --loadpickle      loads the given .pkl-file and plots it

\n

## Performance Profile Plotter
Using this plotter, one can generate performance profiles.
In those, the x-axis represents the factor by which the respective run is worse than the
optimal run, while the y-axis is the corresponding probability.
### Execution

    python3 general/performance_profile.py FILES

with `FILES` being some (space-separated) `.res` or `.out` files in the format as shown in @ref what-files "'What files do I get'".

\htmlinclude performance_profile.html
\n
## General Plotter
The general plotter is able to plot two or more arguments parsed by the `general_parser.py` script in different ways.
You can find all parsed arguments @ref general-args "here".
\n
### Execution

    python3 general/{bubble,plot,twin,time}.py [args] FILES

with `FILES` being either some pickles, or some outfiles, and [args] being as defined below. \n

\n
### Arguments
#### Defining Data to be plotted
> **Times** are the most common arguments to plot (that's the reason for the naming),
> but you can use whatever argument parsed you want. You can find a list of those @subpage visu-args "here".

You can define what data to plot with the `--times` argument:

    -t [TIMES [TIMES ...]], --times [TIMES [TIMES ...]]
                          times to be plotted.

Whereas
- for `bubble.py` and `twin.py`, you may define exactly two times (or none), \n
- for `plot.py` you can define as many times as you like and
- for `time.py` you can also define an arbitrary amount, but the "OTHERS" time will be added automatically.

\n
#### Defining Filters for the Data
For all arguments, you can define filters. For example, if you only want to plot
all those instances that had a Detection Time of between 10 and 20 seconds, you simply
add those bounds to the time argument:

    python3 general/time.py check.gcg.out -A -t "DETECTION TIME" 10 20 "RMP LP TIME"

In this example, we generate all plots (`-A`) for all instances that had a detection time
between 10 and 20 seconds in the `check.gcg.out` file with the arguments `"DETECTION TIME"`
and `"RMP LP TIME"` (without filter).

\n
#### Defining additional Arguments (only for time.py)
For the `time.py` script, you can define which plots you want to generate. The arguments
are as follows:

    -A, --all             create all visualizations
    -B, --bar             create barchart
    -G, --grouped-bar     create barchart, grouped by CPU times
    -P, --plot            create simple plot
    --pie                 create a pie chart plot (will automatically activate --single (see below))

Additionally, for the grouped-bar plot, you may define a number of buckets,
in which the script will automatically sort the (nearly) same number of instances into.

    --buckets BUCKETS     amount of buckets (resolution of grouped_bar plot)

Finally, you can compare two different runs. Though recommended, it is not required for the
runs to have been on the same testset.

    --single              if set, all outfiles will be summed and cumulated in a single plot
    --compare             if set, each outfile will be summed and plotted with all other outfiles


\htmlinclude visualizations/general.html
\htmlinclude visualizations/time.html
\n
## Bounds Plotter
The Bounds Plotter generates a plot showing the development of the primal and dual bound and gap in the root node, as well as the basic variables generated.
#### Execution

    python3 bounds/plotter_bounds.py FILES

with `FILES` being one or more `.out` files.

\htmlinclude visualizations/bounds.html
\n
## Classification and Detection Plotter
This plotter works on a whole testset and makes plots similar to performance profiles, showing the performance of the classifiers and detectors.
#### Execution

    python3 detection/plotter_detection.py FILES

with `FILES` being one or more `.out` files.

\htmlinclude visualizations/detection.html
\n
## Pricing Plotter
The Pricing Plotter generates 7 different plots illustrating the pricing procedure during a single instance's solving process.
When given an outfile with more than one instance, it generates the plots sequentially.
#### Execution

    python3 pricing/plotter_pricing.py FILES --vbcdir VBC

with `FILES` being one `.out` file and `VBC` being the directory where all corresponding `.vbc` files are (per default: `check/results/vbc/`)

\htmlinclude visualizations/pricing.html
\n
## Tree Plotter
The Tree Plotter, just like the Pricing Plotter, needs the `vbc` files to function correctly. It will plot how many nodes were opened on each level.
#### Execution

    python3 tree/plotter_tree.py FILES

with `FILES` being some `.vbc` files.

\htmlinclude visualizations/tree.html
\n

## Comparison Table
A (quite raw) comparison of testruns can be done using this script
in the `general`-folder. This script just puts the statistics of all runs that are
given as arguments into a `.tex`-file and prints it as ASCII on the console.
**Execution**

    ./general/comparison_table.sh run1.res run2.res run3.res ...

with run1, run2, ... being a `.res` file in the format as shown in @ref what-files "'What files do I get'".

\htmlinclude visualizations/table.html
\n

# Test set selection
> A complete guide on this is coming soon.

Test sets can be filtered using the instructions under "General Plotter -> Arguments -> Defining Filters for the Data".

# Custom Visualizations
When creating custom visualizations, one has to know exactly what data is needed to make the visualization. With these arguments in mind, one can then look if they are already parsed. A list of the currently parsed data is located @ref visu-args "here". If so, one of the parsers (<code>parser_general.py</code>, <code>parser_bounds.py</code> or <code>parser_detection.py</code>) can be used, or otherwise, for the <code>plotter_pricing.py</code>, the <code>--savepickle</code> argument of the plotter shall be used each time to parse the runtime data and save it to the pickle, to then read it again for the plotter. Then, the <code>plotter_</code> script can be created which should import the parser that already gets the data needed with a simple <code>import parser_...</code>. Finally, the parser can be used just like in the other plotters.


### Troubleshooting
**Q: Why don't I get any detection times?**\n
A: You probably did not run the test with `MODE=0` (use `.dec` files instead of detecting).
