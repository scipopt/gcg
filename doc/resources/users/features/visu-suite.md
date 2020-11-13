# The Visualization Suite {#visu-suite}

[TOC]

> **Note: The features described on this page are available as a preview only on the branch `visu_suite`.**\n

> On this page, we **present and guide through the main visualization reporting capabilities** of GCG. 
> If you want to generate a specific visualization and no complete report, please visit the guide on @ref visu-manual.
> If you are looking for visualizations of decompositions please visit the guide @ref explore-menu. 

In general, we want to visualize runtime and algorithmic behaviour, either for single test runs (@ref testset-report) 
or for comparisons between different runs (@ref comparison-report). The scripts can automatically conduct test(s), 
generate all visualizations and afterwards compile them into a captioned
report. Furthermore, we give a brief guide on how to use a tool to visualize branch-and-cut algorithmics, graphically
showing how the tree was built during branching.

- @ref testset-report
- @ref comparison-report
- @ref vbc-visu


##### A Remark for CMake Users
Generally, we don't recommend to do testing and experiments using a CMake build, which is due to `make check` not generating 
the usual runtime data, see @ref install-manually for more information). Therefore, automatic test set report generation 
is currently **not supported under CMake builds**. When using the @ref testset-report-manual "manual mode" instead, 
you have to give runtime data generated with a `make` installation using `make test` or with a CMake build using the
`cmake gcg_cluster` target.
 
# Test Set Report {#testset-report}
The test set report **generates a PDF that includes all visualizations**, both for all single instances and for aggregated
statistics for the test set, to give you a compiled and captioned overview of visualizations of different aspects of GCG's algorithmics.\n
The report can be generated in two ways:
1. Execute a **test and generate** the report (continue @ref testset-report-auto "here")
2. Have runtime data **collected already and generate** the report afterwards (continue @ref testset-report-manual "here")

## Automatically generate a Test Set Report {#testset-report-auto}
> The automatic generation is very simple and straightforward and recommended, especially for beginners. \n
> However, **it will always (re-)generate all runtime data**, which might be undesirable, in particular for very big test sets.
> For a different method, see @ref testset-report-manual "here". \n

### Preparation
- For a full report, it is required to have **compiled GCG with `make STATISTICS=true`**. Otherwise, GCG will not print 
out extensive pricing and bounds statistics for the respective visualizations.\n
- All **other requirements** listed on the @ref visu-prerequisites "visualization script page" have to be fulfilled 
(e.g. correctly configured python). 
- The script **will not print any warnings** (e.g. if your python is not configured correctly) 
unless you set `DEBUG=true` in the settings file (see below). 
- You **can create a script settings file**, e.g. `settings.vset` where you can define which plots to generate and which 
arguments to generate them with. Furthermore, if you have a big testset or just want to check if the feature works, 
you can also enable a draft mode that will generate a reduced version of the report. A full list of possible settings can 
be found @subpage report-settings "here".
- To **use the script settings**, execute the test with `VISUSETTINGS=settings.vset`, while having `settings.vset` 
lying in the GCG root directory.

### Generation
When wanting to generate the report automatically, you can work with the usual `make test`. It will support all flags, just
as you are used to when doing tests, and the data generated during the test will remain in the usual place.
To tell GCG to generate the report after the test has been executed, please add the flag `VISU` to your command, i.e. 

    make test VISU=true

This will make the test script call the report generator script afterwards, already with the correct arguments.

### Output
> **An example of the output generated can be found [here](pdf/short.default.testsetreport.pdf).**

As output, the script will create a folder called `testsetreport_<testsetname>_<settingsname>_<timestamp>` inside `check/reports`,
unless defined otherwise inside the script settings file. Inside the report folder, you will find a folder `logs` including 
the out-file, res-file and vbc-folder of your test run. Next, you can find the folder `plots`, where, subdivided into the 
plots types, the visualizations are located. Finally, **the report is inside the directory** as a (LaTeX-generated) PDF, 
already compiled for you and opened automatically.


## Manually generate a Test Set Report {#testset-report-manual}
> The manual generation requires you to have the **runtime data available already**.\n
> You can either give outfile, resfile and vbc folder via the script settings file, else you will be asked for them during runtime.

### Preparation
- For a full report, it is required to have generated your runtime data with both, **a GCG compiled with `make STATISTICS=true`**
and a **test executed with statistics `make test STATISTICS=true`**. Otherwise, GCG will not print out extensive pricing and 
bounds statistics for the respective visualizations. Furthermore, **for detection visualizations, `DETECTIONSTATISTICS=true`** 
must be set during testing.
- All **other requirements** listed on the @ref visu-prerequisites "visualization script page" have to be fulfilled 
(e.g. correctly configured python). 
- The script **will not print any warnings** (e.g. if your python is not configured correctly)
unless you set `DEBUG=true` in the settings file (see below). 
- You **have to create a script settings file**, e.g. `settings.vset`, where you have to define parameters of your
given test run. A full list of possible settings can be found @ref report-settings "here", where it is also described which flags 
_have_ to be defined for the scripts to yield results at all and those that can be defined for the report to have full information.
Note that for `make`, none of the optional flags have to be set - they will be found automatically. 
- To **use the script settings**, execute the test with `VISUSETTINGS=settings.vset`, while having `settings.vset` 
lying in the GCG root directory.

### Generation
In order to generate a **test set report _without_ testing**, you can call the visualization target. This can be done using

    make visu VISUSETTINGS=settings.vset DATADIR=<folder containing .res and .out file and vbc file folder>

The data directory can also be given via the settings file. If you give none, you will be asked for it when the script is called.

### Output
You will get the same output as with the automated report generation. If visualizations are missing, please set the debug flag
```
DEBUG=true
```
inside your script settings file to troubleshoot.

# Comparison Report {#comparison-report}
The comparison report **generates a PDF with a comparison table and all comparing visualizations**, 
to give you a compiled and captioned overview of visualizations of different aspects of your two or more test runs.\n
**Note:** We do not support an automated mode natively, check the @ref compare-versions "version comparison script" for 
a possibility to automatically generate the runtime data.

### Preparation
- For a full report, it is required to have generated your runtime data with both, **a GCG compiled with `make STATISTICS=true`**
and a **test executed with statistics `make test STATISTICS=true`**. Otherwise, GCG will not print out extensive  bounds statistics
 for the respective visualizations. Furthermore, **for detection visualizations, `DETECTIONSTATISTICS=true`** 
must be set during testing.
- All **other requirements** listed on the @ref visu-prerequisites "visualization script page" have to be fulfilled 
(e.g. correctly configured python), and additionally, you have to have the python package `tikzplotlib`.
- The script **will not print any warnings** (e.g. if your python is not configured correctly) 
unless you set `DEBUG=true` in the settings file (see below). 
- You **can create a script settings file**, e.g. `settings.vset` where you can define which plots to generate and which 
arguments to generate them with. Furthermore, if you have a big testset or just want to check if the feature works, 
you can also enable a draft mode that will generate a reduced version of the report. A full list of possible settings can 
be found @ref report-settings "here".
- To **use the script settings**, execute the test with `VISUSETTINGS=settings.vset`, while having `settings.vset` 
lying in the GCG root directory. By default, the script will assume that you conducted the tests with the flag that was set
when you _last compiled this binary_. You can overwrite that by giving `LAST_STATISTICS=true` in the script settings file.

### Generation
In order to generate a **comparison report** (always with already existing runtime data), you can call the visualization target,
but with a data directory that contains multiple runs:

    make visu VISUSETTINGS=settings.vset DATADIR=<folder containing multiple .res and .out files>

The data directory can also be given via the settings file. If you give none, you will be asked for it when the script is called.

### Output
> **An example of the output generated can be found [here](pdf/comparisonreport.pdf).**

As output, the script will create a folder called `comparisonreport_<timestamp>` inside `check/reports`,
unless defined otherwise inside the script settings file. Inside the report folder, you will find a folder `logs` including 
the out-files and res-files of your test runs, as well as a `pickle` folder that includes dataframes with your logs, but parsed. 
Next, you can find the folder `plots`, where, subdivided into the plots types, the _comparing_ visualizations are located. 
Finally, **the report is inside the directory** as a (LaTeX-generated) PDF, 
already compiled for you and opened automatically.

If visualizations are missing, please set the debug flag
```
DEBUG=true
```
inside your script settings file to troubleshoot.


# Tree Visualizations {#vbc-visu}
> Note: The following guide concerns external software. We do not provide warranty nor support for it.

### Installation
In order to generate pictures of the Branch and Bound tree that GCG used during solving, you can use the [vbctool](https://cs.uni-koeln.de/ls-juenger/software/vbctool). Since the executable might have issues with the linking of the libraries, it is suggested to **download the vbctool source code** and additionally the **Motif Framework source code**, both available on the website. Unzip the Motif Framework source code tarball into the `lib/` folder of the vbctool. 
Before starting with the [Build Instructions](https://cs.uni-koeln.de/sites/informatik/projects/vbctool/INSTALL), install the following packages:

    sudo apt-get install libmotif-dev libxext-dev

Then, compile the program (just like explained in the Build Instructions):

    cd lib/GFE
    make
    cd ../GraphInterface/
    make
    cd ../MotifApp/
    make
    cd ../..
    make

Now you can start the program using

    ./vbctool

### Usage
The files you now have to read (File -> Load) are included in the folder `check/results/vbc`.

\image html tree.jpg "A tree." width=80%

In order to generate the tree, click on Emulation -> Start. Before doing that, you can configure the emulation in Emulation -> Setup, where you can also set the time it will need to generate the tree.
- If left on default values, the tree will generate as fast as it generated in GCG during execution, offering you a good insight into how long GCG was 'stuck' in certain nodes.
- If changed, for example to 1 second, it will just generate the tree all at once and you can then save it.

To save the generated tree, just click on File -> Print and it will save a `.ps` file.
