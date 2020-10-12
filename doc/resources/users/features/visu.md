# The Visualization Suite # {#visu}

# Test Set Report {#testset-report}
> **Note: This guide concerns the branch `552-testset-report` and the mentioned features are only available on this branch.**

The test set report **generates a PDF that includes all visualizations**, both for all single instances and for aggregated
statistics for the test set, to give you a compiled and captioned overview of visualizations of  different aspects of GCG's algorithmics.\n
The report can be generated in two ways:
1. Execute a test and directly generate the report (continue @ref testset-report-auto "here")
2. Have runtime data collected already and generate the report afterwards (continue @ref testset-report-manual "here")

## Automatically generate a Test Set Report {#testset-report-auto}
> The automatic generation is very simple and straightforward and recommended, especially for beginners. 
> However, it will always (re-)generate all runtime data, which might be undesirable, in particular for very big test sets.
> For a different method, see @ref testset-report-manual "here".
### Preparation
You do not have to prepare anything except for having **compiled GCG with `make STATISTICS=true`**. This is required because
GCG will only then print out extensive pricing and bounds statistics.\n
Before starting the report generation, you can create a script settings file, e.g. `settings.scset`. In there, you can 
define which plots to generate and which arguments to generate them with. If you have a big testset or just want to check
if the feature works, you can also enable a draft mode that will generate a reduced version of the report. A full
list of arguments can be found @subpage report-args "here". To use the script settings, execute the test with 
`SCRIPTSETTINGS=settings.scset`, with having `settings.scset` lying in the GCG root directory.
### Generation
When wanting to generate the report automatically, you can simply perform a `make test`. It will support all flags, just
as you are used to when doing tests, and the data generated during the test will remain in the usual place. 
Now, to generate the report, please add the flag `VISU` to your command, i.e. 

    make test VISU=true

This will make the test script call the report generator script afterwards, already with the correct arguments.
### Output
As output, the script will create a folder called `report_<testsetname>_<settingsname>_<timestamp>`. Inside this folder,
you will find a folder `logs` including the out-file, res-file and vbc-folder of your test run. Next, you can find the 
folder `plots`, where, subdivided into the plots types, the visualizations are located. Finally, **the report
is inside the directory** as a (LaTeX-generated) PDF, already compiled for you and opened automatically.

## Manually generate a Test Set Report {#testset-report-manual}
> The manual generation requires setting many variables. 
> However, leaving them undefined will only lead to the report having an empty place at the respective table entry for
> most arguments.
### Preparation
First, as above, your (already collected) runtime data has to have been generated with a **GCG version compiled with
`STATISTICS=true`** and also **tested with `STATISTICS=true`** to yield a complete report.\n
Then, in addition, you should know the parameters with which you executed the test. They have to be defined either
via the arguments at the required positions (not recommended) or via the script settings file, as described in the 
section "Preparation" for the automatic report. 
A full list of arguments can be found @ref report-args "here", where it is also described which flags _have_ to be defined for the scripts
to yield results at all and those that can be defined for the report to have full information. As before, to use the script settings, 
execute the test with `SCRIPTSETTINGS=settings.scset`, with having `settings.scset` lying in the GCG root directory.
### Generation
In order to generate a **test set report _without_ testing**, you have to call the script directly. This can be done using

    make test VISU=true SCRIPTSETTINGS=settings.scset

### Output
You will get the same output as with the automated report generation.

# Comparison Report {#comparison-report}
> This feature is coming soon. Stay tuned.

# Tree Visualizations {#visu-tree}
In order to generate pictures of the Branch and Bound tree that GCG used during solving, you can use the [vbctool](https://informatik.uni-koeln.de/ls-juenger/vbctool/). Since the executable might have issues with the linking of the libraries, it is suggested to download the source code (and additionally the Motif Framework). Before you start with the [Build Instructions](https://informatik.uni-koeln.de/fileadmin/projects/vbctool/INSTALL), you have to install a packet:

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

The files you now have to read (File -> Load) are included in the folder `check/results/vbc`.

\image html tree.jpg "A tree." width=80%

In order to generate the tree, click on Emulation -> Start. Before doing that, you can configure the emulation in Emulation -> Setup, where you can also set the time it will need to generate the tree.
- If left on default values, the tree will generate as fast as it generated in GCG during execution, offering you a good insight into how long GCG was 'stuck' in certain nodes.
- If changed, for example to 1 second, it will just generate the tree all at once and you can then save it.

To save the generated tree, just click on File -> Print and it will save a `.ps` file.
