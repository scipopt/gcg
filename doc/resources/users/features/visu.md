# The Visualization Suite # {#visu}

# Instance Runtime Report # {#runtimereport}
From inside GCG, one can generate a runtime report which will export a pdf file containing information about the time distribution and pricing.

> The report needs more information about the runtime than what is usually saved, thus as of now, it is required to rerun the solving on this instance.
> This will probably be fixed in later versions.

In order to make \GCG generate your report, you have to have your instance read into it. Since \GCG will run again on this instance (twice, to be exact), you don't need to have it optimized already. To export the pdf report, simply enter the write menu and type reportinstance just like

    write reportinstance

You will be prompted to give a folder to save the report into. After you gave one, \GCG will run twice (with different modes) on the instance and then start generating the visualizations. Those are written into the `img/` folder. Then, those images will be composed into one single latex file, which can be compiled using the shipped makefile.

> It is scheduled to add descriptions to each visualization in the pdf report.

The pages generated will then look similar to this:
\image html report_instance.png "A page inside the instance report." width=30%
<br/>
All visualizations created inside this report (and more, if a whole testset is given) can also be generated manually. A guide for that can be found under @ref generatevisu.

# Decomposition Report # {#detectionreport}


# Tree Visualizations # {#visu-tree}
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
