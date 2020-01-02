# The Visualization Suite # {#visu}

# Testset Report {#testset-report}
> This feature is still under development. Please be patient.
> Meanwhile, you can @ref generatevisu "use the scripts manually".

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
