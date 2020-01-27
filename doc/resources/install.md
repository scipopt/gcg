Installation {#install}
------------------

> This guide shows you how to install GCG on your **Linux or Mac** computer. Note that
> GCG cannot be installed on a **Windows** machine out of the box (see @ref windows-install "Windows Installation Guide" for more information).
> As a user, you can either
> * @ref install-scipopt including GCG (recommended) or do a
> * @ref install-manually.<br>
>
> If you are a GCG developer and need the current master branch, consult
> * @ref git-install

If you need to install GCG using arguments, for example if you want to enable parallelization, you can find them under @ref makefiles-args and @ref cmake-args respectively.

# Install the SCIP Optimization Suite {#install-scipopt}

The following steps are all executed automatically within
<a href="../scripts/installGCG.sh">this install script</a>, if you prefer it. Note that we don't take any responsibility if you use it.

### Step 1: Prerequisites
Install the required system libraries:

    sudo apt-get install build-essential libreadline-dev libz-dev libgmp3-dev lib32ncurses5-dev


### Step 2: Get the SCIP Optimization Suite
Download the most recent SCIPOptSuite from the [SCIP Page](https://scip.zib.de/index.php#download).
Unzip it and go into the folder, replacing `X.X.X` by the version you downloaded.

    tar xvzf scipoptsuite-X.X.X.tgz
    cd scipoptsuite-X.X.X


### Step 3: Compile via CMake
<i>(recommended for new users)</i><br/>
First, make sure that cmake is installed correctly and up-to-date:

    sudo apt-get install cmake

Then create the build directory and compile the program
(you should execute them as root, if possible, to avoid certain errors):

    mkdir build
    cd build
    cmake ..
    make
    make install  # install scip executable, library and headers

And you're done! To test your installation, you can run a quick check on some instances:

    make check

Note that the execution of `make`-commands, e.g. `make test` is only supported
inside the `build`-folder (in which it requires some more steps to add testsets).
Therefore, if you intend on using `make test` (and not ctest), you should compile
SCIP and GCG via Makefile.

### Step 3 (Alternative): Compile via Makefile
<i>(recommended for developers, if necessary)</i><br/>
#### Creating softlinks

In order to create all necessary links, type

    make links

#### Compilation

In the root folder of the scipoptsuite, you can then compile everything
(you should execute them as root, if possible, to avoid certain errors) with:

    make scipoptsuite
    make gcg
    make install

Arguments to be added to `make gcg` if needed can be found under \ref makefiles-args.

### Step 4: Use GCG!
See our \ref getting-started - Guide for further comments on the usage.
