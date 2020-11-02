# Easy Installation {#easy-install}

> This installation guide will only function for Linux and Mac computers out of the box. 
> If you only have a Windows computer, please see @ref windows-install "Windows Installation Guide" for more information.  

# Install GCG and the SCIP Optimization Suite {#install-scipopt}
### Step 1: Prerequisites
Install the required system libraries:

    sudo apt-get install build-essential libreadline-dev libz-dev libgmp3-dev lib32ncurses5-dev libboost-all-dev


### Step 2: Get the SCIP Optimization Suite
Download the most recent SCIPOptSuite from the [SCIP Page](https://scipopt.org/index.php#download).
Unzip it and go into the folder, replacing `X.X.X` by the version you downloaded.

    tar xvzf scipoptsuite-X.X.X.tgz
    cd scipoptsuite-X.X.X


### Step 3: Compile via CMake
> The installation via CMake is **recommended for new users**. It does
> not offer the testing capabilities that Makefiles offer out-of-the-box.

First, make sure that cmake is installed correctly and up-to-date:

    sudo apt-get install cmake

Then create the build directory and compile the program
(you should execute them as root, if possible, to avoid certain errors):

    mkdir build
    cd build
    cmake ..
    make
    make install  # install scip executable, library and headers

And you're done! To test your installation, you can run a quick check 
(from within the build directory) on some instances:

    make gcg_check

Arguments to be added to `cmake ..` if needed can be found under @ref cmake-args.

Note that the execution of `make`-commands, e.g. `make test` is only supported
inside the `build`-folder (in which it requires some more steps to add testsets).
Therefore, if you intend on using `make test` (and not ctest), you should compile
SCIP and GCG via Makefile.

### Step 3 (Alternative): Compile via Makefile
> The installation using the Makefiles build system is recommended for developers,
> since it offers more testing capabilities (see @ref conduct-experiments).

In the root folder of the SCIP Optimization Suite, you can then compile everything
(you should execute them as root, if possible, to avoid certain errors) with:

    make scipoptlib
    make gcg
    make install

Arguments to be added to `make gcg` if needed can be found under @ref makefiles-args.

### Step 4: Use GCG!
See our @ref getting-started - Guide for further comments on the usage.
