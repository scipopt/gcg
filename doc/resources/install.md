Installation {#install}
------------------
[TOC]
This guide shows you how to install GCG on your machine. You can either
* \ref install-scipopt including GCG (recommended) or
* \ref install-manually.

# Install the SCIP Optimization Suite {#install-scipopt}

The following steps are all executed automatically within
<a href="">this install script</a>, if you prefer it.

### Step 0: Prerequisites
Install the required system libraries:

    sudo apt-get install build-essential libreadline-dev libz-dev libgmp3-dev


### Step 1: Get the SCIP Optimization Suite
Download the archive and unzip it:

    wget https://scip.zib.de/download.php?fname=scipoptsuite-6.0.1.tgz
    tar scipoptsuite-6.0.1.tgz
    cd scipoptsuite-6.0.1.tgz


### Step 2: Compile via CMake
<i>(recommended for new users)</i><br/>
Create the build directory and compile the program:

    mkdir build
    cd build
    cmake ..
    make
    make install    # install scip executable, library and headers

And you're done! To test your installation, you can run a quick check on some instances:

    make check

### Step 2 (Alternative): Compile via Makefile
<i>(recommended only for developers, if necessary)</i><br/>
#### Creating softlinks

In order to create all necessary links, type

    make links

#### Compilation

You can compile GCG with

    make [options]

where the `[options]` are the same options as those from SCIP which is
documented on the official SCIP homepage unde


# Install SCIP, SoPlex, ZIMPL, hMETIS, Bliss, and GCG manually {#install-manually}
