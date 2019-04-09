Installation {#install}
------------------
[TOC]
This guide shows you how to install GCG on your machine. You can either
* \ref install-scipopt including GCG (recommended) or do a
* \ref install-manually.

# Install the SCIP Optimization Suite {#install-scipopt}

The following steps are all executed automatically within
<a href="../installGCG.sh">this install script</a>, if you prefer it.

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
    make install  # install scip executable, library and headers

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
documented on the official SCIP homepage under http://scip.zib.de. For further
options, consult the Makefile.


# Manual Installation {#install-manually}
Note that a manual installation is not recommended to new users. Its only benefit
is the compatibility to your own code, if it was made with makefile instead of
cmake.

### Step 1: Download everything
Installing SCIP, SoPlex, ZIMPL, hMETIS, Bliss, and GCG manually requires <a href="http://www.or.rwth-aachen.de/gcg">GCG</a>, <a href="http://scip.zib.de/">SCIP</a>, <a href="http://soplex.zib.de/">SoPlex</a>, <a href="http://zimpl.zib.de/">ZIMPL</a>, and additionally <a href="http://www.tcs.hut.fi/Software/bliss/">Bliss</a> and <a href="http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview">hMETIS</a>. Let assume that all source code archives were saved within one folder, i.e.

    ls
    bliss-[version].zip  gcg-[version].tgz  hmetis-[version]-linux.tar.gz  scip-[version].tgz  soplex-[version].tgz  zimpl-[version].tgz


### Step 2: Compilation and Links
#### ZIMPL and SoPlex
ZIMPL and SoPlex can be compiled by a simple `make`.

    tar xvfz zimpl-[version].tgz
    cd zimpl-[version]/
    make

    tar xvfz soplex-[version].tgz
    cd soplex-[version]/
    make

#### SCIP
Compiling SCIP requires links to both ZIMPL and SoPlex.

    tar xvfz scip-[version].tgz
    cd scip-[version]/
    make

#### Create Links
When asked for links to SoPlex, we need to set links from `scip-[version]/lib/spxinc` to `soplex-[version]/src/` and `scip-[version]/lib/libsoplex.*` to `soplex-[version]/lib/libsoplex.*`. For the linker, a link to either `libsoplex.*.a` or `libsoplex-*.so` is enough, so in our case we only need to specify the path to `libsoplex.*.a`:

    ../../../soplex-[version]/src/
    ../../../soplex-[version]/lib/libsoplex.linux.x86_64.gnu.opt.a

For ZIMPL, links need to point from `scip-[version]/lib/zimplinc/zimpl` to `zimpl-[version]/src/` and from `scip-[version]/lib/libzimpl.*` to `zimpl-[version]/lib/libzimpl.*`. Again, we may also ignore the link to `libzimpl.*.so`:

    ../../../zimpl-[version]/src/
    ../../../zimpl-[version]/lib/libzimpl.linux.x86_64.gnu.opt.a


Next, compile Bliss:

    unzip bliss-[version].zip
    cd bliss-[version]/
    make

Now, GCG can be compiled. The links to SCIP need to point from `gcg-[version]/lib/scip` to `scip-[version]/`.

    tar xvfz gcg-[version].tgz
    cd gcg-[version]/
    make BLISS=true
    ../../scip-[version]/
    ../../bliss-[version]/
    ../../bliss-[version]/libbliss.a

Finally, set a link from the GCG directory to hMETIS.

    tar xvfz hmetis-[version]-linux.tar.gz
    cd gcg-[version]/
    ln -s ../hmetis-[version]-linux/hmetis hmetis
