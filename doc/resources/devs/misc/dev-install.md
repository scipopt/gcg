# Advanced Installation {#dev-install}

# Installation using Git {#git-install}
## Git prerequisites
First, clone the git repo with

    git clone git@git.or.rwth-aachen.de:gcg/gcg.git [gcg-dir]

You **have** to clone it via ssh, otherwise the submodules won't work.
To initialize the SCIP, SoPlex, bliss and the googletest framework, goto the repository folder "gcg" and run

    git submodule init
    git submodule sync
    git submodule update

## System prerequisites
If you are not working on one of the OR chair's computers make sure that the following libraries are installed:

```or
gcc
gpp
git
cmake
build-essential
libgmp-dev
libreadline-dev
zlib1g-dev
bison
flex
libncurses-dev
```

## Compile SCIP, SOPLEX, BLISS using makefiles
**Note**: Next, ```make``` is used with some arguments that you prefer. Do not use the ```-j``` option on the very first compilation as it is not compatible with the linker. As the linker is not called again once all links are set, using this option on future compilations should be fine.

It should work to compile the depends (scip, soplex, bliss and googletest) with

    make [args] deps

## Links to submodules
Before compilation, you will be asked for some links. Paste the following paths:
 * lib/scip is `../lib/scip-git`
 * lib/include/spxinc is `../../../soplex-git/src/`
 * lib/static/libsoplex.linux.x86_64.gnu.opt.a is `../../../soplex-git/lib/libsoplex.linux.x86_64.gnu.opt.a`

If compiled without flags, this should have been it.<br>
**If you used flags, you might need one of these links:**
 * lib/libsoplex.linux.x86_64.gnu.opt.so is not needed
 * lib/cpxinc is in `PATH_TO_CPLEX_DIR/include/ilcplex/` (PATH_TO_CPLEX_DIR = /opt/cplex/cplex/ on orlabXX)
 * lib/libcplex.a is in `PATH_TO_CPLEX_DIR/lib/x86-64_sles10_4.1/static_pic/libcplex.a`
 * lib/libcplex.linux.x86_64.gnu.so is not needed
 * lib/zimplinc/zimpl is `PATH_TO_ZIMPL_DIR/src` (PATH_TO_ZIMPL_DIR=/opt/scipoptsuite-X.X.X/zimpl-X.X.X/ on orlabXX; compile with argument `ZIMPL=false` if not needed)
 * lib/libzimpl.linux.x86_64.gnu.opt.a is `PATH_TO_ZIMPL_DIR/lib/libzimpl.linux.x86_64.gnu.opt.a`
 * lib/libzimpl.linux.x86_64.gnu.opt.so is not needed

## Compile GCG
Afterwards, compile GCG with

    make [args]

You can set arguments as described in \ref makefiles-args. <br>
**If you used flags, you might need to set these links:**
 * lib/include/bliss is `../bliss-git/`
 * lib/static/libbliss.a is `../bliss-git/libbliss.a`
 * lib/cliquerinc is `cliquer/`
 * lib/libcliquer.a is `cliquer/libcliquer.a`
 * lib/gtest is `googletest/include/`
 * lib/libgtest.a is `googletest/build/libgtest.a`

## Test GCG
You can run a short test with

    make [args] test


## Common Errors
On some distros, including the one used on RWTH cluster, the SCIP link does not work. Do this before compiling:

    cd lib && ln -s scip-git scip && cd ..

You can find all arguments for the installation under \subpage makefiles-args and \subpage cmake-args .

# Manual Installation {#install-manually}
Note that a manual installation is not recommended to new users.

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
