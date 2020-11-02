# Git Installation {#git-install}
# Install GCG using Git
> This guide aims at helping you to install GCG from Git. Note that you have to have access to the repository 
> (it is not public). If you want to let an automated installer do the work for you, you can download it 
> [here](installGCGfromGit.sh). Your SSH key has to be registered correctly for it to work.

## Git prerequisites
Please make sure that you have your ssh key saved in your Git account. Otherwise, you won't be able to
clone or commit anything, since you have not authorized yourself.
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

Compile the depends (scip, soplex, bliss and googletest) with

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
 * lib/zimplinc/zimpl is `PATH_TO_ZIMPL_DIR/src/zimpl` (PATH_TO_ZIMPL_DIR=/opt/scipoptsuite-X.X.X/zimpl/ on orlabXX; compile with argument `ZIMPL=false` if not needed)
 * lib/libzimpl.linux.x86_64.gnu.opt.a is `PATH_TO_ZIMPL_DIR/lib/libzimpl.linux.x86_64.gnu.opt.a`
 * lib/libzimpl.linux.x86_64.gnu.opt.so is not needed

## Compile GCG
Afterward, compile GCG with

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
