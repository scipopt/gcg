# Git Installation {#git-install}
> This page guides you to through the **installation of GCG from Git**. Note that you have to have access to the repository 
> (it is not public). 

# Install GCG using Git
> **If you want to let an automated installer do the work for you, you can download it [here](installGCGfromGit.sh).**\n
> **Note:** Your SSH key has to be registered correctly for it to work. Also, it might be required to make the script
> executable using `sudo chmod +x installGCGfromGit.sh`.

@htmlonly
<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/nsTiXuPW1WE" style="margin:auto; display:block" frameborder="3"  allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
@endhtmlonly

## Automated Installation
We generally recommend to use the [automated installer](installGCGfromGit.sh), since it is updated regularly. You can execute it like that:
```
./installGCGfromGit.sh [additional options]
```
The following options can be set:
```
  -d, --directory  folder to install GCG to (default: 'gcg')
  -s, --system     buildsystem (make/cmake, default: 'make')
  -b, --branch     branch to clone (default: 'master')
  -f, --flags      flags to use for the compilation (for the whole suite)
  -m, --mode       
                   default: clone everything (recommended)
                   fast:    only clone given branch
                   fastest: only clone given branch and newest submodules
```
Make it executable (`sudo` might be required):
```
chmod +x installGCGfromGit.sh
```

## Troubleshooting
In case of compilation errors, please first try this [bugfix script](uploads/260100a47b3b816bc146e974cf14b2a2/bugfixer.sh). Note that GCG has to have been compiled before executing it.
If it does not help, do the following manually:
1. Try to perform a clean make:
```
make clean
make deps
make
```
If that does not help and an `ld` error is returned, you should

2. Remove the links and submodules, since they could be erroneous. Note that you will lose every code changes in the SCIP folder or other submodules located in `lib`.
```
rm -rf lib/
git submodule update
make deps
make
```

3. If that also does not help, consider cloning the repository again, since with the installation script, it will most probably work again.

4. If it does not: Please report that in an issue or in the meeting.

## Manual Installation
### Git prerequisites
Please make sure that you have your ssh key saved in your Git account. Otherwise, you won't be able to
clone or commit anything, since you have not authorized yourself.
First, clone the git repo with

    git clone git@git.or.rwth-aachen.de:gcg/gcg.git [gcg-dir]

You **have** to clone it via ssh, otherwise the submodules won't work.
To initialize the SCIP, SoPlex, bliss and the googletest framework, goto the repository folder "gcg" and run

    git submodule init
    git submodule sync
    git submodule update

### System prerequisites
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
libboost-program-options-dev
```

### Compile SCIP, SOPLEX, BLISS using makefiles
**Note**: Next, ```make``` is used with some arguments that you prefer. Do not use the ```-j``` option on the very first compilation as it is not compatible with the linker. As the linker is not called again once all links are set, using this option on future compilations should be fine.

Compile the depends (scip, soplex, bliss and googletest) with

    make [args] deps

### Links to submodules
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

### Compile GCG
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

### Test GCG
You can run a short test with

    make [args] test


### Common Errors
On some distros, including the one used on RWTH cluster, the SCIP link does not work. Do this before compiling:

    cd lib && ln -s scip-git scip && cd ..
