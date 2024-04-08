# Git Installation {#git-install}
> This page guides you to through the **installation of GCG from Git**. You can use our private repository (if you have access to it) or the public repository on [GitHub](https://github.com/scipopt/gcg) (requires submodules to be cloned/set up manually).
> If you have a Windows computer, please see @subpage windows-install "Windows Installation Guide" for more information.

> **If you want to let an automated installer do the work for you, you can download it [here](installGCGfromGit.sh).**\n
> **Note:** This script uses our private repository. Your SSH key has to be registered correctly for it to work. Also, it might be required to make the script
> executable using `sudo chmod +x installGCGfromGit.sh`.

@htmlonly
<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/nsTiXuPW1WE" style="margin:auto; display:block" frameborder="3"  allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
@endhtmlonly

## Automated Installation
If you have access to our private repository, you can use the [automated installer](installGCGfromGit.sh). You can execute it like that:
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
libboost-program-options-dev
libjansson-dev
```

## Git prerequisites
> The following steps require access to our private repository. You can also use the public one and set up the submodules manually (the folders `lib/soplex-git` and `lib/scip-git`)
Please make sure that you have your ssh key saved in your Git account. Otherwise, you won't be able to
clone or commit anything, since you have not authorized yourself.
First, clone the git repo with

    git clone git@git.or.rwth-aachen.de:gcg/gcg.git [gcg-dir]

You **have** to clone it via ssh, otherwise the submodules won't work.
To initialize the SCIP, SoPlex, and the googletest framework, goto the repository folder "gcg" and run

    git submodule init
    git submodule sync
    git submodule update

## Manual Installation using Makefiles

### Compile SCIP and SOPLEX using makefiles
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

## Manual Installation using CMake
> Note that the available `make`-commands provided by CMake's Makefile differ from the commands provided by the Makefile located in the root directory.

### Build using provided CMake Presets
If your installed CMake version is equal to 3.25 or higher, you can configure and build GCG as well as test the build by calling

    cmake --workflow --preset gcg-linux-release

> You can replace `linux` by `macos` for MacOS builds (or even by `windows`, see @ref msvc).
> By using `--preset gcg-linux-debug` a debug build will be compiled and tested.

### Manual CMake Build
In case of an older CMake version you can also build GCG manually.

#### Build
The following commands will configure and build GCG

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGCG_DEV_BUILD=ON -DZIMPL=OFF -DIPOPT=OFF -DPAPILO=OFF
    cd build
    make gcg -j8

You can specify the number of compile jobs by changing the number appended to the `-j` flag if CPU/RAM usage is too high/low (using just `-j` will not limit the number of jobs).
Moreover, you can use `-DCMAKE_BUILD_TYPE=Debug` instead of `-DCMAKE_BUILD_TYPE=Release` to compile debug builds.

#### Test Build
You can run a short test with

    make gcg_check
