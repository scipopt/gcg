# Easy Installation {#easy-install}
> This page guides you through the **default installation of GCG and the SCIP Optimization Suite**. The steps below will only function for Linux and 
> Mac systems out of the box. At points, it might be required to alter the package management system (here, we use `apt`) according to your system.\n
> **Windows Users:** If you only have a Windows system, please read the @ref windows-install "Windows Installation Guide" for more information.

# Install GCG and the SCIP Optimization Suite {#install-scipopt}
> **If you want to let an automated installer do the work for you, you can download it [here](installSCIPopt.sh).**\n
> **Note:** The script requires you to have the SCIP `.tar` file (the source code) downloaded already 
> (obtainable [here](https://scipopt.org/index.php#download)). Also, it might be required to make the script
> executable, using e.g. `sudo chmod +x installSCIPopt.sh`.

@htmlonly
<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/7E30DoPyCFc" style="margin:auto; display:block" frameborder="3"  allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
@endhtmlonly

### Step 1: Prerequisites
Update your package lists:

    sudo apt-get update

Install the required system libraries:

    sudo apt-get install build-essential libreadline-dev libz-dev libgmp3-dev lib32ncurses5-dev libboost-program-options-dev libblas-dev libjansson-dev

To use all detectors, some more packages are required. See @ref install-optional for more information.

### Step 2: Get the SCIP Optimization Suite
Download the most recent SCIP Optimization Suite **Source Code** from the 
[SCIP Page](https://scipopt.org/index.php#download).
Unzip it and go into the folder, replacing `X.X.X` by the version you downloaded.

    tar xvzf scipoptsuite-X.X.X.tgz
    cd scipoptsuite-X.X.X

### Step 3: Compile via CMake
> The installation via CMake is **recommended for new users**. It does
> not offer all of the capabilities that Makefiles offer out-of-the-box.

First, make sure that cmake is installed correctly and up-to-date:

    sudo apt-get install cmake

Then create the build directory and prepare the build.

    cmake -S. -B./build -DIPOPT=OFF -DPAPILO=OFF

Arguments to be added to `cmake -S. -B./build` if needed can be found under @ref cmake-args. After the build files were prepared, you can start the compilation.

    cd build
    make gcg -j8

You can specify the number of compile jobs by changing the number appended to the `-j` flag if CPU/RAM usage is too high/low (using just `-j` will not limit the number of jobs).
If your compilation results in an error, please check our @ref install-trouble "troubleshooting section".\n
If your compilation succeeded, we recommend to test your installation by running a quick check (from within the 
build directory) on some instances:

    make gcg_check

Note that the available `make`-commands provided by CMake's Makefile differ from the commands provided by the Makefile located in the root directory.

### Step 3 (Alternative): Compile via Makefile
> The installation using the Makefiles build system is recommended for developers,
> since it offers more testing capabilities (see @ref conduct-experiments).

First, make sure that make is installed correctly and up-to-date:

    sudo apt-get install make

Then, from the root folder of the SCIP Optimization Suite, compile everything with

    make
    make gcg

Arguments to be added to `make gcg` if needed can be found under @ref makefiles-args.
If your compilation results in an error, please check our @ref install-trouble "troubleshooting section".\n
If your compilation succeeded, we recommend to test your installation by running a quick check (from within the 
`scipoptsuite-X.X.X/gcg/` directory) on some instances:

    make test

### Step 4: Use GCG!
See our @ref getting-started - Guide for further comments on the usage.

## Troubleshooting {#install-trouble}
In the following, we want to list common mistakes that occur 
during the installation and different solutions (ordered by 
number of successful applications) that have worked in the past.

**Problem**:
```
-- Configuring incomplete, errors occurred!
See also "/home/user/scipoptsuite-7.0.1/build/CMakeFiles/CMakeOutput.log".
See also "/home/user/scipoptsuite-7.0.1/build/CMakeFiles/CMakeError.log".
```
**Possible Solutions**:
- Try executing `cmake ..` again.
- Check if all required packages have been found. Install all
packages listed under `-- The following REQUIRED packages have not been found:`.
- Execute `cmake ..` as root user (`sudo cmake ..`).

\n
**Problem**:
```
/usr/bin/ld: ../../papilo/tbb_release/libtbb.a(task_group_context.o): in function `tbb::internal::basic_tls<unsigned long>::get()':
/home/user/scipoptsuite-7.0.1/papilo/external/tbb/./src/tbb/tls.h:43: undefined reference to `pthread_getspecific'
```

**Possible Solutions**:
- _No solution found yet._