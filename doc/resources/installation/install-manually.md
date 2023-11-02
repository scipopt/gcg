# Manual Installation {#install-manually}
> This page guides you through **installing GCG manually**. It is not necessary for most users.
> Please also check out the compiler and OS @ref compatibility and @ref known-bugs. 
> You can also find all arguments for the installation under @subpage makefiles-args and @subpage cmake-args.

# Install GCG manually
## Prerequisites
Installing SCIP, SoPlex, ZIMPL, hMETIS, Bliss, and GCG manually requires <a href="http://www.or.rwth-aachen.de/gcg">GCG</a>, <a href="https://scipopt.org/">SCIP</a>, <a href="http://soplex.zib.de/">SoPlex</a>, <a href="http://zimpl.zib.de/">ZIMPL</a>, 
and optionally <a href="http://www.tcs.hut.fi/Software/bliss/">Bliss</a>, <a href="https://users.aalto.fi/~pat/cliquer.html">Cliquer</a> and <a href="http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview">hMETIS</a>. In the following, we assume that all source code archives were saved within one folder, i.e.

    ls
    gcg-[version].tgz  scip-[version].tgz  soplex-[version].tgz  zimpl-[version].tgz

Furthermore, the installed packages are just as for the @ref easy-install.
## Main Installation
### ZIMPL and SoPlex
ZIMPL and SoPlex can be compiled by a simple `make`.

    tar xvfz zimpl-[version].tgz
    cd zimpl-[version]/
    make

    tar xvfz soplex-[version].tgz
    cd soplex-[version]/
    make

### SCIP
Compiling SCIP requires links to both ZIMPL and SoPlex.

    tar xvfz scip-[version].tgz
    cd scip-[version]/
    make
    
When asked for links to SoPlex, we need to set links from `scip-[version]/lib/spxinc` to `soplex-[version]/src/` and `scip-[version]/lib/libsoplex.*` to `soplex-[version]/lib/libsoplex.*`. For the linker, a link to either `libsoplex.*.a` or `libsoplex-*.so` is enough, so in our case we only need to specify the path to `libsoplex.*.a`:

    ../../../soplex-[version]/src/
    ../../../soplex-[version]/lib/libsoplex.linux.x86_64.gnu.opt.a

For ZIMPL, links need to point from `scip-[version]/lib/zimplinc/zimpl` to `zimpl-[version]/src/` and from `scip-[version]/lib/libzimpl.*` to `zimpl-[version]/lib/libzimpl.*`. Again, we may also ignore the link to `libzimpl.*.so`:

    ../../../zimpl-[version]/src/
    ../../../zimpl-[version]/lib/libzimpl.linux.x86_64.gnu.opt.a

### GCG
> When compiling GCG, all flags that GCG and SCIP have in common (e.g. `ZIMPL`) will be set to what you set during the compilation of SCIP.
> In general, if you set any flags at all, we recommend to **set them for _both_ SCIP and GCG**. 

Now, GCG can be compiled. The links to SCIP need to point from `gcg-[version]/lib/scip` to `scip-[version]/`.

    tar xvfz gcg-[version].tgz
    cd gcg-[version]/
    make
    ../../scip-[version]/

# Install Optional Packages {#install-optional}
Just as for the required packages, we assume that you have each optional package you want to install lying in
the directory above the GCG root directory for the instructions, i.e.

    ls
    bliss-[version].zip  cliquer-[version].tar.gz gcg/ hmetis-[version]-linux.tar.gz 

## Bliss
> Bliss is used for symmetry detection and required for the @ref det-isomorph. \n
> We recommend using the version shipped with SCIP. Bliss can be obtained
> [here](https://users.aalto.fi/~tjunttil/bliss/) or [here](https://github.com/ds4dm/Bliss) (SCIP's version).

Bliss is bundled with SCIP since version 8.0.1. Hence, no additional download is required anymore. Moreover, Bliss is enabled by default.

If everything went correctly, start GCG and you will be greeted by the line
```
External codes: 
    bliss 0.77           A Tool for Computing Automorphism Groups of Graphs by T. Junttila and P. Kaski (http://www.tcs.hut.fi/Software/bliss/)
```
and you are ready to use the bliss symmetry detector.

## Cliquer
> The Cliquer libary contains routines for clique searching which we use to @ref pricing-solvers "solve pricing problems".\n
> We recommend using version 1.21 or higher. The cliquer source code can be obtained
> [here](https://users.aalto.fi/~pat/cliquer.html).

First, **extract the Cliquer source code** (one folder level above the GCG root directory):

    tar xvfz cliquer-[version].tar.gz

Then, go into the GCG root directory and **link to it** and **compile it**:

    cd gcg/
    ln -s ../../cliquer-[version]/ lib/cliquer-git
    make cliquer

After that, **recompile** using the corresponding flag.

    make deps CLIQUER=true
    make CLIQUER=true
    ../cliquer-git
    ../cliquer-git/libcliquer.a

If everything went correctly, start GCG and you will be greeted by the line
```
External Codes:
    Cliquer              A set of C routines for finding cliques in an arbitrary weighted graph by S. Niskanen and P. Ostergard (https://users.aalto.fi/~pat/cliquer.html)
```
and you will be able to use the Cliquer pricing problem solver. 

## hMETIS
> GCG will only work with hMETIS versions of 2.0 and higher. hMETIS can be obtained
> [here](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/download). It will not be linked,
> but called through a system call.

In our code, hMETIS is used for hypergraph partitioning in three detectors (see @ref detectors).
If you want to compile GCG with hMETIS, you have to set a link to it.

    tar xvfz hmetis[version].tar.gz
    cd gcg-[version]/
    ln -s ../hmetis[version]/Linux-x86_64/hmetis[version] hmetis

Furthermore, depending on your system, it might be required to change the define `HMETIS_EXECUTABLE` 
inside the three detectors (`dec_h*partition.cpp`) source code files, e.g. from "hmetis" to "./hmetis"
and recompile. Alternatively, you can also do `export PATH="$PATH:./"` before calling GCG.

