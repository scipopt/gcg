# Manual Installation {#install-manually}
> This guide aims to help you with installing GCG as a developer. It is not necessary for most users.
> Please also check out the compiler and OS @ref compatibility and @ref known-bugs. 
> You can also find all arguments for the installation under @subpage makefiles-args and @subpage cmake-args.


# Install GCG manually
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
