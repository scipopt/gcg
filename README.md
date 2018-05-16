### README

GCG is an optimization tool based on [SCIP](http://scip.zib.de). This document
is a description of what is needed to get compile GCG.

### Official Homepage

The latest GCG release and the official documentation can be found under
http://www.or.rwth-aachen.de/gcg.

### Requirements

In order to compile GCG, you need

 * SCIP Optimization Suite >= 3.1.0

and optionally

 * bliss >= 0.72 ([Direct download](http://www.tcs.hut.fi/Software/bliss/bliss-0.72.zip))
 * readline-devel
 * zlib-devel
 * CPLEX >= 12.4.1

### Prerequisites

Install the required system libraries and tools for your Linux distribution:

**Ubuntu**

    sudo apt-get install build-essential libreadline-dev libz-dev

**OpenSuSE**

    sudo zypper in -t pattern devel_C_C++ && sudo zypper in readline-devel zlib-devel-static

**CentOs (RHEL, Fedora)**

    ?

### Installation with CMake

In order to create a build directory, type

    mkdir build

enter it via

    cd build

and configure it:

    cmake ..

After a successful configuration, you can compile it:

    make

Installation is done with

    make install

If a dependency such as SCIP is not found by `cmake`, you have to provide its
directory, by adding a `-DSCIP_DIR=/path/to/scip/build/or/install` argument
to the `cmake` call above.

### Installation with shipped Makefile

#### Creating softlinks

In order to create all necessary links, type

    make links

#### Compilation

You can compile GCG with

    make [options]

where the `[options]` are the same options as those from SCIP which is
documented on the official SCIP homepage under http://scip.zib.de. For further
options, consult the Makefile.

### Usage

The GCG binary can be found under `bin/gcg` (relative to where it was built).
An elaborate example can be found in the official documentation
[here](http://www.or.rwth-aachen.de/gcg/EXAMPLE.html).

### Creating the documentation

The official documentation can be found either under
http://www.or.rwth-aachen.de/gcg or it can be locally created using doxygen with

    make doc


### Licensing

GCG is licensed under the GNU Lesser General Public Licence (see [LICENSE](LICENSE)).
The instances in check/instances are distributed with a different license, see
[instance licenses](check/instances/readme) for more details on the license and
origin of these instances.
