# How to use this documentation {#MANUAL}
> **This guide aims to help you in getting started with GCG.** It is recommended
> to read the section "Where to begin?" (see below) completely to avoid unnecessary
> misunderstandings.

# Structure of this Documentation
The GCG documentation and manual is structured into two guides, a misc section
with everything that does not belong into these guides and the C-API
@htmlinclude structureofmanual.html

## Where to begin?
### Installation
In any case, you will have to start with the installation.
- For a new user, we suggest a complete @ref install "SCIP Optimization Suite Installation".
- For a developer or if you need the most recent dev version of GCG, you will have to @ref git-install "install GCG and SCIP using Git".
- <span style="color:grey !important;">If needed, GCG can also be installed @ref install-manually "manually".</span>

### Using and developing with GCG
The following sections are very short descriptions of the longer use cases that
we offer guides for.
\n\n

#### Using GCG as a solver from the *command line*
If you intend to use GCG only as a solver to get optimal solutions for
all your programs without any menu, you can just call GCG from the command line
and advise it to terminate after solving.\n
⇨ @ref u1
\n\n

#### Using GCG as a solver with the *interactive console*
If you intend to use GCG only as a solver to get optimal solutions for
all your programs, you can just @ref start-gcg "start GCG" and you
will be completely satistied with the @ref basic-commands "basic commands".\n
⇨ @ref u2
\n\n

#### Using GCG to gain new *knowledge about the structure* of your model (Explore Menu)
To get to know about how your GCG sees your instance and why, we recommend
to check out the @ref explore-menu "Explore Menu". With it, you can make GCG show the structure that its
@ref detectors "detectors" found. \n
⇨ @ref u3
\n\n

#### Using GCG to gain new *knowledge about how GCG ran* on your instance (Visualization Suite)
To see how GCG ran with your instance, for example how long it detected or how performant
the pricing was, you should check out the @ref visu "Visualization Suite". It illustrates not only
timings, but also algorithmic behaviour.\n
⇨ @ref u4
\n\n

#### Develop new plugins for GCG
If you are not satisfied with the time GCG needs to solve your instance, then
you should consider @ref howtoadd "adding your own plugins" that speed up the solving of your instance.\n
⇨ @ref u5
\n\n
