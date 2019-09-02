# How to use this documentation {#MANUAL}
> **This guide aims to help you in getting started with GCG.** It is recommended
> to read the section "Where to begin?" (see below) completely to avoid unnecessary
> misunderstandings.

## Structure of this Documentation
The GCG documentation and manual is structured into two guides, a misc section
with everything that does not belong into these guides and the C-API
@htmlinclude structureofmanual.html

## Where to begin?
### Installation
In any case, you will have to start with the installation.
- For a new user, we suggest a complete @ref install "SCIP Optimization Suite Installation".
- For a developer or if you need the most recent dev version of GCG, you will have to @ref git-install "install GCG and SCIP using Git".

### Using and developing with GCG
The following sections are very short descriptions of the longer use cases that
we offer guides for.

#### Using GCG as a solver
If you intend to use GCG only as a solver to get optimal solutions for
all your programs, you can just @ref start-gcg "start GCG" and you
will be completely satistied with the @ref basic-commands "basic commands".\n
⇨ @ref u1

#### Using GCG to gain new knowledge about your instance
To get to know about how your GCG sees your instance and why, we recommend
to check out all the built-in features we created. Currently, the biggest ones
are the @ref explore-menu and the @ref visu .\n
⇨ @ref u2

#### Develop new Plugins for GCG
If you are not satisfied with the time GCG needs to solve your instance, then
you should consider implementing plugins that speed up the solving of your instance.\n
⇨ @ref u3
