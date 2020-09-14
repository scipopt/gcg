# Developer's Guide {#devs}

The **first part** of this guide includes a comprehensive description of how things are implemented in SCIP and GCG
as well as the interaction between GCG and SCIP, to make the reader grasp what the differences are when implementing 
own code inside GCG.

@subpage dev-getting-started \n

\n
The **second part** describes what's already implemented in GCG and how it is and further can be used to make things much easier than with SCIP.

@subpage presolving \n
@subpage detection \n
@subpage pricing \n
@subpage cutting \n
@subpage howtouse \n

\n
The **third part** then finally takes you by the hand to implement your own plug-ins for
GCG. First, the howtos show you a rough outline of what's to be done, and then you can
consider some example application projects to get a first basis for your very own project in GCG.

@subpage howtoadd \n
@ref example-projects

\n
For **further reading** (and coding), we then also recommend the following pages as reference:

<a href="modules.html"><b>GCG C API</b></a> \n
@ref CHG \n
<a href="https://scip.zib.de/doc-6.0.2/html/CODE.php"><b>SCIP and GCG Coding Style Guidelines</b></a> \n
