# How to create and use your own decomposition {#own-dec-file}

> If **you know what the underlying structure of your problem** is, there are two common ways
> to let GCG know of it. Both of them are also elaborated in further detail in Use Case @ref u4. \n

## Defining a decomposition
Special file formats allow users to define own decompositions.
GCG can read and write the following file formats (among others):
- @subpage dec-file
- @subpage jdec-file

## Important Remarks
* a decomposition is rejected completely if there are any inconsistencies
* after reading (and and possibly assigning unassigned constraints because of consdefaultmaster, see above) implicit assignments are made 
  * unassigned constraints hitting at least two blocks -> assign to master
  * unassigned variables hitting at least two blocks -> assign to linking
* all constraints of an unassigned variable are master constraints -> variable is master variable;