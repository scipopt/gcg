# How to create and use your own decomposition {#own-dec-file}

[TOC]

> If **you know what the underlying structure of your problem** is, there are two common ways
> to let GCG know of it. Both of them are also elaborated in further detail in Use Case @ref u4. \n

# Defining your own Decomposition
This page will give a short introduction to the `.dec` file format by means of an example
and an extensive list of possible declarations. Note that if GCG has already found the
desired decomposition, you can also write it out using `write selected` (see @ref u3 for more details).
## Introduction
A standard file format used by GCG as well as other solvers such as DECOMP or DIP is the `.dec` file.
Inside it, you can define different things, most importantly, you can 
**fix variables/constraints to a block or the master problem** such that the problem can 
later be decomposed by means of a @ref pricing-process "Dantzig-Wolfe decomposition" as efficiently as possible. 
Partial definitions (i.e. only fixing some variables/constraints) are supported and will be 
@ref detection-process "prioritized in the detection loop". The syntax used inside `.dec` files will be 
explained in the following sections, first by means of an example and then in extensive form.

## Example
Do define a decomposition, you need a given problem instance, whose constraints and variables are referred to by name:
```
\Problem name: /opt/instances/lukas/problems/coloring/coloring.zpl
Minimize
 cost:  + y#1 + y#2
Subject to
 b_use_;1;:
  -100 y#1 + x#2#1 + x#1#1 <= 0
 b_use_;2;:
  -100 y#2 + x#2#1 + x#1#1 <= 0
 b_edge_;1;1;:
  + x#5#1 + x#1#1 <= 1
 b_edge_;2;1;:
  + x#5#2 + x#1#2 <= 1
 b_edge_;1;2;:
  + x#7#1 + x#1#1 <= 1
 b_edge_;2;2;:
  + x#7#2 + x#1#2 <= 1
 m_alco_;1;:
  + x#1#2 + x#1#1 >= 1
 m_alco_;2;:
  + x#2#2 + x#2#1 >= 1
```
Note: This graph coloring problem MIP example is shortened.
Using the names defined in the model, you can then assign constraints to blocks like that:
```
\ Decomposition file for coloring.zpl
PRESOLVED
0
NBLOCKS
2
BLOCK 1
b_use_;1;
b_edge_;1;1;
b_edge_;1;2;
BLOCK 2
b_use_;2;
b_edge_;2;1;
b_edge_;2;2;
MASTERCONSS
m_alco_;1;
m_alco_;2;
```
`PRESOLVED`, `NBLOCKS` and `BLOCK x` are **section key words**. The key words are not case sensitive, 
though commonly written in capital letters and after every keyword, each line provides one value.
In this example, the constraints that each node is colored in one color is used as master constraint 
and each node's constraints concerning edges represent a block, allowing GCG to exploit potential 
symmetry in the problem.

## Section Key Words
Section key words are:
| Section Key Word | Mandatory | Possible Values | Default | Description |
|:---------------------------------------------:|:-:|:------:|:---:|:---:|
| `CONSDEFAULTMASTER`                           |   | Binary | 1   | if set to 1 then (directly after file is read) each unassigned constraint is assigned to the master (needed for backward compatibility) |
| `PRESOLVED`                                   | X | Binary |     | if set to 0 (1) then the decomposition is considered for the unpresolved (presolved) problem |
| `NBLOCKS`                                     | X | N      |     | followed by number of (possibly empty) blocks this decomposition file has information for |
| `BLOCK` (or `BLOCKCONSS` or `BLOCKCONS`)      |   |        |     | followed by block index (starting with 1); each following line contains name of a constraint belonging to this block |
| `MASTERCONSS` (or `MASTERCONS`)               |   |        |     | each following line contains name of a constraint belonging to the master |
| `BLOCKVARS`                                   |   |        |     | followed by block index (starting with 1); each following line contains name of a variable belonging to this block |
| `MASTERVARS` (or `MASTERVAR`)                 |   |        |     | each following line contains name of a master variable; (belongs explicitly only to master constraints) |
| `LINKINGVARS` (or `LINKINGVAR`)               |   |        |     | each following line contains name of a linking variable |

## Important Remarks
* a decomposition is rejected completely if there are any inconsistencies
* after reading (and and possibly assigning unassigned constraints because of consdefaultmaster, see above) implicit assignments are made 
  * unassigned constraints hitting at least two blocks -> assign to master
  * unassigned variables hitting at least two blocks -> assign to linking
* all constraints of an unassigned variable are master constraints -> variable is master variable;