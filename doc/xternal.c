/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal.c
 * @brief  main document page
 * @author   Gerald Gamrath
 * @author   Martin Bergner
 * @author   Christian Puchert
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Generic Column Generation
 *
 * <b>What is GCG?</b>
 *
 * GCG is a generic branch-cut-and-price solver for mixed integer programs. It is based on the branch-and-cut-and-price
 * framework SCIP and is also part of the <a href="http://scip.zib.de">SCIP Optimization Suite</a>.
 *
 * After the standard presolving process of SCIP, GCG performs a Dantzig-Wolfe decomposition of the problem to obtain an
 * extended formulation of the problem. The decomposition is based on a structure either provided by the user or
 * automatically detected by one of the structure detectors included in GCG .
 *
 * During the solving process, GCG manages two SCIP instances, one holding the original problem, the other one
 * representing the reformulated problem. The original instance coordinates the solving process while the other one
 * builds the tree in the same way, transfers branching decisions and bound changes from the original problem and
 * solves the LP relaxation of the extended formulation via columns generation.
 *
 * GCG is developed jointly by <a href="http://www.or.rwth-aachen.de">RWTH Aachen</a> and <a
 * href="http://www.zib.de">Zuse-Institute Berlin</a>
 * and has more than 30,000 lines of C code.
 *
 *
 * <table cellpadding="0px" border="0" width="100%">
 *   <tr>
 *     <td nowrap >
 * <b>How to get started</b>
 * - \ref DOWNLOAD     "Download locations"
 * - \ref INSTALL      "Installation information"
 * - \ref EXAMPLE      "How to get started (example)"
 * - \ref FILEFORMATS  "Input file formats"
 *
 * <table cellpadding="0px" border="0" width="100%">
 *   <tr>
 *     <td nowrap >
 * <b>Further Information</b>
 * - \ref AUTHORS      "Current GCG developers"
 * - \ref CHANGELOG    "Changelog"
 * - \ref RELEASENOTES "Release notes"
 * - \ref LICENSE      "Licensing information"
 *
 * @version  1.0.0
 *
 * <b>Further Documentation</b>
 * - \ref IMPORTANTMETHODS "Important methods for writing GCG plugins"
 * - \ref PUBLICMETHODS "List of callable functions"
 * - \ref PRICINGSOLVER "How to write a custom pricing problem solver"
 * - \ref BRANCH "How to write a custom branching rule"
 * - \ref HEUR "How to write a custom heuristic"
 * - \ref DETECT "How to write a custom structure detector"
 * - \ref FAQ "Frequently asked questions"
 * - \ref PARAMS "Default parameter settings"
 * - \ref bug "Known bugs"
 *     </td>
 *   </tr>
 * </table>
 */

/**@page LICENSE Licensing Information
 *
 * GCG is released under the GNU Lesser General Public License:
 *
 * \verbinclude LICENSE
 */

/**@page EXAMPLE How to get started
 *
 *
 * If you want to use GCG as a solver for your problem, you'll find some starting information on this page
 *
 * You need a compiled GCG binary. The \ref INSTALL "Install section" will guide you through the correct steps. You'll further
 * need a problem you want to solve and you need to know what the structure information for the Dantzig-Wolfe reformulation does
 * look like for your problem.  GCG can read various file formats that SCIP can read, namely MPS, LP, CIP, and many more.
 *
 * If you want to download the complete source code of the GCG, we recommend downloading the complete SCIP optimization suite as
 * it will contain everything you will need in order to produce a working binary. The bin packing instance N1C2W2_O, which will
 * serve as an example in this tutorial, can be found under <code>gcg-[version]/check/instances/bpp/N1C2W2_O.BPP.lg</code>.
 *
 * Now start your binary, without any arguments. The usual place is <code>bin/gcg</code>. This opens the interactive shell, which should look somehow like this:
 *
 * \verbatim
GCG version 1.0.0 [GitHash: v100-0-g2e15c6e]
Copyright (c) 2010-2012 Operations Research, RWTH Aachen University
                        Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)

SCIP version 3.0.0 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 1.6.0.7] [GitHash: 8330fdf]
Copyright (c) 2002-2012 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)

External codes:
  Readline 6.2         GNU library for command line editing (gnu.org/s/readline)
  SoPlex 1.6.0.7       Linear Programming Solver developed at Zuse Institute Berlin (soplex.zib.de) [GitHash: 257a622]
  cppad-20120101.3     Algorithmic Differentiation of C++ algorithms developed by B. Bell (www.coin-or.org/CppAD)
  ZLIB 1.2.5           General purpose compression library by J. Gailly and M. Adler (zlib.net)

user parameter file <gcg.set> not found - using default parameters

GCG>

\endverbatim
 *
 * First of all "help" shows you a list of all available shell commands. Brackets indicate a submenu with further options.
 * \verbatim
GCG> help

 <display>             display information
 <set>                 load/save/change parameters
...
 read                  read a problem
\endverbatim
 *
 * Let's solve a problem with Dantzig-Wolfe reformulation... use "read <path/to/file>" to parse a problem file,
 * "read <path/to/file>" to provide structure information, "optimize" to solve it and "display
 * solution" to show the nonzero variables of the best found solution.

 * \verbatim
GCG> read check/instances/bpp/N1C2W2_O.BPP.lp

read problem <check/instances/bpp/N1C2W2_O.BPP.lp>
============

original problem has 2550 variables (0 bin, 2550 int, 0 impl, 0 cont) and 100 constraints
GCG> read check/instances/bpp/N1C2W2_O.BPP.dec

read problem <check/instances/bpp/N1C2W2_O.BPP.dec>
============

original problem has 2550 variables (0 bin, 2550 int, 0 impl, 0 cont) and 100 constraints

GCG> optimize

presolving:
presolving (0 rounds):
 0 deleted vars, 0 deleted constraints, 0 added constraints, 0 tightened bounds, 0 added holes, 0 changed sides, 0 changed coefficients
 0 implications, 0 cliques
presolved problem has 2550 variables (2550 bin, 0 int, 0 impl, 0 cont) and 100 constraints
    100 constraints of type <linear>
transformed objective value is always integral (scale: 1)
Presolving Time: 0.07
Chosen decomposition with 50 blocks of type arrowhead.


  time | node  | left  |LP iter|MLP iter|LP it/n| mem |mdpt |ovars|mvars|ocons|mcons|mcuts|confs|  dualbound   | primalbound  |  gap
   0.2s|     1 |     0 |     0 |      1 |     - |  16M|   0 |2550 |  10 | 101 | 101 |   0 |   0 | 0.000000e+00 |      --      |    Inf
*r 1.1s|     1 |     0 |     0 |    230 |     - |  17M|   0 |2550 | 500 | 101 | 101 |   0 |   0 | 0.000000e+00 | 5.000000e+01 |    Inf
Starting reduced cost pricing...
*r 1.3s|     1 |     0 |     0 |    255 |     - |  17M|   0 |2550 | 600 | 101 | 101 |   0 |   0 | 0.000000e+00 | 4.600000e+01 |    Inf
...
*r 7.4s|     1 |     0 |     0 |    799 |     - |  22M|   0 |2550 |3350 | 101 | 101 |   0 |   0 | 0.000000e+00 | 3.000000e+01 |    Inf
  time | node  | left  |LP iter|MLP iter|LP it/n| mem |mdpt |ovars|mvars|ocons|mcons|mcuts|confs|  dualbound   | primalbound  |  gap
*r10.2s|     1 |     0 |     0 |   1086 |     - |  24M|   0 |2550 |4600 | 101 | 101 |   0 |   0 | 0.000000e+00 | 3.000000e+01 |    Inf
  10.3s|     1 |     2 |     0 |   1086 |     - |  26M|   0 |2550 |4600 | 101 | 101 |   0 |   0 | 2.900000e+01 | 3.000000e+01 |   3.45%
  10.4s|     1 |     2 |     0 |   1086 |     - |  26M|   0 |2550 |4600 | 101 | 101 |   0 |   0 | 2.900000e+01 | 3.000000e+01 |   3.45%
R 10.4s|    10 |     0 |     0 |   1348 |  29.1 |  26M|   9 |2550 |4600 | 101 | 101 |   0 |   0 | 2.900000e+01 | 2.900000e+01 |   0.00%

SCIP Status        : problem is solved [optimal solution found]
Solving Time (sec) : 10.45
Solving Nodes      : 10
Primal Bound       : +2.90000000000000e+01 (4 solutions)
Dual Bound         : +2.90000000000000e+01
Gap                : 0.00 %

GCG> display solution

objective value:                                   29
x#1#10                                              1 	(obj:0)
x#2#49                                              1 	(obj:0)
x#3#33                                              1 	(obj:0)
...
y#1                                                 1 	(obj:1)
y#2                                                 1 	(obj:1)
...
y#50                                                1 	(obj:1)

GCG>
\endverbatim
 *
 * This tells us the following: After "optimize", GCG would first go into presolving. Since we have specified a structure
 * information for the original problem, GCG will currently disable presolving to not interfere with the decomposition.
 * Each round of presolving will be displayed in a single line, with a short summary at the end. Here, there has been
 * no round. Thus, it is not displayed and presolving is stopped. Afterwards, GCG will print out short information about
 * the currently used decomposition.
 *
 * Then, we see the actual solving process. The second output line and some of the
 * following lines indicate that new incumbent solutions were found by the primal heuristic with display character "r"
 * in the master (indicated by a '*' followed by the character); see, how the "primalbound" column changes from "--"
 * (no feasible solution found so far) to 50 and goes down to 30 in the following lines.
 * The third line indicates that reduced cost pricing is started. Before that, variables are added to the (initially empty)
 * master problem by so-called Farkas pricing to render it feasible. Reduced cost pricing now adds variables with negative
 * reduced costs to the master problem to improve its LP value.
 * After 10.4 seconds, the root node processing is finished. We needed 1086 master LP iterations (MLP iter) to solve the
 * LP relaxation of the master problem to an optimal value of 29, that is why the "dualbound" column changes from 0 to 29.
 * We see that there are now two open nodes in the "left" column. From now on, we will see an output line every hundredth
 * node or whenever a new incumbent is found (e.g. at node 10 in the above output). In our case, the "dualbound" at the
 * root node was optimal, this is why it is not changing anymore. At one point, both primal bound and dualbound will be the
 * same, and the solving process terminates, showing us some wrap-up information.
 *
 * The exact performance varies amongst different architectures, operating systems, and so on. Do not be worried if
 * your installation needs more or less time or nodes to solve.
 *
 * We might want to have some more information now. Which were the primal heuristics that found the solutions? What plugins
 *  were called during the solutions process and how much time did they spend? How did the instance that we were solving
 *  look like?  Information on certain plugin types (e.g., primal heuristics, branching rules, separators) we get by
 *  "display <plugin-type>", information on the solution process, we get by "display statistics", and "display problem"
 *  shows us the current instance.
 *
 * \verbatim
GCG> dis heur

 primal heuristic     c priority freq ofs  description
 ----------------     - -------- ---- ---  -----------
 oneopt               b   -20000    1   0  1-opt heuristic which tries to improve setting of single integer variables
...
 simplerounding       r        0    1   0  simple and fast LP rounding heuristic
 gcgsimplerounding    r        0    1   0  simple and fast LP rounding heuristic on original variables
...
GCG> display statistics

Master Program statistics:
SCIP Status        : solving was interrupted [node limit reached]
Total Time         :      10.17
...
Pricers            :   ExecTime  SetupTime      Calls       Vars
  problem variables:       0.01          -        106       2776
  gcg              :       9.84       0.00        142       4600
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found
...
  simplerounding   :       0.03       0.00         58          8
...
Original Program statistics:
SCIP Status        : problem is solved [optimal solution found]
Total Time         :      10.46
...
  gcgrounding      :       0.00       0.00         10          1
  xpcrossover      :       0.04       0.00          1          1
  xprins           :       0.03       0.00          1          1
  zeroobj          :       0.07       0.00          1          0
Relaxators         :       Time      Calls
  gcg              :      10.28         10
...
GCG>
\endverbatim
 *
 * We see statistics for two different problems: The Dantzig-Wolfe master problem and the original problem.
 * We see that rounding and shifting were the heuristics producing the solutions in the beginning. Rounding is called at
 * every node, shifting only at every tenth level of the tree. The statistics are quite comprehensive, thus, we just
 * explain a few lines here. We get information for all types of plugins and for the overall solving process. Besides
 * others, we see that in 10 calls of the gcgrounding heuristic, 1 solution was found and that the relaxator gcg got called
 * 10 times and took a total of 10.28 seconds to execute. Further on top, we can see that pricing produced 4600 variables
 * in 9.84 seconds.
 *
 * To solve a problem a second time, we have to read it and start the optimization process again. This time, we let GCG
 * discover the decomposition and display it.
 *
 * Lets first see what detectors are available:
 * \verbatim
GCG> dis dete
 detector             priority char  description
 --------------       -------- ----  -----------
 connected                   0    b  Detector for classical and block diagonal problems
\endverbatim
 *
 * We only have the "connected" detector available which claims to detect classical set partitioning master problems.
 * Let's see if that works:
 *
 * \verbatim
GCG> read check/instances/bpp/N1C2W2_O.BPP.lp

read problem <check/instances/bpp/N1C2W2_O.BPP.lp>
============

original problem has 2550 variables (0 bin, 2550 int, 0 impl, 0 cont) and 100 constraints
GCG> detect
Starting detection
Detecting purely block diagonal structure: not found.
Detecting set partitioning master structure: found 50 blocks.
Chosen decomposition with 50 blocks of type bordered.
Detection was successful.
\endverbatim
 *
 * The "connected" detector has found 50 blocks. Let's inspect the decomposition:
 *
 * \verbatim
GCG> display dec
PRESOLVED
0
NBLOCKS
50
BLOCK 1
b_Capacity_1
BLOCK 2
b_Capacity_2
BLOCK 3
b_Capacity_3
BLOCK 4
b_Capacity_4
BLOCK 5
b_Capacity_5
...
MASTERCONSS
m_Allocate_1
m_Allocate_2
m_Allocate_3
m_Allocate_4
m_Allocate_5
...
\endverbatim
 *
 * This tells us the following: The structure was detected from the unpresolved problem and contains 50 blocks. Next,
 * the constraints per block (b_Capacity_1 in block 1) are listed. Finally, all constraints in the master are listed.
 * This DEC format is described in \ref reader_dec.h .
 *
 * Now, we can start playing around with parameters. We have a bin packing example and we know that it can be solved
 * efficiently with discretization, so let us enable this by changing the parameter "relaxing/gcg/discretization".
 *
 * \verbatim
GCG> set
  <branching>           change parameters for branching rules
  ...
  <relaxing>            parameters for <relaxing>
  ...

GCG/set> relax

  <gcg>                 parameters for <gcg>

GCG/set/relaxing> gcg

  aggregation           should identical blocks be aggregated (only for discretization approach)? [TRUE]
  discretization        should discretization (TRUE) or convexification (FALSE) approach be used? [TRUE]
  dispinfos             should additional information about the blocks be displayed? [FALSE]
  enforceproper         should propagated bound changes in the original be enforced in the master (only proper vars)? [TRUE]
  freq                  frequency for calling relaxation handler <gcg> (-1: never, 0: only in root node) [1]
  priority              priority of relaxation handler <gcg> [-1]

GCG/set/relaxing/gcg> discretization
current value: FALSE, new value (TRUE/FALSE): true
relaxing/gcg/discretization = TRUE
GCG> o
  ...
  time | node  | left  |LP iter|MLP iter|LP it/n| mem |mdpt |ovars|mvars|ocons|mcons|mcuts|confs|  dualbound   | primalbound  |  gap
   0.2s|     1 |     0 |     0 |      1 |     - |  16M|   0 |2550 |  10 | 101 | 101 |   0 |   0 | 0.000000e+00 |      --      |    Inf
  time | node  | left  |LP iter|MLP iter|LP it/n| mem |mdpt |ovars|mvars|ocons|mcons|mcuts|confs|  dualbound   | primalbound  |  gap
*  0.8s|     1 |     0 |     0 |    230 |     - |  17M|   0 |2550 | 500 | 101 | 101 |   0 |   0 | 0.000000e+00 | 5.000000e+01 |    Inf
*P 0.8s|     1 |     0 |     0 |    230 |     - |  17M|   0 |2550 | 500 | 101 | 101 |   0 |   0 | 0.000000e+00 | 5.000000e+01 |    Inf
Starting reduced cost pricing...
*r 1.0s|     1 |     0 |     0 |    255 |     - |  17M|   0 |2550 | 600 | 101 | 101 |   0 |   0 | 0.000000e+00 | 4.600000e+01 |    Inf
...
*F 7.4s|     1 |     0 |     0 |   1212 |     - |  24M|   0 |2550 |4600 | 101 | 101 |   0 |   0 | 2.900000e+01 | 2.900000e+01 |   0.00%
   7.4s|     1 |     0 |     0 |   1212 |     - |  24M|   0 |2550 |4600 | 101 | 101 |   0 |   0 | 2.900000e+01 | 2.900000e+01 |   0.00%

SCIP Status        : problem is solved [optimal solution found]
Solving Time (sec) : 7.44
Solving Nodes      : 1
Primal Bound       : +2.90000000000000e+01 (2 solutions)
Dual Bound         : +2.90000000000000e+01
Gap                : 0.00 %
\endverbatim
 *
 * We can navigate through the menus step-by-step and get a list of available options and submenus. Thus, we select
 * "set" to change settings, "relax" to change settings of relaxators, "gcg" for that particular
 * relaxator. Then we see a list of parameters (and sometimes yet another submenu for advanced parameters), and set
 * discretization to TRUE. If we already know the path to a certain setting, we can directly
 * type it (<code>set relaxing gcg discretization TRUE</code>). Note that we do not have to use the full names, but we
 * may use short versions, as long as they are unique.
 *
 * GCG can also write information to files. E.g., we could store the incumbent solution to a file, or output the
 * problem instance in another file format (the LP format is much more human readable than the MPS format, for example).
 * We can also write out the currently used decomposition by saving the problem as a decomposition format (DEC, BLK or REF).
 *
 * \verbatim
GCG> write sol N1C2W2_O.BBP.sol

written solution information to file <N1C2W2_O.BBP.sol>

GCG> write prob "N1C2W2_O.BBP.dec"
written original problem to file <N1C2W2_O.BBP.dec>

GCG> q
...
\endverbatim
 *
 * We hope this tutorial gave you an overview of what is possible using the SCIP interactive shell within GCG. Please also read our
 * \ref FAQ.
 *
 */

/**@page FAQ Frequently Asked Questions (FAQ)
 * \htmlinclude faq.inc
 */

/**@page DOWNLOAD Download Locations
 * @section Downloading GCG
 *
 * GCG can be downloaded from two locations
 * - Standalone from the server at <a href="http://or.rwth-aachen.de/gcg/download"> RWTH Aachen</a>
 * - Complete with SCIP in the SCIPoptSuite from <a href="http://scipoptsuite.zib.de/download">Zuse Institute Berlin</a>
 *
 * Instructions how to compile and build GCG can be found in the \ref INSTALL "Install section".
 *
 * We will not offer precompiled binaries and GCG may not compile on Microsoft Windows. It is developed and tested on GNU/Linux.
 */

/**@page INSTALL Installation information
 * \verbinclude INSTALL
 */

/**@page CHANGELOG Changelog
 *
 * \verbinclude CHANGELOG
 */

/**@page RELEASENOTES Release notes
 *
 * \verbinclude GCG-release-notes-1.0
 *
 */

/**@page AUTHORS GCG Authors
 * \htmlinclude authors.inc
 */


/**@page PARAMS GCG default parameter settings
 *
 * \verbinclude inc/parameters.set
 */


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page FILEFORMATS Input file formats supported by GCG.
 *
 * GCG supports all file formats supported by SCIP to read in problems, solution, etc.
 * E.g., the original problem can be read in as an .lp, .mps, or .cip file.
 *
 * If GCG is not able to automatically detect a structure suitable to perform a Dantzig-Wolfe reformulation, you need
 * to specify the structure yourself and make it available to GCG. There are some file formats for structure information.
 * If in doubt, we recommend the <code>dec</code> format which is documented in \ref reader_dec.h .
 *
 * You can find examples in the <code>check/instances/</code> directory.
 */

/**@defgroup BRANCHINGRULES Branching Rules
 * @brief This page contains a list of all branching rule which are currently available.
 *
 * A detailed description what a branching rule does and how to add a branching rule to GCG can be found
 * \ref BRANCH "here".
 */

/**@defgroup CONSHDLRS  Constraint Handler
 * @brief This page contains a list of all constraint handlers which are currently available.
 *
 * A detailed description what a constraint handler does and how to add a constraint handler to \SCIP can be found
 * in the SCIP documentation.
 */

/**@defgroup DETECTORS Detectors
 * @brief This page contains a list of all detectors which are currently available.
 *
 * A detailed description what a detector does and how to add a detector to GCG can be found
 * \ref DETECT "here".
 */

/**@defgroup DIALOGS Dialogs
 * @brief This page contains a list of all dialogs which are currently available.
 *
 * A detailed description what a dialog does and how to add a dialog to \SCIP can be found
 * n the SCIP documentation.
 */

/**@defgroup DISPLAYS Displays
 * @brief This page contains a list of all displays (output columns)  which are currently available.
 *
 * A detailed description what a display does and how to add a display to \SCIP can be found
 * in the SCIP documentation.
 *
 */

/**@defgroup FILEREADERS File Readers
 * @brief This page contains a list of all file readers which are currently available.
 *
 * A detailed description what a file reader does and how to add a file reader to \SCIP can be found
 * in the SCIP documentation.
 */

/**@defgroup NODESELECTORS Node Selectors
 * @brief This page contains a list of all node selectors which are currently available.
 *
 * A detailed description what a node selector does and how to add a node selector to \SCIP can be found
 * in the SCIP documentation.
 */

/**@defgroup PRICERS Pricers
 * @brief This page contains a list of all pricers which are currently available.
 *
 * Per default there exist no variable pricer. A detailed description what a variable pricer does and how to add a
 * variable pricer to \SCIP can be found in the SCIP documentation.
 */

/**@defgroup PRICINGSOLVERS Pricing solvers
 * @brief This page contains a list of all pricing solvers which are currently available.
 *
 * A detailed description what a pricing solver does and how to add a pricing solver to GCG can be found
 * \ref SOLVER "here".
 */

/**@defgroup PRIMALHEURISTICS Primal Heuristics
 * @brief This page contains a list of all primal heuristics which are currently available.
 *
 * A detailed description what a primal heuristic does and how to add a primal heuristic to \SCIP can be found
 * \ref HEUR "here".
 */

/**@defgroup PUBLICMETHODS Public Methods
 *
 * This page lists headers containing branch-and-price specific public methods provided by GCG.
 *
 * All of the headers listed below include functions that are allowed to be called by external users. Besides those
 * functions it is also valid to call methods that are listed in one of the headers of the (default) GCG plugins; in
 * particular, this holds for relax_gcg.h and pricer_gcg.h.
 *
 */

/**@defgroup RELAXATORS Relaxators
 * @brief This page contains a list of all relaxators which are currently available.
 */

/**@defgroup SEPARATORS Separators
 * @brief This page contains a list of all separators  which are currently available.
 *
 * A detailed description what a separator does and how to add a separator to \SCIP can be found
 * in the SCIP documentation.
 */

/**@defgroup TYPEDEFINITIONS Type Definitions
 * This page lists headers which contain type definitions of callback methods.
 *
 * All headers below include the descriptions of callback methods of
 * certain plugins. For more detail see the corresponding header.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page HEUR How to add primal heuristics
 *
 * For general information on how to add your own primal heuristics to GCG, first check the SCIP documentation.
 * However, one has to take into account some peculiarities when implementing heuristics that are included
 * in the original SCIP instance, i.e. work on the original variables.
 *
 * @section HEUR_LPSOLUTIONS Access to LP feasible solutions (on the original variables)
 *
 * Many MIP heuristics make use of an LP feasible solution. In SCIP, such a solution is obtained by solving the LP relaxation.
 * In GCG however, no LP relaxation is solved by default. A linearly feasible solution on the original variables comes from the
 * GCG relaxator plugin; it is a solution of the master LP that has been translated back into the original variables. To access
 * it, one should use the method GCGrelaxGetCurrentOrigSol() in relax_gcg.h.
 * Its fractional variables can be accessed by SCIPgetExternBranchCands() (rather than SCIPgetLPBranchCands() which is used
 * in SCIP heuristics).
 * \n
 * Note also that heuristics using LP solutions should use another timing than SCIP heuristics. Heuristics that are called after
 * solving a node's relaxation typically have the timing SCIP_HEURTIMING_AFTERLPNODE.
 * By default, no LPs are solved on the original problem. A heuristic relying on a linearly feasible solution should therefore
 * have the timing SCIP_HEURTIMING_AFTERNODE to ensure that the heuristic is called at all. One then must ensure that the node's
 * relaxation has indeed been solved to optimality and that the relaxation solution is valid. This can be done by placing
 * \code
 *  if( !SCIPisRelaxSolValid(scip) )
 *     return SCIP_OKAY;
 *
 *  if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) > SCIP_STAGE_SOLVING || SCIPgetLPSolstat(GCGrelaxGetMasterprob(scip)) != SCIP_LPSOLSTAT_OPTIMAL )
 *     return SCIP_OKAY;
 * \endcode
 * at the beginning of the HEUREXEC callback.
 *
 * @section HEUR_DIVING Diving on original variables
 *
 * A common class of heuristics are diving heuristics; they solve LPs with modified bounds to perform a depth-first
 * search on the Branch-and-bound tree. For this purpose, a probing mode and a diving mode have been implemented in SCIP,
 * which can be invoked by SCIPstartProbing(scip) and SCIPstartDive(scip), respectively. In this mode, temporary bound changes
 * on variables can be made, and modified LPs can be solved.
 * \n
 * In GCG, a special probing mode has been implemented for the original instance. This mode serves for performing changes on
 * the original instance but using the master LP instead of the original LP. It is invoked by GCGrelaxStartProbing() and terminated
 * by GCGrelaxEndProbing() and features the methods GCGrelaxPerformProbing() and GCGrelaxPerformProbingWithPricing(), which will
 * propagate any bound changes on the original instance to the extended instance and solve the resulting modified master LP, either
 * without or with pricing new variables in. See e.g. heur_gcgcoefdiving.c for an example on how to use them.
 *
 * @section HEURCOPY The HEURCOPY callback
 *
 * The HEURCOPY callback is executed when a SCIP instance is copied, e.g. to
 * solve a sub-SCIP. By
 * defining this callback as
 * <code>NULL</code> the user disables the execution of the specified
 * heuristic for all copied SCIP instances. This may deteriorate the performance
 * of primal heuristics using sub-SCIPs.
 * \n
 * For heuristics that are included in the original instance and make use of the extended instance as well (in
 * particular, most of the heur_gcg* and heur_xp* plugins), this callback should be set to null. This is because
 * sub-SCIPs are solved by SCIP rather than GCG and therefore do not know any master problem; including a GCG
 * specific heuristic into them would cause errors.
 *
 */

 /*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DETECT How to add structure detectors
 *
 * Structure detectors are used to detect or enforce a structure suitable for Dantzig-Wolfe Reformulation (DWR).
 * \n
 * A complete list of all detectors/enforcers contained in this release can be found \ref DETECTORS "here".
 *
 * In the following, we explain how the user can add its own structure enforcement plugin.
 * Take the basic detector (dec_connected.c) as an example.
 * As all other default plugins, it is written in C. There is currently no C++ wrapper available.
 *
 * Additional documentation for the callback methods of structure detectors, in particular for their input parameters,
 * can be found in the file type_detector.h.
 *
 * Here is what you have to do to implement a detector:
 * -# Copy the template files src/dec_xyz.c and src/dec_xyz.h into files named "dec_mydetector.c"
 *    and "dec_mydetector.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "mydetector".
 * -# Adjust the properties of the detector (see \ref DEC_PROPERTIES).
 * -# Define the detector data (see \ref DEC_DATA). This is optional.
 * -# Implement the interface methods (see \ref DEC_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref DEC_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref DEC_ADDITIONALCALLBACKS). This is optional.
 *
 *
 * @section DEC_PROPERTIES Properties of a Detector
 *
 * At the top of the new file "dec_mydetector.c", you can find the detector properties.
 * These are given as compiler defines.
 * The properties you have to set have the following meaning:
 *
 * \par DEC_DETECTORNAME: the name of the detector.
 * This name is used in the interactive shell to address the detector.
 * Additionally, if you are searching for a detector with SCIPfindDetector(), this name is looked up.
 * Names have to be unique: no two detectors may have the same name.
 *
 * \par DEC_DESC: the description of the detector.
 * This string is printed as description of the detector in the interactive shell.
 *
 * \par DEC_PRIORITY: the priority of the detector.
 * At the begin of the solving process, the detectors are called in a predefined order, which is given by the priorities
 * of the detectors.
 * The detectors are called in the order of decreasing priority.
 * \n
 * The priority of the detector should be set according to the complexity of the detection algorithm and the quality of the decomposition:
 * detectors that provide fast algorithms that usually have a good decomposition (i.e., provide good dual bound) should have a high
 * priority. An easy way to list the priorities of all detectors is "display detectors" in the interactive shell of GCG.
 *
 * \par DEC_DECCHAR: Display character of the detector.
 * The unique display character for the detector. It can be used to quickly distinguish similar structures by the detector which found them.
 *
 * \par DEC_ENABLED: Flag to indicate whether the detector should be enabled by default.
 * Disabled detectors are not started.
 *
 * @section DEC_DATA Detector Data
 *
 * Below the header "Data structures" you can find the struct "struct DEC_DetectorData".
 * In this data structure, you can store the data of your detector. For example, you should store the adjustable parameters
 * of the detector in this data structure.
 * \n
 * Defining detector data is optional. You can leave this struct empty.
 *
 *
 * @section DEC_INTERFACE Interface Methods
 *
 * At the bottom of "dec_mydetector.c", you can find the interface method SCIPincludeDetectionMyDetector(),
 * which also appears in "dec_mydetector.h".
 * \n
 * This method has to be adjusted only slightly.
 * It is responsible for notifying GCG (and especially cons_decomp.c) of the presence of the detector by calling the method
 * DECincludeDetector().
 * SCIPincludeDetectionMyDetector() is called by the user, if he wants to include the detector,
 * i.e., if he wants to use the detector in his application.
 *
 * If you are using detector data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_DetectorData afterwards. For freeing the
 * detector data, see \ref DETECTOREXIT.
 *
 * You may also add user parameters for your detector, see the parameters documentation of \SCIP for how to add user parameters and
 * the method SCIPincludeDetectionBorderheur() in dec_connected.c for an example.
 *
 *
 * @section DEC_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Detector
 *
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain
 * an operational algorithm. Detector plugins have only one fundamental callback method, namely the DETECTSTRUCTURE method.
 * This method has to be implemented for every detector; the other callback methods are optional.
 *
 * Additional documentation to the callback methods, in particular to their input parameters,
 * can be found in type_detector.h.
 *
 * @subsection DETECTSTRUCTURE
 *
 * The DETECTSTRUCTURE callback is called during the detection loop and should perform the actual detection.
 * It should inspect the problem instance at hand and deduct some structure from the constraint matrix.
 * It needs to store the structure information in DEC_DECOMP and needs to allocate the array where to store the
 * information.
 *
 * Typical methods called by a detector are, for example, SCIPgetVars(), SCIPGetConss(), DECdecompSetNBlocks(),
 * DECdecompSet*(), etc. .
 *
 * @section DEC_ADDITIONALCALLBACKS Additional Callback Methods of a Detector
 *
 * @subsection DETECTORINIT
 *
 * The INITDETECTOR callback is executed after the problem was transformed.
 * The detector may, e.g., use this call to initialize his detector data.
 * The difference between the original and the transformed problem is explained in
 * "What is this thing with the original and the transformed problem about?" on \ref FAQ.
 *
 * @subsection DETECTOREXIT
 *
 * If you are using detection data (see \ref DEC_DATA and \ref DEC_INTERFACE), you have to implement this method in order to free the detection data.
 * This can be done by the following procedure:
 * \code
 * static
 * DEC_DECL_EXITDETECTOR(decExitMydetector)
 * {
 *    DEC_DETECTORDATA* detectordata;
 *
 *    detectordata = DECdetectorGetData(detector);
 *    assert(detectordata != NULL);
 *
 *    SCIPfreeMemory(scip, &detectordata);
 *
 *    DECdetectorSetData(detector, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you have allocated memory for fields in your detector data, remember to free this memory
 * before freeing the detector data itself.
 * The DETECTOREXIT callback is executed before the solution process is started.
 * In this method, the detector should free all resources that have been allocated for the detection process in \ref DETECTORINIT.
 *
 */

 /*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page PRICINGSOLVER How to add pricing problem solvers
 *
 * Pricing problem solvers are called by the pricer to solve a single pricing problem either heuristically or to optimality
 * and return a set of solutions.
 * \n
 * A complete list of all pricing problem solvers contained in this release can be found \ref PRICINGSOLVERS "here".
 *
 * In the following, we explain how the user can add his own pricing problem solver.
 * Take the generic MIP pricing problem solver (src/solver_mip.c) as an example.
 * As all other default plugins, it is written in C. There is currently no C++ wrapper available.
 *
 * Additional documentation for the callback methods of pricing problem solvers, in particular for their input parameters,
 * can be found in the file type_solver.h.
 *
 * Here is what you have to do to implement a pricing problem solver:
 * -# Copy the template files src/solver_xyz.c and src/solver_xyz.h into files named "solver_mysolver.c"
 *    and "solver_mysolver.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "mysolver".
 * -# Adjust the properties of the solver (see \ref SOLVER_PROPERTIES).
 * -# Define the solver data (see \ref SOLVER_DATA). This is optional.
 * -# Implement the interface methods (see \ref SOLVER_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref SOLVER_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref SOLVER_ADDITIONALCALLBACKS). This is optional.
 *
 *
 * @section SOLVER_PROPERTIES Properties of a pricing problem solver
 *
 * At the top of the new file "solver_mysolver.c", you can find the solver properties.
 * These are given as compiler defines.
 * The properties you have to set have the following meaning:
 *
 * \par SOLVER_NAME: the name of the solver.
 * This name is used in the interactive shell to address the solver.
 * Names have to be unique: no two solvers may have the same name.
 *
 * \par SOLVER_DESC: the description of the solver.
 * This string is printed as description of the solver in the interactive shell.
 *
 * \par SOLVER_PRIORITY: the priority of the solver.
 * Whenever a pricing problem has to be solver, the solvers are called in a predefined order, until one of the solvers solved the
 * pricing problem. This order is given by the priorities of the solvers; they are called in the order of decreasing priority.
 * \n
 * The priority of the solver should be set according to the range of problems it can solve. The more specialized a solver is,
 * the higher its priority should be in order to be called before more general solvers that also cover this type of problem.
 * The default pricing solver that solves the pricing problem as a MIP using SCIP has priority 0, so all other pricing solvers
 * should have positive priority.
 * An easy way to list the priorities and descriptions of all solvers is "display solvers" in the interactive shell of GCG.
 *
 * @section SOLVER_DATA Solver Data
 *
 * Below the header "Data structures" you can find the struct "struct GCG_SolverData".
 * In this data structure, you can store the data of your solver. For example, you should store the adjustable parameters
 * of the solver in this data structure and also arrays to store the solutions that you return after solving a pricing problem.
 * \n
 * Defining solver data is optional. You can leave this struct empty.
 *
 *
 * @section SOLVER_INTERFACE Interface Methods
 *
 * At the bottom of "solver_mysolver.c", you can find the interface method GCGincludeSolverMysolver(),
 * which also appears in "solver_mysolver.h".
 * \n
 * This method has to be adjusted only slightly.
 * It is responsible for notifying GCG (and especially the pricer in the master SCIP instance) of the presence of the solver by
 * calling the method GCGpricerIncludeSolver().
 * SCIPincludeSolverMysolver() is should be called by the user, if he wants to include the solver,
 * i.e., if he wants to use the solver in his application. This can done for example by adding
 * \code SCIP_CALL( GCGincludeSolverMysolver(scip) );\endcode to masterplugins.c
 *
 * If you are using solver data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &solverdata) );
 * \endcode
 * You can also initialize the fields in struct SCIP_SolverData afterwards. For freeing the
 * solver data, see \ref SOLVEREXIT. Alternatively, you can also initialize and free the fields in the solver data
 * in the SOLVERINITSOL and SOLVEREXITSOL callback methods, respectively. For an example see, solver_mip.c.
 *
 * You may also add user parameters for your solver, see the SCIP documentation for how to add user parameters and
 * the method SCIPincludeSolverMip() in src/solver_mip.c for an example.
 *
 *
 * @section SOLVER_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Solver
 *
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain
 * an operational algorithm. Pricing problem solvers do not have fundamental callback methods, but they should
 * implement at least of the SOLVERSOLVE and SOVLERSOLVEHEUR methods.
 *
 * Additional documentation to the callback methods, in particular to their input parameters,
 * can be found in type_solver.h.
 *
 * @section SOLVER_ADDITIONALCALLBACKS Additional Callback Methods of a Solver
 *
 * @subsection SOLVERSOLVE
 *
 * The SOLVERSOLVE callback is called when the variable pricer in the master SCIP instance wants to solve a pricing problem
 * to optimality.
 * It is given a SCIP instance representing the pricing problem that should be solved and should check whether the pricing
 * problem has a structure that can be handled by the solver. If so, it should solve the pricing problem to optimality
 * and return the optimal objective value of the problem (to be stored in the lowerbound pointer) along with a set of
 * primal feasible solutions (including the optimal one).
 *
 * The solutions should be returned by setting some given pointers:
 * solvars should store the number of solutions returned, nsolvars should point to an array storing for each solution
 * the number of variables with non-zero value, and solisray should point to an array storing for each solution whether it
 * represents an extreme point of the pricing problem or an extreme ray. Furthermore, solvars and solvals should point to
 * arrays storing for each solution an array of variables and corresponding solution values, respectively (variables with
 * value 0 can be omitted). Therefore, a pricing problem solver should manage these arrays in its own data structure,
 * fill them after solving a pricing problem and set the return pointers to these arrays. Have a look at solver_mip for
 * an example of how this can be done.
 *
 * Last, the callback should adjust the given result pointer to SCIP_STATUS_OPTIMAL if the problem was solved to optimality,
 * to SCIP_STATUS_UNBOUNDED if the problem was solved and is unbounded, or SCIP_STATUS_UNKNOWN if the solver was not applicable
 * to the pricing problem or if the solving was stopped.
 *
 * @subsection SOLVERSOLVEHEUR
 *
 * The SOLVERSOLVEHEUR callback is called during heuristical pricing when the variable pricer in the master SCIP instance
 * wants to solve a pricing problem heuristically.
 * It has the same input and return parameters as the SOLVERSOLVE callback. It does not need to solve the pricing problem to
 * optimality, but should try to construct good feasible solutions using fast methods. Nevertheless, it can return a lower bound
 * for the optimal solution value of the problem, if it computes one.
 *
 * @subsection SOLVERFREE
 *
 * If you are using solver data, you have to implement this method in order to free the solver data.
 * This can be done by the following procedure:
 * \code
 * static
 * GCG_DECL_SOLVERFREE(solverFreeMip)
 * {
 *    GCG_SOLVERDATA* solverdata;
 *
 *    assert(scip != NULL);
 *    assert(solver != NULL);
 *
 *    solverdata = GCGpricerGetSolverdata(scip, solver);
 *    assert(solverdata != NULL);
 *
 *    SCIPfreeMemory(scip, &solverdata);
 *
 *    GCGpricerSetSolverdata(scip, solver, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you have allocated memory for fields in your solver data, remember to free this memory
 * before freeing the pricer data itself.
 *
 * @subsection SOLVERINIT
 *
 * The SOLVERINIT callback is executed after the problem is transformed.
 * The pricing problem solver may, e.g., use this call to replace the original constraints stored in its solver data by transformed
 * constraints, or to initialize other elements of his solver data.
 *
 * @subsection SOLVEREXIT
 *
 * The SOLVEREXIT callback is executed before the transformed problem is freed.
 * In this method the pricing problem solver should free all resources that have been allocated for the solving process in SOLVERINIT.
 *
 * @subsection SOLVERINITSOL
 *
 * The SOLVERINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to begin.
 * The pricing problem solver may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection SOLVEREXITSOL
 *
 * The SOLVEREXITSOL callback is executed before the branch-and-bound process is freed.
 * The pricing problem solver should use this call to clean up its branch-and-bound data, which was allocated in SOLVERINITSOL.
 *
 */

 /*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page BRANCH How to add branching rules to GCG
 *
 * Branching in the branch-cut-and-price context is a bit more complicated than in a branch-and-cut algorithm.
 * GCG manages two SCIP instances, one for the original instance and one for the reformulated instance, that are solved jointly.
 * The original SCIP instance coordinates the solving process, while the reformulated instance builds the tree in the same
 * way, with each node representing the Dantzig-Wolfe reformulated subproblem of the related node in the original instance.
 *
 * Therefore, if you want to implement a new branching rule, it has to be added to the original SCIP instance, but it gets some
 * more callback methods, e.g., to transfer the branching decision to the corresponding nodes of the master problem and
 * enforce these changes during propagation.
 * \n
 *
 * The linking of the nodes of the two trees is done via two constraint handler, cons_origbranch and cons_masterbranch,
 * that add local constraints to the nodes whoch know about the corresponding constraint in the other SCIP instance and by
 * this also the corresponding node. Therefore, each branchrule in GCG has to add one of the origbranch constraints to each
 * child node it creates. This origbranch constraint also stores a branching data, as described \ref BRANCHDATA "below", that
 * can be used to store information about the branching decisions.
 *
 * \code
 *  SCIP_NODE* childleft;
 *  SCIP_NODE* childright;
 *  GCG_BRANCHDATA* branchleftdata;
 *  GCG_BRANCHDATA* branchrightdata;
 *  SCIP_CONS* origbranchleft;
 *  SCIP_CONS* origbranchright;
 *
 *  /* create the branching data for the children */
 *  SCIP_CALL( SCIPallocMemory(scip, &(branchleftdata)) );
 *  SCIP_CALL( SCIPallocMemory(scip, &(branchrightdata)) );
 *
 *  /* add branching decision to the subproblems and store them in the branching data */
 *  (...)
 *
 *  /* create the origbranch constraints */
 *  SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchleft, "left", childleft,
 *        GCGconsOrigbranchGetActiveCons(scip), branchrule, branchleftdata) );
 *  SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchright, "right", childright,
 *        GCGconsOrigbranchGetActiveCons(scip), branchrule, branchrightdata) );
 *
 *  /* add constraints to nodes */
 *  SCIP_CALL( SCIPaddConsNode(scip, childleft, origbranchleft, NULL) );
 *  SCIP_CALL( SCIPaddConsNode(scip, childright, origbranchright, NULL) );
 *
 *  /* release constraints */
 *  SCIP_CALL( SCIPreleaseCons(scip, &origbranchleft) );
 *  SCIP_CALL( SCIPreleaseCons(scip, &origbranchright) );
 * \endcode
 *
 * A complete list of all branching rules contained in this release can be found \ref BRANCHINGRULES "here".
 *
 * In the following, we explain how the user can add his own pricing problem solver.
 * Take the branching rule that branches on original variables (src/branch_orig.c) as an example.
 * As all other default plugins, it is written in C. There is currently no C++ wrapper available.
 *
 * In the following, we will only explain the additional callback methods which a branching rule in GCG can implement.
 * For a general introduction about how to write a branching rule in SCIP and a description of the default callback methods
 * of a branching rule in SCIP, we refer to the SCIP documentation and the "How to add constraint handlers" there.
 *
 * Additional documentation for the GCG-specific callback methods of branching rules, in particular for their input parameters,
 * can be found in the file type_branchgcg.h.
 *
 * @section BRANCHDATA GCG specific branching data
 *
 * Below the header "Data structures" you can find the new struct "struct GCG_BranchData".
 * In this data structure, each branching rule can store information about the branching decisions applied to a newly created node.
 * Later, this information is used to transfer the branching decisions to the corresponding nodes in the reformulated SCIP instance.
 * \n
 * Defining branching data is optional. You can leave this struct empty.
 *
 * @section BRANCH_CALLBACKS GCG-specific callbacks of branching rules
 *
 * @subsection BRANCHACTIVEMASTER
 *
 * The BRANCHACTIVEMASTER callback is called whenever a node in the master problem is activated.
 * It should transfer the branching decisions that were made at the corresponding node in the original SCIP instance to the
 * current node of the reformulated SCIP instance. Therefore, it gets branching data stored at the corresponding node in the
 * original SCIP instance as a parameter. It should then for example reformulate a banching constraint to the
 * reformulated problem space and add it locally to the current node or update the pricing problem according to the
 * branching decision.
 *
 * @subsection BRANCHDEACTIVEMASTER
 *
 * The BRANCHDEACTIVEMASTER callback is called whenever a node in the master problem is deactivated.
 * It should undo all changes that were performed in the BRANCHACTIVEMASTER callback, e.g., remove branching restrictions
 * from the pricing problems.
 *
 * @subsection BRANCHPROPMASTER
 *
 * The BRANCHPROPMASTER callback is called whenever a node in the master problem is propagated by the masterbranch
 * constraint handler.
 * It should perform propagation according to the branching changes done at the related original node, which are stored in
 * the given banching data. In particular, it should fix all variable locally to 0, that contradict the branching decision
 * and represent an extreme point or extreme ray which is not feasible for the current pricing problem.
 *
 * @subsection BRANCHMASTERSOLVED
 *
 * The BRANCHMASTERSOLVED callback is called whenever the relaxation of a node is solved for the first time. It can be used
 * to collect statistical information about branching decision, e.g., update pseudocost values.
 *
 * @subsection BRANCHDATADELETE
 *
 * The BRANCHDATADELETE callback is called when a branch-and-bound node is freed and should free the branching data.
 *
 */

 /*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DEC_DECOMP How to store the structure information
 *
 * struct_decomp.h is responsible for storing structure information. The memory has to be allocated by caller and is freed
 * later. You have to use getter and setter functions ins pub_decomp.h to fill this structure.
 *
 * Very quick until more elaborate, these are the relevant fields of the structure and its basic usage:
 *  - <code>subscipcons   </code>
 *   - an array of array of constraints in each block
 *   - Usage: <code>subscipcons[b][c]</code> is constraint <em>c</em> in block <em>b</em>
 *  - <code>nsubscipconss </code>
 *   - an array of the number of constraints in each block
 *   - Usage: <code>nsubscipcons[b]</code> gives you the number of constraints of block <em>b</em>
 *  - <code>subscipvars   </code>
 *   - an array of arrays of variables in each block
 *   - Usage: <code>subscipvars[b][v]</code> is variable <em>v</em> in block <em>b</em>
 *  - <code>nsubscipvars  </code>
 *   - an array of the number of variables in each block
 *   - Usage: <code>nsubscipvars[b]</code> gives you the number of variables of block <em>b</em>
 *  - <code>nblocks       </code>
 *   - number of blocks/pricing problems in the matrix resp. reformulation
 *  - <code>type          </code>
 *   - Type of the decomposition (DEC_DECTYPE_ARROWHEAD is the most general)
 *   - Current supported types: DEC_DECTYPE_DIAGONAL, DEC_DECTYPE_BORDERED, DEC_DECTYPE_ARROWHEAD, DEC_DECTYPE_UNKNOWN
 *  - <code>constoblock   </code>
 *   - SCIP_HASHMAP linking constraints to blocks
 *   - Usage: <code>SCIPhashmapGetImage(constoblock, cons)</code> returns <em>b+1</em>, where <em>b</em> is the block of
 *     constraint cons. This map is somehow inverse of the subscipconss array.
 *  - <code>vartoblock    </code>
 *   - SCIP_HASHMAP linking variables to blocks
 *   - Usage: <code>SCIPhashmapGetImage(vartoblock, var)</code> returns <em>b+1</em>, where <em>b</em> is the block of
 *     variable var. This map is somehow inverse of the subscipcvars array.
 *  - <code>linkingvars   </code>
 *   - array of linking variables (to be in the master)
 *   - Usage: <code>linkingvars[v]</code> is linking variable number <em>v</em>
 *  - <code>nlinkingvars  </code>
 *   - number of linking variables
 *  - <code>linkingconss  </code>
 *   - array of linking constraints (to be in the master)
 *   - Usage: <code>linkingcons[c]</code> is linking constraint number <em>c</em>
 *  - <code>nlinkingconss </code>
 *   - number of linking constraints
 */
