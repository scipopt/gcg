/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal.c
 * @brief  main document page
 * @author   Gerald Gamrath
 * @author   Martin Bergner
 * @author   Christian Puchert
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Generic Column Genration
 * @version  0.9.0.1
 * @author   Gerald Gamrath
 * @author   Martin Bergner
 * @author   Christian Puchert
 *
 * <table cellpadding="0px" border="0" width="100%">
 *   <tr>
 *     <td nowrap >
 * <b>GCG Information</b>
 *
 * - \ref CHANGELOG    "Change log"
 * - \ref RELEASENOTES "Release notes"
 * - \ref INSTALL      "Installation information"
 *
 * <b>SCIP Authors</b>
 * - <a class="el" href="AUTHORS.html#main">Current main developers</a>
 * - <a class="el" href="AUTHORS.html#further">Further developers</a>
 * 
 * <b>Further Documentation</b>
 * - \ref DETECT "How to write a custom structure detection"
 * - \ref HEUR "How to write custom heuristics"
 * 
 *     </td>
 *     <td valign="bottom" width="200">
 *       \image html scippy.png
 *     </td>
 *   </tr>
 * </table>
 */

/**@page CHANGELOG CHANGELOG
 *
 * \verbinclude CHANGELOG
 */

/**@page INSTALL INSTALL
 *
 * \verbinclude INSTALL
 */

/**@page RELEASENOTES Release notes
 *
 * \verbinclude GCG-release-notes-1.0
 *
 */

/**@page AUTHORS SCIP Authors
 * \htmlinclude authors.inc  
 */

/**@page INSTALL Installation information
 * \verbinclude INSTALL  
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CODE Coding style guidelines
 *
 * We follow the following coding style guidlines and we recommend to use it in your code.
 *
 * - Indentation is 3 spaces. No tabs anywhere in the code.
 * - Always only one declaration in a line.
 * - Braces are on a new line and not indented.
 * - Spaces around all operators.
 * - No spaces between control structure keywords like "if", "for", "while", "switch" and the corresponding brackets.
 * - No spaces between a function name and the parenthesis in both the definition and function calls.
 * - Use assert() to show preconditions for the parameters, invariants and postconditions.
 * - All global functions start with "SCIP". In the usual naming scheme this is followed by the object and a method name
 *   like in SCIPlpAddRow(). Functions return TRUE or FALSE should be named like SCIPisFeasEQ().
 * - Make all functions that are not used outside the module 'static'. Naming should start with a lower case letter.
 * - Variable names should be all lower case.
 * - For each structure there is a typedef with the name in all upper case.
 * - Defines should be named all upper case.
 * - Document functions, parameters, and variables in a doxygen conform way.
 *
 * As an example have a look at tree.c .
 */

/**@defgroup TYPEDEFINITIONS Type Definitions
 * This page lists headers which contain type definitions of callback methods.
 *
 * All headers below include the descriptions of callback methods of
 * certain plugins. For more detail see the corresponding header.
 */

/**@defgroup BRANCHINGRULES Branching Rules 
 * @brief This page contains a list of all branching rule which are currently available. 
 *
 * A detailed description what a branching rule does and how to add a branching rule to \SCIP can be found 
 * \ref BRANCH "here".
 */

/**@defgroup CONSHDLRS  Constraint Handler 
 * @brief This page contains a list of all constraint handlers which are currently available.
 *
 * A detailed description what a constraint handler does and how to add a constraint handler to \SCIP can be found 
 * \ref CONS "here".
 */

/**@defgroup DIALOGS Dialogs 
 * @brief This page contains a list of all dialogs which are currently available.
 *
 * A detailed description what a dialog does and how to add a dialog to \SCIP can be found 
 * \ref DIALOG "here".
 */

/**@defgroup DISPLAYS Displays 
 * @brief This page contains a list of all displays (output columns)  which are currently available.
 *
 * A detailed description what a display does and how to add a display to \SCIP can be found 
 * \ref DISP "here".
 *
 */

/**@defgroup FILEREADERS File Readers 
 * @brief This page contains a list of all file readers which are currently available.
 *
 * A detailed description what a file reader does and how to add a file reader to \SCIP can be found 
 * \ref READER "here".
 */
 
/**@defgroup NODESELECTORS Node Selectors
 * @brief This page contains a list of all node selectors which are currently available.
 *
 * A detailed description what a node selector does and how to add a node selector to \SCIP can be found 
 * \ref NODESEL "here".
 */

/**@defgroup PRESOLVERS Presolvers
 * @brief This page contains a list of all presolvers which are currently available.
 *
 * A detailed description what a presolver does and how to add a presolver to \SCIP can be found 
 * \ref PRESOL "here".
 */

/**@defgroup PRICERS Pricers
 * @brief This page contains a list of all pricers which are currently available.
 *
 * Per default there exist no variable pricer. A detailed description what a variable pricer does and how to add a
 * variable pricer to \SCIP can be found \ref PRICER "here".
 */

/**@defgroup PRIMALHEURISTICS Primal Heuristics
 * @brief This page contains a list of all primal heuristics which are currently available.
 *
 * A detailed description what a primal heuristic does and how to add a primal heuristic to \SCIP can be found 
 * \ref HEUR "here".
 */


/**@defgroup RELAXATORS Relaxators
 * @brief This page contains a list of all relaxators which are currently available.
 */

/**@defgroup SEPARATORS Separators
 * @brief This page contains a list of all separators  which are currently available.
 *
 * A detailed description what a separator does and how to add a separator to \SCIP can be found 
 * \ref SEPA "here".
 */

/**@defgroup DETECTORS  Structure Enforcers 
 * @brief This page contains a list of all constraint handlers which are currently available.
 *
 * A detailed description what a constraint handler does and how to add a constraint handler to \SCIP can be found 
 * \ref DETECT "here".
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page HEUR How to add primal heuristics
 *
 * Feasible solutions can be found in two different ways during the traversal of the branch-and-bound tree. On the one 
 * hand, the solution of a node's relaxation may be feasible with respect to the constraints (including the integrality). 
 * On the other hand, feasible solutions can be discovered by primal heuristics.  
 * \n
 * A complete list of all primal heuristics contained in this release can be found \ref PRIMALHEURISTICS "here".
 *
 * In the following, we explain how the user can add an own primal heuristic.
 * Take the simple and fast LP rounding heuristic (src/scip/heur_simplerounding.c) as an example.
 * The idea of simple rounding is to iterate over all fractional variables of an LP solution and round them down,
 * if the variables appears only with nonnegative coefficients in the system Ax <= b and round them up if
 * the variables appears only with nonpositive coefficients.
 * If one of both conditions applies for each of the fractional variables, this will give a feasible solution.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the ObjHeur wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_HEUR... callback methods.
 *
 * Additional documentation for the callback methods of a primal heuristic can be found in the file type_heur.h.
 *
 * Here is what you have to do to implement a primal heuristic:
 * -# Copy the template files src/scip/heur_xyz.c and src/scip/heur_xyz.h into files named "heur_myheuristic.c"
 *    and "heur_myheuristic.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "myheuristic".
 * -# Adjust the properties of the primal heuristic (see \ref HEUR_PROPERTIES).
 * -# Define the primal heuristic data (see \ref HEUR_DATA). This is optional.
 * -# Implement the interface methods (see \ref HEUR_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref HEUR_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref HEUR_ADDITIONALCALLBACKS). This is optional.
 *
 *
 * @section HEUR_PROPERTIES Properties of a Primal Heuristic
 *
 * At the top of the new file "heur_myheuristic.c" you can find the primal heuristic properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the primal heuristic properties by calling the constructor
 * of the abstract base class ObjHeur from within your constructor.
 * Of course, all of them are of relevance, but the most important ones for controlling the performance
 * usually are HEUR_FREQ and HEUR_TIMING.
 * The properties you have to set have the following meaning:
 *
 * \par HEUR_NAME: the name of the primal heuristic.
 * This name is used in the interactive shell to address the primal heuristic.
 * Additionally, if you are searching for a primal heuristic with SCIPfindHeur(), this name is looked up.
 * Names have to be unique: no two primal heuristics may have the same name.
 *
 * \par HEUR_DESC: the description of the primal heuristic.
 * This string is printed as description of the primal heuristic in the interactive shell when you call "display heuristics".
 *
 * \par HEUR_DISPCHAR: the display character of the primal heuristic.
 * In the interactive shell, this character is printed in the first column of a status information row, if the primal 
 * heuristic found the feasible solution belonging to the primal bound. Note that a star '*' stands for an integral 
 * LP-relaxation.
 * In order to avoid confusion, display characters should be unique: no two primal heuristics should have the same display character.
 * You can get a list of all primal heuristics along with their display characters by entering "display heuristics" in the 
 * SCIP interactive shell.
 *
 * \par HEUR_PRIORITY: the priority of the primal heuristic.
 * At each of the different entry points of the primal heuristics during the solving process (see HEUR_TIMING), they are 
 * called in decreasing order of their priority. 
 * \n
 * The priority of a primal heuristic should be set according to the complexity of the heuristic and the likelihood to find
 * feasible solutions: primal heuristics that provide fast algorithms that often succeed in finding a feasible solution should have
 * a high priority (like simple rounding). In addition, the interaction between different types of primal heuristics should be taken into account.
 * For example, improvement heuristics, which try to generate improved solutions by inspecting one or more of the feasible
 * solutions that have already been found, should have a small priority (like Crossover which by default needs at least 3 feasible solutions).
 *
 * \par HEUR_FREQ: the default frequency for executing the primal heuristic.
 * The frequency together with the frequency offset (see HEUR_FREQOFS) defines the depth levels at which the execution
 * method of the primal heuristic \ref HEUREXEC is called. For example, a frequency of 7 together with a frequence offset 
 * of 5 means, that the \ref HEUREXEC callback is executed for subproblems that are in depth 5, 12, 19, ... of the branching tree. A 
 * frequency of 0 together with a frequence offset of 3 means, that the execution method is only called at those nodes that are in
 * depth level 3 (i.e., at most for \f$2^3 = 8\f$ nodes if binary branching is applied).
 * Typical cases are: A frequency of 0 and an offset of 0 which means that
 * the heuristic is only called at the root node and a frequency of -1 which disables the heuristic.
 * \n
 * The frequency can be adjusted by the user. The property of the primal heuristic only defines the default value of the 
 * frequency. If you want to have a more flexible control of when to execute the primal heuristic, you have to assign
 * a frequency of 1 and implement a check at the beginning of your execution method whether you really want to search for feasible
 * solutions or not. If you do not want to execute the method, set the result code to SCIP_DIDNOTRUN.
 *
 * \par HEUR_FREQOFS: the frequency offset for executing the primal heuristic.
 * The frequency offset defines the depth of the branching tree at which the primal heuristic is executed for the first 
 * time. For example, a frequency of 7 (see HEUR_FREQ) together with a frequency offset of 10 means, that the 
 * callback is executed for subproblems that are in depth 10, 17, 24, ... of the branching tree. In particular, assigning 
 * different offset values to heuristics of the same type, like diving heuristics, can be useful for evenly spreading the 
 * application of these heuristics across the branch-and-bound tree.
 * Note that if the frequency is equal to 1, the heuristic is applied for all nodes with depth level larger or equal to
 * the frequency offset.
 *
 * \par HEUR_MAXDEPTH: the maximal depth level for executing the primal heuristic.
 * This parameter denotes the maximal depth level in the branching tree up to which the execution method of the primal 
 * heuristic is called. Use -1 for no limit (a usual case). 
 *
 * \par HEUR_TIMING: the execution timing of the primal heuristic.
 * Primal heuristics have different entry points during the solving process and the execution timing parameter defines the
 * entry point at which the primal heuristic is executed first. 
 * \n
 * The primal heuristic can be called first:
 * - before the processing of the node starts (SCIP_HEURTIMING_BEFORENODE)
 * - after each LP solving during the cut-and-price loop (SCIP_HEURTIMING_DURINGLPLOOP) 
 * - after the cut-and-price loop was finished (SCIP_HEURTIMING_AFTERLPLOOP) 
 * - after the processing of a node <em>with solved LP</em>  was finished (SCIP_HEURTIMING_AFTERLPNODE)
 * - after the processing of a node <em>without solved LP</em> was finished (SCIP_HEURTIMING_AFTERPSEUDONODE)
 * - after the processing of the last node in the current plunge was finished, <em>and only if the LP was solved for 
 *   this node</em> (SCIP_HEURTIMING_AFTERLPPLUNGE) 
 * - after the processing of the last node in the current plunge was finished, <em>and only if the LP was not solved 
 *   for this node</em> (SCIP_HEURTIMING_AFTERPSEUDOPLUNGE).
 * \par
 * A plunge is the successive solving of child and sibling nodes in the search tree.
 * The flags listed above can be combined to call the heuristic at multiple times by concatenating them with a bitwise OR. 
 * Two useful combinations are already predefined:
 * - after the processing of a node was finished (SCIP_HEURTIMING_AFTERNODE; combines SCIP_HEURTIMING_AFTERLPNODE and
 *   SCIP_HEURTIMING_AFTERPSEUDONODE)
 * - after the processing of the last node in the current plunge was finished (SCIP_HEURTIMING_AFTERPLUNGE; combines
 *   SCIP_HEURTIMING_AFTERLPPLUNGE and SCIP_HEURTIMING_AFTERPSEUDOPLUNGE)
 * \par
 * Calling a primal heuristic "before the processing of the node starts" is particularly useful for heuristics 
 * that do not need to access the LP solution of the current node. If such a heuristic finds a feasible solution, the 
 * leaves of the branching tree exceeding the new primal bound are pruned. It may happen that even the current node can 
 * be cut off without solving the LP relaxation. Combinatorial heuristics, like the farthest insert heuristic for the TSP 
 * (see examples/TSP/src/HeurFarthestInsert.cpp), are often applicable at this point.
 * \n
 * Very fast primal heuristics that require an LP solution can also be called "after each LP solving during the 
 * cut-and-price loop". Rounding heuristics, like the simple and fast LP rounding heuristic 
 * (src/scip/heur_simplerounding.c), belong to this group of primal heuristics. 
 * \n
 * Most heuristics, however, are called either after a node was completely processed 
 * (e.g. expensive rounding heuristics like RENS), or even only after a full plunge was finished (e.g., diving heuristics). 
 *
 * \par HEUR_USESSUBSCIP: Does the heuristic use a secondary SCIP instance? 
 * Some heuristics and separators solve MIPs or SAT problems and use a secondary SCIP instance therefor. Examples are
 * Large Neighborhood Search heuristics such as RINS and Local Branching or the CGMIP separator. To avoid recursion,
 * these plugins usually deactivate all other plugins that solve MIPs. If a heuristic uses a secondary SCIP instance,
 * this parameter has to be TRUE and it is recommended to call SCIPsetSubscipsOff() for the secondary SCIP instance.
 *
 * Computational experiments indicate that for the overall performance of a MIP solver, it is important to evenly 
 * spread the application of the heuristics across the branch-and-bound tree. Thus, the assignment of the parameters 
 * HEUR_FREQ, HEUR_FREQOFS, and HEUR_TIMING should contribute to this aim.
 *
 * Note that all diving heuristics in the SCIP distribution (see, e.g., src/scip/heur_guideddiving.c) check whether other diving
 * heuristics have already been called at the current node. This can be done by comparing SCIPgetLastDivenode(scip) and
 * SCIPgetNNodes(scip). If the two are equal, and if the current node is not the root node (SCIPgetDepth(scip) > 0), diving
 * heuristics should be delayed by returning the result code 'SCIP_DELAYED'. This is an additional contribution to the goal of
 * not calling multiple similar heuristics at the same node.
 *
 *
 * @section HEUR_DATA Primal Heuristic Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_HeurData".
 * In this data structure, you can store the data of your primal heuristic. For example, you should store the adjustable 
 * parameters of the primal heuristic or a working solution in this data structure. 
 * If you are using C++, you can add primal heuristic data as usual as object variables to your class.
 * \n
 * Defining primal heuristic data is optional. You can leave the struct empty.
 *
 *
 * @section HEUR_INTERFACE Interface Methods
 *
 * At the bottom of "heur_myheuristic.c" you can find the interface method SCIPincludeHeurMyheuristic(), which also 
 * appears in "heur_myheuristic.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the primal heuristic by calling the method SCIPincludeHeur().
 * It is called by the user, if he wants to include the primal heuristic, i.e., if he wants to use the primal heuristic
 * in his application.
 *
 * If you are using primal heuristic data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_HeurData afterwards.
 *
 * You may also add user parameters for your primal heuristic, see the method SCIPincludeHeurFeaspump() in 
 * src/scip/heur_oneopt.c for an example where a single Boolean parameter is added.
 *
 * 
 * @section HEUR_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Primal Heuristic
 *
 * Primal heuristic plugins have only one fundamental callback method, namely the HEUREXEC method.
 * This method has to be implemented for every primal heuristic; the other callback methods are optional.
 * In the C++ wrapper class ObjHeur, the scip_exec() method (which corresponds to the HEUREXEC callback) is a virtual
 * abstract member function. You have to implement it in order to be able to construct an object of your primal heuristic 
 * class.
 *
 * Additional documentation to the callback methods can be found in type_heur.h.
 *
 * @subsection HEUREXEC
 *
 * The HEUREXEC callback is called at different positions during the node processing loop, see HEUR_TIMING. It should
 * search for feasible solutions and add them to the solution pool. For creating a new feasible solution, the 
 * methods SCIPcreateSol() and SCIPsetSolVal() can be used. Afterwards, the solution can be added to the storage by 
 * calling the method SCIPtrySolFree() (or SCIPtrySol() and SCIPfreeSol()).
 *
 * The HEUREXEC callback gets a SCIP pointer, a pointer to the heuristic itself, the current point in the 
 * solve loop and a result pointer as input (see type_heur.h).
 *
 * The heuristic has to set the result pointer appropriately!
 * Therefore it has the following options:
 *  - finding at least one feasible solution (result SCIP_FOUNDSOL)
 *  - stating that the primal heuristic searched, but did not find a feasible solution (result SCIP_DIDNOTFIND)
 *  - stating that the primal heuristic was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the primal heuristic was skipped, but should be called again (result SCIP_DELAYED).
 *
 *
 * @section HEUR_ADDITIONALCALLBACKS Additional Callback Methods of a Primal Heuristic
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 *
 * @subsection HEURFREE
 *
 * If you are using primal heuristic data, you have to implement this method in order to free the primal heuristic data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_HEURFREE(heurFreeMyheuristic)
 * {
 *    SCIP_HEURDATA* heurdata;
 *  
 *    heurdata = SCIPheurGetData(heur);
 *    assert(heurdata != NULL);
 *
 *    SCIPfreeMemory(scip, &heurdata);
 *
 *    SCIPheurSetData(heur, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you have allocated memory for fields in your primal heuristic data, remember to free this memory 
 * before freeing the primal heuristic data itself. 
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection HEURINIT
 *
 * The HEURINIT callback is executed after the problem was transformed.
 * The primal heuristic may, e.g., use this call to initialize his primal heuristic data.
 *
 * @subsection HEURCOPY
 *
 * The HEURCOPY callback is executed when a SCIP instance is copied, e.g. to 
 * solve a sub-SCIP. By
 * defining this callback as
 * <code>NULL</code> the user disables the execution of the specified 
 * heuristic for all copied SCIP instances. This may deteriorate the performance 
 * of primal heuristics using sub-SCIPs.
 *
 * @subsection HEUREXIT
 *
 * The HEUREXIT callback is executed before the transformed problem is freed.
 * In this method, the primal heuristic should free all resources that have been allocated for the solving process in 
 * HEURINIT.
 *
 * @subsection HEURINITSOL
 *
 * The HEURINITSOL callback is executed when the presolving was finished and the branch and bound process is about to 
 * begin. The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 * @subsection HEUREXITSOL
 *
 * The HEUREXITSOL callback is executed before the branch and bound process is freed. The primal heuristic should use this
 * call to clean up its branch and bound data, which was allocated in HEURINITSOL.
 */

/**@page FAQ Frequently Asked Questions (FAQ)
 * \htmlinclude faq.inc  
 */

 /*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DETECT How to add structure detectors
 *
 * Structure detectors are used to detect or enforce a structure suitable for Dantzig-Wolfe Reformulation (DWR). 
 * \n 
 * A complete list of all detectors/enforcers contained in this release can be found \ref DETECTORS "here".
 *
 * In the following, we explain how the user can add its own structure enforcement plugin.
 * Take the bordered block diagonal detector (src/dec_borderheur.c) as an example.
 * As all other default plugins, it is written in C. There is currently no C++ wrapper available.
 *
 * Additional documentation for the callback methods of a structure detection, in particular for their input parameters, 
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
 * \par DEC_NAME: the name of the detector.
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
 * priority. An easy way to list the priorities of all detectors "display detectors" (TODO) in the interactive shell of GCG.
 * 
 *
 * @section DEC_DATA Detector Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct DEC_DetectorData".
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
 * This method has only to be adjusted slightly.
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
 * detector data, see \ref DETECTORFREE.
 *
 * You may also add user parameters for your detector, see \ref PARAM for how to add user parameters and 
 * the method SCIPincludeDetectionBorderheur() in src/dec_borderheur.c for an example.
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
 * The DETECTSTRUCTURE callback is called detection loop and should perform the actual detection.
 * It should inspect the problem instance at hand and deduct some structure from the constraint matrix.
 *
 * Typical methods called by a detector are, for example, … .
 *
 * @subsection GETPRIORITY
 * 
 * The GETPRIORITY callback is called by cons_decomp in order to get the priority of the detector.
 *
 * @subsection SETSTRUCTDECOMP
 * 
 * The SETSTRUCTDECOMP callback is called in order to notify the detector where to store the structure 
 * information for GCG. See \ref DEC_DECOMP "\"Storing structure information\"" for an explanation of what data needs to be stored there.
 *
 * @section DEC_ADDITIONALCALLBACKS Additional Callback Methods of a Detector
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 *
 * @subsection DECINIT
 *
 * The INITDETECTOR callback is executed after the problem was transformed.
 * The detector may, e.g., use this call to initialize his detector data.
 * The difference between the original and the transformed problem is explained in 
 * "What is this thing with the original and the transformed problem about?" on \ref FAQ. 
 *
 * @subsection DECEXIT
 *
 * If you are using detection data (see \ref DEC_DATA and \ref DEC_INTERFACE), you have to implement this method in order to free the detection data.
 * This can be done by the following procedure: (TODO)
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
 * In this method, the detector should free all resources that have been allocated for the solving process in PRESOLINIT.
 *
 */

 /*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DEC_DECOMP Storing structure information
 *
 * struct_decomp.h is responsible for storing structure information. The memory has to be allocated by caller and is freed
 * later
 * 
 * Very quick until more elaborate, these are the relavant fields of the structure:
 *  - subscipconss - an array of array of constraints in each block
 *  - nsubscipconss - an array of the number of constraints in each block
 *  - subscipvars - an array of arrays of variables in each block 
 *  - nsubscipvars - an array of the number of constraints in each blocks
 *  - nblocks - number of blocks
 *  - type - Type of the decomposition (DEC_STAIRCASE is the most general)
 *  - constoblock - SCIP_HASHMAP linking constraints to blocks
 *  - vartoblock - SCIP_HASHMAP linking variables to blocks
 *  - linkingvars - array of linking variables (to be in the master)
 *  - nlinkingvars - number of linking variables
 *  - linkingconss - array of linking constraints (to be in the master)
 *  - nlinkingconss - number of linking constraints
 */

