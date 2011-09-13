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
 * @version  0.9
 * @author   Gerald Gamrath
 * @author   Martin Bergner
 * @author   Christian Puchert
 *
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
 * \ref CONS "here".
 */

 /*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page PRESOL How to add presolvers
 *
 * Presolvers are used to reduce the size of the model by removing irrelevant information like redundant constraints,
 * to strengthen the LP relaxation by exploiting integrality information, and to extract useful information in the 
 * presolving step.
 * Constraint based presolving is done in the CONSPRESOL callback methods of the constraint handlers, see \ref CONSPRESOL.
 * The presolver plugins complement the constraint based presolving by additional, usually optimality based, presolving
 * reductions. 
 * \n 
 * A complete list of all presolvers contained in this release can be found \ref PRESOLVERS "here".
 *
 * In the following, we explain how the user can add its own presolver.
 * Take the dual fixing presolver (src/scip/presol_dualfix.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the ObjPresol wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_PRESOL... callback methods.
 *
 * Additional documentation for the callback methods of a presolver, in particular for their input parameters, 
 * can be found in the file type_presol.h.
 *
 * Here is what you have to do to implement a presolver:
 * -# Copy the template files src/scip/presol_xyz.c and src/scip/presol_xyz.h into files named "presol_mypresolver.c"
 *    and "presol_mypresolver.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "mypresolver".
 * -# Adjust the properties of the presolver (see \ref PRESOL_PROPERTIES).
 * -# Define the presolver data (see \ref PRESOL_DATA). This is optional.
 * -# Implement the interface methods (see \ref PRESOL_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref PRESOL_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref PRESOL_ADDITIONALCALLBACKS). This is optional.
 *
 * 
 * @section PRESOL_PROPERTIES Properties of a Presolver
 *
 * At the top of the new file "presol_mypresolver.c", you can find the presolver properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the presolver properties by calling the constructor
 * of the abstract base class ObjPresol from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par PRESOL_NAME: the name of the presolver.
 * This name is used in the interactive shell to address the presolver.
 * Additionally, if you are searching for a presolver with SCIPfindPresol(), this name is looked up.
 * Names have to be unique: no two presolvers may have the same name.
 *
 * \par PRESOL_DESC: the description of the presolver.
 * This string is printed as description of the presolver in the interactive shell.
 *
 * \par PRESOL_PRIORITY: the priority of the presolver.
 * In each presolving round, the presolvers and presolving methods of the constraint handlers are called in
 * a predefined order, which is given by the priorities of the presolvers and the check priorities of the
 * constraint handlers, see \ref CONS_PROPERTIES.
 * First, the presolvers with non-negative priority are called in the order of decreasing priority.
 * Next, the presolving methods of the different constraint handlers are called in the order of decreasing check
 * priority.
 * Finally, the presolvers with negative priority are called in the order of decreasing priority.
 * \n
 * The priority of the presolver should be set according to the complexity of the presolving algorithm and the impact of the reduction:
 * presolvers that provide fast algorithms that usually have a high impact (i.e., remove lots of variables or tighten 
 * bounds of many variables) should have a high priority. An easy way to list the 
 * priorities of all presolvers and constraint handlers is to type "display presolvers" and "display conshdlrs" in 
 * the interactive shell of SCIP.
 *
 * \par PRESOL_MAXROUNDS: the default maximal number of rounds the presolver participates in.
 * The presolving is conducted in rounds: the presolvers and presolving methods of the constraint handlers
 * are called iteratively until no more reductions have been found or some other abort criterion applies.
 * The "maxrounds" parameter of a presolver imposes a limit on the number of presolving rounds in which the
 * presolver is called. The PRESOL_MAXROUNDS property specifies the default value for this parameter.
 * A value of -1 represents an unlimited number of rounds.
 *
 * \par PRESOL_DELAY: the default for whether the presolver should be delayed, if other presolvers found reductions.
 * If the presolver is marked to be delayed, it is only executed if no other presolvers found a reduction during the current
 * presolving round.
 * If the presolver is very expensive, you may want to mark it to be delayed until all cheap presolving methods have been executed.
 * 
 *
 * @section PRESOL_DATA Presolver Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_PresolData".
 * In this data structure, you can store the data of your presolver. For example, you should store the adjustable parameters
 * of the presolver in this data structure.  
 * If you are using C++, you can add presolver data as usual as object variables to your class.
 * \n
 * Defining presolver data is optional. You can leave this struct empty.
 *
 *
 * @section PRESOL_INTERFACE Interface Methods
 *
 * At the bottom of "presol_mypresolver.c", you can find the interface method SCIPincludePresolMypresolver(), 
 * which also appears in "presol_mypresolver.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the presolver by calling the method
 * SCIPincludePresol(). 
 * SCIPincludePresolMypresolver() is called by the user, if he wants to include the presolver, 
 * i.e., if he wants to use the presolver in his application.
 *
 * If you are using presolver data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &presoldata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_PresolData afterwards. For freeing the 
 * presolver data, see \ref PRESOLFREE.
 *
 * You may also add user parameters for your presolver, see \ref PARAM for how to add user parameters and 
 * the method SCIPincludePresolProbing() in src/scip/presol_probing.c for an example.
 *
 * 
 * @section PRESOL_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Presolver
 *
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain 
 * an operational algorithm. Presolver plugins have only one fundamental callback method, namely the PRESOLEXEC method.
 * This method has to be implemented for every presolver; the other callback methods are optional.
 * In the C++ wrapper class ObjPresol, the scip_exec() method (which corresponds to the PRESOLEXEC callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your presolver class.
 *
 * Additional documentation to the callback methods, in particular to their input parameters, 
 * can be found in type_presol.h.
 *
 * @subsection PRESOLEXEC
 *
 * The PRESOLEXEC callback is called inside the presolving loop and should perform the actual presolving reductions.
 * It should inspect the problem instance at hand and simplify it by tightening bounds of variables, aggregating or fixing
 * variables, changing the type of variables, modifying the graph that represents the instance of your application, and
 * the like.
 *
 * Typical methods called by a presolver are, for example, SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), SCIPtightenVarLb(),
 * and SCIPtightenVarUb().
 *
 *
 * @section PRESOL_ADDITIONALCALLBACKS Additional Callback Methods of a Presolver
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 *
 * @subsection PRESOLFREE
 *
 * If you are using presolver data (see \ref PRESOL_DATA and \ref PRESOL_INTERFACE), you have to implement this method in order to free the presolver data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_PRESOLFREE(presolFreeMypresolver)
 * {
 *    SCIP_PRESOLDATA* presoldata;
 *  
 *    presoldata = SCIPpresolGetData(presol);
 *    assert(presoldata != NULL);
 *
 *    SCIPfreeMemory(scip, &presoldata);
 *
 *    SCIPpresolSetData(presol, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you have allocated memory for fields in your presolver data, remember to free this memory 
 * before freeing the presolver data itself. 
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection PRESOLINIT
 *
 * The PRESOLINIT callback is executed after the problem was transformed.
 * The presolver may, e.g., use this call to initialize his presolver data.
 * The difference between the original and the transformed problem is explained in 
 * "What is this thing with the original and the transformed problem about?" on \ref FAQ. 
 *
 * @subsection PRESOLCOPY
 *
 * The PRESOLCOPY callback is executed when a SCIP instance is copied, e.g. to 
 * solve a sub-SCIP. By
 * defining this callback as
 * <code>NULL</code> the user disables the execution of the specified 
 * presolver for all copied SCIP instances. This may deteriorate the performance 
 * of primal heuristics using sub-SCIPs.
 *
 * @subsection PRESOLEXIT
 *
 * The PRESOLEXIT callback is executed before the transformed problem is freed.
 * In this method, the presolver should free all resources that have been allocated for the solving process in PRESOLINIT.
 *
 * @subsection PRESOLINITPRE
 *
 * The PRESOLINITPRE callback is executed when the presolving is about to begin.
 * The presolver may use this call to initialize its presolving data which only need to exist during the presolving stage.
 *
 * @subsection PRESOLEXITPRE
 *
 * The PRESOLEXITPRE callback is executed after presolving has been finished and before the branch and bound process begins.
 * The presolver should use this call to clean up its presolving data, which was allocated in PRESOLINITPRE.
 */

