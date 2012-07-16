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

/**@mainpage Generic Column Genration
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
 * and has more than 40,000 lines of C code.
 *
 *
 * <table cellpadding="0px" border="0" width="100%">
 *   <tr>
 *     <td nowrap >
 * <b>How to get started</b>
 *
 * - \ref INSTALL      "Installation information"
 * - \ref FILEFORMATS  "Input file formats"
 *
 * <table cellpadding="0px" border="0" width="100%">
 *   <tr>
 *     <td nowrap >
 * <b>Further Information</b>
 * - \ref AUTHORS      "Current GCG developers"
 * - \ref CHANGELOG    "Change log"
 * - \ref RELEASENOTES "Release notes"
 *
 * @version  1.0
 *
 * <b>Further Documentation</b>
 * - \ref DETECT "How to write a custom structure detection"
 * - \ref HEUR "How to write custom heuristics"
 * - \ref PUBLICMETHODS "List of callable functions"
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

/**@page AUTHORS GCG Authors
 * \htmlinclude authors.inc
 */

/**@page INSTALL Installation information
 * \verbinclude INSTALL
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page FILEFORMATS Input file formats supported by GCG.
 *
 * GCG supports all file formats supported by SCIP to read in problems, solution, etc.
 * E.g., the original problem can be read in as an .lp, .mps, or .cip file.
 *
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

/**@defgroup DETECTORS Detectors
 * @brief This page contains a list of all detectors which are currently available.
 *
 * A detailed description what a detector does and how to add a detector to GCG can be found
 * \ref DETECT "here".
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
 * /* do not execute the heuristic on invalid relaxation solutions
 *  * (which is the case if the node has been cut off)
 *  */
 *  if( !SCIPisRelaxSolValid(scip) )
 *     return SCIP_OKAY;
 *
 *  /* only call heuristic, if an optimal LP solution is at hand */
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
 * Take the basic detector (src/dec_connected.c) as an example.
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
 * priority. An easy way to list the priorities of all detectors "display detectors" in the interactive shell of GCG.
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
 * It needs to store the structure information in DECDECOMP.
 *
 * Typical methods called by a detector are, for example, SCIPgetVars(), SCIPGetConss, DECdecompSetNBlocks() .
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
 * In this method, the detector should free all resources that have been allocated for the solving process in \ref DETECTORINIT.
 *
 */

 /*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DEC_DECOMP Storing structure information
 *
 * struct_decomp.h is responsible for storing structure information. The memory has to be allocated by caller and is freed
 * later
 *
 * Very quick until more elaborate, these are the relavant fields of the structure:
 *  - subscipconss - an array of array of constraints in each block - array[blocknr][constraintid]
 *  - nsubscipconss - an array of the number of constraints in each block
 *  - subscipvars - an array of arrays of variables in each block - array[blocknr][varid]
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

