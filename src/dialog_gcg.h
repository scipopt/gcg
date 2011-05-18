/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   dialog_gcg.h
 * @brief  gcg user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DIALOG_GCG_H__
#define __SCIP_DIALOG_GCG_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** dialog execution method for the display statistics command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetMaster);

/** standard menu dialog execution method, that displays it's help screen if the remaining command line is empty */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecMenu);

/** standard menu dialog execution method, that doesn't display it's help screen */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecMenuLazy);

/** dialog execution method for the checksol command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecChecksol);

/** dialog execution method for the conflictgraph command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecConflictgraph);

/** dialog execution method for the display branching command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayBranching);

/** dialog execution method for the display conflict command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayConflict);

/** dialog execution method for the display conshdlrs command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayConshdlrs);

/** dialog execution method for the display displaycols command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayDisplaycols);

/** dialog execution method for the display heuristics command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayHeuristics);

/** dialog execution method for the display memory command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayMemory);

/** dialog execution method for the display nodeselectors command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayNodeselectors);

/** dialog execution method for the display parameters command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayParameters);

/** dialog execution method for the display presolvers command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayPresolvers);

/** dialog execution method for the display problem command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayProblem);

/** dialog execution method for the display propagators command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayPropagators);

/** dialog execution method for the display readers command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayReaders);

/** dialog execution method for the display separators command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplaySeparators);

/** dialog execution method for the display solution command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplaySolution);

/** dialog execution method for the display statistics command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayStatistics);

/** dialog execution method for the display transproblem command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayTransproblem);

/** dialog execution method for the display value command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayValue);

/** dialog execution method for the display varbranchstatistics command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayVarbranchstatistics);

/** dialog execution method for the display transsolution command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayTranssolution);

/** dialog execution method for the help command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecHelp);

/** dialog execution method for the free command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecFree);

/** dialog execution method for the newstart command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecNewstart);

/** dialog execution method for the optimize command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecOptimize);

/** dialog execution method for the presolve command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecPresolve);

/** dialog execution method for the quit command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecQuit);

/** dialog execution method for the read command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecRead);

/** dialog execution method for the set default command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetDefault);

/** dialog execution method for the set load command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetLoad);

/** dialog execution method for the set save command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetSave);

/** dialog execution method for the set diffsave command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetDiffsave);

/** dialog execution method for the set parameter command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetParam);

/** dialog description method for the set parameter command */
extern
SCIP_DECL_DIALOGDESC(GCGdialogDescSetParam);

/** dialog execution method for the set branching direction command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetBranchingDirection);

/** dialog execution method for the set branching priority command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetBranchingPriority);

/** dialog execution method for the set limits objective command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetLimitsObjective);

/** creates a root dialog */
extern
SCIP_RETCODE GCGcreateRootDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         root                /**< pointer to store the root dialog */
   );

/** includes or updates the gcg dialog menus in SCIP */
extern
SCIP_RETCODE SCIPincludeDialogGcg(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes or updates the "set" menu for each available parameter setting */
extern
SCIP_RETCODE SCIPincludeDialogGcgSet(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
