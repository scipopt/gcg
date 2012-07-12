/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dialog_gcg.h
 * @brief  gcg user interface dialog
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @ingroup DIALOGS
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
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayStatistics);

/** dialog execution method for the display detectors command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayDetectors);

/** dialog execution method for the master command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetMaster);

/** dialog execution method for the detect command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecDetect);

/** dialog execution method for the detect command */
extern
SCIP_DECL_DIALOGEXEC(GCGdialogExecOptimize);

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

#ifdef __cplusplus
}
#endif

#endif
