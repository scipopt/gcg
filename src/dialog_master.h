/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dialog_master.h
 * @ingroup DIALOGS
 * @brief  user interface dialog for master problem
 * @author Tobias Achterberg
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DIALOG_MASTER_H__
#define __SCIP_DIALOG_MASTER_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** dialog execution method telling that a command is not available */
extern
SCIP_DECL_DIALOGEXEC(GCGmasterDialogExecNotAvailable);

/** creates a root dialog */
extern
SCIP_RETCODE GCGcreateRootMasterDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         root                /**< pointer to store the root dialog */
   );

/** includes or updates the master dialog menus in GCG */
extern
SCIP_RETCODE SCIPincludeDialogMaster(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
