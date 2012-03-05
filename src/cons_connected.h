/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_connected.h
 * @brief  constraint handler for connected constraints
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_CONNECTED_H__
#define __SCIP_CONS_CONNECTED_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for connected constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrConnected(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a connected constraint */
extern
SCIP_RETCODE SCIPcreateConsConnected(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name                /**< name of constraint */
   );

/** returns whether a block diagonal structure was found */
extern
SCIP_Bool SCIPisMatrixBlockDiagonal(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
