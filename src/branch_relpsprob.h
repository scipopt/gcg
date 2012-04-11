/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_relpsprob.h
 * @brief  reliable pseudo costs branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_RELPSPROB_H__
#define __SCIP_BRANCH_RELPSPROB_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the reliable pseudo cost braching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchruleRelpsprob(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** execution reliability pseudo cost probing branching with the given branching candidates */
extern
SCIP_RETCODE SCIPgetRelpsprobBranchVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            branchcands,        /**< brancing candidates */
   SCIP_Real*            branchcandssol,     /**< solution value for the branching candidates */
   SCIP_Real*            branchcandsfrac,    /**< fractional part of the branching candidates, zero possible */
   int                   nbranchcands,       /**< number of branching candidates */
   int                   nvars,              /**< number of variables to be watched by bdchgdata */
   SCIP_RESULT*          result,             /**< pointer to the result of the execution */
   SCIP_VAR**            branchvar           /**< pointer to the variable to branch on */
   );

#ifdef __cplusplus
}
#endif

#endif
