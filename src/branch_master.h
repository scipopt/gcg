/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   branch_master.h
 * @brief  branching rule for master problem
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_MASTER_H__
#define __SCIP_BRANCH_MASTER_H__


#include "scip/scip.h"


/** creates the most infeasible LP braching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchruleMaster(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
