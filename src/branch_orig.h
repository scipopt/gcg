/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_orig.h
 * @brief  branching rule for original problem in gcg
 * @author Gerald Gamrath
 * @ingroup BRANCHINGRULES
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_ORIG_H__
#define __SCIP_BRANCH_ORIG_H__


#include "scip/scip.h"


/** creates the branching on original variable branching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchruleOrig(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
