/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nodesel_master.h
 * @brief  node selector for depth first search
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NODESEL_MASTER_H__
#define __SCIP_NODESEL_MASTER_H__


#include "scip/scip.h"


/** creates the node selector for depth first search and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeNodeselMaster(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
void GCGnodeselMasterSetOrigscip(
   SCIP*                 scip,
   SCIP*                 origscip
   );

#endif
