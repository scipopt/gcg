/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_greedycolsel.h
 * @brief  greedy column selection primal heuristic
 * @author Christian Puchert
 * @ingroup PRIMALHEURISTICS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_GREEDYCOLSEL_H__
#define __SCIP_HEUR_GREEDYCOLSEL_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the greedy column selection primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurGreedycolsel(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
