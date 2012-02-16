/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_colgenfeaspump.h
 * @brief  column generation based feasibility pump primal heuristic
 * @author Christian Puchert
 * @ingroup PRIMALHEURISTICS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_COLGENFEASPUMP_H__
#define __SCIP_HEUR_COLGENFEASPUMP_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the colgenfeaspump primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurColgenfeaspump(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
