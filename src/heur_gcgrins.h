/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_gcgrins.h
 * @brief  GCG RINS primal heuristic
 * @author Timo Berthold
 * @author Christian Puchert
 * @ingroup PRIMALHEURISTICS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_GCGRINS_H__
#define __SCIP_HEUR_GCGRINS_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates RINS primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurGcgrins(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
