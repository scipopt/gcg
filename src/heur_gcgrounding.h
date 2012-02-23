/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_gcgrounding.h
 * @brief  LP gcgrounding heuristic that tries to recover from intermediate infeasibilities
 * @author Tobias Achterberg
 * @author Christian Puchert
 * @ingroup PRIMALHEURISTICS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_GCGROUNDING_H__
#define __SCIP_HEUR_GCGROUNDING_H__


#include "scip/scip.h"


/** creates the GCG rounding heuristic with infeasibility recovering and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurGcgrounding(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
