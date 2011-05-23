/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident ""

/**@file   heur_gcgshifting.h
 * @brief  LP gcgrounding heuristic that tries to recover from intermediate infeasibilities and shifts continuous variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_GCGSHIFTING_H__
#define __SCIP_HEUR_GCGSHIFTING_H__


#include "scip/scip.h"


/** creates the GCG shifting heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurGcgshifting(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
