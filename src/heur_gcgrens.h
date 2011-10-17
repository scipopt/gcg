/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_gcgrens.h
 * @brief  GCG RENS primal heuristic
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_GCGRENS_H__
#define __SCIP_HEUR_GCGRENS_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates RENS primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurGcgrens(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** main procedure of the RENS heuristic, creates and solves a subMIP */
SCIP_RETCODE SCIPapplyGcgrens(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< result data structure                                          */
   SCIP_Real             minfixingrate,      /**< minimum percentage of integer variables that have to be fixed  */
   SCIP_Real             minimprove,         /**< factor by which RENS should at least improve the incumbent     */
   SCIP_Longint          maxnodes,           /**< maximum number of  nodes for the subproblem                    */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem                    */
   SCIP_Bool             binarybounds,       /**< should general integers get binary bounds [floor(.),ceil(.)]?  */
   SCIP_Bool             uselprows,          /**< should subproblem be created out of the rows in the LP rows?   */
   SCIP_Bool             usegcg              /**< should the subproblem be solved with GCG as well?              */
   );

#ifdef __cplusplus
}
#endif

#endif
