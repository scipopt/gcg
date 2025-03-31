/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_masterlinesdiving.h
 * @ingroup DIVINGHEURISTICS
 * @brief  master LP diving heuristic that fixes variables with a large difference to their root solution
 * @author Tobias Achterberg
 *
 * Diving heuristic: Iteratively fixes some fractional variable and resolves the LP-relaxation, thereby simulating a
 * depth-first-search in the tree. Line search diving chooses the variable with the greatest difference of its root LP
 * solution and the current LP solution, hence, the variable that developed most.  It is fixed to the next integer in
 * the direction it developed. One-level backtracking is applied: If the LP gets infeasible, the last fixing is undone,
 * and the opposite fixing is tried. If this is infeasible, too, the procedure aborts.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_MASTERLINESDIVING_H__
#define __SCIP_HEUR_MASTERLINESDIVING_H__



#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the masterlinesdiving primal heuristic and includes it in GCG */
GCG_EXPORT
SCIP_RETCODE GCGincludeHeurMasterlinesdiving(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
