/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
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

/**@file   heur_ipcolgen.h
 * @brief  LNS heuristic that generates columns that are more likely to be in integer feasible solutions
 * @author Stephen J. Maher
 * @ingroup PRIMALHEURISTICS
 *
 * This heuristic is detailed in the paper:
 *
 * Maher, S. J. & Roennberg, E. Integer programming column generation: Accelerating branch-and-price using a novel
 * pricing scheme for finding high-quality solutions in set covering, packing, and partitioning problems. Mathematical
 * Programming Computation, under review.
 *
 * The heuristic is based on the destroy-and-repair framework for large neighbourhood search heuristics. The main
 * components of this heuristic are:
 * - A partial solution: this is found by taking the current incumbent and then fixing a subset of active variables to
 *   zero. This constitutes the destroy method of the heuristic.
 * - A modified pricing problem for generating repair columns
 * - A repair sub-MIP that is solved to find improving solutions given the set of repair columns.
 *
 * The pricing problem is modified by scaling the dual value in the objective function and adding dynamic penalties to
 * selected set covering, packing or partitioning constraints. The penalties are computed based on the activities of the
 * set covering, packing and partitioning constraints.
 *
 * During the procedure for generating repair columns, the dynamic penalties are updated with respect to the generated
 * columns and the activities in the penalised constraints. The update depends on the type of constraint. The update
 * procedure makes used of a pricing callback plugin. This plugin is executed before and after each pricing round to
 * update the dynamic penalties.
 *
 * The repair sub-MIP is constructed with the column from the restricted master problem and all repair columns. The
 * columns active in the partial solution are fixed in the repair problem. The repair problem is then solved to find a
 * configuration of the repair columns to repair the partial solution.
 *
 * This heuristic is very computationally expensive. It is only useful for problems where there is a large optimality
 * gap and it is difficult to find good primal feasible solutions. In this regard, a parameter has been added
 * (mininitialgap) that is used to set the minimum gap for calling the IPColGen heuristic. The first call to the
 * IPColGen heuristic occurs once the tail-off effect has been observed. This is computed by comparing the most recent
 * LP solution to the average LP solution from previous \f$n\f$ iterations. If the relative difference is less than
 * 0.01%, then the heuristic may be called (depending on other parameter settings).
 *
 * One of the most sensitive parameters, with respect to the performance of the heuristic, is the number of iterations
 * of the destroy-and-repair procedure without improvement in the incumbent solution. This is set to 3 by default.
 * Experiments conducted for the paper found that any deviation from this value results in a significant degradation in
 * the performance.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_HEUR_IPCOLGEN_H__
#define GCG_HEUR_IPCOLGEN_H__

#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

enum IPC_PenaltyType
{
   PENALTYTYPE_BIGM      = 1,
   PENALTYTYPE_SETCOVER  = 2,
   PENALTYTYPE_SETPACK   = 3
};
typedef enum IPC_PenaltyType IPC_PENALTYTYPE;

/** creates IP Column Generation primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurIPcolgen(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
