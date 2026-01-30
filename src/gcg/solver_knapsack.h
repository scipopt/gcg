/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   solver_knapsack.h
 * @brief  knapsack solver for pricing problems
 * @author Gerald Gamrath
 * @ingroup PRICINGSOLVERS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_SOLVER_KNAPSACK_H__
#define GCG_SOLVER_KNAPSACK_H__

#include "scip/scip.h"

#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the knapsack solver for pricing problems and includes it in GCG */
GCG_EXPORT
SCIP_RETCODE GCGincludeSolverKnapsack(
   GCG*                  gcg                 /**< GCG data structure */
   );

GCG_EXPORT
SCIP_RETCODE GCGsolverKnapsackSolveKnapsack(
   SCIP_Bool             exactly,            /**< should the pricing problem be solved to optimality or heuristically? */
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   SCIP_Real*            solval,             /**< pointer to store the solution objective value */
   GCG_PRICINGSTATUS*    status,             /**< pointer to store pricing problem status */
   SCIP_VAR***           solvars,            /**< pointer to store solution vars (will be allocated as buffer array and must be freed afterwards) */
   SCIP_Real**           solvals,            /**< pointer to store solution vals (will be allocated as buffer array and must be freed afterwards) */
   int*                  nsolvars            /**< pointer to store the number of solution vars */
   );

#ifdef __cplusplus
}
#endif

#endif
