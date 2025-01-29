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
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
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
#include "def.h"
#include "type_pricingstatus.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the knapsack solver for pricing problems and includes it in GCG */
GCG_EXPORT
SCIP_RETCODE GCGincludeSolverKnapsack(
   SCIP*                 scip                /**< SCIP data structure */
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
