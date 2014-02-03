/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2014 Operations Research, RWTH Aachen University       */
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

/**@file   solver_knapsack.c
 * @brief  knapsack solver for pricing problems
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "solver_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"
#include "type_solver.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"

#define SOLVER_NAME          "knapsack"
#define SOLVER_DESC          "knapsack solver for pricing problems"
#define SOLVER_PRIORITY      -100

/** knapsack pricing solver needs no solverdata */
/* struct GCG_SolverData {}; */


/*
 * Local methods
 */

/** solve the pricing problem as a knapsack problem, either exactly or approximately */
static
SCIP_RETCODE solveKnapsack(
   SCIP_Bool             exactly,            /**< should the pricing problem be solved to optimality or heuristically? */
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   GCG_SOLVER*           solver,             /**< solver data structure */
   int                   probnr,             /**< problem number */
   SCIP_Real*            lowerbound,         /**< pointer to store lower bound */
   SCIP_SOL**            sols,               /**< array of solutions */
   SCIP_Bool*            solisray,           /**< array indicating whether solution is a ray */
   int                   maxsols,            /**< size of preallocated array */
   int*                  nsols,              /**< pointer to store number of solutions */
   SCIP_STATUS*          result              /**< pointer to store pricing problem status */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR** consvars;
   int nconsvars;
   SCIP_Real* consvals;
   SCIP_VAR** solvars;
   SCIP_Real* solvals;
   int nsolvars;
   SCIP_VAR** pricingprobvars;
   int npricingprobvars;
   int nconss;

   int                   nitems;
   SCIP_Longint*         weights;
   SCIP_Real*            profits;
   SCIP_Longint          capacity;
   int*                  items;
   int*                  solitems;
   int                   nsolitems;
   int*                  nonsolitems;
   int                   nnonsolitems;
   SCIP_Real             solval;
   SCIP_Bool             success;

   int i;
   int j;
   int k;

   /* check preconditions */
   assert(pricingprob != NULL);
   assert(solver != NULL);
   assert(lowerbound != NULL);
   assert(sols != NULL);
   assert(solisray != NULL);
   assert(nsols != NULL);
   assert(result != NULL);

   pricingprobvars = SCIPgetVars(pricingprob);
   npricingprobvars = SCIPgetNVars(pricingprob);

   /* check prerequisites: the pricing problem can be solved as a knapsack problem only if
    * - all variables are nonnegative integer variables
    * - there is only one constraint, which has infinite lhs and integer rhs
    */
   if( SCIPgetNBinVars(pricingprob) + SCIPgetNIntVars(pricingprob) < npricingprobvars )
   {
      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }
   for( i = SCIPgetNBinVars(pricingprob); i < SCIPgetNBinVars(pricingprob) + SCIPgetNIntVars(pricingprob); ++i )
   {
      if( SCIPisNegative(pricingprob, SCIPvarGetLbLocal(pricingprobvars[i])) )
      {
         *result = SCIP_STATUS_UNKNOWN;
         return SCIP_OKAY;
      }
   }

   nconss = SCIPgetNConss(pricingprob);
   if( nconss != 1 )
   {
      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   cons = SCIPgetConss(pricingprob)[0];
   assert(cons != NULL);

   if( !SCIPisIntegral(pricingprob, SCIPgetRhsLinear(pricingprob, cons)) ||
      !SCIPisInfinity(pricingprob, - SCIPgetLhsLinear(pricingprob, cons)) )
   {
      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   consvars = SCIPgetVarsLinear(pricingprob, cons);
   nconsvars = SCIPgetNVarsLinear(pricingprob, cons);
   consvals = SCIPgetValsLinear(pricingprob, cons);

   for( i = 0; i < nconsvars; i++ )
   {
      if( !SCIPisIntegral(pricingprob, consvals[i]) )
      {
         *result = SCIP_STATUS_UNKNOWN;
         return SCIP_OKAY;
      }
   }

   capacity = (SCIP_Longint)SCIPfloor(pricingprob, SCIPgetRhsLinear(pricingprob, cons));
   nitems = 0;
   for( i = 0; i < npricingprobvars; i++ )
      nitems += (int)(SCIPvarGetUbLocal(pricingprobvars[i]) - SCIPvarGetLbLocal(pricingprobvars[i]) + 0.5);

   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvars, npricingprobvars) );
   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvals, npricingprobvars) );

   SCIP_CALL( SCIPallocBufferArray(pricingprob, &items, nitems) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &weights, nitems) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &profits, nitems) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &solitems, nitems) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &nonsolitems, nitems) );

   BMSclearMemoryArray(weights, nitems);

   k = 0;
   for( i = 0; i < npricingprobvars; i++ )
   {
      assert(!SCIPisInfinity(pricingprob, SCIPvarGetUbLocal(pricingprobvars[i])));
      for( j = 0; j < (int)(SCIPvarGetUbLocal(pricingprobvars[i]) - SCIPvarGetLbLocal(pricingprobvars[i]) + 0.5); ++j )
      {
         items[k] = i;
         profits[k] = - SCIPvarGetObj(pricingprobvars[i]);
         k++;
      }
   }
   assert(k == nitems);

   for( i = 0; i < nconsvars; i++ )
   {
      assert(SCIPisIntegral(pricingprob, consvals[i]));

      if( SCIPisEQ(pricingprob, SCIPvarGetUbLocal(consvars[i]), 0.0) )
         continue;
      if( SCIPisGE(pricingprob, SCIPvarGetLbLocal(consvars[i]), 1.0) )
      {
         capacity -= (SCIP_Longint)SCIPfloor(pricingprob, SCIPvarGetLbLocal(consvars[i]))
            * (SCIP_Longint)SCIPfloor(pricingprob, consvals[i]);
      }

      j = 0;
      for( k = 0; k < nitems && j < (int)(SCIPvarGetUbLocal(consvars[i]) - SCIPvarGetLbLocal(consvars[i]) + 0.5); k++ )
      {
         if( pricingprobvars[items[k]] == consvars[i] )
         {
            if( SCIPisPositive(pricingprob, consvals[i]) )
            {
               weights[k] = (SCIP_Longint)SCIPfloor(pricingprob, consvals[i]);
            }
            else
            {
               capacity -= (SCIP_Longint)SCIPfloor(pricingprob, consvals[i]);
               weights[k] = (SCIP_Longint)SCIPfloor(pricingprob, -1.0*consvals[i]);
               profits[k] *= -1.0;
            }
            ++j;
         }
      }
      assert(j == (int)(SCIPvarGetUbLocal(consvars[i]) - SCIPvarGetLbLocal(consvars[i]) + 0.5));
   }

   success = TRUE;

   /* problem is infeasible */
   if( capacity < 0 )
   {
      *result = SCIP_STATUS_INFEASIBLE;
      goto TERMINATE;
   }

   /* solve knapsack problem, all result pointers are needed! */
   if( exactly )
   {
      SCIP_CALL( SCIPsolveKnapsackExactly(pricingprob, nitems, weights, profits, capacity, items, solitems,
         nonsolitems, &nsolitems, &nnonsolitems, &solval, &success ));
   }
   else
   {
      SCIP_CALL( SCIPsolveKnapsackApproximately(pricingprob, nitems, weights, profits, capacity, items, solitems,
         nonsolitems, &nsolitems, &nnonsolitems, &solval ));
   }

   assert(success);
   SCIPdebugMessage("knapsack solved, solval = %g\n", solval);

   nsolvars = 0;
   solisray[0] = FALSE;

   for( i = 0; i < nsolitems; i++ )
   {
      if( !SCIPisNegative(pricingprob, consvals[solitems[i]]) )
      {
         for( j = 0; j < nsolvars; ++j )
            if( solvars[j] == pricingprobvars[solitems[i]] )
               break;

         if( j == nsolvars )
         {
            solvars[j] = pricingprobvars[solitems[i]];
            solvals[j] = 1.0;
            ++nsolvars;
         }
         else
            solvals[j] += 1.0;
      }
   }

   for( i = 0; i < nnonsolitems; i++ )
   {
      if( SCIPisNegative(pricingprob, consvals[nonsolitems[i]]) )
      {
         for( j = 0; j < nsolvars; ++j )
            if( solvars[j] == pricingprobvars[nonsolitems[i]] )
               break;

         if( j == nsolvars )
         {
            solvars[j] = pricingprobvars[nonsolitems[i]];
            solvals[j] = 1.0;
            ++nsolvars;
         }
         else
            solvals[j] += 1.0;
      }
   }

   for( i = 0; i < npricingprobvars; i++ )
   {
      if( SCIPisGE(pricingprob, SCIPvarGetLbLocal(pricingprobvars[i]), 1.0) )
      {
         for( j = 0; j < nsolvars; ++j )
            if( solvars[j] == pricingprobvars[i] )
               break;

         if( j == nsolvars )
         {
            solvars[j] = pricingprobvars[i];
            solvals[j] = SCIPfloor(pricingprob, SCIPvarGetLbLocal(pricingprobvars[i]));
            ++nsolvars;
         }
         else
            solvals[j] += SCIPfloor(pricingprob, SCIPvarGetLbLocal(pricingprobvars[i]));
      }
   }

   SCIP_CALL( SCIPcreateSol(pricingprob, &sols[0], NULL) );
   SCIP_CALL( SCIPsetSolVals(pricingprob, sols[0], nsolvars, solvars, solvals) );

   *nsols = 1;

   *lowerbound = solval;

   *result = SCIP_STATUS_OPTIMAL;

 TERMINATE:
   SCIPfreeBufferArray(pricingprob, &nonsolitems);
   SCIPfreeBufferArray(pricingprob, &solitems);
   SCIPfreeBufferArray(pricingprob, &profits);
   SCIPfreeBufferArray(pricingprob, &weights);
   SCIPfreeBufferArray(pricingprob, &items);
   SCIPfreeMemoryArray(pricingprob, &solvars);
   SCIPfreeMemoryArray(pricingprob, &solvals);

   return SCIP_OKAY;
}

/*
 * Callback methods for pricing problem solver
 */

#define solverFreeKnapsack NULL
#define solverInitsolKnapsack NULL
#define solverExitsolKnapsack NULL
#define solverInitKnapsack NULL
#define solverExitKnapsack NULL

/** exact solving method for knapsack solver */
static
GCG_DECL_SOLVERSOLVE(solverSolveKnapsack)
{  /*lint --e{715}*/

   /* solve the knapsack problem exactly */
   SCIP_CALL( solveKnapsack(TRUE, pricingprob, solver, probnr, lowerbound, sols, solisray, maxsols, nsols, result) );

   return SCIP_OKAY;
}


/** heuristic solving method of knapsack solver */
static
GCG_DECL_SOLVERSOLVEHEUR(solverSolveHeurKnapsack)
{  /*lint --e{715}*/

   /* solve the knapsack problem approximately */
   SCIP_CALL( solveKnapsack(FALSE, pricingprob, solver, probnr, lowerbound, sols, solisray, maxsols, nsols, result) );

   return SCIP_OKAY;
}


/** creates the knapsack solver for pricing problems and includes it in GCG */
SCIP_RETCODE GCGincludeSolverKnapsack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY, solverSolveKnapsack,
         solverSolveHeurKnapsack, solverFreeKnapsack, solverInitKnapsack, solverExitKnapsack,
         solverInitsolKnapsack, solverExitsolKnapsack, NULL) );

   return SCIP_OKAY;
}
