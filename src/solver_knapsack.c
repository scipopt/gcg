/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2012 Operations Research, RWTH Aachen University       */
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

/** knapsack pricing solverdata */
struct GCG_SolverData
{
   SCIP*                 origprob;           /**< original problem */
};


/*
 * Callback methods for pricing problem solver
 */

/** free method of knapsack solver */
static
GCG_DECL_SOLVERFREE(solverFreeKnapsack)
{
   GCG_SOLVERDATA* solverdata;

   assert(scip != NULL);
   assert(solver != NULL);

   solverdata = GCGpricerGetSolverdata(scip, solver);
   assert(solverdata != NULL);

   SCIPfreeMemory(scip, &solverdata);

   GCGpricerSetSolverdata(scip, solver, NULL);

   return SCIP_OKAY;
}

/** initialization method of knapsack solver */
static
GCG_DECL_SOLVERINITSOL(solverInitsolKnapsack)
{
   return SCIP_OKAY;
}

static
GCG_DECL_SOLVEREXITSOL(solverExitsolKnapsack)
{
   return SCIP_OKAY;
}

#define solverInitKnapsack NULL
#define solverExitKnapsack NULL

/** exact solving method for knapsack solver */
static
GCG_DECL_SOLVERSOLVE(solverSolveKnapsack)
{  /*lint --e{715}*/
   GCG_SOLVERDATA* solverdata;
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
   SCIP_Bool success;
   int i;
   int k;

   assert(pricingprob != NULL);
   assert(scip != NULL);
   assert(result != NULL);
   assert(solver != NULL);

   solverdata = GCGpricerGetSolverdata(scip, solver);
   assert(solverdata != NULL);

   pricingprobvars = SCIPgetVars(pricingprob);
   npricingprobvars = SCIPgetNVars(pricingprob);

   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvars, npricingprobvars) );
   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvals, npricingprobvars) );

   nconss = SCIPgetNConss(pricingprob);
   if( nconss != 1 )
   {
      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   cons = SCIPgetConss(pricingprob)[0];
   assert(cons != NULL);

   nitems = 0;

   for( i = 0; i < npricingprobvars; i++ )
   {
      if( SCIPvarGetUbLocal(pricingprobvars[i]) > SCIPvarGetLbLocal(pricingprobvars[i]) + 0.5 )
         nitems++;
   }

   if( !SCIPisIntegral(scip, SCIPgetRhsLinear(pricingprob, cons)) ||
      !SCIPisInfinity(scip, - SCIPgetLhsLinear(pricingprob, cons)) )
   {
      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   capacity = (SCIP_Longint)SCIPfloor(scip, SCIPgetRhsLinear(pricingprob, cons));
   consvars = SCIPgetVarsLinear(pricingprob, cons);
   nconsvars = SCIPgetNVarsLinear(pricingprob, cons);
   consvals = SCIPgetValsLinear(pricingprob, cons);

   for( i = 0; i < nconsvars; i++ )
   {
      if( !SCIPisIntegral(scip, consvals[i]) )
      {
         *result = SCIP_STATUS_UNKNOWN;
         return SCIP_OKAY;
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &items, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &weights, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &profits, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solitems, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsolitems, nitems) );

   BMSclearMemoryArray(weights, nitems);

   k = 0;
   for( i = 0; i < npricingprobvars; i++ )
   {
      if( SCIPvarGetUbLocal(pricingprobvars[i]) > SCIPvarGetLbLocal(pricingprobvars[i]) + 0.5 )
      {
         items[k] = i;
         profits[k] = - SCIPvarGetObj(pricingprobvars[i]);
         k++;
      }
   }
   assert(k == nitems);

   for( i = 0; i < nconsvars; i++ )
   {
      assert(SCIPisIntegral(scip, consvals[i]));

      if( SCIPisEQ(scip, SCIPvarGetUbLocal(consvars[i]), 0.0) )
         continue;
      if( SCIPisEQ(scip, SCIPvarGetLbLocal(consvars[i]), 1.0) )
      {
         capacity -= (SCIP_Longint)SCIPfloor(scip, consvals[i]);
         continue;
      }
      for( k = 0; k < nitems; k++ )
      {
         if( pricingprobvars[items[k]] == consvars[i] )
         {
            if( SCIPisPositive(scip, consvals[i]) )
            {
               weights[k] = (SCIP_Longint)SCIPfloor(scip, consvals[i]);
               break;
            }
            else
            {
               capacity -= (SCIP_Longint)SCIPfloor(scip, consvals[i]);
               weights[k] = (SCIP_Longint)SCIPfloor(scip, -1.0*consvals[i]);
               profits[k] *= -1.0;

               break;
            }
         }
      }
      assert(k < nitems);
   }

   /* solve knapsack problem exactly, all result pointers are needed! */
   SCIP_CALL( SCIPsolveKnapsackExactly(pricingprob, nitems, weights, profits, capacity, items, solitems,
         nonsolitems, &nsolitems, &nnonsolitems, &solval, &success ));

   assert(success);
   SCIPdebugMessage("knapsack solved, solval = %g\n", solval);

   nsolvars = 0;
   solisray[0] = FALSE;

   for( i = 0; i < nsolitems; i++ )
   {
      if( !SCIPisNegative(scip, consvals[solitems[i]]) )
      {
         solvars[nsolvars] = pricingprobvars[solitems[i]];
         solvals[nsolvars] = 1;
         nsolvars++;
      }
   }

   for( i = 0; i < nnonsolitems; i++ )
   {
      if( SCIPisNegative(scip, consvals[nonsolitems[i]]) )
      {
         solvars[nsolvars] = pricingprobvars[nonsolitems[i]];
         solvals[nsolvars] = 1;
         nsolvars++;
      }
   }

   for( i = 0; i < npricingprobvars; i++ )
   {
      if( SCIPvarGetLbLocal(pricingprobvars[i]) > 0.5 )
      {
         solvars[nsolvars] = pricingprobvars[i];
         solvals[nsolvars] = 1;
         nsolvars++;
      }
   }

   SCIP_CALL( SCIPcreateSol(pricingprob, &sols[0], NULL) );
   SCIP_CALL( SCIPsetSolVals(pricingprob, sols[0], nsolvars, solvars, solvals) );

   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &profits);
   SCIPfreeBufferArray(scip, &weights);
   SCIPfreeBufferArray(scip, &items);
   SCIPfreeMemoryArray(pricingprob, &solvars);
   SCIPfreeMemoryArray(pricingprob, &solvals);

   *nsols = 1;

   *lowerbound = solval;

   *result = SCIP_STATUS_OPTIMAL;

   return SCIP_OKAY;
}


/** heuristic solving method of knapsack solver */
static
GCG_DECL_SOLVERSOLVEHEUR(solverSolveHeurKnapsack)
{  /*lint --e{715}*/

   GCG_SOLVERDATA* solverdata;
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

   int i;
   int k;

   assert(pricingprob != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   solverdata = GCGpricerGetSolverdata(scip, solver);
   assert(solverdata != NULL);

   pricingprobvars = SCIPgetVars(pricingprob);
   npricingprobvars = SCIPgetNVars(pricingprob);

   nconss = SCIPgetNConss(pricingprob);
   if( nconss != 1 )
   {
      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   cons = SCIPgetConss(pricingprob)[0];
   assert(cons != NULL);

   nitems = 0;

   for( i = 0; i < npricingprobvars; i++ )
   {
      if( SCIPvarGetUbLocal(pricingprobvars[i]) > SCIPvarGetLbLocal(pricingprobvars[i]) + 0.5 )
         nitems++;
   }

   if( !SCIPisIntegral(scip, SCIPgetRhsLinear(pricingprob, cons)) ||
      !SCIPisInfinity(scip, - SCIPgetLhsLinear(pricingprob, cons)) )
   {
      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   capacity = (SCIP_Longint)SCIPfloor(scip, SCIPgetRhsLinear(pricingprob, cons));
   consvars = SCIPgetVarsLinear(pricingprob, cons);
   nconsvars = SCIPgetNVarsLinear(pricingprob, cons);
   consvals = SCIPgetValsLinear(pricingprob, cons);

   for( i = 0; i < nconsvars; i++ )
   {
      if( !SCIPisIntegral(scip, consvals[i]) )
      {
         *result = SCIP_STATUS_UNKNOWN;
         return SCIP_OKAY;
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &items, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &weights, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &profits, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solitems, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsolitems, nitems) );

   BMSclearMemoryArray(weights, nitems);

   k = 0;
   for( i = 0; i < npricingprobvars; i++ )
   {
      if( SCIPvarGetUbLocal(pricingprobvars[i]) > SCIPvarGetLbLocal(pricingprobvars[i]) + 0.5 )
      {
         items[k] = i;
         profits[k] = - SCIPvarGetObj(pricingprobvars[i]);
         k++;
      }
   }
   assert(k == nitems);

   for( i = 0; i < nconsvars; i++ )
   {
      assert(SCIPisIntegral(scip, consvals[i]));

      if( SCIPisEQ(scip, SCIPvarGetUbLocal(consvars[i]), 0.0) )
         continue;
      if( SCIPisEQ(scip, SCIPvarGetLbLocal(consvars[i]), 1.0) )
      {
         capacity -= (SCIP_Longint)SCIPfloor(scip, consvals[i]);
         continue;
      }
      for( k = 0; k < nitems; k++ )
      {
         if( pricingprobvars[items[k]] == consvars[i] )
         {
            if( SCIPisPositive(scip, consvals[i]) )
            {
               weights[k] = (SCIP_Longint)SCIPfloor(scip, consvals[i]);
               break;
            }
            else
            {
               capacity -= (SCIP_Longint)SCIPfloor(scip, consvals[i]);
               weights[k] = (SCIP_Longint)SCIPfloor(scip, -1.0*consvals[i]);
               profits[k] *= -1.0;

               break;
            }
         }
      }
      assert(k < nitems);
   }

   /* solve knapsack problem exactly, all result pointers are needed! */
   SCIP_CALL( SCIPsolveKnapsackApproximately(pricingprob, nitems, weights, profits, capacity, items, solitems,
         nonsolitems, &nsolitems, &nnonsolitems, &solval ));

   SCIPdebugMessage("knapsack solved, solval = %g\n", solval);

   nsolvars = 0;
   solisray[0] = FALSE;

   for( i = 0; i < nsolitems; i++ )
   {
      if( !SCIPisNegative(scip, consvals[solitems[i]]) )
      {
         solvars[nsolvars] = pricingprobvars[solitems[i]];
         solvals[nsolvars] = 1;
         nsolvars++;
      }
   }

   for( i = 0; i < nnonsolitems; i++ )
   {
      if( SCIPisNegative(scip, consvals[nonsolitems[i]]) )
      {
         solvars[nsolvars] = pricingprobvars[nonsolitems[i]];
         solvals[nsolvars] = 1;
         nsolvars++;
      }
   }

   for( i = 0; i < npricingprobvars; i++ )
   {
      if( SCIPvarGetLbLocal(pricingprobvars[i]) > 0.5 )
      {
         solvars[nsolvars] = pricingprobvars[i];
         solvals[nsolvars] = 1;
         nsolvars++;
      }
   }

   SCIP_CALL( SCIPcreateSol(pricingprob, &sols[0], NULL) );
   SCIP_CALL( SCIPsetSolVals(pricingprob, sols[0], nsolvars, solvars, solvals) );

   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &profits);
   SCIPfreeBufferArray(scip, &weights);
   SCIPfreeBufferArray(scip, &items);
   SCIPfreeMemoryArray(pricingprob, &solvars);
   SCIPfreeMemoryArray(pricingprob, &solvals);

   *nsols = 1;

   *lowerbound = -SCIPinfinity(scip);

   *result = SCIP_STATUS_OPTIMAL;

   return SCIP_OKAY;
}


/** creates the knapsack solver for pricing problems and includes it in GCG */
SCIP_RETCODE GCGincludeSolverKnapsack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   GCG_SOLVERDATA* data;

   SCIP_CALL( SCIPallocMemory( scip, &data) );
   data->origprob = GCGpricerGetOrigprob(scip);

   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY, solverSolveKnapsack,
         solverSolveHeurKnapsack, solverFreeKnapsack, solverInitKnapsack, solverExitKnapsack,
         solverInitsolKnapsack, solverExitsolKnapsack, data) );

   return SCIP_OKAY;
}
