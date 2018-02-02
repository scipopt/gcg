/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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
#include "pub_gcgcol.h"

#define SOLVER_NAME          "knapsack"
#define SOLVER_DESC          "knapsack solver for pricing problems"
#define SOLVER_PRIORITY      200

#define SOLVER_ENABLED       TRUE  /**< indicates whether the solver should be enabled */

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
   GCG_COL**             cols,               /**< array of columns corresponding to solutions */
   int                   maxcols,            /**< size of preallocated array */
   int*                  ncols,              /**< pointer to store number of columns */
   SCIP_STATUS*          result              /**< pointer to store pricing problem status */
   )
{ /*lint -e715 */
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
   SCIP_Real*            ubs;
   SCIP_Longint          capacity;
   SCIP_Longint          prelcapacity;
   int*                  items;
   int*                  solitems;
   int                   nsolitems;
   int*                  nonsolitems;
   int                   nnonsolitems;
   SCIP_Real             solval;
   SCIP_Bool             success;
   SCIP_Bool             inferbounds;
   int i;
   int j;
   int k;

   /* check preconditions */
   assert(pricingprob != NULL);
   assert(solver != NULL);
   assert(lowerbound != NULL);
   assert(cols != NULL);
   assert(ncols != NULL);
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
   prelcapacity = capacity;

   inferbounds = FALSE;
   for( i = 0; i < nconsvars; i++ )
   {
      if( SCIPisInfinity(pricingprob, SCIPvarGetUbLocal(consvars[i])) )
      {
         inferbounds = TRUE;
      }
      if( SCIPisNegative(pricingprob, consvals[i]) )
      {
         /* Handle the cases where the transformation is not clear:
          *
          * a column with infinite upper bound (capacity not deducible) or
          * a column column with negative weight and negative cost (should we add it?)
          */
         if( SCIPisInfinity(pricingprob, SCIPvarGetUbLocal(consvars[i])) )
         {
            *result = SCIP_STATUS_UNKNOWN;
            return SCIP_OKAY;
         }
         else if( SCIPisNegative(pricingprob, SCIPvarGetObj(consvars[i])) )
         {
            *result = SCIP_STATUS_UNKNOWN;
            return SCIP_OKAY;
         }

         prelcapacity -= (SCIP_Longint)SCIPfloor(pricingprob, consvals[i]* SCIPvarGetUbLocal(consvars[i]));
      }
   }
   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &ubs, npricingprobvars) );

   for( i = 0; i < nconsvars; i++ )
   {
      if( inferbounds && SCIPisInfinity(pricingprob, SCIPvarGetUbLocal(consvars[i])) )
      {
         SCIP_Real newbound = SCIPfloor(pricingprob, ABS((SCIP_Real)prelcapacity/consvals[i]));
         SCIPdebugMessage("newbound: %.2f/%.2f = %.2f\n", (SCIP_Real)prelcapacity, consvals[i], newbound);
         ubs[i] = newbound;
      }
      else
         ubs[i] = SCIPvarGetUbLocal(consvars[i]);

   }

   nitems = 0;
   for( i = 0; i < nconsvars; i++ )
   {
      assert(!SCIPisInfinity(pricingprob, ubs[i]));
      SCIPdebugMessage("%d: %d+%d\n",i, nitems,  (int)(ubs[i] - SCIPvarGetLbLocal(consvars[i]) + 0.5));
      nitems += (int)(ubs[i] - SCIPvarGetLbLocal(consvars[i]) + 0.5);
   }

   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvars, npricingprobvars) );
   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvals, npricingprobvars) );

   SCIP_CALL( SCIPallocBufferArray(pricingprob, &items, nitems) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &weights, nitems) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &profits, nitems) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &solitems, nitems) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &nonsolitems, nitems) );

   BMSclearMemoryArray(weights, nitems);

   k = 0;
   for( i = 0; i < nconsvars; i++ )
   {
      assert(!SCIPisInfinity(pricingprob, ubs[i]));
      for( j = 0; j < (int)(ubs[i] - SCIPvarGetLbLocal(consvars[i]) + 0.5); ++j )
      {
         items[k] = i;
         profits[k] = - SCIPvarGetObj(consvars[i]);
         SCIPdebugMessage("%d: <%s> %d\n",k, SCIPvarGetName(consvars[i]), i);

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
   }

   for( k = 0; k < nitems; k++ )
   {
      i = items[k];
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

   if( ! success )
   {
      SCIPwarningMessage(pricingprob, "Knapsack solver could not solve pricing problem!");
      goto TERMINATE;
   }

   SCIPdebugMessage("knapsack solved, solval = %g\n", solval);

   nsolvars = 0;

   for( i = 0; i < nsolitems; i++ )
   {
      if( !SCIPisNegative(pricingprob, consvals[solitems[i]]) )
      {
         for( j = 0; j < nsolvars; ++j )
            if( solvars[j] == consvars[solitems[i]] )
               break;

         if( j == nsolvars )
         {
            solvars[j] = consvars[solitems[i]];
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
            if( solvars[j] == consvars[nonsolitems[i]] )
               break;

         if( j == nsolvars )
         {
            solvars[j] = consvars[nonsolitems[i]];
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

   SCIP_CALL( GCGcreateGcgCol(pricingprob, &cols[0], probnr, solvars, solvals, nsolvars, FALSE, SCIPinfinity(pricingprob)) );

   *ncols = 1;

   solval = 0.0;

   for( i = 0; i < nsolvars; ++i )
   {
      solval += solvals[i] * SCIPvarGetObj(solvars[i]);
   }

   *lowerbound = solval;

   *result = SCIP_STATUS_OPTIMAL;

 TERMINATE:
   SCIPfreeBufferArray(pricingprob, &nonsolitems);
   SCIPfreeBufferArray(pricingprob, &solitems);
   SCIPfreeBufferArray(pricingprob, &profits);
   SCIPfreeBufferArray(pricingprob, &weights);
   SCIPfreeBufferArray(pricingprob, &items);
   SCIPfreeMemoryArray(pricingprob, &ubs);
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
   SCIP_CALL( solveKnapsack(TRUE, pricingprob, solver, probnr, lowerbound, cols, maxcols, ncols, result) );

   return SCIP_OKAY;
}


/** heuristic solving method of knapsack solver */
static
GCG_DECL_SOLVERSOLVEHEUR(solverSolveHeurKnapsack)
{  /*lint --e{715}*/

   /* solve the knapsack problem approximately */
   SCIP_CALL( solveKnapsack(FALSE, pricingprob, solver, probnr, lowerbound, cols, maxcols, ncols, result) );

   return SCIP_OKAY;
}


/** creates the knapsack solver for pricing problems and includes it in GCG */
SCIP_RETCODE GCGincludeSolverKnapsack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY, SOLVER_ENABLED, solverSolveKnapsack,
         solverSolveHeurKnapsack, solverFreeKnapsack, solverInitKnapsack, solverExitKnapsack,
         solverInitsolKnapsack, solverExitsolKnapsack, NULL) );

   return SCIP_OKAY;
}
