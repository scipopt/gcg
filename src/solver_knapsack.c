/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"
//#define SCIP_DEBUG
/**@file   solver_knapsack.c
 * @brief  knapsack solver for pricing problems
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "solver_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"
#include "type_solver.h"
#include "pricer_gcg.h"



#define SOLVER_NAME          "knapsack"
#define SOLVER_DESC          "knapsack solver for pricing problems"
#define SOLVER_PRIORITY      100



/** branching data for branching decisions */
struct GCG_SolverData
{
};


/*
 * Callback methods for pricing problem solver
 */

static
GCG_DECL_SOLVERSOLVE(solverSolveKnapsack)
{
   SCIP_CONS* cons;
   SCIP_VAR** consvars;
   int nconsvars;
   SCIP_Real* consvals;
   
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
   SCIP_Bool             applyable;      

   int i;
   int k;

   SCIP_SOL* sol;

   SCIP_Bool stored;

   assert(pricingprob != NULL);
   assert(scip != NULL);
   assert(result != NULL);
   assert(nsols != NULL);
   assert(sols != NULL);
   assert(*nsols > 0);

   //printf("solver knapsack\n");

   pricingprobvars = SCIPgetVars(pricingprob);
   npricingprobvars = SCIPgetNVars(pricingprob);

   nconss = SCIPgetNConss(pricingprob);
   if ( nconss != 1 )
   {
      //printf("%d conss in problem, abort!\n", nconss);

      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   //cons = SCIPconshdlrGetConss(linearhandler)[0];
   cons = SCIPgetConss(pricingprob)[0];
   assert(cons != NULL);
   
   nitems = 0;
   
   for( i = 0; i < npricingprobvars; i++ )
   {
      if( SCIPvarGetUbLocal(pricingprobvars[i]) > SCIPvarGetLbLocal(pricingprobvars[i]) + 0.5 )
         nitems++;
   }

   applyable = TRUE;
   if( !SCIPisIntegral(scip, SCIPgetRhsLinear(pricingprob, cons)) || 
      !SCIPisInfinity(scip, - SCIPgetLhsLinear(pricingprob, cons)) )
   {
      //printf("wrong structure, abort!\n");

      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   capacity = SCIPfloor(scip, SCIPgetRhsLinear(pricingprob, cons));
   consvars = SCIPgetVarsLinear(pricingprob, cons);
   nconsvars = SCIPgetNVarsLinear(pricingprob, cons);
   consvals = SCIPgetValsLinear(pricingprob, cons);

   for( i = 0; i < nconsvars; i++ )
   {
      if( !SCIPisIntegral(scip, consvals[i]) )
      {
         //printf("wrong structure, abort!\n");
         
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
         capacity -= SCIPfloor(scip, consvals[i]);
         continue;
      }
      for( k = 0; k < nitems; k++ )
      {
         if( pricingprobvars[items[k]] == consvars[i] )
         {
            if( SCIPisPositive(scip, consvals[i]) )
            {
               weights[k] = SCIPfloor(scip, consvals[i]);
               break;
            }
            else
            {
               capacity -= SCIPfloor(scip, consvals[i]);
               weights[k] = SCIPfloor(scip, -1.0*consvals[i]);
               profits[k] *= -1.0; 

               break;
            }
         }
      }
      assert(k < nitems);
   }




   /* solve knapsack problem exactly, all result pointers are needed! */
   SCIP_CALL( SCIPsolveKnapsackExactly(pricingprob, nitems, weights, profits, capacity, items, solitems, 
         nonsolitems, &nsolitems, &nnonsolitems, &solval ));

   //printf("knapsack solved, solval = %g\n", solval);

   SCIP_CALL( SCIPtransformProb(pricingprob) );

   SCIP_CALL( SCIPcreateSol( pricingprob, &sol, NULL) );

   for( i = 0; i < nsolitems; i++ )
   {
      if( !SCIPisNegative(scip, consvals[solitems[i]]) )
      {
         SCIP_CALL( SCIPsetSolVal(pricingprob, sol, pricingprobvars[solitems[i]], 1) );
      }
   }

   for( i = 0; i < nnonsolitems; i++ )
   {
      if( SCIPisNegative(scip, consvals[nonsolitems[i]]) )
      {
         SCIP_CALL( SCIPsetSolVal(pricingprob, sol, pricingprobvars[nonsolitems[i]], 1) );
      }
   }

   for( i = 0; i < npricingprobvars; i++ )
   {
      if( SCIPvarGetLbLocal(pricingprobvars[i]) > 0.5 )
      {
         SCIP_CALL( SCIPsetSolVal(pricingprob, sol, pricingprobvars[i], 1) );
      }
   }


   SCIP_CALL( SCIPaddSolFree(pricingprob, &sol, &stored) );
   assert(stored);

   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &profits);
   SCIPfreeBufferArray(scip, &weights);
   SCIPfreeBufferArray(scip, &items);

   *nsols = SCIPgetNSols(pricingprob);
   *sols = SCIPgetSols(pricingprob);

#ifndef NDEBUG
   SCIP_Bool feasible;
   for( i = 0; i < *nsols; i++ )
   {
      SCIP_CALL( SCIPcheckSolOrig(pricingprob, (*sols)[i], &feasible, TRUE, TRUE) );
      //SCIP_CALL( SCIPprintSol(pricingprob, (*sols)[i], NULL, FALSE) );
      assert(feasible);
   }
#endif

   *result = SCIP_STATUS_OPTIMAL;

   return SCIP_OKAY;
}


/** creates the most infeasible LP braching rule and includes it in SCIP */
SCIP_RETCODE GCGincludeSolverKnapsack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{   
   GCG_SOLVERDATA* data;

   data = NULL;
   
   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY, solverSolveKnapsack, data) );

   return SCIP_OKAY;
}
