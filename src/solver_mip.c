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
/**@file   solver_mip.c
 * @brief  mip solver for pricing problems
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "solver_mip.h"
#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"
#include "type_solver.h"
#include "pricer_gcg.h"



#define SOLVER_NAME          "mip"
#define SOLVER_DESC          "mip solver for pricing problems"
#define SOLVER_PRIORITY      1000



/** branching data for branching decisions */
struct GCG_SolverData
{


};


/*
 * Callback methods for pricing problem solver
 */

static
GCG_DECL_SOLVERSOLVE(solverSolveMip)
{
#if 0

#ifdef DEBUG_PRICING_ALL_OUTPUT
   //if( pricetype == GCG_PRICETYPE_REDCOST )
   {
      char probname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "pricingmip_%d_%d_vars.lp", prob, SCIPgetNVars(scip));
      SCIP_CALL( SCIPwriteOrigProblem(pricingprob, probname, NULL, FALSE) );
      
      SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", SCIP_VERBLEVEL_HIGH) );
      
      SCIP_CALL( SCIPwriteParams(pricingprob, "pricing.set", TRUE, TRUE) );
      
   }
#endif

   /* start clock measuring the transformation effort */
   SCIP_CALL( SCIPstartClock(scip, solverdata->transformclock) );
   if( !solverdata->useheurpricing )
   {
      SCIP_CALL( SCIPtransformProb(pricingprob) );
   }
   SCIP_CALL( SCIPstopClock(scip, solverdata->transformclock) );

   /* start clock measuring the presolving effort */
   if( pricetype == GCG_PRICETYPE_REDCOST )
      SCIP_CALL( SCIPstartClock(scip, solverdata->redcostpresolveclock) );
   else
      SCIP_CALL( SCIPstartClock(scip, solverdata->farkaspresolveclock) );

   /* presolve the pricing submip */
   if( !solverdata->useheurpricing )
   {
      SCIP_CALL( SCIPpresolve(pricingprob) );
   }

   /* stop clock measuring the presolving effort, start clock measuring the solving effort */
   if( pricetype == GCG_PRICETYPE_REDCOST )
   {
      SCIP_CALL( SCIPstopClock(scip, solverdata->redcostpresolveclock) );
      SCIP_CALL( SCIPstartClock(scip, solverdata->redcostsolveclock) );
   }
   else
   {
      SCIP_CALL( SCIPstopClock(scip, solverdata->farkaspresolveclock) );
      SCIP_CALL( SCIPstartClock(scip, solverdata->farkassolveclock) );
   }

   /* solve the pricing submip */
   SCIP_CALL( SCIPsolve(pricingprob) );

   /* stop clock measuring the solving effort */
   if( pricetype == GCG_PRICETYPE_REDCOST )
      SCIP_CALL( SCIPstopClock(scip, solverdata->redcostsolveclock) );
   else
      SCIP_CALL( SCIPstopClock(scip, solverdata->farkassolveclock) );

   solverdata->solvedsubmipsoptimal++;

   SCIP_CALL( SCIPstopClock(scip, solverdata->subsolveclock) );
   SCIP_CALL( SCIPstartClock(scip, solverdata->owneffortclock) );

#ifdef DEBUG_PRICING_ALL_OUTPUT
   if( pricetype == GCG_PRICETYPE_REDCOST )
   {
      SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", 0) );
      SCIP_CALL( SCIPprintStatistics(pricingprob, NULL) );
   }
#endif

   /* so far, the pricing problem should be solved to optimality */
   assert( SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_GAPLIMIT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT 
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFEASIBLE
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT );
   /** @todo handle userinterrupt: set result pointer (and lowerbound), handle other solution states */

   if( SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT )
   {
      *nsols = 0;
      *sols = NULL;
      bestredcostvalid = FALSE;
      if( result != NULL ) 
         *result = SCIP_DIDNOTRUN;
   } 
   else
   {
      *nsols = SCIPgetNSols(pricingprob);
      *sols = SCIPgetSols(pricingprob);
      //printf("Pricingprob %d has found %d sols!\n", prob, nsols);
   }
#endif

   //printf("solver MIP\n");

   SCIP_CALL( SCIPtransformProb(pricingprob) );
   SCIP_CALL( SCIPpresolve(pricingprob) );
   /* solve the pricing submip */
   SCIP_CALL( SCIPsolve(pricingprob) );
   
   /* so far, the pricing problem should be solved to optimality */
   assert( SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_GAPLIMIT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT 
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFEASIBLE
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT );
   /** @todo handle userinterrupt: set result pointer (and lowerbound), handle other solution states */
   
   if( SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT )
   {
      *nsols = 0;
      *sols = NULL;
      *result = SCIP_DIDNOTRUN;
   } 
   else
   {
      *nsols = SCIPgetNSols(pricingprob);
      *sols = SCIPgetSols(pricingprob);
      *result = SCIP_OKAY;
      //printf("Pricingprob %d has found %d sols!\n", prob, nsols);
   }

#ifndef NDEBUG
   SCIP_Bool feasible;
   int i;
   for( i = 0; i < *nsols; i++ )
   {
      SCIP_CALL( SCIPcheckSolOrig(pricingprob, (*sols)[i], &feasible, TRUE, TRUE) );
      assert(feasible);
   }
#endif

#ifdef DEBUG_PRICING_ALL_OUTPUT
   if( pricetype == GCG_PRICETYPE_REDCOST )
   {
      SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", 0) );
      SCIP_CALL( SCIPprintStatistics(pricingprob, NULL) );
   }
#endif


   //SCIP_CALL( SCIPfreeTransform(pricingprob) );

   return SCIP_OKAY;
}


/** creates the most infeasible LP braching rule and includes it in SCIP */
SCIP_RETCODE GCGincludeSolverMip(
   SCIP*                 scip                /**< SCIP data structure */
   )
{   
   GCG_SOLVERDATA* data;

   data = NULL;
   
   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY, solverSolveMip, data) );

   return SCIP_OKAY;
}
