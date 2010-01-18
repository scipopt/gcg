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
#define SOLVER_PRIORITY      0



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

#if 0

   /* start clock measuring the presolving effort */
   if( pricetype == GCG_PRICETYPE_REDCOST )
      SCIP_CALL( SCIPstartClock(scip, solverdata->redcostpresolveclock) );
   else
      SCIP_CALL( SCIPstartClock(scip, solverdata->farkaspresolveclock) );


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

   /* stop clock measuring the solving effort */
   if( pricetype == GCG_PRICETYPE_REDCOST )
      SCIP_CALL( SCIPstopClock(scip, solverdata->redcostsolveclock) );
   else
      SCIP_CALL( SCIPstopClock(scip, solverdata->farkassolveclock) );

#endif

   SCIP_CALL( SCIPtransformProb(pricingprob) );

   /* presolve the pricing submip */
   if( SCIPgetStage(pricingprob) < SCIP_STAGE_PRESOLVING )
   {
      SCIP_CALL( SCIPpresolve(pricingprob) );
   }

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
      *result = SCIP_STATUS_UNKNOWN;
   } 
   else
   {
      *result = SCIP_STATUS_OPTIMAL;
      //printf("Pricingprob %d has found %d sols!\n", prob, nsols);
   }

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


static
GCG_DECL_SOLVERSOLVEHEUR(solverSolveHeurMip)
{

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

#if 0

   /* start clock measuring the presolving effort */
   if( pricetype == GCG_PRICETYPE_REDCOST )
      SCIP_CALL( SCIPstartClock(scip, solverdata->redcostpresolveclock) );
   else
      SCIP_CALL( SCIPstartClock(scip, solverdata->farkaspresolveclock) );


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

   /* stop clock measuring the solving effort */
   if( pricetype == GCG_PRICETYPE_REDCOST )
      SCIP_CALL( SCIPstopClock(scip, solverdata->redcostsolveclock) );
   else
      SCIP_CALL( SCIPstopClock(scip, solverdata->farkassolveclock) );

#endif

   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/stallnodes", 100) );
   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/nodes", 1000) );
   SCIP_CALL( SCIPsetRealParam(pricingprob, "limits/gap", 0.2) );
   //SCIP_CALL( SCIPsetIntParam(pricingprob, "limits/bestsol", 5) ); 


   SCIP_CALL( SCIPtransformProb(pricingprob) );

   /* presolve the pricing submip */
   if( SCIPgetStage(pricingprob) < SCIP_STAGE_PRESOLVING )
   {
      SCIP_CALL( SCIPpresolve(pricingprob) );
   }

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
      *result = SCIP_STATUS_UNKNOWN;
   } 
   else
   {
      *result = SCIP_STATUS_OPTIMAL;
      //printf("Pricingprob %d has found %d sols!\n", prob, nsols);
   }

#ifdef DEBUG_PRICING_ALL_OUTPUT
   if( pricetype == GCG_PRICETYPE_REDCOST )
   {
      SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", 0) );
      SCIP_CALL( SCIPprintStatistics(pricingprob, NULL) );
   }
#endif

   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/stallnodes", -1) );
   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/nodes", -1) ); 
   SCIP_CALL( SCIPsetRealParam(pricingprob, "limits/gap", 0.0) ); 
   SCIP_CALL( SCIPsetIntParam(pricingprob, "limits/bestsol", -1) ); 


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
   
   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY, 
         solverSolveMip, solverSolveHeurMip, data) );

   return SCIP_OKAY;
}
