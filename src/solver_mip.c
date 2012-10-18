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

/**@file   solver_mip.c
 * @brief  pricing solver solving the pricing problem as a sub-MIP, using SCIP
 * @author Gerald Gamrath
 * @author Martin Bergner
 *
 * \bug
 * there are some issues with unstable or infinite pricing problem solutions
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define DEBUG_PRICING_ALL_OUTPUT*/

#include <assert.h>

#include "solver_mip.h"
#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"
#include "type_solver.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"

/*
  #define EXPERIMENTALUNBOUNDED
*/

#define SOLVER_NAME          "mip"
#define SOLVER_DESC          "mip solver for pricing problems"
#define SOLVER_PRIORITY      0

#define DEFAULT_CHECKSOLS    TRUE


/** branching data for branching decisions */
struct GCG_SolverData
{
   SCIP*                 origprob;           /**< original SCIP data structure */

   SCIP_Bool             checksols;          /**< should solutions be checked extensively */
};

/** checks whether the given solution is equal to one of the former solutions in the sols array */
static
SCIP_RETCODE checkSolNew(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   SCIP_SOL**            sols,               /**< array of solutions */
   int                   idx,                /**< index of the solution */
   SCIP_Bool*            isnew               /**< pointer to store whether the solution is new */
   )
{
   SCIP_VAR** probvars;
   int nprobvars;
   SCIP_Real* newvals;

   int s;
   int i;

   assert(scip != NULL);
   assert(pricingprob != NULL);
   assert(sols != NULL);
   assert(sols[idx] != NULL);
   assert(isnew != NULL);

   assert(GCGpricerGetOrigprob(scip) != NULL);

   probvars = SCIPgetVars(pricingprob);
   nprobvars = SCIPgetNVars(pricingprob);

   *isnew = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &newvals, nprobvars) );

   SCIP_CALL( SCIPgetSolVals(pricingprob, sols[idx], nprobvars, probvars, newvals) );

   for( s = 0; s < idx && *isnew == TRUE; s++ )
   {
      assert(sols[s] != NULL);
      /** @todo ensure that the solutions are sorted  */
      /*assert(SCIPisLE(scip, SCIPgetSolOrigObj(pricingprob, sols[s]), SCIPgetSolOrigObj(pricingprob, sols[idx]))
       *|| ABS(SCIPgetSolOrigObj(pricingprob, sols[s])) > 1e+15 * SCIPepsilon(pricingprob));*/
      if( !SCIPisEQ(scip, SCIPgetSolOrigObj(pricingprob, sols[s]), SCIPgetSolOrigObj(pricingprob, sols[idx])) )
         continue;

      if( SCIPsolGetOrigin(sols[s]) != SCIP_SOLORIGIN_ORIGINAL && SCIPsolGetOrigin(sols[idx]) != SCIP_SOLORIGIN_ORIGINAL )
         continue;

      for( i = 0; i < nprobvars; i++ )
         if( !SCIPisEQ(scip, SCIPgetSolVal(pricingprob, sols[s], probvars[i]), newvals[i]) )
            break;

      if( i == nprobvars )
         *isnew = FALSE;
   }

   SCIPfreeBufferArray(scip, &newvals);

   return SCIP_OKAY;
}

static
SCIP_RETCODE solveProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   int                   probnr,             /**< problem number */
   GCG_SOLVERDATA*       solverdata,         /**< solver data structure */
   SCIP_SOL**            sols,               /**< array of solutions */
   SCIP_Bool*            solisray,           /**< array indicating whether solution is a ray */
   int                   maxsols,            /**< size of preallocated array */
   int*                  nsols,              /**< pointer to store number of solutions */
   SCIP_Real*            lowerbound,         /**< pointer to store lower bound */
   SCIP_STATUS*          result              /**< pointer to store pricing problem status */
   )
{
   SCIP_SOL** probsols;
   int nprobsols;
   SCIP_VAR** probvars;
   int nprobvars;
   SCIP_Bool newsol;

   int s;
   int i;

   SCIP_Real memlimit;

   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memlimit) );
   if( !SCIPisInfinity(scip, memlimit) )
   {
      memlimit -= SCIPgetMemUsed(scip)/1048576.0 + GCGgetPricingprobsMemUsed(solverdata->origprob) - SCIPgetMemUsed(pricingprob)/1048576.0;
      if( memlimit < 0 )
         memlimit = 0.0;
      SCIP_CALL( SCIPsetRealParam(pricingprob, "limits/memory", memlimit) );
   }

   /* solve the pricing submip */
   SCIP_CALL( SCIPsolve(pricingprob) );

   /* all SCIP statuses handled so far */
   assert(SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_GAPLIMIT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFEASIBLE
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_UNBOUNDED
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFORUNBD
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_STALLNODELIMIT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_MEMLIMIT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_UNKNOWN);

   if( SCIPgetStatus(pricingprob) == SCIP_STATUS_UNBOUNDED || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFORUNBD )
   {
      /* the pricing problem was declared to be (infeasible or) unbounded, but SCIP did not compute a primal ray;
       * this occurs when presolving detected (infeasibility or) unboundedness; since we need a primal ray to create
       * the corresponding variable, we disable presolving and resolve the problem to get the primal ray out of the LP
       */
      if( !SCIPhasPrimalRay(pricingprob) )
      {
         SCIP_CALL( SCIPfreeTransform(pricingprob) );

         SCIP_CALL( SCIPsetIntParam(pricingprob, "presolving/maxrounds", 0) );
         SCIP_CALL( SCIPtransformProb(pricingprob) );

         /* solve the pricing submip */
         SCIP_CALL( SCIPsolve(pricingprob) );
      }
      assert(SCIPhasPrimalRay(pricingprob)
         || SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT
         || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT );
   }

   SCIPdebugMessage("MIP pricing solver: status = %d\n", SCIPgetStatus(pricingprob));

   /* the pricing problem was declared to be (infeasible or) unbounded and we should have a primal ray at hand,
    * so copy the primal ray into the solution structure and mark it to be a primal ray
    */
   if( SCIPgetStatus(pricingprob) == SCIP_STATUS_UNBOUNDED || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFORUNBD )
   {
      SCIP_VAR** solvars;
      SCIP_Real* solvals;
      int nsolvars;

      assert(SCIPhasPrimalRay(pricingprob));

      probvars  = SCIPgetOrigVars(pricingprob);
      nprobvars = SCIPgetNOrigVars(pricingprob);

      nsolvars = 0;

      SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvars, nprobvars) );
      SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvals, nprobvars) );
      /* store the primal ray values */
      for( i = 0; i < nprobvars; i++ )
      {
         if( SCIPisZero(scip, SCIPgetPrimalRayVal(pricingprob, probvars[i])) )
            continue;

         assert(!SCIPisInfinity(scip, SCIPgetPrimalRayVal(pricingprob, probvars[i])));
         assert(!SCIPisInfinity(scip, -SCIPgetPrimalRayVal(pricingprob, probvars[i])));

         solvars[nsolvars] = probvars[i];
         solvals[nsolvars] = SCIPgetPrimalRayVal(pricingprob, probvars[i]);
         nsolvars++;

         SCIPdebugMessage("%s: %g (obj = %g)\n", SCIPvarGetName(probvars[i]), SCIPgetPrimalRayVal(pricingprob, probvars[i]), SCIPvarGetObj(probvars[i]));
      }
      solisray[0] = TRUE;

      SCIP_CALL( SCIPcreateSol(pricingprob, &sols[0], NULL) );
      SCIP_CALL( SCIPsetSolVals(pricingprob, sols[0], nsolvars, solvars, solvals) );
      *nsols = 1;
      *result = SCIP_STATUS_UNBOUNDED;

      SCIPfreeMemoryArray(pricingprob, &solvars);
      SCIPfreeMemoryArray(pricingprob, &solvals);

      SCIPdebugMessage("pricingproblem has an unbounded ray!\n");
   }
   /* the solving process was interrupted, so we have no solutions and set the status pointer accordingly */
   else if( SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT || SCIPgetStatus(pricingprob) == SCIP_STATUS_MEMLIMIT || SCIPgetStatus(pricingprob) == SCIP_STATUS_UNKNOWN )
   {
      *nsols = 0;
      *result = SCIPgetStatus(pricingprob);
   }
   /* the pricing problem was solved to optimality, copy all solutions found into the solution arrays */
   else
   {
      /* get variables of the pricing problem */
      probvars = SCIPgetOrigVars(pricingprob);
      nprobvars = SCIPgetNOrigVars(pricingprob);

      nprobsols = SCIPgetNSols(pricingprob);
      probsols = SCIPgetSols(pricingprob);

      *nsols = 0;

      for( s = 0; s < nprobsols && s < maxsols; s++ )
      {
         SCIP_Bool feasible;

         if( SCIPisInfinity(pricingprob, -SCIPgetSolOrigObj(pricingprob, probsols[s])) )
         {
           SCIPdebugMessage("unbounded solution\n");
         }

         SCIP_CALL( SCIPcheckSolOrig(pricingprob, probsols[s], &feasible, FALSE, FALSE) );

         if( !feasible )
         {
            SCIPwarningMessage(scip, "solution of pricing problem %d not feasible:\n", probnr);
            SCIP_CALL( SCIPcheckSolOrig(pricingprob, probsols[s], &feasible, TRUE, TRUE) );
         }

         /* check whether the solution is equal to one of the previous solutions */
         if( solverdata->checksols )
         {
            SCIP_CALL( checkSolNew(scip, pricingprob, probsols, s, &newsol) );

            if( !newsol )
               continue;
         }

         solisray[*nsols] = FALSE;
         sols[*nsols] = probsols[s];
         *nsols = *nsols + 1;
      }

      if( SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL )
         *lowerbound = SCIPgetDualbound(pricingprob);

      *result = SCIP_STATUS_OPTIMAL;
      SCIPdebugMessage("pricingproblem found %d sols!\n", *nsols);
   }
   return SCIP_OKAY;
}


/*
 * Callback methods for pricing problem solver
 */

/** destructor of pricing solver to free user data (called when SCIP is exiting) */
static
GCG_DECL_SOLVERFREE(solverFreeMip)
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

#define solverInitsolMip NULL
#define solverExitsolMip NULL
#define solverInitMip NULL
#define solverExitMip NULL

/** solving method for pricing solver which solves the pricing problem to optimality */
static
GCG_DECL_SOLVERSOLVE(solverSolveMip)
{  /*lint --e{715}*/
   GCG_SOLVERDATA* solverdata;

#ifdef DEBUG_PRICING_ALL_OUTPUT
   SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", SCIP_VERBLEVEL_HIGH) );
   SCIP_CALL( SCIPwriteParams(pricingprob, "pricing.set", TRUE, TRUE) );
#endif

   solverdata = GCGpricerGetSolverdata(scip, solver);
   assert(solverdata != NULL);

   *lowerbound = -SCIPinfinity(scip);

   SCIP_CALL( solveProblem(scip, pricingprob, probnr, solverdata, sols, solisray, maxsols, nsols, lowerbound, result) );

#ifdef DEBUG_PRICING_ALL_OUTPUT
   SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", 0) );
   SCIP_CALL( SCIPprintStatistics(pricingprob, NULL) );
#endif

   return SCIP_OKAY;
}

/** heuristic solving method of mip solver */
static
GCG_DECL_SOLVERSOLVEHEUR(solverSolveHeurMip)
{  /*lint --e{715}*/
   GCG_SOLVERDATA* solverdata;

#ifdef DEBUG_PRICING_ALL_OUTPUT
   SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", SCIP_VERBLEVEL_HIGH) );
#endif

   solverdata = GCGpricerGetSolverdata(scip, solver);
   assert(solverdata != NULL);

   *lowerbound = -SCIPinfinity(scip);

   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/stallnodes", 100LL) );
   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/nodes", 1000LL) );
   SCIP_CALL( SCIPsetRealParam(pricingprob, "limits/gap", 0.2) );
   /*SCIP_CALL( SCIPsetIntParam(pricingprob, "limits/bestsol", 5) );*/ /* TODO: do we want a solution limit? */

   SCIP_CALL( solveProblem(scip, pricingprob, probnr, solverdata, sols, solisray, maxsols, nsols, lowerbound, result) );

#ifdef DEBUG_PRICING_ALL_OUTPUT
   SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", 0) );
   SCIP_CALL( SCIPprintStatistics(pricingprob, NULL) );
#endif

   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/stallnodes", -1LL) );
   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/nodes", -1LL) );
   SCIP_CALL( SCIPsetRealParam(pricingprob, "limits/gap", 0.0) );
   SCIP_CALL( SCIPsetIntParam(pricingprob, "limits/bestsol", -1) );

   return SCIP_OKAY;
}

/** creates the mip solver for pricing problems and includes it in GCG */
SCIP_RETCODE GCGincludeSolverMip(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   GCG_SOLVERDATA* data;

   SCIP_CALL( SCIPallocMemory( scip, &data) );
   data->origprob = GCGpricerGetOrigprob(scip);

   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY,
         solverSolveMip, solverSolveHeurMip, solverFreeMip, solverInitMip, solverExitMip,
         solverInitsolMip, solverExitsolMip, data) );

   SCIP_CALL( SCIPaddBoolParam(data->origprob, "pricingsolver/mip/checksols",
         "should solutions of the pricing MIPs be checked for duplicity?",
         &data->checksols, TRUE, DEFAULT_CHECKSOLS, NULL, NULL) );


   return SCIP_OKAY;
}
