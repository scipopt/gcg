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

/**@file   solver_mip.c
 * @brief  pricing solver solving the pricing problem as a sub-MIP, using SCIP
 * @author Gerald Gamrath
 * @author Martin Bergner
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define DEBUG_PRICING_ALL_OUTPUT */

#include <assert.h>
#include <string.h>

#include "solver_mip.h"
#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"
#include "gcg.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "scip/scipdefplugins.h"
#include "pub_gcgcol.h"

#define SOLVER_NAME          "mip"
#define SOLVER_DESC          "mip solver for pricing problems"
#define SOLVER_PRIORITY      0

#define SOLVER_ENABLED      TRUE  /**< indicates whether the solver should be enabled */

#define DEFAULT_CHECKSOLS    TRUE
#define DEFAULT_SETTINGSFILE "-"

/** branching data for branching decisions */
struct GCG_SolverData
{
   SCIP_Bool             checksols;          /**< should solutions be checked extensively */
   char*                 settingsfile;       /**< settings file to be applied in pricing problems */
};

/** extracts ray from pricing problem */
static
SCIP_RETCODE createSolutionFromRay(
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   int                   probnr,             /**< problem number */
   GCG_COL**             newcol              /**< column pointer to store new column */
)
{
   SCIP_VAR** probvars;
   SCIP_VAR** solvars;
   SCIP_Real* solvals;
   int nprobvars;
   int nsolvars;
   int i;

   assert(pricingprob != NULL);
   assert(newcol != NULL);
   assert(SCIPhasPrimalRay(pricingprob));

   probvars = SCIPgetOrigVars(pricingprob);
   nprobvars = SCIPgetNOrigVars(pricingprob);
   nsolvars = 0;
   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvars, nprobvars) );
   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvals, nprobvars) );

   /* store the primal ray values */
   for ( i = 0; i < nprobvars; i++ )
   {
      if ( SCIPisZero(pricingprob, SCIPgetPrimalRayVal(pricingprob, probvars[i])) )
         continue;

      assert(!SCIPisInfinity(pricingprob, SCIPgetPrimalRayVal(pricingprob, probvars[i])));
      assert(!SCIPisInfinity(pricingprob, -SCIPgetPrimalRayVal(pricingprob, probvars[i])));

      solvars[nsolvars] = probvars[i];
      solvals[nsolvars] = SCIPgetPrimalRayVal(pricingprob, probvars[i]);
      assert(!SCIPisInfinity(pricingprob, solvals[nsolvars]));
      /* todo: check if/ensure that this value is integral! */
      nsolvars++;

      SCIPdebugMessage("%s: %g (obj = %g)\n", SCIPvarGetName(probvars[i]), SCIPgetPrimalRayVal(pricingprob, probvars[i]), SCIPvarGetObj(probvars[i]));
   }

   SCIP_CALL( SCIPfreeSolve(pricingprob, TRUE) );
   SCIP_CALL( SCIPtransformProb(pricingprob) );

   /* todo: is it okay to write pricingprob into column structure? */
   SCIP_CALL( GCGcreateGcgCol(pricingprob, newcol, probnr, solvars, solvals, nsolvars, TRUE, SCIPinfinity(pricingprob)) );

   SCIPfreeMemoryArray(pricingprob, &solvars);
   SCIPfreeMemoryArray(pricingprob, &solvals);

   SCIPdebugMessage("pricingproblem has an unbounded ray!\n");

   return SCIP_OKAY;
}

/** solves the pricing problem again without presolving */
static
SCIP_RETCODE resolvePricingWithoutPresolving(
   SCIP*                 pricingprob         /**< pricing problem to solve again */
   )
{
   assert(pricingprob != NULL);
   SCIP_CALL( SCIPfreeTransform(pricingprob) );

   SCIP_CALL( SCIPsetIntParam(pricingprob, "presolving/maxrounds", 0) );
   SCIP_CALL( SCIPtransformProb(pricingprob) );
   SCIP_CALL( SCIPsolve(pricingprob) );
   SCIP_CALL( SCIPsetIntParam(pricingprob, "presolving/maxrounds", -1) );

   return SCIP_OKAY;
}

/** checks whether the given solution is equal to one of the former solutions in the sols array */
static
SCIP_RETCODE checkSolNew(
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

   assert(pricingprob != NULL);
   assert(sols != NULL);
   assert(sols[idx] != NULL);
   assert(isnew != NULL);

   probvars = SCIPgetVars(pricingprob);
   nprobvars = SCIPgetNVars(pricingprob);

   *isnew = TRUE;

   SCIP_CALL( SCIPallocBufferArray(pricingprob, &newvals, nprobvars) );

   SCIP_CALL( SCIPgetSolVals(pricingprob, sols[idx], nprobvars, probvars, newvals) );

   for( s = 0; s < idx && *isnew == TRUE; s++ )
   {
      assert(sols[s] != NULL);
      /** @todo ensure that the solutions are sorted  */
      if( (!SCIPisInfinity(pricingprob, -SCIPgetSolOrigObj(pricingprob, sols[s])) && !SCIPisInfinity(pricingprob, -SCIPgetSolOrigObj(pricingprob, sols[idx])) )
        && !SCIPisEQ(pricingprob, SCIPgetSolOrigObj(pricingprob, sols[s]), SCIPgetSolOrigObj(pricingprob, sols[idx])) )
         continue;

      if( (SCIPisInfinity(pricingprob, -SCIPgetSolOrigObj(pricingprob, sols[s])) && !SCIPisInfinity(pricingprob, -SCIPgetSolOrigObj(pricingprob, sols[idx])) )
       ||(!SCIPisInfinity(pricingprob, -SCIPgetSolOrigObj(pricingprob, sols[s])) &&  SCIPisInfinity(pricingprob, -SCIPgetSolOrigObj(pricingprob, sols[idx]))) )
         continue;


      if( SCIPsolGetOrigin(sols[s]) != SCIP_SOLORIGIN_ORIGINAL && SCIPsolGetOrigin(sols[idx]) != SCIP_SOLORIGIN_ORIGINAL )
         continue;

      for( i = 0; i < nprobvars; i++ )
         if( !SCIPisEQ(pricingprob, SCIPgetSolVal(pricingprob, sols[s], probvars[i]), newvals[i]) )
            break;

      if( i == nprobvars )
         *isnew = FALSE;
   }

   SCIPfreeBufferArray(pricingprob, &newvals);

   return SCIP_OKAY;
}

/** returns whether the solution is unbounded or not */
static
SCIP_Bool problemIsUnbounded(
   SCIP*                 pricingprob         /**< pricing problem SCIP data structure */
   )

{
   return SCIPgetStatus(pricingprob) == SCIP_STATUS_UNBOUNDED || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFORUNBD;
}

/** returns whether the solution process was aborted */
static
SCIP_Bool problemIsInterrupted(
   SCIP*                 pricingprob         /**< pricing problem SCIP data structure */
   )
{
   return SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT ||
          SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT ||
          SCIPgetStatus(pricingprob) == SCIP_STATUS_MEMLIMIT ||
          SCIPgetStatus(pricingprob) == SCIP_STATUS_STALLNODELIMIT ||
          SCIPgetStatus(pricingprob) == SCIP_STATUS_UNKNOWN;
}

/** returns whether the solution process finished */
static
SCIP_Bool problemIsFeasible(
   SCIP*                 pricingprob         /**< pricing problem SCIP data structure */
   )
{
   return SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL || SCIPgetStatus(pricingprob) == SCIP_STATUS_GAPLIMIT;
}

/** returns whether the solution has an infinite value for at least one variable */
static
SCIP_Bool problemHasUnboundedSolution(
   SCIP*                 pricingprob         /**< pricing problem SCIP data structure */
   )
{
   int i;
   SCIP_SOL* sol;
   sol = SCIPgetBestSol(pricingprob);
   assert(SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL || SCIPgetStatus(pricingprob) == SCIP_STATUS_GAPLIMIT );
   /* assert(!SCIPisInfinity(pricingprob, SCIPsolGetOrigObj(sol))); */

   if( sol == NULL )
      return FALSE;

   for( i = 0; i < SCIPgetNOrigVars(pricingprob); ++i )
   {
      if( SCIPisInfinity(pricingprob, SCIPgetSolVal(pricingprob, sol, SCIPgetOrigVars(pricingprob)[i])) )
      {
         return TRUE;
      }
   }

   return FALSE;
}

/** removes wrong solutions from returned solutions */
static
SCIP_RETCODE filterInfiniteColumns(
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   GCG_COL**             cols,               /**< array of columns corresponding to solutions */
   int                   *ncols              /**< number of columns */
   )
{
   int s;
   int i;

   assert(pricingprob != NULL);
   assert(cols != NULL);
   assert(ncols != NULL);

   /*lint --e{850} s is modified in the loop */
   for( s = 0; s < *ncols; ++s )
   {
      SCIP_Real* vals;
      int nvals;

      vals = GCGcolGetVals(cols[s]);
      nvals = GCGcolGetNVars(cols[s]);

      for( i = 0; i < nvals; ++i )
      {
         if( SCIPisInfinity(pricingprob, vals[i]) )
         {
            if( s == 0 )
            {
               SCIPwarningMessage(pricingprob, "Removing solution with infinite value.\n");
            }

            GCGfreeGcgCol(&cols[s]);

            cols[s] = cols[*ncols-1];
            --(*ncols);
            --s;
            break;
         }
      }
   }
   assert(*ncols >= 0);

   return SCIP_OKAY;
}

/** solves the given pricing problem as a sub-SCIP */
static
SCIP_RETCODE solveProblem(
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   int                   probnr,             /**< problem number */
   GCG_SOLVERDATA*       solverdata,         /**< solver data structure */
   GCG_COL**             cols,               /**< array of columns corresponding to solutions */
   int                   maxcols,            /**< size of preallocated array */
   int*                  ncols,              /**< pointer to store number of columns */
   SCIP_Real*            lowerbound,         /**< pointer to store lower bound */
   SCIP_STATUS*          status              /**< pointer to store pricing problem status */
   )
{
   SCIP_Bool newsol;

   int s;

   SCIP_RETCODE retcode;

   assert(pricingprob != NULL);
   assert(probnr >= 0);
   assert(solverdata != NULL);
   assert(cols != NULL);
   assert(maxcols > 0);
   assert(ncols != NULL);
   assert(lowerbound != NULL);
   assert(status != NULL);

   /* solve the pricing submip */
   retcode = SCIPsolve(pricingprob);

   if( retcode != SCIP_OKAY )
   {
      SCIPwarningMessage(pricingprob, "Encountered non recoverable issues solving pricingproblem, ignoring problem\n");
      *status = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }
   SCIPdebugMessage("MIP pricing solver: status = %d\n", SCIPgetStatus(pricingprob));

   /* all SCIP statuses handled so far */
   assert(SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_GAPLIMIT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFEASIBLE
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_UNBOUNDED
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFORUNBD
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_STALLNODELIMIT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_MEMLIMIT);
   /* @todo: can UNKNOWN happen, too? */

   /* the pricing problem was declared to be (infeasible or) unbounded and we should have a primal ray at hand,
    * so copy the primal ray into the solution structure and mark it to be a primal ray
    */

   if( SCIPgetStatus(pricingprob) == SCIP_STATUS_INFEASIBLE )
   {
      SCIPdebugMessage("Pricing is infeasible, abort immediately.\n");
      *status = SCIP_STATUS_INFEASIBLE;
      return SCIP_OKAY;
   }

   if( problemIsUnbounded(pricingprob) && !SCIPhasPrimalRay(pricingprob) )
   {
      SCIP_CALL( resolvePricingWithoutPresolving(pricingprob) );
   }

   if( problemIsUnbounded(pricingprob) )
   {
      SCIP_CALL( createSolutionFromRay(pricingprob, probnr, &cols[0]) );

      *ncols = 1;
      *status = SCIP_STATUS_UNBOUNDED;
   }
   /* the solving process was interrupted, so we have no solutions and set the status pointer accordingly */
   else if( problemIsInterrupted(pricingprob) )
   {
      *ncols = 0;
      *status = SCIP_STATUS_UNKNOWN;
   }
   /* the pricing problem has an unbounded solution but finite solution, resolve it and extract a finite solution if necessary */
   else if( problemHasUnboundedSolution(pricingprob) )
   {
      SCIP_SOL* sol = NULL;
      SCIP_Bool success = FALSE;

      SCIPdebugMessage("solution has infinite values, create a copy with finite values\n");

      SCIP_CALL( SCIPcreateFiniteSolCopy(pricingprob, &sol, SCIPgetBestSol(pricingprob), &success) );
      assert(success);
      assert(sol != NULL);

      SCIP_CALL( GCGcreateGcgColFromSol(pricingprob, &cols[0], probnr, sol, FALSE, SCIPinfinity(pricingprob)) );
      *ncols = 1;
      *status = SCIP_STATUS_OPTIMAL;
   }
   /* the pricing problem was solved to optimality, copy all solutions found into the solution arrays */
   else if( problemIsFeasible(pricingprob) )
   {
      SCIP_SOL** probsols = SCIPgetSols(pricingprob);
      int nprobsols = SCIPgetNSols(pricingprob);;

      *ncols = 0;
      *status = SCIP_STATUS_OPTIMAL;

      for( s = 0; s < nprobsols && *ncols < maxcols; s++ )
      {
         SCIP_Bool feasible;
         assert(probsols[s] != NULL);
         SCIP_CALL( SCIPcheckSolOrig(pricingprob, probsols[s], &feasible, FALSE, FALSE) );

         if( !feasible )
         {
            SCIPwarningMessage(pricingprob, "solution of pricing problem %d not feasible:\n", probnr);
            SCIP_CALL( SCIPcheckSolOrig(pricingprob, probsols[s], &feasible, TRUE, TRUE) );

            /* if the best solution is not feasible, we set the status to UNKNOWN */
            if( s == 0 )
            {
               *status = SCIP_STATUS_UNKNOWN;
            }
         }

         /* check whether the solution is equal to one of the previous solutions */
         if( solverdata->checksols )
         {
            SCIP_CALL( checkSolNew(pricingprob, probsols, s, &newsol) );

            if( !newsol )
               continue;
         }

         SCIP_CALL( GCGcreateGcgColFromSol(pricingprob, &(cols[*ncols]), probnr, probsols[s], FALSE, SCIPinfinity(pricingprob)) );
         ++(*ncols);
      }

      if( (SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL) || (SCIPgetStatus(pricingprob) == SCIP_STATUS_GAPLIMIT) )
         *lowerbound = SCIPgetDualbound(pricingprob);

      SCIPdebugMessage("pricingproblem found %d sols, lowerbound = %.4g!\n", *ncols, *lowerbound);
   }
   else
   {
      SCIPerrorMessage("invalid status of pricing problem %d: %d", probnr, SCIPgetStatus(pricingprob));
      *status = SCIP_STATUS_UNKNOWN;
   }

   if( *ncols > 0 )
   {
      SCIP_CALL( filterInfiniteColumns(pricingprob, cols, ncols) );

      if( *ncols == 0 )
      {
         *status = SCIP_STATUS_UNKNOWN;
      }
   }

   assert(*ncols >= 0);

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

   solverdata = GCGsolverGetSolverdata(solver);
   assert(solverdata != NULL);

   SCIPfreeMemory(scip, &solverdata);

   GCGsolverSetSolverdata(solver, NULL);

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

   solverdata = GCGsolverGetSolverdata(solver);
   assert(solverdata != NULL);

   if( strcmp(solverdata->settingsfile, "-") != 0 )
   {
      SCIP_CALL( SCIPreadParams(pricingprob, solverdata->settingsfile) );
   }

#ifdef DEBUG_PRICING_ALL_OUTPUT
   SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", SCIP_VERBLEVEL_HIGH) );
   SCIP_CALL( SCIPwriteParams(pricingprob, "pricing.set", TRUE, TRUE) );
#endif

   *lowerbound = -SCIPinfinity(pricingprob);
   SCIPdebugMessage("solving pricing %d (pointer: %p)\n", probnr, (void*)pricingprob);
   SCIP_CALL( solveProblem(pricingprob, probnr, solverdata, cols, maxcols, ncols, lowerbound, result) );

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

   solverdata = GCGsolverGetSolverdata(solver);
   assert(solverdata != NULL);

   *lowerbound = -SCIPinfinity(pricingprob);

   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/stallnodes", 100LL) );
   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/nodes", 1000LL) );
   SCIP_CALL( SCIPsetRealParam(pricingprob, "limits/gap", 0.2) );
   /*SCIP_CALL( SCIPsetIntParam(pricingprob, "limits/bestsol", 5) );*/ /* TODO: do we want a solution limit? */

   SCIP_CALL( solveProblem(pricingprob, probnr, solverdata, cols, maxcols, ncols, lowerbound, result) );

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

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->settingsfile = NULL;

   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY, SOLVER_ENABLED,
         solverSolveMip, solverSolveHeurMip, solverFreeMip, solverInitMip, solverExitMip,
         solverInitsolMip, solverExitsolMip, data) );

   SCIP_CALL( SCIPaddBoolParam(GCGmasterGetOrigprob(scip), "pricingsolver/mip/checksols",
         "should solutions of the pricing MIPs be checked for duplicity?",
         &data->checksols, TRUE, DEFAULT_CHECKSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(GCGmasterGetOrigprob(scip), "pricingsolver/mip/settingsfile",
         "settings file for pricing problems",
         &data->settingsfile, TRUE, DEFAULT_SETTINGSFILE, NULL, NULL) );


   return SCIP_OKAY;
}
