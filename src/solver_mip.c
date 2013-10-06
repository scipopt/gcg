/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
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
/*
  #define EXPERIMENTALUNBOUNDED
*/

#define SOLVER_NAME          "mip"
#define SOLVER_DESC          "mip solver for pricing problems"
#define SOLVER_PRIORITY      0

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
   SCIP_SOL**            newsol              /**< solution pointer to store new solution */
)
{
   SCIP_VAR** probvars;
   SCIP_VAR** solvars;
   SCIP_Real* solvals;
   int nprobvars;
   int nsolvars;
   int i;

   assert(pricingprob != NULL);
   assert(newsol != NULL);
   assert(SCIPhasPrimalRay(pricingprob));

   probvars = SCIPgetOrigVars(pricingprob);
   nprobvars = SCIPgetNOrigVars(pricingprob);
   nsolvars = 0;
   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvars, nprobvars));
   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvals, nprobvars));

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
      nsolvars++;

      SCIPdebugMessage("%s: %g (obj = %g)\n", SCIPvarGetName(probvars[i]), SCIPgetPrimalRayVal(pricingprob, probvars[i]), SCIPvarGetObj(probvars[i]));
   }

   SCIP_CALL( SCIPfreeSolve(pricingprob, TRUE));
   SCIP_CALL( SCIPtransformProb(pricingprob));
   SCIP_CALL( SCIPcreateSol(pricingprob, newsol, NULL));
   SCIP_CALL( SCIPsetSolVals(pricingprob, *newsol, nsolvars, solvars, solvals));

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

   return SCIP_OKAY;
}
static
SCIP_RETCODE inferFiniteSolution(
   SCIP*                 pricingprob,        /**< pricing problem to solve again */
   SCIP_SOL*             oldsol,             /**< old solution with infinite solution values */
   SCIP_SOL**            newsol,             /**< new solution with finite solution values, or NULL if not possible */
   SCIP_Bool*            solisray,           /**< whether the new solution is a ray */
   SCIP_Bool*            freesol             /**< whether the new solution has to be freed */
   )
{

   SCIP* newscip;
   SCIP_HASHMAP* varmap;
   SCIP_Bool valid;
   SCIP_Real* solvals;
   int i;
   SCIP_Bool isinfinity;

   assert(pricingprob != NULL);
   assert(oldsol != NULL);
   assert(newsol != NULL);
   assert(solisray != NULL);
   assert(freesol != NULL);

   *newsol = NULL;
   *freesol = FALSE;
   SCIPdebugMessage("Trying to infer finite value!\n");

   SCIP_CALL( resolvePricingWithoutPresolving(pricingprob) );

   if( SCIPhasPrimalRay(pricingprob) )
   {
      SCIP_CALL(createSolutionFromRay(pricingprob, newsol) );

      *freesol = TRUE;
      *solisray = TRUE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &solvals, SCIPgetNOrigVars(pricingprob)) );

   *solisray = FALSE;
   oldsol = SCIPgetBestSol(pricingprob);
   SCIP_CALL( SCIPgetSolVals(pricingprob, oldsol, SCIPgetNOrigVars(pricingprob), SCIPgetOrigVars(pricingprob), solvals) );

   isinfinity = FALSE;
   for( i = 0; i < SCIPgetNOrigVars(pricingprob); ++i )
   {
      if( SCIPisInfinity(pricingprob, solvals[i]) )
      {
         isinfinity = TRUE;
         break;
      }
   }

   if( !isinfinity )
   {
      SCIPdebugMessage("No infinite value found after resolving \n");
      SCIPdebug(SCIP_CALL( SCIPprintSol(pricingprob, oldsol, NULL, FALSE) ));
      *newsol = oldsol;
      *freesol = FALSE;
      SCIPfreeMemoryArray(newscip, &solvals);
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Infinite value found\n");
   SCIPdebug(SCIP_CALL( SCIPprintSol(pricingprob, oldsol, NULL, FALSE) ));

   SCIP_CALL( SCIPcreate(&newscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(newscip));
   SCIP_CALL( SCIPcreateProb(newscip, "finite", NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(newscip), SCIPgetNOrigVars(pricingprob)) );
   SCIP_CALL( SCIPsetIntParam(newscip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
   SCIP_CALL( SCIPsetPresolving(newscip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetBoolParam(newscip, "presolving/donotmultaggr", TRUE) );

   /* make sure that we create a copy of the original problem */
   SCIP_CALL( SCIPfreeTransform(pricingprob) );

   /* copy variables to temporary problem */
   for( i = 0; i < SCIPgetNOrigVars(pricingprob); ++i)
   {
      SCIP_VAR* copyvar;
      SCIP_CALL( SCIPgetVarCopy(pricingprob, newscip, SCIPgetOrigVars(pricingprob)[i], &copyvar, varmap, NULL, TRUE, &valid) );
      assert(valid);
      SCIP_CALL( SCIPaddVar(newscip, copyvar) );
   }

   /* copy constraints to temporary problem */
   for( i = 0; i < SCIPgetNOrigConss(pricingprob); ++i)
   {
      SCIP_CONS* copycons;
      SCIP_CALL( SCIPgetConsCopy(pricingprob, newscip, SCIPgetOrigConss(pricingprob)[i], &copycons, SCIPconsGetHdlr(SCIPgetOrigConss(pricingprob)[i]), varmap, NULL, NULL, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, &valid) );
      assert(valid);
      SCIP_CALL( SCIPaddCons(newscip, copycons) );
   }

   /* roll original problem back to transformed stage */
   SCIP_CALL( SCIPtransformProb(pricingprob) );

   /* change objective function to 0 in temporary problem and create absolute value for infinite variables */
   for( i = 0; i < SCIPgetNOrigVars(pricingprob); ++i )
   {
      SCIP_VAR* newvar;
      SCIP_VAR* origvar;
      origvar = SCIPgetOrigVars(pricingprob)[i];
      SCIPdebugMessage("Origvar <%s>\n", SCIPvarGetName(origvar));
      newvar = SCIPhashmapGetImage(varmap, origvar);
      assert(newvar != NULL);

      SCIP_CALL( SCIPchgVarObj(newscip, newvar, 0.0) );
      if( SCIPisInfinity(pricingprob, solvals[i]) )
      {
         SCIP_VAR* plus;
         SCIP_VAR* minus;
         SCIP_CONS* cons;
         char name[SCIP_MAXSTRLEN];
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_p", SCIPvarGetName(origvar));

         assert(valid);
         SCIP_CALL( SCIPcreateVarBasic(newscip, &plus, name, 0.0, SCIPinfinity(newscip), 0.0001, SCIP_VARTYPE_CONTINUOUS));
         SCIP_CALL( SCIPaddVar(newscip, plus) );

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_m", SCIPvarGetName(origvar));
         SCIP_CALL( SCIPcreateVarBasic(newscip, &minus, name, 0.0, SCIPinfinity(newscip), 0.0001, SCIP_VARTYPE_CONTINUOUS));
         SCIP_CALL( SCIPaddVar(newscip, minus) );

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "abs_%s_cons", SCIPvarGetName(origvar));
         SCIP_CALL( SCIPcreateConsBasicLinear(newscip, &cons, name, 0, NULL, NULL, 0.0, 0.0 ) );
         SCIP_CALL( SCIPaddCoefLinear(newscip, cons, newvar, -1.0) );
         SCIP_CALL( SCIPaddCoefLinear(newscip, cons, plus, 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(newscip, cons, minus, -1.0) );
         SCIP_CALL( SCIPaddCons(newscip, cons) );
         assert(SCIPisZero(pricingprob, SCIPvarGetObj(origvar)));
      }
      else
      {
         SCIP_Bool infeasible;
         SCIPdebugMessage("Fixing <%s> to %.3g\n", SCIPvarGetName(newvar), solvals[i]);
         SCIP_CALL( SCIPfixVar(newscip, newvar, solvals[i], &infeasible, &valid) );
         assert(valid);
         assert(!infeasible);
      }
   }

   /* solve problem */
   SCIP_CALL( SCIPsolve(newscip) );

   /* if the new solution is optimal, copy new solution to old pricing problem*/
   if( SCIPgetStatus(newscip) == SCIP_STATUS_OPTIMAL )
   {
      SCIPdebug(SCIP_CALL( SCIPprintSol(newscip, SCIPgetBestSol(newscip), NULL, FALSE) ));
      *freesol = TRUE;
      SCIP_CALL( SCIPcreateSol(pricingprob, newsol, NULL) );
      SCIP_CALL( SCIPreallocMemoryArray(pricingprob, &solvals, SCIPgetNOrigVars(pricingprob)) );
      for( i = 0; i < SCIPgetNOrigVars(pricingprob); ++i )
      {
         SCIP_VAR* newvar;
         SCIP_VAR* origvar;
         origvar = SCIPgetOrigVars(pricingprob)[i];
         newvar = SCIPhashmapGetImage(varmap, origvar);
         solvals[i] = SCIPgetSolVal(newscip, SCIPgetBestSol(newscip), newvar);
      }
      SCIP_CALL( SCIPsetSolVals(pricingprob, *newsol, SCIPgetNOrigVars(pricingprob), SCIPgetOrigVars(pricingprob), solvals) );
      SCIPdebug(SCIP_CALL( SCIPprintSol(pricingprob, *newsol, NULL, TRUE) ));
   }

   SCIPhashmapFree(&varmap);
   SCIPfreeMemoryArray(pricingprob, &solvals);
   assert(oldsol != NULL);
   assert(SCIPgetStatus(newscip) == SCIP_STATUS_OPTIMAL || *newsol != NULL);

   SCIPfree(&newscip);

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
      if( !SCIPisEQ(pricingprob, SCIPgetSolOrigObj(pricingprob, sols[s]), SCIPgetSolOrigObj(pricingprob, sols[idx])) )
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
   if( SCIPgetStatus(pricingprob) == SCIP_STATUS_UNBOUNDED || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFORUNBD )
      return TRUE;
   else
      return SCIPisInfinity(pricingprob, -SCIPgetSolOrigObj(pricingprob, SCIPgetBestSol(pricingprob))) ||
            SCIPisLT(pricingprob, SCIPinfinity(pricingprob), -SCIPgetSolOrigObj(pricingprob, SCIPgetBestSol(pricingprob)));
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
SCIP_RETCODE filterInfiniteSolutions(
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   SCIP_SOL**            sols,               /**< array of solutions */
   int                   *nsols              /**< number of solutions */
   )
{
   int s;
   int i;
   SCIP_VAR** origvars;
   int norigvars;

   origvars = SCIPgetOrigVars(pricingprob);
   norigvars = SCIPgetNOrigVars(pricingprob);
   for( s = 0; s < *nsols; ++s)
   {
      for( i = 0; i < norigvars; ++i )
      {
         if( SCIPisInfinity(pricingprob, SCIPgetSolVal(pricingprob, sols[s], origvars[i])) )
         {
            sols[s] = sols[*nsols-1];
            --(*nsols);
            --s;
         }
      }
   }
   assert(*nsols >= 0);

   return SCIP_OKAY;
}

static
SCIP_RETCODE solveProblem(
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   int                   probnr,             /**< problem number */
   GCG_SOLVERDATA*       solverdata,         /**< solver data structure */
   SCIP_SOL**            sols,               /**< array of solutions */
   SCIP_Bool*            solisray,           /**< array indicating whether solution is a ray */
   int                   maxsols,            /**< size of preallocated array */
   int*                  nsols,              /**< pointer to store number of solutions */
   SCIP_Real*            lowerbound,         /**< pointer to store lower bound */
   SCIP_STATUS*          status              /**< pointer to store pricing problem status */
   )
{
   SCIP_SOL** probsols;
   int nprobsols;
   SCIP_Bool newsol;

   int s;
   int i;

   SCIP_RETCODE retcode;
   /* solve the pricing submip */
   retcode = SCIPsolve(pricingprob);

   if( retcode != SCIP_OKAY )
   {
      SCIPwarningMessage(pricingprob, "Encountered non recoverable issues solving pricingproblem, ignoring problem\n");
      *status = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }
   SCIPdebugMessage("MIP pricing solver: status = %d\n", SCIPgetStatus(pricingprob));

   *status = SCIPgetStatus(pricingprob);

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

   if( problemIsUnbounded(pricingprob) )
   {
      if( !SCIPhasPrimalRay(pricingprob) )
      {
         SCIP_CALL( resolvePricingWithoutPresolving(pricingprob) );
      }

      SCIP_CALL( createSolutionFromRay(pricingprob, &sols[0]) );

      *nsols = 1;
      solisray[0] = TRUE;
      *status = SCIP_STATUS_UNBOUNDED;
   }
   /* the solving process was interrupted, so we have no solutions and set the status pointer accordingly */
   else if( problemIsInterrupted(pricingprob) )
   {
      *nsols = 0;
      *status = SCIPgetStatus(pricingprob);
   }
   /* the pricing problem has an unbounded solution but finite solution, resolve it and extract a finite solution if necessary */
   else if( problemHasUnboundedSolution(pricingprob) )
   {
      SCIP_SOL* sol = NULL;
      SCIP_Bool isray;
      SCIP_Bool freesol;
      *status = SCIP_STATUS_UNKNOWN;
      SCIP_CALL( inferFiniteSolution(pricingprob, SCIPgetBestSol(pricingprob), &sol, &isray, &freesol) );
      if( sol != NULL )
      {
         *nsols = 1;
         sols[0] = sol;
         solisray[0] = isray;
         SCIPdebugMessage("# changing solution\n");
         *status = isray ? SCIP_STATUS_UNBOUNDED : SCIP_STATUS_OPTIMAL;
      }
   }
   /* the pricing problem was solved to optimality, copy all solutions found into the solution arrays */
   else if( problemIsFeasible(pricingprob) )
   {
      nprobsols = SCIPgetNSols(pricingprob);
      probsols = SCIPgetSols(pricingprob);

      *nsols = 0;

      for( s = 0; s < nprobsols && s < maxsols; s++ )
      {
         SCIP_Bool feasible;
         assert(probsols[s] != NULL);
         SCIP_CALL( SCIPcheckSolOrig(pricingprob, probsols[s], &feasible, FALSE, FALSE) );

         if( !feasible )
         {
            SCIPwarningMessage(pricingprob, "solution of pricing problem %d not feasible:\n", probnr);
            SCIP_CALL( SCIPcheckSolOrig(pricingprob, probsols[s], &feasible, TRUE, TRUE) );
            if( *status != SCIP_STATUS_OPTIMAL )
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

         solisray[*nsols] = FALSE;
         SCIP_CALL( SCIPcreateSolCopy(pricingprob, &sols[*nsols], probsols[s]) );
         *nsols = *nsols + 1;
      }

      if( SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL || SCIPgetStatus(pricingprob) == SCIP_STATUS_GAPLIMIT )
         *lowerbound = SCIPgetDualbound(pricingprob);

      *status = SCIP_STATUS_OPTIMAL;
      SCIPdebugMessage("pricingproblem found %d sols!\n", *nsols);
   }

   if( *nsols > 0)
   {
      SCIP_CALL( filterInfiniteSolutions(pricingprob, sols, nsols) );
      if( nsols == 0)
      {
         *status = SCIP_STATUS_UNKNOWN;
      }
   }

   assert(*nsols >= 0);

   for( s = 0; s < *nsols; ++s)
   {
      for( i = 0; i < SCIPgetNOrigVars(pricingprob); ++i )
      {
         assert( !SCIPisInfinity(pricingprob, SCIPgetSolVal(pricingprob, sols[s], SCIPgetOrigVars(pricingprob)[i])) );
      }
      for( i = 0; i < SCIPgetNVars(pricingprob); ++i )
      {
         assert( !SCIPisInfinity(pricingprob, SCIPgetSolVal(pricingprob, sols[s], SCIPgetVars(pricingprob)[i])) );
      }
      SCIPdebugMessage("Solution has no infinite values.\n");
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
   SCIP_CALL( solveProblem(pricingprob, probnr, solverdata, sols, solisray, maxsols, nsols, lowerbound, result) );

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

   SCIP_CALL( solveProblem(pricingprob, probnr, solverdata, sols, solisray, maxsols, nsols, lowerbound, result) );

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

   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY,
         solverSolveMip, solverSolveHeurMip, solverFreeMip, solverInitMip, solverExitMip,
         solverInitsolMip, solverExitsolMip, data) );

   SCIP_CALL( SCIPaddBoolParam(GCGpricerGetOrigprob(scip), "pricingsolver/mip/checksols",
         "should solutions of the pricing MIPs be checked for duplicity?",
         &data->checksols, TRUE, DEFAULT_CHECKSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(GCGpricerGetOrigprob(scip), "pricingsolver/mip/settingsfile",
         "settings file for pricing problems",
         &data->settingsfile, TRUE, DEFAULT_SETTINGSFILE, NULL, NULL) );


   return SCIP_OKAY;
}
