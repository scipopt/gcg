/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file   solver_gcg.c
 * @brief  gcg solver for pricing problem
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

// #define SCIP_DEBUG
// #define DEBUG_PRICING_ALL_OUTPUT
#define SUBGCG_DEBUG_ITER -1

#include "scip/scip.h"
#include "gcg.h"
#include "objpricer_gcg.h"
#include "solver_gcg.h"
#include "gcgplugins.h"
#include "struct_solver.h"
#include "solver_mip.h"
#include "cons_decomp.h"
#include "cons_decomp.hpp"
#include "pub_solver.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "cons_origbranch.h"
#include "class_partialdecomp.h"
#include "scip/scip_timing.h"
#include "struct_decomp.h"
#include "class_detprobdata.h"

#define SOLVER_NAME          "gcg"
#define SOLVER_DESC          "gcg solver for pricing problems"
#define SOLVER_PRIORITY      100

#define SOLVER_ENABLED      TRUE  /**< indicates whether the solver should be enabled */
#define SOLVER_HEU_ENABLED  FALSE  /**< indicates whether the heuristic solving method of the solver should be enabled */

#define DEFAULT_MAX_RECURSION_DEPTH  0
#define DEFAULT_CHECKSOLS            TRUE    /**< should solutions be checked extensively */
#define DEFAULT_STARTNODELIMIT       1000LL  /**< start node limit for heuristic pricing */
#define DEFAULT_STARTSTALLNODELIMIT  100LL   /**< start stalling node limit for heuristic pricing */
#define DEFAULT_STARTGAPLIMIT        0.2     /**< start gap limit for heuristic pricing */
#define DEFAULT_STARTSOLLIMIT        10      /**< start solution limit for heuristic pricing */
#define DEFAULT_NODELIMITFAC         1.25    /**< factor by which to increase node limit for heuristic pricing */
#define DEFAULT_STALLNODELIMITFAC    1.25    /**< factor by which to increase stalling node limit for heuristic pricing */
#define DEFAULT_GAPLIMITFAC          0.8     /**< factor by which to decrease gap limit for heuristic pricing */
#define DEFAULT_SOLLIMITFAC          1.5     /**< factor by which to increase solution limit for heuristic pricing */
#define DEFAULT_SETTINGSFILE         "-"     /**< settings file to be applied in pricing problems */

/*
 * Data structures
 */

/** branching data for branching decisions */
struct GCG_SolverData
{
   SCIP*                   origprob;            /**< original problem SCIP instance */
   SCIP*                   masterprob;          /**< master problem SCIP instance */
   SCIP**                  pricingprobs;        /**< array storing the SCIP instances for all pricing problems */
   SCIP_HASHMAP**          varmaps;
   int                     depth;
   int                     maxdepth;
   int                     npricingprobs;
   int*                    nbasicpricingconss;  /**< array storing the basic number of constraints of the pricing problems */

   SCIP_Longint            startnodelimit;      /**< start node limit for heuristic pricing */
   SCIP_Longint            startstallnodelimit; /**< start stalling node limit for heuristic pricing */
   SCIP_Real               startgaplimit;       /**< start gap limit for heuristic pricing */
   int                     startsollimit;       /**< start solution limit for heuristic pricing */
   SCIP_Real               nodelimitfac;        /**< factor by which to increase node limit for heuristic pricing */
   SCIP_Real               stallnodelimitfac;   /**< factor by which to increase stalling node limit for heuristic pricing */
   SCIP_Real               gaplimitfac;         /**< factor by which to decrease gap limit for heuristic pricing */
   SCIP_Real               sollimitfac;         /**< factor by which to increase solution limit for heuristic pricing */
   char*                   settingsfile;        /**< settings file to be applied in pricing problems */

   SCIP_Longint*           curnodelimit;        /**< current node limit per pricing problem */
   SCIP_Longint*           curstallnodelimit;   /**< current stalling node limit per pricing problem */
   SCIP_Real*              curgaplimit;         /**< current gap limit per pricing problem */
   int*                    cursollimit;         /**< current solution limit per pricing problem */


   SCIP_Clock* inittime;
   SCIP_Clock* updatetime;
   SCIP_Clock* solvingtime;
   SCIP_Clock* postprocessingtime;
   int64_t count;

   SCIP_Bool               checksols;           /**< should solutions be checked extensively */
};


/*
 * Local methods
 */

static
void solverGcgSetDepth(GCG_SOLVER* solver, int depth)
{
   GCG_SOLVERDATA* solverdata;
   solverdata = GCGsolverGetData(solver);
   solverdata->depth = depth;
}

static
SCIP_RETCODE adjustSettings(SCIP* pricingprob, SCIP* subgcg)
{
   SCIP_Real infinity;
   SCIP_Real epsilon;
   SCIP_Real sumepsilon;
   SCIP_Real feastol;
   SCIP_Real lpfeastolfactor;
   SCIP_Real dualfeastol;

   SCIP_CALL( SCIPsetBoolParam(subgcg, "conflict/useprop", FALSE) );
   SCIP_CALL( SCIPsetCharParam(subgcg, "conflict/useinflp", 'o') );
   SCIP_CALL( SCIPsetCharParam(subgcg, "conflict/useboundlp", 'o') );
   SCIP_CALL( SCIPsetBoolParam(subgcg, "conflict/usesb", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subgcg, "conflict/usepseudo", FALSE) );

   //SCIP_CALL( SCIPsetBoolParam(subgcg, "misc/usevartable", FALSE) );
   //SCIP_CALL( SCIPsetBoolParam(subgcg, "misc/useconstable", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subgcg, "misc/usesmalltables", TRUE) );

   SCIP_CALL( SCIPsetBoolParam(subgcg, "constraints/linear/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subgcg, "constraints/setppc/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subgcg, "constraints/logicor/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subgcg, "constraints/linear/presolusehashing", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subgcg, "constraints/setppc/presolusehashing", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subgcg, "constraints/logicor/presolusehashing", FALSE) );

   SCIP_CALL( SCIPsetIntParam(subgcg, "propagating/dualfix/maxprerounds", 0) );
   SCIP_CALL( SCIPfixParam(subgcg, "propagating/dualfix/maxprerounds") );

   SCIP_CALL( SCIPsetIntParam(subgcg, "limits/maxorigsol", 0) );
   SCIP_CALL( SCIPfixParam(subgcg, "limits/maxorigsol") );

   SCIP_CALL( SCIPsetBoolParam(subgcg, "presolving/donotmultaggr", TRUE) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subgcg, "misc/catchctrlc", FALSE) );

   SCIP_CALL( SCIPsetBoolParam(subgcg, "misc/calcintegral", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subgcg, "misc/finitesolutionstore", TRUE) );

   // TODO: disable dualreds?
   //SCIP_CALL( SCIPsetBoolParam(subgcg, "misc/allowstrongdualreds", FALSE) );
   //SCIP_CALL( SCIPsetBoolParam(subgcg, "misc/allowweakdualreds", FALSE) );

   SCIP_CALL( SCIPgetRealParam(pricingprob, "numerics/infinity", &infinity) );
   SCIP_CALL( SCIPgetRealParam(pricingprob, "numerics/epsilon", &epsilon) );
   SCIP_CALL( SCIPgetRealParam(pricingprob, "numerics/sumepsilon", &sumepsilon) );
   SCIP_CALL( SCIPgetRealParam(pricingprob, "numerics/feastol", &feastol) );
   SCIP_CALL( SCIPgetRealParam(pricingprob, "numerics/lpfeastolfactor", &lpfeastolfactor) );
   SCIP_CALL( SCIPgetRealParam(pricingprob, "numerics/dualfeastol", &dualfeastol) );

   SCIP_CALL( SCIPsetRealParam(subgcg, "numerics/infinity", infinity) );
   SCIP_CALL( SCIPsetRealParam(subgcg, "numerics/epsilon", epsilon) );
   SCIP_CALL( SCIPsetRealParam(subgcg, "numerics/sumepsilon", sumepsilon) );
   SCIP_CALL( SCIPsetRealParam(subgcg, "numerics/feastol", feastol) );
   SCIP_CALL( SCIPsetRealParam(subgcg, "numerics/lpfeastolfactor", lpfeastolfactor) );
   SCIP_CALL( SCIPsetRealParam(subgcg, "numerics/dualfeastol", dualfeastol) );

   // we disable presolving by default (can be enabled using a settings file)
   SCIP_CALL( SCIPsetIntParam(subgcg, "presolving/maxrounds", 0) );

   return SCIP_OKAY;
}

static
SCIP_Bool buildProblem(
   GCG_SOLVERDATA*       solverdata,         /**< solver data structure */
   SCIP*                 pricingprob,        /**< pricing problem */
   int                   probnr              /**< problem number */
)
{
   SCIP* subgcg = NULL;
   SCIP_Bool valid;
   ObjPricerGcg* pricer;
   int nsolvers;
   int i;
   char name[SCIP_MAXSTRLEN];
   GCG_SOLVER* childsolver;
   SCIP_HASHMAP* varmap;
   SCIP_RESULT decompresult = SCIP_DIDNOTRUN;
#ifdef SCIP_DEBUG
   char fname[SCIP_MAXSTRLEN];

   SCIPsnprintf(fname, SCIP_MAXSTRLEN, "%s_subgcg.lp", SCIPgetProbName(pricingprob));
   SCIPwriteOrigProblem(pricingprob, fname, NULL, FALSE);
#endif

   SCIPstartClock(solverdata->origprob, solverdata->inittime);
   SCIP_CALL( SCIPcreate(&subgcg) );
   solverdata->pricingprobs[probnr] = subgcg;

#ifndef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(subgcg, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
#if SCIP_VERSION > 210
   SCIP_CALL( SCIPsetBoolParam(subgcg, "misc/printreason", FALSE) );
#endif
#endif

   SCIP_CALL( SCIPincludeGcgPlugins(subgcg, (SCIP_Bool)(solverdata->depth + 1 < solverdata->maxdepth)) );
   (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_subgcg", SCIPgetProbName(pricingprob) );
   // TODO: check memory leak
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subgcg), SCIPgetNOrigVars(pricingprob)) );
   solverdata->varmaps[probnr] = varmap;
   solverdata->nbasicpricingconss[probnr] = SCIPgetNOrigConss(pricingprob);

   SCIP_CALL( adjustSettings(pricingprob, subgcg) );

   if( strcmp(solverdata->settingsfile, "-") != 0 )
   {
      SCIP_CALL( SCIPreadParams(subgcg, solverdata->settingsfile) );
   }

   // TODO: setting needed parameters (also user specified ones)

   SCIPcreateProb(subgcg, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
   SCIP_CALL( SCIPcopyOrigVars(pricingprob, subgcg, varmap, NULL, NULL, NULL, 0) );
   SCIP_CALL( SCIPcopyOrigConss(pricingprob, subgcg, varmap, NULL, TRUE, &valid) );
   assert(valid);

   pricer = dynamic_cast<ObjPricerGcg*>(SCIPfindObjPricer(GCGgetMasterprob(subgcg), "gcg"));
   assert(pricer != NULL);

   nsolvers = pricer->getNumSolvers();

   for (i = 0; i < nsolvers; ++i)
   {
      childsolver = pricer->getSolvers()[i];
      if (strcmp(childsolver->name, SOLVER_NAME) == 0)
      {
          solverGcgSetDepth(childsolver, solverdata->depth + 1);
          break;
      }
   }

   SCIPstopClock(solverdata->origprob, solverdata->inittime);

   SCIPdebugMessage("SUBGCG Problem %i built, stage: %i\n", probnr, SCIPgetStage(subgcg));

   SCIPdebugMessage("SUBGCG Detecting structure of problem %i\n", probnr);
   SCIPstartClock(solverdata->origprob, solverdata->inittime);
   GCG_DECOMP* decomp = GCGgetStructDecomp(solverdata->origprob);
   assert(decomp);
   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(solverdata->origprob, decomp->partialdecid);
   if( partialdec && partialdec->isNested() )
   {
      gcg::BLOCK_STRUCTURE* blockstructure = partialdec->getBlockStructure(probnr);
      if( blockstructure )
      {
         gcg::DETPROBDATA* detprobdata = GCGconshdlrDecompGetDetprobdataOrig(subgcg);
         assert(subgcg == detprobdata->getScip());
         gcg::PARTIALDECOMP* newpartialdec = blockstructure->createPartialdec(partialdec, detprobdata, probnr);
      }
   }
   else
   {
      SCIP_CALL(GCGdetectStructure(subgcg, &decompresult));
      if (decompresult != SCIP_SUCCESS)
      {
         SCIPwarningMessage(pricingprob, "No decomposition found!\n");
         SCIP_CALL( SCIPfree(&subgcg) );
         solverdata->pricingprobs[probnr] = NULL;
         return SCIP_OKAY;
      }
   }

   SCIPstopClock(solverdata->origprob, solverdata->inittime);

   SCIPdebugMessage("SUBGCG Problem %i structure detected, stage: %i\n", probnr, SCIPgetStage(subgcg));

   return TRUE;
}

/** updates bounds and objective coefficients of variables in the given pricing problem */
static
SCIP_RETCODE updateVars(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_SOLVERDATA*       solverdata,         /**< solver data structure */
   SCIP*                 pricingprob,        /**< pricing problem */
   int                   probnr,             /**< problem number */
   SCIP_Bool             varobjschanged,     /**< have the objective coefficients changed? */
   SCIP_Bool             varbndschanged      /**< have the lower and upper bounds changed? */
)
{
   SCIP_VAR** vars;
   int nvars;
   int npricingvars;
   SCIP* subgcg;
   SCIP_HASHMAP* varmap;
   SCIP_RETCODE retval;
   int i;

   subgcg = solverdata->pricingprobs[probnr];
   varmap = solverdata->varmaps[probnr];
   vars = SCIPgetOrigVars(pricingprob);
   nvars = SCIPgetNOrigVars(pricingprob);
   npricingvars = SCIPgetNOrigVars(subgcg);

   assert(npricingvars == nvars);

   retval = SCIP_OKAY;

   /* get new bounds and objective coefficients of variables */
   for (i = 0; i < nvars; i++)
   {
      SCIP_VAR* origvar;
      SCIP_VAR* var;
      SCIP_VAR* suborigvar;
      SCIP_VAR* subvar;

      origvar = vars[i];
      suborigvar = static_cast<SCIP_VAR*>(SCIPhashmapGetImage(varmap, origvar));

      if( SCIPgetStage(pricingprob) >= SCIP_STAGE_TRANSFORMED && !SCIPvarIsTransformed(origvar) )
         var = SCIPvarGetTransVar(origvar);
      else
         var = origvar;

      if( SCIPgetStage(pricingprob) >= SCIP_STAGE_TRANSFORMED && SCIPvarIsTransformed(suborigvar) )
         subvar = SCIPvarGetTransVar(suborigvar);
      else
         subvar = suborigvar;

      if (varbndschanged)
      {
         SCIPchgVarLb(subgcg, subvar, SCIPvarGetLbGlobal(var));
         SCIPchgVarUb(subgcg, subvar, SCIPvarGetUbGlobal(var));
      }

      if (varobjschanged)
         SCIPchgVarObj(subgcg, suborigvar, SCIPvarGetObj(origvar));
   }

   return retval;
}

/** updates branching constraints in the given pricing problem */
static
SCIP_RETCODE updateBranchingConss(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_SOLVERDATA*       solverdata,         /**< solver data structure */
   SCIP*                 pricingprob,        /**< pricing problem */
   int                   probnr              /**< problem number */
)
{
   SCIP_CONS** conss;
   SCIP_CONS** sconss;
   int nconss;
   int nbasicpricingconss;
   int nsubconss;
   int nnewconss = 0;
   SCIP* subgcg;
   SCIP_HASHMAP* varmap;
   int c;

   subgcg = solverdata->pricingprobs[probnr];
   varmap = solverdata->varmaps[probnr];
   conss = SCIPgetOrigConss(pricingprob);
   sconss = SCIPgetOrigConss(subgcg);
   nconss = SCIPgetNOrigConss(pricingprob);
   nbasicpricingconss = solverdata->nbasicpricingconss[probnr];
   nsubconss = SCIPgetNOrigConss(subgcg);

   for (c = nbasicpricingconss; c < nsubconss; ++c)
   {
      SCIPdelCons(subgcg, sconss[c]);
   }

   nnewconss = nconss - nbasicpricingconss;

   if (nnewconss == 0)
      return SCIP_OKAY;

   for (c = nbasicpricingconss; c < nconss; ++c)
   {
      SCIP_Bool valid;
      SCIP_CONS* newcons;
      SCIP_CONS* cons = conss[c];

      SCIP_CALL( SCIPgetConsCopy(pricingprob, subgcg, cons, &newcons, SCIPconsGetHdlr(cons), varmap, NULL, NULL,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), FALSE, SCIPconsIsModifiable(cons),
         SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), FALSE, TRUE, &valid) );

      if (newcons != NULL && valid)
      {
         SCIP_CALL( SCIPaddCons(subgcg, newcons) );

         // TODO: SCIPconsIsConflict(cons) ?

         SCIP_CALL( SCIPreleaseCons(subgcg, &newcons) );
      }
      else
      {
         SCIPerrorMessage("Could not copy constraint %s (conshdlr: %s)!\n", SCIPconsGetName(cons), SCIPconshdlrGetName(SCIPconsGetHdlr(cons)));
         return SCIP_ERROR;
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE solveProblem(
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   SCIP*                 subgcg,
   int                   probnr,             /**< problem number */
   GCG_SOLVERDATA*       solverdata,         /**< solver data structure */
   SCIP_Real*            lowerbound,         /**< pointer to store lower bound */
   GCG_PRICINGSTATUS*    status              /**< pointer to store pricing problem status */
)
{
   GCG_COL* col;
   SCIP_RETCODE retcode;
   SCIP_Real timelimit;
   SCIP_HASHMAP* varmap = solverdata->varmaps[probnr];

   solverdata->count++;

   // set time limit
   SCIP_CALL( SCIPgetRealParam(pricingprob, "limits/time", &timelimit) );
   SCIP_CALL( SCIPsetRealParam(subgcg, "limits/time", timelimit) );

   if( SCIPgetStage(subgcg) == SCIP_STAGE_PROBLEM )
   {
#ifdef SCIP_DEBUG
      if( SUBGCG_DEBUG_ITER >= 0 && solverdata->count == SUBGCG_DEBUG_ITER )
      {
         SCIPwriteOrigProblem(pricingprob, "pricingprob.lp", NULL, FALSE);
         SCIPwriteOrigProblem(subgcg, "subgcg.lp", NULL, FALSE);
         SCIPwriteOrigProblem(subgcg, "subgcg.dec", NULL, FALSE);
         SCIPwriteParams(subgcg, "params.txt", FALSE, FALSE);
      }
#endif
      SCIPstartClock(solverdata->origprob, solverdata->inittime);
      SCIP_CALL( SCIPpresolve(subgcg) );
      GCGconshdlrDecompTranslateNBestOrigPartialdecs(subgcg, 1, TRUE);
      SCIPstopClock(solverdata->origprob, solverdata->inittime);
#ifdef SCIP_DEBUG
      if( SUBGCG_DEBUG_ITER >= 0 && solverdata->count == SUBGCG_DEBUG_ITER )
      {
         SCIPwriteTransProblem(subgcg, "subgcg_p.lp", NULL, FALSE);
         SCIPwriteTransProblem(subgcg, "subgcg_p.dec", NULL, FALSE);
      }
#endif
   }

   SCIPstartClock(solverdata->origprob, solverdata->solvingtime);
   retcode = SCIPsolve(subgcg);
   SCIPstopClock(solverdata->origprob, solverdata->solvingtime);
   SCIPdebugMessage("Problem %i solved\n", probnr);

   //solverdata->iters += GCGmasterGetPricingSimplexIters(GCGgetMasterprob(subgcg));
   //solverdata->iters += SCIPgetNLPIterations(GCGgetMasterprob(subgcg));

   SCIPstartClock(solverdata->origprob, solverdata->postprocessingtime);

   if (retcode != SCIP_OKAY)
   {
      SCIPwarningMessage(pricingprob, "Pricing problem %d terminated with retcode = %d, ignoring\n", probnr, retcode);
      return SCIP_OKAY;
   }
   SCIPdebugMessage("  -> status = %d\n", SCIPgetStatus(subgcg));
   SCIPdebugMessage("  -> nsols = %d\n", SCIPgetNSols(subgcg));

   *status = getPricingstatus(subgcg);
   SCIPdebugMessage("GCG Solver: Pricingstatus after solve: %u\n", *status);

   switch (*status)
   {
   case GCG_PRICINGSTATUS_INFEASIBLE:
      SCIPdebugMessage("  -> infeasible.\n");
      break;

      /* The pricing problem was declared to be unbounded and we should have a primal ray at hand,
      * so copy the primal ray into the solution structure and mark it to be a primal ray
      */
   case GCG_PRICINGSTATUS_UNBOUNDED:
      if (!SCIPhasPrimalRay(subgcg))
      {
         GCGconshdlrDecompFreeOrigOnExit(subgcg, FALSE);
         SCIP_CALL( SCIPfreeTransform(subgcg));
         GCGconshdlrDecompFreeOrigOnExit(subgcg, TRUE);

         SCIP_CALL( SCIPsetIntParam(subgcg, "presolving/maxrounds", 0) );
         SCIP_CALL( SCIPtransformProb(subgcg) );
         SCIP_CALL( SCIPsolve(subgcg) );
         SCIP_CALL( SCIPsetIntParam(subgcg, "presolving/maxrounds", -1) );
      }

      SCIPdebugMessage("  -> unbounded, creating column from ray\n");
      SCIP_CALL( createColumnFromRay(pricingprob, subgcg, varmap, probnr, &col) );
      SCIP_CALL( GCGpricerAddCol(solverdata->masterprob, col) );
      break;

      /* If the pricing problem is neither infeasible nor unbounded, try to extract feasible columns */
   case GCG_PRICINGSTATUS_UNKNOWN:
   case GCG_PRICINGSTATUS_SOLVERLIMIT:
   case GCG_PRICINGSTATUS_OPTIMAL:
      assert(SCIPgetNSols(subgcg) > 0
         || (SCIPgetStatus(subgcg) != SCIP_STATUS_OPTIMAL
            && SCIPgetStatus(subgcg) != SCIP_STATUS_GAPLIMIT
            && SCIPgetStatus(subgcg) != SCIP_STATUS_SOLLIMIT));

      /* Transform at most maxcols many solutions from the pricing problem into columns */
      SCIP_CALL( getColumnsFromPricingprob(solverdata->masterprob, pricingprob, subgcg, varmap, probnr, solverdata->checksols) );

      *lowerbound = SCIPgetDualbound(subgcg);

      SCIPdebugMessage("  -> lowerbound = %.4g\n", *lowerbound);
      break;

   default:
      SCIPerrorMessage("Pricing problem %d has invalid status: %d\n", probnr, SCIPgetStatus(subgcg));
      break;
   }

   SCIPstopClock(solverdata->origprob, solverdata->postprocessingtime);
   SCIPdebugMessage("Postprocessing of problem %d finished\n", probnr);

   return SCIP_OKAY;
}

/*
 * Callback methods of propagator
 */

/** destructor of pricing solver to free user data (called when SCIP is exiting) */
static
GCG_DECL_SOLVERFREE(solverFreeGcg)
{  /*lint --e{715}*/
   GCG_SOLVERDATA* solverdata;

   assert(scip != NULL);
   assert(solver != NULL);

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   SCIPfreeMemory(scip, &solverdata);

   GCGsolverSetData(solver, NULL);

   return SCIP_OKAY;
}

static
SCIP_RETCODE freeBlockMemory(
   SCIP* scip,
   GCG_SOLVERDATA* solverdata
   )
{
   int i;
   int npricingprobs = solverdata->npricingprobs;

   if (solverdata->pricingprobs == NULL)
      return SCIP_OKAY;

   for (i = 0; i < npricingprobs; ++i)
   {
      if (GCGisPricingprobRelevant(solverdata->origprob, i))
      {
         if (solverdata->pricingprobs[i] != NULL)
         {
            GCGconshdlrDecompFreeDetprobdata(solverdata->pricingprobs[i]);
            SCIPhashmapFree(&(solverdata->varmaps[i]));
            SCIP_CALL( SCIPfree(&(solverdata->pricingprobs[i])) );
            solverdata->pricingprobs[i] = NULL;
         }
         else
            break;
      }
   }

   SCIPfreeBlockMemoryArray(scip, &(solverdata->pricingprobs), npricingprobs);
   solverdata->pricingprobs = NULL;
   SCIPfreeBlockMemoryArray(scip, &(solverdata->varmaps), npricingprobs);
   solverdata->varmaps = NULL;
   SCIPfreeBlockMemoryArray(scip, &(solverdata->nbasicpricingconss), npricingprobs);
   solverdata->nbasicpricingconss = NULL;
   SCIPfreeBlockMemoryArray(scip, &(solverdata->curnodelimit), npricingprobs);
   solverdata->curnodelimit = NULL;
   SCIPfreeBlockMemoryArray(scip, &(solverdata->curgaplimit), npricingprobs);
   solverdata->curgaplimit = NULL;
   SCIPfreeBlockMemoryArray(scip, &(solverdata->cursollimit), npricingprobs);
   solverdata->cursollimit = NULL;
   SCIPfreeBlockMemoryArray(scip, &(solverdata->curstallnodelimit), npricingprobs);
   solverdata->curstallnodelimit = NULL;

   SCIPfreeClock(solverdata->origprob, &solverdata->inittime);
   SCIPfreeClock(solverdata->origprob, &solverdata->updatetime);
   SCIPfreeClock(solverdata->origprob, &solverdata->solvingtime);
   SCIPfreeClock(solverdata->origprob, &solverdata->postprocessingtime);

   return SCIP_OKAY;
}

/** solving process initialization method of pricing solver (called when branch and bound process is about to begin) */
static
GCG_DECL_SOLVERINITSOL(solverInitsolGcg)
{
   GCG_SOLVERDATA* solverdata;
   int npricingprobs;

   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   if (solverdata->depth >= solverdata->maxdepth)
   {
      SCIPdebugMessage("GCG Solver is disabled (depth %i)!\n", solverdata->depth);
      return SCIP_OKAY;
   }

   solverdata->npricingprobs = GCGgetNPricingprobs(solverdata->origprob);
   npricingprobs = solverdata->npricingprobs;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->pricingprobs), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->varmaps), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->nbasicpricingconss), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->curnodelimit), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->curgaplimit), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->cursollimit), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->curstallnodelimit), npricingprobs) );

   SCIPcreateClock(scip, &solverdata->inittime);
   SCIPcreateClock(scip, &solverdata->updatetime);
   SCIPcreateClock(scip, &solverdata->solvingtime);
   SCIPcreateClock(scip, &solverdata->postprocessingtime);

   for (i = 0; i < npricingprobs; ++i)
   {
      solverdata->pricingprobs[i] = NULL;
      solverdata->nbasicpricingconss[i] = 0;
      solverdata->varmaps[i] = NULL;

      if (GCGisPricingprobRelevant(solverdata->origprob, i))
      {
         if (!buildProblem(solverdata, GCGgetPricingprob(solverdata->origprob, i), i))
         {
            SCIP_CALL( freeBlockMemory(scip, solverdata) );
            return SCIP_OKAY;
         }
         else
         {
            solverdata->curnodelimit[i] = solverdata->startnodelimit;
            solverdata->curgaplimit[i] = solverdata->startgaplimit;
            solverdata->cursollimit[i] = solverdata->startsollimit;
            solverdata->curstallnodelimit[i] = solverdata->startstallnodelimit;
         }
      }
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of pricing solver (called before branch and bound process data is freed) */
static
GCG_DECL_SOLVEREXITSOL(solverExitsolGcg)
{
   GCG_SOLVERDATA* solverdata;

   assert(scip != NULL);
   assert(solver != NULL);

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   if( solverdata->depth == 0 && solverdata->inittime )
      SCIPinfoMessage(
         GCGgetOriginalprob(scip),
         NULL,
         "GCG Solver: Init: %.2f, Update: %.2f, Solving: %.2f, Postprocessing: %.2f, Iters: %li\n",
         SCIPgetClockTime(solverdata->origprob, solverdata->inittime),
         SCIPgetClockTime(solverdata->origprob, solverdata->updatetime),
         SCIPgetClockTime(solverdata->origprob, solverdata->solvingtime),
         SCIPgetClockTime(solverdata->origprob, solverdata->postprocessingtime),
         solverdata->count
      );

   SCIP_CALL( freeBlockMemory(scip, solverdata) );
   return SCIP_OKAY;
}


/** initialization method of pricing solver (called after problem was transformed and solver is active) */
#define solverInitGcg NULL

/** deinitialization method of pricing solver (called before transformed problem is freed and solver is active) */
#define solverExitGcg NULL

/** solving method for pricing solver which solves the pricing problem to optimality */
static
GCG_DECL_SOLVERSOLVE(solverSolveGcg)
{  /*lint --e{715}*/
   SCIP* subgcg = NULL;
   GCG_SOLVERDATA* solverdata;

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   if (solverdata->pricingprobs == NULL)
   {
      *status = GCG_PRICINGSTATUS_NOTAPPLICABLE;
      return SCIP_OKAY;
   }

   *lowerbound = -SCIPinfinity(pricingprob);

   SCIPdebugMessage("GCG Solver %li: solve start, probnr: %i, status: %u\n", solverdata->count, probnr, *status);

   subgcg = solverdata->pricingprobs[probnr];

#ifdef DEBUG_PRICING_ALL_OUTPUT
   SCIP_CALL( SCIPsetIntParam(subgcg, "display/verblevel", SCIP_VERBLEVEL_HIGH) );
#endif

   SCIP_CALL( SCIPsetLongintParam(subgcg, "limits/stallnodes", -1LL) );
   SCIP_CALL( SCIPsetLongintParam(subgcg, "limits/nodes", -1LL) );
   SCIP_CALL( SCIPsetRealParam(subgcg, "limits/gap", 0.0) );
   SCIP_CALL( SCIPsetIntParam(subgcg, "limits/solutions", -1) );

   SCIPdebugMessage("Solving pricing problem %d (pointer: %p)\n", probnr, (void*)pricingprob);

   SCIP_CALL( solveProblem(pricingprob, subgcg, probnr, solverdata, lowerbound, status) );

#ifdef DEBUG_PRICING_ALL_OUTPUT
   SCIP_CALL( SCIPsetIntParam(subgcg, "display/verblevel", 0) );
   SCIP_CALL( SCIPprintStatistics(subgcg, NULL) );
#endif

   SCIPdebugMessage("GCG Solver: solve finished, probnr: %i, status: %u\n", probnr, *status);

   return SCIP_OKAY;
}


static
GCG_DECL_SOLVERSOLVEHEUR(solverSolveHeurGcg)
{  /*lint --e{715}*/
   SCIP* subgcg = NULL;
   GCG_SOLVERDATA* solverdata;
   int heurpricingiters;

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   if (solverdata->pricingprobs == NULL)
   {
      *status = GCG_PRICINGSTATUS_NOTAPPLICABLE;
      return SCIP_OKAY;
   }

   *lowerbound = -SCIPinfinity(pricingprob);

   SCIPdebugMessage("GCG Solver %li: solveHeur start, probnr: %i, status: %u\n", solverdata->count, probnr, *status);

   subgcg = solverdata->pricingprobs[probnr];

#ifdef DEBUG_PRICING_ALL_OUTPUT
   SCIP_CALL( SCIPsetIntParam(subgcg, "display/verblevel", SCIP_VERBLEVEL_HIGH) );
#endif


   /* setup heuristic solver parameters */
   if (SCIPgetStage(subgcg) == SCIP_STAGE_PROBLEM)
   {
      solverdata->curnodelimit[probnr] = solverdata->startnodelimit;
      solverdata->curstallnodelimit[probnr] = solverdata->startstallnodelimit;
      solverdata->curgaplimit[probnr] = solverdata->startgaplimit;
      solverdata->cursollimit[probnr] = solverdata->startsollimit;
   }
   else
   {
      switch( SCIPgetStatus(subgcg) )
      {
      case SCIP_STATUS_NODELIMIT:
         if (solverdata->nodelimitfac > 1.0)
         {
            solverdata->curnodelimit[probnr] *= solverdata->nodelimitfac;
            break;
         }
      case SCIP_STATUS_STALLNODELIMIT:
         if (solverdata->stallnodelimitfac > 1.0)
         {
            solverdata->curstallnodelimit[probnr] *= solverdata->stallnodelimitfac;
            break;
         }
      case SCIP_STATUS_GAPLIMIT:
         if (solverdata->gaplimitfac < 1.0)
         {
            solverdata->curgaplimit[probnr] *= solverdata->gaplimitfac;
            break;
         }
      case SCIP_STATUS_SOLLIMIT:
         if (solverdata->sollimitfac > 1.0)
         {
            solverdata->cursollimit[probnr] *= solverdata->sollimitfac;
            break;
         }
      default:
         *status = GCG_PRICINGSTATUS_UNKNOWN;
         SCIPwarningMessage(pricingprob, "GCG solver: cancelled with status %u\n", SCIPgetStatus(subgcg));
         return SCIP_OKAY;
      }
   }
   SCIP_CALL( SCIPsetLongintParam(subgcg, "limits/nodes", solverdata->curnodelimit[probnr]) );
   SCIP_CALL( SCIPsetLongintParam(subgcg, "limits/stallnodes", solverdata->curstallnodelimit[probnr]) );
   SCIP_CALL( SCIPsetRealParam(subgcg, "limits/gap", solverdata->curgaplimit[probnr]) );
   SCIP_CALL( SCIPsetIntParam(subgcg, "limits/solutions", solverdata->cursollimit[probnr]) );

   SCIP_CALL( SCIPgetIntParam(solverdata->origprob, "pricing/masterpricer/heurpricingiters", &heurpricingiters) );
   SCIP_CALL( SCIPsetIntParam(subgcg, "pricing/masterpricer/heurpricingiters", heurpricingiters) );

   /* solve the pricing problem */
   SCIPdebugMessage("Solving pricing problem %d heuristically (pointer: %p)\n", probnr, (void*)pricingprob);

   SCIP_CALL( solveProblem(pricingprob, subgcg, probnr, solverdata, lowerbound, status) );

#ifdef DEBUG_PRICING_ALL_OUTPUT
   SCIP_CALL( SCIPsetIntParam(subgcg, "display/verblevel", 0) );
#endif

   SCIPdebugMessage("GCG Solver: solveHeur finished, probnr: %i, status: %u\n", probnr, *status);

   return SCIP_OKAY;
}

/** update method for pricing solver, used to update solver specific pricing problem data */
static
GCG_DECL_SOLVERUPDATE(solverUpdateGcg)
{
   GCG_SOLVERDATA* solverdata;

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   if (solverdata->pricingprobs == NULL || solverdata->pricingprobs[probnr] == NULL)
      return SCIP_OKAY;

   SCIPdebugMessage("GCG solver -- update data for problem %d: varobjschanged = %u, varbndschanged = %u, consschanged = %u\n",
      probnr, varobjschanged, varbndschanged, consschanged);

   SCIPstartClock(solverdata->origprob, solverdata->updatetime);

   GCGconshdlrDecompFreeOrigOnExit(solverdata->pricingprobs[probnr], FALSE);
   SCIPfreeTransform(solverdata->pricingprobs[probnr]);
   GCGconshdlrDecompFreeOrigOnExit(solverdata->pricingprobs[probnr], TRUE);

   /* update pricing problem information */
   SCIP_CALL( updateVars(solverdata->masterprob, solverdata, pricingprob, probnr, varobjschanged, varbndschanged) );
   if (consschanged)
   {
      SCIP_CALL( updateBranchingConss(solverdata->masterprob, solverdata, pricingprob, probnr) );
   }

   /* reset heuristic pricing limits */
   solverdata->curnodelimit[probnr] = solverdata->startnodelimit;
   solverdata->curgaplimit[probnr] = solverdata->startgaplimit;
   solverdata->cursollimit[probnr] = solverdata->startsollimit;
   solverdata->curstallnodelimit[probnr] = solverdata->startstallnodelimit;

   SCIPstopClock(solverdata->origprob, solverdata->updatetime);
   SCIPdebugMessage("Updated problem %i\n", probnr);

   return SCIP_OKAY;
}

/** creates the most infeasible LP braching rule and includes it in SCIP */
SCIP_RETCODE GCGincludeSolverGcg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   GCG_SOLVERDATA* solverdata = NULL;
   SCIP* origprob = GCGmasterGetOrigprob(scip);

   /* create gcg solver data */
   SCIP_CALL( SCIPallocMemory(scip, &solverdata) );
   solverdata->depth = 0;
   solverdata->npricingprobs = 0;
   solverdata->inittime = NULL;
   solverdata->updatetime = NULL;
   solverdata->solvingtime = NULL;
   solverdata->postprocessingtime = NULL;
   solverdata->count = 0;
   solverdata->origprob = origprob;
   solverdata->masterprob = scip;
   solverdata->pricingprobs = NULL;
   solverdata->varmaps = NULL;
   solverdata->nbasicpricingconss = NULL;
   solverdata->settingsfile = NULL;

   /* include pricing problem solver */
   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY, SOLVER_HEU_ENABLED, SOLVER_ENABLED,
      solverUpdateGcg, solverSolveGcg, solverSolveHeurGcg, solverFreeGcg, solverInitGcg, solverExitGcg,
         solverInitsolGcg, solverExitsolGcg, solverdata) );

   /* add gcg propagator parameters */
   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/recursivedecomp/maxdepth",
         "maximal recursive decomposition depth",
         &solverdata->maxdepth, FALSE, DEFAULT_MAX_RECURSION_DEPTH, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "pricingsolver/gcg/checksols",
      "should solutions of the pricing MIPs be checked for duplicity?",
      &solverdata->checksols, TRUE, DEFAULT_CHECKSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(origprob, "pricingsolver/gcg/startnodelimit",
      "start node limit for heuristic pricing",
      &solverdata->startnodelimit, TRUE, DEFAULT_STARTNODELIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(origprob, "pricingsolver/gcg/startstallnodelimit",
      "start stalling node limit for heuristic pricing",
      &solverdata->startstallnodelimit, TRUE, DEFAULT_STARTSTALLNODELIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/gcg/startgaplimit",
      "start gap limit for heuristic pricing",
      &solverdata->startgaplimit, TRUE, DEFAULT_STARTGAPLIMIT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricingsolver/gcg/startsollimit",
      "start solution limit for heuristic pricing",
      &solverdata->startsollimit, TRUE, DEFAULT_STARTSOLLIMIT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/gcg/nodelimitfac",
      "factor by which to increase node limit for heuristic pricing",
      &solverdata->nodelimitfac, TRUE, DEFAULT_NODELIMITFAC, 1.0, SCIPinfinity(origprob), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/gcg/stallnodelimitfac",
      "factor by which to increase stalling node limit for heuristic pricing",
      &solverdata->stallnodelimitfac, TRUE, DEFAULT_STALLNODELIMITFAC, 1.0, SCIPinfinity(origprob), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/gcg/gaplimitfac",
      "factor by which to decrease gap limit for heuristic pricing",
      &solverdata->gaplimitfac, TRUE, DEFAULT_GAPLIMITFAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/gcg/sollimitfac",
      "factor by which to increase solution limit for heuristic pricing",
      &solverdata->sollimitfac, TRUE, DEFAULT_SOLLIMITFAC, 1.0, SCIPinfinity(origprob), NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(origprob, "pricingsolver/gcg/settingsfile",
      "settings file for pricing problems",
      &solverdata->settingsfile, TRUE, DEFAULT_SETTINGSFILE, NULL, NULL) );

   return SCIP_OKAY;
}
