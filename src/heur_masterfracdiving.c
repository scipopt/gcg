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

/**@file   heur_masterfracdiving.c
 * @brief  master LP diving heuristic that chooses fixings w.r.t. the fractionalities
 * @author Tobias Achterberg
 * @author Christian Puchert
 */
//#define SCIP_DEBUG
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_masterfracdiving.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"


#define HEUR_NAME             "masterfracdiving"
#define HEUR_DESC             "master LP diving heuristic that chooses fixings w.r.t. the fractionalities"
#define HEUR_DISPCHAR         'f'
#define HEUR_PRIORITY         -1003000
#define HEUR_FREQ             -1
//#define HEUR_FREQOFS          3
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.05 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXPRICEQUOT       0.10 /**< maximal fraction of pricing rounds compared to node pricing rounds */
#define DEFAULT_MAXPRICEOFS          10 /**< additional number of allowed pricing rounds (-1: no limit) */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */

#define MINLPITER                 10000 /**< minimal number of LP iterations allowed in each LP solving call */


/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             minreldepth;        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth;        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   int                   maxlpiterofs;       /**< additional number of allowed LP iterations */
   SCIP_Real             maxpricequot;       /**< maximal fraction of pricing rounds compared to node pricing rounds */
   int                   maxpriceofs;        /**< additional number of allowed pricing rounds (-1: no limit) */
   SCIP_Real             maxdiveubquot;      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveavgquot;     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Bool             backtrack;          /**< use one level of backtracking if infeasibility is encountered? */
   SCIP_Longint          nlpiterations;      /**< LP iterations used in this heuristic */
   int                   npricerounds;       /**< pricing rounds used in this heuristic */
   int                   nsuccess;           /**< number of runs that produced at least one feasible solution */
};


/*
 * local methods
 */


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#define heurCopyMasterfracdiving NULL

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeMasterfracdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitMasterfracdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* initialize data */
   heurdata->nlpiterations = 0;
   heurdata->nsuccess = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitMasterfracdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecMasterfracdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP* origprob;
   SCIP_HEURDATA* heurdata;
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_VAR* var;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   SCIP_Real searchubbound;
   SCIP_Real searchavgbound;
   SCIP_Real searchbound;
   SCIP_Real objval;
   SCIP_Real oldobjval;
   SCIP_Real obj;
   SCIP_Real objgain;
   SCIP_Real bestobjgain;
   SCIP_Real frac;
   SCIP_Real bestfrac;
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   SCIP_Bool bestcandroundup;
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   SCIP_Bool roundup;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;
   SCIP_Bool backtracked;
   SCIP_Bool origfeas;
   SCIP_Longint ncalls;
   SCIP_Longint nsolsfound;
   SCIP_Longint nlpiterations;
   SCIP_Longint maxnlpiterations;
   int maxpricerounds;
   int npricerounds;
   int totalpricerounds;
   int nlpcands;
   int startnlpcands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int bestcand;
   int c;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   *result = SCIP_DELAYED;

   SCIPdebugMessage("called Masterfracdiving heuristic\n");

   /* only call heuristic, if an optimal LP solution is at hand */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("not executing Masterfracdiving heuristic: master LP not solved to optimality\n");
      return SCIP_OKAY;
   }

   /* only call heuristic, if the LP solution is basic (which allows fast resolve in diving) */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only try to dive, if we are in the correct part of the tree, given by minreldepth and maxreldepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   maxdepth = MAX(maxdepth, 30);
   if( depth < heurdata->minreldepth*maxdepth || depth > heurdata->maxreldepth*maxdepth )
      return SCIP_OKAY;

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNNodeLPIterations(scip);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = 10*SCIPheurGetNBestSolsFound(heur) + heurdata->nsuccess;
   maxnlpiterations = (SCIP_Longint)((1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * nlpiterations);
   maxnlpiterations += heurdata->maxlpiterofs;

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;

   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + MINLPITER);

   /** @todo limit number of pricing rounds, play with parameters */
   if( heurdata->maxpriceofs > -1 )
   {
      npricerounds = SCIPgetNPriceRounds(scip);
      SCIPdebugMessage("masterfracdiving - pricing rounds at this node: %d\n", npricerounds);
      maxpricerounds = (int)((1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxpricequot * npricerounds);
      maxpricerounds += heurdata->maxpriceofs;
   }
   else
      maxpricerounds = -1;

   SCIPdebugMessage("Maximum number of LP iters and price rounds: %"SCIP_LONGINT_FORMAT", %d\n", maxnlpiterations, maxpricerounds);

   /* get fractional variables that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );

   /* don't try to dive, if there are no fractional variables */
   if( nlpcands == 0 )
      return SCIP_OKAY;

   /* calculate the objective search bound */
   if( SCIPgetNSolsFound(scip) == 0 )
   {
      if( heurdata->maxdiveubquotnosol > 0.0 )
         searchubbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveubquotnosol * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if( heurdata->maxdiveavgquotnosol > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveavgquotnosol * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   else
   {
      if( heurdata->maxdiveubquot > 0.0 )
         searchubbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveubquot * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if( heurdata->maxdiveavgquot > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveavgquot * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   searchbound = MIN(searchubbound, searchavgbound);
   if( SCIPisObjIntegral(scip) )
      searchbound = SCIPceil(scip, searchbound);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   maxdivedepth = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;


   *result = SCIP_DIDNOTFIND;

   /* start diving */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   SCIPenableVarHistory(scip);

   /* get LP objective value*/
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   objval = SCIPgetLPObjval(scip);

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") executing fracdiving heuristic: depth=%d, %d fractionals, dualbound=%g, searchbound=%g\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), SCIPretransformObj(scip, searchbound));

   /* dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   lperror = FALSE;
   cutoff = FALSE;
   origfeas = FALSE;
   divedepth = 0;
   npricerounds = 0;
   totalpricerounds = 0;
   bestcandmayrounddown = FALSE;
   bestcandmayroundup = FALSE;
   startnlpcands = nlpcands;
   while( !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && (divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations && objval < searchbound))
      && !SCIPisStopped(scip) )
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      divedepth++;

      /* choose variable fixing:
       * - prefer variables that may not be rounded without destroying LP feasibility:
       *   - of these variables, round least fractional variable in corresponding direction
       * - if all remaining fractional variables may be rounded without destroying LP feasibility:
       *   - round variable with least increasing objective value
       */
      bestcand = -1;
      bestobjgain = SCIPinfinity(scip);
      bestfrac = SCIP_INVALID;
      bestcandmayrounddown = TRUE;
      bestcandmayroundup = TRUE;
      bestcandroundup = FALSE;
      for( c = 0; c < nlpcands; ++c )
      {
         var = lpcands[c];
         mayrounddown = SCIPvarMayRoundDown(var);
         mayroundup = SCIPvarMayRoundUp(var);
         frac = lpcandsfrac[c];
         obj = SCIPvarGetObj(var);
         if( mayrounddown || mayroundup )
         {
            /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
            if( bestcandmayrounddown || bestcandmayroundup )
            {
               /* choose rounding direction:
                * - if variable may be rounded in both directions, round corresponding to the fractionality
                * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
                *   the current fractional solution
                */
               if( mayrounddown && mayroundup )
                  roundup = (frac > 0.5);
               else
                  roundup = mayrounddown;

               if( roundup )
               {
                  frac = 1.0 - frac;
                  objgain = frac*obj;
               }
               else
                  objgain = -frac*obj;

               /* penalize too small fractions */
               if( frac < 0.01 )
                  objgain *= 1000.0;

               /* prefer decisions on binary variables */
               if( !SCIPvarIsBinary(var) )
                  objgain *= 1000.0;

               /* check, if candidate is new best candidate */
               if( SCIPisLT(scip, objgain, bestobjgain) || (SCIPisEQ(scip, objgain, bestobjgain) && frac < bestfrac) )
               {
                  bestcand = c;
                  bestobjgain = objgain;
                  bestfrac = frac;
                  bestcandmayrounddown = mayrounddown;
                  bestcandmayroundup = mayroundup;
                  bestcandroundup = roundup;
               }
            }
         }
         else
         {
            /* the candidate may not be rounded */
            if( frac < 0.5 )
               roundup = FALSE;
            else
            {
               roundup = TRUE;
               frac = 1.0 - frac;
            }

            /* penalize too small fractions */
            if( frac < 0.01 )
               frac += 10.0;

            /* prefer decisions on binary variables */
            if( !SCIPvarIsBinary(var) )
               frac *= 1000.0;

            /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
            if( bestcandmayrounddown || bestcandmayroundup || frac < bestfrac )
            {
               bestcand = c;
               bestfrac = frac;
               bestcandmayrounddown = FALSE;
               bestcandmayroundup = FALSE;
               bestcandroundup = roundup;
            }
            assert(bestfrac < SCIP_INVALID);
         }
      }
      assert(bestcand != -1);

      /* if all candidates are roundable, try to round the solution */
      if( bestcandmayrounddown || bestcandmayroundup )
      {
         SCIP_Bool success;

         /* create solution from diving LP and try to round it */
         SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
         SCIP_CALL( SCIProundSol(scip, heurdata->sol, &success) );

         if( success )
         {
            SCIPdebugMessage("fracdiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

            /* try to add solution to SCIP */
            SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, &success) );

            /* check, if solution was feasible and good enough */
            if( success )
            {
               SCIPdebugMessage(" -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;
            }
         }
      }

      var = lpcands[bestcand];

      backtracked = FALSE;
      do
      {
         /* if the variable is already fixed or if the solution value is outside the domain, numerical troubles may have
          * occured or variable was fixed by propagation while backtracking => Abort diving!
          */
         if( SCIPvarGetLbLocal(var) >= SCIPvarGetUbLocal(var) - 0.5 )
         {
            SCIPdebugMessage("Selected variable <%s> already fixed to [%g,%g] (solval: %.9f), diving aborted \n",
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), lpcandssol[bestcand]);
            cutoff = TRUE;
            break;
         }
         if( SCIPisFeasLT(scip, lpcandssol[bestcand], SCIPvarGetLbLocal(var)) || SCIPisFeasGT(scip, lpcandssol[bestcand], SCIPvarGetUbLocal(var)) )
         {
            SCIPdebugMessage("selected variable's <%s> solution value is outside the domain [%g,%g] (solval: %.9f), diving aborted\n",
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), lpcandssol[bestcand]);
            assert(backtracked);
            break;
         }

         /* apply rounding of best candidate */
         if( bestcandroundup == !backtracked )
         {
            /* round variable up */
            SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", pricerounds %d/%d: var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations, totalpricerounds, maxpricerounds,
               SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
               lpcandssol[bestcand], SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
               SCIPfeasCeil(scip, lpcandssol[bestcand]), SCIPvarGetUbLocal(var));
            SCIP_CALL( SCIPchgVarLbProbing(scip, var, SCIPfeasCeil(scip, lpcandssol[bestcand])) );
         }
         else
         {
            /* round variable down */
            SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", pricerounds %d/%d: var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations, totalpricerounds, maxpricerounds,
               SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
               lpcandssol[bestcand], SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
               SCIPvarGetLbLocal(var), SCIPfeasFloor(scip, lpcandssol[bestcand]));
            SCIP_CALL( SCIPchgVarUbProbing(scip, var, SCIPfeasFloor(scip, lpcandssol[bestcand])) );
         }

         /* apply domain propagation */
         SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );
         if( !cutoff )
         {
            /* resolve the diving LP */
            /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
             * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
             */
#ifdef NDEBUG
            SCIP_RETCODE retstat;
            nlpiterations = SCIPgetNLPIterations(scip);
            npricerounds = SCIPgetNPriceRounds(scip);
            if( maxpricerounds == 0 )
               retstat = SCIPsolveProbingLP(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror);
            else
               retstat = SCIPsolveProbingLPWithPricing(scip, FALSE, TRUE, maxpricerounds == -1 ? -1 : maxpricerounds - totalpricerounds, &lperror);
            if( retstat != SCIP_OKAY )
            {
               SCIPwarningMessage(scip, "Error while solving LP in Masterfracdiving heuristic; LP solve terminated with code <%d>\n",retstat);
            }
#else
            nlpiterations = SCIPgetNLPIterations(scip);
            npricerounds = SCIPgetNPriceRounds(scip);
            if( maxpricerounds == 0 )
               SCIP_CALL( SCIPsolveProbingLP(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror) );
            else
               SCIP_CALL( SCIPsolveProbingLPWithPricing(scip, FALSE, TRUE, maxpricerounds == -1 ? -1 : maxpricerounds - totalpricerounds, &lperror) );
#endif

            if( lperror )
               break;

            /* update iteration count */
            heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;
            heurdata->npricerounds += SCIPgetNPriceRounds(scip) - npricerounds;
            totalpricerounds += SCIPgetNPriceRounds(scip) - npricerounds;

            /* get LP solution status, objective value, and fractional variables, that should be integral */
            lpsolstat = SCIPgetLPSolstat(scip);
            cutoff = (lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT || lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
         }

         /* perform backtracking if a cutoff was detected */
         if( cutoff && !backtracked && heurdata->backtrack )
         {
            SCIPdebugMessage("  *** cutoff detected at level %d - backtracking\n", SCIPgetProbingDepth(scip));
            SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );
            SCIP_CALL( SCIPnewProbingNode(scip) );
            backtracked = TRUE;
         }
         else
            backtracked = FALSE;
      }
      while( backtracked );

      if( !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new objective value */
         oldobjval = objval;
         objval = SCIPgetLPObjval(scip);

         /* update pseudo cost values */
         if( SCIPisGT(scip, objval, oldobjval) )
         {
            if( bestcandroundup )
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, var, 1.0-bestfrac,
                     objval - oldobjval, 1.0) );
            }
            else
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, var, 0.0-bestfrac,
                     objval - oldobjval, 1.0) );
            }
         }

         /* get new fractional variables */
         SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );

         /* update original LP solution */
         SCIP_CALL( GCGrelaxUpdateCurrentSol(origprob, &origfeas) );
         if( origfeas )
         {
            SCIPdebugMessage("   -> found feasible original solution\n");
         }
      }
      SCIPdebugMessage("   -> lpsolstat=%d, objval=%g/%g, nfrac=%d\n", lpsolstat, objval, searchbound, nlpcands);
   }

   /* check if a solution has been found */
   if( nlpcands == 0 && !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_Bool success;

      /* create solution from diving LP */
      SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
      SCIPdebugMessage("fracdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* end diving */
   SCIP_CALL( SCIPendProbing(scip) );

   if( *result == SCIP_FOUNDSOL )
      heurdata->nsuccess++;

   SCIPdebugMessage("fracdiving heuristic finished\n");

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the masterfracdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurMasterfracdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Masterfracdiving primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecMasterfracdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeMasterfracdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitMasterfracdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitMasterfracdiving) );

   /* masterfracdiving heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/masterfracdiving/minreldepth",
         "minimal relative depth to start diving",
         &heurdata->minreldepth, TRUE, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/masterfracdiving/maxreldepth",
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, TRUE, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/masterfracdiving/maxlpiterquot",
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, FALSE, DEFAULT_MAXLPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/masterfracdiving/maxlpiterofs",
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, FALSE, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/masterfracdiving/maxpricequot",
         "maximal fraction of pricing rounds compared to node pricing rounds",
         &heurdata->maxpricequot, FALSE, DEFAULT_MAXPRICEQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/masterfracdiving/maxpriceofs",
         "additional number of allowed pricing rounds (-1: no limit)",
         &heurdata->maxpriceofs, FALSE, DEFAULT_MAXPRICEOFS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/masterfracdiving/maxdiveubquot",
         "maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveubquot, TRUE, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/masterfracdiving/maxdiveavgquot",
         "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveavgquot, TRUE, DEFAULT_MAXDIVEAVGQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/masterfracdiving/maxdiveubquotnosol",
         "maximal UBQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveubquotnosol, TRUE, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/masterfracdiving/maxdiveavgquotnosol",
         "maximal AVGQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveavgquotnosol, TRUE, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/masterfracdiving/backtrack",
         "use one level of backtracking if infeasibility is encountered?",
         &heurdata->backtrack, FALSE, DEFAULT_BACKTRACK, NULL, NULL) );

   return SCIP_OKAY;
}

