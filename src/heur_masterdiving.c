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

/**@file   heur_masterdiving.c
 * @brief  LP diving heuristic on the master variables
 * @author Tobias Achterberg
 * @author Christian Puchert
 */
#define SCIP_STATISTIC
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_masterdiving.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"


#define HEUR_NAME             "masterdiving"
#define HEUR_DESC             "LP diving heuristic on the master variables"
#define HEUR_DISPCHAR         'm'
#define HEUR_PRIORITY         -1000600
#define HEUR_FREQ             1
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
#define DEFAULT_MAXPRICEROUNDS       30 /**< maximal number of allowed pricing rounds (-1: no limit) */
#define DEFAULT_USEFARKASONLY      TRUE /**< perform pricing only if infeasibility is encountered */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use backtracking (discrepancy search) if infeasibility is encountered? */
#define DEFAULT_MAXDISCREPANCY        2 /**< maximal discrepancy in limited discrepancy search */
#define DEFAULT_MAXDISCDEPTH          3 /**< maximal depth until which a limited discrepancy search is performed */
#define DEFAULT_VARSELRULE          'v' /**< which variable selection should be used? ('c'oefficient, 'f'ractionality,
                                         *   'l'inesearch, 'p'scost, 'v'eclen; '*': alternate between rules) */

#define ALLOWEDRULES              "cflpv" /**< possible variable selection rules */
#define MINLPITER                 10000 /**< minimal number of LP iterations allowed in each LP solving call */

#ifdef SCIP_STATISTIC
#define EVENTHDLR_NAME         "masterdiving"
#define EVENTHDLR_DESC         "event handler for masterdiving solution statistics"
#endif


/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             minreldepth;        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth;        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   int                   maxlpiterofs;       /**< additional number of allowed LP iterations */
   int                   maxpricerounds;     /**< maximal number of allowed pricing rounds (-1: no limit) */
   SCIP_Bool             usefarkasonly;      /**< perform pricing only if infeasibility is encountered */
   SCIP_Real             maxdiveubquot;      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveavgquot;     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Bool             backtrack;          /**< use one level of backtracking if infeasibility is encountered? */
   int                   maxdiscrepancy;     /**< maximal discrepancy in limited discrepancy search */
   int                   maxdiscdepth;       /**< maximal depth until which a limited discrepancy search is performed */
   char                  varselrule;         /**< which variable selection should be used? ('c'oefficient, 'f'ractionality,
                                              *   'l'inesearch, 'p'scost, 'v'eclen; '*': alternate between rules) */
   char                  currentrule;        /**< variable selection rule that is to be used at the next call */
   SCIP_Longint          nlpiterations;      /**< LP iterations used in this heuristic */
   int                   npricerounds;       /**< pricing rounds used in this heuristic */
   int                   nsuccess;           /**< number of runs that produced at least one feasible solution */

#ifdef SCIP_STATISTIC
   SCIP_Longint*         ncalls;             /**< number of calls per diving strategy                               */
   SCIP_Longint*         nsols;              /**< number of solutions                                               */
   SCIP_Longint*         nimpsols;           /**< number of improving solutions                                     */
   SCIP_Longint*         ndivesols;          /**< number of integral diving LP solutions                            */
   SCIP_Longint*         nimpdivesols;       /**< number of improving integral diving LP solutions                  */
   SCIP_Longint*         nroundsols;         /**< number of integral solutions that have been obtained by rounding  */
   SCIP_Longint*         nimproundsols;      /**< number of improving integral solutions obtained by rounding       */
   SCIP_Longint*         ndives;             /**< number of dives                                                   */
   SCIP_Longint*         nrulelpiters;       /**< number of diving LP iterations (per diving rule)                  */
   SCIP_Longint*         nrulepricerds;      /**< number of pricing rounds (per diving rule)                        */
   SCIP_Real*            bestprimalbds;      /**< objective value of best solution found by this heuristic          */
   SCIP_Bool*            bestsolrounded;     /**< was the best solution obtained by rounding?                       */
#endif
};

#ifdef SCIP_STATISTIC
/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_Bool             heurisrunning;      /**< is the masterdiving heuristic currently running? */
};
#endif


/*
 * local methods
 */

/** finds best candidate variable w.r.t. locking numbers:
 * - prefer variables that may not be rounded without destroying LP feasibility:
 *   - of these variables, round variable with least number of locks in corresponding direction
 * - if all remaining fractional variables may be rounded without destroying LP feasibility:
 *   - round variable with least number of locks in opposite of its feasible rounding direction
 * - binary variables are preferred
 */
static
SCIP_RETCODE chooseCoefVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            lpcands,            /**< array of LP fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of LP fractional variables solution values */
   SCIP_Real*            lpcandsfrac,        /**< array of LP fractional variables fractionalities */
   int                   nlpcands,           /**< number of LP fractional variables */
   int*                  tabulist,           /**< array of variables that must not be chosen */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround    /**< pointer to store whether best candidate is trivially roundable */
   )
{
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   int bestnviolrows;             /* number of violated rows for best candidate */
   SCIP_Real bestcandfrac;        /* fractionality of best candidate */
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(lpcands != NULL);
   assert(lpcandsfrac != NULL);
   assert(lpcandssol != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);

   bestcandmayrounddown = TRUE;
   bestcandmayroundup = TRUE;
   bestnviolrows = INT_MAX;
   bestcandfrac = SCIP_INVALID;

   /* get best candidate */
   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_VAR* var;

      int nviolrows;

      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;
      SCIP_Real frac;

      int i;

      var = lpcands[c];
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      frac = lpcandsfrac[c];

      /* if the variable is on the tabu list, do not choose it */
      for( i = 0; i < heurdata->maxdiscrepancy; ++i )
         if( tabulist[i] == SCIPvarGetProbindex(var) )
            break;
      if( i < heurdata->maxdiscrepancy )
         continue;

      if( mayrounddown || mayroundup )
      {
         /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
         if( bestcandmayrounddown || bestcandmayroundup )
         {
            frac = 1.0 - frac;
            nviolrows = SCIPvarGetNLocksUp(var);

            /* penalize too small fractions */
            if( frac < 0.01 )
               nviolrows *= 100;

            /* prefer decisions on binary variables */
            if( !SCIPvarIsBinary(var) )
               nviolrows *= 1000;

            /* check, if candidate is new best candidate */
            assert( (0.0 < frac && frac < 1.0) || SCIPvarIsBinary(var) );
            if( nviolrows + frac < bestnviolrows + bestcandfrac )
            {
               *bestcand = c;
               bestnviolrows = nviolrows;
               bestcandfrac = frac;
               bestcandmayrounddown = mayrounddown;
               bestcandmayroundup = mayroundup;
            }
         }
      }
      else
      {
         /* the candidate may not be rounded */
         frac = 1.0 - frac;
         nviolrows = SCIPvarGetNLocksUp(var);

         /* penalize too small fractions */
         if( frac < 0.01 )
            nviolrows *= 100;

         /* prefer decisions on binary variables */
         if( !SCIPvarIsBinary(var) )
            nviolrows *= 100;

         /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
         assert((0.0 < frac && frac < 1.0) || SCIPvarIsBinary(var));
         if( bestcandmayrounddown || bestcandmayroundup || nviolrows + frac < bestnviolrows + bestcandfrac )
         {
            *bestcand = c;
            bestnviolrows = nviolrows;
            bestcandfrac = frac;
            bestcandmayrounddown = FALSE;
            bestcandmayroundup = FALSE;
         }
         assert(bestcandfrac < SCIP_INVALID);
      }
   }

   *bestcandmayround = bestcandmayroundup || bestcandmayrounddown;

   return SCIP_OKAY;
}

/** finds best candidate variable w.r.t. fractionality:
 * - prefer variables that may not be rounded without destroying LP feasibility:
 *   - of these variables, round least fractional variable in corresponding direction
 * - if all remaining fractional variables may be rounded without destroying LP feasibility:
 *   - round variable with least increasing objective value
 * - binary variables are preferred
 */
static
SCIP_RETCODE chooseFracVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            lpcands,            /**< array of LP fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of LP fractional variables solution values */
   SCIP_Real*            lpcandsfrac,        /**< array of LP fractional variables fractionalities */
   int                   nlpcands,           /**< number of LP fractional variables */
   int*                  tabulist,           /**< array of variables that must not be chosen */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround    /**< pointer to store whether best candidate is trivially roundable */
   )
{
   SCIP_Real bestobjgain;
   SCIP_Real bestfrac;
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(lpcandsfrac != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);

   bestcandmayrounddown = TRUE;
   bestcandmayroundup = TRUE;
   bestobjgain = SCIPinfinity(scip);
   bestfrac = SCIP_INVALID;

   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_VAR* var;
      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;
      SCIP_Real frac;
      SCIP_Real obj;
      SCIP_Real objgain;

      int i;

      var = lpcands[c];

      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      frac = lpcandsfrac[c];
      obj = SCIPvarGetObj(var);

      /* if the variable is on the tabu list, do not choose it */
      for( i = 0; i < heurdata->maxdiscrepancy; ++i )
         if( tabulist[i] == SCIPvarGetProbindex(var) )
            break;
      if( i < heurdata->maxdiscrepancy )
         continue;

      if( mayrounddown || mayroundup )
      {
         /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
         if( bestcandmayrounddown || bestcandmayroundup )
         {
            objgain = (1.0-frac)*obj;

            /* penalize too small fractions */
            if( ABS(1.0 - frac) < 0.01 )
               objgain *= 1000.0;

            /* prefer decisions on binary variables */
            if( !SCIPvarIsBinary(var) )
               objgain *= 1000.0;

            /* check, if candidate is new best candidate */
            if( SCIPisLT(scip, objgain, bestobjgain) || (SCIPisEQ(scip, objgain, bestobjgain) && frac > bestfrac) )
            {
               *bestcand = c;
               bestobjgain = objgain;
               bestfrac = frac;
               bestcandmayrounddown = mayrounddown;
               bestcandmayroundup = mayroundup;
            }
         }
      }
      else
      {
         /* penalize too small fractions */
         if( ABS(1.0-frac) < 0.01 )
            frac += 10.0;

         /* prefer decisions on binary variables */
         if( !SCIPvarIsBinary(var) )
            frac *= 1000.0;

         /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
         if( bestcandmayrounddown || bestcandmayroundup || frac > bestfrac )
         {
            *bestcand = c;
            bestfrac = frac;
            bestcandmayrounddown = FALSE;
            bestcandmayroundup = FALSE;
         }
         assert(bestfrac < SCIP_INVALID);
      }
   }

   *bestcandmayround = bestcandmayroundup || bestcandmayrounddown;

   return SCIP_OKAY;
}

/** finds best candidate variable w.r.t. the incumbent solution:
 * - prefer variables that may not be rounded without destroying LP feasibility:
 *   - of these variables, round a variable to its value in direction of incumbent solution, and choose the
 *     variable that is closest to its rounded value
 * - if all remaining fractional variables may be rounded without destroying LP feasibility:
 *   - round variable in direction that destroys LP feasibility (other direction is checked by SCIProundSol())
 *   - round variable with least increasing objective value
 * - binary variables are preferred
 */
static
SCIP_RETCODE chooseGuidedVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            lpcands,            /**< array of LP fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of LP fractional variables solution values */
   SCIP_Real*            lpcandsfrac,        /**< array of LP fractional variables fractionalities */
   int                   nlpcands,           /**< number of LP fractional variables */
   int*                  tabulist,           /**< array of variables that must not be chosen */
   SCIP_SOL*             bestsol,            /**< incumbent solution */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround,   /**< pointer to store whether best candidate is trivially roundable */
   SCIP_Bool*            bestcandroundup     /**< pointer to store whether best candidate should be rounded up */
   )
{
   SCIP_Real bestobjgain;
   SCIP_Real bestfrac;
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(lpcands != NULL);
   assert(lpcandsfrac != NULL);
   assert(lpcandssol != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);
   assert(bestcandroundup != NULL);

   bestcandmayrounddown = TRUE;
   bestcandmayroundup = TRUE;
   bestobjgain = SCIPinfinity(scip);
   bestfrac = SCIP_INVALID;

   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_VAR* var;
      SCIP_Real bestsolval;
      SCIP_Real solval;
      SCIP_Real obj;
      SCIP_Real frac;
      SCIP_Real objgain;

      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;
      SCIP_Bool roundup;

      int i;

      var = lpcands[c];
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      solval = lpcandssol[c];
      frac = lpcandsfrac[c];
      obj = SCIPvarGetObj(var);
      bestsolval = SCIPgetSolVal(scip, bestsol, var);

      /* if the variable is on the tabu list, do not choose it */
      for( i = 0; i < heurdata->maxdiscrepancy; ++i )
         if( tabulist[i] == SCIPvarGetProbindex(var) )
            break;
      if( i < heurdata->maxdiscrepancy )
         continue;

      /* select default rounding direction */
      roundup = (solval < bestsolval);

      if( mayrounddown || mayroundup )
      {
         /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
         if( bestcandmayrounddown || bestcandmayroundup )
         {
            /* choose rounding direction:
             * - if variable may be rounded in both directions, round corresponding to its value in incumbent solution
             * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
             *   the current fractional solution with SCIProundSol()
             */
            if( !mayrounddown || !mayroundup )
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
               *bestcand = c;
               bestobjgain = objgain;
               bestfrac = frac;
               bestcandmayrounddown = mayrounddown;
               bestcandmayroundup = mayroundup;
               *bestcandroundup = roundup;
            }
         }
      }
      else
      {
         /* the candidate may not be rounded */
         if( roundup )
            frac = 1.0 - frac;

         /* penalize too small fractions */
         if( frac < 0.01 )
            frac += 10.0;

         /* prefer decisions on binary variables */
         if( !SCIPvarIsBinary(var) )
            frac *= 1000.0;

         /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
         if( bestcandmayrounddown || bestcandmayroundup || frac < bestfrac )
         {
            *bestcand = c;
            bestfrac = frac;
            bestcandmayrounddown = FALSE;
            bestcandmayroundup = FALSE;
            *bestcandroundup = roundup;
         }
      }
   }

   *bestcandmayround = bestcandmayroundup || bestcandmayrounddown;

   return SCIP_OKAY;
}

/** finds best candidate variable w.r.t. the root LP solution:
 * - in the projected space of fractional variables, extend the line segment connecting the root solution and
 *   the current LP solution up to the point, where one of the fractional variables becomes integral
 * - round this variable to the integral value
 */
static
SCIP_RETCODE chooseLinesearchVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            lpcands,            /**< array of LP fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of LP fractional variables solution values */
   SCIP_Real*            lpcandsfrac,        /**< array of LP fractional variables fractionalities */
   int                   nlpcands,           /**< number of LP fractional variables */
   int*                  tabulist,           /**< array of variables that must not be chosen */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround    /**< pointer to store whether best candidate is trivially roundable */
   )
{
   SCIP_Real bestdistquot;
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(lpcandsfrac != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);

   *bestcandmayround = TRUE;
   bestdistquot = SCIPinfinity(scip);

   /* get best candidate */
   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_VAR* var;
      SCIP_Real solval;
      SCIP_Real rootsolval;
      SCIP_Real distquot;

      int i;

      var = lpcands[c];
      solval = lpcandssol[c];
      rootsolval = SCIPvarGetRootSol(var);

      /* if the variable is on the tabu list, do not choose it */
      for( i = 0; i < heurdata->maxdiscrepancy; ++i )
         if( tabulist[i] == SCIPvarGetProbindex(var) )
            break;
      if( i < heurdata->maxdiscrepancy )
         continue;

      if( SCIPisGT(scip, solval, rootsolval) )
      {
         distquot = (SCIPfeasCeil(scip, solval) - solval) / (solval - rootsolval);

         /* avoid roundable candidates */
         if( SCIPvarMayRoundUp(var) )
            distquot *= 1000.0;
      }
      else
         distquot = SCIPinfinity(scip);

      /* check whether the variable is roundable */
      *bestcandmayround = *bestcandmayround && (SCIPvarMayRoundDown(var) || SCIPvarMayRoundUp(var));

      /* check, if candidate is new best candidate */
      if( distquot < bestdistquot )
      {
         *bestcand = c;
         bestdistquot = distquot;
      }
   }

   return SCIP_OKAY;
}

/** calculates the pseudocost score for a given variable w.r.t. a given solution value */
static
void calcPscostQuot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             primsol,            /**< primal solution of variable */
   SCIP_Real             frac,               /**< fractionality of variable */
   SCIP_Real*            pscostquot          /**< pointer to store pseudo cost quotient */
   )
{
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   assert(pscostquot != NULL);
   assert(SCIPisEQ(scip, frac, primsol - SCIPfeasFloor(scip, primsol)));

   /* bound fractions to not prefer variables that are nearly integral */
   frac = MAX(frac, 0.1);
   frac = MIN(frac, 0.9);

   /* get pseudo cost quotient */
   pscostdown = SCIPgetVarPseudocostVal(scip, var, 0.0-frac);
   pscostup = SCIPgetVarPseudocostVal(scip, var, 1.0-frac);
   assert(pscostdown >= 0.0 && pscostup >= 0.0);

   /* calculate pseudo cost quotient */
   *pscostquot = sqrt(frac) * (1.0+pscostdown) / (1.0+pscostup);

   /* reward or punish variables:
    *  - a variable which has moved downwards from its root LP value should not be rounded up,
    *    hence its score is decreased; the same is done for variables which are near to their
    *    rounded down values
    *  - on the other hand, increase the scores of variables that have moved upwards from
    *    their root LP value or which are near to their rounded up values
    */
   if( primsol < SCIPvarGetRootSol(var) - 0.4 )
      (*pscostquot) /= 100.0;
   else if( primsol > SCIPvarGetRootSol(var) + 0.4 )
      (*pscostquot) *= 100.0;
   else if( frac < 0.3 )
      (*pscostquot) /= 100.0;
   else if( frac > 0.7 )
      (*pscostquot) *= 100.0;

   /* prefer decisions on binary variables */
   if( SCIPvarIsBinary(var) )
      (*pscostquot) *= 1000.0;
}

/** finds best candidate variable w.r.t. pseudo costs:
 * - prefer variables that may not be rounded without destroying LP feasibility:
 *   - of these variables, round variable with largest rel. difference of pseudo cost values in corresponding
 *     direction
 * - if all remaining fractional variables may be rounded without destroying LP feasibility:
 *   - round variable in the objective value direction
 * - binary variables are preferred
 */
static
SCIP_RETCODE choosePscostVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            lpcands,            /**< array of LP fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of LP fractional variables solution values */
   SCIP_Real*            lpcandsfrac,        /**< array of LP fractional variables fractionalities */
   int                   nlpcands,           /**< number of LP fractional variables */
   int*                  tabulist,           /**< array of variables that must not be chosen */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround    /**< pointer to store whether best candidate is trivially roundable */
   )
{
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   SCIP_Real bestpscostquot;
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(lpcands != NULL);
   assert(lpcandsfrac != NULL);
   assert(lpcandssol != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);

   bestcandmayrounddown = TRUE;
   bestcandmayroundup = TRUE;
   bestpscostquot = -1.0;

   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_VAR* var;
      SCIP_Real primsol;

      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;
      SCIP_Real frac;
      SCIP_Real pscostquot;

      int i;

      var = lpcands[c];
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      primsol = lpcandssol[c];
      frac = lpcandsfrac[c];
      pscostquot = SCIP_INVALID;

      /* if the variable is on the tabu list, do not choose it */
      for( i = 0; i < heurdata->maxdiscrepancy; ++i )
         if( tabulist[i] == SCIPvarGetProbindex(var) )
            break;
      if( i < heurdata->maxdiscrepancy )
         continue;

      if( mayrounddown || mayroundup )
      {
         /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
         if( bestcandmayrounddown || bestcandmayroundup )
         {
            /* calculate pseudo cost */
            calcPscostQuot(scip, var, primsol, frac, &pscostquot);
            assert(!SCIPisInfinity(scip,ABS(pscostquot)));

            /* check, if candidate is new best candidate */
            if( pscostquot > bestpscostquot )
            {
               *bestcand = c;
               bestpscostquot = pscostquot;
               bestcandmayrounddown = mayrounddown;
               bestcandmayroundup = mayroundup;
            }
         }
      }
      else
      {
         /* the candidate may not be rounded: calculate pseudo cost */
         calcPscostQuot(scip, var, primsol, frac, &pscostquot);
         assert(!SCIPisInfinity(scip, ABS(pscostquot)));

         /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
         if( bestcandmayrounddown || bestcandmayroundup || pscostquot > bestpscostquot )
         {
            *bestcand = c;
            bestpscostquot = pscostquot;
            bestcandmayrounddown = FALSE;
            bestcandmayroundup = FALSE;
         }
      }
   }

   *bestcandmayround = bestcandmayroundup || bestcandmayrounddown;

   return SCIP_OKAY;
}

/** finds best candidate variable w.r.t. vector length:
 * - round variables in direction where objective value gets worse; for zero objective coefficient, round upwards
 * - round variable with least objective value deficit per row the variable appears in
 *   (we want to "fix" as many rows as possible with the least damage to the objective function)
 */
static
SCIP_RETCODE chooseVeclenVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            lpcands,            /**< array of NLP fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of NLP fractional variables solution values */
   SCIP_Real*            lpcandsfrac,        /**< array of NLP fractional variables fractionalities */
   int                   nlpcands,           /**< number of NLP fractional variables */
   int*                  tabulist,           /**< array of variables that must not be chosen */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround    /**< pointer to store whether best candidate is trivially roundable */
   )
{
   SCIP_Real bestscore;
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(lpcands != NULL);
   assert(lpcandsfrac != NULL);
   assert(lpcandssol != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);

   *bestcandmayround = TRUE;
   bestscore = SCIP_REAL_MAX;

   /* get best candidate */
   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_VAR* var;

      SCIP_Real obj;
      SCIP_Real objdelta;
      SCIP_Real frac;
      SCIP_Real score;
      int colveclen;

      int i;

      var = lpcands[c];

      /* if the variable is on the tabu list, do not choose it */
      for( i = 0; i < heurdata->maxdiscrepancy; ++i )
         if( tabulist[i] == SCIPvarGetProbindex(var) )
            break;
      if( i < heurdata->maxdiscrepancy )
         continue;

      frac = lpcandsfrac[c];
      obj = SCIPvarGetObj(var);
      objdelta = (1.0 - frac) * obj;

      colveclen = (SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN ? SCIPcolGetNNonz(SCIPvarGetCol(var)) : 0);

      /* check whether the variable is roundable */
      *bestcandmayround = *bestcandmayround && (SCIPvarMayRoundDown(var) || SCIPvarMayRoundUp(var));

      /* smaller score is better */
      score = (objdelta + SCIPsumepsilon(scip))/((SCIP_Real)colveclen+1.0);

      /* penalize negative scores (i.e. improvements in the objective) */
      if( score <= 0.0 )
         score *= 100.0;

      /* prefer decisions on binary variables */
      if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
         score *= 1000.0;

      /* check, if candidate is new best candidate */
      if( score < bestscore )
      {
         *bestcand = c;
         bestscore = score;
      }
   }

   return SCIP_OKAY;
}

/** finds the best candidate variable for diving */
static
SCIP_RETCODE chooseVariable(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            lpcands,            /**< array of LP fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of LP fractional variables solution values */
   SCIP_Real*            lpcandsfrac,        /**< array of LP fractional variables fractionalities */
   int                   nlpcands,           /**< number of LP fractional variables */
   int*                  tabulist,           /**< array of variables that must not be chosen */
   SCIP_SOL*             bestsol,            /**< incumbent solution */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround,   /**< pointer to store whether best candidate is trivially roundable */
   SCIP_Bool*            bestcandroundup     /**< pointer to store whether best candidate should be rounded up */
   )
{
   switch( heurdata->currentrule )
   {
   case 'c':
      SCIP_CALL( chooseCoefVar(scip, heurdata, lpcands, lpcandssol, lpcandsfrac, nlpcands,
            tabulist, bestcand, bestcandmayround) );
      break;
   case 'f':
      SCIP_CALL( chooseFracVar(scip, heurdata, lpcands, lpcandssol, lpcandsfrac, nlpcands,
            tabulist, bestcand, bestcandmayround) );
      break;
   case 'g':
      SCIP_CALL( chooseGuidedVar(scip, heurdata, lpcands, lpcandssol, lpcandsfrac, nlpcands,
            tabulist, bestsol, bestcand, bestcandmayround, bestcandroundup) );
      break;
   case 'l':
      SCIP_CALL( chooseLinesearchVar(scip, heurdata, lpcands, lpcandssol, lpcandsfrac, nlpcands,
            tabulist, bestcand, bestcandmayround) );
      break;
   case 'p':
      SCIP_CALL( choosePscostVar(scip, heurdata, lpcands, lpcandssol, lpcandsfrac, nlpcands,
            tabulist, bestcand, bestcandmayround) );
      break;
   case 'v':
      SCIP_CALL( chooseVeclenVar(scip, heurdata, lpcands, lpcandssol, lpcandsfrac, nlpcands,
            tabulist, bestcand, bestcandmayround) );
      break;
   default:
      SCIPerrorMessage("invalid variable selection rule\n");
      return SCIP_INVALIDDATA;
   }

   *bestcandroundup = TRUE;
   return SCIP_OKAY;
}

/** gets to the variable selection rule for the next call of this heuristic */
static
SCIP_RETCODE getNextRule(
      SCIP*                 scip,               /**< original SCIP data structure */
      SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
      )
{
   if( heurdata->varselrule == '*' )
   {
      const char* rules = ALLOWEDRULES;
      int nrules = strlen(rules);
      int i;

      assert(nrules > 0);

      for( i = 0; i < nrules; ++i )
         if( rules[i] == heurdata->currentrule )
            break;

      assert(i < nrules);

      if( i == nrules-1 )
         heurdata->currentrule = rules[0];
      else
         heurdata->currentrule = rules[i+1];
   }
#ifndef NDEBUG
   else
   {
      assert(heurdata->currentrule == heurdata->varselrule);
   }
#endif

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeMasterdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
#ifdef SCIP_STATISTIC
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;
#endif
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

#ifdef SCIP_STATISTIC
   /* free event handler data */
   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   SCIPfreeMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(eventhdlr, NULL);
#endif

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitMasterdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   const char* rules;
   int nrules;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   if ( nrules == 0 )
   {
      SCIPerrorMessage("no valid variable selection rule found!\n");
      return SCIP_INVALIDDATA;
   }

   /* initialize data */
   if( heurdata->varselrule == '*' )
      heurdata->currentrule = rules[0];
   else
      heurdata->currentrule = heurdata->varselrule;

   heurdata->nlpiterations = 0;
   heurdata->nsuccess = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitMasterdiving) /*lint --e{715}*/
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


#ifdef SCIP_STATISTIC
/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolMasterdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   const char* rules;
   int nrules;

   int i;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   /* initialize statistical data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->ncalls, nrules) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->nsols, nrules) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->nimpsols, nrules) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->ndivesols, nrules) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->nimpdivesols, nrules) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->nroundsols, nrules) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->nimproundsols, nrules) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->ndives, nrules) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->nrulelpiters, nrules) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->nrulepricerds, nrules) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->bestprimalbds, nrules) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->bestsolrounded, nrules) );

   for( i = 0; i < nrules; ++i )
   {
      heurdata->ncalls[i] = 0;
      heurdata->nsols[i] = 0;
      heurdata->nimpsols[i] = 0;
      heurdata->ndivesols[i] = 0;
      heurdata->nimpdivesols[i] = 0;
      heurdata->nroundsols[i] = 0;
      heurdata->nimproundsols[i] = 0;
      heurdata->ndives[i] = 0;
      heurdata->nrulelpiters[i] = 0;
      heurdata->nrulepricerds[i] = 0;
      heurdata->bestprimalbds[i] = SCIPinfinity(scip);
      heurdata->bestsolrounded[i] = FALSE;
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolMasterdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   const char* rules;
   int nrules;

   int i;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   /* print detailed statistics */
   SCIPstatisticPrintf("Master Diving Heuristics   :      Calls       Sols  Improving   DiveSols  Improving  RoundSols  Improving      Dives   LP iters  Price rds    BestPrimal Rounded?\n");
   for( i = 0; i < nrules; ++i )
   {
      SCIPstatisticPrintf("%c                          : %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT,
         rules[i], heurdata->ncalls[i], heurdata->nsols[i], heurdata->nimpsols[i], heurdata->ndivesols[i], heurdata->nimpdivesols[i], heurdata->nroundsols[i], heurdata->nimproundsols[i], heurdata->ndives[i], heurdata->nrulelpiters[i], heurdata->nrulepricerds[i]);
      if( SCIPisInfinity(scip, heurdata->bestprimalbds[i]) )
         SCIPstatisticPrintf("      infinity");
      else
         SCIPstatisticPrintf(" %13.6e", heurdata->bestprimalbds[i]);
      SCIPstatisticPrintf(heurdata->bestsolrounded[i] ? "      yes\n" : "       no\n");
   }
   SCIPstatisticPrintf("\n");

   SCIPfreeMemoryArray(scip, &heurdata->bestsolrounded);
   SCIPfreeMemoryArray(scip, &heurdata->bestprimalbds);
   SCIPfreeMemoryArray(scip, &heurdata->nrulepricerds);
   SCIPfreeMemoryArray(scip, &heurdata->nrulelpiters);
   SCIPfreeMemoryArray(scip, &heurdata->ndives);
   SCIPfreeMemoryArray(scip, &heurdata->nimproundsols);
   SCIPfreeMemoryArray(scip, &heurdata->nroundsols);
   SCIPfreeMemoryArray(scip, &heurdata->nimpdivesols);
   SCIPfreeMemoryArray(scip, &heurdata->ndivesols);
   SCIPfreeMemoryArray(scip, &heurdata->nimpsols);
   SCIPfreeMemoryArray(scip, &heurdata->nsols);
   SCIPfreeMemoryArray(scip, &heurdata->ncalls);

   return SCIP_OKAY;
}
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecMasterdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP* origprob;
#ifdef SCIP_STATISTIC
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
#endif
   SCIP_HEURDATA* heurdata;
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_SOL* bestsol;
   SCIP_VAR* var;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int* discrepancies;
   int* selectedvars;
   int* tabulist;
   const char* rules;
   SCIP_Real searchubbound;
   SCIP_Real searchavgbound;
   SCIP_Real searchbound;
   SCIP_Real objval;
   SCIP_Real oldobjval;
   SCIP_Real bestfrac;
   SCIP_Bool bestcandmayround;
   SCIP_Bool bestcandroundup;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;
   SCIP_Bool backtracked;
   SCIP_Bool farkaspricing;
   SCIP_Bool origfeas;
   SCIP_Longint ncalls;
   SCIP_Longint nsolsfound;
   SCIP_Longint nlpiterations;         /* lp iterations performed in one single diving loop */
   SCIP_Longint maxnlpiterations;
   int npricerounds;                   /* pricing rounds performed in one single diving loop */
   int totalpricerounds;               /* pricing rounds performed in one call of the heuristic */
   int nlpcands;
   int startnlpcands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int bestcand;
   int nrules;
   int ruleindex;

#ifdef SCIP_STATISTIC
   /* variable declarations for additional statistics */
   int ndives;                         /* diving loops performed in one call of the heuristic */
   SCIP_Longint totallpiters;          /* lp iterations performed in one call of the heuristic */
#endif

   int i;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

#ifdef SCIP_STATISTIC
   /* get the masterdiving event handler and its data */
   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
#endif

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
//      SCIPdebugMessage("not executing Masterdiving heuristic: master LP not solved to optimality\n");
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

   /* get fractional variables that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );

   /* don't try to dive, if there are no fractional variables */
   if( nlpcands == 0 )
      return SCIP_OKAY;

   /* store a copy of the best solution, if guided diving should be used */
   bestsol = NULL;
   if( heurdata->currentrule == 'g' )
   {
      /* do not perform guided diving if no feasible solutions exist or if the best solution lives
       * in the original variable space (we then cannot use it since it might violate the global
       * bounds of the current problem)
       */
      if( SCIPgetNSols(scip) == 0 || SCIPsolIsOriginal(SCIPgetBestSol(scip)) )
         if( heurdata->varselrule == '*' )
         {
            SCIP_CALL( getNextRule(scip, heurdata) );
         }
         else
            return SCIP_OKAY;
      else
         SCIP_CALL( SCIPcreateSolCopy(scip, &bestsol, SCIPgetBestSol(scip)) );
   }

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

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &discrepancies, heurdata->maxdiscdepth) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tabulist, heurdata->maxdiscrepancy) );
   SCIP_CALL( SCIPallocBufferArray(scip, &selectedvars, heurdata->maxdiscdepth) );

   /* initialize arrays */
   for( i = 0; i < heurdata->maxdiscdepth; ++i )
      discrepancies[i] = 0;
   for( i = 0; i < heurdata->maxdiscrepancy; ++i )
      tabulist[i] = -1;

   /* get the index of the current variable selection rule in the ALLOWEDRULES string */
   for( ruleindex = 0; ruleindex < nrules; ++ruleindex )
      if( rules[ruleindex] == heurdata->currentrule )
         break;
   assert(ruleindex < nrules);


   *result = SCIP_DIDNOTFIND;

#ifdef SCIP_STATISTIC
   eventhdlrdata->heurisrunning = TRUE;
   ++heurdata->ncalls[ruleindex];
#endif

   /* start diving */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   SCIPenableVarHistory(scip);

   /* get LP objective value*/
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   objval = SCIPgetLPObjval(scip);

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") executing masterdiving heuristic: depth=%d, %d fractionals, dualbound=%g, searchbound=%g, divingrule=%c\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), SCIPretransformObj(scip, searchbound), heurdata->currentrule);

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
   startnlpcands = nlpcands;

#ifdef SCIP_STATISTIC
   ndives = 0;
   totallpiters = 0;
#endif

   while( !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && (divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations && objval < searchbound))
      && (divedepth >= heurdata->maxdiscdepth || discrepancies[divedepth] <= heurdata->maxdiscrepancy)
      && !SCIPisStopped(scip) )
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      divedepth++;

#ifdef SCIP_STATISTIC
      ++heurdata->ndives[ruleindex];
      ++ndives;
#endif

      bestcand = -1;
      bestfrac = SCIP_INVALID;
      bestcandmayround = TRUE;
      bestcandroundup = FALSE;

      /* choose a variable to dive on */
      SCIP_CALL( chooseVariable(scip, heurdata, lpcands, lpcandssol, lpcandsfrac, nlpcands, tabulist,
         bestsol, &bestcand, &bestcandmayround, &bestcandroundup) );

      /* if no variable could be chosen, abort diving */
      if( bestcand == -1 )
      {
         SCIPdebugMessage("No variable for diving could be selected, diving aborted\n");
         break;
      }

      assert(bestcand >= 0);
      var = lpcands[bestcand];
      bestfrac = lpcandsfrac[bestcand];

      /* memorize selected variables up to the maximal depth for discrepancy search */
      if( divedepth-1 < heurdata->maxdiscdepth )
         selectedvars[divedepth-1] = SCIPvarGetProbindex(var);

      /* if all candidates are roundable, try to round the solution */
      if( bestcandmayround )
      {
         SCIP_Bool success;

         /* create solution from diving LP and try to round it */
         SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
         SCIP_CALL( SCIProundSol(scip, heurdata->sol, &success) );

         if( success )
         {
            SCIPdebugMessage("masterdiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

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

      /* round variable up */
      SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", pricerounds %d/%d: var <%s>, round=%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
         divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations, totalpricerounds, heurdata->maxpricerounds,
         SCIPvarGetName(var), bestcandmayround,
         lpcandssol[bestcand], SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
         SCIPfeasCeil(scip, lpcandssol[bestcand]), SCIPvarGetUbLocal(var));
      SCIP_CALL( SCIPchgVarLbProbing(scip, var, SCIPfeasCeil(scip, lpcandssol[bestcand])) );

      backtracked = FALSE;
      farkaspricing = FALSE;
      do
      {
         /* apply domain propagation */
         SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );

         if( !cutoff || backtracked || farkaspricing )
         {
            /* resolve the diving LP */
            /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
             * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
             */
#ifdef NDEBUG
            SCIP_RETCODE retstat;
            nlpiterations = SCIPgetNLPIterations(scip);
            npricerounds = SCIPgetNPriceRounds(scip);
            if( (!heurdata->usefarkasonly || farkaspricing)
               && (heurdata->maxpricerounds == -1 || totalpricerounds < heurdata->maxpricerounds) )
               retstat = SCIPsolveProbingLPWithPricing(scip, FALSE, TRUE, heurdata->maxpricerounds == -1 ? -1 : heurdata->maxpricerounds - totalpricerounds, &lperror);

            else
               retstat = SCIPsolveProbingLP(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror);
            if( retstat != SCIP_OKAY )
            {
               SCIPwarningMessage(scip, "Error while solving LP in Masterdiving heuristic; LP solve terminated with code <%d>\n",retstat);
            }
#else
            nlpiterations = SCIPgetNLPIterations(scip);
            npricerounds = SCIPgetNPriceRounds(scip);
            if( (!heurdata->usefarkasonly || farkaspricing)
               && (heurdata->maxpricerounds == -1 || totalpricerounds < heurdata->maxpricerounds) )
               SCIP_CALL( SCIPsolveProbingLPWithPricing(scip, FALSE, TRUE, heurdata->maxpricerounds == -1 ? -1 : heurdata->maxpricerounds - totalpricerounds, &lperror) );
            else
               SCIP_CALL( SCIPsolveProbingLP(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror) );
#endif

            if( lperror )
               break;

            /* update iteration counts */
            heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;
            heurdata->npricerounds += SCIPgetNPriceRounds(scip) - npricerounds;
            totalpricerounds += SCIPgetNPriceRounds(scip) - npricerounds;
#ifdef SCIP_STATISTIC
            SCIPstatistic( totallpiters += SCIPgetNLPIterations(scip) - nlpiterations );

            /* update summarized statistics */
            heurdata->nrulelpiters[ruleindex] += SCIPgetNLPIterations(scip) - nlpiterations;
            heurdata->nrulepricerds[ruleindex] += SCIPgetNPriceRounds(scip) - npricerounds;
#endif

            /* get LP solution status */
            lpsolstat = SCIPgetLPSolstat(scip);
            cutoff = (lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT || lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
         }

         /* if infeasibility is encountered, perform Farkas pricing
          * in order to reach feasibility again */
         if( lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE && heurdata->usefarkasonly
            && !farkaspricing && (heurdata->maxpricerounds == -1 || totalpricerounds < heurdata->maxpricerounds)
            && !backtracked )
         {
            SCIPdebugMessage("  *** infeasibility detected at level %d - perform Farkas pricing\n", SCIPgetProbingDepth(scip));
            farkaspricing = TRUE;
         }
         else
            farkaspricing = FALSE;

         /* perform backtracking if a cutoff was detected */
         if( cutoff && !backtracked && heurdata->backtrack && !farkaspricing )
         {
            SCIPdebugMessage("  *** cutoff detected at level %d - backtracking\n", SCIPgetProbingDepth(scip));

            /* go back until the search can differ from the previous search tree */
            do
            {
               SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );
               --divedepth;
            }
            while( divedepth > 0
               && (divedepth >= heurdata->maxdiscdepth
                  || discrepancies[divedepth] >= heurdata->maxdiscrepancy) );

            assert(divedepth < heurdata->maxdiscdepth);

            /* add variable selected previously at this depth to the tabu list */
            tabulist[discrepancies[divedepth]] = selectedvars[divedepth];

            ++discrepancies[divedepth];
            for( i = divedepth + 1; i < heurdata->maxdiscdepth; ++i )
               discrepancies[i] = discrepancies[divedepth];

            backtracked = TRUE;
         }
         else
            backtracked = FALSE;
      }
      while( backtracked || farkaspricing );

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
            SCIPdebugMessage("   -> LP solution is feasible in the original problem\n");
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
      SCIPdebugMessage("masterdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

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

#ifdef SCIP_STATISTIC
   eventhdlrdata->heurisrunning = FALSE;
#endif

   if( *result == SCIP_FOUNDSOL )
      heurdata->nsuccess++;

#ifdef SCIP_STATISTIC
   SCIPstatisticPrintf("Masterdiving statistic: rule %c, %3d diveloops, %"SCIP_LONGINT_FORMAT" lp iterations, %5d pricing rounds\n",
      heurdata->currentrule, ndives, totallpiters, totalpricerounds);
#endif

   /* free memory */
   SCIPfreeBufferArray(scip, &selectedvars);
   SCIPfreeBufferArray(scip, &tabulist);
   SCIPfreeBufferArray(scip, &discrepancies);

   /* free copied best solution */
   if( heurdata->currentrule == 'g' )
   {
      assert(bestsol != NULL);
      SCIP_CALL( SCIPfreeSol(scip, &bestsol) );
   }
   else
      assert(bestsol == NULL);

   SCIP_CALL( getNextRule(scip, heurdata) );

   SCIPdebugMessage("masterdiving heuristic finished\n");

   return SCIP_OKAY;
}


#ifdef SCIP_STATISTIC
/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitMasterdiving)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(eventhdlr != NULL);

   /* get event handler data */
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->heurisrunning = FALSE;

   /* notify GCG that this event should catch the SOLFOUND event */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitMasterdiving)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);

   /* notify GCG that this event should drop the SOLFOUND event */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}


/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecMasterdiving)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;
   SCIP_SOL* sol;
   SCIP_HEUR* solheur;
   const char* rules;
   int nrules;
   int ruleindex;

   assert(eventhdlr != NULL);

   /* get event handler data */
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* get masterdiving primal heuristic */
   heur = SCIPfindHeur(scip, HEUR_NAME);
   assert(heur != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   /* get new primal solution */
   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   /* get heuristic that found the solution */
   solheur = SCIPgetSolHeur(scip, sol);

   /* get the index of the current variable selection rule in the ALLOWEDRULES string */
   for( ruleindex = 0; ruleindex < nrules; ++ruleindex )
      if( rules[ruleindex] == heurdata->currentrule )
         break;
   assert(ruleindex < nrules);

   /* if the heuristic is currently running, update its solution statistics */
   if( eventhdlrdata->heurisrunning )
   {
      ++heurdata->nsols[ruleindex];
      if( SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND )
         ++heurdata->nimpsols[ruleindex];

      if( solheur != NULL && strcmp(SCIPheurGetName(solheur), "simplerounding") == 0 )
      {
         ++heurdata->nroundsols[ruleindex];
         if( SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND )
            ++heurdata->nimproundsols[ruleindex];
         if( SCIPgetSolTransObj(scip, sol) < heurdata->bestprimalbds[ruleindex] )
         {
            heurdata->bestprimalbds[ruleindex] = SCIPgetSolTransObj(scip, sol);
            heurdata->bestsolrounded[ruleindex] = TRUE;
         }
      }
      else if( solheur == heur )
      {
         ++heurdata->ndivesols[ruleindex];
         if( SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND )
            ++heurdata->nimpdivesols[ruleindex];
         if( SCIPgetSolTransObj(scip, sol) < heurdata->bestprimalbds[ruleindex] )
         {
            heurdata->bestprimalbds[ruleindex] = SCIPgetSolTransObj(scip, sol);
            heurdata->bestsolrounded[ruleindex] = FALSE;
         }
      }
   }

   return SCIP_OKAY;
}
#endif


/*
 * heuristic specific interface methods
 */

/** creates the masterdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurMasterdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;
#ifdef SCIP_STATISTIC
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   eventhdlr = NULL;
#endif

   /* create Masterdiving primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecMasterdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeMasterdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitMasterdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitMasterdiving) );
#ifdef SCIP_STATISTIC
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolMasterdiving) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolMasterdiving) );

   /* create masterdiving eventhandler data */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );

   /* include event handler */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecMasterdiving, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitMasterdiving) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitMasterdiving) );
#endif

   /* masterdiving heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/"HEUR_NAME"/minreldepth",
         "minimal relative depth to start diving",
         &heurdata->minreldepth, TRUE, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/"HEUR_NAME"/maxreldepth",
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, TRUE, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/"HEUR_NAME"/maxlpiterquot",
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, FALSE, DEFAULT_MAXLPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/"HEUR_NAME"/maxlpiterofs",
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, FALSE, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/"HEUR_NAME"/maxpricerounds",
         "maximal number of allowed pricing rounds (-1: no limit)",
         &heurdata->maxpricerounds, FALSE, DEFAULT_MAXPRICEROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/"HEUR_NAME"/usefarkasonly",
         "perform pricing only if infeasibility is encountered",
         &heurdata->usefarkasonly, FALSE, DEFAULT_USEFARKASONLY, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/"HEUR_NAME"/maxdiveubquot",
         "maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveubquot, TRUE, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/"HEUR_NAME"/maxdiveavgquot",
         "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveavgquot, TRUE, DEFAULT_MAXDIVEAVGQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/"HEUR_NAME"/maxdiveubquotnosol",
         "maximal UBQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveubquotnosol, TRUE, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/"HEUR_NAME"/maxdiveavgquotnosol",
         "maximal AVGQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveavgquotnosol, TRUE, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/"HEUR_NAME"/backtrack",
         "use one level of backtracking if infeasibility is encountered?",
         &heurdata->backtrack, FALSE, DEFAULT_BACKTRACK, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/"HEUR_NAME"/maxdiscrepancy",
         "maximal discrepancy in limited discrepancy search",
         &heurdata->maxdiscrepancy, FALSE, DEFAULT_MAXDISCREPANCY, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/"HEUR_NAME"/maxdiscdepth",
         "maximal depth until which a limited discrepancy search is performed",
         &heurdata->maxdiscdepth, FALSE, DEFAULT_MAXDISCDEPTH, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip,
         "heuristics/"HEUR_NAME"/varselrule",
         "which variable selection should be used? ('c'oefficient, 'f'ractionality, 'l'inesearch, 'p'scost, 'v'eclen; '*': alternate between rules)",
         &heurdata->varselrule, FALSE, DEFAULT_VARSELRULE, ALLOWEDRULES"*", NULL, NULL) );

   return SCIP_OKAY;
}
