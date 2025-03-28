/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_mastercoefdiving.c
 * @brief  master LP diving heuristic that chooses fixings w.r.t. the matrix coefficients
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "gcg/heur_mastercoefdiving.h"
#include "gcg/heur_masterdiving.h"


#define HEUR_NAME             "mastercoefdiving"
#define HEUR_DESC             "master LP diving heuristic that chooses fixings w.r.t. the matrix coefficients"
#define HEUR_DISPCHAR         'c'
#define HEUR_PRIORITY         -1001000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          1
#define HEUR_MAXDEPTH         -1


/*
 * Callback methods
 */

/** variable selection method of diving heuristic;
 * finds best candidate variable w.r.t. locking numbers:
 * - prefer variables that may not be rounded without destroying LP feasibility:
 *   - of these variables, round variable with least number of locks in corresponding direction
 * - if all remaining fractional variables may be rounded without destroying LP feasibility:
 *   - round variable with least number of locks in opposite of its feasible rounding direction
 * - binary variables are preferred
 */
static
GCG_DECL_MASTER_DIVINGSELECTVAR(heurSelectVarMastercoefdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   int bestnviolrows;             /* number of violated rows for best candidate */
   SCIP_Real bestcandfrac;        /* fractionality of best candidate */
   int c;
   SCIP* masterprob = GCGgetMasterprob(gcg);

   /* check preconditions */
   assert(masterprob != NULL);
   assert(heur != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);

   /* get fractional variables that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(masterprob, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL, NULL) );
   assert(lpcands != NULL);
   assert(lpcandsfrac != NULL);
   assert(lpcandssol != NULL);

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
      for( i = 0; i < tabulistsize; ++i )
         if( tabulist[i] == var )
            break;
      if( i < tabulistsize )
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
               *bestcand = var;
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
            *bestcand = var;
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


/*
 * heuristic specific interface methods
 */

/** creates the mastercoefdiving heuristic and includes it in GCG */
SCIP_RETCODE GCGincludeHeurMastercoefdiving(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_HEUR* heur;

   /* include diving heuristic */
   SCIP_CALL( GCGincludeDivingHeurMaster(gcg, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, NULL, NULL, NULL, NULL, NULL, NULL, NULL, heurSelectVarMastercoefdiving, NULL) );

   assert(heur != NULL);

   return SCIP_OKAY;
}

