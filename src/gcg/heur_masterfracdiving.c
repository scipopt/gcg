/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
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

/**@file   heur_masterfracdiving.c
 * @brief  master LP diving heuristic that chooses fixings w.r.t. the fractionalities
 * @author Tobias Achterberg
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "gcg/heur_masterfracdiving.h"
#include "gcg/heur_masterdiving.h"


#define HEUR_NAME             "masterfracdiving"
#define HEUR_DESC             "master LP diving heuristic that chooses fixings w.r.t. the fractionalities"
#define HEUR_DISPCHAR         'f'
#define HEUR_PRIORITY         -1003000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          3
#define HEUR_MAXDEPTH         -1


/*
 * Callback methods
 */

/** variable selection method of diving heuristic;
 * finds best candidate variable w.r.t. fractionality:
 * - prefer variables that may not be rounded without destroying LP feasibility:
 *   - of these variables, round least fractional variable in corresponding direction
 * - if all remaining fractional variables may be rounded without destroying LP feasibility:
 *   - round variable with least increasing objective value
 * - binary variables are preferred
 */
static
GCG_DECL_MASTER_DIVINGSELECTVAR(heurSelectVarMasterfracdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   SCIP_Real bestobjgain;
   SCIP_Real bestfrac;
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
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
   bestobjgain = SCIPinfinity(masterprob);
   bestfrac = SCIP_INVALID;

   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_VAR* var;
      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;
      SCIP_Real frac;
      SCIP_Real obj;

      int i;

      var = lpcands[c];

      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      frac = lpcandsfrac[c];
      obj = SCIPvarGetObj(var);

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
            SCIP_Real objgain;

            objgain = (1.0-frac)*obj;

            /* penalize too small fractions */
            if( ABS(1.0 - frac) < 0.01 )
               objgain *= 1000.0;

            /* prefer decisions on binary variables */
            if( !SCIPvarIsBinary(var) )
               objgain *= 1000.0;

            /* check, if candidate is new best candidate */
            if( SCIPisLT(masterprob, objgain, bestobjgain) || (SCIPisEQ(masterprob, objgain, bestobjgain) && frac > bestfrac) )
            {
               *bestcand = var;
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
            *bestcand = var;
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


/*
 * heuristic specific interface methods
 */

/** creates the masterfracdiving heuristic and includes it in GCG */
SCIP_RETCODE GCGincludeHeurMasterfracdiving(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_HEUR* heur;

   /* include diving heuristic */
   SCIP_CALL( GCGincludeDivingHeurMaster(gcg, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, NULL, NULL, NULL, NULL, NULL, NULL, NULL, heurSelectVarMasterfracdiving, NULL) );

   assert(heur != NULL);

   return SCIP_OKAY;
}

