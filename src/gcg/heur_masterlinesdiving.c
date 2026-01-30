/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   heur_masterlinesdiving.c
 * @brief  LP diving heuristic that fixes variables with a large difference to their root solution
 * @author Tobias Achterberg
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "gcg/heur_masterlinesdiving.h"
#include "gcg/heur_masterdiving.h"


#define HEUR_NAME             "masterlinesdiving"
#define HEUR_DESC             "master LP diving heuristic that chooses fixings following the line from root solution to current solution"
#define HEUR_DISPCHAR         'l'
#define HEUR_PRIORITY         -1006000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          6
#define HEUR_MAXDEPTH         -1


/*
 * Callback methods
 */

/** variable selection method of diving heuristic;
 * finds best candidate variable w.r.t. the root LP solution:
 * - in the projected space of fractional variables, extend the line segment connecting the root solution and
 *   the current LP solution up to the point, where one of the fractional variables becomes integral
 * - round this variable to the integral value
 */
static
GCG_DECL_MASTER_DIVINGSELECTVAR(heurSelectVarMasterlinesdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   SCIP_Real bestdistquot;
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

   *bestcandmayround = TRUE;
   bestdistquot = SCIPinfinity(masterprob);

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
      for( i = 0; i < tabulistsize; ++i )
         if( tabulist[i] == var )
            break;
      if( i < tabulistsize )
         continue;

      if( SCIPisGT(masterprob, solval, rootsolval) )
      {
         distquot = (SCIPfeasCeil(masterprob, solval) - solval) / (solval - rootsolval);

         /* avoid roundable candidates */
         if( SCIPvarMayRoundUp(var) )
            distquot *= 1000.0;
      }
      else
         distquot = SCIPinfinity(masterprob);

      /* check whether the variable is roundable */
      *bestcandmayround = *bestcandmayround && (SCIPvarMayRoundDown(var) || SCIPvarMayRoundUp(var));

      /* check, if candidate is new best candidate */
      if( distquot < bestdistquot )
      {
         *bestcand = var;
         bestdistquot = distquot;
      }
   }

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the masterlinesdiving heuristic and includes it in GCG */
SCIP_RETCODE GCGincludeHeurMasterlinesdiving(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_HEUR* heur;

   /* include diving heuristic */
   SCIP_CALL( GCGincludeDivingHeurMaster(gcg, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, NULL, NULL, NULL, NULL, NULL, NULL, NULL, heurSelectVarMasterlinesdiving, NULL) );

   assert(heur != NULL);

   return SCIP_OKAY;
}
