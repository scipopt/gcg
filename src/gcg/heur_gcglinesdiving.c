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

/**@file   heur_gcglinesdiving.c
 * @brief  LP diving heuristic that fixes variables with a large difference to their root solution
 * @author Tobias Achterberg
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "gcg/heur_gcglinesdiving.h"
#include "gcg/heur_origdiving.h"
#include "gcg/gcg.h"
#include "gcg/pricer_gcg.h"


#define HEUR_NAME             "gcglinesdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings following the line from root solution to current solution"
#define HEUR_DISPCHAR         'l'
#define HEUR_PRIORITY         -1006000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          6
#define HEUR_MAXDEPTH         -1


/*
 * Default diving rule specific parameter settings
 */



/* locally defined diving heuristic data */
struct GCG_DivingData
{
   SCIP_SOL*             rootsol;            /**< relaxation solution at the root node */
   SCIP_Bool             firstrun;           /**< is the heuristic running for the first time? */
};


/*
 * Local methods
 */

/** get relaxation solution of root node (in original variables) */
static
SCIP_RETCODE getRootRelaxSol(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_SOL**            rootsol             /**< pointer to store root relaxation solution */
   )
{
   SCIP* masterprob;
   SCIP_SOL* masterrootsol;
   SCIP_VAR** mastervars;
   int nmastervars;
   int i;

   /* get master problem */
   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);

   /* allocate memory for master root LP solution */
   SCIP_CALL( SCIPcreateSol(masterprob, &masterrootsol, NULL) );

   /* get master variable information */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   /* store root LP values in working master solution */
   for( i = 0; i < nmastervars; i++ )
      SCIP_CALL( SCIPsetSolVal(masterprob, masterrootsol, mastervars[i], SCIPvarGetRootSol(mastervars[i])) );

   /* calculate original root LP solution */
   SCIP_CALL( GCGtransformMastersolToOrigsol(gcg, masterrootsol, rootsol, TRUE, NULL) );

   /* free memory */
   SCIP_CALL( SCIPfreeSol(masterprob, &masterrootsol) );

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** destructor of diving heuristic to free user data (called when GCG is exiting) */
static
GCG_DECL_DIVINGFREE(heurFreeGcglinesdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   GCG_DIVINGDATA* divingdata;

   assert(heur != NULL);

   /* free diving rule specific data */
   divingdata = GCGheurGetDivingDataOrig(heur);
   assert(divingdata != NULL);
   SCIPfreeMemory(GCGgetOrigprob(gcg), &divingdata);
   GCGheurSetDivingDataOrig(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of diving heuristic (called after problem was transformed) */
static
GCG_DECL_DIVINGINIT(heurInitGcglinesdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   GCG_DIVINGDATA* divingdata;

   assert(heur != NULL);

   /* get diving data */
   divingdata = GCGheurGetDivingDataOrig(heur);
   assert(divingdata != NULL);

   /* initialize data */
   divingdata->firstrun = TRUE;
   divingdata->rootsol = NULL;

   return SCIP_OKAY;
}


/** deinitialization method of diving heuristic (called before transformed problem is freed) */
static
GCG_DECL_DIVINGEXIT(heurExitGcglinesdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   GCG_DIVINGDATA* divingdata;
   SCIP* origprob = GCGgetOrigprob(gcg);

   assert(heur != NULL);

   /* get diving data */
   divingdata = GCGheurGetDivingDataOrig(heur);
   assert(divingdata != NULL);

   assert(divingdata->firstrun == TRUE || divingdata->rootsol != NULL);

   /* free root relaxation solution */
   if( divingdata->rootsol != NULL )
      SCIP_CALL( SCIPfreeSol(origprob, &divingdata->rootsol) );

   return SCIP_OKAY;
}


/** execution initialization method of diving heuristic (called when execution of diving heuristic is about to begin) */
static
GCG_DECL_DIVINGINITEXEC(heurInitexecGcglinesdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   GCG_DIVINGDATA* divingdata;

   assert(heur != NULL);

   /* get diving data */
   divingdata = GCGheurGetDivingDataOrig(heur);
   assert(divingdata != NULL);

   /* if the heuristic is running for the first time, the root relaxation solution needs to be stored */
   if( divingdata->firstrun )
   {
      assert(divingdata->rootsol == NULL);
      SCIP_CALL( getRootRelaxSol(gcg, &divingdata->rootsol) );
      assert(divingdata->rootsol != NULL);
      divingdata->firstrun = FALSE;
   }

   return SCIP_OKAY;
}


/** variable selection method of diving heuristic;
 * finds best candidate variable w.r.t. the root LP solution:
 * - in the projected space of fractional variables, extend the line segment connecting the root solution and
 *   the current LP solution up to the point, where one of the fractional variables becomes integral
 * - round this variable to the integral value
 */
static
GCG_DECL_DIVINGSELECTVAR(heurSelectVarGcglinesdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   GCG_DIVINGDATA* divingdata;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   int nlpcands;
   SCIP_Real bestdistquot;
   int c;
   SCIP* origprob = GCGgetOrigprob(gcg);

   /* check preconditions */
   assert(origprob != NULL);
   assert(heur != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);
   assert(bestcandroundup != NULL);

   /* get diving data */
   divingdata = GCGheurGetDivingDataOrig(heur);
   assert(divingdata != NULL);
   assert(divingdata->rootsol != NULL);

   /* get fractional variables that should be integral */
   SCIP_CALL( SCIPgetExternBranchCands(origprob, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL, NULL, NULL) );
   assert(lpcands != NULL);
   assert(lpcandssol != NULL);

   *bestcandmayround = TRUE;
   bestdistquot = SCIPinfinity(origprob);

   /* get best candidate */
   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_VAR* var;
      SCIP_Bool roundup;
      SCIP_Real distquot;
      SCIP_Real solval;
      SCIP_Real rootsolval;

      int i;

      var = lpcands[c];

      /* if the variable is on the tabu list, do not choose it */
       for( i = 0; i < tabulistsize; ++i )
          if( tabulist[i] == var )
             break;
       if( i < tabulistsize )
          continue;

      solval = lpcandssol[c];
      rootsolval = SCIPgetSolVal(origprob, divingdata->rootsol, var);

      /* calculate distance to integral value divided by distance to root solution */
      if( SCIPisLT(origprob, solval, rootsolval) )
      {
         roundup = FALSE;
         distquot = (solval - SCIPfeasFloor(origprob, solval)) / (rootsolval - solval);

         /* avoid roundable candidates */
         if( SCIPvarMayRoundDown(var) )
            distquot *= 1000.0;
      }
      else if( SCIPisGT(origprob, solval, rootsolval) )
      {
         roundup = TRUE;
         distquot = (SCIPfeasCeil(origprob, solval) - solval) / (solval - rootsolval);

         /* avoid roundable candidates */
         if( SCIPvarMayRoundUp(var) )
            distquot *= 1000.0;
      }
      else
      {
         roundup = FALSE;
         distquot = SCIPinfinity(origprob);
      }

      /* check, if candidate is new best candidate */
      if( distquot < bestdistquot )
      {
         *bestcand = var;
         *bestcandmayround = SCIPvarMayRoundDown(var) || SCIPvarMayRoundUp(var);
         *bestcandroundup = roundup;
         bestdistquot = distquot;
      }
   }

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the gcglinesdiving heuristic and includes it in GCG */
SCIP_RETCODE GCGincludeHeurGcglinesdiving(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_HEUR* heur;
   GCG_DIVINGDATA* divingdata;

   /* create gcglinesdiving primal heuristic data */
   SCIP_CALL( SCIPallocMemory(origprob, &divingdata) );

   /* include diving heuristic */
   SCIP_CALL( GCGincludeDivingHeurOrig(gcg, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, heurFreeGcglinesdiving, heurInitGcglinesdiving, heurExitGcglinesdiving, NULL, NULL,
         heurInitexecGcglinesdiving, NULL, heurSelectVarGcglinesdiving, divingdata) );

   assert(heur != NULL);

   return SCIP_OKAY;
}
