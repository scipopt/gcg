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

/**@file   gcgheur.c
 * @brief  public methods for GCG heuristics
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/gcg.h"
#include "gcg/pub_gcgheur.h"


/** resets the parameters to their default value */
static
SCIP_RETCODE setOrigHeuristicsDefault(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   /* set specific parameters for LNS heuristics */
   SCIP_CALL( SCIPresetParam(scip, "heuristics/gcgrens/nodesofs") );
   SCIP_CALL( SCIPresetParam(scip, "heuristics/gcgrens/minfixingrate") );
   SCIP_CALL( SCIPresetParam(scip, "heuristics/gcgrins/nodesofs") );
   SCIP_CALL( SCIPresetParam(scip, "heuristics/gcgrins/minfixingrate") );
   SCIP_CALL( SCIPresetParam(scip, "heuristics/xpcrossover/nodesofs") );
   SCIP_CALL( SCIPresetParam(scip, "heuristics/xpcrossover/minfixingrate") );
   SCIP_CALL( SCIPresetParam(scip, "heuristics/xprins/nodesofs") );
   SCIP_CALL( SCIPresetParam(scip, "heuristics/xprins/minfixingrate") );

   return SCIP_OKAY;
}

/** sets the parameters to aggressive values */
static
SCIP_RETCODE setOrigHeuristicsAggressive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   /* set specific parameters for GCG RENS heuristic, if the heuristic is included */
#ifndef NDEBUG
   if( SCIPfindHeur(scip, "gcgrens") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetLongintParam(scip, "heuristics/gcgrens/nodesofs", (SCIP_Longint)2000) );
      SCIP_CALL( SCIPsetRealParam(scip, "heuristics/gcgrens/minfixingrate", 0.3) );
   }

   /* set specific parameters for GCG RINS heuristic, if the heuristic is included */
#ifndef NDEBUG
   if( SCIPfindHeur(scip, "gcgrins") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetIntParam(scip, "heuristics/gcgrins/nodesofs", 2000) );
      SCIP_CALL( SCIPsetRealParam(scip, "heuristics/gcgrins/minfixingrate", 0.3) );
   }

   /* set specific parameters for XP Crossover heuristic, if the heuristic is included */
#ifndef NDEBUG
   if( SCIPfindHeur(scip, "xpcrossover") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetLongintParam(scip, "heuristics/xpcrossover/nodesofs", (SCIP_Longint)2000) );
      SCIP_CALL( SCIPsetRealParam(scip, "heuristics/xpcrossover/minfixingrate", 0.3) );
   }

   /* set specific parameters for XP RINS heuristic, if the heuristic is included */
#ifndef NDEBUG
   if( SCIPfindHeur(scip, "xprins") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetLongintParam(scip, "heuristics/xprins/nodesofs", (SCIP_Longint)2000) );
      SCIP_CALL( SCIPsetRealParam(scip, "heuristics/xprins/minfixingrate", 0.3) );
   }

   return SCIP_OKAY;
}

/** sets the parameters to fast values */
static
SCIP_RETCODE setOrigHeuristicsFast(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   int i;

#define NEXPENSIVEHEURS 11
   static const char* const expensiveheurs[NEXPENSIVEHEURS] = {
      "gcgcoefdiving",
      "gcgfeaspump",
      "gcgfracdiving",
      "gcgguideddiving",
      "gcglinesdiving",
      "gcgpscostdiving",
      "gcgrens",
      "gcgrins",
      "gcgveclendiving",
      "xpcrossover",
      "xprins"
   };

   assert(scip != NULL);

   SCIP_CALL( setOrigHeuristicsDefault(scip) );

   /* explicitly turn off expensive heuristics, if included */
   for( i = 0; i < NEXPENSIVEHEURS; ++i )
   {
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/freq", expensiveheurs[i]);
      SCIP_CALL( SCIPsetIntParam(scip, paramname, -1) );
   }

   return SCIP_OKAY;
}

/** sets heuristic parameters values to
 *
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all heuristic parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for heuristic is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the heuristic are called more aggregative
 *  - SCIP_PARAMSETTING_OFF which turns off all heuristics
 */
SCIP_RETCODE GCGsetHeuristics(
   GCG*                  gcg,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting        /**< parameter settings */
   )
{
   SCIP* origprob;
   assert(paramsetting == SCIP_PARAMSETTING_DEFAULT || paramsetting == SCIP_PARAMSETTING_FAST
      || paramsetting == SCIP_PARAMSETTING_AGGRESSIVE || paramsetting == SCIP_PARAMSETTING_OFF);

   origprob = GCGgetOrigprob(gcg);

   switch( paramsetting )
   {
   case SCIP_PARAMSETTING_AGGRESSIVE:
      SCIP_CALL( setOrigHeuristicsAggressive(origprob) );
      break;
   case SCIP_PARAMSETTING_OFF:
      break;
   case SCIP_PARAMSETTING_FAST:
      SCIP_CALL( setOrigHeuristicsFast(origprob) );
      break;
   case SCIP_PARAMSETTING_DEFAULT:
      SCIP_CALL( setOrigHeuristicsDefault(origprob) );
      break;
   default:
      SCIPerrorMessage("The given paramsetting is invalid!\n");
      break;
   }

   return SCIP_OKAY;
}
