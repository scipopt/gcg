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

/**@file   event_origdiving.c
 * @brief  eventhdlr for origdiving event
 * @author Christian Puchert
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "event_origdiving.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"

#define EVENTHDLR_NAME         "origdiving"
#define EVENTHDLR_DESC         "event handler for original diving solution statistic"

#define DEFAULT_PRINTSTATISTICS FALSE       /**< shall additional statistics about original diving  heuristics be printed? */

#define ALLOWEDRULES           "cfglpv"     /**< possible variable selection rules */


/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_Bool             printstatistics;    /**< shall additional statistics about this heuristic be printed?      */
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
};

/*
 * Callback methods of event handler
 */

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeOrigdiving)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);

   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitOrigdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);

   /* notify GCG that this event should catch the SOLFOUND event */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitOrigdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);

   /* notify GCG that this event should drop the SOLFOUND event */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolOrigdiving)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   const char* rules;
   int nrules;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   /* get event handler data */
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   /* initialize statistical data */
   if( eventhdlrdata->printstatistics )
   {
      int i;

      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->ncalls, nrules) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->nsols, nrules) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->nimpsols, nrules) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->ndivesols, nrules) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->nimpdivesols, nrules) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->nroundsols, nrules) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->nimproundsols, nrules) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->ndives, nrules) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->nrulelpiters, nrules) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->nrulepricerds, nrules) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->bestprimalbds, nrules) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->bestsolrounded, nrules) );

      for( i = 0; i < nrules; ++i )
      {
         eventhdlrdata->ncalls[i] = 0;
         eventhdlrdata->nsols[i] = 0;
         eventhdlrdata->nimpsols[i] = 0;
         eventhdlrdata->ndivesols[i] = 0;
         eventhdlrdata->nimpdivesols[i] = 0;
         eventhdlrdata->nroundsols[i] = 0;
         eventhdlrdata->nimproundsols[i] = 0;
         eventhdlrdata->ndives[i] = 0;
         eventhdlrdata->nrulelpiters[i] = 0;
         eventhdlrdata->nrulepricerds[i] = 0;
         eventhdlrdata->bestprimalbds[i] = SCIPinfinity(scip);
         eventhdlrdata->bestsolrounded[i] = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolOrigdiving)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   const char* rules;
   int nrules;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   /* get event handler data */
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   /* print detailed statistics */
   if( eventhdlrdata->printstatistics )
   {
      int i;

      SCIPinfoMessage(scip, NULL, "Original Diving Heuristics :      Calls       Sols  Improving   DiveSols  Improving  RoundSols  Improving      Dives   LP iters  Price rds    BestPrimal Rounded?\n");
      for( i = 0; i < nrules; ++i )
      {
         SCIPinfoMessage(scip, NULL, "%c                          : %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT,
            rules[i], eventhdlrdata->ncalls[i], eventhdlrdata->nsols[i], eventhdlrdata->nimpsols[i], eventhdlrdata->ndivesols[i], eventhdlrdata->nimpdivesols[i], eventhdlrdata->nroundsols[i], eventhdlrdata->nimproundsols[i], eventhdlrdata->ndives[i], eventhdlrdata->nrulelpiters[i], eventhdlrdata->nrulepricerds[i]);
         if( SCIPisInfinity(scip, eventhdlrdata->bestprimalbds[i]) )
            SCIPinfoMessage(scip, NULL, "      infinity");
         else
            SCIPinfoMessage(scip, NULL, " %13.6e", eventhdlrdata->bestprimalbds[i]);
         SCIPinfoMessage(scip, NULL, eventhdlrdata->bestsolrounded[i] ? "      yes\n" : "       no\n");
      }
      SCIPinfoMessage(scip, NULL, "\n");

      SCIPfreeMemoryArray(scip, &eventhdlrdata->bestsolrounded);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->bestprimalbds);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->nrulepricerds);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->nrulelpiters);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->ndives);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->nimproundsols);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->nroundsols);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->nimpdivesols);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->ndivesols);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->nimpsols);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->nsols);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->ncalls);
   }

   return SCIP_OKAY;
}

/** execution method of event handler; captures the event that simplerounding finds a feasible solution during diving */
static
SCIP_DECL_EVENTEXEC(eventExecOrigdiving)
{  /*lint --e{715}*/
   SCIP* origprob;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_SOL* sol;
   SCIP_HEUR* probingheur;
   SCIP_HEUR* solheur;
   const char* rules;
   char probingchar;
   int nrules;
   int ruleindex;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   /* get event handler data */
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   /* get new primal solution */
   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   /* check if the original problem is currently in probing mode */
   probingheur = GCGrelaxGetProbingheur(origprob);
   probingchar = probingheur != NULL ? SCIPheurGetDispchar(probingheur) : ' ';

   /* this event is irrelevant if we are not in probing */
   if( probingchar == ' ' )
      return SCIP_OKAY;

   /* check if the heuristic is one of the diving heuristics */
   for( ruleindex = 0; ruleindex < nrules; ++ruleindex )
      if( rules[ruleindex] == probingchar )
         break;

   /* this event is irrelevant if probing was not invoked by one of our diving heuristics */
   if( ruleindex == nrules )
      return SCIP_OKAY;

   /* get the heuristic that found the solution */
   solheur = SCIPgetSolHeur(scip, sol);

   /* update the solution statistics */
   if( solheur != NULL && strcmp(SCIPheurGetName(solheur), "simplerounding") == 0 )
   {
      ++eventhdlrdata->nsols[ruleindex];
      ++eventhdlrdata->nroundsols[ruleindex];

      if( SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND )
      {
         ++eventhdlrdata->nimpsols[ruleindex];
         ++eventhdlrdata->nimproundsols[ruleindex];
      }

      if( SCIPgetSolTransObj(scip, sol) < eventhdlrdata->bestprimalbds[ruleindex] )
      {
         eventhdlrdata->bestprimalbds[ruleindex] = SCIPgetSolTransObj(scip, sol);
         eventhdlrdata->bestsolrounded[ruleindex] = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** informs the event handler that a diving heuristic has been called */
SCIP_RETCODE GCGeventOrigdivingCalled(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur                /**< diving heuristic */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   const char* rules;
   int nrules;
   char rule;
   int ruleindex;

   assert(scip != NULL);
   assert(heur != NULL);

   /* get origdiving event handler */
   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /* get event handler data */
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   /* check if the heuristic is one of the diving heuristics */
   rule = SCIPheurGetDispchar(heur);
   for( ruleindex = 0; ruleindex < nrules; ++ruleindex )
      if( rules[ruleindex] == rule )
         break;

   if( ruleindex == nrules )
      return SCIP_OKAY;

   ++eventhdlrdata->ncalls[ruleindex];

   return SCIP_OKAY;
}

/** informs the event handler that a diving heuristic has found a new solution */
SCIP_RETCODE GCGeventOrigdivingNewDivingsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< new solution */
   )
{
   SCIP* origprob;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_HEUR* solheur;
   const char* rules;
   char solchar;
   int nrules;
   int ruleindex;

   assert(scip != NULL);
   assert(sol != NULL);

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   /* get origdiving event handler */
   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /* get event handler data */
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   /* get the heuristic that found the solution */
   solheur = SCIPgetSolHeur(origprob, sol);
   assert(solheur != NULL);
   assert(solheur == GCGrelaxGetProbingheur(origprob));
   solchar = SCIPheurGetDispchar(solheur);

   /* check if the heuristic is one of the diving heuristics */
   for( ruleindex = 0; ruleindex < nrules; ++ruleindex )
      if( rules[ruleindex] == solchar )
         break;

   /* this event is irrelevant if probing was not invoked by one of the diving heuristics */
   if( ruleindex == nrules )
      return SCIP_OKAY;

   ++eventhdlrdata->nsols[ruleindex];
   ++eventhdlrdata->ndivesols[ruleindex];

   /* @todo: I need some better way to check whether the diving solution is improving */
   if( SCIPgetSolTransObj(origprob,sol) == SCIPgetSolTransObj(origprob, SCIPgetBestSol(origprob)) )
   {
      ++eventhdlrdata->nimpsols[ruleindex];
      ++eventhdlrdata->nimpdivesols[ruleindex];
   }

   if( SCIPgetSolTransObj(origprob, sol) < eventhdlrdata->bestprimalbds[ruleindex] )
   {
      eventhdlrdata->bestprimalbds[ruleindex] = SCIPgetSolTransObj(origprob, sol);
      eventhdlrdata->bestsolrounded[ruleindex] = FALSE;
   }

   return SCIP_OKAY;
}

/** updates diving loop statistics of a diving heuristic */
SCIP_RETCODE GCGeventOrigdivingDiveround(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur                /**< diving heuristic */   
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   const char* rules;
   int nrules;
   char rule;
   int ruleindex;

   assert(scip != NULL);
   assert(heur != NULL);

   /* get origdiving event handler */
   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /* get event handler data */
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   /* check if the heuristic is one of the diving heuristics */
   rule = SCIPheurGetDispchar(heur);
   for( ruleindex = 0; ruleindex < nrules; ++ruleindex )
      if( rules[ruleindex] == rule )
         break;

   if( ruleindex == nrules )
      return SCIP_OKAY;

   ++eventhdlrdata->ndives[ruleindex];

   return SCIP_OKAY;
}

/** updates LP statistics of a diving heuristic */
SCIP_RETCODE GCGeventOrigdivingUpdateLPstats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< diving heuristic */   
   SCIP_Longint          nlpiters,           /**< number of new LP iterations */
   int                   npricerounds        /**< number of new pricing rounds */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   const char* rules;
   int nrules;
   char rule;
   int ruleindex;

   assert(scip != NULL);
   assert(heur != NULL);

   /* get origdiving event handler */
   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /* get event handler data */
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* get possible variable selection rules */
   rules = ALLOWEDRULES;
   nrules = strlen(rules);

   /* check if the heuristic is one of the diving heuristics */
   rule = SCIPheurGetDispchar(heur);
   for( ruleindex = 0; ruleindex < nrules; ++ruleindex )
      if( rules[ruleindex] == rule )
         break;

   if( ruleindex == nrules )
      return SCIP_OKAY;

   eventhdlrdata->nrulelpiters[ruleindex] += nlpiters;
   eventhdlrdata->nrulepricerds[ruleindex] += npricerounds;

   return SCIP_OKAY;
}

/** creates event handler for origdiving event */
SCIP_RETCODE SCIPincludeEventHdlrOrigdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP* origprob;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   /* create origdiving event handler data */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );

   eventhdlr = NULL;

   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecOrigdiving, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeOrigdiving) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitOrigdiving) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitOrigdiving) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolOrigdiving) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolOrigdiving) );

   /* add origdiving event handler parameters */
   SCIP_CALL( SCIPaddBoolParam(origprob, "heuristics/"EVENTHDLR_NAME"/printstatistics",
      "shall additional statistics about original diving heuristics be printed?",
      &eventhdlrdata->printstatistics, TRUE, DEFAULT_PRINTSTATISTICS, NULL, NULL) );

   return SCIP_OKAY;
}
