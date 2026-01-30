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

/**@file   event_relaxsol.c
 * @brief  eventhandler to update the relaxation solution in the original problem when the master LP has been solved
 * @author Christian Puchert
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "gcg/event_relaxsol.h"
#include "gcg/relax_gcg.h"
#include "gcg/pricer_gcg.h"
#include "gcg/gcg.h"
#include "gcg/event_mastersol.h"

#define EVENTHDLR_NAME         "relaxsol"
#define EVENTHDLR_DESC         "eventhandler to update the relaxation solution in the original problem when the master LP has been solved"


/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   GCG*                  gcg;                /**< GCG data structure */
   SCIP_Bool             triggered;          /**< flag to indicate whether event has been triggered */
};


/*
 * Callback methods of event handler
 */

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeRelaxsol)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitRelaxsol)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* notify SCIP that your event handler wants to react on the event type lp solved and solution found */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_LPSOLVED | SCIP_EVENTTYPE_SOLFOUND, eventhdlr, NULL, NULL) );
   eventhdlrdata->triggered = FALSE;

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitRelaxsol)
{  /*lint --e{715}*/

   /* notify SCIP that your event handler wants to drop the event type lp solved and solution found */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_LPSOLVED | SCIP_EVENTTYPE_SOLFOUND, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecRelaxsol)
{  /*lint --e{715}*/
   SCIP* origprob;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_Bool violatesvarbnds;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* get original problem */
   origprob = GCGgetOrigprob(eventhdlrdata->gcg);
   assert(origprob != NULL);

   /* Only transfer the master solution if it is an LP solution or if it is a feasible solution that
    * comes from a master heuristic; otherwise it is assumed to already come from the original problem
    */
   if( (SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND) && SCIPsolGetHeur(SCIPeventGetSol(event)) == NULL
      && GCGeventhdlrMastersolIsTriggered(eventhdlrdata->gcg) )
      return SCIP_OKAY;

   eventhdlrdata->triggered = TRUE;

   if( SCIPeventGetType(event) & SCIP_EVENTTYPE_LPSOLVED )
   {
      SCIPdebugMessage("Transferring master LP solution to the original problem\n");
      SCIP_CALL( GCGrelaxUpdateCurrentSol(eventhdlrdata->gcg) );
   }
   else if( SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND )
   {
      SCIP_SOL* sol = SCIPeventGetSol(event);
      SCIP_SOL* origsol;
      SCIP_Bool stored;
      SCIP_Bool foundbyheur = SCIPsolGetHeur(sol) != NULL;

      SCIPdebugMessage("Master feasible solution found by <%s> -- transferring to original problem\n",
         foundbyheur ? SCIPheurGetName(SCIPsolGetHeur(sol)) : "relaxation");

      /* transform the master solution to the original variable space */
      SCIP_CALL( GCGtransformMastersolToOrigsol(eventhdlrdata->gcg, sol, &origsol, foundbyheur, &violatesvarbnds) );
      assert(!violatesvarbnds || !GCGmasterIsSolValid(eventhdlrdata->gcg, sol));

      SCIP_CALL( SCIPtrySolFree(origprob, &origsol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );
      SCIPdebugMessage("  ->%s stored\n", stored ? "" : " not");
   }

   eventhdlrdata->triggered = FALSE;

   return SCIP_OKAY;
}

/** creates event handler for relaxsol event */
SCIP_RETCODE GCGincludeEventHdlrRelaxsol(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP* masterprob = GCGgetMasterprob(gcg);

   eventhdlr = NULL;

   SCIP_CALL( SCIPallocMemory(masterprob, &eventhdlrdata) );
   assert(eventhdlrdata != NULL);
   eventhdlrdata->gcg = gcg;

   /* include event handler into GCG */
   SCIP_CALL( SCIPincludeEventhdlrBasic(masterprob, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecRelaxsol, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrFree(masterprob, eventhdlr, eventFreeRelaxsol) );
   SCIP_CALL( SCIPsetEventhdlrInit(masterprob, eventhdlr, eventInitRelaxsol) );
   SCIP_CALL( SCIPsetEventhdlrExit(masterprob, eventhdlr, eventExitRelaxsol) );

   return SCIP_OKAY;
}

/** return whether event has been triggered */
SCIP_Bool GCGeventhdlrRelaxsolIsTriggered(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(gcg != NULL);

   /* the relaxation solution event handler is not included if BENDERS or ORIGINAL mode is used. As such, it will
    * never be triggered. In this case, it will always return FALSE.
    */
   if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_BENDERS || GCGgetDecompositionMode(gcg) == GCG_DECMODE_ORIGINAL )
      return FALSE;

   eventhdlr = SCIPfindEventhdlr(GCGgetMasterprob(gcg), EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   return eventhdlrdata->triggered;
}
