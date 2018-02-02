/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file   event_mastersol.c
 * @brief  eventhdlr to transfer solutions found in the original problem to the master problem
 * @author Christian Puchert
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "event_mastersol.h"
#include "pricer_gcg.h"
#include "gcg.h"
#include "relax_gcg.h"

#define EVENTHDLR_NAME         "mastersol"
#define EVENTHDLR_DESC         "event handler to to transfer solutions found in the original problem to the master problem"


/*
 * Callback methods of event handler
 */

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitMastersol)
{  /*lint --e{715}*/

   /* notify SCIP that your event handler wants to react on the event types best solution found and node solved */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitMastersol)
{  /*lint --e{715}*/

   /* notify SCIP that your event handler wants to drop the event types best solution found and node solved */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecMastersol)
{  /*lint --e{715}*/
   SCIP* masterprob;
   SCIP_SOL* sol;
   SCIP_Bool discretization;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* get new primal solution */
   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   /* get master problem */
   masterprob = GCGgetMasterprob(scip);
   assert(masterprob != NULL);

   /* get discretization parameter */
   SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/gcg/discretization", &discretization) );

   /* transfer solution to the master problem if it was found by a heuristic in the original problem
    * or if discretization is used
    */
   if( SCIPgetStage(scip) > SCIP_STAGE_TRANSFORMED && SCIPgetStage(masterprob) > SCIP_STAGE_TRANSFORMED &&
      (SCIPsolGetHeur(sol) != NULL || (discretization && SCIPgetStage(masterprob) != SCIP_STAGE_SOLVED)) )
   {
      SCIPdebugMessage("Original feasible solution found by <%s> -- transferring to master problem\n",
         SCIPsolGetHeur(sol) == NULL ? "relaxation" : SCIPheurGetName(SCIPsolGetHeur(sol)));
      SCIP_CALL( GCGmasterTransOrigSolToMasterVars(masterprob, sol, NULL) );
   }

   return SCIP_OKAY;
}

/** creates event handler for mastersol event */
SCIP_RETCODE SCIPincludeEventHdlrMastersol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   eventhdlr = NULL;

   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecMastersol, NULL) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitMastersol) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitMastersol) );

   return SCIP_OKAY;
}
