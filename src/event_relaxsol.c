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

/**@file   event_relaxsol.c
 * @brief  eventhandler to update the relaxation solution in the original problem when the master LP has been solved
 * @author Christian Puchert
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "event_relaxsol.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"

#define EVENTHDLR_NAME         "relaxsol"
#define EVENTHDLR_DESC         "eventhandler to update the relaxation solution in the original problem when the master LP has been solved"


/*
 * Callback methods of event handler
 */

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitRelaxsol)
{  /*lint --e{715}*/

   /* notify SCIP that your event handler wants to react on the event type lp solved and solution found */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_LPSOLVED | SCIP_EVENTTYPE_SOLFOUND, eventhdlr, NULL, NULL) );

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

   /* get original problem */
   origprob = GCGmasterGetOrigprob(scip);
   assert(origprob != NULL);

   /* Only transfer the master solution if it is an LP solution or if it is a feasible solution that
    * comes from a master heuristic; otherwise it is assumed to already come from the original problem
    */
   if( (SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND) && SCIPsolGetHeur(SCIPeventGetSol(event)) == NULL )
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   if( SCIPeventGetType(event) & SCIP_EVENTTYPE_LPSOLVED )
   {
      SCIPdebugMessage("Transferring master LP solution to the original problem\n");
   }
   else if( SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND )
   {
      SCIP_SOL* sol = SCIPeventGetSol(event);
      SCIPdebugMessage("Master feasible solution found by <%s> -- transferring to original problem\n",
         SCIPsolGetHeur(sol) == NULL ? "relaxation" : SCIPheurGetName(SCIPsolGetHeur(sol)));
   }
#endif

   SCIP_CALL( GCGrelaxUpdateCurrentSol(origprob) );

   return SCIP_OKAY;
}

/** creates event handler for relaxsol event */
SCIP_RETCODE SCIPincludeEventHdlrRelaxsol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   eventhdlr = NULL;

   /* include event handler into GCG */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecRelaxsol, NULL) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitRelaxsol) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitRelaxsol) );

   return SCIP_OKAY;
}
