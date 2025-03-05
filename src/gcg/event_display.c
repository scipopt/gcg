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

/**@file   event_display.c
 * @brief  eventhdlr to disable the master display after the root node
 * @author Martin Bergner
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/event_display.h"
#include "gcg/gcg.h"
#include <string.h>

#define EVENTHDLR_NAME         "display"
#define EVENTHDLR_DESC         "event handler to disable the master display after the root node"


/*
 * Callback methods of event handler
 */


/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecDisplay)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* solved node does not have to be the root node (can happen when solving was paused and resumed) */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** activates the eventhandler in SCIP */
SCIP_RETCODE GCGactivateEventHdlrDisplay(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(scip != NULL);
   assert(GCGisMaster(scip));

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates event handler for display event */
SCIP_RETCODE GCGincludeEventHdlrDisplay(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_EVENTHDLR* eventhdlr = NULL;
   SCIP* masterprob = GCGgetMasterprob(gcg);

   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));

   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(masterprob, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecDisplay, NULL) );
   assert(eventhdlr != NULL);

   return SCIP_OKAY;
}
