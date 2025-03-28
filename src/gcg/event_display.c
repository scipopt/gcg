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
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_EVENTHDLR* eventhdlr;

   assert(gcg != NULL);
   scip = GCGgetMasterprob(gcg);
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
