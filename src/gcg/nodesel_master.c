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

/**@file   nodesel_master.c
 * 
 * @brief  node selector for coordination of master and original formulation
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "gcg/nodesel_master.h"
#include "gcg/cons_origbranch.h"
#include "gcg/cons_masterbranch.h"
#include "gcg/pricer_gcg.h"

#define NODESEL_NAME             "master"
#define NODESEL_DESC             "orig master coordination"
#define NODESEL_STDPRIORITY           0
#define NODESEL_MEMSAVEPRIORITY  100000
#define EVENTHDLR_NAME           "masterfocusnode"
#define EVENTHDLR_DESC           "event handler to transfer dual bounds from orig to master"


/** node selector data */
struct SCIP_NodeselData
{
   GCG*                 gcg;                /**< GCG data structure */
   SCIP_Longint         lastorignodenumber;
};

struct SCIP_EventhdlrData
{
   GCG*                  gcg;                /**< GCG data structure */
};

/*
 * Callback methods
 */

#define nodeselCopyMaster NULL

/** destructor of node selector to free user data (called when SCIP is exiting) */
static
SCIP_DECL_NODESELFREE(nodeselFreeMaster)
{
   SCIP_NODESELDATA* nodeseldata;

   assert(nodesel != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   SCIPfreeMemory(scip, &nodeseldata);

   return SCIP_OKAY;
}


/** initialization method of node selector (called after problem was transformed) */
#define nodeselInitMaster NULL


/** deinitialization method of node selector (called before transformed problem is freed) */
#define nodeselExitMaster NULL


/** solving process initialization method of node selector (called when branch and bound process is about to begin) */
static
SCIP_DECL_NODESELINITSOL(nodeselInitsolMaster)
{
   SCIP_NODESELDATA* nodeseldata;
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   nodeseldata->lastorignodenumber = -1LL;

   return SCIP_OKAY;
}


/** solving process deinitialization method of node selector (called before branch and bound process data is freed) */
#define nodeselExitsolMaster NULL


/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectMaster)
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODE** nodes;
   SCIP* origscip;
   int nnodes;
   SCIP_NODE* orignode;
   SCIP_Longint orignodenumber;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   origscip = GCGgetOrigprob(nodeseldata->gcg);

   *selnode = NULL;

   orignode = SCIPgetCurrentNode(origscip);
   orignodenumber = SCIPnodeGetNumber(orignode);

   if( orignodenumber != nodeseldata->lastorignodenumber )
   {
      SCIP_CONS* origcons = GCGconsOrigbranchGetActiveCons(nodeseldata->gcg);
      SCIP_CONS* parentorigcons = GCGconsOrigbranchGetParentcons(origcons);

      nodeseldata->lastorignodenumber = orignodenumber;

      /* check whether the current node is the root node and has no parent */
      if( parentorigcons == NULL )
      {
         if( SCIPgetRootNode(origscip) == orignode )
         {
            assert((GCGconsOrigbranchGetNode(origcons) == SCIPgetRootNode(origscip)) || (GCGconsOrigbranchGetNode(origcons) == NULL) );
            assert(GCGconsOrigbranchGetMastercons(origcons) != NULL);
            assert((GCGconsMasterbranchGetNode(GCGconsOrigbranchGetMastercons(origcons)) == SCIPgetRootNode(scip)) || (GCGconsMasterbranchGetNode(GCGconsOrigbranchGetMastercons(origcons)) == NULL));

            *selnode = SCIPgetRootNode(scip);
            SCIPdebugMessage("selected root node in the master program\n");
            SCIP_CALL( GCGrestoreLimitSettings(nodeseldata->gcg) );
         }
         else
         {
            /* solving was interrupted and SCIP created a dummy node in the original problem */
            assert(SCIPgetNChildren(scip) < 2);
            if( SCIPgetNChildren(scip) == 0 )
            {
               /* create and select dummy node */
               SCIP_NODE* child;
               SCIPcreateChild(scip, &child, 0.0, SCIPgetLocalTransEstimate(scip));
               *selnode = child;
               SCIPdebugMessage("created and selected dummy node in the master program\n");
            }
            else
            {
               /* select dummy node */
               SCIP_NODE** children;
               int nchildren;
               SCIPgetChildren(scip, &children, &nchildren);
               assert(nchildren > 0);
               *selnode = children[0];
               SCIPdebugMessage("selected dummy node in the master program\n");
            }
         }
      }
      else
      {
         assert(GCGconsOrigbranchGetMastercons(parentorigcons) != NULL);
         *selnode = GCGconsMasterbranchGetNode(GCGconsOrigbranchGetMastercons(origcons));

         assert(SCIPnodeGetDepth(GCGconsMasterbranchGetNode(GCGconsOrigbranchGetMastercons(parentorigcons))) == SCIPnodeGetDepth(GCGconsOrigbranchGetNode(parentorigcons)));
         assert( *selnode != NULL );
      }

      if( *selnode == NULL )
      {
         SCIPerrorMessage("nodesel_master could not find a node corresponding to the current original node!\n");
      }
      assert(*selnode != NULL);
   }
   else
   {
      SCIPdebugMessage("select random node\n");

      if( SCIPgetNChildren(scip) > 0 )
      {
         SCIP_CALL( SCIPgetChildren(scip, &nodes, &nnodes) );
         *selnode = nodes[0];
      }
      else if( SCIPgetNSiblings(scip) > 0 )
      {
         SCIP_CALL( SCIPgetSiblings(scip, &nodes, &nnodes) );
         *selnode = nodes[0];
      }
      else if( SCIPgetNLeaves(scip) > 0 )
      {
         SCIP_CALL( SCIPgetLeaves(scip, &nodes, &nnodes) );
         *selnode = nodes[0];
      }
   }

#ifndef NDEBUG
   GCGconsOrigbranchCheckConsistency(nodeseldata->gcg);
   GCGconsMasterbranchCheckConsistency(nodeseldata->gcg);
#endif

   return SCIP_OKAY;
}


/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompMaster)
{
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   if( SCIPnodeGetNumber(node1) < SCIPnodeGetNumber(node2) )
      return 1;
   else
      return -1;
}


/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitFocusnode)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);

   /* catch SCIP_EVENTTYPE_NODEFOCUSED events */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}


/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitFocusnode)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);

   /* stop event handling */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** destructor of event handler */
static
SCIP_DECL_EVENTFREE(eventFreeFocusnode)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}


/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecFocusnode)
{  /*lint --e{715}*/
   SCIP* origscip;
   SCIP_NODE* focusnode;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   origscip = GCGgetOrigprob(eventhdlrdata->gcg);
   assert(origscip != NULL);
   focusnode = SCIPeventGetNode(event);
   assert(focusnode != NULL);

   /* set the dual bound to the lower bound of the corresponding original node */
   SCIP_CALL( SCIPupdateNodeDualbound(scip, focusnode, SCIPgetNodeLowerbound(origscip, SCIPgetCurrentNode(origscip))) );
   assert(
      (GCGconsOrigbranchGetNode(GCGconsOrigbranchGetActiveCons(eventhdlrdata->gcg)) == SCIPgetRootNode(origscip)
         && SCIPnodeGetDepth(SCIPgetCurrentNode(origscip)) == SCIPnodeGetDepth(focusnode))
      || focusnode == GCGconsMasterbranchGetNode(GCGconsOrigbranchGetMastercons(GCGconsOrigbranchGetActiveCons(eventhdlrdata->gcg))));

   return SCIP_OKAY;
}


/*
 * master specific interface methods
 */

/** creates the node selector for depth first search and includes it in SCIP */
SCIP_RETCODE GCGincludeNodeselMaster(
   GCG*                  gcg                 /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP* scip = GCGgetMasterprob(gcg);

   /* create master node selector data */
   SCIP_CALL( SCIPallocMemory(scip, &nodeseldata) );
   nodeseldata->gcg =  gcg;
   nodeseldata->lastorignodenumber = -1LL;

   /* include node selector */
   SCIP_CALL( SCIPincludeNodesel(scip, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
      nodeselCopyMaster, nodeselFreeMaster, nodeselInitMaster, nodeselExitMaster,
      nodeselInitsolMaster, nodeselExitsolMaster, nodeselSelectMaster, nodeselCompMaster,
      nodeseldata) );

   /* include event handler into SCIP */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   eventhdlrdata->gcg =  gcg;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecFocusnode, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitFocusnode) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitFocusnode) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeFocusnode) );

   return SCIP_OKAY;
}
