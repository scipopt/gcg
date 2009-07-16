/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"
#define SCIP_DEBUG
/**@file   nodesel_master.c
 * @ingroup NODESELECTORS
 * @brief  node selector for coordination of master and original formulation
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "nodesel_master.h"


#define NODESEL_NAME             "master"
#define NODESEL_DESC             "depth first search"
#define NODESEL_STDPRIORITY           0
#define NODESEL_MEMSAVEPRIORITY  100000


/** node selector data */
struct SCIP_NodeselData
{
   SCIP* origscip;
   int lastorignodenumber;
};

/*
 * Callback methods
 */

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
#define nodeselInitsolMaster NULL


/** solving process deinitialization method of node selector (called before branch and bound process data is freed) */
#define nodeselExitsolMaster NULL


/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectMaster)
{  
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODE** nodes;
   int nnodes;
   int i;
   int orignodenumber;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   *selnode = NULL;

   orignodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(nodeseldata->origscip));
   
   if ( orignodenumber != nodeseldata->lastorignodenumber )
   {
      nodeseldata->lastorignodenumber = orignodenumber;

      SCIPdebugMessage("nleaves = %d, nsibling = %d, nchildren = %d\n", SCIPgetNLeaves(scip), SCIPgetNSiblings(scip), SCIPgetNChildren(scip));

      if ( SCIPgetNChildren(scip) > 0 )
      {
         SCIP_CALL( SCIPgetChildren(scip, &nodes, &nnodes) );

         for ( i = 0; i < nnodes; i++ )
         {
            if ( SCIPnodeGetNumber(nodes[i]) == orignodenumber )
            {
               *selnode = nodes[i];
               SCIPdebugMessage("select node (child) with number %d\n", orignodenumber);
               break;
            }
         }

      }

      if ( *selnode == NULL && SCIPgetNSiblings(scip) > 0 )
      {
         SCIP_CALL( SCIPgetSiblings(scip, &nodes, &nnodes) );

         for ( i = 0; i < nnodes; i++ )
         {
            if ( SCIPnodeGetNumber(nodes[i]) == orignodenumber )
            {
               *selnode = nodes[i];
               SCIPdebugMessage("select node (sibling) with number %d\n", orignodenumber);
               break;
            }
         }

      }

      if ( *selnode == NULL && SCIPgetNLeaves(scip) > 0 )
      {
         SCIP_CALL( SCIPgetLeaves(scip, &nodes, &nnodes) );

         for ( i = 0; i < nnodes; i++ )
         {
            if ( SCIPnodeGetNumber(nodes[i]) == orignodenumber )
            {
               *selnode = nodes[i];
               SCIPdebugMessage("select node (leaf) with number %d\n", orignodenumber);
               break;
            }
         }

      }



      if ( *selnode == NULL )
      {
         printf("nodesel_master could not find a node with node number %d!\n", orignodenumber);
      }
      assert(*selnode != NULL);
   }
   else
   {
      SCIPdebugMessage("select random node\n");

      if ( SCIPgetNChildren(scip) > 0 )
      {
         SCIP_CALL( SCIPgetChildren(scip, &nodes, &nnodes) );
         *selnode = nodes[0];
      }
      else if ( SCIPgetNSiblings(scip) > 0 )
      {
         SCIP_CALL( SCIPgetSiblings(scip, &nodes, &nnodes) );
         *selnode = nodes[0];
      }
      else if ( SCIPgetNLeaves(scip) > 0 )
      {
         SCIP_CALL( SCIPgetLeaves(scip, &nodes, &nnodes) );
         *selnode = nodes[0];
      }
   }

   return SCIP_OKAY;
}


/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompMaster)
{  
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   printf("nodeselcomp master!\n");

   if( SCIPnodeGetNumber(node1) < SCIPnodeGetNumber(node2) )
      return 1;
   else 
      return -1;

}





/*
 * master specific interface methods
 */

/** creates the node selector for depth first search and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselMaster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;

   /* create master node selector data */
   SCIP_CALL( SCIPallocMemory(scip, &nodeseldata) );

   /* include node selector */
   SCIP_CALL( SCIPincludeNodesel(scip, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
         nodeselFreeMaster, nodeselInitMaster, nodeselExitMaster, 
         nodeselInitsolMaster, nodeselExitsolMaster, nodeselSelectMaster, nodeselCompMaster,
         nodeseldata) );

   return SCIP_OKAY;
}



void GCGnodeselMasterSetOrigscip(
   SCIP*                 scip,
   SCIP*                 origscip
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODESELDATA* nodeseldata;
   
   assert(scip != NULL);
   assert(origscip != NULL);

   nodesel = SCIPfindNodesel(scip, NODESEL_NAME);
   assert(nodesel != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   nodeseldata->origscip = origscip;
   nodeseldata->lastorignodenumber = -1;
}
