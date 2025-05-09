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

/**@file   cons_origbranch.c
 * 
 * @brief  constraint handler for storing the branching decisions at each node of the tree
 * @author Gerald Gamrath
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */
/* #define CHECKCONSISTENCY */
#include <assert.h>
#include <string.h>

#include "gcg/gcg.h"
#include "gcg/branch_generic.h"
#include "gcg/cons_masterbranch.h"
#include "gcg/cons_origbranch.h"

#include "scip/nodesel_bfs.h"
#include "scip/nodesel_dfs.h"
#include "scip/nodesel_estimate.h"
#include "scip/nodesel_hybridestim.h"
#include "scip/nodesel_restartdfs.h"
#include "scip/branch_allfullstrong.h"
#include "scip/branch_fullstrong.h"
#include "scip/branch_inference.h"
#include "scip/branch_mostinf.h"
#include "scip/branch_leastinf.h"
#include "scip/branch_pscost.h"
#include "scip/branch_random.h"
#include "scip/branch_relpscost.h"
#include "scip/cons_linear.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "origbranch"
#define CONSHDLR_DESC          "store branching decision at nodes of the tree constraint handler"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              * propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */


/** constraint data for branch orig constraints */
struct SCIP_ConsData
{
   SCIP_NODE*            node;               /**< the node at which the cons is sticking */
   SCIP_CONS*            parentcons;         /**< the origbranch constraint of the parent node */
   SCIP_CONS**           childconss;         /**< array of the masterbranch constraints of child nodes */
   int                   nchildconss;        /**< number of the masterbranch constraints of child nodes */
   int                   maxchildconss;      /**< capacity of childconss */
   SCIP_CONS*            probingtmpcons;     /**< pointer to save the last child in the childconss array if it is overwritten in probing mode */
   SCIP_CONS*            mastercons;         /**< the masterbranch constraint of the corresponding node
                                              *   in the master program */
   GCG_BRANCHDATA*       branchdata;         /**< branching data stored by the branching rule containing information
                                              *   about the branching restrictions */
   SCIP_BRANCHRULE*      branchrule;         /**< branching rule that created the corresponding node and imposed
                                              *   branching restrictions */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   GCG*                  gcg;                /**< GCG data structure */
   SCIP_CONS**           stack;              /**< stack for storing active constraints */
   int                   nstack;             /**< number of elements on the stack */
   int                   maxstacksize;       /**< maximum size of the stack */
   SCIP_CONS*            rootcons;           /**< constraint in the root node */
};


/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   SCIPdebugMessage("freeing branch orig constraint handler (nconss = %d)\n", SCIPconshdlrGetNConss(conshdlr));

   /* free constraint handler storage */
   assert(conshdlrData->stack == NULL);
   if( conshdlrData->rootcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conshdlrData->rootcons) );
   }

   SCIPfreeMemory(scip, &conshdlrData);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   /* prepare stack */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrData->stack, conshdlrData->maxstacksize) );
   assert(conshdlrData->nstack >= 0);

   /* check consistency */
   if( conshdlrData->rootcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conshdlrData->rootcons) );
      conshdlrData->rootcons = NULL;
      --(conshdlrData->nstack);
   }

   GCGconsOrigbranchCheckConsistency(conshdlrData->gcg);

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   GCG_BRANCHDATA* branchdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->nstack <= 1);

   SCIPdebugMessage("exiting solution process branch orig constraint handler (nconss = %d)\n", SCIPconshdlrGetNConss(conshdlr));

   /* check for root */
   if( conshdlrdata->rootcons != NULL )
   {
      branchdata = GCGconsOrigbranchGetBranchdata(conshdlrdata->rootcons);

      SCIPfreeMemoryNull(scip, &branchdata );

      SCIP_CALL( SCIPreleaseCons(scip, &conshdlrdata->rootcons) );
      conshdlrdata->rootcons = NULL;
   }

   /* free stack */
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->stack, conshdlrdata->maxstacksize);
   conshdlrdata->stack = NULL;

   return SCIP_OKAY;
}

/** exit method of constraint handler (called before problem is free transformed) */
static
SCIP_DECL_CONSEXIT(consExitOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(conshdlr != NULL);
   assert(scip != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMessage("exiting transformed branch orig constraint handler (nconss = %d)\n", SCIPconshdlrGetNConss(conshdlr));

   if( conshdlrdata->rootcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conshdlrdata->rootcons) );
      conshdlrdata->rootcons = NULL;
   }
   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSDATA* parentdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(*consdata != NULL);

   SCIPdebugMessage("Deleting branch orig constraint: <%s>.\n", SCIPconsGetName(cons));

   /* set the origcons pointer of the corresponding mastercons to NULL */
   if( (*consdata)->mastercons != NULL )
   {
      assert(GCGconsMasterbranchGetOrigcons((*consdata)->mastercons) == cons);
      GCGconsMasterbranchSetOrigcons((*consdata)->mastercons, NULL);
   }

   /* set the pointer in the parent constraint to NULL */
   if( (*consdata)->parentcons != NULL )
   {
#ifndef NDEBUG
      SCIP_Bool childdeleted = FALSE;
#endif
      parentdata = SCIPconsGetData((*consdata)->parentcons);

      if( SCIPinProbing(scip) )
         parentdata->probingtmpcons = NULL;

      for( i = 0; i < parentdata->nchildconss; ++i )
      {
         if( parentdata->childconss[i] == cons )
         {

            parentdata->childconss[i] = parentdata->childconss[parentdata->nchildconss-1];
            parentdata->childconss[parentdata->nchildconss-1] = NULL;
#ifndef NDEBUG
            childdeleted = TRUE;
#endif
            (parentdata->nchildconss) -= 1;
            break;
         }
      }
      assert(childdeleted || SCIPinProbing(scip));
   }

   /* no child nodes may exist */
   for( i = 0; i < (*consdata)->nchildconss; ++i )
      assert((*consdata)->childconss[i] == NULL);

   /* allow the correspondig branchrule to delete the branch data */
   if( (*consdata)->branchdata != NULL && (*consdata)->branchrule != NULL )
   {
      SCIP_Bool force = ((*consdata)->mastercons == NULL);
      SCIP_CALL( GCGrelaxBranchDataDelete(GCGorigGetGcg(scip), (*consdata)->branchrule, &(*consdata)->branchdata, TRUE, force) );
      if( (*consdata)->mastercons != NULL && (*consdata)->branchdata == NULL )
         GCGconsMasterbranchSetBranchdata((*consdata)->mastercons, NULL);
   }

   (*consdata)->branchdata = NULL;

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->childconss, (*consdata)->maxchildconss);
   (*consdata)->childconss = NULL;
   (*consdata)->nchildconss = 0;
   (*consdata)->maxchildconss = 0;

   /* free constraint data */
   SCIPfreeBlockMemory(scip, consdata);
   *consdata = NULL;

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   assert(SCIPconsGetData(cons) != NULL);

   if( SCIPconsGetData(cons)->node == NULL )
      SCIPconsGetData(cons)->node = SCIPgetRootNode(scip);

   SCIPdebugMessage("Activating branch orig constraint: <%s>[stack size: %d].\n", SCIPconsGetName(cons),
      conshdlrData->nstack+1);

   /* put constraint on the stack */
   if( conshdlrData->nstack >= conshdlrData->maxstacksize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, conshdlrData->nstack+1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(conshdlrData->stack), conshdlrData->maxstacksize, newsize) );
      SCIPdebugMessage("reallocating Memory for stack! %d --> %d\n", conshdlrData->maxstacksize, newsize);
      conshdlrData->maxstacksize = newsize;
   }

   /* put constraint on the stack */
   assert(conshdlrData->stack != NULL);
   conshdlrData->stack[conshdlrData->nstack] = cons;
   ++(conshdlrData->nstack);

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL || conshdlrData->nstack <= 1);
   assert(conshdlrData->nstack <= 1 || cons == conshdlrData->stack[conshdlrData->nstack-1]);

   assert(SCIPconsGetData(cons) != NULL);

   SCIPdebugMessage("Deactivating branch orig constraint: <%s> [stack size: %d].\n",
      SCIPconsGetName(cons), conshdlrData->nstack-1);

   /* remove constraint from the stack */
   if( conshdlrData->nstack > 0 )
      --(conshdlrData->nstack);

   return SCIP_OKAY;
}

/** lp solution enforcement method */
static
SCIP_DECL_CONSENFOLP(consEnfolpOrigbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** lp solution enforcement method */
static
SCIP_DECL_CONSENFORELAX(consEnforeOrigbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** pseudo solution enforcement method */
static
SCIP_DECL_CONSENFOPS(consEnfopsOrigbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** solution check method */
static
SCIP_DECL_CONSCHECK(consCheckOrigbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** variable lock method */
static
SCIP_DECL_CONSLOCK(consLockOrigbranch)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/*
 * interface methods
 */


/** creates the handler for origbranch constraints and includes it in SCIP */
SCIP_RETCODE GCGincludeConshdlrOrigbranch(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSHDLR* conshdlr;
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   SCIP_CALL( SCIPallocMemory(origprob, &conshdlrData) );
   conshdlrData->gcg = gcg;
   conshdlrData->stack = NULL;
   conshdlrData->nstack = 0;
   conshdlrData->maxstacksize = 25;
   conshdlrData->rootcons = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(origprob, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpOrigbranch, consEnfopsOrigbranch, consCheckOrigbranch,
         consLockOrigbranch, conshdlrData) );
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPsetConshdlrEnforelax(origprob, conshdlr, consEnforeOrigbranch) );
   SCIP_CALL( SCIPsetConshdlrFree(origprob, conshdlr, consFreeOrigbranch) );
   SCIP_CALL( SCIPsetConshdlrExit(origprob, conshdlr, consExitOrigbranch) );
   SCIP_CALL( SCIPsetConshdlrInitsol(origprob, conshdlr, consInitsolOrigbranch) );
   SCIP_CALL( SCIPsetConshdlrExitsol(origprob, conshdlr, consExitsolOrigbranch) );
   SCIP_CALL( SCIPsetConshdlrDelete(origprob, conshdlr, consDeleteOrigbranch) );
   SCIP_CALL( SCIPsetConshdlrActive(origprob, conshdlr, consActiveOrigbranch) );
   SCIP_CALL( SCIPsetConshdlrDeactive(origprob, conshdlr, consDeactiveOrigbranch) );

   return SCIP_OKAY;
}


/** creates and captures a origbranch constraint */
SCIP_RETCODE GCGcreateConsOrigbranch(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_NODE*            node,               /**< the node to which this origbranch constraint belongs */
   SCIP_CONS*            parentcons,         /**< origbranch constraint associated with the father node */
   SCIP_BRANCHRULE*      branchrule,         /**< the branching rule that created the b&b node the constraint belongs to */
   GCG_BRANCHDATA*       branchdata          /**< branching data storing information about the branching restrictions at the
                                              *   corresponding node */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(gcg != NULL);
   assert((parentcons == NULL) == (node == NULL));

   scip = GCGgetOrigprob(gcg);
   /* find the origbranch constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   /* initialize the fields in the constraint data */
   consdata->node = node;
   consdata->parentcons = parentcons;
   consdata->childconss = NULL;
   consdata->nchildconss = 0;
   consdata->maxchildconss = 0;
   consdata->probingtmpcons = NULL;
   consdata->mastercons = NULL;
   consdata->branchdata = branchdata;
   consdata->branchrule = branchrule;

   SCIPdebugMessage("Creating branch orig constraint: <%s> (nconss = %d).\n", name, SCIPconshdlrGetNConss(conshdlr));

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         TRUE, FALSE, FALSE, FALSE, TRUE) );

   /* store the pointer to the new constraint in the origbranch constraint of the parent node */
   if( parentcons != NULL )
   {
      SCIP_CONSDATA* parentdata;

      parentdata = SCIPconsGetData(parentcons);
      assert(parentdata != NULL);

      if( SCIPinProbing(scip) )
      {
         parentdata->probingtmpcons = *cons;
      }
      else
      {
         int newsize;
         ++parentdata->nchildconss;
         newsize = SCIPcalcMemGrowSize(scip, parentdata->nchildconss);
         if( parentdata->nchildconss == 1 )
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(parentdata->childconss), newsize) );
            parentdata->childconss[0] = NULL;
         }
         else
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(parentdata->childconss), parentdata->maxchildconss, newsize) );
            parentdata->childconss[parentdata->nchildconss - 1] = NULL;
         }
         parentdata->maxchildconss = newsize;
         parentdata->childconss[parentdata->nchildconss-1] = *cons;
      }
   }

   return SCIP_OKAY;
}

/** returns the branch orig constraint of the current node, only needs the pointer to scip */
SCIP_CONS* GCGconsOrigbranchGetActiveCons(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(gcg != NULL);
   scip = GCGgetOrigprob(gcg);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);
   assert(conshdlrData->nstack > 0);

   return conshdlrData->stack[conshdlrData->nstack-1];
}

/** returns the stack and the number of elements on it */
void GCGconsOrigbranchGetStack(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(gcg != NULL);
   scip = GCGgetOrigprob(gcg);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   *stack = conshdlrData->stack;
   *nstackelements = conshdlrData->nstack;
}

/** sets the branching data for a given origbranch constraint */
void GCGconsOrigbranchSetBranchdata(
   SCIP_CONS*            cons,               /**< origbranch constraint for which the branching data is requested */
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->branchdata = branchdata;
}

/** returns the branching data for a given origbranch constraint */
GCG_BRANCHDATA* GCGconsOrigbranchGetBranchdata(
   SCIP_CONS*            cons                /**< origbranch constraint for which the branching data is requested */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->branchdata;
}

/** returns the branching rule for a given origbranch constraint */
SCIP_BRANCHRULE* GCGconsOrigbranchGetBranchrule(
   SCIP_CONS*            cons                /**< origbranch constraint for which the branchrule is requested */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->branchrule;
}

/** returns the node in the B&B tree at which the given origbranch constraint is sticking */
SCIP_NODE* GCGconsOrigbranchGetNode(
   SCIP_CONS*            cons                /**< origbranch constraint for which the corresponding node is requested */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->node;
}

/** returns the origbranch constraint of the B&B father of the node at which the
  * given origbranch constraint is sticking
  */
SCIP_CONS* GCGconsOrigbranchGetParentcons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the origbranch constraint of
                                              *   the father node is requested */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->parentcons;
}

/** returns the number of origbranch constraints of the children of the node at which the
  * given origbranch constraint is sticking
  */
int GCGconsOrigbranchGetNChildconss(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nchildconss;
}

/** returns an origbranch constraint of a child of the node at which the
  * given origbranch constraint is sticking
  */
SCIP_CONS* GCGconsOrigbranchGetChildcons(
   SCIP_CONS*            cons,               /**< constraint */
   int                   childnr             /**< number of child */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->childconss != NULL);
   assert(consdata->nchildconss > childnr);

   return consdata->childconss[childnr];
}

/** sets the masterbranch constraint of the node in the master program corresponding to the node
  * at which the given origbranchbranch constraint is sticking
  */
void GCGconsOrigbranchSetMastercons(
   SCIP_CONS*            cons,               /**< origbranch constraint for which the masterbranch constraint should be set */
   SCIP_CONS*            mastercons          /**< masterbranch constraint corresponding to the given origbranch constraint */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->mastercons = mastercons;
}

/** returns the masterbranch constraint of the node in the master program corresponding to the node
  * at which the given origbranchbranch constraint is sticking
  */
SCIP_CONS* GCGconsOrigbranchGetMastercons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the corresponding masterbranch
                                              *   constraint is requested */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->mastercons;
}

/** adds initial constraint to root node */
SCIP_RETCODE GCGconsOrigbranchAddRootCons(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONS** conss;
   int nconss;
   int i;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNConss(conshdlr);
   assert(nconss <= 1);
   conss = SCIPconshdlrGetConss(conshdlr);
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( SCIPdelCons(scip, conss[i]) );
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(SCIPconshdlrGetNConss(conshdlr) == 0);
   if( conshdlrdata->rootcons == NULL )
   {
      SCIP_CALL( GCGcreateConsOrigbranch(gcg, &cons, "root-origbranch", NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddConsNode(scip, SCIPgetRootNode(scip), cons, SCIPgetRootNode(scip)) );
      conshdlrdata->rootcons = cons;
   }

   /* check consistency */
   GCGconsOrigbranchCheckConsistency(gcg);

   return SCIP_OKAY;
}

/** checks the consistency of the origbranch constraints in the problem */
void GCGconsOrigbranchCheckConsistency(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
#ifdef CHECKCONSISTENCY
   SCIP* scip;
   SCIP_CONSHDLR*     conshdlr;

#ifndef NDEBUG
   SCIP_CONS** conss;
   int nconss;
   int i;
#endif

   assert(gcg != NULL);
   scip = GCGgetOrigprob(gcg);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("origbranch constraint handler not found\n");
      return;
   }
#ifndef NDEBUG
   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);

   for( i = 0; i < nconss; i++ )
   {
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      assert(consdata->node != NULL);
      assert((consdata->parentcons == NULL) == (SCIPnodeGetDepth(consdata->node) == 0));
      /** todo assertions for all nchildcons */
      if( consdata->nchildconss == 2 )
         assert(consdata->parentcons == NULL || SCIPconsGetData(consdata->parentcons)->childconss[0] == conss[i]
            || SCIPconsGetData(consdata->parentcons)->childconss[1] == conss[i]
            || ( SCIPinProbing(scip) && SCIPconsGetData(consdata->parentcons)->probingtmpcons == conss[i]));
      if( consdata->nchildconss == 2 )
         assert(consdata->childconss[0] == NULL || SCIPconsGetData(consdata->childconss[0])->parentcons == conss[i]);
      if( consdata->nchildconss == 2 )
         assert(consdata->childconss[1] == NULL || SCIPconsGetData(consdata->childconss[1])->parentcons == conss[i]);
      assert(consdata->probingtmpcons == NULL || SCIPinProbing(scip));
      assert(consdata->probingtmpcons == NULL || SCIPconsGetData(consdata->probingtmpcons)->parentcons == conss[i]);
      assert(consdata->mastercons == NULL || GCGconsMasterbranchGetOrigcons(consdata->mastercons) == conss[i]);
   }
#endif
#endif
}
