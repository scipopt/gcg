/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"
//#define SCIP_DEBUG
/**@file   cons_origbranch.c
 * @brief  constraint handler for storing the branching decisions at each node of the tree
 * @author Gerald Gamrath
 */


#include <assert.h>
#include <string.h>

#include "scip/type_cons.h"
#include "scip/cons_linear.h"
#include "cons_origbranch.h"
#include "relax_gcg.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "origbranch"
#define CONSHDLR_DESC          "store branching decision at nodes of the tree constraint handler"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              * propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/** constraint data for branch orig constraints */
struct SCIP_ConsData
{
   SCIP_NODE*         node;                  /* the node at which the cons is sticking */
   SCIP_CONS*         parentcons;            /* the origbranch constraint of the parent node */
   SCIP_CONS*         child1cons;            /* the origbranch constraint of the first child node */
   SCIP_CONS*         child2cons;            /* the origbranch constraint of the second child node */
   SCIP_CONS*         mastercons;            /* the masterbranch constraint of the corresponding node 
                                              * in the master program */
   GCG_BRANCHDATA*    branchdata;
   SCIP_BRANCHRULE*   branchrule;

   SCIP_VAR**         propvars;              /**< original variable for which the propagation found domain reductions */
   SCIP_BOUNDTYPE*    propboundtypes;        /**< type of the domain new bound found by propagation */
   SCIP_Real*         propbounds;            /**< new lower/upper bound of the propagated original variable */
   int                npropbounds;
   int                maxpropbounds;

};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_CONS**        stack;                 /**< stack for storing active constraints */
   int                nstack;                /**< number of elements on the stack */
   int                maxstacksize;          /**< maximum size of the stack */
};


/*
 * Local methods
 */


/*
 * Callback methods of constraint handler
 */


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeOrigbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   SCIPdebugMessage("freeing branch orig constraint handler\n");

   /* free constraint handler storage */
   assert(conshdlrData->stack == NULL);
   SCIPfreeMemory(scip, &conshdlrData);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolOrigbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONS* cons;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   /* prepare stack */
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrData->stack, conshdlrData->maxstacksize) );
   conshdlrData->nstack = 0;

   assert(SCIPgetRootNode(scip) != NULL);

   SCIP_CALL( GCGcreateConsOrigbranch(scip, &cons, "root-origbranch", SCIPgetRootNode(scip), NULL, NULL, NULL) );

   SCIP_CALL( SCIPaddConsNode(scip, SCIPgetRootNode(scip), cons, SCIPgetRootNode(scip)) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   GCGconsOrigbranchCheckConsistency(scip);

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolOrigbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->nstack == 1);
   SCIPdebugMessage("exiting branch orig constraint handler\n");

   /* free stack */
   SCIPfreeMemoryArray(scip, &conshdlrData->stack);

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteOrigbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA* consdata2;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(*consdata != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);

   SCIPdebugMessage("Deleting branch orig constraint: <%s>.\n", SCIPconsGetName(cons));

   /* set the mastercons pointer of the corresponding origcons to NULL */
   if ( (*consdata)->mastercons != NULL )
      GCGconsMasterbranchSetOrigcons((*consdata)->mastercons, NULL);
   /* set the pointer in parents and children to NULL */
   if ( (*consdata)->parentcons != NULL )
   {
      consdata2 = SCIPconsGetData((*consdata)->parentcons);
      if ( consdata2->child1cons == cons )
      {
         consdata2->child1cons = NULL;
      }
      else
      {
         assert(consdata2->child2cons == cons);
         consdata2->child2cons = NULL;
      }
   }
   assert((*consdata)->child1cons == NULL);
   assert((*consdata)->child2cons == NULL);

   /* delete branchdata, if no mastercons is linked, which would still need the branchdata */
   if ( (*consdata)->mastercons == NULL && (*consdata)->branchdata != NULL )
   {
      SCIP_CALL( GCGrelaxBranchDataDelete(scip, (*consdata)->branchrule, &(*consdata)->branchdata) );
   }

   /* free constraint data */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveOrigbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("Activating branch orig constraint: <%s>[stack size: %d].\n", SCIPconsGetName(cons),
      conshdlrData->nstack+1);


   /* put constraint on the stack */
   if ( conshdlrData->nstack >= conshdlrData->maxstacksize )
   {
      SCIPreallocMemoryArray(scip, &(conshdlrData->stack), 2*(conshdlrData->maxstacksize));
      conshdlrData->maxstacksize = 2*(conshdlrData->maxstacksize);
      SCIPdebugMessage("reallocating Memory for stack! %d --> %d\n", conshdlrData->maxstacksize/2, conshdlrData->maxstacksize);
   }
   conshdlrData->stack[conshdlrData->nstack] = cons;
   ++(conshdlrData->nstack);

   return SCIP_OKAY;
}



/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveOrigbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL || conshdlrData->nstack == 1);
   assert(conshdlrData->nstack > 0);
   assert(conshdlrData->nstack == 1 || cons == conshdlrData->stack[conshdlrData->nstack-1]);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("Deactivating branch orig constraint: <%s> [stack size: %d].\n", 
      SCIPconsGetName(cons), conshdlrData->nstack-1);


   /* remove constraint from the stack */
   --(conshdlrData->nstack);

   return SCIP_OKAY;
}



/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropOrigbranch)
{
   return SCIP_OKAY;
}

static
SCIP_DECL_CONSENFOLP(consEnfolpOrigbranch)
{
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSENFOPS(consEnfopsOrigbranch)
{
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSCHECK(consCheckOrigbranch)
{
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSLOCK(consLockOrigbranch)
{
   return SCIP_OKAY;
}

/* define not used callbacks as NULL */
#define consPresolOrigbranch NULL
#define consRespropOrigbranch NULL
#define consInitOrigbranch NULL
#define consExitOrigbranch NULL
#define consInitpreOrigbranch NULL
#define consExitpreOrigbranch NULL
#define consTransOrigbranch NULL
#define consInitlpOrigbranch NULL
#define consSepalpOrigbranch NULL
#define consSepasolOrigbranch NULL
#define consEnableOrigbranch NULL
#define consDisableOrigbranch NULL
#define consPrintOrigbranch NULL
#define consCopyOrigbranch NULL
#define consParseOrigbranch NULL

/*
 * interface methods
 */


/** creates the handler for origbranch constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrOrigbranch(
   SCIP*                 scip                /**< SCIP data structure */
		       )
{
   SCIP_CONSHDLRDATA* conshdlrData;

   SCIPdebugMessage("Including branch orig constraint handler.\n");

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrData) );
   conshdlrData->stack = NULL;
   conshdlrData->nstack = 0;
   conshdlrData->maxstacksize = 25;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeOrigbranch, consInitOrigbranch, consExitOrigbranch,
         consInitpreOrigbranch, consExitpreOrigbranch, consInitsolOrigbranch, consExitsolOrigbranch,
         consDeleteOrigbranch, consTransOrigbranch, consInitlpOrigbranch,
         consSepalpOrigbranch, consSepasolOrigbranch, consEnfolpOrigbranch, consEnfopsOrigbranch, consCheckOrigbranch,
         consPropOrigbranch, consPresolOrigbranch, consRespropOrigbranch, consLockOrigbranch,
         consActiveOrigbranch, consDeactiveOrigbranch,
         consEnableOrigbranch, consDisableOrigbranch,
         consPrintOrigbranch, consCopyOrigbranch, consParseOrigbranch, 
         conshdlrData) );

   return SCIP_OKAY;
}


/** creates and captures a origbranch constraint*/
SCIP_RETCODE GCGcreateConsOrigbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */          
   SCIP_NODE*            node,
   SCIP_CONS*            parentcons,
   SCIP_BRANCHRULE*      branchrule,
   GCG_BRANCHDATA*       branchdata
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert((parentcons == NULL) == (SCIPnodeGetDepth(node) == 0));

   /* find the origbranch constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("origbranch constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->parentcons = parentcons;
   consdata->node = node;
   consdata->child1cons = NULL;
   consdata->child2cons = NULL;
   consdata->mastercons = NULL;
   consdata->branchrule = branchrule;
   consdata->branchdata = branchdata;
   consdata->npropbounds = 0;
   consdata->maxpropbounds = 0;
   consdata->propvars = NULL;
   consdata->propboundtypes = NULL;
   consdata->propbounds = NULL;


   SCIPdebugMessage("Creating branch orig constraint: <%s>.\n", name);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         TRUE, FALSE, FALSE, FALSE, TRUE) );

   if ( parentcons != NULL )
   {
      SCIP_CONSDATA* parentdata;

      parentdata = SCIPconsGetData(parentcons);
      assert(parentdata != NULL);

      if ( parentdata->child1cons == NULL )
      {
         parentdata->child1cons = *cons;
      }
      else
      {
         assert(parentdata->child2cons == NULL);
         parentdata->child2cons = *cons;
      }
   }


   return SCIP_OKAY;
}




/* ----------------------------------- external methods -------------------------- */

/** returns the branch orig constraint of the current node, only needs the pointer to scip */
SCIP_CONS* GCGconsOrigbranchGetActiveCons(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("origbranch constraint handler not found\n");
      return NULL;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);
   assert(conshdlrData->nstack > 0);

   return conshdlrData->stack[conshdlrData->nstack-1];
}


/** returns the stack and the number of elements on it */
void GCGconsOrigbranchGetStack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("origbranch constraint handler not found\n");
      return;
   }   
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   *stack = conshdlrData->stack;
   *nstackelements = conshdlrData->nstack;   

}

/** returns the branching data for a given origbranch constraint */
GCG_BRANCHDATA* GCGconsOrigbranchGetBranchdata(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->branchdata;
}

/** returns the branchrule for a given origbranch constraint */
SCIP_BRANCHRULE* GCGconsOrigbranchGetBranchrule(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->branchrule;
}

/** returns the node in the B&B tree at which the given origbranch constraint is sticking */
SCIP_NODE* GCGconsOrigbranchGetNode(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->node;
}

/** returns the origbranch constraint of the B&B father of the node at which the 
    given origbranch constraint is sticking */
SCIP_CONS* GCGconsOrigbranchGetParentcons(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->parentcons;
}

/** returns the origbranch constraint of the first child of the node at which the 
    given origbranch constraint is sticking */
SCIP_CONS* GCGconsOrigbranchGetChild1cons(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->child1cons;
}

/** returns the origbranch constraint of the second child of the node at which the 
    given origbranch constraint is sticking */
SCIP_CONS* GCGconsOrigbranchGetChild2cons(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->child2cons;
}

/** sets the masterbranch constraint of the node in the master program corresponding to the node 
    at which the given origbranchbranch constraint is sticking */
void GCGconsOrigbranchSetMastercons(
   SCIP_CONS*            cons,
   SCIP_CONS*            mastercons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->mastercons == NULL || mastercons == NULL);

   consdata->mastercons = mastercons;
}

/** returns the masterbranch constraint of the node in the master program corresponding to the node 
    at which the given origbranchbranch constraint is sticking */
SCIP_CONS* GCGconsOrigbranchGetMastercons(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->mastercons != NULL);

   return consdata->mastercons;
}

/** checks the consistency of the origbranch constraints in the problem */
void GCGconsOrigbranchCheckConsistency(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONS** conss;
   SCIP_CONSDATA* consdata;
   int nconss;
   int i;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("origbranch constraint handler not found\n");
      return;
   }   

   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);

   for ( i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      assert(consdata->node != NULL);
      assert((consdata->parentcons == NULL) == (SCIPnodeGetDepth(consdata->node) == 0));
      assert(consdata->parentcons == NULL || SCIPconsGetData(consdata->parentcons)->child1cons == conss[i]
         || SCIPconsGetData(consdata->parentcons)->child2cons == conss[i]);
      assert(consdata->child1cons == NULL || SCIPconsGetData(consdata->child1cons)->parentcons == conss[i]);
      assert(consdata->child2cons == NULL || SCIPconsGetData(consdata->child2cons)->parentcons == conss[i]);
      assert(consdata->mastercons == NULL || 
         GCGconsMasterbranchGetOrigcons(consdata->mastercons) == conss[i]);
   }

   //SCIPdebugMessage("checked consistency of %d origbranch constraints, all ok!\n", nconss);
}


SCIP_RETCODE GCGconsOrigbranchAddPropBoundChg(
   SCIP*                 scip,
   SCIP_CONS*            cons, 
   SCIP_VAR*             var,
   SCIP_BOUNDTYPE        boundtype,
   SCIP_Real             newbound
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->npropbounds >= consdata->maxpropbounds )
   {
      consdata->maxpropbounds = consdata->npropbounds+5;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(consdata->propvars), consdata->maxpropbounds) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(consdata->propboundtypes), consdata->maxpropbounds) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(consdata->propbounds), consdata->maxpropbounds) );
   }
   
   SCIPdebugMessage("Bound change stored at branch orig constraint: <%s>.\n", SCIPconsGetName(cons));


   consdata->propvars[consdata->npropbounds] = var;
   consdata->propboundtypes[consdata->npropbounds] = boundtype;
   consdata->propbounds[consdata->npropbounds] = newbound;
   consdata->npropbounds++;
 
   return SCIP_OKAY;
}

SCIP_RETCODE GCGconsOrigbranchGetPropBoundChgs(
   SCIP*                 scip,
   SCIP_CONS*            cons, 
   SCIP_VAR***           vars,
   SCIP_BOUNDTYPE**      boundtypes,
   SCIP_Real**           newbounds,
   int*                  npropbounds
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *vars = consdata->propvars;
   *boundtypes = consdata->propboundtypes;
   *newbounds = consdata->propbounds;
   *npropbounds = consdata->npropbounds;

   consdata->npropbounds = 0;
 
   return SCIP_OKAY;
}


