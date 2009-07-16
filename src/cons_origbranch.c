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
#include "probdata_gcg.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "origbranch"
#define CONSHDLR_DESC          "store branching decision at nodes of the tree constraint handler"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
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
   SCIP_CONS*         branchcons;            /* constraint in the original problem that forces the 
                                              * branching decision */
   SCIP_VAR*          origvar;               /* original variable on which the branching is done */
   GCG_CONSSENSE      conssense;             /* sense of the branching on the original variable: 
                                                greater-equal (GCG_CONSSENSE_GE) or smaller-equal (GCG_CONSSENSE_LE) */
   SCIP_Real          val;                   /* new lower/upper bound of the original variable */

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
 * Callback methods
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
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   /* prepare stack */
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrData->stack, conshdlrData->maxstacksize) );
   conshdlrData->nstack = 0;

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
   assert(conshdlrData->nstack == 0);
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

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(*consdata != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   
   SCIPdebugMessage("Deleting branch orig constraint: <%s> %s %s %f.\n", SCIPconsGetName(cons), SCIPvarGetName((*consdata)->origvar), 
      ( (*consdata)->conssense == GCG_CONSSENSE_GE ? ">=" : "<=" ), (*consdata)->val);

   /* free constraint data */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
#define consEnfolpOrigbranch NULL


/** constraint enforcing method of constraint handler for pseudo solutions */
#define consEnfopsOrigbranch NULL

/** feasibility check method of constraint handler for integral solutions */
#define consCheckOrigbranch NULL

/** variable rounding lock method of constraint handler */
#define consLockOrigbranch NULL



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

   SCIPdebugMessage("Activating branch orig constraint: <%s> %s %s %f [stack size: %d].\n", SCIPconsGetName(cons), SCIPvarGetName(consdata->origvar), 
      ( consdata->conssense == GCG_CONSSENSE_GE ? ">=" : "<=" ), consdata->val, conshdlrData->nstack+1);

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
   assert(conshdlrData->stack != NULL);
   assert(conshdlrData->nstack > 0);
   assert(cons == conshdlrData->stack[conshdlrData->nstack-1]);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("Deactivating branch orig constraint: <%s> %s %s %f [stack size: %d].\n", SCIPconsGetName(cons), SCIPvarGetName(consdata->origvar), 
      ( consdata->conssense == GCG_CONSSENSE_GE ? ">=" : "<=" ), consdata->val, conshdlrData->nstack-1);

   /* remove constraint from the stack */
   --conshdlrData->nstack;

   return SCIP_OKAY;
}



/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropOrigbranch)
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
   SCIP_CONS*            branchcons,         /**< linear constraint in the original problem */
   SCIP_VAR*             origvar,
   GCG_CONSSENSE         conssense,
   SCIP_Real             val
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);

   /* find the origbranch constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("origbranch constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->branchcons = branchcons;

   consdata->origvar = origvar;
   consdata->conssense = conssense;
   consdata->val = val;

   SCIPdebugMessage("Creating branch orig constraint: <%s>: %s %s %f.\n", name, SCIPvarGetName(origvar), 
      ( conssense == GCG_CONSSENSE_GE ? ">=" : "<=" ), val);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         TRUE, FALSE, FALSE, FALSE, TRUE) );

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
   conshdlr = SCIPfindConshdlr(scip, "origbranch");
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
   conshdlr = SCIPfindConshdlr(scip, "origbranch");
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

/** returns the branch orig constraint of the current node, only needs the pointer to scip */
SCIP_VAR* GCGconsOrigbranchGetOrigvar(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);

   return consdata->origvar;
}

/** returns the branch orig constraint of the current node, only needs the pointer to scip */
GCG_CONSSENSE GCGconsOrigbranchGetConssense(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);

   return consdata->conssense;
}

/** returns the branch orig constraint of the current node, only needs the pointer to scip */
SCIP_Real GCGconsOrigbranchGetVal(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);

   return consdata->val;
}
