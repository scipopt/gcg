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
/**@file   cons_branchOrig.c
 * @brief  constraint handler for storing the graph at each node of the tree
 * @author Gerald Gamrath
 *
 */

#include <assert.h>
#include <string.h>

#include "scip/type_cons.h"
#include "scip/cons_linear.h"
#include "cons_branchOrig.h"
#include "probdata_gcg.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "branchOrig"
#define CONSHDLR_DESC          "storing graph at nodes of the tree constraint handler"
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

#define CONSNAMELEN 50

/** constraint data for storing graph constraints */
struct SCIP_ConsData
{
   int                propagatedvars;        /* number of Vars that existed, the last time, the related node was propagated,
                                                used to determine whether the constraint should be repropagated*/
   SCIP_VAR*          origvar;
   GCG_CONSSENSE      sense;
   SCIP_Real          val;
   SCIP_CONS*         pricingcons;
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
SCIP_DECL_CONSFREE(consFreeBranchOrig)
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
SCIP_DECL_CONSINITSOL(consInitsolBranchOrig)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONS*         cons;
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
SCIP_DECL_CONSEXITSOL(consExitsolBranchOrig)
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
SCIP_DECL_CONSDELETE(consDeleteBranchOrig)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(*consdata != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   
   SCIPdebugMessage("Deleting branch orig constraint: <%s>.\n", SCIPconsGetName(cons));

   /* free constraint data */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBranchOrig)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBranchOrig)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckBranchOrig)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockBranchOrig)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   SCIPdebugMessage("Locking method for branch orig constraint: <%s>.\n", SCIPconsGetName(cons));

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveBranchOrig)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA*     consdata;
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("Activating branch orig constraint: <%s> [stack size: %d].\n", SCIPconsGetName(cons), conshdlrData->nstack+1);

   /* put constraint on the stack */
   if ( conshdlrData->nstack >= conshdlrData->maxstacksize )
   {
      SCIPreallocMemoryArray(scip, &(conshdlrData->stack), 2*(conshdlrData->maxstacksize));
      conshdlrData->maxstacksize = 2*(conshdlrData->maxstacksize);
      SCIPdebugMessage("reallocating Memory for stack! %d --> %d\n", conshdlrData->maxstacksize/2, conshdlrData->maxstacksize);
   }
   conshdlrData->stack[conshdlrData->nstack] = cons;
   ++(conshdlrData->nstack);

   /* enable corresponding constraitn in the pricing problem */
   vardata = SCIPvarGetData(consdata->origvar);
   assert(vardata != NULL);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata->blocknr >= 0 && vardata->blocknr < GCGprobGetNPricingprobs(scip));
   assert(vardata->data.origvardata.pricingvar != NULL);

   SCIP_CALL( SCIPdisableCons(GCGprobGetPricingprob(scip, vardata->blocknr), consdata->pricingcons) );

   return SCIP_OKAY;
}



/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveBranchOrig)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA*     consdata;
   SCIP_VARDATA* vardata;

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

   /* disable corresponding constraitn in the pricing problem */
   vardata = SCIPvarGetData(consdata->origvar);
   assert(vardata != NULL);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata->blocknr >= 0 && vardata->blocknr < GCGprobGetNPricingprobs(scip));
   assert(vardata->data.origvardata.pricingvar != NULL);

   SCIP_CALL( SCIPdisableCons(GCGprobGetPricingprob(scip, vardata->blocknr), consdata->pricingcons) );

   SCIPdebugMessage("Deactivating branch orig constraint: <%s> [stack size: %d].\n", SCIPconsGetName(cons), conshdlrData->nstack-1);

   /* remove constraint from the stack */
   --conshdlrData->nstack;

   return SCIP_OKAY;
}



/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropBranchOrig)
{
   SCIP_CONSHDLRDATA* conshdlrData;  
   SCIP_CONS*         cons;
   SCIP_CONSDATA*     consdata;
   SCIP_VAR*          var;
   int                i;
   int                propcount;

   assert(conshdlr != NULL); 
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   *result = SCIP_DIDNOTFIND;
   propcount = 0;

   /* the constraint data of the cons related to the current node */
   cons = conshdlrData->stack[conshdlrData->nstack-1];
   consdata = SCIPconsGetData(cons);

   
   
   SCIPdebugMessage( "Starting propagation of branch orig constraint <%s> %d conss, %d useful.\n", SCIPconsGetName(cons), nconss, nusefulconss);
  
   SCIPdebugMessage( "Finished propagation of branch orig constraint <%s>, %d vars fixed.\n", SCIPconsGetName(cons), propcount);

   consdata = SCIPconsGetData(GCGconsGetActiveBranchOrigCons(scip));

   return SCIP_OKAY;
}

/* define not used callbacks as NULL */
#define consPresolBranchOrig NULL
#define consRespropBranchOrig NULL
#define consInitBranchOrig NULL
#define consExitBranchOrig NULL
#define consInitpreBranchOrig NULL
#define consExitpreBranchOrig NULL
#define consTransBranchOrig NULL
#define consInitlpBranchOrig NULL
#define consSepalpBranchOrig NULL
#define consSepasolBranchOrig NULL
#define consEnableBranchOrig NULL
#define consDisableBranchOrig NULL
#define consPrintBranchOrig NULL


/*
 * interface methods
 */


/** creates the handler for branchOrig constraints and includes it in SCIP */
SCIP_RETCODE GCGincludeConshdlrBranchOrig(
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
         consFreeBranchOrig, consInitBranchOrig, consExitBranchOrig,
         consInitpreBranchOrig, consExitpreBranchOrig, consInitsolBranchOrig, consExitsolBranchOrig,
         consDeleteBranchOrig, consTransBranchOrig, consInitlpBranchOrig,
         consSepalpBranchOrig, consSepasolBranchOrig, consEnfolpBranchOrig, consEnfopsBranchOrig, consCheckBranchOrig,
         consPropBranchOrig, consPresolBranchOrig, consRespropBranchOrig, consLockBranchOrig,
         consActiveBranchOrig, consDeactiveBranchOrig,
         consEnableBranchOrig, consDisableBranchOrig,
         consPrintBranchOrig,
         conshdlrData) );

   return SCIP_OKAY;
}


/** creates and captures a branchOrig constraint*/
SCIP_RETCODE GCGcreateConsBranchOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   SCIP_VAR*             origvar,
   GCG_CONSSENSE         sense,
   SCIP_Real             val
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_VARDATA* vardata;
   SCIP_VAR* pricingvar;
   char* consname;

   assert(scip != NULL);

   /* find the branchOrig constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("branchOrig constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   vardata = SCIPvarGetData(origvar);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata->blocknr >= 0 && vardata->blocknr < GCGprobGetNPricingprobs(scip));
   assert(vardata->data.origvardata.pricingvar != NULL);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->origvar = origvar;
   consdata->sense = sense;
   consdata->val = val;

   (void) SCIPsnprintf(consname, CONSNAMELEN, "%s %s %f", SCIPvarGetName(origvar), (sense == GCG_CONSSENSE_GE ? ">=" : "<="), val);


   /* create constraint in the pricing problem */
   SCIP_CALL( SCIPcreateConsLinear(GCGprobGetPricingprob(scip, vardata->blocknr), &(consdata->pricingcons), consname, 
         0, NULL, NULL, (sense == GCG_CONSSENSE_GE ? val : SCIPinfinity(scip)), (sense == GCG_CONSSENSE_GE ? SCIPinfinity(scip) : val),
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCoefLinear(GCGprobGetPricingprob(scip, vardata->blocknr), consdata->pricingcons, vardata->data.origvardata.pricingvar, 1) );
   SCIP_CALL( SCIPaddCons(GCGprobGetPricingprob(scip, vardata->blocknr), consdata->pricingcons) );
   SCIP_CALL( SCIPdisableCons(GCGprobGetPricingprob(scip, vardata->blocknr), consdata->pricingcons) );

   SCIPdebugMessage("Creating branch orig constraint: <%s>.\n", consname);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, consname, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
         TRUE, FALSE, FALSE, FALSE, TRUE) );

   
   return SCIP_OKAY;
}




/* ----------------------------------- external methods -------------------------- */

/** returns the branch orig constraint of the current node, needs the pointer to the constraint handler */
SCIP_CONS* GCGconsGetActiveBranchOrigConsFromHandler(
   SCIP_CONSHDLR*        conshdlr            /**< constaint handler for branchOrig constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(conshdlr != NULL);
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);
   assert(conshdlrData->nstack > 0);

   return conshdlrData->stack[conshdlrData->nstack-1];
}


/** returns the branch orig constraint of the current node, only needs the pointer to scip */
SCIP_CONS* GCGconsGetActiveBranchOrigCons(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, "branchOrig");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("branchOrig constraint handler not found\n");
      return NULL;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);
   assert(conshdlrData->nstack > 0);

   return conshdlrData->stack[conshdlrData->nstack-1];
}


/** returns the stack and the number of elements on it */
void GCGconsGetStack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, "branchOrig");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("branchOrig constraint handler not found\n");
      return;
   }   
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   *stack = conshdlrData->stack;
   *nstackelements = conshdlrData->nstack;   
}


