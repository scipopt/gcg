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
/**@file   cons_masterbranch.c
 * @brief  constraint handler for storing the branching decisions at each node of the tree
 * @author Gerald Gamrath
 *
 */

#include <assert.h>
#include <string.h>

#include "scip/type_cons.h"
#include "scip/cons_linear.h"
#include "cons_masterbranch.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "cons_origbranch.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "masterbranch"
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
   int                propagatedvars;        /* number of Vars that existed, the last time, the related node was propagated,
                                                used to determine whether the constraint should be repropagated */
   SCIP_VAR*          origvar;               /* original variable on which the branching is done */
   GCG_CONSSENSE      conssense;             /* sense of the branching on the original variable: 
                                                greater-equal (GCG_CONSSENSE_GE) or smaller-equal (GCG_CONSSENSE_LE) */
   SCIP_Real          val;                   /* new lower/upper bound of the original variable */
   SCIP_CONS*         pricingcons;
   SCIP_Bool          created;
   SCIP_NODE*         node;                  /* the node at which the cons is sticking */
   SCIP_CONS*         parentcons;            /* the masterbranch constraint of the parent node */
   SCIP_CONS*         child1cons;            /* the masterbranch constraint of the first child node */
   SCIP_CONS*         child2cons;            /* the masterbranch constraint of the second child node */
   SCIP_CONS*         origcons;              /* the corresponding origbranch cons in the original program */
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

#if 0
static
void printVar(
   SCIP_VAR*             var
   )
{
   SCIP_VARDATA* vardata;
   int i;

   vardata = SCIPvarGetData(var);

   printf("name = %s, ", SCIPvarGetName(var));
   if ( vardata->data.mastervardata.origvals[0] > 10000 )
   {
      printf("vals = (%s: inf", SCIPvarGetName(vardata->data.mastervardata.origvars[0]));
   }
   else
   {
      printf("vals = (%s: %f", SCIPvarGetName(vardata->data.mastervardata.origvars[0]), 
         vardata->data.mastervardata.origvals[0]);
   }
   for ( i = 1; i < vardata->data.mastervardata.norigvars; i++ )
   {
      if ( vardata->data.mastervardata.origvals[i] > 10000 )
      {
         printf("; %s: inf", SCIPvarGetName(vardata->data.mastervardata.origvars[i]));
      }
      else
      {
         printf("; %s: %f", SCIPvarGetName(vardata->data.mastervardata.origvars[i]), 
            vardata->data.mastervardata.origvals[i]);
      }
   }
   if ( SCIPvarGetUbLocal(var) < 10000 )
      printf("), ub = %f", SCIPvarGetUbLocal(var));
   else
      printf("), ub = inf");
   printf("\n");
}
#endif

/*
 * Callback methods
 */


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   SCIPdebugMessage("freeing masterbranch constraint handler\n");

   /* free constraint handler storage */
   assert(conshdlrData->stack == NULL);
   SCIPfreeMemory(scip, &conshdlrData);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolMasterbranch)
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

   SCIPdebugMessage("consInitsolMasterbranch()\n");

   assert(SCIPgetRootNode(scip) != NULL);

   SCIP_CALL( GCGcreateConsMasterbranch(scip, &cons, SCIPgetRootNode(scip), NULL) );

   SCIP_CALL( SCIPaddConsNode(scip, SCIPgetRootNode(scip), cons, SCIPgetRootNode(scip)) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->nstack == 1);
   SCIPdebugMessage("exiting masterbranch constraint handler\n");

   /* free stack */
   SCIPfreeMemoryArray(scip, &conshdlrData->stack);

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(*consdata != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);

   SCIPdebugMessage("Deleting masterbranch constraint: <%s>.\n", SCIPconsGetName(cons));

   /* free constraint data */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA*     consdata;
   SCIP_VARDATA* vardata;
   SCIP_CONS* origcons;
   SCIP* origscip;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->node != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   if ( !consdata->created )
   {
      origcons = GCGconsOrigbranchGetActiveCons(origscip);
      assert(origcons != NULL);

      consdata->origvar = GCGconsOrigbranchGetOrigvar(origcons);
      consdata->val = GCGconsOrigbranchGetVal(origcons);
      consdata->conssense = GCGconsOrigbranchGetConssense(origcons);
      consdata->origcons = origcons;
      GCGconsOrigbranchSetMastercons(origcons, scip, cons);

      assert(SCIPgetCurrentNode(scip) == consdata->node || consdata->node == SCIPgetRootNode(scip));
      assert((SCIPgetNNodesLeft(scip)+SCIPgetNNodes(scip) == 1) == (consdata->node == SCIPgetRootNode(scip)));
      assert(SCIPnodeGetDepth(GCGconsOrigbranchGetNode(consdata->origcons)) == SCIPnodeGetDepth(consdata->node));
      assert((consdata->node == SCIPgetRootNode(scip)) == (SCIPnodeGetDepth(consdata->node) == 0));
      assert(consdata->parentcons != NULL || SCIPnodeGetDepth(consdata->node) == 0);
      assert(consdata->parentcons == NULL || 
         SCIPconsGetData(consdata->parentcons)->origcons == GCGconsOrigbranchGetParentcons(consdata->origcons));
      
      /* round current value of the branching-variable */
      if ( consdata->conssense == GCG_CONSSENSE_GE )
      {
         consdata->val = SCIPceil(scip, consdata->val);
      }
      else
      {
         consdata->val = SCIPfloor(scip, consdata->val);
      }

      consdata->created = TRUE;
   }

   /* put constraint on the stack */
   if ( conshdlrData->nstack >= conshdlrData->maxstacksize )
   {
      SCIPreallocMemoryArray(scip, &(conshdlrData->stack), 2*(conshdlrData->maxstacksize));
      conshdlrData->maxstacksize = 2*(conshdlrData->maxstacksize);
      SCIPdebugMessage("reallocating Memory for stack! %d --> %d\n", conshdlrData->maxstacksize/2, conshdlrData->maxstacksize);
   }
   conshdlrData->stack[conshdlrData->nstack] = cons;
   (conshdlrData->nstack)++;

   if ( consdata->conssense == GCG_CONSSENSE_NONE )
   {
      SCIPdebugMessage("Activating masterbranch constraint at root: <%s> [stack size: %d].\n", 
         SCIPconsGetName(cons), conshdlrData->nstack);

      return SCIP_OKAY;
   }

   SCIPdebugMessage("Activating masterbranch constraint: <%s> %s %s %f [stack size: %d].\n", SCIPconsGetName(cons), 
      SCIPvarGetName(consdata->origvar), (consdata->conssense == GCG_CONSSENSE_GE ? ">=" : "<="), consdata->val, conshdlrData->nstack);

   /* create corresponding constraint in the pricing problem */
   vardata = SCIPvarGetData(consdata->origvar);
   assert(vardata != NULL);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata->blocknr >= 0 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
   assert(vardata->data.origvardata.pricingvar != NULL);

   /* create constraint in the pricing problem */
   SCIP_CALL( SCIPcreateConsLinear(GCGrelaxGetPricingprob(origscip, vardata->blocknr), &(consdata->pricingcons), 
         SCIPconsGetName(cons), 0, NULL, NULL, 
         (consdata->conssense == GCG_CONSSENSE_GE ? consdata->val : -1.0 * SCIPinfinity(scip)), 
         (consdata->conssense == GCG_CONSSENSE_GE ? SCIPinfinity(scip) : consdata->val), 
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCoefLinear(GCGrelaxGetPricingprob(origscip, vardata->blocknr), consdata->pricingcons, vardata->data.origvardata.pricingvar, 1) );
   SCIP_CALL( SCIPaddCons(GCGrelaxGetPricingprob(origscip, vardata->blocknr), consdata->pricingcons) );

#if 0
   SCIP_VAR** vars;
   vars = SCIPgetVars(scip);
  
   
   if ( SCIPgetDepth(scip) == 1 )
   {
      for ( i = 0; i < SCIPgetNVars(scip); i++)
      {
         printVar(vars[i]);
      }
   }
#endif

   return SCIP_OKAY;
}



/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA*     consdata;
   SCIP_VARDATA* vardata;
   SCIP* origscip;

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
   assert(consdata->created);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   if ( consdata->conssense != GCG_CONSSENSE_NONE )
   {
      /* disable corresponding constraint in the pricing problem */
      vardata = SCIPvarGetData(consdata->origvar);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->blocknr >= 0 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
      assert(vardata->data.origvardata.pricingvar != NULL);
      
      SCIP_CALL( SCIPdelCons(GCGrelaxGetPricingprob(origscip, vardata->blocknr), consdata->pricingcons) );

      SCIPdebugMessage("Deactivating masterbranch constraint: <%s> %s %s %f [stack size: %d].\n", SCIPconsGetName(cons), 
         SCIPvarGetName(consdata->origvar), (consdata->conssense == GCG_CONSSENSE_GE ? ">=" : "<="), consdata->val, conshdlrData->nstack-1);
   }
   else
   {
      SCIPdebugMessage("Deactivating masterbranch constraint at root: <%s> [stack size: %d].\n", SCIPconsGetName(cons), 
         conshdlrData->nstack-1);
   }

   /* remove constraint from the stack */
   (conshdlrData->nstack)--;

   return SCIP_OKAY;
}



/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;  
   SCIP_CONS*         cons;
   SCIP_CONSDATA*     consdata;
   SCIP_VAR**         vars;
   SCIP_VARDATA*      vardata;
   SCIP* origscip;
   int                i;
   int                j;
   int                c;
   int                propcount;
   

   assert(conshdlr != NULL); 
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   *result = SCIP_DIDNOTFIND;
   propcount = 0;

   for ( c = 0; c < conshdlrData->nstack; c++ )
   {
      /* the constraint data of the cons related to the current node */
      cons = conshdlrData->stack[c];
      consdata = SCIPconsGetData(cons);

      if ( consdata->conssense == GCG_CONSSENSE_NONE )
         continue;
   
      SCIPdebugMessage("Starting propagation of masterbranch constraint: <%s> %s %s %f.\n", SCIPconsGetName(cons), 
         SCIPvarGetName(consdata->origvar), (consdata->conssense == GCG_CONSSENSE_GE ? ">=" : "<="), consdata->val);

      vars = SCIPgetVars(scip);
      
      for ( i = 0; i < SCIPgetNVars(scip); i++)
      {
         //printVar(vars[i]);
         
         if ( !SCIPisFeasZero(scip, SCIPvarGetUbLocal(vars[i])) && SCIPvarGetData(vars[i]) != NULL)
         {
            vardata = SCIPvarGetData(vars[i]);
            assert(vardata->vartype == GCG_VARTYPE_MASTER);
            assert(vardata->blocknr >= -1 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
            assert(vardata->data.mastervardata.norigvars > 0);
            assert(vardata->data.mastervardata.origvals != NULL);
            assert(vardata->data.mastervardata.origvars != NULL);
            
            for ( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
            {
               if ( vardata->data.mastervardata.origvars[j] == consdata->origvar )
               {
                  if ( consdata->conssense == GCG_CONSSENSE_GE && 
                     SCIPisFeasLT(scip, vardata->data.mastervardata.origvals[j], consdata->val) )
                  {
                     SCIPchgVarUb(scip, vars[i], 0.0);
                     propcount++;
                     //printf("var %s: upper bound set to 0 (ge)\n", SCIPvarGetName(vars[i]));
                     //printVar(vars[i]);
                     break;
                  }
                  if ( consdata->conssense == GCG_CONSSENSE_LE && 
                     SCIPisFeasGT(scip, vardata->data.mastervardata.origvals[j], consdata->val) )
                  {
                     SCIPchgVarUb(scip, vars[i], 0.0);
                     propcount++;
                     //printf("var %s: upper bound set to 0 (le)\n", SCIPvarGetName(vars[i]));
                     //printVar(vars[i]);
                     break;
                  }
                  
               }
               
            }
            /*
            if ( j == vardata->data.mastervardata.norigvars && consdata->conssense == GCG_CONSSENSE_GE )
            {
               SCIPchgVarUb(scip, vars[i], 0.0);
               //printf("var %s: upper bound set to 0 (ge)\n", SCIPvarGetName(vars[i]));
               //printVar(vars[i]);
               
            }
            */
         }
      }
      
      SCIPdebugMessage("Finished propagation of masterbranch constraint: <%s> %s %s %f, %d vars fixed.\n", SCIPconsGetName(cons), 
         SCIPvarGetName(consdata->origvar), (consdata->conssense == GCG_CONSSENSE_GE ? ">=" : "<="), consdata->val, propcount);

   }
   //consdata = SCIPconsGetData(GCGconsGetActiveMasterbranchCons(scip));

   return SCIP_OKAY;
}

/* define not used callbacks as NULL */
#define consEnfolpMasterbranch NULL
#define consEnfopsMasterbranch NULL
#define consCheckMasterbranch NULL
#define consLockMasterbranch NULL
#define consPresolMasterbranch NULL
#define consRespropMasterbranch NULL
#define consInitMasterbranch NULL
#define consExitMasterbranch NULL
#define consInitpreMasterbranch NULL
#define consExitpreMasterbranch NULL
#define consTransMasterbranch NULL
#define consInitlpMasterbranch NULL
#define consSepalpMasterbranch NULL
#define consSepasolMasterbranch NULL
#define consEnableMasterbranch NULL
#define consDisableMasterbranch NULL
#define consPrintMasterbranch NULL
#define consCopyMasterbranch NULL
#define consParseMasterbranch NULL


/*
 * interface methods
 */


/** creates the handler for masterbranch constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrMasterbranch(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrData;

   SCIPdebugMessage("Including masterbranch constraint handler.\n");

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrData) );
   conshdlrData->stack = NULL;
   conshdlrData->nstack = 0;
   conshdlrData->maxstacksize = 25;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeMasterbranch, consInitMasterbranch, consExitMasterbranch,
         consInitpreMasterbranch, consExitpreMasterbranch, consInitsolMasterbranch, consExitsolMasterbranch,
         consDeleteMasterbranch, consTransMasterbranch, consInitlpMasterbranch,
         consSepalpMasterbranch, consSepasolMasterbranch, consEnfolpMasterbranch, consEnfopsMasterbranch, consCheckMasterbranch,
         consPropMasterbranch, consPresolMasterbranch, consRespropMasterbranch, consLockMasterbranch,
         consActiveMasterbranch, consDeactiveMasterbranch,
         consEnableMasterbranch, consDisableMasterbranch,
         consPrintMasterbranch, consCopyMasterbranch, consParseMasterbranch, 
         conshdlrData) );

   return SCIP_OKAY;
}


/** creates and captures a masterbranch constraint*/
SCIP_RETCODE GCGcreateConsMasterbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   SCIP_NODE*            node,
   SCIP_CONS*            parentcons
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(node != NULL);
   assert((parentcons == NULL) == (SCIPnodeGetDepth(node) == 0));


   /* find the masterbranch constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->node = node;
   consdata->parentcons = parentcons;
   consdata->child1cons = NULL;
   consdata->child2cons = NULL;
   consdata->created = FALSE;
   consdata->origcons = NULL;

   SCIPdebugMessage("Creating masterbranch constraint.\n");

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, "masterbranch", conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
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

/** returns the masterbranch constraint of the current node */
SCIP_CONS* GCGconsMasterbranchGetActiveCons(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return NULL;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);
   assert(conshdlrData->nstack > 0);

   return conshdlrData->stack[conshdlrData->nstack-1];
}


/** returns the stack and the number of elements on it */
void GCGconsMasterbranchGetStack(
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
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return;
   }   
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   *stack = conshdlrData->stack;
   *nstackelements = conshdlrData->nstack;

}


/** returns the original variable for a given masterbranch constraint */
SCIP_VAR* GCGconsMasterbranchGetOrigvar(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);

   return consdata->origvar;
}

/** returns the conssense for a given masterbranch constraint */
GCG_CONSSENSE GCGconsMasterbranchGetConssense(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);

   return consdata->conssense;
}

/** returns the new bound for a given masterbranch constraint */
SCIP_Real GCGconsMasterbranchGetVal(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);

   return consdata->val;
}

/** returns the node in the B&B tree at which the given masterbranch constraint is sticking */
SCIP_NODE* GCGconsMasterbranchGetNode(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);

   return consdata->node;
}

/** returns the masterbranch constraint of the B&B father of the node at which the 
    given masterbranch constraint is sticking */
SCIP_CONS* GCGconsMasterbranchGetParentcons(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);

   return consdata->parentcons;
}

/** returns the masterbranch constraint of the first child of the node at which the 
    given masterbranch constraint is sticking */
SCIP_CONS* GCGconsMasterbranchGetChild1cons(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);

   return consdata->child1cons;
}

/** returns the masterbranch constraint of the second child of the node at which the 
    given masterbranch constraint is sticking */
SCIP_CONS* GCGconsMasterbranchGetChild2cons(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);

   return consdata->child2cons;
}

/** returns the origbranch constraint of the node in the original program corresponding to the node 
    at which the given masterbranch constraint is sticking */
SCIP_CONS* GCGconsMasterbranchGetOrigcons(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);

   return consdata->origcons;
}

/** checks the consistency of the masterbranch constraints in the problem */
void GCGconsMasterbranchCheckConsistency(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONS** conss;
   SCIP_CONSDATA* consdata;
   int nconss;
   int i;

   if ( scip == NULL )
      return;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      assert(0);
      return;
   }   

   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);

   for ( i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata->node != NULL);
      assert((consdata->parentcons == NULL) == (SCIPnodeGetDepth(consdata->node) == 0));
      assert((consdata->origcons == NULL) == !consdata->created);
      assert(consdata->parentcons == NULL || SCIPconsGetData(consdata->parentcons)->child1cons == conss[i]
         || SCIPconsGetData(consdata->parentcons)->child2cons == conss[i]);
      assert(consdata->child1cons == NULL || SCIPconsGetData(consdata->child1cons)->parentcons == conss[i]);
      assert(consdata->child2cons == NULL || SCIPconsGetData(consdata->child2cons)->parentcons == conss[i]);
      assert(consdata->origcons == NULL || 
         GCGconsOrigbranchGetMastercons(consdata->origcons) == conss[i]);
   }

   //SCIPdebugMessage("checked consistency of %d masterbranch constraints, all ok!\n", nconss);
}
