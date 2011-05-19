/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"
//#define SCIP_DEBUG
//#define CHECKPROPAGATEDVARS
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
#include "cons_origbranch.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"

#include "struct_vardata.h"


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


#define EVENTHDLR_NAME         "origvarbound"
#define EVENTHDLR_DESC         "event handler for origvarbound event"


/** constraint data for masterbranch constraints */
struct SCIP_ConsData
{
   int                propagatedvars;        /**< number of Vars that existed, the last time, the related node was propagated,
                                                  used to determine whether the constraint should be repropagated */
   SCIP_Bool          needprop;              /**< should the constraint be propagated? */
   SCIP_Bool          created;
   SCIP_NODE*         node;                  /**< the node at which the cons is sticking */
   SCIP_CONS*         parentcons;            /**< the masterbranch constraint of the parent node */
   SCIP_CONS*         child1cons;            /**< the masterbranch constraint of the first child node */
   SCIP_CONS*         child2cons;            /**< the masterbranch constraint of the second child node */
   SCIP_CONS*         probingtmpcons;        /**< pointer to save the second child if the child2cons pointer is overwritten in probing mode */
   SCIP_CONS*         origcons;              /**< the corresponding origbranch cons in the original program */

   GCG_BRANCHDATA*    branchdata;            /**< branching data stored by the branching rule at the corresponding origcons constraint 
                                              *   containing information about the branching restrictions */
   SCIP_BRANCHRULE*   branchrule;            /**< branching rule that created the corresponding node in the original problem and imposed 
                                              *   branching restrictions */
   
   SCIP_VAR**         boundchgvars;          /**< variables of bound changes stored at the current node */
   SCIP_Real*         newbounds;             /**< new bounds for the bound changes stored at the current node */
   SCIP_Real*         oldbounds;             /**< old bounds for the bound changes stored at the current node */
   SCIP_BOUNDTYPE*    boundtypes;            /**< types of the bound changes stored at the current node */
   int*               nboundchangestreated;  /**< number of bound changes of the nodes on the way from the current node to 
                                              *   the root node that are treated so far */
   int                nboundchanges;         /**< number of bound changes */
   int                nbranchingchanges;     /**< number of bound changes due to branching (<= nboundchanges) */
   int                nactivated;            /**< number of times the constraint was activated so far */
   char*              name;                  /** name of the constraint */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_CONS**        stack;                 /**< stack for storing active constraints */
   int                nstack;                /**< number of elements on the stack */
   int                maxstacksize;          /**< maximum size of the stack */
   SCIP_VAR**         pendingvars;           /**< variables corresponding to pending bound changes (global bound changes) */
   SCIP_BOUNDTYPE*    pendingbndtypes;       /**< types of the pending bound changes (global bound changes) */
   SCIP_Real*         pendingnewbnds;        /**< new bounds corresponding to pending bound changes (global bound changes) */
   SCIP_Real*         pendingoldbnds;        /**< old bounds corresponding to pending bound changes (global bound changes) */
   int                npendingbnds;          /**< number of pending bound changes (global bound changes) */ 
   SCIP_Bool          pendingbndsactivated;  /**< were pending bound changes already activated? */
   int                maxpendingbnds;        /**< size of the array corresponding to pending bound changes */
   SCIP_Bool          enforceproper;         /**< should proper variables be enforced? */
};


/** event handler data */
struct SCIP_EventhdlrData
{
};

/*
 * Local methods
 */

/* local method to add a global bound change to the pending bound changes array */
static
SCIP_RETCODE GCGconsMasterbranchAddPendingBndChg(
   SCIP*                 scip, 
   SCIP_VAR*             var,
   SCIP_BOUNDTYPE        boundtype,
   SCIP_Real             oldbound,
   SCIP_Real             newbound
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }   
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert((conshdlrData->npendingbnds > 0) || conshdlrData->pendingbndsactivated);

   /* realloc memory if needed */
   if( conshdlrData->npendingbnds >= conshdlrData->maxpendingbnds )
   {
      conshdlrData->maxpendingbnds = conshdlrData->npendingbnds+5;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrData->pendingvars), conshdlrData->maxpendingbnds) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrData->pendingbndtypes), conshdlrData->maxpendingbnds) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrData->pendingoldbnds), conshdlrData->maxpendingbnds) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrData->pendingnewbnds), conshdlrData->maxpendingbnds) );
   }

   /* store pending bound change */
   conshdlrData->pendingvars[conshdlrData->npendingbnds] = var;
   conshdlrData->pendingbndtypes[conshdlrData->npendingbnds] = boundtype;
   conshdlrData->pendingoldbnds[conshdlrData->npendingbnds] = oldbound;
   conshdlrData->pendingnewbnds[conshdlrData->npendingbnds] = newbound;
   conshdlrData->npendingbnds++;
   conshdlrData->pendingbndsactivated = FALSE;
 
   return SCIP_OKAY;
}

#ifdef CHECKPROPAGATEDVARS
/** method to check whether all master variables that violate the bounds in the current original problem are fixed to 0 */
static
SCIP_Bool checkVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< data structure of the masterbranch constraint handler */
   SCIP_Bool             printall            /**< should all violations be printed or only the first one? */
   )
{
   SCIP_CONSHDLRDATA* conshdlrData;  
   SCIP_CONS*         cons;
   SCIP_CONSDATA*     consdata;
   SCIP_VAR**         vars;
   int                nvars;
   SCIP_VARDATA*      vardata;
   SCIP*              origscip;
   int                i;
   int                j;
   int                c;

   assert(conshdlr != NULL); 
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   vars = SCIPgetVars(scip);
   assert(vars != NULL);
   nvars = SCIPgetNVars(scip);

   printf("checkVars()\n");

   /* first of all, check whether variables not fixed to 0 are really valid for the current node */
   /* iterate over all constraints */
   for( c = 0; c < conshdlrData->nstack; c++ )
   {
      cons = conshdlrData->stack[c];
      consdata = SCIPconsGetData(cons);

      if( consdata->branchrule == NULL )
         continue;
           
      /* iterate over all vars and check whether they violate the current cons */
      for( i = 0; i < nvars; i++)
      {
         if( !SCIPisFeasZero(scip, SCIPvarGetUbLocal(vars[i])) )
         {
            vardata = SCIPvarGetData(vars[i]);
            assert(vardata != NULL);
            assert(vardata->vartype == GCG_VARTYPE_MASTER);
            assert(vardata->blocknr >= -1 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
            assert(vardata->data.mastervardata.norigvars >= 0);
            assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
            assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);
            
            for( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
            {
               if( vardata->data.mastervardata.origvars[j] == consdata->origvar )
               {
                  if( consdata->conssense == GCG_CONSSENSE_GE && 
                     SCIPisFeasLT(scip, vardata->data.mastervardata.origvals[j], consdata->val) )
                  {
                     printf("var %s: upper bound should be fixed to 0 because of cons %s [c=%d], but it is not!\n", SCIPvarGetName(vars[i]), SCIPconsGetName(cons), c);
                     printf("--> Reason: origvars[j] = %s >= origvals[j] = %g violated!\n",
                        SCIPvarGetName(vardata->data.mastervardata.origvars[j]), vardata->data.mastervardata.origvals[j]);
                     if( !printall )
                        return FALSE;
                  }
                  if( consdata->conssense == GCG_CONSSENSE_LE && 
                     SCIPisFeasGT(scip, vardata->data.mastervardata.origvals[j], consdata->val) )
                  {
                     printf("var %s: upper bound should be fixed to 0 because of cons %s [c=%d], but it is not!\n", SCIPvarGetName(vars[i]), SCIPconsGetName(cons), c);
                     printf("--> Reason: origvars[j] = %s <= origvals[j] = %g violated!\n",
                        SCIPvarGetName(vardata->data.mastervardata.origvars[j]), vardata->data.mastervardata.origvals[j]);                       
                     if( !printall )
                        return FALSE;
                  }
                  
               }
               
            }
         }
      }
   }
   /* now check for all variables fixed to 0, whether there is a reason for this fixing active at the current node */

   return TRUE;   
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


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   SCIPdebugMessage("consInitMasterbranch()\n");

   /* prepare stack */
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrData->stack, conshdlrData->maxstacksize) );
   conshdlrData->nstack = 0;

   /* prepare pending bound changes */
   conshdlrData->npendingbnds = 0;
   conshdlrData->maxpendingbnds = 5;
   conshdlrData->pendingbndsactivated = TRUE;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrData->pendingvars), conshdlrData->maxpendingbnds) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrData->pendingbndtypes), conshdlrData->maxpendingbnds) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrData->pendingoldbnds), conshdlrData->maxpendingbnds) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrData->pendingnewbnds), conshdlrData->maxpendingbnds) );

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

   SCIPdebugMessage("consInitsolMasterbranch()\n");

   /* create masterbranch constraint for the root node */
   assert(SCIPgetRootNode(scip) != NULL);
   SCIP_CALL( GCGcreateConsMasterbranch(scip, &cons, SCIPgetRootNode(scip), NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, SCIPgetRootNode(scip), cons, SCIPgetRootNode(scip)) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitMasterbranch)
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
   SCIPfreeMemoryArray(scip, &(conshdlrData->stack));
   SCIPfreeMemoryArray(scip, &(conshdlrData->pendingvars));
   SCIPfreeMemoryArray(scip, &(conshdlrData->pendingbndtypes));
   SCIPfreeMemoryArray(scip, &(conshdlrData->pendingoldbnds));
   SCIPfreeMemoryArray(scip, &(conshdlrData->pendingnewbnds));

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteMasterbranch)
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

   SCIPdebugMessage("Deleting masterbranch constraint: <%s>.\n", (*consdata)->name);

   /* set the mastercons pointer of the corresponding origcons to NULL */
   if( (*consdata)->origcons != NULL )
      GCGconsOrigbranchSetMastercons((*consdata)->origcons, NULL);
   /* set the pointer in the parent node to NULL */
   if( (*consdata)->parentcons != NULL )
   {
      consdata2 = SCIPconsGetData((*consdata)->parentcons);
      if( consdata2->child1cons == cons )
      {
         consdata2->child1cons = NULL;
      }
      else
      {
         if( consdata2->child2cons == cons )
         {
            consdata2->child2cons = NULL;

            if( SCIPinProbing(scip) )
            {
               consdata2->child2cons = consdata2->probingtmpcons;
               consdata2->probingtmpcons = NULL;
            }
         }
         else
         {
            assert(SCIPinProbing(scip));
            assert(consdata2->probingtmpcons == cons);
            assert(SCIPisLE(scip, SCIPgetCutoffbound(scip), SCIPgetNodeLowerbound(scip, (*consdata)->node)));

            consdata2->probingtmpcons = NULL;
         }
      }
   }
   /* the node should not have children anymore */
   assert((*consdata)->child1cons == NULL);
   assert((*consdata)->child2cons == NULL);

   /* delete branchdata, if the corresponding origcons was already deleted, otherwise, it will be deleted by the 
    * corresponding origbranch constraint */
   if( (*consdata)->origcons == NULL && (*consdata)->branchdata != NULL )
   {
      SCIP_CALL( GCGrelaxBranchDataDelete(GCGpricerGetOrigprob(scip), (*consdata)->branchrule, &(*consdata)->branchdata) );
   }

   /* delete array with bound changes */
   if( (*consdata)->nboundchanges > 0 )
   {
      SCIPfreeMemoryArray(scip, &(*consdata)->oldbounds);
      SCIPfreeMemoryArray(scip, &(*consdata)->newbounds);
      SCIPfreeMemoryArray(scip, &(*consdata)->boundtypes);
      SCIPfreeMemoryArray(scip, &(*consdata)->boundchgvars);
   }

   if( (*consdata)->nboundchangestreated != NULL )
   {
      SCIPfreeMemoryArray(scip, &(*consdata)->nboundchangestreated);
   }

   /* free constraint data */
   if( (*consdata)->name != NULL )
   {
      BMSfreeBlockMemoryArray(SCIPblkmem(scip), &(*consdata)->name, strlen((*consdata)->name)+1);
   }
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveMasterbranch)
{
   SCIP* origscip;
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* origcons;

   SCIP_CONSDATA* parentdata;
   SCIP_VARDATA* vardata;
   int i;
   int j;

   SCIP_DOMCHG* domchg;
   SCIP_BOUNDCHG* boundchg;

   SCIP_CONSDATA* stackconsdata;


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

   consdata->nactivated++;

   /* if the node is activated the first time, we first have to setup the constraint data */
   if( !consdata->created )
   {
      origcons = GCGconsOrigbranchGetActiveCons(origscip);
      assert(origcons != NULL);

      consdata->origcons = origcons;
      consdata->branchrule = GCGconsOrigbranchGetBranchrule(origcons);
      consdata->branchdata = GCGconsOrigbranchGetBranchdata(origcons);
      GCGconsOrigbranchSetMastercons(origcons, cons);

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(SCIPblkmem(scip), &consdata->name, SCIPconsGetName(consdata->origcons), 
            strlen(SCIPconsGetName(consdata->origcons))+1) );

      assert(SCIPgetCurrentNode(scip) == consdata->node || consdata->node == SCIPgetRootNode(scip));
      assert((SCIPgetNNodesLeft(scip)+SCIPgetNNodes(scip) == 1) == (consdata->node == SCIPgetRootNode(scip)));
      assert(SCIPnodeGetDepth(GCGconsOrigbranchGetNode(consdata->origcons)) == SCIPnodeGetDepth(consdata->node));
      assert(consdata->parentcons != NULL || SCIPnodeGetDepth(consdata->node) == 0);
      assert(consdata->parentcons == NULL || 
         SCIPconsGetData(consdata->parentcons)->origcons == GCGconsOrigbranchGetParentcons(consdata->origcons));

      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->nboundchangestreated, conshdlrData->nstack+1) );
      
      domchg = SCIPnodeGetDomchg(GCGconsOrigbranchGetNode(origcons));
      consdata->nboundchanges = SCIPdomchgGetNBoundchgs(domchg);
      consdata->nboundchangestreated[conshdlrData->nstack] = consdata->nboundchanges;

      if( consdata->nboundchanges > 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->boundchgvars, consdata->nboundchanges) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->boundtypes, consdata->nboundchanges) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->newbounds, consdata->nboundchanges) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->oldbounds, consdata->nboundchanges) );
      }

      consdata->nbranchingchanges = 0;

      for( i = 0; i < consdata->nboundchanges; i++ )
      {
         boundchg = SCIPdomchgGetBoundchg(domchg, i);

         consdata->boundchgvars[i] = SCIPboundchgGetVar(boundchg);
         consdata->newbounds[i] = SCIPboundchgGetNewbound(boundchg);
         consdata->boundtypes[i] = SCIPboundchgGetBoundtype(boundchg);

         if( SCIPboundchgGetBoundchgtype(boundchg) == SCIP_BOUNDCHGTYPE_BRANCHING )
         {
            consdata->nbranchingchanges++;
            assert(consdata->nbranchingchanges == i+1);
         }
      }

      consdata->created = TRUE;
      consdata->needprop = TRUE;

      assert((consdata->parentcons == NULL) == (conshdlrData->nstack == 0));
      if( consdata->parentcons != NULL )
      {
         parentdata = SCIPconsGetData(consdata->parentcons);
      
         assert(consdata->parentcons == conshdlrData->stack[conshdlrData->nstack-1]);
         assert(SCIPconsGetData(conshdlrData->stack[0])->parentcons == NULL);

         /* check whether bound changes were added in nodes on the path to the current node after first activation */
         for( i = 1; i < conshdlrData->nstack; i++ )
         {
            stackconsdata = SCIPconsGetData(conshdlrData->stack[i]);
            domchg = SCIPnodeGetDomchg(GCGconsOrigbranchGetNode(stackconsdata->origcons));

            assert(SCIPdomchgGetNBoundchgs(domchg) >= parentdata->nboundchangestreated[i]);

            if( SCIPdomchgGetNBoundchgs(domchg) != parentdata->nboundchangestreated[i] )
            {
               int diff;

               diff = SCIPdomchgGetNBoundchgs(domchg) - parentdata->nboundchangestreated[i];

               SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->boundchgvars, consdata->nboundchanges + diff) );
               SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->boundtypes, consdata->nboundchanges + diff) );
               SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->newbounds, consdata->nboundchanges + diff) );
               SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->oldbounds, consdata->nboundchanges + diff) );
         
               /* add bound changes to the boundchanges array */
               for( j = 0; j < SCIPdomchgGetNBoundchgs(domchg); j++ )
               {
                  boundchg = SCIPdomchgGetBoundchg(domchg, j);
                  if( j < stackconsdata->nboundchangestreated[i] )
                  { 
                     assert(stackconsdata->boundchgvars[j] == SCIPboundchgGetVar(boundchg)
                        && stackconsdata->newbounds[j] == SCIPboundchgGetNewbound(boundchg)
                        && stackconsdata->boundtypes[j] == SCIPboundchgGetBoundtype(boundchg));
                     continue;
                  }
                  if( j < parentdata->nboundchangestreated[i] )
                     continue;
                  consdata->boundchgvars[consdata->nboundchanges + j - parentdata->nboundchangestreated[i]] 
                     = SCIPboundchgGetVar(boundchg);
                  consdata->newbounds[consdata->nboundchanges + j - parentdata->nboundchangestreated[i]] 
                     = SCIPboundchgGetNewbound(boundchg);
                  consdata->boundtypes[consdata->nboundchanges + j - parentdata->nboundchangestreated[i]] 
                     = SCIPboundchgGetBoundtype(boundchg);

               }
               consdata->nboundchanges += diff;
            }
            consdata->nboundchangestreated[i] = SCIPdomchgGetNBoundchgs(domchg);
         }
      }

   }

   /* the node has to be repropagated if new variables were created after the node was left the last time 
    * or if new bound changes on directly transferred variables were found */
   assert(GCGpricerGetNPricedvars(scip) >= consdata->propagatedvars);
   if( GCGpricerGetNPricedvars(scip) > consdata->propagatedvars || GCGconsOrigbranchGetNPropBoundChgs(origscip, consdata->origcons) > 0 )
   {
      consdata->needprop = TRUE;
      SCIP_CALL( SCIPrepropagateNode(scip, consdata->node) );
   }

   if( consdata->nboundchanges - consdata->nboundchangestreated[conshdlrData->nstack] > 0 )
      SCIPdebugMessage("added %d boundchanges from previous nodes!\n", consdata->nboundchanges - consdata->nboundchangestreated[conshdlrData->nstack]);

   /* put constraint on the stack */
   if( conshdlrData->nstack >= conshdlrData->maxstacksize )
   {
      SCIPreallocMemoryArray(scip, &(conshdlrData->stack), 2*(conshdlrData->maxstacksize));
      conshdlrData->maxstacksize = 2*(conshdlrData->maxstacksize);
      SCIPdebugMessage("reallocating Memory for stack! %d --> %d\n", conshdlrData->maxstacksize/2, conshdlrData->maxstacksize);
   }
   conshdlrData->stack[conshdlrData->nstack] = cons;
   (conshdlrData->nstack)++;

   SCIPdebugMessage("Activating masterbranch constraint: <%s> [stack size: %d], needprop = %d.\n", 
      consdata->name, conshdlrData->nstack, consdata->needprop);

   /* apply global bound changes in the original problem to the master problem */
   if( !conshdlrData->pendingbndsactivated )
   {
      assert(conshdlrData->npendingbnds > 0);
      for( i = 0; i < conshdlrData->npendingbnds; i++ )
      {
         vardata = SCIPvarGetData(conshdlrData->pendingvars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_MASTER || vardata->vartype == GCG_VARTYPE_PRICING);
         
         if( vardata->vartype == GCG_VARTYPE_MASTER )
         {
            if( conshdlrData->pendingbndtypes[i] == SCIP_BOUNDTYPE_LOWER )
            {
               if( !SCIPisEQ(scip, SCIPvarGetLbGlobal(conshdlrData->pendingvars[i]), conshdlrData->pendingnewbnds[i]) )
               {
                  SCIPdebugMessage("Global lower bound of var <%s> set to %g\n", SCIPvarGetName(conshdlrData->pendingvars[i]),
                     conshdlrData->pendingnewbnds[i]);
                  SCIP_CALL( SCIPchgVarLbGlobal(scip, conshdlrData->pendingvars[i], conshdlrData->pendingnewbnds[i]) );
               }
            }
            else
            {
               if( !SCIPisEQ(scip, SCIPvarGetUbGlobal(conshdlrData->pendingvars[i]), conshdlrData->pendingnewbnds[i]) )
               {
                  SCIPdebugMessage("Global upper bound of var <%s> set to %g\n", SCIPvarGetName(conshdlrData->pendingvars[i]),
                     conshdlrData->pendingnewbnds[i]);
                  SCIP_CALL( SCIPchgVarUbGlobal(scip, conshdlrData->pendingvars[i], conshdlrData->pendingnewbnds[i]) );
               }
            }
         }
         else
         {
            /* this is a global boundchange on a variable that belongs to a block, 
             * we have to adjust the bound of the corresponding variable in the pricing problem */
            if( conshdlrData->pendingbndtypes[i] == SCIP_BOUNDTYPE_LOWER )
            {
               SCIP_CALL( SCIPchgVarLb(GCGrelaxGetPricingprob(origscip, vardata->blocknr), 
                     conshdlrData->pendingvars[i], conshdlrData->pendingnewbnds[i]) );
            }
            else
            {
               SCIP_CALL( SCIPchgVarUbGlobal(GCGrelaxGetPricingprob(origscip, vardata->blocknr), 
                     conshdlrData->pendingvars[i], conshdlrData->pendingnewbnds[i]) );
            }
         }
      }
      conshdlrData->pendingbndsactivated = TRUE;
   }
   

   /* apply local bound changes in the original problem to the pricing problems */
   for( i = 0; i < consdata->nboundchanges; i++ )
   {
         vardata = SCIPvarGetData(consdata->boundchgvars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         assert(vardata->blocknr >= -1 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
         assert(vardata->blocknr == -1 || vardata->data.origvardata.pricingvar != NULL);


         /* if variable belongs to no block, skip it here because the bound changes are treated in the propagation */
         if( vardata->blocknr == -1 )
            continue;

         assert(vardata->data.origvardata.pricingvar != NULL);

         /* set corresponding bound in the pricing problem */
         /* lower bound was changed */
         if( consdata->boundtypes[i] == SCIP_BOUNDTYPE_LOWER )
         {
            consdata->oldbounds[i] = SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar);

            SCIP_CALL( SCIPchgVarLb(GCGrelaxGetPricingprob(origscip, vardata->blocknr),
                  vardata->data.origvardata.pricingvar, consdata->newbounds[i]) );
            SCIPdebugMessage("tightened lower bound of var %s from %g to %g\n", 
               SCIPvarGetName(vardata->data.origvardata.pricingvar), consdata->oldbounds[i],
               consdata->newbounds[i]);
         }
         /* upper bound was changed */
         else
         {
            consdata->oldbounds[i] = SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar);

            SCIP_CALL( SCIPchgVarUb(GCGrelaxGetPricingprob(origscip, vardata->blocknr),
                  vardata->data.origvardata.pricingvar, consdata->newbounds[i]) );
            SCIPdebugMessage("tightened upper bound of var %s from %g to %g\n", 
               SCIPvarGetName(vardata->data.origvardata.pricingvar), consdata->oldbounds[i],
               consdata->newbounds[i]);
         }
   }

   /* call branching specific activation method */
   if( consdata->branchrule != NULL )
   {
      SCIP_CALL( GCGrelaxBranchActiveMaster(origscip, consdata->branchrule, consdata->branchdata) );
   }

   return SCIP_OKAY;
}



/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA* consdata;
   SCIP_VARDATA* vardata;
   SCIP* origscip;
   int i;

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

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      consdata->propagatedvars = GCGpricerGetNPricedvars(scip);

   /* remove constraint from the stack */
   (conshdlrData->nstack)--;


   SCIPdebugMessage("Deactivating masterbranch constraint: <%s> [stack size: %d].\n", 
      consdata->name, conshdlrData->nstack);

   /* undo local bound changes in the original problem to the pricing problems */
   for( i = consdata->nboundchanges - 1; i >= 0; i-- )
   {
         vardata = SCIPvarGetData(consdata->boundchgvars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         assert(vardata->blocknr >= -1 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
         assert(vardata->blocknr == -1 || vardata->data.origvardata.pricingvar != NULL);

         /* if variable belongs to no block, local bound in master was set, is reset automatically */
         if( vardata->blocknr == -1 )
            continue;

         assert(vardata->data.origvardata.pricingvar != NULL);
         assert(GCGrelaxGetPricingprob(origscip, vardata->blocknr) != NULL);

         /* reset corresponding bound in the pricing problem */
         /* lower bound was changed */
         if( consdata->boundtypes[i] == SCIP_BOUNDTYPE_LOWER )
         {
            assert(SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar) == consdata->newbounds[i]);
            if( SCIPvarGetLbGlobal(consdata->boundchgvars[i]) == consdata->newbounds[i] )
               continue;

            SCIP_CALL( SCIPchgVarLb(GCGrelaxGetPricingprob(origscip, vardata->blocknr),
                  vardata->data.origvardata.pricingvar, consdata->oldbounds[i]) );
            SCIPdebugMessage("relaxed lower bound of var %s from %g to %g\n", 
               SCIPvarGetName(vardata->data.origvardata.pricingvar), consdata->newbounds[i],
               consdata->oldbounds[i]);
         }
         /* upper bound was changed */
         else
         {
            assert(SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar) == consdata->newbounds[i]);
            if( SCIPvarGetUbGlobal(consdata->boundchgvars[i]) == consdata->newbounds[i] )
               continue;

            SCIP_CALL( SCIPchgVarUb(GCGrelaxGetPricingprob(origscip, vardata->blocknr),
                  vardata->data.origvardata.pricingvar, consdata->oldbounds[i]) );
            SCIPdebugMessage("relaxed upper bound of var %s from %g to %g\n", 
               SCIPvarGetName(vardata->data.origvardata.pricingvar), consdata->newbounds[i],
               consdata->oldbounds[i]);
         }
   }

   /* call branching specific deactivation method */
   if( consdata->branchrule != NULL )
   {
      SCIP_CALL( GCGrelaxBranchDeactiveMaster(GCGpricerGetOrigprob(scip), consdata->branchrule, consdata->branchdata) );
   }

   return SCIP_OKAY;
}



/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropMasterbranch)
{
   SCIP* origscip;
   SCIP_CONSHDLRDATA* conshdlrData;  
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_VARDATA* vardata;
   SCIP_Real val;

   int propcount;
   int i;
   int j;
   int k;
   SCIP_Bool fixed;

   SCIP_VARDATA* boundchgvardata;
   SCIP_VAR** vars;
   int nvars;
   int nboundchanges;

   SCIP_VAR**         propvars;              /**< original variable for which the propagation found domain reductions */
   SCIP_BOUNDTYPE*    propboundtypes;        /**< type of the domain new bound found by propagation */
   SCIP_Real*         propbounds;            /**< new lower/upper bound of the propagated original variable */
   int                npropbounds;
   SCIP_VAR* mastervar;


   assert(conshdlr != NULL); 
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   *result = SCIP_DIDNOTRUN;

   /* the constraint data of the cons related to the current node */
   cons = conshdlrData->stack[conshdlrData->nstack-1];
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !consdata->needprop && GCGconsOrigbranchGetNPropBoundChgs(origscip, consdata->origcons) == 0)
   {
#ifdef CHECKPROPAGATEDVARS
      SCIP_Bool consistent;
      consistent = checkVars(scip, conshdlr, TRUE);
      assert(consistent);
#endif

      SCIPdebugMessage("No propagation of masterbranch constraint needed: <%s>, stack size = %d.\n", 
         consdata->name, conshdlrData->nstack);
      
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   vars = GCGpricerGetPricedvars(scip);
   nvars = GCGpricerGetNPricedvars(scip);

   SCIPdebugMessage("Starting propagation of masterbranch constraint: <%s>, stack size = %d, newvars = %d, npendingbnds = %d, npropbounds = %d.\n", 
      consdata->name, conshdlrData->nstack, nvars - consdata->propagatedvars, conshdlrData->npendingbnds, GCGconsOrigbranchGetNPropBoundChgs(origscip, consdata->origcons));

   *result = SCIP_DIDNOTFIND;

   propcount = 0;


   /* propagate all bound changes or only the branching bound changes, depending on the setting for the enforcement of proper variables */
   nboundchanges = (conshdlrData->enforceproper ? consdata->nboundchanges : consdata->nbranchingchanges);

   assert((conshdlrData->npendingbnds > 0) || conshdlrData->pendingbndsactivated);
   
   /* iterate over all master variables and apply global bound changes */
   if( conshdlrData->npendingbnds > 0 && conshdlrData->pendingbndsactivated )
   {
      for( i = 0; i < nvars; i++)
      {
         vardata = SCIPvarGetData(vars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_MASTER);
         assert(vardata->blocknr >= -1 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
         assert(vardata->data.mastervardata.norigvars >= 0);
         assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
         assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);
         assert(vardata->blocknr != -1 || vardata->data.mastervardata.norigvars == 2 );

         fixed = FALSE;

         /* only look at master variables not globally fixed to zero that belong to a block */
         if( !SCIPisFeasZero(scip, SCIPvarGetUbGlobal(vars[i])) && vardata->blocknr != -1 )
         {
            /* iterate over global bound changes that were not yet checked for the master variables */
            for( k = 0; k < conshdlrData->npendingbnds && !fixed; k++ )
            {
               assert(conshdlrData->pendingvars[k] != NULL);
               boundchgvardata = SCIPvarGetData(conshdlrData->pendingvars[k]);
               assert(boundchgvardata != NULL);
               assert(boundchgvardata->vartype != GCG_VARTYPE_ORIGINAL);
               assert(boundchgvardata->blocknr >= -1 && boundchgvardata->blocknr < GCGrelaxGetNPricingprobs(origscip));

               /* the boundchage was performed on a variable in another block, continue */
               if( boundchgvardata->blocknr != vardata->blocknr )
                  continue;

               assert(boundchgvardata->blocknr != -1);
               assert(boundchgvardata->data.pricingvardata.origvars[0] != NULL);

               /* val is the value of the branching variable in the current mastervar,
                * we set it to 0.0, since variables with 0 coefficient are not stored in the origvars array,
                * if we do not find the branching variable in this array, it has value 0.0 */
               val = 0.0;

               for( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
               {
                  assert(SCIPvarGetData(vardata->data.mastervardata.origvars[j])->blocknr == vardata->blocknr);
                  
                  /* check whether the original variable contained in the master variable equals the variable 
                   * on which the current branching was performed */
                  if( vardata->data.mastervardata.origvars[j] == boundchgvardata->data.pricingvardata.origvars[0] )
                  {
                     val = vardata->data.mastervardata.origvals[j];
                     break;
                  }
               }

               /* if the variable contains a part of the branching variable that violates the bound, 
                * fix the master variable to 0 */
               
               /* branching imposes new lower bound */
               if( conshdlrData->pendingbndtypes[k] == SCIP_BOUNDTYPE_LOWER && 
                  SCIPisFeasLT(scip, val, conshdlrData->pendingnewbnds[k]) )
               {
                  SCIPchgVarUbGlobal(scip, vars[i], 0.0);
                  propcount++;
                  fixed = TRUE;
                  break;
               }
               /* branching imposes new upper bound */
               if( conshdlrData->pendingbndtypes[k] == SCIP_BOUNDTYPE_UPPER && 
                  SCIPisFeasGT(scip, val, conshdlrData->pendingnewbnds[k]) )
               {
                  SCIPchgVarUbGlobal(scip, vars[i], 0.0);
                  propcount++;
                  fixed = TRUE;
                  break;
               }
            }
         }
      }
      conshdlrData->npendingbnds = 0;

      SCIPdebugMessage("Finished handling of pending global bound changes: %d changed bounds\n", propcount);
   }

   /* iterate over all master variables created after the current node was left the last time */
   for( i = consdata->propagatedvars; i < nvars; i++)
   {
      vardata = SCIPvarGetData(vars[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_MASTER);
      assert(vardata->blocknr >= -1 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
      assert(vardata->data.mastervardata.norigvars >= 0);
      assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
      assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);
      assert(vardata->blocknr != -1 || vardata->data.mastervardata.norigvars == 2 );

      fixed = FALSE;

      /* only look at variables not already fixed to 0 or that belong to no block */
      if( (!SCIPisFeasZero(scip, SCIPvarGetUbLocal(vars[i]))) || vardata->blocknr == -1 )
      {
         /* the variable was copied from original to master */
         if( vardata->blocknr == -1 )
         {
            /* iterate over bound changes performed at the current node's equivalent in the original tree */
            for( k = 0; k < nboundchanges; k++ )
            {
               assert(SCIPisFeasEQ(scip, vardata->data.mastervardata.origvals[0], 1.0));
               assert(SCIPisFeasEQ(scip, vardata->data.mastervardata.origvals[1], 0.0));

               if( vardata->data.mastervardata.origvars[0] == consdata->boundchgvars[k] )
               {
                  /* branching imposes new lower bound */
                  if( consdata->boundtypes[k] == SCIP_BOUNDTYPE_LOWER 
                     && SCIPisGT(scip, consdata->newbounds[k], SCIPvarGetLbLocal(vars[i])) )
                  {
                     SCIPchgVarLb(scip, vars[i], consdata->newbounds[k]);
                     propcount++;
                  }
                  /* branching imposes new upper bound */
                  if( consdata->boundtypes[k] == SCIP_BOUNDTYPE_UPPER
                     && SCIPisLT(scip, consdata->newbounds[k], SCIPvarGetUbLocal(vars[i])) )
                  {
                     SCIPchgVarUb(scip, vars[i], consdata->newbounds[k]);
                     propcount++;
                  }

               }
            }
         }
         else
         {
            /* iterate over bound changes performed at the current node's equivalent in the original tree */
            for( k = 0; k < nboundchanges && !fixed; k++ )
            {
               boundchgvardata = SCIPvarGetData(consdata->boundchgvars[k]);
               assert(boundchgvardata != NULL);
               assert(boundchgvardata->vartype == GCG_VARTYPE_ORIGINAL);
               assert(boundchgvardata->blocknr >= -1 && boundchgvardata->blocknr < GCGrelaxGetNPricingprobs(origscip));

               /* the boundchage was performed on a variable in another block, continue */
               if( boundchgvardata->blocknr != vardata->blocknr )
                  continue;

               assert(boundchgvardata->blocknr != -1);
               
               /* val is the value of the branching variable in the current mastervar,
                * we set it to 0.0, since variables with 0 coefficient are not stored in the origvars array,
                * if we do not find the branching variable in this array, it has value 0.0 */
               val = 0.0;

               /* iterate over all original variables contained in the current master variable */
               for( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
               {
                  assert(SCIPvarGetData(vardata->data.mastervardata.origvars[j])->blocknr == vardata->blocknr);
                  
                  /* check whether the original variable contained in the master variable equals the variable 
                   * on which the current branching was performed */
                  if( vardata->data.mastervardata.origvars[j] == consdata->boundchgvars[k] )
                  {
                     val = vardata->data.mastervardata.origvals[j];
                     break;
                  }
               }

               /* if the variable contains a part of the branching variable that violates the bound, 
                * fix the master variable to 0 */
               
               /* branching imposes new lower bound */
               if( consdata->boundtypes[k] == SCIP_BOUNDTYPE_LOWER && 
                  SCIPisFeasLT(scip, val, consdata->newbounds[k]) )
               {
                  SCIPchgVarUb(scip, vars[i], 0.0);
                  propcount++;
                  fixed = TRUE;
                  break;
               }
               /* branching imposes new upper bound */
               if( consdata->boundtypes[k] == SCIP_BOUNDTYPE_UPPER && 
                  SCIPisFeasGT(scip, val, consdata->newbounds[k]) )
               {
                  SCIPchgVarUb(scip, vars[i], 0.0);
                  propcount++;
                  fixed = TRUE;
                  break;
               }
            }
         }
      }
   }
   SCIPdebugMessage("Finished propagation of newly created variables: %d changed bounds\n", propcount);

   /* get local bound changes on variables directly transferred to the master problem and apply them */
   SCIP_CALL( GCGconsOrigbranchGetPropBoundChgs(origscip, consdata->origcons, &propvars, &propboundtypes,
         &propbounds, &npropbounds) );
   for( i = 0; i < npropbounds; i++ )
   {
      vardata = SCIPvarGetData(propvars[i]);
      assert(vardata != NULL);
      assert(vardata->blocknr == -1);
      assert(vardata->data.origvardata.nmastervars == 1);
      assert(vardata->data.origvardata.mastervals[0] == 1);
      assert(vardata->data.origvardata.mastervars[0] != NULL);

      mastervar = vardata->data.origvardata.mastervars[0];

      if( propboundtypes[i] == SCIP_BOUNDTYPE_LOWER )
      {
         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(mastervar), propbounds[i]) )
         {
            SCIP_CALL( SCIPchgVarLb(scip, mastervar, propbounds[i]) );
            propcount++;
            SCIPdebugMessage("changed lb of var %s locally to %g\n", SCIPvarGetName(propvars[i]), propbounds[i]);
         }
      }
      else 
      {
         if( !SCIPisEQ(scip, SCIPvarGetUbLocal(mastervar), propbounds[i]) )
         {
            SCIP_CALL( SCIPchgVarUb(scip, mastervar, propbounds[i]) );
            propcount++;
            SCIPdebugMessage("changed ub of var %s locally to %g\n", SCIPvarGetName(propvars[i]), propbounds[i]);
         }
      }
   }

   SCIPdebugMessage("Finished propagation of %d stored propagated bounds: %d vars fixed.\n", npropbounds, propcount);

   /* call branching rule specific propagation method */     
   if( consdata->branchrule != NULL )
   {
      /* TODO: count number of propagations */
      SCIP_CALL( GCGrelaxBranchPropMaster(GCGpricerGetOrigprob(scip), consdata->branchrule, consdata->branchdata, result) );
   }

   /*SCIPdebugMessage("Finished external propagation: %d vars fixed.\n", propcount);*/

   if( *result != SCIP_CUTOFF )
      if( propcount > 0 )
         *result = SCIP_REDUCEDDOM;

   consdata->needprop = FALSE;
   consdata->propagatedvars = GCGpricerGetNPricedvars(scip);

#ifdef CHECKPROPAGATEDVARS
   {
      SCIP_Bool consistent;
      consistent = checkVars(scip, conshdlr, TRUE);
      assert(consistent);
   }
#endif

   return SCIP_OKAY;
}


static
SCIP_DECL_CONSENFOLP(consEnfolpMasterbranch)
{
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSENFOPS(consEnfopsMasterbranch)
{
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSCHECK(consCheckMasterbranch)
{
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSLOCK(consLockMasterbranch)
{
   return SCIP_OKAY;
}


/* define not used callbacks as NULL */
#define conshdlrCopyMasterbranch NULL
#define consPresolMasterbranch NULL
#define consRespropMasterbranch NULL
#define consExitsolMasterbranch NULL
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
 * Callback methods of event handler
 */

#define eventCopyOrigvarbound NULL
#define eventFreeOrigvarbound NULL
#define eventExitOrigvarbound NULL

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitOrigvarbound)
{  
   return SCIP_OKAY;
}



/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolOrigvarbound)
{  
   SCIP_VAR** vars;
   int nvars;
   int i;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, vars[i], SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_BOUNDCHANGED, 
            eventhdlr, NULL, NULL) );
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_EVENTEXITSOL(eventExitsolOrigvarbound)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of origvarbound event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventExitsolOrigvarbound NULL
#endif

/** frees specific event data */
#if 0
static
SCIP_DECL_EVENTDELETE(eventDeleteOrigvarbound)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of origvarbound event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventDeleteOrigvarbound NULL
#endif

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecOrigvarbound)
{  
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR* var;
   SCIP_Real oldbound;
   SCIP_Real newbound;

   SCIP_VARDATA* vardata;

   eventtype = SCIPeventGetType(event);
   var = SCIPeventGetVar(event);
   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);

   SCIPdebugMessage("eventexec: eventtype = %d, var = %s, oldbound = %f, newbound = %f\n", eventtype, SCIPvarGetName(var), oldbound, newbound);
   //printf("eventexec: eventtype = %d, var = %s, oldbound = %f, newbound = %f, diff = %g\n", eventtype, SCIPvarGetName(var), oldbound, newbound, oldbound-newbound);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   
   if( vardata->blocknr != -1 && GCGrelaxIsPricingprobRelevant(scip, vardata->blocknr) )
   {
      if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
      {
         SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip), 
               vardata->data.origvardata.pricingvar, SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );
      }
      if( (eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0 )
      {
         SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip),
               vardata->data.origvardata.pricingvar, SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );
      }
   }

   if( vardata->blocknr == -1 && SCIPgetStage(GCGrelaxGetMasterprob(scip)) >= SCIP_STAGE_SOLVING )
   {
      assert(vardata->data.origvardata.nmastervars == 1);
      assert(vardata->data.origvardata.mastervals[0] == 1);
      assert(vardata->data.origvardata.mastervars[0] != NULL);
      
      if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
      {
         assert(SCIPvarGetLbGlobal(vardata->data.origvardata.mastervars[0]) == oldbound);
         SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip), 
               vardata->data.origvardata.mastervars[0], SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );
         //printf("-> saved change of lb of var %s to %g\n", SCIPvarGetName(vardata->data.origvardata.mastervars[0]), newbound);
      }
      if( (eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0 )
      {
         assert(SCIPvarGetUbGlobal(vardata->data.origvardata.mastervars[0]) == oldbound);
         SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip), 
               vardata->data.origvardata.mastervars[0], SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );
         //printf("-> saved change of ub of var %s to %g\n", SCIPvarGetName(vardata->data.origvardata.mastervars[0]), newbound);
      }
      if( (eventtype & SCIP_EVENTTYPE_LBTIGHTENED) != 0 )
      {
         SCIP_CALL( GCGconsOrigbranchAddPropBoundChg(scip, GCGconsOrigbranchGetActiveCons(scip), var,
               SCIP_BOUNDTYPE_LOWER, newbound) );
      }
      if( (eventtype & SCIP_EVENTTYPE_UBTIGHTENED) != 0 )
      {
         SCIP_CALL( GCGconsOrigbranchAddPropBoundChg(scip, GCGconsOrigbranchGetActiveCons(scip), var,
               SCIP_BOUNDTYPE_UPPER, newbound) );
      }
      
   }

   return SCIP_OKAY;
}


/*
 * interface methods
 */


/** creates the handler for masterbranch constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrMasterbranch(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

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
         conshdlrCopyMasterbranch, consFreeMasterbranch, consInitMasterbranch, consExitMasterbranch,
         consInitpreMasterbranch, consExitpreMasterbranch, consInitsolMasterbranch, consExitsolMasterbranch,
         consDeleteMasterbranch, consTransMasterbranch, consInitlpMasterbranch,
         consSepalpMasterbranch, consSepasolMasterbranch, consEnfolpMasterbranch, consEnfopsMasterbranch, consCheckMasterbranch,
         consPropMasterbranch, consPresolMasterbranch, consRespropMasterbranch, consLockMasterbranch,
         consActiveMasterbranch, consDeactiveMasterbranch,
         consEnableMasterbranch, consDisableMasterbranch,
         consPrintMasterbranch, consCopyMasterbranch, consParseMasterbranch, 
         conshdlrData) );

   /* create event handler data */
   eventhdlrdata = NULL;

   /* include event handler into original SCIP */
   SCIP_CALL( SCIPincludeEventhdlr(GCGpricerGetOrigprob(scip), EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventCopyOrigvarbound, eventFreeOrigvarbound, eventInitOrigvarbound, eventExitOrigvarbound, 
         eventInitsolOrigvarbound, eventExitsolOrigvarbound, eventDeleteOrigvarbound, eventExecOrigvarbound,
         eventhdlrdata) );

   SCIP_CALL( SCIPaddBoolParam(GCGpricerGetOrigprob(scip), "relaxing/gcg/enforceproper",
         "should propagated bound changes in the original be enforced in the master (only proper vars)?",
         &conshdlrData->enforceproper, FALSE, TRUE, NULL, NULL) );



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
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->propagatedvars = 0;
   consdata->needprop = TRUE;

   consdata->node = node;
   consdata->parentcons = parentcons;
   consdata->child1cons = NULL;
   consdata->child2cons = NULL;
   consdata->probingtmpcons = NULL;
   consdata->created = FALSE;
   consdata->origcons = NULL;
   consdata->name = NULL;

   consdata->branchrule = NULL;
   consdata->branchdata = NULL;

   consdata->boundchgvars = NULL;
   consdata->boundtypes = NULL;
   consdata->newbounds = NULL;
   consdata->oldbounds = NULL;
   consdata->nboundchangestreated = NULL;
   consdata->nboundchanges = 0;
   consdata->nactivated = 0;
   

   SCIPdebugMessage("Creating masterbranch constraint.\n");

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, "masterbranch", conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
         TRUE, FALSE, FALSE, FALSE, TRUE) );

   if( parentcons != NULL )
   {      
      SCIP_CONSDATA* parentdata;

      parentdata = SCIPconsGetData(parentcons);
      assert(parentdata != NULL);

      if( parentdata->child1cons == NULL )
      {
         parentdata->child1cons = *cons;
      }
      else
      {
         assert(parentdata->child2cons == NULL || SCIPinProbing(scip));

         /* store the second child in case we are in probing and have to overwrite it */
         if( SCIPinProbing(scip) )
         {
            assert(parentdata->probingtmpcons == NULL);
            parentdata->probingtmpcons = parentdata->child2cons;
         }

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
   if( conshdlr == NULL )
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
   if( conshdlr == NULL )
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

/** returns the number of elements on the stack */
int GCGconsMasterbranchGetNStackelements(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return -1;
   }   
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   return conshdlrData->nstack;
}

/** returns the branching data for a given masterbranch constraint */
GCG_BRANCHDATA* GCGconsMasterbranchGetBranchdata(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->branchdata;
}

/** returns the node in the B&B tree at which the given masterbranch constraint is sticking */
SCIP_NODE* GCGconsMasterbranchGetNode(
   SCIP_CONS*            cons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

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
   assert(consdata != NULL);

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
   assert(consdata != NULL);

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
   assert(consdata != NULL);

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
   assert(consdata != NULL);

   return consdata->origcons;
}

/** sets the origbranch constraint of the node in the master program corresponding to the node 
    at which the given masterbranchbranch constraint is sticking */
void GCGconsMasterbranchSetOrigcons(
   SCIP_CONS*            cons,
   SCIP_CONS*            origcons
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->origcons == NULL || origcons == NULL);

   consdata->origcons = origcons;
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

   if( scip == NULL )
      return;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      assert(0);
      return;
   }   

   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);

   for( i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      assert(consdata->node != NULL);
      assert((consdata->parentcons == NULL) == (SCIPnodeGetDepth(consdata->node) == 0));
      assert(consdata->origcons == NULL || consdata->created);
      assert(consdata->parentcons == NULL || SCIPconsGetData(consdata->parentcons)->child1cons == conss[i]
         || SCIPconsGetData(consdata->parentcons)->child2cons == conss[i]
         || ( SCIPinProbing(scip) && SCIPconsGetData(consdata->parentcons)->probingtmpcons == conss[i]));
      assert(consdata->child1cons == NULL || SCIPconsGetData(consdata->child1cons)->parentcons == conss[i]);
      assert(consdata->child2cons == NULL || SCIPconsGetData(consdata->child2cons)->parentcons == conss[i]);
      assert(consdata->probingtmpcons == NULL || SCIPinProbing(scip));
      assert(consdata->probingtmpcons == NULL || SCIPconsGetData(consdata->probingtmpcons)->parentcons == conss[i]);
      assert(consdata->origcons == NULL || 
         GCGconsOrigbranchGetMastercons(consdata->origcons) == conss[i]);
   }

   SCIPdebugMessage("checked consistency of %d masterbranch constraints, all ok!\n", nconss);
}
