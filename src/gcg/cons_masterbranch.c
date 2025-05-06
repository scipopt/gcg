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

/**@file   cons_masterbranch.c
 * 
 * @brief  constraint handler for storing the branching decisions at each node of the tree
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @author Christian Puchert
 * @author Marcel Schmickerath
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#define MASTERSEP_DEBUG
#include <assert.h>
#include <string.h>

#include "gcg/gcg.h"
#include "gcg/branch_generic.h"
#include "gcg/cons_masterbranch.h"
#include "gcg/cons_origbranch.h"
#include "gcg/relax_gcg.h"
#include "gcg/pricer_gcg.h"
#include "gcg/pub_colpool.h"
#include "gcg/gcgvarhistory.h"
#include "gcg/event_sepacuts.h"
#include "gcg/mastersepacut.h"
#include "gcg/struct_sepagcg.h"


/*#define CHECKPROPAGATEDVARS*/

/* constraint handler properties */
#define CONSHDLR_NAME          "masterbranch"
#define CONSHDLR_DESC          "store branching decision at nodes of the tree constraint handler"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         * propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */
#define CONSHDLR_PROPTIMING SCIP_PROPTIMING_ALWAYS

#define EVENTHDLR_NAME         "origvarbound"
#define EVENTHDLR_DESC         "event handler for bound changes on original variables"


/** constraint data for masterbranch constraints */
struct SCIP_ConsData
{
   char*                 name;               /**< name of the constraint */
   int                   npropvars;          /**< number of variables that existed the last time the related node was propagated,
                                                  used to determine whether the constraint should be repropagated */
   SCIP_Bool             needprop;           /**< should the constraint be propagated? */
   SCIP_NODE*            node;               /**< the node at which the constraint is sticking */
   int                   nactivated;         /**< number of times the constraint has been activated so far */

   SCIP_CONS*            parentcons;         /**< the masterbranch constraint of the parent node */
   SCIP_CONS**           childconss;         /**< array of the masterbranch constraints of child nodes */
   int                   nchildconss;        /**< number of the masterbranch constraints of child nodes */
   int                   maxchildconss;     /**< capacity of childconss array */
   SCIP_CONS*            probingtmpcons;     /**< pointer to save the last child in the childconss array if it is overwritten in probing mode */
   SCIP_CONS*            origcons;           /**< the corresponding origbranch cons in the original program */

   GCG_BRANCHDATA*       branchdata;         /**< branching data stored by the branching rule at the corresponding origcons constraint
                                              *   containing information about the branching restrictions */
   SCIP_BRANCHRULE*      branchrule;         /**< branching rule that created the corresponding node in the original problem and imposed
                                              *   branching restrictions */

   /* pointer to the last variable that we have seen, any newer variables are unseen */
   GCG_VARHISTORY*       knownvarhistory;    /**< pointer to the history of priced variables */

   /* local bound changes on original variables that belong to a unique block */
   SCIP_VAR**            localbndvars;       /**< original variables of bound changes stored at the current node */
   SCIP_BOUNDTYPE*       localbndtypes;      /**< types of the bound changes stored at the current node */
   SCIP_Real*            localnewbnds;       /**< new bounds for the bound changes stored at the current node */
   SCIP_Real*            localoldbnds;       /**< old bounds for the bound changes stored at the current node */

   int*                  nlocalbndchgstreated; /**< number of bound changes of the nodes on the way from the current node to
                                                *   the root node that have been treated so far */
   int                   maxlocalbndchgstreated; /**< capacity of nlocalbndchgstreated */
   int                   nlocalbndchgs;      /**< number of bound changes */
   int                   maxlocalbndchgs;    /**< capacity of corresponding arrays */
   int                   nbranchingchgs;     /**< number of bound changes due to branching (<= nlocalbndchgs) */

   /* local bound changes on original variables that have been directly copied to the master problem */
   SCIP_VAR**            copiedvars;         /**< original variables on which local bounds were changed */
   GCG_BOUNDTYPE*        copiedvarbndtypes;  /**< types of the new local bounds of the coped original variables */
   SCIP_Real*            copiedvarbnds;      /**< new lower/upper bounds of the coped original variables */
   int                   ncopiedvarbnds;     /**< number of new local bounds stored */
   int                   maxcopiedvarbnds;   /**< size of copiedvars, copiedvarbndtypes, and copiedvarbnds arrays */

   /* constraints that enforce the branching restrictions on the original problem */
   SCIP_CONS**           origbranchconss;    /**< constraints in the original problem that enforce the branching decision */
   int                   norigbranchconss;   /**< number of constraints in the original problem that enforce the branching decision */
   int                   maxorigbranchconss;/**< capacity of origbranchconss array */

   /* information needed to update cuts generated by master separators */
   GCG_EXTENDEDMASTERCONSDATA** addedcuts;        /**< master separator cuts which were created in the node defined by this branch */
   int                        naddedcuts;         /**< number of master separator cuts which were added in this node */
   int                        firstnewcut;        /**< index of first cut from this node in active cuts */
   SCIP_Bool                  addedcutsinit;      /**< indicates if above data structures have been initialized */
   SCIP_Bool                  nodestoredcuts;     /**< indicates if any cuts were stored at this node */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   GCG*                  gcg;                /**< GCG data structure */
   /* active masterbranch constraints on the path from the root node to the current node */
   SCIP_CONS**           stack;              /**< stack for storing active constraints */
   int                   nstack;             /**< number of elements on the stack */
   int                   maxstacksize;       /**< maximum size of the stack */

   /* global bound changes on the original problem */
   SCIP_VAR**            pendingvars;        /**< pricing variables or master variable copies corresponding to pending bound changes (global bound changes) */
   SCIP_BOUNDTYPE*       pendingbndtypes;    /**< types of the pending bound changes (global bound changes) */
   SCIP_Real*            pendingnewbnds;     /**< new bounds corresponding to pending bound changes (global bound changes) */
   SCIP_HASHMAP*         pendingvarmaplb;    /**< maps vars to indices of the pendingvars array for global lower bound changes */
   SCIP_HASHMAP*         pendingvarmapub;    /**< maps vars to indices of the pendingvars array for global upper bound changes */
   int                   npendingbnds;       /**< number of pending bound changes (global bound changes) */
   SCIP_Bool             pendingbndsactivated; /**< were pending bound changes already activated? */
   int                   maxpendingbnds;     /**< size of the array corresponding to pending bound changes */
   SCIP_Bool             enforceproper;      /**< should proper variables be enforced? */

   /* information needed by applyLocalBndchgsToPricedMastervars */
   SCIP_VAR***           collectedbndvars;   /**< collected original variables of bound changes stored for each block */
   SCIP_Real**           collectedlbnds;     /**< collected upper bound changes for each block */
   SCIP_Real**           collectedubnds;     /**< collected upper bound changes for each block */
   int                   maxblocknum;        /**< maximal number of blocks that can be handled by the collected* arrays */
   int*                  ncollectedbndvars;  /**< number of collected bndvars per block */
   int*                  maxcollectedbndvars;/**< capacity of the collected* arrays per block */
   int**                 linkingvaridxs;     /**< position of linking vars in collectedbndvars */
   int                   maxlinkingvaridxs;  /**< capacity of linkingvaridxs */

   /* handler which manages the master separator cuts */
   SCIP_EVENTHDLR*       eventhdlr;
};

struct SCIP_EventhdlrData
{
   GCG*                  gcg;                /**< GCG data structure */
};

/*
 * Local methods
 */

/**< updates a master separator cut with all the variables it 'missed' while being inactive */
static
SCIP_RETCODE addMissedVariables(
   GCG*                          gcg,              /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*   mastersepacut     /**< master separator cut */
   )
{
   SCIP* scip;
   GCG_VARHISTORY*         varhistory;
   GCG_SEPA*               sepa;

   assert(gcg != NULL);
   assert(mastersepacut != NULL);

   scip = GCGgetMasterprob(gcg);
   sepa = GCGmastersepacutGetSeparator(GCGextendedmasterconsGetSepamastercut(mastersepacut));
   assert(sepa != NULL);
   varhistory = GCGmastersepacutGetVarHistory(GCGextendedmasterconsGetSepamastercut(mastersepacut));
   assert(varhistory != NULL);

   while( GCGvarhistoryHasNext(varhistory) )
   {
      SCIP_VAR*  mastervar;
      SCIP_VAR** origvars;
      SCIP_VAR** pricingvars;
      SCIP_Real* pricingvals;
      SCIP_Real  coef;
      int        blocknr;
      int        npricingvars;
      int        nnonzeros;
      int        j;

      SCIP_CALL( GCGvarhistoryNext(scip, &varhistory) );
      SCIP_CALL( GCGvarhistoryGetVar(varhistory, &mastervar) );

      assert(mastervar != NULL);
      if( SCIPvarIsDeleted(mastervar) )
      {
         SCIPdebugMessage("Skipping deleted Variable <%s>!\n", SCIPvarGetName(mastervar));
         continue;
      }
      assert(GCGvarIsMaster(mastervar));

      /* get the pricing variables corresponding to the original variables which define the master variable */
      npricingvars = GCGmasterVarGetNOrigvars(mastervar);
      origvars = GCGmasterVarGetOrigvars(mastervar);
      pricingvals = GCGmasterVarGetOrigvals(mastervar);
      SCIPallocBufferArray(scip, &pricingvars, npricingvars);
      nnonzeros = 0;
      for( j = 0; j < npricingvars; j++ )
      {
         pricingvars[j] = GCGoriginalVarGetPricingVar(origvars[j]);
         assert(pricingvars[j] != NULL);
         if( pricingvals[j] != 0.0 )
            nnonzeros++;
      }

      /* compute the coefficient for this master variable */
      coef = 0.0;
      if( (npricingvars > 0)  && (nnonzeros > 0) )
      {
         blocknr = GCGvarGetBlock(pricingvars[0]);
         SCIP_CALL( sepa->gcgsepagetvarcoefficient(gcg, sepa, mastersepacut, pricingvars, pricingvals, npricingvars,
                                                   blocknr, &coef) );
      }

      /* add variable with its coefficient to the cut */
      if( !SCIPisZero(scip, coef) )
      {
         SCIP_ROW* mastercutrow;

         mastercutrow = GCGextendedmasterconsGetRow(mastersepacut);
         assert(mastercutrow != NULL);
         SCIP_CALL( SCIPaddVarToRow(scip, mastercutrow, mastervar, coef) );
      }

      SCIPfreeBufferArrayNull(scip, &pricingvars);
   }

   return SCIP_OKAY;
}


/** remove the separator mastercuts generated and applied at this node from activecuts */
static
SCIP_RETCODE removeStoredCutsFromActiveCuts(
   GCG*                 gcg,            /**< GCG data structure */
   SCIP_CONSDATA*       consdata,       /**< data of current constraint */
   SCIP_CONSHDLRDATA*   conshdlrdata    /**< masterbranch constraint handler data */
   )
{
   assert(gcg != NULL);
   assert(consdata != NULL);

   SCIP_CALL( GCGsepacutShrinkActiveCuts(gcg, consdata->firstnewcut, conshdlrdata->eventhdlr) );


   return SCIP_OKAY;
}

/** add the separator mastercuts generated and applied at this node to active cuts */
static
SCIP_RETCODE addStoredCutsToActiveCuts(
   GCG*               gcg,                /**< GCG data structure */
   SCIP_CONSDATA*     consdata,           /**< consdata of current constraint */
   SCIP_CONSHDLRDATA* conshdlrdata        /**< constraint handler data */
   )
{
   int nactivecuts;
   int j;

   assert(gcg != NULL);
   assert(consdata != NULL);
   assert(consdata->addedcutsinit);
   assert(conshdlrdata != NULL);

#ifndef MASTERSEP_DEBUG
   SCIPinfoMessage(scip, NULL, "add stored cuts: node %lli of type: %i\n", SCIPnodeGetNumber(consdata->node), SCIPnodeGetType(consdata->node));
#endif

   nactivecuts = GCGsepacutGetNActiveCuts(gcg, conshdlrdata->eventhdlr);
   assert(consdata->firstnewcut == nactivecuts);

   /* store the current number of activecuts */
   consdata->firstnewcut = nactivecuts;

   /* if this node did not store any cuts, we do nothing */
   if( !consdata->nodestoredcuts )
      return SCIP_OKAY;

   assert(SCIPnodeGetType(consdata->node) == SCIP_NODETYPE_FORK
          || SCIPnodeGetType(consdata->node) == SCIP_NODETYPE_PSEUDOFORK
          || SCIPnodeGetType(consdata->node) == SCIP_NODETYPE_SUBROOT);

   for( j = 0; j < consdata->naddedcuts; j++ )
   {
      /* we update the cut to include all the master variables which were generated while the cut was inactive */
      SCIP_CALL( addMissedVariables(gcg, consdata->addedcuts[j]) );
      /* we add this cut to activecuts */
      SCIP_CALL( GCGsepacutAddCutToActiveCuts(gcg, consdata->addedcuts[j], conshdlrdata->eventhdlr) );

#ifndef MASTERSEP_DEBUG
      SCIP_ROW* mastercutrow = NULL;

      SCIP_CALL( GCGextendedmasterconsGetRow(consdata->addedcuts[j], &mastercutrow) );
      assert(mastercutrow != NULL);
      SCIPinfoMessage(GCGgetMasterprob(gcg), NULL, "add to active cuts: update row %s from node %lli and add it to active cuts\n",
                            SCIProwGetName(mastercutrow), SCIPnodeGetNumber(consdata->node));
#endif

   }

   return SCIP_OKAY;
}

/** stores the separator mastercuts generated and applied at this node to the branchdata */
static
SCIP_RETCODE initializeAddedCuts(
   GCG*                 gcg,             /**< GCG data structure */
   SCIP_CONSDATA*       consdata,        /**< constraint data of current constraint */
   SCIP_CONSHDLRDATA*   conshdlrdata     /**< constraint handler data */
   )
{
   GCG_EXTENDEDMASTERCONSDATA** activecuts;
   int nactivecuts;
   int j;

   assert(gcg != NULL);
   assert(consdata != NULL);
   assert(!consdata->addedcutsinit);


#ifndef MASTERSEP_DEBUG
   SCIPinfoMessage(scip, NULL, "init added cuts: remove and free new inactive rows from node %lli of type %i\n",
                      SCIPnodeGetNumber(consdata->node), SCIPnodeGetType(consdata->node));
#endif
   /* cleanup (remove AND free) the rows generated at node which have already been removed from the LP */
   SCIP_CALL( GCGsepacutRemoveNewInactiveRows(gcg, consdata->firstnewcut, conshdlrdata->eventhdlr) );

   /* the only type of nodes which can store rows */
   if( !(SCIPnodeGetType(consdata->node) == SCIP_NODETYPE_FORK
         || SCIPnodeGetType(consdata->node) == SCIP_NODETYPE_PSEUDOFORK
         || SCIPnodeGetType(consdata->node) == SCIP_NODETYPE_SUBROOT) )
   {
      consdata->nodestoredcuts = FALSE;
      consdata->addedcutsinit = TRUE;
      return SCIP_OKAY;
   }

   nactivecuts = GCGsepacutGetNActiveCuts(gcg, conshdlrdata->eventhdlr);
   activecuts = GCGsepacutGetActiveCuts(gcg, conshdlrdata->eventhdlr);
   consdata->naddedcuts = 0;
   assert(consdata->firstnewcut <= nactivecuts);

   /* no cuts were applied at this node --> nothing to store */
   if( consdata->firstnewcut == nactivecuts )
   {
#ifndef MASTERSEP_DEBUG
      SCIPinfoMessage(scip, NULL, "init consdata: no cuts to store\n");
#endif
      consdata->nodestoredcuts = FALSE;
      consdata->addedcutsinit = TRUE;
      return SCIP_OKAY;
   }
   consdata->naddedcuts = nactivecuts - consdata->firstnewcut;
#ifndef MASTERSEP_DEBUG
   SCIPinfoMessage(scip, NULL, "init consdata: store %i cuts for (pseudo)fork/subroot\n",
                         consdata->naddedcuts);
#endif

   /* copy, store and capture the cuts */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(GCGgetMasterprob(gcg), &(consdata->addedcuts),
                                            &(activecuts[consdata->firstnewcut]), consdata->naddedcuts) );
   for( j = 0; j < consdata->naddedcuts; j++ )
   {
      SCIP_CALL( GCGcaptureMasterSepaCut(consdata->addedcuts[j]) );
   }
   consdata->nodestoredcuts = TRUE;
   consdata->addedcutsinit = TRUE;

   return SCIP_OKAY;
}

/** initialize the consdata data structure */
static
SCIP_RETCODE initializeConsdata(
   GCG*                  gcg,                /**< GCG data structure*/
   SCIP_CONS*            cons                /**< constraint for which the consdata is created */
   )
{
#ifdef SCIP_DEBUG
   SCIP_CONS* origcons_parent;
   SCIP_CONS* parent_origcons;
#endif

   SCIP* masterprob;
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* origcons;

   SCIP_DOMCHG* domchg;
   SCIP_BOUNDCHG* boundchg;
   SCIP_VAR* boundchgvar;

   int i;

   assert(gcg != NULL);

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));
   assert(cons != NULL);

   /* get constraint handler */
   conshdlr = SCIPfindConshdlr(masterprob, CONSHDLR_NAME);
   assert(conshdlr != NULL);
   assert(SCIPconsGetHdlr(cons) == conshdlr);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get corresponding origbranch constraint in the original problem */
   origcons = GCGconsOrigbranchGetActiveCons(gcg);
   assert(origcons != NULL);

   if( consdata->origcons == NULL ) /*rootnode?*/
   {
      SCIPdebugMessage("set root origcons\n");
      consdata->origcons = origcons;
      GCGconsOrigbranchSetMastercons(origcons, cons);
   }
   else if( consdata->origcons != origcons )
   {
      // todo: Check this case.
      SCIPdebugMessage("B&B trees could be out of sync\n");
   }

   /* @fixme: Why should anything else happen? */
   if( GCGconsOrigbranchGetNChildconss(origcons) == 0 )
   {
      consdata->nchildconss = 0;
      consdata->childconss = NULL;
   }

   /*GCGconsOrigbranchSetMastercons(origcons, cons);*/


#ifdef SCIP_DEBUG
   if( consdata->parentcons != NULL )
      parent_origcons = SCIPconsGetData(consdata->parentcons)->origcons;
   else
      parent_origcons = NULL;

   if( consdata->origcons != NULL )
      origcons_parent = GCGconsOrigbranchGetParentcons(consdata->origcons);
   else
      origcons_parent = NULL;

   SCIPdebugMessage("cons: %s, origcons: %s, parent: %s => %s\n", SCIPconsGetName(cons), consdata->origcons == NULL? "NULL" : SCIPconsGetName( consdata->origcons ),
      parent_origcons == NULL? "NULL" : SCIPconsGetName(parent_origcons), origcons_parent == NULL? "NULL" : SCIPconsGetName(origcons_parent) );
#endif

   assert(SCIPgetCurrentNode(masterprob) == consdata->node || consdata->node == SCIPgetRootNode(masterprob));
/*    assert((SCIPgetNNodesLeft(masterprob)+SCIPgetNNodes(masterprob) == 1) == (consdata->node == SCIPgetRootNode(masterprob))); */
   assert(SCIPnodeGetDepth(GCGconsOrigbranchGetNode(consdata->origcons)) == SCIPnodeGetDepth(consdata->node));
   assert(consdata->parentcons != NULL || SCIPnodeGetDepth(consdata->node) == 0);

   assert(consdata->parentcons == NULL ||
      SCIPconsGetData(consdata->parentcons)->origcons == GCGconsOrigbranchGetParentcons(consdata->origcons));

   consdata->maxlocalbndchgstreated = SCIPcalcMemGrowSize(masterprob, conshdlrdata->nstack+1);
   SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &consdata->nlocalbndchgstreated, (size_t)consdata->maxlocalbndchgstreated) );

   /* get all bound changes at the corresponding node in the original problem */

   domchg = SCIPnodeGetDomchg(GCGconsOrigbranchGetNode(origcons));
   consdata->nlocalbndchgs = SCIPdomchgGetNBoundchgs(domchg);
   consdata->nlocalbndchgstreated[conshdlrdata->nstack] = consdata->nlocalbndchgs;

   if( consdata->nlocalbndchgs > 0 )
   {
      consdata->maxlocalbndchgs = SCIPcalcMemGrowSize(masterprob, consdata->nlocalbndchgs);
      SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &consdata->localbndvars, consdata->maxlocalbndchgs) );
      SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &consdata->localbndtypes, consdata->maxlocalbndchgs) );
      SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &consdata->localnewbnds, consdata->maxlocalbndchgs) );
      SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &consdata->localoldbnds, consdata->maxlocalbndchgs) );
   }
   else
   {
      consdata->maxlocalbndchgs = 0;
   }

   consdata->nbranchingchgs = 0;

   for( i = 0; i < consdata->nlocalbndchgs; ++i )
   {
      boundchg = SCIPdomchgGetBoundchg(domchg, i);

      consdata->localbndvars[i] = SCIPboundchgGetVar(boundchg);
      consdata->localnewbnds[i] = SCIPboundchgGetNewbound(boundchg);
      consdata->localbndtypes[i] = SCIPboundchgGetBoundtype(boundchg);

      if( SCIPboundchgGetBoundchgtype(boundchg) == SCIP_BOUNDCHGTYPE_BRANCHING )
      {
         consdata->nbranchingchgs++;
         assert(consdata->nbranchingchgs == i+1);
      }
   }

   consdata->needprop = TRUE;

   assert((consdata->parentcons == NULL) == (conshdlrdata->nstack == 0));
   if( consdata->parentcons != NULL )
   {
      SCIP_CONSDATA* parentdata = SCIPconsGetData(consdata->parentcons);

      assert(consdata->parentcons == conshdlrdata->stack[conshdlrdata->nstack-1]);
      assert(SCIPconsGetData(conshdlrdata->stack[0])->parentcons == NULL);

      /* check whether bound changes were added in nodes on the path
       * to the current node after activation of the parent node
       */
      for( i = 1; i < conshdlrdata->nstack; ++i )
      {
         int ndomboundchgs;
         SCIP_CONSDATA* stackconsdata = SCIPconsGetData(conshdlrdata->stack[i]);
         domchg = SCIPnodeGetDomchg(GCGconsOrigbranchGetNode(stackconsdata->origcons));
         ndomboundchgs = SCIPdomchgGetNBoundchgs(domchg);

         assert(ndomboundchgs >= parentdata->nlocalbndchgstreated[i]);

         if( ndomboundchgs != parentdata->nlocalbndchgstreated[i] )
         {
            int diff;
            int j;

            diff = ndomboundchgs - parentdata->nlocalbndchgstreated[i];

            if( consdata->maxlocalbndchgs < consdata->nlocalbndchgs + diff )
            {
               int newmaxsize = SCIPcalcMemGrowSize(masterprob, consdata->nlocalbndchgs + diff);
               SCIP_CALL(SCIPreallocBlockMemoryArray(masterprob, &consdata->localbndvars, consdata->maxlocalbndchgs,
                                                     (size_t) newmaxsize));
               SCIP_CALL(SCIPreallocBlockMemoryArray(masterprob, &consdata->localbndtypes, consdata->maxlocalbndchgs,
                                                     (size_t) newmaxsize));
               SCIP_CALL(SCIPreallocBlockMemoryArray(masterprob, &consdata->localnewbnds, consdata->maxlocalbndchgs,
                                                     (size_t) newmaxsize));
               SCIP_CALL(SCIPreallocBlockMemoryArray(masterprob, &consdata->localoldbnds, consdata->maxlocalbndchgs,
                                                     (size_t) newmaxsize));
               consdata->maxlocalbndchgs = newmaxsize;
            }

            /* add bound changes to the boundchanges array */
            for( j = 0; j < ndomboundchgs; ++j )
            {
               int bndchgindex;
               SCIP_BOUNDTYPE boundchgtype;
               SCIP_Real boundchgnewbound;

               boundchg = SCIPdomchgGetBoundchg(domchg, j);
               boundchgvar = SCIPboundchgGetVar(boundchg);
               boundchgtype = SCIPboundchgGetBoundtype(boundchg);
               boundchgnewbound = SCIPboundchgGetNewbound(boundchg);

               if( j < stackconsdata->nlocalbndchgstreated[i] )
               {
                  assert(stackconsdata->localbndvars[j] == boundchgvar
                     && SCIPisEQ(masterprob, stackconsdata->localnewbnds[j], boundchgnewbound)
                     && stackconsdata->localbndtypes[j] == boundchgtype);
                  continue;
               }
               if( j < parentdata->nlocalbndchgstreated[i] )
                  continue;

               bndchgindex = consdata->nlocalbndchgs + j - parentdata->nlocalbndchgstreated[i];

               consdata->localbndvars[bndchgindex] = boundchgvar;
               consdata->localnewbnds[bndchgindex] = boundchgnewbound;
               consdata->localbndtypes[bndchgindex] = boundchgtype;
            }

            consdata->nlocalbndchgs += diff;
         }

         consdata->nlocalbndchgstreated[i] = ndomboundchgs;
      }
       /* store the master cuts (generated and) applied at parent node in parent constraint data, when
        *   - parent master cuts have not been stored yet, because the corresponding constraint was not deactivated */
       if( !parentdata->addedcutsinit && SCIPnodeGetType(consdata->node) == SCIP_NODETYPE_FOCUSNODE )
          SCIP_CALL( initializeAddedCuts(gcg, parentdata, conshdlrdata) );

   }

   /* store the number of activecuts at activation */
   consdata->firstnewcut = GCGsepacutGetNActiveCuts(gcg, conshdlrdata->eventhdlr);
#ifndef MASTERSEP_DEBUG
   SCIPinfoMessage(scip, NULL, "init cons data: firstnewcut %i in node %lli\n", consdata->firstnewcut, SCIPnodeGetNumber(consdata->node));
#endif

   return SCIP_OKAY;
}

/** add a global bound change on the original problem to the pending bound changes array */
static
SCIP_RETCODE addPendingBndChg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable on which the bound change is applied (corresponding master variable copy or pricing variable) */
   SCIP_BOUNDTYPE        boundtype,          /**< type of the bound (lower or upper) */
   SCIP_Real             oldbound,           /**< previous bound value */
   SCIP_Real             newbound            /**< new bound value */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_HASHMAP* pendingvarmap;
   int idx;

   assert(scip != NULL);

   /* get constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert((conshdlrdata->npendingbnds > 0) || conshdlrdata->pendingbndsactivated);

   pendingvarmap = (boundtype == SCIP_BOUNDTYPE_LOWER) ? conshdlrdata->pendingvarmaplb : conshdlrdata->pendingvarmapub;
   idx = SCIPhashmapGetImageInt(pendingvarmap, var);
   if( idx != INT_MAX )
   {
      assert(conshdlrdata->pendingvars[idx] == var);
      assert(conshdlrdata->pendingbndtypes[idx] == boundtype);
      conshdlrdata->pendingnewbnds[idx] = newbound;
   }
   else
   {
      /* reallocate memory if needed */
      if( conshdlrdata->npendingbnds >= conshdlrdata->maxpendingbnds )
      {
         int newsize = SCIPcalcMemGrowSize(scip, conshdlrdata->npendingbnds+5);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(conshdlrdata->pendingvars), conshdlrdata->maxpendingbnds, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(conshdlrdata->pendingbndtypes), conshdlrdata->maxpendingbnds, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(conshdlrdata->pendingnewbnds), conshdlrdata->maxpendingbnds, newsize) );
         conshdlrdata->maxpendingbnds = newsize;
      }

      /* store pending bound change */
      conshdlrdata->pendingvars[conshdlrdata->npendingbnds] = var;
      conshdlrdata->pendingbndtypes[conshdlrdata->npendingbnds] = boundtype;
      conshdlrdata->pendingnewbnds[conshdlrdata->npendingbnds] = newbound;
      SCIP_CALL( SCIPhashmapInsertInt(pendingvarmap, var, conshdlrdata->npendingbnds) );
      conshdlrdata->npendingbnds++;
      conshdlrdata->pendingbndsactivated = FALSE;
   }

   return SCIP_OKAY;
}

/** For a given global bound change on a pricing variable, check if the global bounds on all corresponding original variables are still the same
 *
 *  @return TRUE if the variable is in a relevant block AND all variables identical to it have the same bounds
 */
static
SCIP_Bool checkAggregatedGlobalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             bndvar,             /**< pricing variable whose new global bound is to be checked */
   SCIP_BOUNDTYPE        bndtype,            /**< type of the new global bound */
   SCIP_Real             bound               /**< new global bound */
   )
{
   SCIP_VAR** identvars;
   int nidentvars;
   SCIP_Bool identical;

   assert(GCGvarIsPricing(bndvar));

   /* get all identical variables */
   identvars = GCGpricingVarGetOrigvars(bndvar);
   nidentvars = GCGpricingVarGetNOrigvars(bndvar);

   identical = TRUE;

   /* If the variable was not aggregated, there is nothing to check */
   if( nidentvars > 1 )
   {
      int i;

      /* Check if the bounds of all identical variables are equal to the one of the representative */
      for( i = 0; i < nidentvars; ++i )
      {
         SCIP_Real identbound = bndtype == SCIP_BOUNDTYPE_UPPER ? SCIPvarGetUbGlobal(identvars[i]) : SCIPvarGetLbGlobal(identvars[i]);
         if( !SCIPisEQ(scip, identbound, bound) )
         {
            SCIPwarningMessage(scip, "Var <%s> has new global %s bound %g, but identical var <%s> has %g -- don't know how to handle!\n",
               SCIPvarGetName(bndvar), bndtype == SCIP_BOUNDTYPE_UPPER ? "upper" : "lower",
                  bound, SCIPvarGetName(identvars[i]), identbound);
            identical = FALSE;
         }
      }
   }

   return identical;
}

/** apply global bound changes on original problem variables either
 *  to their copies in the master problem and/or to the corresponding pricing problem variables
 */
static
SCIP_RETCODE applyGlobalBndchgsToPricingprobs(
   GCG*                  gcg,                /* GCG data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /* constraint handler data */
   )
{
   SCIP* masterprob;
   SCIP* origscip;
   int i;

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));

   assert(conshdlrdata != NULL);

   /* get original problem */
   origscip = GCGgetOrigprob(gcg);
   assert(origscip != NULL);

   if( !conshdlrdata->pendingbndsactivated )
   {
      assert(conshdlrdata->npendingbnds > 0);
      for( i = 0; i < conshdlrdata->npendingbnds; i++ )
      {
         /* This should not have an effect on linking variables */
         assert(GCGvarIsMaster(conshdlrdata->pendingvars[i]) || GCGvarIsPricing(conshdlrdata->pendingvars[i]));

         if( GCGvarIsMaster(conshdlrdata->pendingvars[i]) )
         {
            if( conshdlrdata->pendingbndtypes[i] == SCIP_BOUNDTYPE_LOWER )
            {
               if( SCIPisLT(masterprob, SCIPvarGetLbGlobal(conshdlrdata->pendingvars[i]), conshdlrdata->pendingnewbnds[i]) )
               {
                  SCIPdebugMessage("Global lower bound of mastervar <%s> set to %g\n", SCIPvarGetName(conshdlrdata->pendingvars[i]),
                     conshdlrdata->pendingnewbnds[i]);
                  SCIP_CALL( SCIPchgVarLbGlobal(masterprob, conshdlrdata->pendingvars[i], conshdlrdata->pendingnewbnds[i]) );
               }
            }
            else
            {
               if( SCIPisGT(masterprob, SCIPvarGetUbGlobal(conshdlrdata->pendingvars[i]), conshdlrdata->pendingnewbnds[i]) )
               {
                  SCIPdebugMessage("Global upper bound of mastervar <%s> set to %g\n", SCIPvarGetName(conshdlrdata->pendingvars[i]),
                     conshdlrdata->pendingnewbnds[i]);
                  SCIP_CALL( SCIPchgVarUbGlobal(masterprob, conshdlrdata->pendingvars[i], conshdlrdata->pendingnewbnds[i]) );
               }
            }
         }
         else
         {
            /* this is a global boundchange on a variable that belongs to a block,
             * we have to adjust the bound of the corresponding variable in the pricing problem
             */

            /* check if all identical variables have the same global bound */
            if( !checkAggregatedGlobalBounds(origscip, conshdlrdata->pendingvars[i], conshdlrdata->pendingbndtypes[i], conshdlrdata->pendingnewbnds[i]) )
               continue;

            if( conshdlrdata->pendingbndtypes[i] == SCIP_BOUNDTYPE_LOWER )
            {
               SCIP_CALL( SCIPchgVarLbGlobal(GCGgetPricingprob(gcg, GCGvarGetBlock(conshdlrdata->pendingvars[i]) ),
                     conshdlrdata->pendingvars[i], conshdlrdata->pendingnewbnds[i]) );
            }
            else
            {
               SCIP_CALL( SCIPchgVarUbGlobal(GCGgetPricingprob(gcg, GCGvarGetBlock(conshdlrdata->pendingvars[i]) ),
                     conshdlrdata->pendingvars[i], conshdlrdata->pendingnewbnds[i]) );
            }
         }
      }
      conshdlrdata->pendingbndsactivated = TRUE;
   }

   return SCIP_OKAY;
}

/** apply global bound changes on original problem variables to the master problem */
static
SCIP_RETCODE applyGlobalBndchgsToPricedMastervars(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   int*                  propcount           /**< number of applied bound changes */
   )
{
   SCIP* masterprob;
   SCIP_VAR** vars;
   int nvars;
   int i;
   int k;

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));

   /* get master problem variables generated during pricing */
   vars = GCGmasterGetPricedvars(gcg);
   nvars = GCGmasterGetNPricedvars(gcg);

   assert((conshdlrdata->npendingbnds > 0) || conshdlrdata->pendingbndsactivated);

   /* iterate over all master variables and apply global bound changes */
   if( conshdlrdata->npendingbnds > 0 && conshdlrdata->pendingbndsactivated )
   {
      for( i = 0; i < nvars; i++ )
      {
         SCIP_Bool ismastervarrelevant;
         SCIP_VAR** origvars;
         SCIP_HASHMAP* origvals;
         int norigvars;
         int blocknr;
         blocknr = GCGvarGetBlock(vars[i]);

         /* get the original variables that are contained in the master variable */
         assert(GCGvarIsMaster(vars[i]));
         norigvars = GCGmasterVarGetNOrigvars(vars[i]);
         origvars = GCGmasterVarGetOrigvars(vars[i]);
         origvals = GCGmasterVarGetOrigvalmap(vars[i]);

         assert(blocknr < GCGgetNPricingprobs(conshdlrdata->gcg));
         assert(norigvars >= 0);
         assert(origvars != NULL || norigvars == 0);

         /* only look at master variables not globally fixed to zero that belong to a block */
         ismastervarrelevant = !SCIPisFeasZero(masterprob, SCIPvarGetUbGlobal(vars[i]));
         ismastervarrelevant = ismastervarrelevant && (norigvars > 0);
         ismastervarrelevant = ismastervarrelevant && (blocknr >= 0 || GCGmasterVarIsLinking(vars[i])); /*lint !e613*/
         if( !ismastervarrelevant )
            continue;

         /* iterate over global bound changes on the original variables
          * that have not yet been checked for the master variables
          */
         for( k = 0; k < conshdlrdata->npendingbnds; k++ )
         {
            int bndchgblocknr;
            SCIP_VAR** bndchgorigvars;
            SCIP_Real val;

            assert(!GCGvarIsOriginal(conshdlrdata->pendingvars[k]));

            bndchgblocknr = GCGvarGetBlock(conshdlrdata->pendingvars[k]);
            if( GCGvarIsMaster(conshdlrdata->pendingvars[k]) )
            {
               assert(bndchgblocknr == -1);
               bndchgorigvars = GCGmasterVarGetOrigvars(conshdlrdata->pendingvars[k]);
            }
            else if( GCGvarIsPricing(conshdlrdata->pendingvars[k]) )
               bndchgorigvars = GCGpricingVarGetOrigvars(conshdlrdata->pendingvars[k]);
            else
            {
               SCIPerrorMessage("Variable %s is not pricing nor master.\n", SCIPvarGetName(conshdlrdata->pendingvars[k]));
               assert(GCGvarIsMaster(conshdlrdata->pendingvars[k]) || GCGvarIsPricing(conshdlrdata->pendingvars[k]));
               bndchgorigvars = NULL;
            }
            assert(bndchgblocknr < GCGgetNPricingprobs(conshdlrdata->gcg));
            assert(bndchgorigvars != NULL);
            assert(origvars != NULL);

            /* The bound change is only relevant for the master variable if either
             *  - the bound change was performed in the same block as the master variable, or
             *  - the master variable is a copied linking variable and the bound change was performed
             *    in one of the blocks that the variable is linking
             */
            if( (bndchgblocknr != blocknr )
               && !(GCGmasterVarIsLinking(vars[i]) && GCGisLinkingVarInBlock(origvars[0], bndchgblocknr)) )
               continue;

            assert(bndchgorigvars[0] != NULL);

            val = SCIPhashmapGetImageReal(origvals, bndchgorigvars[0]);
            /* variables belong to the same block -> set origval to 0.0 if not in map */
            if( val == SCIP_INVALID )
               val = 0.0;

            /* if the variable contains a part of the branching variable that violates the bound,
             * fix the master variable to 0
             */
            /* @todo: This is the wrong way to treat bound changes on original variable copies in the master problem;
             * I think they have already been treated during constraint activation
             */

            assert(GCGvarGetBlock(bndchgorigvars[0]) == blocknr || (GCGisLinkingVarInBlock(bndchgorigvars[0], blocknr))); /*lint !e613*/

            /* new lower bound */
            if( conshdlrdata->pendingbndtypes[k] == SCIP_BOUNDTYPE_LOWER &&
               SCIPisFeasLT(masterprob, val, conshdlrdata->pendingnewbnds[k]) )
            {
               SCIP_CALL( SCIPchgVarUbGlobal(masterprob, vars[i], 0.0) );
               ++(*propcount);
               break;
            }
            /* new upper bound */
            if( conshdlrdata->pendingbndtypes[k] == SCIP_BOUNDTYPE_UPPER &&
               SCIPisFeasGT(masterprob, val, conshdlrdata->pendingnewbnds[k]) )
            {
               SCIP_CALL( SCIPchgVarUbGlobal(masterprob, vars[i], 0.0) );
               ++(*propcount);
               break;
            }
         }
      }
      conshdlrdata->pendingbndsactivated = TRUE;
      conshdlrdata->npendingbnds = 0;
      SCIPhashmapRemoveAll(conshdlrdata->pendingvarmaplb);
      SCIPhashmapRemoveAll(conshdlrdata->pendingvarmapub);

      SCIPdebugMessage("Finished handling of pending global bound changes: %d changed bounds\n", *propcount);
   }

   return SCIP_OKAY;
}

/** reset bound changes on pricing variables (called when a node is deactivated) */
static
SCIP_RETCODE resetPricingVarBound(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_VAR*             pricingvar,         /**< variable on which the bound change should be performed */
   SCIP_CONSDATA*        consdata,           /**< constraint data structure */
   int                   i,                  /**< index of the information in the constraint data structure */
   int                   blocknr             /**< number of the pricing problem */
   )
{
   SCIP* masterprob = GCGgetMasterprob(gcg);

   assert(masterprob != NULL);
   assert(pricingvar != NULL);
   assert(consdata != NULL);
   assert(consdata->nactivated >= 1);
   assert(consdata->localbndvars != NULL);
   assert(consdata->localbndtypes != NULL);
   assert(consdata->localnewbnds != NULL);
   assert(consdata->localoldbnds != NULL);

   assert(blocknr >= 0 && blocknr < GCGgetNPricingprobs(gcg));

   /* lower bound was changed */
   if( consdata->localbndtypes[i] == SCIP_BOUNDTYPE_LOWER )
   {
      if( GCGgetNIdenticalBlocks(gcg, blocknr) > 1 || GCGgetNIdenticalBlocks(gcg, blocknr) == 0 )
         return SCIP_OKAY;

      assert(SCIPisGE(masterprob, SCIPvarGetLbLocal(pricingvar), consdata->localnewbnds[i])
         || SCIPisLE(masterprob, SCIPvarGetLbLocal(pricingvar), SCIPvarGetLbGlobal(consdata->localbndvars[i])));

      if( SCIPisEQ(masterprob, SCIPvarGetLbGlobal(consdata->localbndvars[i]), consdata->localnewbnds[i]) )
         return SCIP_OKAY;

      if( SCIPisGT(masterprob, consdata->localoldbnds[i], consdata->localnewbnds[i]) )
         return SCIP_OKAY;

      if( SCIPisGT(masterprob, SCIPvarGetLbGlobal(consdata->localbndvars[i]), consdata->localoldbnds[i]) )
      {
         SCIP_CALL( SCIPchgVarLb(GCGgetPricingprob(gcg, blocknr), pricingvar, SCIPvarGetLbGlobal(consdata->localbndvars[i])) );
         SCIPdebugMessage("relaxed lower bound of pricing var %s from %g to global bound %g (%s)\n",
            SCIPvarGetName(pricingvar), consdata->localnewbnds[i], SCIPvarGetLbGlobal(consdata->localbndvars[i]), consdata->name);
      }
      else
      {
         SCIP_CALL( SCIPchgVarLb(GCGgetPricingprob(gcg, blocknr), pricingvar, consdata->localoldbnds[i]) );
         SCIPdebugMessage("relaxed lower bound of pricing var %s from %g to %g (%s)\n",
            SCIPvarGetName(pricingvar), consdata->localnewbnds[i], consdata->localoldbnds[i], consdata->name);
      }
   }
   /* upper bound was changed */
   else
   {
      if( GCGgetNIdenticalBlocks(gcg, blocknr) > 1 || GCGgetNIdenticalBlocks(gcg, blocknr) == 0 )
         return SCIP_OKAY;

      assert(SCIPisLE(masterprob, SCIPvarGetUbLocal(pricingvar), consdata->localnewbnds[i])
         || SCIPisGE(masterprob, SCIPvarGetUbLocal(pricingvar), SCIPvarGetUbGlobal(consdata->localbndvars[i])));

      if( SCIPisEQ(masterprob, SCIPvarGetUbGlobal(consdata->localbndvars[i]), consdata->localnewbnds[i]) )
         return SCIP_OKAY;

      if( SCIPisLT(masterprob, consdata->localoldbnds[i], consdata->localnewbnds[i]) )
         return SCIP_OKAY;

      if( SCIPisLT(masterprob, SCIPvarGetUbGlobal(consdata->localbndvars[i]), consdata->localoldbnds[i]) )
      {
         SCIP_CALL( SCIPchgVarUb(GCGgetPricingprob(gcg, blocknr), pricingvar, SCIPvarGetUbGlobal(consdata->localbndvars[i])) );
         SCIPdebugMessage("relaxed upper bound of pricing var %s from %g to global bound %g (%s)\n",
            SCIPvarGetName(pricingvar), consdata->localnewbnds[i], SCIPvarGetUbGlobal(consdata->localbndvars[i]), consdata->name);
      }
      else
      {
         SCIP_CALL( SCIPchgVarUb(GCGgetPricingprob(gcg, blocknr), pricingvar, consdata->localoldbnds[i]) );
         SCIPdebugMessage("relaxed upper bound of pricing var %s from %g to %g (%s)\n",
            SCIPvarGetName(pricingvar), consdata->localnewbnds[i], consdata->localoldbnds[i], consdata->name);
      }
   }

   return SCIP_OKAY;
}

/** perform bound changes on pricing variables (called when a node is activated) */
static
SCIP_RETCODE tightenPricingVarBound(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_VAR*             pricingvar,         /**< variable on which the bound change should be performed */
   SCIP_CONSDATA*        consdata,           /**< constraint data structure */
   int                   i,                  /**< index of the information in the constraint data structure */
   int                   blocknr             /**< number of the pricing problem */
   )
{
   SCIP* masterprob;

   assert(gcg != NULL);
   assert(pricingvar != NULL);
   assert(consdata != NULL);
   assert(consdata->nactivated >= 1);
   assert(consdata->localbndvars != NULL);
   assert(consdata->localbndtypes != NULL);
   assert(consdata->localnewbnds != NULL);
   assert(consdata->localoldbnds != NULL);

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);
   assert(blocknr >= 0 && blocknr < GCGgetNPricingprobs(gcg));

   /* lower bound was changed */
   if( consdata->localbndtypes[i] == SCIP_BOUNDTYPE_LOWER )
   {
      consdata->localoldbnds[i] = SCIPvarGetLbLocal(pricingvar);

      if( SCIPisGT(masterprob, consdata->localnewbnds[i], consdata->localoldbnds[i]) )
      {
         SCIP_CALL( SCIPchgVarLb(GCGgetPricingprob(gcg, blocknr), pricingvar, consdata->localnewbnds[i]) );
         SCIPdebugMessage("tightened lower bound of var %s from %g to %g\n",
            SCIPvarGetName(pricingvar), consdata->localoldbnds[i], consdata->localnewbnds[i]);
      }
   }
   /* upper bound was changed */
   else
   {
      assert(consdata->localbndtypes[i] == SCIP_BOUNDTYPE_UPPER);

      consdata->localoldbnds[i] = SCIPvarGetUbLocal(pricingvar);

      if( SCIPisLT(masterprob, consdata->localnewbnds[i], consdata->localoldbnds[i]) )
      {
         SCIP_CALL( SCIPchgVarUb(GCGgetPricingprob(gcg, blocknr), pricingvar, consdata->localnewbnds[i]) );
         SCIPdebugMessage("tightened upper bound of var %s from %g to %g\n",
            SCIPvarGetName(pricingvar), consdata->localoldbnds[i], consdata->localnewbnds[i]);
      }
   }

   return SCIP_OKAY;
}

/** For a given local bound change on an original variable, check if the bounds on the variables identical to it are the same
 *
 *  @note If the variable is represented by another one, we check only the representative;
 *        otherwise, we check all variables identical to it
 *
 *  @return TRUE if the variable is in a relevant block AND all variables identical to it have the same bounds
 */
static
SCIP_Bool checkAggregatedLocalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            bndvars,            /**< all variables whose local bounds were changed */
   SCIP_Real*            bounds,             /**< corresponding new local bounds */
   int                   nbndvars,           /**< number of all variables whose local bounds were changed */
   SCIP_VAR*             bndvar,             /**< original variable whose local bound was changed and which is to be checked */
   SCIP_BOUNDTYPE        bndtype,            /**< type of the new local bound */
   SCIP_VAR*             pricingvar          /**< pricing variable corresponding to the original variable */
   )
{
   SCIP_VAR** identvars;
   int nidentvars;

   assert(GCGvarIsOriginal(bndvar));
   assert(GCGvarIsPricing(pricingvar));

   /* get variables with which the original variable was aggregated */
   identvars = GCGpricingVarGetOrigvars(pricingvar);
   nidentvars = GCGpricingVarGetNOrigvars(pricingvar);

   /* First case: The variable is not represented by another one - check the bounds of all variables it represents */
   if( identvars[0] == bndvar )
   {
      SCIP_Bool identical = TRUE;

      /* If the variable was not aggregated, there is nothing to check */
      if( nidentvars > 1 )
      {
         int i;
         int j;
         SCIP_Real* identbounds;  /* most recent bounds of all identical variables */

         SCIP_CALL( SCIPallocBufferArray(scip, &identbounds, nidentvars) );
         for( j = 0; j < nidentvars; ++j )
            identbounds[j] = SCIP_INVALID;

         /* For all variables on which a bound was changed *and* which are identical to the current variable,
          * get the most recent bound
          */
         for( i = 0; i < nbndvars; ++i )
         {
            assert(GCGvarIsOriginal(bndvars[i]));

            if( GCGvarGetBlock(bndvars[i]) < 0 )
               continue;

            if( GCGpricingVarGetOrigvars(GCGoriginalVarGetPricingVar(bndvars[i]))[0] == identvars[0] )
               for( j = 0; j < nidentvars; ++j )
                  if( identvars[j] == bndvars[i] )
                     identbounds[j] = bounds[i];
         }

         /* Check if the bounds of all identical variables are equal to the one of the representative */
         for( j = 1; j < nidentvars; ++j )
         {
            if( !SCIPisEQ(scip, identbounds[j], identbounds[0]) )
            {
               SCIPwarningMessage(scip, "Var <%s> has new local %s bound %g, but identical var <%s> has %g -- don't know how to handle!\n",
                  SCIPvarGetName(bndvar), bndtype == SCIP_BOUNDTYPE_UPPER ? "upper" : "lower",
                     identbounds[0], SCIPvarGetName(identvars[j]), identbounds[j]);
               identical = FALSE;
            }
         }

         SCIPfreeBufferArray(scip, &identbounds);
      }

      return identical;
   }

   /* Second case: The variable is represented by another one due to aggregation; check if its representative has the same bound */
   else
   {
      int i;
      SCIP_Real reprbound;
      SCIP_Real bound;

      /* Get the most recent bound for the bound change variable as well as for its representative */
      reprbound = SCIP_INVALID;
      bound = SCIP_INVALID;
      for( i = 0; i < nbndvars; ++i )
      {
         assert(GCGvarIsOriginal(bndvars[i]));

         if( bndvars[i] == identvars[0] )
            reprbound = bounds[i];
         else if( bndvars[i] == bndvar )
            bound = bounds[i];
      }

      /* Check if the bounds are equal */
      if( !SCIPisEQ(scip, bound, reprbound) )
      {
         SCIPwarningMessage(scip, "Var <%s> has new local %s bound %g, but representative <%s> has %g -- don't know how to handle!\n",
            SCIPvarGetName(bndvar), bndtype == SCIP_BOUNDTYPE_UPPER ? "upper" : "lower",
            bound, SCIPvarGetName(identvars[0]), reprbound);
      }

      /* Since the block is not relevant, there is no corresponding pricing variable */
      return FALSE;
   }
}

/** apply local bound changes in the original problem to the pricing problems */
static
SCIP_RETCODE applyLocalBndchgsToPricingprobs(
   GCG*                  gcg,                /* GCG data structure */
   SCIP_CONS*            cons                /* current masterbranch constraint */
   )
{
   SCIP* masterprob;
   SCIP_CONSDATA* consdata;
   int i;

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* iterate over all local bound changes in the original problem */
   for( i = 0; i < consdata->nlocalbndchgs; i++ )
   {
      int blocknr;
      assert(GCGvarIsOriginal(consdata->localbndvars[i]));
      blocknr = GCGvarGetBlock(consdata->localbndvars[i]);
      assert(blocknr < GCGgetNPricingprobs(gcg));

      /* if variable belongs to no block, skip it here because the bound changes are treated in the propagation */
      if( blocknr == -1 )
         continue;

      else if( blocknr >= 0 )
      {
         if( checkAggregatedLocalBounds(masterprob, consdata->localbndvars, consdata->localnewbnds, consdata->nlocalbndchgs, consdata->localbndvars[i],
         consdata->localbndtypes[i], GCGoriginalVarGetPricingVar(consdata->localbndvars[i])) )
         {
            SCIPdebugMessage("adjusting bound of pricing var <%s>\n", SCIPvarGetName(consdata->localbndvars[i]));
            /* set corresponding bound in the pricing problem */
            SCIP_CALL( tightenPricingVarBound(gcg, GCGoriginalVarGetPricingVar(consdata->localbndvars[i]), consdata, i, blocknr) );
         }
      }

      else if( blocknr == -2 )
      {
         int j;
         int npricingprobs;
         SCIP_VAR** pricingvars;
         SCIP_Bool aggregated;

         npricingprobs = GCGgetNPricingprobs(gcg);
         pricingvars = GCGlinkingVarGetPricingVars(consdata->localbndvars[i]);
         aggregated = FALSE;
         /* check the blocks in which the linking variable appears */
         for( j = 0; j < npricingprobs; ++j )
         {
            if( pricingvars[j] == NULL )
               continue;

            if( !checkAggregatedLocalBounds(masterprob, consdata->localbndvars, consdata->localnewbnds, consdata->nlocalbndchgs,
               consdata->localbndvars[i], consdata->localbndtypes[i], pricingvars[j]) )
               aggregated = TRUE;
         }
         if( aggregated )
            continue;

         SCIPdebugMessage("adjusting bound of linking pricing var <%s>\n", SCIPvarGetName(consdata->localbndvars[i]));

         /* set corresponding bound in the pricing problem */
         for( j = 0; j < npricingprobs; ++j )
         {
            if( pricingvars[j] == NULL )
               continue;

            SCIP_CALL( tightenPricingVarBound(gcg, pricingvars[j], consdata, i, j) );
         }
      }

      else
      {
         SCIPerrorMessage("blocknr = %d is not valid! This is a serious error!", GCGvarGetBlock(consdata->localbndvars[i]));
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** undo local bound changes in the original problem to the pricing problems */
static
SCIP_RETCODE undoLocalBndchgsToPricingprobs(
   GCG*                  gcg,                /* GCG data structure */
   SCIP_CONS*            cons                /* current masterbranch constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* iterate over all local bound changes in the original problem */
   for( i = consdata->nlocalbndchgs - 1; i >= 0; i-- )
   {
      int blocknr;
      blocknr = GCGvarGetBlock(consdata->localbndvars[i]);
      assert(GCGvarIsOriginal(consdata->localbndvars[i]));
      assert(blocknr < GCGgetNPricingprobs(gcg));

      /* if variable belongs to no block, local bound in master was set, is reset automatically */
      if( blocknr == -1 )
         continue;

      else if( blocknr >= 0 )
      {
         assert(GCGgetPricingprob(gcg, GCGvarGetBlock(consdata->localbndvars[i])) != NULL);

         /* reset corresponding bound in the pricing problem */

         SCIP_CALL( resetPricingVarBound(gcg,
               GCGoriginalVarGetPricingVar(consdata->localbndvars[i]), consdata, i, blocknr));
      }
      else if( blocknr == -2 )
      {
         int j;
         SCIP_VAR** pricingvars;
         int npricingprobs;

         /* if the variable is linking, we have to perform the same step as above for every existing block*/
         assert(GCGoriginalVarIsLinking(consdata->localbndvars[i]));
         pricingvars = GCGlinkingVarGetPricingVars(consdata->localbndvars[i]);
         npricingprobs = GCGgetNPricingprobs(gcg);

         /* reset corresponding bound in the pricing problem */
         /* lower bound was changed */
         for( j = 0; j < npricingprobs; ++j )
         {
            assert(GCGgetPricingprob(gcg, j) != NULL);
            if( pricingvars[j] == NULL )
               continue;

            assert(GCGgetPricingprob(gcg, j) != NULL);

            /* reset corresponding bound in the pricing problem */
            SCIP_CALL( resetPricingVarBound(gcg, pricingvars[j], consdata, i, j) );
         }
      }
      else
      {
         SCIPerrorMessage("blocknr = %d is not valid! This is a serious error!", blocknr);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE ensureCollectedBndvarsSize(
   SCIP* scip,
   SCIP_CONSHDLRDATA *conshdlrdata,
   int blocknr,
   int minsize
   )
{
   int oldsize = conshdlrdata->maxcollectedbndvars[blocknr];
   if( oldsize < minsize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, minsize);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(conshdlrdata->collectedbndvars[blocknr]), oldsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(conshdlrdata->collectedlbnds[blocknr]), oldsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(conshdlrdata->collectedubnds[blocknr]), oldsize, newsize) );
      conshdlrdata->maxcollectedbndvars[blocknr] = newsize;
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE ensureLinkingvarIndxsSize(
   SCIP* scip,
   SCIP_CONSHDLRDATA *conshdlrdata,
   int minsize
)
{
   int oldsize = conshdlrdata->maxlinkingvaridxs;
   if( oldsize < minsize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, minsize);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(conshdlrdata->linkingvaridxs), oldsize, newsize) );
      BMSclearMemoryArray(&conshdlrdata->linkingvaridxs[oldsize], (newsize - oldsize));
      conshdlrdata->maxlinkingvaridxs = newsize;
   }
   return SCIP_OKAY;
}

/** apply local bound changes on the original variables on newly generated master variables */
static
SCIP_RETCODE applyLocalBndchgsToPricedMastervars(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS*            cons,               /**< current masterbranch constraint */
   int*                  propcount           /**< number of applied bound changes */
   )
{
   SCIP* origprob;
   SCIP* masterprob;
   SCIP_CONSDATA* consdata;
   SCIP_CONSDATA* curconsdata;
   SCIP_VAR** vars;
   int nvars;

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));

   origprob = GCGgetOrigprob(gcg);

   assert(conshdlrdata != NULL);

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get master problem variables generated during pricing */
   vars = GCGmasterGetPricedvars(gcg);
   nvars = GCGmasterGetNPricedvars(gcg);

   if( consdata->npropvars < nvars )
   {
      SCIP_CONS *curcons;
      SCIP_HASHMAP *origvar2idx;
      int nbndvars;
      int nlocalbndchgs;
      int i;
      int j;
      int k;
      SCIP_VAR*** collectedbndvars;
      SCIP_Real** collectedlbnds;
      SCIP_Real** collectedubnds;
      int* ncollectedbndvars;
      int hashmapsize;
      int nblocks;
      int blocknr;
      SCIP_VAR *bndvar;
      int idx;
      SCIP_Real bnd;
      SCIP_Bool islinking;
      int* linkingvaridxs;
      int nlinkingvars;

      nblocks = GCGgetNPricingprobs(gcg);
      nlinkingvars = 0;

      collectedbndvars = conshdlrdata->collectedbndvars;
      collectedlbnds = conshdlrdata->collectedlbnds;
      collectedubnds = conshdlrdata->collectedubnds;
      ncollectedbndvars = conshdlrdata->ncollectedbndvars;

      BMSclearMemoryArray(ncollectedbndvars, conshdlrdata->maxblocknum);

      /* propagate local bound changes on the path from the current node to the root */
      curcons = cons;
      hashmapsize = 0;
      while ( curcons != NULL )
      {
         curconsdata = SCIPconsGetData(curcons);
         assert(curconsdata != NULL);
         //if( curconsdata->npropvars == nvars )
         //   break;
         hashmapsize += (conshdlrdata->enforceproper ? curconsdata->nlocalbndchgs : curconsdata->nbranchingchgs);
         curcons = curconsdata->parentcons;
      }
      if( hashmapsize > SCIPgetNVars(origprob))
         hashmapsize = SCIPgetNVars(origprob);

      SCIP_CALL(SCIPhashmapCreate(&origvar2idx, SCIPblkmem(masterprob), hashmapsize));

      curcons = cons;
      nbndvars = 0;
      while ( curcons != NULL)
      {
         curconsdata = SCIPconsGetData(curcons);
         assert(curconsdata != NULL);
         //if( curconsdata->npropvars == nvars )
         //   break;

         /* propagate all bound changes or only the branching bound changes, depending on the setting for the enforcement of proper variables */
         nlocalbndchgs = (conshdlrdata->enforceproper ? curconsdata->nlocalbndchgs : curconsdata->nbranchingchgs);

         /* iterate over bound changes performed at the current node's equivalent in the original tree */
         for ( k = 0; k < nlocalbndchgs; ++k )
         {
            linkingvaridxs = NULL;
            bndvar = curconsdata->localbndvars[k];
            islinking = GCGoriginalVarIsLinking(bndvar);

            /* get the block the original variable is in */
            blocknr = GCGvarGetBlock(bndvar);
            assert(GCGvarIsOriginal(bndvar));
            assert(blocknr < nblocks);

            if( blocknr >= 0 || islinking)
            {
               idx = SCIPhashmapGetImageInt(origvar2idx, bndvar);
               bnd = curconsdata->localnewbnds[k];

               if( islinking )
               {
                  blocknr = 0;
                  if( idx < INT_MAX )
                  {
                     linkingvaridxs = conshdlrdata->linkingvaridxs[idx];
                  }
                  else
                  {
                     ensureLinkingvarIndxsSize(masterprob, conshdlrdata, nlinkingvars + 1);
                     linkingvaridxs = conshdlrdata->linkingvaridxs[nlinkingvars];
                     if( linkingvaridxs == NULL )
                     {
                        SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &linkingvaridxs, conshdlrdata->maxblocknum) );
                        conshdlrdata->linkingvaridxs[nlinkingvars] = linkingvaridxs;
                     }
                     SCIPhashmapInsertInt(origvar2idx, bndvar, nlinkingvars);
                     ++nlinkingvars;

                     for( i = 0; i < nblocks; ++i )
                        linkingvaridxs[i] = INT_MAX;
                  }
               }

               do
               {
                  assert(blocknr >= 0 && blocknr < GCGgetNPricingprobs(gcg));

                  if( islinking )
                  {
                     if( !GCGisLinkingVarInBlock(bndvar, blocknr) )
                        continue;

                     idx = linkingvaridxs[blocknr];
                  }

                  if( curconsdata->localbndtypes[k] == SCIP_BOUNDTYPE_LOWER )
                  {
                     if( idx < INT_MAX )
                     {
                        if( collectedlbnds[blocknr][idx] == SCIP_INVALID || collectedlbnds[blocknr][idx] < bnd )
                           collectedlbnds[blocknr][idx] = bnd;
                     }
                     else
                     {
                        SCIP_CALL( ensureCollectedBndvarsSize(masterprob, conshdlrdata, blocknr, ncollectedbndvars[blocknr] + 1) );
                        if ( islinking )
                           linkingvaridxs[blocknr] = ncollectedbndvars[blocknr];
                        else
                           SCIPhashmapInsertInt(origvar2idx, bndvar, ncollectedbndvars[blocknr]);
                        collectedbndvars[blocknr][ncollectedbndvars[blocknr]] = bndvar;
                        collectedlbnds[blocknr][ncollectedbndvars[blocknr]] = bnd;
                        collectedubnds[blocknr][ncollectedbndvars[blocknr]] = SCIP_INVALID;
                        ++ncollectedbndvars[blocknr];
                        ++nbndvars;
                     }
                  }
                  else
                  {
                     if( idx < INT_MAX )
                     {
                        if( collectedubnds[blocknr][idx] == SCIP_INVALID || collectedubnds[blocknr][idx] > bnd )
                           collectedubnds[blocknr][idx] = bnd;
                     }
                     else
                     {
                        SCIP_CALL( ensureCollectedBndvarsSize(masterprob, conshdlrdata, blocknr, ncollectedbndvars[blocknr] + 1) );
                        if ( islinking )
                           linkingvaridxs[blocknr] = ncollectedbndvars[blocknr];
                        else
                           SCIPhashmapInsertInt(origvar2idx, bndvar, ncollectedbndvars[blocknr]);
                        collectedbndvars[blocknr][ncollectedbndvars[blocknr]] = bndvar;
                        collectedubnds[blocknr][ncollectedbndvars[blocknr]] = bnd;
                        collectedlbnds[blocknr][ncollectedbndvars[blocknr]] = SCIP_INVALID;
                        ++ncollectedbndvars[blocknr];
                        ++nbndvars;
                     }
                  }
               }
               while( islinking && (++blocknr) < nblocks );
            }
         }
         /* proceed with the parent node */
         curcons = curconsdata->parentcons;
      }

      if( nbndvars > 0 )
      {
         SCIP_Real lbnd;
         SCIP_Real ubnd;
         SCIP_HASHMAP* origvals;
         SCIP_Real origval;

         /* iterate over all master variables created after the current node was left the last time */
         for ( i = consdata->npropvars; i < nvars; i++ )
         {
            assert(GCGvarIsMaster(vars[i]));
            blocknr = GCGvarGetBlock(vars[i]);
            assert(blocknr >= 0 && blocknr < GCGgetNPricingprobs(gcg));

            /** @todo check if this really works with linking variables */

            /* only look at variables not already fixed to 0 or that belong to no block */
            if( SCIPisFeasZero(masterprob, SCIPvarGetUbLocal(vars[i])) )
               continue;

            origvals = GCGmasterVarGetOrigvalmap(vars[i]);

            /* iterate over all original variables whose bound was changed */
            for ( j = 0; j < ncollectedbndvars[blocknr]; j++ )
            {
               bndvar = collectedbndvars[blocknr][j];

               assert(GCGvarGetBlock(bndvar) == blocknr || (GCGoriginalVarIsLinking(bndvar) && GCGisLinkingVarInBlock(bndvar, blocknr)));

               origval = SCIPhashmapGetImageReal(origvals, bndvar);
               /* variables belong to the same block -> set origval to 0.0 if not in map */
               if( origval == SCIP_INVALID )
                  origval = 0.0;
               lbnd = collectedlbnds[blocknr][j];
               ubnd = collectedubnds[blocknr][j];

               /* branching imposes new bound */
               if( (lbnd != SCIP_INVALID && SCIPisFeasLT(masterprob, origval, lbnd))
                   || (ubnd != SCIP_INVALID && SCIPisFeasGT(masterprob, origval, ubnd)) )
               {
                  SCIPdebugMessage("Changing upper bound of var %s\n", SCIPvarGetName(vars[i]));
                  SCIP_CALL(SCIPchgVarUb(masterprob, vars[i], 0.0));
                  ++(*propcount);
                  break;
               }
            }
         }
      }
      SCIPhashmapFree(&origvar2idx);
   }

   SCIPdebugMessage("Finished propagation of newly created variables: %d changed bounds\n", *propcount);

   return SCIP_OKAY;
}

/** apply local bound changes on original variables that have been directly copied to the master problem */
static
SCIP_RETCODE applyLocalBndchgsToCopiedMastervars(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS*            cons,               /**< current masterbranch constraint */
   int*                  propcount           /**< number of applied bound changes */
   )
{
   SCIP* masterprob;
   SCIP_CONSDATA* consdata;
   int i;

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* apply local bound changes */
   for( i = 0; i < consdata->nlocalbndchgs; i++ )
   {
      SCIP_VAR* mastervar;

      assert(GCGvarIsOriginal(consdata->localbndvars[i]));

      /** @todo this might lead to an error with linking variables ? */
      if( GCGvarGetBlock(consdata->localbndvars[i]) >= 0 )
         continue;

      assert(GCGoriginalVarGetNMastervars(consdata->localbndvars[i]) >= 1);

      mastervar = GCGoriginalVarGetMastervars(consdata->localbndvars[i])[0];
      assert(GCGvarGetBlock(mastervar) == -1);

      if( consdata->localbndtypes[i] == SCIP_BOUNDTYPE_LOWER )
      {
         if( SCIPisLT(masterprob, SCIPvarGetLbLocal(mastervar), consdata->localnewbnds[i]) )
         {
            SCIP_CALL( SCIPchgVarLb(masterprob, mastervar, consdata->localnewbnds[i]) );
            ++(*propcount);
            SCIPdebugMessage("changed lb of copied original var %s locally to %g\n", SCIPvarGetName(consdata->localbndvars[i]), consdata->localnewbnds[i]);
         }
      }
      else
      {
         if( SCIPisGT(masterprob, SCIPvarGetUbLocal(mastervar), consdata->localnewbnds[i]) )
         {
            SCIP_CALL( SCIPchgVarUb(masterprob, mastervar, consdata->localnewbnds[i]) );
            ++(*propcount);
            SCIPdebugMessage("changed ub of copied original var %s locally to %g\n", SCIPvarGetName(consdata->localbndvars[i]), consdata->localnewbnds[i]);
         }
      }
   }

   SCIPdebugMessage("Finished propagation of bounds of copied original variables: %d bounds changed.\n", *propcount);

   return SCIP_OKAY;
}

/** forward the seen history */
static
SCIP_RETCODE forwardUpdateSeenHistory(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   SCIP* masterprob;
   assert(gcg != NULL);
   masterprob = GCGgetMasterprob(gcg);
   assert(consdata != NULL);

   if( consdata->branchrule == NULL )
   {
      SCIP_CALL( GCGvarhistoryJumpToLatest(masterprob, &(consdata->knownvarhistory)) );
   }
   else
   {
      SCIP_VAR* var = NULL;
      while( GCGvarhistoryHasNext(consdata->knownvarhistory) )
      {
         SCIP_CALL( GCGvarhistoryNext(masterprob, &consdata->knownvarhistory) );
         SCIP_CALL( GCGvarhistoryGetVar(consdata->knownvarhistory, &var) );
         assert(var != NULL);
         if( SCIPvarIsDeleted(var) )
         {
            SCIPdebugMessage("Skipping deleted Variable <%s>!\n", SCIPvarGetName(var));
            continue;
         }
         SCIP_CALL( GCGrelaxBranchNewCol(gcg, consdata->branchrule, consdata->branchdata, var) );
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods
 */


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free constraint handler storage */
   assert(conshdlrdata->stack == NULL);
   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitMasterbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* prepare stack */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->stack, conshdlrdata->maxstacksize) );
   conshdlrdata->nstack = 0;

   /* prepare pending bound changes */
   conshdlrdata->npendingbnds = 0;
   conshdlrdata->maxpendingbnds = SCIPcalcMemGrowSize(scip, 1);
   conshdlrdata->pendingbndsactivated = TRUE;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->pendingvars), conshdlrdata->maxpendingbnds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->pendingbndtypes), conshdlrdata->maxpendingbnds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->pendingnewbnds), conshdlrdata->maxpendingbnds) );
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->pendingvarmaplb, SCIPblkmem(scip), conshdlrdata->maxpendingbnds) );
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->pendingvarmapub, SCIPblkmem(scip), conshdlrdata->maxpendingbnds) );

   conshdlrdata->maxblocknum = SCIPcalcMemGrowSize(scip, GCGgetNPricingprobs(conshdlrdata->gcg));
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->collectedbndvars),  conshdlrdata->maxblocknum) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->collectedlbnds),  conshdlrdata->maxblocknum) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->collectedubnds),  conshdlrdata->maxblocknum) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->ncollectedbndvars),  conshdlrdata->maxblocknum) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->maxcollectedbndvars),  conshdlrdata->maxblocknum) );
   BMSclearMemoryArray(conshdlrdata->collectedbndvars, conshdlrdata->maxblocknum);
   BMSclearMemoryArray(conshdlrdata->collectedlbnds, conshdlrdata->maxblocknum);
   BMSclearMemoryArray(conshdlrdata->collectedubnds, conshdlrdata->maxblocknum);
   BMSclearMemoryArray(conshdlrdata->ncollectedbndvars, conshdlrdata->maxblocknum);
   BMSclearMemoryArray(conshdlrdata->maxcollectedbndvars, conshdlrdata->maxblocknum);

   conshdlrdata->maxlinkingvaridxs = 0;
   conshdlrdata->linkingvaridxs = NULL;

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolMasterbranch)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create masterbranch constraint for the root node */
   SCIP_CALL( GCGcreateConsMasterbranch(conshdlrdata->gcg, &cons, "root-masterbranch", NULL,  NULL, NULL, NULL, NULL, 0, 0) );
   GCGconsOrigbranchSetMastercons(GCGconsOrigbranchGetActiveCons(conshdlrdata->gcg), cons);

   conshdlrdata->nstack = 1;
   conshdlrdata->stack[0] = cons;
   conshdlrdata->eventhdlr = SCIPfindEventhdlr(scip, "mastersepacut");

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXITSOL(consExitSolMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* we release all the separator mastercuts stored in the data of all the constraints which have not been deleted yet:
    * normally we release the mastercuts via/in consDeleteMasterbranch, but some constraints are only deleted after SCIP_STAGE_EXITSOLVE
    * and at that point, we cannot free rows anymore */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[i]);
      if( consdata->addedcutsinit )
      {
         if( consdata->nodestoredcuts )
         {
            for( j = 0; j < consdata->naddedcuts; j++ )
            {
               SCIP_CALL( GCGreleaseMasterSepaCut(conshdlrdata->gcg, &(consdata->addedcuts[j])) );
            }
            SCIPfreeBlockMemoryArray(scip, &(consdata->addedcuts), consdata->naddedcuts);
            consdata->naddedcuts = 0;
         }
      }
      consdata->nodestoredcuts = FALSE;
      consdata->addedcutsinit = TRUE;
   }

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitMasterbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->nstack == 1 || SCIPgetNNodes(scip) == 0);

   /* free stack */
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->stack), conshdlrdata->maxstacksize);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->pendingvars), conshdlrdata->maxpendingbnds);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->pendingbndtypes), conshdlrdata->maxpendingbnds);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->pendingnewbnds), conshdlrdata->maxpendingbnds);
   SCIPhashmapFree(&conshdlrdata->pendingvarmaplb);
   SCIPhashmapFree(&conshdlrdata->pendingvarmapub);

   for( i = 0; i < conshdlrdata->maxblocknum; ++i )
   {
      if( conshdlrdata->maxcollectedbndvars[i] > 0)
      {
         SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->collectedbndvars[i]), conshdlrdata->maxcollectedbndvars[i]);
         SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->collectedlbnds[i]), conshdlrdata->maxcollectedbndvars[i]);
         SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->collectedubnds[i]), conshdlrdata->maxcollectedbndvars[i]);
      }
   }
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->collectedbndvars), conshdlrdata->maxblocknum);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->collectedlbnds), conshdlrdata->maxblocknum);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->collectedubnds), conshdlrdata->maxblocknum);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->ncollectedbndvars), conshdlrdata->maxblocknum);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->maxcollectedbndvars), conshdlrdata->maxblocknum);

   for( i = 0; i < conshdlrdata->maxlinkingvaridxs; ++i )
   {
      if( conshdlrdata->linkingvaridxs[i] != NULL)
         SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->linkingvaridxs[i]), conshdlrdata->maxblocknum);
      else
         break;
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(conshdlrdata->linkingvaridxs), conshdlrdata->maxlinkingvaridxs);

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->stack != NULL);

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert((consdata->node == NULL) == (consdata->parentcons == NULL));

   /* @fixme: This is a hack */
   if( consdata->node == NULL )
   {
      SCIPdebugMessage("root node not present in masterconsdata!\n");
      consdata->node = SCIPgetRootNode(scip);
   }

   assert(consdata->node != NULL);

   SCIPdebugMessage("Activating branch master constraint: <%s>[stack size: %d].\n", SCIPconsGetName(cons), conshdlrdata->nstack+1);
   /* If the node is activated the first time, we have to initialize the constraint data first */
   if( consdata->nactivated == 0 )
   {
      SCIPdebugPrintf("for the first time\n");
      SCIP_CALL( initializeConsdata(conshdlrdata->gcg, cons) );
   }
   else
      SCIPdebugPrintf("\n");

   consdata->nactivated++;

   /* The node has to be repropagated if new variables were created after the node was left the last time
    * or if new bound changes on directly transferred variables were found
    */
   assert(GCGmasterGetNPricedvars(conshdlrdata->gcg) >= consdata->npropvars);
   if( GCGmasterGetNPricedvars(conshdlrdata->gcg) > consdata->npropvars || consdata->ncopiedvarbnds > 0 )
   {
      consdata->needprop = TRUE;
      SCIP_CALL( SCIPrepropagateNode(scip, consdata->node) );
   }

   if( consdata->nlocalbndchgs - consdata->nlocalbndchgstreated[conshdlrdata->nstack] > 0 )
      SCIPdebugMessage("added %d boundchanges from previous nodes!\n", consdata->nlocalbndchgs - consdata->nlocalbndchgstreated[conshdlrdata->nstack]);

   /* put constraint on the stack */
   if( conshdlrdata->nstack >= conshdlrdata->maxstacksize )
   {
      int newsize = SCIPcalcMemGrowSize(scip,  conshdlrdata->nstack+1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(conshdlrdata->stack), conshdlrdata->maxstacksize, newsize) );

      SCIPdebugMessage("reallocating Memory for stack! %d --> %d\n", conshdlrdata->maxstacksize, newsize);
      conshdlrdata->maxstacksize = newsize;
   }

   conshdlrdata->stack[conshdlrdata->nstack] = cons;
   (conshdlrdata->nstack)++;

   SCIPdebugMessage("Activating masterbranch constraint: <%s> [stack size: %d], needprop = %u.\n",
      consdata->name, conshdlrdata->nstack, consdata->needprop);

   /* apply global bound changes in the original problem to the pricing problems */
   SCIP_CALL( applyGlobalBndchgsToPricingprobs(conshdlrdata->gcg, conshdlrdata) );

   /* apply local bound changes in the original problem to the pricing problems */
   SCIP_CALL( applyLocalBndchgsToPricingprobs(conshdlrdata->gcg, cons) );

   /* call branching specific activation method */
   if( consdata->branchrule != NULL )
   {
      SCIP_CALL( GCGrelaxBranchActiveMaster(conshdlrdata->gcg, consdata->branchrule, consdata->branchdata) );
   }

   /* forward history of node we are activating */
   forwardUpdateSeenHistory(conshdlrdata->gcg, consdata);
   /* forward history of possible ancestor nodes (all active) */
   SCIP_CONS* parentcons = consdata->parentcons;
   SCIP_CONSDATA* parentconsdata;
   while( parentcons != NULL )
   {
      parentconsdata = SCIPconsGetData(parentcons);
      GCGvarhistoryJumpToLatest(scip, &parentconsdata->knownvarhistory);
      parentcons = parentconsdata->parentcons;
   }
#ifndef MASTERSEP_DEBUG
   SCIPinfoMessage(scip, NULL, "Activation of node %lli of type %i with cons %s\n", SCIPnodeGetNumber(consdata->node),
                   SCIPnodeGetType(consdata->node), consdata->name);
#endif
   if( consdata->addedcutsinit )
   {
      assert(SCIPnodeGetType(consdata->node) != SCIP_NODETYPE_FOCUSNODE);
      assert(consdata->nactivated >= 1);
      SCIP_CALL( addStoredCutsToActiveCuts(conshdlrdata->gcg, consdata, conshdlrdata) );
   }

   /* if tree is currently probing, we do not clear generated cuts
    * - the cuts in generated cuts may not have been separated yet, as separation store gets switched when probing */
   if( SCIPnodeGetType(consdata->node) == SCIP_NODETYPE_FOCUSNODE )
      SCIP_CALL( GCGsepacutClearGeneratedCuts(conshdlrdata->gcg, conshdlrdata->eventhdlr) ); // potentially: move to pricing

   return SCIP_OKAY;
}

/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->stack != NULL || conshdlrdata->nstack == 1);
   assert(conshdlrdata->nstack > 0);
   assert(conshdlrdata->nstack == 1 || cons == conshdlrdata->stack[conshdlrdata->nstack-1]);

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nactivated >= 1);

   if( !conshdlrdata->pendingbndsactivated )
   {
      SCIPdebugMessage("We need repropagation\n");
      consdata->needprop = TRUE;
   }

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      consdata->npropvars = GCGmasterGetNPricedvars(conshdlrdata->gcg);

   /* remove constraint from the stack */
   (conshdlrdata->nstack)--;

   SCIPdebugMessage("Deactivating masterbranch constraint: <%s> [stack size: %d].\n",
      consdata->name, conshdlrdata->nstack);

   /* undo local bound changes in the original problem to the pricing problems */
   SCIP_CALL( undoLocalBndchgsToPricingprobs(conshdlrdata->gcg, cons) );

   /* call branching specific deactivation method */
   if( consdata->branchrule != NULL )
   {
      SCIP_CALL( GCGrelaxBranchDeactiveMaster(conshdlrdata->gcg, consdata->branchrule, consdata->branchdata) );
   }

   GCGvarhistoryJumpToLatest(scip, &consdata->knownvarhistory);
#ifndef MASTERSEP_DEBUG
   SCIPinfoMessage(scip, NULL, "Deactivation of node %lli of type %i with cons %s\n", SCIPnodeGetNumber(consdata->node),
                   SCIPnodeGetType(consdata->node), consdata->name);
#endif
   /* node is deactivated without any of its children being activated:
    *    - store the separator mastercuts which were (generated and) applied in this node in the data of its defining branch */
   if( !consdata->addedcutsinit )
      SCIP_CALL( initializeAddedCuts(conshdlrdata->gcg, consdata, conshdlrdata) );

   /* remove all the mastercuts (generated and) applied at this node from activecuts */
   SCIP_CALL( removeStoredCutsFromActiveCuts(conshdlrdata->gcg, consdata, conshdlrdata) );

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteMasterbranch)
{
   SCIP* origscip;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* parentconsdata;
   SCIP_CONSDATA** childconsdatas;
   SCIP_CONS** childconss;
   int nchildconss;
   int i;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   assert(cons != NULL);
   assert(consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(*consdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* get original problem */
   origscip = GCGgetOrigprob(conshdlrdata->gcg);
   assert(origscip != NULL);

   SCIPdebugMessage("Deleting masterbranch constraint: <%s>.\n", (*consdata)->name);

   /* remove original branching constraints if not yet done
    * (might happen if node is cut off before branching decisions are transferred to the original problem)
    */
   SCIP_CALL( GCGconsMasterbranchReleaseOrigbranchConss(conshdlrdata->gcg, cons) );

   /* free arrays with local bound changes on copied original variables */
   if( (*consdata)->maxcopiedvarbnds > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->copiedvars), (*consdata)->maxcopiedvarbnds);
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->copiedvarbndtypes), (*consdata)->maxcopiedvarbnds);
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->copiedvarbnds), (*consdata)->maxcopiedvarbnds);
   }

   /* free arrays with local bound changes on original variables belonging to a unique block */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->nlocalbndchgstreated, (*consdata)->maxlocalbndchgstreated);
   if( (*consdata)->maxlocalbndchgs > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->localoldbnds, (*consdata)->maxlocalbndchgs);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->localnewbnds, (*consdata)->maxlocalbndchgs);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->localbndtypes, (*consdata)->maxlocalbndchgs);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->localbndvars, (*consdata)->maxlocalbndchgs);
   }

   assert((*consdata)->origcons == NULL || GCGconsOrigbranchGetMastercons((*consdata)->origcons) == cons);

   /* allow the correspondig branchrule to delete the branch data */
   if( (*consdata)->branchdata != NULL && (*consdata)->branchrule != NULL )
   {
      SCIP_Bool force = ((*consdata)->origcons == NULL);
      SCIP_CALL( GCGrelaxBranchDataDelete(conshdlrdata->gcg, (*consdata)->branchrule, &(*consdata)->branchdata, FALSE, force) );
      if( (*consdata)->origcons != NULL && (*consdata)->branchdata == NULL )
         GCGconsOrigbranchSetBranchdata((*consdata)->origcons, NULL);
   }

   (*consdata)->branchdata = NULL;

   /* set the mastercons pointer of the corresponding origcons to NULL */
   if( (*consdata)->origcons != NULL )
   {
      GCGconsOrigbranchSetMastercons((*consdata)->origcons, NULL);
   }

   assert((*consdata)->knownvarhistory != NULL);
   SCIP_CALL( GCGvarhistoryFreeReference(scip, &(*consdata)->knownvarhistory) );

   /* remove branching constraints at child nodes */
   nchildconss = (*consdata)->nchildconss;
   if( nchildconss > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &childconsdatas, nchildconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &childconss, nchildconss) );

      for( i = 0; i < nchildconss; ++i )
      {
         if( (*consdata)->childconss != NULL && (*consdata)->childconss[i] != NULL )
         {
            childconsdatas[i] = SCIPconsGetData((*consdata)->childconss[i]);
            childconss[i] = (*consdata)->childconss[i];
         }
         else
         {
            childconsdatas[i] = NULL;
            childconss[i] = NULL;
         }
      }

      /* delete childnodes */
      for( i = 0; i < nchildconss; ++i )
      {
         SCIPdebugMessage("Deleting %d childnodes\n", nchildconss);

         if( childconss[i] != NULL )
         {
            /*SCIP_CALL( consDeleteMasterbranch(scip, conshdlr, childcons[i], &childconsdatas[i]) );*/
            SCIP_CALL( SCIPreleaseCons(scip, &childconss[i]) );
            childconss[i] = NULL;
         }
      }

      SCIPfreeBufferArray(scip, &childconss);
      SCIPfreeBufferArray(scip, &childconsdatas);
   }
   assert((*consdata)->nchildconss == 0);

   /* set the pointer in the parent node to NULL */
   if( (*consdata)->parentcons != NULL )
   {
      SCIP_Bool isinprobing;
#ifndef NDEBUG
      SCIP_Bool childdeleted = FALSE;
#endif

      parentconsdata = SCIPconsGetData((*consdata)->parentcons);

      isinprobing = (SCIPgetStage(scip) <= SCIP_STAGE_SOLVING && SCIPinProbing(scip)) || (SCIPgetStage(origscip) <= SCIP_STAGE_SOLVING && SCIPinProbing(origscip));
      if( isinprobing )
         parentconsdata->probingtmpcons = NULL;

      for( i = 0; i < parentconsdata->nchildconss; ++i )
      {
         if( parentconsdata->childconss[i] == cons )
         {
            parentconsdata->childconss[i] = parentconsdata->childconss[parentconsdata->nchildconss-1];

            parentconsdata->childconss[parentconsdata->nchildconss-1] = NULL;
            parentconsdata->nchildconss--;
#ifndef NDEBUG
            childdeleted = TRUE;
#endif
            break;
         }
      }
      assert(childdeleted || isinprobing);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->childconss, (*consdata)->maxchildconss);

   if( (*consdata)->name != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->name, strlen((*consdata)->name)+1);
   }

   /* we delete all the data associated with the mastercuts from the constraint data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if(  (*consdata)->addedcutsinit )
   {
#ifndef MASTERSEP_DEBUG
      SCIPinfoMessage(scip, NULL, "delete master: node stored cuts %i\n", (*consdata)->nodestoredcuts);
#endif
      if( (*consdata)->nodestoredcuts )
      {
#ifndef MASTERSEP_DEBUG
         SCIPinfoMessage(scip, NULL, "delete master: release %i cuts and free addedcuts\n", i, (*consdata)->naddedcuts);
#endif
         for( j = 0; j < (*consdata)->naddedcuts; j++ )
         {
            SCIP_CALL( GCGreleaseMasterSepaCut(conshdlrdata->gcg, &((*consdata)->addedcuts[j])) );
         }
         SCIPfreeBlockMemoryArray(scip, &((*consdata)->addedcuts), (*consdata)->naddedcuts);
         (*consdata)->naddedcuts = 0;
      }
   }

   SCIPfreeBlockMemoryNull(scip, consdata);
   *consdata = NULL;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropMasterbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;

   int propcount;

   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->stack != NULL);

   *result = SCIP_DIDNOTRUN;

   /* the constraint data of the cons related to the current node */
   cons = conshdlrdata->stack[conshdlrdata->nstack-1];
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( conshdlrdata->npendingbnds == 0 && !consdata->needprop && consdata->ncopiedvarbnds == 0 )
   {
      SCIPdebugMessage("No propagation of masterbranch constraint needed: <%s>, stack size = %d.\n",
         consdata->name, conshdlrdata->nstack);

      return SCIP_OKAY;
   }

   SCIPdebugMessage("Starting propagation of masterbranch constraint: <%s>, stack size = %d, newvars = %d, npendingbnds = %d, npropbounds = %d.\n",
      consdata->name, conshdlrdata->nstack, GCGmasterGetNPricedvars(conshdlrdata->gcg) - consdata->npropvars, conshdlrdata->npendingbnds, consdata->ncopiedvarbnds);

   propcount = 0;

   if( conshdlrdata->npendingbnds > 0 )
   {
      *result = SCIP_DIDNOTFIND;

      if( !conshdlrdata->pendingbndsactivated )
      {
         /* apply global bound changes in the original problem to the pricing problems */
         SCIP_CALL(applyGlobalBndchgsToPricingprobs(conshdlrdata->gcg, conshdlrdata));
      }

      /* apply global bound changes on original problem variables to the master problem */
      SCIP_CALL( applyGlobalBndchgsToPricedMastervars(conshdlrdata->gcg, conshdlrdata, &propcount) );

      GCGcolpoolPropagateGlobalBounds(GCGgetColpool(conshdlrdata->gcg));
   }

   if( consdata->needprop || consdata->ncopiedvarbnds != 0 )
   {
      *result = SCIP_DIDNOTFIND;

      /* apply local bound changes on the original variables on newly generated master variables */
      SCIP_CALL( applyLocalBndchgsToPricedMastervars(conshdlrdata->gcg, conshdlrdata, cons, &propcount) );

      /* apply local bound changes on original variables that have been directly copied to the master problem */
      SCIP_CALL( applyLocalBndchgsToCopiedMastervars(conshdlrdata->gcg, cons, &propcount) );

      /* call branching rule specific propagation method */
      if ( consdata->branchrule != NULL )
      {
         /** @todo count number of propagations */
         SCIP_CALL( GCGrelaxBranchPropMaster(conshdlrdata->gcg, consdata->branchrule, consdata->branchdata, result) );
      }

      consdata->needprop = FALSE;
      consdata->npropvars = GCGmasterGetNPricedvars(conshdlrdata->gcg);
   }

   if( *result != SCIP_CUTOFF && propcount > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** enforcement method for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpMasterbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** enforcement method for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsMasterbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** check method for solutions */
static
SCIP_DECL_CONSCHECK(consCheckMasterbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** variable lock method */
static
SCIP_DECL_CONSLOCK(consLockMasterbranch)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/*
 * Callback methods of event handler
 */

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolOrigvarbound)
{  /*lint --e{715}*/
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

/** eventhdlr destructor */
static
SCIP_DECL_EVENTFREE(eventFreeOrigvarbound)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecOrigvarbound)
{  /*lint --e{715}*/
   SCIP* masterscip;
   SCIP_EVENTTYPE eventtype;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_VAR* var;
   SCIP_Real oldbound;
   SCIP_Real newbound;
   int blocknr;
#ifdef SCIP_DEBUG
   SCIP_Bool handled = FALSE;
#endif
   SCIP_VAR** mastervars;
#ifndef NDEBUG
   int nmastervars;
   SCIP_Real* mastervals;
#endif

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* get master problem */
   masterscip = GCGgetMasterprob(eventhdlrdata->gcg);
   assert(masterscip != NULL);

   /* get event data */
   eventtype = SCIPeventGetType(event);
   var = SCIPeventGetVar(event);
   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);

   SCIPdebugMessage("eventexec: eventtype = 0x%x, var = %s, oldbound = %f, newbound = %f\n", (unsigned int) eventtype, SCIPvarGetName(var), oldbound, newbound);

   if( GCGgetDecompositionMode(eventhdlrdata->gcg) != GCG_DECMODE_DANTZIGWOLFE || !GCGrelaxIsInitialized(eventhdlrdata->gcg) )
   {
//      assert(SCIPvarGetData(var) == NULL);
      SCIPdebugMessage("Ignoring since in presolving / propagating.\n");
      return SCIP_OKAY;
   }

   if( !SCIPisTransformed(masterscip) )
   {
      GCGinitializeMasterProblemSolve(eventhdlrdata->gcg);
   }

   assert(GCGvarIsOriginal(var));
   blocknr = GCGvarGetBlock(var);

   mastervars = GCGoriginalVarGetMastervars(var);
#ifndef NDEBUG
   nmastervars = GCGoriginalVarGetNMastervars(var);
   mastervals = GCGoriginalVarGetMastervals(var);
#endif

   assert(SCIPgetStage(masterscip) >= SCIP_STAGE_TRANSFORMED);

   /* A global bound change might turn the current relaxation solution invalid */
   if( SCIPisRelaxSolValid(scip)
      && (((eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 && SCIPisFeasLT(scip, SCIPgetRelaxSolVal(scip, var), newbound))
      || ((eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0 && SCIPisFeasGT(scip, SCIPgetRelaxSolVal(scip, var), newbound))) )
   {
      SCIP_CALL( SCIPmarkRelaxSolInvalid(scip) );
   }

   /* deal with variables present in the pricing */
   if( blocknr >= 0 && GCGisPricingprobRelevant(eventhdlrdata->gcg, blocknr) )
   {
      SCIPdebugMessage("Pricing var!\n");
      if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( addPendingBndChg(masterscip, GCGoriginalVarGetPricingVar(var), SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );
      }
      if( (eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( addPendingBndChg(masterscip, GCGoriginalVarGetPricingVar(var), SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );
      }
   }
   /* deal with variables appearing in the master only */
   if( blocknr == -1 )
   {
      assert(nmastervars == 1);
      assert(mastervals[0] == 1);
      assert(mastervars[0] != NULL);

      SCIPdebugMessage("Master var!\n");
      if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( addPendingBndChg(masterscip, mastervars[0], SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );
      }
      if( (eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( addPendingBndChg(masterscip, mastervars[0], SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );
      }
   }
   /* deal with linking variables */
   if( blocknr == -2 )
   {
      int npricingprobs;
      SCIP_VAR** pricingvars;
      int i;

      SCIPdebugMessage("Linking var!\n");
      pricingvars = GCGlinkingVarGetPricingVars(var);
      npricingprobs = GCGgetNPricingprobs(eventhdlrdata->gcg);

      assert(nmastervars >= 1);
      assert(mastervals[0] == 1);
      assert(mastervars[0] != NULL);
      assert(GCGvarGetBlock(mastervars[0]) == -1);

      if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
      {
#ifdef SCIP_DEBUG
            handled = TRUE;
#endif
         /* add the bound change in the master */
         SCIP_CALL( addPendingBndChg(masterscip, mastervars[0], SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );

         /* add the bound change to the pricing problems */
         for( i = 0; i < npricingprobs; ++i )
         {
            if( pricingvars[i] == NULL )
               continue;
#ifdef SCIP_DEBUG
            handled = TRUE;
#endif
            SCIP_CALL( addPendingBndChg(masterscip, pricingvars[i], SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );
         }
      }
      if( (eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0 )
      {
#ifdef SCIP_DEBUG
            handled = TRUE;
#endif
         /* add the bound change in the master */
         SCIP_CALL( addPendingBndChg(masterscip, mastervars[0], SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );

         /* add the bound change to the pricing problems */
         for( i = 0; i < npricingprobs; ++i )
         {
            if( pricingvars[i] == NULL )
               continue;
#ifdef SCIP_DEBUG
            handled = TRUE;
#endif
            SCIP_CALL( addPendingBndChg(masterscip, pricingvars[i], SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );
         }

      }

   }
#ifdef SCIP_DEBUG
   if( !handled )
   {
      SCIPdebugMessage("Effectively ignoring this change\n");
   }
#endif

   return SCIP_OKAY;
}


/*
 * interface methods
 */


/** creates the handler for masterbranch constraints and includes it in SCIP */
SCIP_RETCODE GCGincludeConshdlrMasterbranch(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* masterprob;
   SCIP* origscip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);

   /* get original problem */
   origscip = GCGgetOrigprob(gcg);
   assert(origscip != NULL);

   SCIP_CALL( SCIPallocMemory(masterprob, &conshdlrdata) );
   conshdlrdata->gcg = gcg;
   conshdlrdata->stack = NULL;
   conshdlrdata->nstack = 0;
   conshdlrdata->maxstacksize = 25;
   conshdlrdata->eventhdlr = NULL;
   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(masterprob, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpMasterbranch, consEnfopsMasterbranch, consCheckMasterbranch, consLockMasterbranch,
         conshdlrdata) );
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPsetConshdlrFree(masterprob, conshdlr, consFreeMasterbranch) );
   SCIP_CALL( SCIPsetConshdlrInit(masterprob, conshdlr, consInitMasterbranch) );
   SCIP_CALL( SCIPsetConshdlrExit(masterprob, conshdlr, consExitMasterbranch) );
   SCIP_CALL( SCIPsetConshdlrInitsol(masterprob, conshdlr, consInitsolMasterbranch) );
   SCIP_CALL( SCIPsetConshdlrExitsol(masterprob, conshdlr, consExitSolMasterbranch) );
   SCIP_CALL( SCIPsetConshdlrDelete(masterprob, conshdlr, consDeleteMasterbranch) );
   SCIP_CALL( SCIPsetConshdlrActive(masterprob, conshdlr, consActiveMasterbranch) );
   SCIP_CALL( SCIPsetConshdlrDeactive(masterprob, conshdlr, consDeactiveMasterbranch) );
   SCIP_CALL( SCIPsetConshdlrProp(masterprob, conshdlr, consPropMasterbranch, CONSHDLR_PROPFREQ,
         CONSHDLR_DELAYPROP, CONSHDLR_PROPTIMING) );

   /* create event handler data */
   SCIP_CALL( SCIPallocBlockMemory(origscip, &eventhdlrdata) );
   eventhdlrdata->gcg = gcg;

   /* include event handler into original SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(origscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecOrigvarbound, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrInitsol(origscip, eventhdlr, eventInitsolOrigvarbound) );
   SCIP_CALL( SCIPsetEventhdlrFree(origscip, eventhdlr, eventFreeOrigvarbound) );

   SCIP_CALL( SCIPaddBoolParam(origscip, "relaxing/gcg/enforceproper",
         "should propagated bound changes in the original be enforced in the master (only proper vars)?",
         &conshdlrdata->enforceproper, FALSE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a masterbranch constraint */
SCIP_RETCODE GCGcreateConsMasterbranch(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of the constraint */
   SCIP_NODE*            node,               /**< node at which the constraint should be created */
   SCIP_CONS*            parentcons,         /**< parent constraint */
   SCIP_BRANCHRULE*      branchrule,         /**< pointer to the branching rule */
   GCG_BRANCHDATA*       branchdata,         /**< branching data */
   SCIP_CONS**           origbranchconss,    /**< original constraints enforcing the branching decision */
   int                   norigbranchconss,   /**< number of original constraints */
   int                   maxorigbranchconss  /**< capacity of origbranchconss */
   )
{
   SCIP* masterprob;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);
   assert(node != NULL || parentcons == NULL);
   if( node != NULL )
      assert((parentcons == NULL) == (SCIPnodeGetDepth(node) == 0));
   else
      assert(parentcons == NULL);

   assert(SCIPgetStage(masterprob) <= SCIP_STAGE_SOLVING);

   /* find the masterbranch constraint handler */
   conshdlr = SCIPfindConshdlr(masterprob, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(masterprob, &consdata) );

   SCIP_CALL( SCIPduplicateBlockMemoryArray(masterprob, &(consdata->name), name, strlen(name)+1) );

   consdata->npropvars = 0;
   consdata->needprop = TRUE;
   consdata->node = node;
   consdata->nactivated = 0;

   consdata->parentcons = parentcons;
   consdata->childconss = NULL;
   consdata->nchildconss = 0;
   consdata->maxchildconss = 0;
   consdata->probingtmpcons = NULL;
   consdata->origcons = NULL;

   consdata->branchdata = branchdata;
   consdata->branchrule = branchrule;

   consdata->knownvarhistory = NULL;
   SCIP_CALL( GCGvarhistoryCopyReference(masterprob, &consdata->knownvarhistory, GCGgetCurrentVarhistoryReference(gcg)) );

   consdata->localbndvars = NULL;
   consdata->localbndtypes = NULL;
   consdata->localnewbnds = NULL;
   consdata->localoldbnds = NULL;
   consdata->nlocalbndchgstreated = NULL;
   consdata->maxlocalbndchgstreated = 0;
   consdata->nlocalbndchgs = 0;
   consdata->maxlocalbndchgs = 0;
   consdata->nbranchingchgs = 0;

   consdata->copiedvars = NULL;
   consdata->copiedvarbndtypes = NULL;
   consdata->copiedvarbnds = NULL;
   consdata->ncopiedvarbnds = 0;
   consdata->maxcopiedvarbnds = 0;

   consdata->maxorigbranchconss = maxorigbranchconss;
   consdata->origbranchconss = origbranchconss;
   consdata->norigbranchconss = norigbranchconss;

   /* initialize structures for managing separator generated cuts for master problem*/
   consdata->addedcuts = NULL;
   consdata->naddedcuts = 0;
   consdata->firstnewcut = 0;
   consdata->nodestoredcuts = FALSE;
   consdata->addedcutsinit = FALSE;

   SCIPdebugMessage("Creating masterbranch constraint with parent %p.\n", (void*) parentcons);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(masterprob, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
         TRUE, FALSE, FALSE, FALSE, TRUE) );

   /* add the new masterbranch constraint to the parent node's masterbranch constraint data
    * (unless the current node is the root node)
    */
   if( parentcons != NULL )
   {
      SCIP_CONSDATA* parentdata;

      parentdata = SCIPconsGetData(parentcons);
      assert(parentdata != NULL);

      if( SCIPinProbing(masterprob) || SCIPinProbing(GCGgetOrigprob(gcg)) )
      {
         parentdata->probingtmpcons = *cons;
      }
      else
      {
         int newmaxsize;
         ++parentdata->nchildconss;
         newmaxsize = SCIPcalcMemGrowSize(masterprob, parentdata->nchildconss);
         if( parentdata->nchildconss == 1 )
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &(parentdata->childconss), newmaxsize) );
            parentdata->childconss[0] = NULL;
         }
         else
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &(parentdata->childconss), parentdata->maxchildconss, newmaxsize) );
            parentdata->childconss[parentdata->nchildconss - 1] = NULL;
         }
         parentdata->maxchildconss = newmaxsize;

         assert(parentdata->childconss[parentdata->nchildconss - 1] == NULL);
         parentdata->childconss[parentdata->nchildconss - 1] = *cons;

         // stash limit settings until branching is applied to the original problem
         SCIP_CALL( GCGstashLimitSettings(gcg) );
      }
   }

   return SCIP_OKAY;
}

/** returns the name of the constraint */
char* GCGconsMasterbranchGetName(
   SCIP_CONS*            cons                /**< masterbranch constraint for which the data is requested */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->name;
}

/** returns the node in the B&B tree at which the given masterbranch constraint is sticking */
SCIP_NODE* GCGconsMasterbranchGetNode(
   SCIP_CONS*            cons                /**< constraint pointer */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->node;
}

/** returns the masterbranch constraint of the B&B father of the node at which the
  * given masterbranch constraint is sticking
  */
SCIP_CONS* GCGconsMasterbranchGetParentcons(
   SCIP_CONS*            cons                /**< constraint pointer */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->parentcons;
}

/** returns the number of masterbranch constraints of the children of the node at which the
  * given masterbranch constraint is sticking
  */
int GCGconsMasterbranchGetNChildconss(
   SCIP_CONS*            cons                /**< constraint pointer */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nchildconss;
}

/** returns a masterbranch constraint of a child of the node at which the
  * given masterbranch constraint is sticking
  */
SCIP_CONS* GCGconsMasterbranchGetChildcons(
   SCIP_CONS*            cons,                /**< constraint pointer */
   int                   childnr              /**< index of the child node */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->childconss != NULL);

   assert(consdata->nchildconss > childnr);

   return consdata->childconss[childnr];
}

/** returns the origbranch constraint of the node in the original program corresponding to the node
  * at which the given masterbranch constraint is sticking
  */
SCIP_CONS* GCGconsMasterbranchGetOrigcons(
   SCIP_CONS*            cons                /**< constraint pointer */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->origcons;
}

/** sets the origbranch constraint of the node in the master program corresponding to the node
  * at which the given masterbranchbranch constraint is sticking
  */
void GCGconsMasterbranchSetOrigcons(
   SCIP_CONS*            cons,               /**< constraint pointer */
   SCIP_CONS*            origcons            /**< original branching constraint */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->origcons == NULL || origcons == NULL);

   consdata->origcons = origcons;
}

/** returns the branching data for a given masterbranch constraint */
GCG_BRANCHDATA* GCGconsMasterbranchGetBranchdata(
   SCIP_CONS*            cons                /**< constraint pointer */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->branchdata;
}

void GCGconsMasterbranchSetBranchdata(
   SCIP_CONS*            cons,               /**< masterbranch constraint for which the branching data is requested */
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->branchdata = branchdata;
}

/** returns the branching rule of the constraint */
SCIP_BRANCHRULE* GCGconsMasterbranchGetBranchrule(
   SCIP_CONS*            cons                /**< masterbranch constraint for which the data is requested */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->branchrule;
}

/** adds a bound change on an original variable that was directly copied to the master problem */
SCIP_RETCODE GCGconsMasterbranchAddCopiedVarBndchg(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS*            cons,               /**< masterbranch constraint to which the bound change is added */
   SCIP_VAR*             var,                /**< variable on which the bound change was performed */
   GCG_BOUNDTYPE         boundtype,          /**< bound type of the bound change */
   SCIP_Real             newbound            /**< new bound of the variable after the bound change */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP* scip;

   assert(gcg != NULL);
   
   scip = GCGgetMasterprob(gcg);
   assert(GCGisMaster(scip));
   assert(cons != NULL);
   assert(var != NULL);

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* realloc the arrays, if needed */
   if( consdata->ncopiedvarbnds >= consdata->maxcopiedvarbnds )
   {
      int newsize = SCIPcalcMemGrowSize(scip, consdata->ncopiedvarbnds+1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->copiedvars), consdata->maxcopiedvarbnds, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->copiedvarbndtypes), consdata->maxcopiedvarbnds, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->copiedvarbnds), consdata->maxcopiedvarbnds, newsize) );
      consdata->maxcopiedvarbnds = newsize;
   }

   SCIPdebugMessage("Bound change on copied original variable stored at masterbranch constraint: <%s>.\n", SCIPconsGetName(cons));

   /* store the new bound change */
   consdata->copiedvars[consdata->ncopiedvarbnds] = var;
   consdata->copiedvarbndtypes[consdata->ncopiedvarbnds] = boundtype;
   consdata->copiedvarbnds[consdata->ncopiedvarbnds] = newbound;
   consdata->ncopiedvarbnds++;

   /* mark the corresponding master node to be repropagated */
   SCIP_CALL( SCIPrepropagateNode(scip, consdata->node) );

   return SCIP_OKAY;
}

/** returns the constraints in the original problem that enforce the branching decision */
SCIP_CONS** GCGconsMasterbranchGetOrigbranchConss(
   SCIP_CONS*            cons                /**< masterbranch constraint for which the data is requested */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->origbranchconss;
}

/** returns the number of constraints in the original problem that enforce the branching decision */
int GCGconsMasterbranchGetNOrigbranchConss(
   SCIP_CONS*            cons                /**< masterbranch constraint for which the data is requested */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->norigbranchconss;
}

/** releases the constraints in the original problem that enforce the branching decision
 *  and frees the array holding the constraints
 */
SCIP_RETCODE GCGconsMasterbranchReleaseOrigbranchConss(
   GCG*                  gcg,                /**< GCG instance */
   SCIP_CONS*            cons                /**< masterbranch constraint for which the data is freed */
   )
{
   SCIP* origscip;
   SCIP* masterscip;
   SCIP_CONSDATA* consdata;

   masterscip = GCGgetMasterprob(gcg);
   origscip = GCGgetOrigprob(gcg);
   assert(GCGisMaster(masterscip));
   assert(GCGisOriginal(origscip));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->norigbranchconss > 0 )
   {
      int c;

      for( c = consdata->norigbranchconss - 1; c >= 0; --c )
      {
         SCIP_CALL( SCIPreleaseCons(origscip, &consdata->origbranchconss[c]) );
      }
      SCIPfreeBlockMemoryArray(masterscip, &consdata->origbranchconss, consdata->maxorigbranchconss);
      consdata->origbranchconss = NULL;
      consdata->norigbranchconss = 0;
      consdata->maxorigbranchconss = 0;
   }

   return SCIP_OKAY;
}

/** returns the masterbranch constraint of the current node */
SCIP_CONS* GCGconsMasterbranchGetActiveCons(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(gcg != NULL);

   scip = GCGgetMasterprob(gcg);
   assert(GCGisMaster(scip));

   /* get constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->stack != NULL);

   if( conshdlrdata->nstack == 0 || SCIPgetStage(scip) > SCIP_STAGE_SOLVING )
      return NULL;

   assert(conshdlrdata->stack[conshdlrdata->nstack-1] != NULL);
   return conshdlrdata->stack[conshdlrdata->nstack-1];
}

/** returns the stack and the number of elements on it */
void GCGconsMasterbranchGetStack(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(gcg != NULL);
   scip = GCGgetMasterprob(gcg);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->stack != NULL);

   *stack = conshdlrdata->stack;
   *nstackelements = conshdlrdata->nstack;
}

/** returns the number of elements on the stack */
int GCGconsMasterbranchGetNStackelements(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(gcg != NULL);
   scip = GCGgetMasterprob(gcg);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->stack != NULL);

   return conshdlrdata->nstack;
}

/** adds initial constraint to root node */
SCIP_RETCODE GCGconsMasterbranchAddRootCons(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;

   assert(gcg != NULL);
   scip = GCGgetMasterprob(gcg);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   cons = conshdlrdata->stack[0];
   conshdlrdata->stack[0] = NULL;
   assert(conshdlrdata->nstack == 1);
   conshdlrdata->nstack = 0;

   SCIP_CALL( SCIPaddConsNode(scip, SCIPgetRootNode(scip), cons, SCIPgetRootNode(scip)) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/** check whether the node was generated by generic branching */
SCIP_Bool GCGcurrentNodeIsGeneric(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_CONS* mastercons;
   SCIP_BRANCHRULE* branchrule;

   assert(gcg != NULL);

   mastercons = GCGconsMasterbranchGetActiveCons(gcg);

   /* @todo: Why might mastercons be NULL? */
   if( mastercons == NULL || SCIPnodeGetDepth(GCGconsMasterbranchGetNode(mastercons)) == 0 )
      return FALSE;

   branchrule = GCGconsMasterbranchGetBranchrule(mastercons);

   if( branchrule == NULL || strcmp(SCIPbranchruleGetName(branchrule), "generic") != 0 )
      return FALSE;

   return TRUE;
}

/** checks the consistency of the masterbranch constraints in the problem */
void GCGconsMasterbranchCheckConsistency(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR*     conshdlr;
   int nconss;
#ifndef NDEBUG
   int i;
   SCIP_CONS** conss;
#endif

   assert(gcg != NULL);
   scip = GCGgetMasterprob(gcg);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNConss(conshdlr);
#ifndef NDEBUG
   conss = SCIPconshdlrGetConss(conshdlr);

   for( i = 0; i < nconss; i++ )
   {
      SCIP_CONSDATA* consdata;
#ifdef SCIP_DEBUG
      SCIP_CONS* parent_origcons;
      SCIP_CONS* origcons_parent;
#endif

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      assert(consdata->probingtmpcons == NULL || SCIPinProbing(scip));
      assert(consdata->probingtmpcons == NULL || SCIPconsGetData(consdata->probingtmpcons)->parentcons == conss[i]);
      assert(consdata->origcons == NULL || GCGconsOrigbranchGetMastercons(consdata->origcons) == conss[i]);

#ifdef SCIP_DEBUG
      if( consdata->parentcons != NULL )
         parent_origcons = SCIPconsGetData(consdata->parentcons)->origcons;
      else
         parent_origcons = NULL;

      if( consdata->origcons != NULL )
         origcons_parent = GCGconsOrigbranchGetParentcons(consdata->origcons);
      else
         origcons_parent = NULL;

      SCIPdebugMessage("cons: %s (node %p), origcons: %s, parent %s: %s => %s\n",
         SCIPconsGetName(conss[i]),
         (void*) consdata->node,
         consdata->origcons == NULL? "NULL" : SCIPconsGetName(consdata->origcons),
         consdata->parentcons == NULL? "NULL" : SCIPconsGetName(consdata->parentcons),
         parent_origcons == NULL? "NULL" :  SCIPconsGetName(parent_origcons),
         origcons_parent == NULL? "NULL" : SCIPconsGetName(origcons_parent) );
#endif
   }
#endif

   SCIPdebugMessage("checked consistency of %d masterbranch constraints, all ok!\n", nconss);
}
