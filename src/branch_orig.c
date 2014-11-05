/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2014 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_orig.c
 * @brief  branching rule for the original problem in GCG
 * @author Gerald Gamrath
 * @author Marcel Schmickerath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */
#include <assert.h>
#include <string.h>

#include "branch_orig.h"
#include "gcg.h"
#include "branch_relpsprob.h"
#include "cons_integralorig.h"
#include "cons_masterbranch.h"
#include "cons_origbranch.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "type_branchgcg.h"

#include "scip/cons_linear.h"


#define BRANCHRULE_NAME          "orig"
#define BRANCHRULE_DESC          "branching for the original program in generic column generation"
#define BRANCHRULE_PRIORITY      100
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define DEFAULT_ENFORCEBYCONS FALSE
#define DEFAULT_MOSTFRAC      FALSE
#define DEFAULT_USEPSEUDO     TRUE
#define DEFAULT_USEPSSTRONG   FALSE

/** branching data for branching decisions */
struct GCG_BranchData
{
   SCIP_VAR*             origvar;            /**< original variable on which the branching is done */
   GCG_BOUNDTYPE         boundtype;          /**< type of the new bound of original variable */
   SCIP_Real             newbound;           /**< new lower/upper bound of the original variable */
   SCIP_Real             oldbound;           /**< old lower/upper bound of the pricing variable */
   SCIP_Real             oldvalue;           /**< old value of the original variable */
   SCIP_Real             olddualbound;       /**< dual bound before the branching was performed */
   SCIP_CONS*            cons;               /**< constraint that enforces the branching restriction in the original
                                              *   problem, or NULL if this is done by variable bounds */
};

/** updates extern branching candidates before branching */
SCIP_RETCODE updateExternBranchcandsForMasterbranch(
   SCIP*                 scip               /**< SCIP data structure */
)
{
   SCIP_VAR** origvars;
   int norigvars;
   int i;

   assert(GCGisOriginal(scip));

   origvars = SCIPgetVars(scip);
   norigvars = SCIPgetNVars(scip);
   assert(origvars != NULL);

   SCIPclearExternBranchCands(scip);

   /* store branching candidates */
   for( i = 0; i < norigvars; i++ )
   {
      if( SCIPvarGetType(origvars[i]) <= SCIP_VARTYPE_INTEGER && !SCIPisFeasIntegral(scip, SCIPgetRelaxSolVal(scip, origvars[i])) )
      {
         assert(!SCIPisEQ(scip, SCIPvarGetLbLocal(origvars[i]), SCIPvarGetUbLocal(origvars[i])));

         SCIP_CALL( SCIPaddExternBranchCand(scip, origvars[i], SCIPgetRelaxSolVal(scip,
            origvars[i]) - SCIPfloor(scip, SCIPgetRelaxSolVal(scip, origvars[i])),
            SCIPgetRelaxSolVal(scip, origvars[i])) );
      }
   }
   SCIPdebugMessage("updated relaxation branching candidates\n");

   return SCIP_OKAY;
}

/** branches on a given variable */
static
SCIP_RETCODE branchVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< pointer to the original variable branching rule */
   SCIP_VAR*             branchvar,          /**< variable to branch on */
   SCIP_Real             solval              /**< value of the variable in the current solution */
   )
{
   /* data for b&b child creation */
   SCIP* masterscip;
   SCIP_NODE* child1;
   SCIP_NODE* child2;
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;
   SCIP_CONS** origbranchconss1;
   SCIP_CONS** origbranchconss2;
   GCG_BRANCHDATA* branchupdata;
   GCG_BRANCHDATA* branchdowndata;
   SCIP_Bool propagatebndchg;
   char upname[SCIP_MAXSTRLEN];
   char downname[SCIP_MAXSTRLEN];
   int norigbranchconss;

   /* parameter data */
   SCIP_Bool enforcebycons;

   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(branchvar != NULL);

   /* get master problem */
   masterscip = GCGgetMasterprob(scip);
   assert(masterscip != NULL);

   propagatebndchg = FALSE;
   origbranchconss1 = NULL;
   origbranchconss2 = NULL;
   norigbranchconss = 0;

   /* for cons_masterbranch */

   /*  create two child-nodes of the current node in the b&b-tree and add the masterbranch constraints */
   SCIP_CALL( SCIPcreateChild(masterscip, &child1, 0.0, SCIPgetLocalTransEstimate(masterscip)) );
   SCIP_CALL( SCIPcreateChild(masterscip, &child2, 0.0, SCIPgetLocalTransEstimate(masterscip)) );
   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &cons1, child1, GCGconsMasterbranchGetActiveCons(masterscip)) );
   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &cons2, child2, GCGconsMasterbranchGetActiveCons(masterscip)) );
   SCIP_CALL( SCIPaddConsNode(masterscip, child1, cons1, NULL) );
   SCIP_CALL( SCIPaddConsNode(masterscip, child2, cons2, NULL) );

   /* get values of parameters */
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/enforcebycons", &enforcebycons) );

   SCIPdebugMessage("Branching on var %s with value %g in current solution\n", SCIPvarGetName(branchvar), solval);

   /** @todo use block memory here */
   /* create the branch data for the children and assign the values */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchupdata) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchdowndata) );

   branchupdata->origvar = branchvar;
   branchupdata->oldvalue = solval;
   branchupdata->olddualbound = SCIPgetLocalLowerbound(GCGgetMasterprob(scip));
   branchupdata->boundtype = GCG_BOUNDTYPE_LOWER;
   branchupdata->newbound = SCIPceil(scip, solval);
   branchupdata->oldbound = SCIPvarGetLbLocal(branchvar);
   branchupdata->cons = cons1;

   branchdowndata->origvar = branchvar;
   branchdowndata->oldvalue = solval;
   branchdowndata->olddualbound = SCIPgetLocalLowerbound(GCGgetMasterprob(scip));
   branchdowndata->boundtype = GCG_BOUNDTYPE_UPPER;
   branchdowndata->newbound = SCIPfloor(scip, solval);
   branchdowndata->oldbound = SCIPvarGetUbLocal(branchvar);
   branchdowndata->cons = cons2;


   (void) SCIPsnprintf(upname, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchupdata->origvar),
      ">=", branchupdata->newbound);
   (void) SCIPsnprintf(downname, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchdowndata->origvar),
      "<=", branchdowndata->newbound);

   /* enforce branching decision by a constraint rather than by bound changes */
   if( enforcebycons )
   {
      /* enforce new bounds by linear constraints */
      SCIP_CONS* consup;
      SCIP_CONS* consdown;

      norigbranchconss = 1;

      SCIPdebugMessage("enforced by cons\n");

      SCIP_CALL( GCGconsOrigbranchCreateOrigconsArray(masterscip, &origbranchconss1, norigbranchconss) );
      SCIP_CALL( GCGconsOrigbranchCreateOrigconsArray(masterscip, &origbranchconss2, norigbranchconss) );

      /* create corresponding constraints */
      SCIP_CALL( SCIPcreateConsLinear(scip, &consup, upname, 0, NULL, NULL,
            SCIPceil(scip, solval), SCIPinfinity(scip),
            TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &consdown, downname, 0, NULL, NULL,
            -1.0 * SCIPinfinity(scip), SCIPfloor(scip, solval),
            TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPaddCoefLinear(scip, consup, branchvar, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, consdown, branchvar, 1.0) );

      origbranchconss1[0] = consup;
      origbranchconss2[0] = consdown;

      branchupdata->cons = consup;
      branchdowndata->cons = consdown;
   }
   else
   {
      branchupdata->cons = NULL;
      branchdowndata->cons = NULL;
   }

   /* store bound change of variables that were directly transferred to the master problem */
   if( !enforcebycons && GCGvarGetBlock(branchvar) == -1 )
      propagatebndchg = TRUE;

   SCIP_CALL( GCGconsMasterbranchSetOrigConsData(masterscip, cons1, upname, branchrule,
         branchupdata, origbranchconss1, norigbranchconss, branchvar,
         GCG_BOUNDTYPE_LOWER, branchupdata->newbound, propagatebndchg) );
   SCIP_CALL( GCGconsMasterbranchSetOrigConsData(masterscip, cons2, downname, branchrule,
         branchdowndata, origbranchconss2, norigbranchconss, branchvar,
         GCG_BOUNDTYPE_UPPER, branchdowndata->newbound, propagatebndchg) );

   return SCIP_OKAY;
}


/** branching method for relaxation solutions */
static
SCIP_RETCODE branchExtern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< pointer to the original variable branching rule */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching call */
   )
{
   SCIP* masterscip;
   int i;

   /* parameter data */
   SCIP_Bool mostfrac;
   SCIP_Bool usepseudocosts;
   SCIP_Bool usepsstrong;

   /* branching candidates */
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandsscore;
   SCIP_Real* branchcandssol;
   int nbranchcands;
   int npriobranchcands;

   /* values for choosing the variable to branch on */
   SCIP_VAR* branchvar;
   SCIP_Real solval;
   SCIP_Real maxfrac;
   SCIP_Real frac;
   SCIP_Real maxpsscore;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPisRelaxSolValid(scip));

   *result = SCIP_DIDNOTRUN;

   /* get master problem */
   masterscip = GCGgetMasterprob(scip);
   assert(masterscip != NULL);

   /* get values of parameters */
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/mostfrac", &mostfrac) );
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/usepseudocosts", &usepseudocosts) );
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/usepsstrong", &usepsstrong) );

   /* get the branching candidates */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &branchcands, &branchcandssol, &branchcandsscore, &nbranchcands,
         &npriobranchcands, NULL, NULL, NULL) );

   branchvar = NULL;
   solval = 0.0;

   maxfrac = 0.0;
   maxpsscore = -1.0;

   if( usepsstrong )
   {
      SCIP_CALL( SCIPgetRelpsprobBranchVar(masterscip, branchcands, branchcandssol, branchcandsscore, npriobranchcands,
            npriobranchcands, result, &branchvar) );
      assert(branchvar != NULL || *result == SCIP_CUTOFF);
      assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF);

      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;

      solval = SCIPgetRelaxSolVal(scip, branchvar);
   }

   /* branch on an integer variable belonging to a unique block with fractional value */
   if( branchvar == NULL )
      for( i = 0; i < npriobranchcands; i++ )
      {
         assert(GCGvarIsOriginal(branchcands[i]));

         /* variable belongs to no block */
         if( GCGvarGetBlock(branchcands[i]) == -1 )
            continue;

         /* block is not unique (non-linking variables) */
         if( !GCGvarIsLinking(branchcands[i]) && GCGgetNIdenticalBlocks(scip, GCGvarGetBlock(branchcands[i])) != 1 )
            continue;

         /* check that blocks of linking variable are unique */
         if( GCGvarIsLinking(branchcands[i]) )
         {
            int nvarblocks;
            int* varblocks;
            SCIP_Bool unique;
            int j;

            nvarblocks = GCGlinkingVarGetNBlocks(branchcands[i]);
            SCIP_CALL( SCIPallocBufferArray(scip, &varblocks, nvarblocks) );
            SCIP_CALL( GCGlinkingVarGetBlocks(branchcands[i], nvarblocks, varblocks) );

            unique = TRUE;
            for( j = 0; j < nvarblocks; ++j )
               if( GCGgetNIdenticalBlocks(scip, varblocks[j]) != 1 )
                  unique = FALSE;

            SCIPfreeBufferArray(scip, &varblocks);

            if( !unique )
               continue;
         }

         /* use pseudocost variable selection rule */
         if( usepseudocosts )
         {
            /* select the variable, if its pseudocost are higher than the ones of the currently saved variable */
            if( SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]) > maxpsscore )
            {
               branchvar = branchcands[i];
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               maxpsscore = SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]);
            }
         }
         /* use most fractional variable selection rule */
         else
         {
            /* compute the fractionality */
            frac = MIN( branchcandsscore[i], 1.0 - branchcandsscore[i] );
            assert(frac > 0);

            /* fractionality is higher than that of the current highest fractionality */
            if( frac >= maxfrac )
            {
               SCIPdebugMessage("Var %s has fractional value in current solution: %f\n", SCIPvarGetName(branchcands[i]), branchcandssol[i]);
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               branchvar = branchcands[i];
               /* if we do not look for the most fractional variable, but for the first fractional variable,
                * we can stop here since we found a variable to branch on */
               if( !mostfrac )
                  break;
            }
         }
      }

   /* we did not find a variable to branch on so far, so we look for an integer variable that belongs to no block
    * but was directly transferred to the master problem and which has fractional value in the current solution */
   if( branchvar == NULL )
   {
      for( i = 0; i < npriobranchcands; i++ )
      {
         assert(GCGvarIsOriginal(branchcands[i]));

         /* continue if variable belongs to a block */
         if( GCGvarGetBlock(branchcands[i]) != -1 )
            continue;

         /* use pseudocost variable selection rule */
         if( usepseudocosts )
         {
            if( SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]) > maxpsscore )
            {
               branchvar = branchcands[i];
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               maxpsscore = SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]);
            }
         }
         /* use most fractional variable selection rule */
         else
         {
            /* compute fractionality */
            frac = MIN( branchcandsscore[i], 1.0 - branchcandsscore[i] );
            assert(frac > 0);

            if( frac >= maxfrac )
            {
               SCIPdebugMessage("Var %s has fractional value in current solution: %f\n",
                  SCIPvarGetName(branchcands[i]), branchcandssol[i]);
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               branchvar = branchcands[i];
               /* if we do not look for the most fractional variable, but for the first fractional variable,
                * we stop here since we found a variable to branch on */
               if( !mostfrac )
                  break;
            }
         }
      }
   }

   if( branchvar == NULL )
   {
      SCIPdebugMessage("Original branching rule could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(branchvar != NULL);

   SCIP_CALL( branchVar(scip, branchrule, branchvar, solval) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/*
 * Callback methods for enforcing branching constraints
 */

#define branchDeactiveMasterOrig NULL
#define branchPropMasterOrig NULL

/** callback activation method */
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterOrig)
{
   SCIP* origscip;
   SCIP_CONS* mastercons;

   assert(scip != NULL);
   assert(branchdata != NULL);

   /* branching restrictions are enforced by variable bounds, this is done automatically, so we can abort here */
   if( branchdata->cons == NULL )
      return SCIP_OKAY;

   assert(branchdata->origvar != NULL);

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   SCIPdebugMessage("branchActiveMasterOrig: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == GCG_BOUNDTYPE_LOWER ? ">=" : "<=" ), branchdata->newbound);

   /* transform constraint to the master variable space */
   SCIP_CALL( GCGrelaxTransOrigToMasterCons(origscip, branchdata->cons, &mastercons) );
   assert(mastercons != NULL);

   /* add constraint to the master problem */
   SCIP_CALL( SCIPaddConsNode(scip, SCIPgetCurrentNode(scip), mastercons, NULL) );

   /* constraint was added locally to the node where it is needed, so we do not need to care about this
    * at the next activation of the node and can set the constraint pointer to NULL */
   SCIP_CALL( SCIPreleaseCons(scip, &branchdata->cons) );
   branchdata->cons = NULL;

   return SCIP_OKAY;
}

/** callback solved method */
static
GCG_DECL_BRANCHMASTERSOLVED(branchMasterSolvedOrig)
{
   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->origvar != NULL);

   SCIPdebugMessage("branchMasterSolvedOrig: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == GCG_BOUNDTYPE_LOWER ? ">=" : "<=" ), branchdata->newbound);

   if( !SCIPisInfinity(scip, newlowerbound) && SCIPgetStage(GCGgetMasterprob(scip)) == SCIP_STAGE_SOLVING
      && SCIPisRelaxSolValid(GCGgetMasterprob(scip)) )
   {
      SCIP_CALL( SCIPupdateVarPseudocost(scip, branchdata->origvar,
            SCIPgetRelaxSolVal(scip, branchdata->origvar) - branchdata->oldvalue,
            newlowerbound - branchdata->olddualbound, 1.0) );
   }

   return SCIP_OKAY;
}

/** callback deletion method for branching data */
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteOrig)
{
   assert(scip != NULL);
   assert(branchdata != NULL);

   if( *branchdata == NULL )
      return SCIP_OKAY;

   SCIPdebugMessage("branchDataDeleteOrig: %s %s %f\n", SCIPvarGetName((*branchdata)->origvar),
      ( (*branchdata)->boundtype == GCG_BOUNDTYPE_LOWER ? ">=" : "<=" ), (*branchdata)->newbound);

   /* release constraint */
   if( (*branchdata)->cons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*branchdata)->cons) );
      (*branchdata)->cons = NULL;
   }

   SCIPfreeBlockMemoryNull(scip, branchdata);
   *branchdata = NULL;

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpOrig)
{  /*lint --e{715}*/
   SCIP* origscip;
   SCIP_Bool feasible;

   SCIPdebugMessage("Execlp method of orig branching\n");

   /* get original problem */
   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   feasible = FALSE;

   if( GCGcurrentNodeIsGeneric(scip) )
   {
      SCIPdebugMessage("Not executing orig branching, node was branched by generic branchrule\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Update current sol.\n");
   SCIP_CALL( GCGrelaxUpdateCurrentSol(origscip, &feasible) );

   /* if the transferred master solution is feasible, the current node is solved to optimality and can be pruned */
   if( feasible )
   {
      *result = SCIP_CUTOFF;
      SCIPdebugMessage("solution was feasible, node can be cut off!");
   }

   if( SCIPgetNExternBranchCands(origscip) > 0 )
   {
      assert(SCIPisRelaxSolValid(origscip));
      SCIP_CALL( branchExtern(origscip, branchrule, result) );
   }

   return SCIP_OKAY;
}

/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextOrig)
{  /*lint --e{715}*/
   SCIP* origscip;
   SCIP_Bool feasible;

   SCIPdebugMessage("Execext method of orig branching\n");

   /* get original problem */
   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   feasible = FALSE;

   if( GCGcurrentNodeIsGeneric(scip) )
   {
      SCIPdebugMessage("Not executing orig branching, node was branched by generic branchrule\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Update current sol.\n");
   SCIP_CALL( GCGrelaxUpdateCurrentSol(origscip, &feasible) );

   /* if the transferred master solution is feasible, the current node is solved to optimality and can be pruned */
   if( feasible )
   {
      *result = SCIP_CUTOFF;
      SCIPdebugMessage("solution was feasible, node can be cut off!");
   }
   SCIP_CALL( branchExtern(origscip, branchrule, result) );

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitOrig)
{
   SCIP* origprob;

   origprob = GCGmasterGetOrigprob(scip);
   assert(branchrule != NULL);
   assert(origprob != NULL);

   SCIP_CALL( GCGrelaxIncludeBranchrule( origprob, branchrule, branchActiveMasterOrig,
         branchDeactiveMasterOrig, branchPropMasterOrig, branchMasterSolvedOrig, branchDataDeleteOrig) );

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsOrig)
{  /*lint --e{715}*/
   int i;
   SCIP* origscip;

   /* branching candidates */
   SCIP_VAR** branchcands;
   int nbranchcands;
   int npriobranchcands;

   /* values for choosing the variable to branch on */
   SCIP_VAR* branchvar;
   SCIP_Real solval;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /* get original problem */
   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   if( GCGcurrentNodeIsGeneric(scip) )
   {
      SCIPdebugMessage("Not executing orig branching, node was branched by generic branchrule\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Execps method of orig branching\n");

   *result = SCIP_DIDNOTRUN;
   if( SCIPgetStage(scip) > SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* get the branching candidates */
   SCIP_CALL( SCIPgetPseudoBranchCands(origscip, &branchcands, &nbranchcands, &npriobranchcands) );

   branchvar = NULL;
   solval = 0.0;

   /* branch on an integer variable belonging to a unique block with fractional value */
   for( i = 0; i < npriobranchcands; i++ )
   {
      assert(GCGvarIsOriginal(branchcands[i]));

      /* variable belongs to no block or the block is not unique */
      if( GCGvarGetBlock(branchcands[i]) <= -1 || GCGgetNIdenticalBlocks(origscip, GCGvarGetBlock(branchcands[i])) != 1 )
         continue;

      branchvar = branchcands[i];
      assert(SCIPvarGetUbLocal(branchvar) - SCIPvarGetLbLocal(branchvar) > 0.8);
      solval = SCIPvarGetLbLocal(branchvar) + 0.5;
   }

   /* we did not find a variable to branch on so far, so we look for an unfixed linking variable or an integer variable
    * that belongs to no block but was directly transferred to the master problem
    */
   if( branchvar == NULL )
   {
      for( i = 0; i < npriobranchcands; i++ )
      {
         assert(GCGvarIsOriginal(branchcands[i]));

         /* continue if variable belongs to a block */
         if( GCGvarGetBlock(branchcands[i]) > -1 )
            continue;

         /* check that blocks of linking variable are unique */
         if( GCGvarIsLinking(branchcands[i]) )
         {
            int nvarblocks;
            int* varblocks;
            SCIP_Bool unique;
            int j;

            nvarblocks = GCGlinkingVarGetNBlocks(branchcands[i]);
            SCIP_CALL( SCIPallocBufferArray(origscip, &varblocks, nvarblocks) );
            SCIP_CALL( GCGlinkingVarGetBlocks(branchcands[i], nvarblocks, varblocks) );

            unique = TRUE;
            for( j = 0; j < nvarblocks; ++j )
               if( GCGgetNIdenticalBlocks(origscip, varblocks[j]) != 1 )
                  unique = FALSE;

            SCIPfreeBufferArray(origscip, &varblocks);

            if( !unique )
               continue;
         }

         branchvar = branchcands[i];
         assert(SCIPvarGetUbLocal(branchvar) - SCIPvarGetLbLocal(branchvar) > 0.8);
         solval = SCIPvarGetLbLocal(branchvar) + 0.5;
      }
   }

   if( branchvar == NULL )
   {
      SCIPdebugMessage("Original branching rule could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(branchvar != NULL);

   SCIP_CALL( branchVar(origscip, branchrule, branchvar, solval) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}



/*
 * branching specific interface methods
 */

/** creates the branching on original variable branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleOrig(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP* origscip;
   SCIP_BRANCHRULE* branchrule;

   /* get original problem */
   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
            BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, NULL) );
   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitOrig) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpOrig) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextOrig) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsOrig) );

   /* add original variable branching rule parameters */

   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/orig/enforcebycons",
         "should bounds on variables be enforced by constraints(TRUE) or by bounds(FALSE)",
         NULL, FALSE, DEFAULT_ENFORCEBYCONS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/orig/mostfrac",
         "should branching be performed on the most fractional variable instead of the first variable?",
         NULL, FALSE, DEFAULT_MOSTFRAC, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/orig/usepseudocosts",
         "should pseudocosts be used to determine the variable on which the branching is performed?",
         NULL, FALSE, DEFAULT_USEPSEUDO, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/orig/usepsstrong",
         "should strong branching with propagation be used to determine the variable on which the branching is performed?",
         NULL, FALSE, DEFAULT_USEPSSTRONG, NULL, NULL) );

   /* notify cons_integralorig about the original variable branching rule */
   SCIP_CALL( GCGconsIntegralorigAddBranchrule(scip, branchrule) );

   return SCIP_OKAY;
}
