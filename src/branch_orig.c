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
//#define SCIP_DEBUG
/**@file   branch_orig.c
 * @ingroup BRANCHINGRULES
 * @brief  branching rule for original problem in gcg
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_orig.h"
#include "relax_gcg.h"
#include "cons_origbranch.h"
#include "type_branchgcg.h"
#include "scip/type_lp.h"

#define BRANCHRULE_NAME          "orig"
#define BRANCHRULE_DESC          "branching for the original program in generic column generation"
#define BRANCHRULE_PRIORITY      100
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define DEFAULT_ENFORCEBYCONS FALSE
#define DEFAULT_MOSTFRAC      FALSE
#define DEFAULT_USEPSEUDO     TRUE


/** branching data for branching decisions */
struct GCG_BranchData
{
   SCIP_VAR*          origvar;               /**< original variable on which the branching is done */
   SCIP_BOUNDTYPE     boundtype;             /**< type of the new bound of original variable */
   SCIP_Real          newbound;              /**< new lower/upper bound of the original variable */
   SCIP_Real          oldbound;              /**< old lower/upper bound of the pricing variable */
   SCIP_Real          oldvalue;              /**< old value of the original variable */
   SCIP_Real          olddualbound;          /**< dual bound before the branching was performed */
};


/*
 * Callback methods for enforcing branching constraints
 */

#if 0
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterOrig)
{
   SCIP* origscip;
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(branchdata != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   assert(branchdata->origvar != NULL);

   SCIPdebugMessage("branchActiveMasterOrig: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=" ), branchdata->newbound);

   /* get vardata*/
   vardata = SCIPvarGetData(branchdata->origvar);
   assert(vardata != NULL);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata->blocknr >= 0 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
   assert(vardata->data.origvardata.pricingvar != NULL);

   /* set corresponding bound in the pricing problem */
   /* lower bound was changed */
   if ( branchdata->boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      //branchdata->oldbound = SCIPvarGetLbOriginal(vardata->data.origvardata.pricingvar);
      assert(branchdata->oldbound < branchdata->newbound);
      assert(SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar) == branchdata->newbound);
      if ( SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar) == branchdata->newbound )
         SCIP_CALL( SCIPchgVarLb(GCGrelaxGetPricingprob(origscip, vardata->blocknr),
               vardata->data.origvardata.pricingvar, branchdata->newbound) );
   }
   /* upper bound was changed */
   else
   {
      //branchdata->oldbound = SCIPvarGetUbOriginal(vardata->data.origvardata.pricingvar);
      assert(branchdata->oldbound > branchdata->newbound);
      assert(SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar) == branchdata->newbound);
      if ( SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar) != branchdata->newbound )
         SCIP_CALL( SCIPchgVarUb(GCGrelaxGetPricingprob(origscip, vardata->blocknr),
               vardata->data.origvardata.pricingvar, branchdata->newbound) );
   }

   return SCIP_OKAY;
}
#else
#define branchActiveMasterOrig NULL
#endif

#if 0
static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterOrig)
{
   SCIP* origscip;
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(branchdata != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   assert(branchdata->origvar != NULL);

   SCIPdebugMessage("branchDeactiveMasterOrig: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=" ), branchdata->newbound);

   /* get vardata*/
   vardata = SCIPvarGetData(branchdata->origvar);
   assert(vardata != NULL);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata->blocknr >= 0 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
   assert(vardata->data.origvardata.pricingvar != NULL);

   /* reset corresponding bound in the pricing problem */
   /* lower bound was changed */
   if ( branchdata->boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      assert(SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar) == branchdata->oldbound);
      //SCIP_CALL( SCIPchgVarLb(GCGrelaxGetPricingprob(origscip, vardata->blocknr), vardata->data.origvardata.pricingvar, branchdata->oldbound) );
   }
   /* upper bound was changed */
   else
   {
      assert(SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar) == branchdata->oldbound);
      //SCIP_CALL( SCIPchgVarUb(GCGrelaxGetPricingprob(origscip, vardata->blocknr), vardata->data.origvardata.pricingvar, branchdata->oldbound) );
   }

   return SCIP_OKAY;
}
#else
#define branchDeactiveMasterOrig NULL
#endif

#if 0
static
GCG_DECL_BRANCHPROPMASTER(branchPropMasterOrig)
{
   SCIP* origscip;
   SCIP_VARDATA* vardata;
   SCIP_VAR** vars;
   int nvars;
   int propcount;
   int i;
   int j;

   assert(scip != NULL);
   assert(branchdata != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   assert(branchdata->origvar != NULL);

   SCIPdebugMessage("branchPropMasterOrig: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=" ), branchdata->newbound);

   *result = SCIP_DIDNOTFIND;

   propcount = 0;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
          
   /* iterate over all master variables */
   for ( i = 0; i < nvars; i++)
   {
      /* only look at variables not fixed to 0 */
      if ( !SCIPisZero(scip, SCIPvarGetUbLocal(vars[i])) )
      {
         vardata = SCIPvarGetData(vars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_MASTER);
         assert(vardata->blocknr >= -1 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
         assert(vardata->data.mastervardata.norigvars >= 0);
         assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
         assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);


         /* iterate over all original variables contained in the current master variable */
         for ( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
         {
            /* check whether the original variable contained in the master variable equals the variable 
             * on which the current branching was performed and if so, fix the master variable to 0,
             * if the master variable contains a part of the branching variable that violates the bound */
            if ( vardata->data.mastervardata.origvars[j] == branchdata->origvar )
            {
               /* branching imposes new lower bound */
               if ( branchdata->boundtype == SCIP_BOUNDTYPE_LOWER && 
                  SCIPisLT(scip, vardata->data.mastervardata.origvals[j], branchdata->newbound) )
               {
                  SCIPchgVarUb(scip, vars[i], 0.0);
                  propcount++;
                  break;
               }
               /* branching imposes new upper bound */
               if ( branchdata->boundtype == SCIP_BOUNDTYPE_UPPER && 
                  SCIPisGT(scip, vardata->data.mastervardata.origvals[j], branchdata->newbound) )
               {
                  SCIPchgVarUb(scip, vars[i], 0.0);
                  propcount++;
                  break;
               }
               
            }
               
         }
      }
   }
      
   SCIPdebugMessage("Finished propagation of branching decision constraint: %s %s %f, %d vars fixed.\n",
      SCIPvarGetName(branchdata->origvar), (branchdata->boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<="), 
      branchdata->newbound, propcount);

   if ( propcount > 0 )
   {
      *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}
#else
#define branchPropMasterOrig NULL
#endif

static
GCG_DECL_BRANCHMASTERSOLVED(branchMasterSolved)
{
   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->origvar != NULL);

   SCIPdebugMessage("branchMasterSolved: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=" ), branchdata->newbound);

   if( !SCIPisInfinity(scip, newlowerbound) )
   {
      SCIP_CALL( SCIPupdateVarPseudocost(scip, branchdata->origvar, 
            SCIPgetRelaxSolVal(scip, branchdata->origvar) - branchdata->oldvalue, 
            newlowerbound - branchdata->olddualbound, 1.0) );
   }

   //printf("Old lowerbound = %g, current lowerbound = %g, pseudocost of var %s = %g / %g\n",
   //   branchdata->olddualbound, newlowerbound, SCIPvarGetName(branchdata->origvar),
   //   SCIPgetVarPseudocost(scip, branchdata->origvar, 1.0), SCIPgetVarPseudocost(scip, branchdata->origvar, -1.0));

      
   SCIPdebugMessage("Finished branchMasterSolved: %s %s %f.\n",
      SCIPvarGetName(branchdata->origvar), (branchdata->boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<="), 
      branchdata->newbound);

   return SCIP_OKAY;
}

static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteOrig)
{
   assert(scip != NULL);
   assert(branchdata != NULL);
   
   SCIPdebugMessage("branchDataDeleteOrig: %s %s %f\n", SCIPvarGetName((*branchdata)->origvar),
      ( (*branchdata)->boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=" ), (*branchdata)->newbound);

   SCIPfreeMemory(scip, branchdata);

   *branchdata = NULL;

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpOrig)
{  
   SCIPdebugMessage("Execlp method of orig branching\n");
   printf("Execlp method of orig branching\n");

   return SCIP_OKAY;
}

/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECREL(branchExecrelOrig)
{
   SCIP_SOL* currentsol;
   SCIP_VAR** vars;
   int nvars;
   int nbinvars;
   int nintvars;
   int i;

   SCIP_VAR* branchvar;
   SCIP_Real solval;

   SCIP_Bool enforcebycons;
   SCIP_Bool mostfrac;
   SCIP_Bool usepseudocosts;


   SCIP_NODE* childup;
   SCIP_NODE* childdown;

   SCIP_CONS* origbranchup;
   SCIP_CONS* origbranchdown;

   GCG_BRANCHDATA* branchupdata;
   GCG_BRANCHDATA* branchdowndata;

   SCIP_VARDATA* vardata;

   char upname[SCIP_MAXSTRLEN];
   char downname[SCIP_MAXSTRLEN];

   SCIP_VAR** branchcands;
   SCIP_Real* branchcandsscore;
   SCIP_Real* branchcandssol;
   int nbranchcands;
   int npriobranchcands;

   SCIP_Real maxfrac;
   SCIP_Real frac;
   SCIP_Real maxpsscore;


   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPisRelaxSolValid(scip));

   SCIPdebugMessage("Execps method of orig branching\n");

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/enforcebycons", &enforcebycons) );
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/enforcebycons", &mostfrac) );
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/usepseudocosts", &usepseudocosts) );

   /* get current sol */
   currentsol = GCGrelaxGetCurrentOrigSol(scip);

   /* get the variables of the original problem and the numbers of variable types */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   SCIP_CALL( SCIPgetRelaxBranchCands(scip, &branchcands, &branchcandssol, &branchcandsscore, &nbranchcands,
         &npriobranchcands, NULL, NULL, NULL) );

   branchvar = NULL;
   solval = 0.0;

   maxfrac = 0.0;
   maxpsscore = -1.0;

   for ( i = 0; i < npriobranchcands; i++ )
   {
      vardata = SCIPvarGetData(branchcands[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->blocknr >= -1 && vardata->blocknr < GCGrelaxGetNPricingprobs(scip));
      
      if ( vardata->blocknr == -1 || GCGrelaxGetNIdenticalBlocks(scip, vardata->blocknr) != 1 )
         continue;

      frac = MIN( branchcandsscore[i], 1.0 - branchcandsscore[i] );
      assert(frac > 0);

      if( usepseudocosts )
      {
         if( SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]) > maxpsscore )
         {
            branchvar = branchcands[i];
            solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
            maxpsscore = SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]);
         }
      }
      else
      {
         if( frac >= maxfrac  )
         {
            SCIPdebugMessage("Var %s has fractional value in current solution: %f\n", SCIPvarGetName(vars[i]), branchcandssol[i]);
            solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
            branchvar = branchcands[i];
            if ( !mostfrac )
               break;
         }
      }
   }

   if (branchvar == NULL)
   {
      for ( i = 0; i < npriobranchcands; i++ )
      {
         vardata = SCIPvarGetData(branchcands[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         assert(vardata->blocknr >= -1 && vardata->blocknr < GCGrelaxGetNPricingprobs(scip));
         
         if ( vardata->blocknr != -1 )
            continue;
         
         frac = MIN( branchcandsscore[i], 1.0 - branchcandsscore[i] );
         assert(frac > 0);
         
         if( usepseudocosts )
         {
            if( SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]) > maxpsscore )
            {
               branchvar = branchcands[i];
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               maxpsscore = SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]);
            }
         }
         else
         {
            if( frac >= maxfrac  )
            {
               SCIPdebugMessage("Var %s has fractional value in current solution: %f\n", SCIPvarGetName(vars[i]), branchcandssol[i]);
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               branchvar = branchcands[i];
               if ( !mostfrac )
                  break;
            }
         }
      }
   }

#if 0
   if ( branchvar == NULL )
   {
      SCIPdebugMessage("All vars have integral value in current solution, branch on a non-fixed var\n");


      /* search for an integer variable which is not fixed */
      for ( i = 0; i < nbinvars + nintvars; i++ )
      {
         assert(SCIPvarGetType(vars[i]) == (i < nbinvars ? SCIP_VARTYPE_BINARY : SCIP_VARTYPE_INTEGER));

         if ( SCIPvarGetUbLocal(vars[i]) - SCIPvarGetLbLocal(vars[i]) < 0.5 )
         {
            continue;
         }

         vardata = SCIPvarGetData(vars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         assert(vardata->blocknr >= 0 && vardata->blocknr < GCGrelaxGetNPricingprobs(scip));
      
         if ( GCGrelaxGetNIdenticalBlocks(scip, vardata->blocknr) != 1 )
            continue;

         //solval = SCIPgetSolVal(scip, currentsol, vars[i]);
         solval = SCIPgetRelaxSolVal(scip, vars[i]);
         assert(SCIPisIntegral(scip, solval));
         
         branchvar = vars[i];

         if ( SCIPvarGetUbLocal(vars[i]) - solval < 0.5 )
         {
            solval = solval - 0.5;
         }
         else
         {
            solval = solval + 0.5;
         }
         break;
      }
   }
#endif
   if ( branchvar == NULL )
   {
      SCIPdebugMessage("Original branching rule could not find a variable to branch on!\n");
      //printf("Original branching rule could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(branchvar != NULL);

   SCIPdebugMessage("Branching on var %s with value %g in current solution\n", SCIPvarGetName(branchvar), solval);

   /* create the b&b-tree child-nodes of the current node */
   SCIP_CALL( SCIPcreateChild(scip, &childup, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &childdown, 0.0, SCIPgetLocalTransEstimate(scip)) );

   SCIP_CALL( SCIPallocMemory(scip, &(branchupdata)) );
   SCIP_CALL( SCIPallocMemory(scip, &(branchdowndata)) );

   if ( enforcebycons )
   {
      /* enforce new bounds by linear constraints */
      SCIP_CONS* consup;
      SCIP_CONS* consdown;

      /* create corresponding constraints */
      SCIP_CALL( SCIPcreateConsLinear(scip, &consup, "branch_up", 0, NULL, NULL, 
            SCIPceil(scip, solval), SCIPinfinity(scip), 
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &consdown, "branch_down", 0, NULL, NULL, 
            -1.0 * SCIPinfinity(scip), SCIPfloor(scip, solval),  
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPaddCoefLinear(scip, consup, branchvar, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, consdown, branchvar, 1.0) );

      /* add constraints to nodes */
      SCIP_CALL( SCIPaddConsNode(scip, childup, consup, NULL) );
      SCIP_CALL( SCIPaddConsNode(scip, childdown, consdown, NULL) );

      /* release constraints */
      SCIP_CALL( SCIPreleaseCons(scip, &consup) );
      SCIP_CALL( SCIPreleaseCons(scip, &consdown) );
   }
   else
   {
      /* enforce new bounds by setting variable bounds */
      SCIP_CALL( SCIPchgVarUbNode(scip, childdown, branchvar, solval) );
      SCIP_CALL( SCIPchgVarLbNode(scip, childup, branchvar, solval) );
   }

   branchupdata->origvar = branchvar;
   branchupdata->oldvalue = solval;
   branchupdata->olddualbound = SCIPgetLocalLowerbound(GCGrelaxGetMasterprob(scip));
   branchupdata->boundtype = SCIP_BOUNDTYPE_LOWER;
   branchupdata->newbound = SCIPceil(scip, solval);
   branchupdata->oldbound = SCIPvarGetLbLocal(branchvar);

   branchdowndata->origvar = branchvar;
   branchdowndata->oldvalue = solval;
   branchdowndata->olddualbound = SCIPgetLocalLowerbound(GCGrelaxGetMasterprob(scip));
   branchdowndata->boundtype = SCIP_BOUNDTYPE_UPPER;
   branchdowndata->newbound = SCIPfloor(scip, solval);
   branchdowndata->oldbound = SCIPvarGetUbLocal(branchvar);

   (void) SCIPsnprintf(upname, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchupdata->origvar), 
      ">=", branchupdata->newbound);
   (void) SCIPsnprintf(downname, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchdowndata->origvar), 
      "<=", branchdowndata->newbound);

   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchup, upname, childup, 
         GCGconsOrigbranchGetActiveCons(scip), branchrule, branchupdata) );
   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchdown, downname, childdown, 
         GCGconsOrigbranchGetActiveCons(scip), branchrule, branchdowndata) );

   /* add constraints to nodes */
   SCIP_CALL( SCIPaddConsNode(scip, childup, origbranchup, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdown, origbranchdown, NULL) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &origbranchup) );
   SCIP_CALL( SCIPreleaseCons(scip, &origbranchdown) );
      
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitOrig)
{  
   assert(branchrule != NULL);

   SCIP_CALL( GCGrelaxIncludeBranchrule(scip, branchrule, branchActiveMasterOrig, 
         branchDeactiveMasterOrig, branchPropMasterOrig, branchMasterSolved, branchDataDeleteOrig) );
   
   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsOrig)
{  
   SCIPdebugMessage("Execps method of orig branching\n");
   assert(0);

   return SCIP_OKAY;
}

/* define not used callback as NULL*/
#define branchFreeOrig NULL
#define branchExitOrig NULL
#define branchInitsolOrig NULL
#define branchExitsolOrig NULL


/*
 * branching specific interface methods
 */

/** creates the most infeasible LP braching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleOrig(
   SCIP*                 scip                /**< SCIP data structure */
   )
{   
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchFreeOrig, branchInitOrig, branchExitOrig, branchInitsolOrig, branchExitsolOrig, 
         branchExeclpOrig, branchExecrelOrig, branchExecpsOrig,
         NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/orig/enforcebycons",
         "should bounds on variables be enforced by constraints(TRUE) or by bounds(FALSE)",
         NULL, FALSE, DEFAULT_ENFORCEBYCONS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/orig/mostfrac",
         "should branching be performed on the most fractional variable instead of the first variable?",
         NULL, FALSE, DEFAULT_MOSTFRAC, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/orig/usepseudocosts",
         "should pseudocosts be used to determine the variable on which the branching is performed?",
         NULL, FALSE, DEFAULT_USEPSEUDO, NULL, NULL) );


   return SCIP_OKAY;
}
