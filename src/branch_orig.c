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

#define BRANCHRULE_NAME          "orig"
#define BRANCHRULE_DESC          "branching for the original program in generic column generation"
#define BRANCHRULE_PRIORITY      100
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0


/** branching data for branching decisions */
struct GCG_BranchData
{
   SCIP_VAR*          origvar;               /**< original variable on which the branching is done */
   GCG_CONSSENSE      conssense;             /**< sense of the branching on the original variable: 
                                              *   greater-equal (GCG_CONSSENSE_GE) or smaller-equal (GCG_CONSSENSE_LE) */
   SCIP_Real          val;                   /**< new lower/upper bound of the original variable */
   SCIP_Real          oldbound;              /**< old lower/upper bound of the pricing variable */
};


/*
 * Callback methods for enforcing branching constraints
 */

static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterOrig)
{
   SCIP* origscip;
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(branchdata != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   printf("branchActiveMasterOrig\n");

   /* get vardata*/
   vardata = SCIPvarGetData(branchdata->origvar);
   assert(vardata != NULL);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata->blocknr >= 0 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
   assert(vardata->data.origvardata.pricingvar != NULL);

   /* set corresponding bound in the pricing problem */
   /* lower bound was changed */
   if ( branchdata->conssense == GCG_CONSSENSE_GE )
   {
      branchdata->oldbound = SCIPvarGetLbOriginal(vardata->data.origvardata.pricingvar);
      SCIP_CALL( SCIPchgVarLb(GCGrelaxGetPricingprob(origscip, vardata->blocknr),
            vardata->data.origvardata.pricingvar, branchdata->val) );
   }
   /* upper bound was changed */
   else
   {
      branchdata->oldbound = SCIPvarGetUbOriginal(vardata->data.origvardata.pricingvar);
      SCIP_CALL( SCIPchgVarUb(GCGrelaxGetPricingprob(origscip, vardata->blocknr),
            vardata->data.origvardata.pricingvar, branchdata->val) );
   }
   
   return SCIP_OKAY;
}

static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterOrig)
{
   SCIP* origscip;
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(branchdata != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   /* get vardata*/
   vardata = SCIPvarGetData(branchdata->origvar);
   assert(vardata != NULL);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata->blocknr >= 0 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
   assert(vardata->data.origvardata.pricingvar != NULL);

   /* reset corresponding bound in the pricing problem */
   /* lower bound was changed */
   if ( branchdata->conssense == GCG_CONSSENSE_GE )
   {
      SCIP_CALL( SCIPchgVarLb(GCGrelaxGetPricingprob(origscip, vardata->blocknr), vardata->data.origvardata.pricingvar, branchdata->oldbound) );
   }
   /* upper bound was changed */
   else
   {
      SCIP_CALL( SCIPchgVarUb(GCGrelaxGetPricingprob(origscip, vardata->blocknr), vardata->data.origvardata.pricingvar, branchdata->oldbound) );
   }
   
   return SCIP_OKAY;
}

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

   *result = SCIP_DIDNOTFIND;

   propcount = 0;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
          
   for ( i = 0; i < nvars; i++)
   {
      if ( !SCIPisFeasZero(scip, SCIPvarGetUbLocal(vars[i])) )
      {
         vardata = SCIPvarGetData(vars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_MASTER);
         assert(vardata->blocknr >= -1 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
         assert(vardata->data.mastervardata.norigvars > 0);
         assert(vardata->data.mastervardata.origvals != NULL);
         assert(vardata->data.mastervardata.origvars != NULL);
            
         for ( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
         {
               if ( vardata->data.mastervardata.origvars[j] == branchdata->origvar )
               {
                  if ( branchdata->conssense == GCG_CONSSENSE_GE && 
                     SCIPisFeasLT(scip, vardata->data.mastervardata.origvals[j], branchdata->val) )
                  {
                     SCIPchgVarUb(scip, vars[i], 0.0);
                     propcount++;
                     break;
                  }
                  if ( branchdata->conssense == GCG_CONSSENSE_LE && 
                     SCIPisFeasGT(scip, vardata->data.mastervardata.origvals[j], branchdata->val) )
                  {
                     SCIPchgVarUb(scip, vars[i], 0.0);
                     propcount++;
                     break;
                  }
                  
               }
               
         }
      }
   }
      
   SCIPdebugMessage("Finished propagation of masterbranch constraint: %s %s %f, %d vars fixed.\n",
      SCIPvarGetName(branchdata->origvar), (branchdata->conssense == GCG_CONSSENSE_GE ? ">=" : "<="), branchdata->val, propcount);

   if ( propcount > 0 )
   {
      *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}

#define branchDataDeleteOrig NULL


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


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsOrig)
{
   SCIP_SOL* currentsol;
   SCIP_VAR** vars;
   int nvars;
   int nbinvars;
   int nintvars;
   int i;

   SCIP_NODE* childup;
   SCIP_NODE* childdown;
   SCIP_CONS* consup;
   SCIP_CONS* consdown;

   SCIP_CONS* origbranchup;
   SCIP_CONS* origbranchdown;

   GCG_BRANCHDATA* branchupdata;
   GCG_BRANCHDATA* branchdowndata;

   SCIP_VARDATA* vardata;


   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of orig branching\n");

   *result = SCIP_DIDNOTRUN;

   /* get current sol */
   currentsol = GCGrelaxGetCurrentOrigSol(scip);

   /* get the variables of the original problem and the numbers of variable types */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* search for an integer variable with fractional value */
   for ( i = 0; i < nbinvars + nintvars; i++ )
   {
      assert(SCIPvarGetType(vars[i]) == (i < nbinvars ? SCIP_VARTYPE_BINARY : SCIP_VARTYPE_INTEGER));

      vardata = SCIPvarGetData(vars[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->blocknr >= 0 && vardata->blocknr < GCGrelaxGetNPricingprobs(scip));
      
      if ( GCGrelaxGetNIdenticalBlocks(scip, vardata->blocknr) != 1 )
      {
         continue;
      }
      
      if ( !SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, currentsol, vars[i])))
      {
         SCIPdebugMessage("Var %s has fractional value in current solution: %f\n", SCIPvarGetName(vars[i]), SCIPgetSolVal(scip, currentsol, vars[i]));
         break;
      }
   }

   if ( i == nbinvars + nintvars )
   {
      printf("Original branching rule could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(i < nbinvars + nintvars);

   /* create the b&b-tree child-nodes of the current node */
   SCIP_CALL( SCIPcreateChild(scip, &childup, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &childdown, 0.0, SCIPgetLocalTransEstimate(scip)) );

   /* create corresponding constraints */
   SCIP_CALL( SCIPcreateConsLinear(scip, &consup, "branch_up", 0, NULL, NULL, 
         SCIPceil(scip, SCIPgetSolVal(scip, currentsol, vars[i])), SCIPinfinity(scip), 
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcreateConsLinear(scip, &consdown, "branch_down", 0, NULL, NULL, 
         -1.0 * SCIPinfinity(scip), SCIPfloor(scip, SCIPgetSolVal(scip, currentsol, vars[i])),  
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddCoefLinear(scip, consup, vars[i], 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, consdown, vars[i], 1.0) );

   SCIP_CALL( SCIPallocMemory(scip, &(branchupdata)) );
   SCIP_CALL( SCIPallocMemory(scip, &(branchdowndata)) );

   branchupdata->origvar = vars[i];
   branchupdata->conssense = GCG_CONSSENSE_GE;
   branchupdata->val = SCIPceil(scip, SCIPgetSolVal(scip, currentsol, vars[i]));
   branchupdata->oldbound = 0.0;

   branchdowndata->origvar = vars[i];
   branchdowndata->conssense = GCG_CONSSENSE_GE;
   branchdowndata->val = SCIPfloor(scip, SCIPgetSolVal(scip, currentsol, vars[i]));
   branchdowndata->oldbound = 0.0;

   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchup, "branchup", consup, vars[i], GCG_CONSSENSE_GE, 
         SCIPceil(scip, SCIPgetSolVal(scip, currentsol, vars[i])), childup, 
         GCGconsOrigbranchGetActiveCons(scip), branchrule, branchupdata) );
   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchdown, "branchdown", consdown, vars[i], GCG_CONSSENSE_LE, 
         SCIPfloor(scip, SCIPgetSolVal(scip, currentsol, vars[i])), childdown, 
         GCGconsOrigbranchGetActiveCons(scip), branchrule, branchdowndata) );

   /* add constraints to nodes */
   SCIP_CALL( SCIPaddConsNode(scip, childup, consup, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdown, consdown, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childup, origbranchup, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdown, origbranchdown, NULL) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &consup) );
   SCIP_CALL( SCIPreleaseCons(scip, &consdown) );
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
         branchDeactiveMasterOrig, branchPropMasterOrig, branchDataDeleteOrig) );
   
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
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create inference branching rule data */
   branchruledata = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchFreeOrig, branchInitOrig, branchExitOrig, branchInitsolOrig, branchExitsolOrig, 
         branchExeclpOrig, branchExecpsOrig,
         branchruledata) );

   return SCIP_OKAY;
}
