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
/**@file   pricer_gcg.c
 * @ingroup PRICERS
 * @brief  pricer for generic column generation, solves the pricing problem as a MIP
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "pricer_gcg.h"
#include "scip/cons_linear.h"
#include "scip/scip.h"
#include "sepa_master.h"
#include <stdio.h>
#include <stdlib.h>


#define PRICER_NAME            "gcg"
#define PRICER_DESC            "pricer for gcg"
#define PRICER_PRIORITY        5000000
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

#define EPS 0.0001

#define DEFAULT_MAXVARSROUNDFARKAS 1
#define DEFAULT_MAXVARSROUNDREDCOST INT_MAX


/*
 * Data structures
 */


/** variable pricer data */
struct SCIP_PricerData
{
   int pricingprobnr;              /* number of the pricingproblem, that was used last time to find a new var */
   int npricingprobs;              /* number of pricing problems */
   SCIP** pricingprobs;            /* pointers to the pricing problems */
   SCIP_Real* redcost;             /* array of reduced cost for the constraints of the master problem */
   SCIP_Real* redcostconv;         /* array of reduced cost for the convexity constraints */
   SCIP* origprob;                 /* the original program */
   SCIP_Real* solvals;             /* solution values of variables in the pricing problems */
   int* nvarsprob;                 /* number of variables created by the pricing probs */
   SCIP_CLOCK* subsolveclock;
   SCIP_CLOCK* redcostclock;
   SCIP_CLOCK* redcostsolveclock;
   SCIP_CLOCK* redcostpresolveclock;
   SCIP_CLOCK* farkasclock;
   SCIP_CLOCK* farkassolveclock;
   SCIP_CLOCK* farkaspresolveclock;
   SCIP_CLOCK* owneffortclock;
   SCIP_CLOCK* freeclock;
   SCIP_CLOCK* transformclock;
   int solvedsubmips;
   int calls;
   int farkascalls;
   int redcostcalls;
   int initcalls;
   
   /* vartype of created master variables */
   SCIP_VARTYPE vartype;
   int maxvarsroundfarkas;
   int maxvarsroundredcost;

   SCIP_HASHMAP* mapcons2idx;

};



/*
 * Vardata methods
 */

static
SCIP_DECL_VARDELTRANS(gcgvardeltrans)
{
   assert((*vardata)->vartype == GCG_VARTYPE_MASTER);
   SCIPfreeBlockMemoryArray(scip, &((*vardata)->data.mastervardata.origvals), (*vardata)->data.mastervardata.norigvars);
   SCIPfreeBlockMemoryArray(scip, &((*vardata)->data.mastervardata.origvars), (*vardata)->data.mastervardata.norigvars);
   
   SCIPfreeBlockMemory(scip, vardata);

   return SCIP_OKAY;
}


/*
 * Local methods
 */

static
SCIP_RETCODE checkNewVar(
   SCIP*                 scip,
   SCIP_VAR*             newvar,
   SCIP_Real             redcost,
   SCIP_Real             redcostconv
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;
   int i;

   SCIP_VARDATA* newvardata;
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(newvar != NULL);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* compare newvar with all existing variables */
   for ( v = 0; v < nvars; v++ )
   {
      assert(vars[v] != NULL);

      /* vars[v] is the new variable itself */
      if ( vars[v] == newvar )
         continue;

      /* vars[v] has a different objective function value, may not be equal to newvar */
      if ( !SCIPisEQ(scip, SCIPvarGetObj(vars[v]), SCIPvarGetObj(newvar)) )
         continue;

      vardata = SCIPvarGetData(vars[v]);
      newvardata = SCIPvarGetData(newvar);
      assert(vardata != NULL);
      assert(newvardata != NULL);

      /* vars[v] belongs to a different block, may not be equal to newvar */
      if ( vardata->blocknr != newvardata->blocknr )
         continue;

      assert(vardata->vartype == GCG_VARTYPE_MASTER);
      assert(newvardata->vartype == GCG_VARTYPE_MASTER);
      assert(vardata->data.mastervardata.norigvars = newvardata->data.mastervardata.norigvars);

      /* compare the parts of the original variables contained in vars[i] and newvar */
      for ( i = 0; i < vardata->data.mastervardata.norigvars; i++ )
      {
         assert(vardata->data.mastervardata.origvars[i] == newvardata->data.mastervardata.origvars[i]);

         if ( !SCIPisEQ(scip, vardata->data.mastervardata.origvals[i], newvardata->data.mastervardata.origvals[i]) )
            break;
      }

      if ( i == vardata->data.mastervardata.norigvars )
      {
         printf("var %s is equal to var %s! solval = %f, redcost = %f, lpsolstat = %d, redcostcons = %f\n", 
            SCIPvarGetName(newvar), SCIPvarGetName(vars[v]), SCIPgetSolVal(scip, NULL, vars[v]), 
            redcost, SCIPgetLPSolstat(scip), redcostconv);
      }
   }

   return SCIP_OKAY;
}

#ifndef NDEBUG
static
SCIP_RETCODE checkVarBounds(
   SCIP*                 scip
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;
   SCIP_VARDATA* vardata;
   SCIP* origscip;

   assert(scip != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   SCIP_CALL( SCIPgetVarsData(origscip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* check whether the corresponding pricing MIP has the same bound for the variable */
   for ( v = 0; v < nvars; v++ )
   {
      assert(vars[v] != NULL);

      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->data.origvardata.pricingvar != NULL);

      if ( !GCGrelaxIsPricingprobRelevant(origscip, vardata->blocknr) || GCGrelaxGetNIdenticalBlocks(origscip, vardata->blocknr) != 1 ) 
         continue;

      if ( SCIPvarGetUbLocal(vars[v]) != SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar) )
      {
         printf("var %s: orig upper bound = %f, pricing upper bound = %f!\n", SCIPvarGetName(vars[v]),
            SCIPvarGetUbLocal(vars[v]), SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar) );
      }
      if ( SCIPvarGetLbLocal(vars[v]) != SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar) )
      {
         printf("var %s: orig lower bound = %f, pricing lower bound = %f!\n", SCIPvarGetName(vars[v]),
            SCIPvarGetLbLocal(vars[v]), SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar) );
      }

      assert(SCIPvarGetUbLocal(vars[v]) == SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar));
      assert(SCIPvarGetLbLocal(vars[v]) == SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar));
   }

   return SCIP_OKAY;
}
#endif

/* informs an original variable, that a variable in the master problem was created, 
 * that contains a part of the original variable.
 * Saves this information in the original variable's data */
SCIP_RETCODE GCGpricerAddMasterVarToOrigVar(
   SCIP*                 scip,
   SCIP_VAR*             origvar,
   SCIP_VAR*             var,
   SCIP_Real             val
   )
{
   SCIP_VARDATA* vardata;
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(origvar != NULL);
   assert(var != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);   

   vardata = SCIPvarGetData(origvar);
   assert(vardata != NULL);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata->data.origvardata.mastervars != NULL);
   assert(vardata->data.origvardata.mastervals != NULL);
   assert(vardata->data.origvardata.nmastervars >= 0);
   assert(vardata->data.origvardata.maxmastervars >= vardata->data.origvardata.nmastervars);

   /* realloc mastervars array of the original variable, if needed */
   if ( vardata->data.origvardata.maxmastervars == vardata->data.origvardata.nmastervars )
   {
      SCIP_CALL( SCIPreallocMemoryArray(pricerdata->origprob, &(vardata->data.origvardata.mastervars),
            2*vardata->data.origvardata.maxmastervars) );
      SCIP_CALL( SCIPreallocMemoryArray(pricerdata->origprob, &(vardata->data.origvardata.mastervals),
            2*vardata->data.origvardata.maxmastervars) );
      SCIPdebugMessage("mastervars array of var %s resized from %d to %d\n", SCIPvarGetName(origvar), 
         vardata->data.origvardata.maxmastervars, 2*vardata->data.origvardata.maxmastervars);
      vardata->data.origvardata.maxmastervars = 2*vardata->data.origvardata.maxmastervars;
   }
   /* add information to the original variable's vardata */
   vardata->data.origvardata.mastervars[vardata->data.origvardata.nmastervars] = var;
   vardata->data.origvardata.mastervals[vardata->data.origvardata.nmastervars] = val;
   vardata->data.origvardata.nmastervars++;

   return SCIP_OKAY;
}


/* performs the pricing routine, gets the type of pricing that should be done: farkas or redcost pricing */
static
SCIP_RETCODE performPricing(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_PRICER*          pricer,             /* the pricer */
   GCG_PRICETYPE         pricetype           /* type of the pricing */
   )
{
   SCIP_PRICERDATA* pricerdata;            /* the data of the pricer */
   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_CONS** origconss;
   char varname[SCIP_MAXSTRLEN];
   SCIP* origprob;
   SCIP_Real* mastercoefs;

   int i;
   int j;
   int k;
   int l;
   int prob;

   int nfoundvars;

   SCIP_VAR** consvars;
   int nconsvars;
   SCIP_Real* consvals;
   SCIP_VAR** probvars;
   int nprobvars;

   SCIP_VARDATA* vardata;
   SCIP_VARDATA* newvardata;

   int nsols;
   SCIP_SOL** sols;

   SCIP_VAR* newvar;

   SCIP_Real objcoeff;
   SCIP_Real conscoeff;

   SCIP_ROW** mastercuts;
   int nmastercuts;
   SCIP_ROW** origcuts;
   int norigcuts;
   SCIP_Real redcost;
   SCIP_COL** cols;
   SCIP_VAR* var;

   SCIP_Real* tmpconvredcost;
   int* permu;

   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   if ( pricetype == GCG_PRICETYPE_REDCOST )
      pricerdata->redcostcalls++;
   if ( pricetype == GCG_PRICETYPE_FARKAS )
      pricerdata->farkascalls++;

   SCIP_CALL( SCIPstartClock(scip, pricerdata->owneffortclock) );
   pricerdata->calls++;
   nfoundvars = 0;

#ifndef NDEBUG
   SCIP_CALL( checkVarBounds(scip) );
#endif

   if ( pricetype == GCG_PRICETYPE_REDCOST || SCIPgetNVars(scip) % 50 == 0 )
      printf("nvars = %d, current lowerbound = %g\n", SCIPgetNVars(scip), SCIPgetLPObjval(scip));

   /* get the constraints of the master problem and the corresponding constraints in the original problem */
   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);
   origconss = GCGrelaxGetLinearOrigMasterConss(origprob);

   SCIP_CALL( SCIPallocBufferArray(scip, &mastercoefs, nmasterconss) );

   /* set objective value of all variables in the pricing problems to 0 (for farkas pricing) /
    * to the original objective of the variable (for redcost pricing) */
   for ( i = 0; i < pricerdata->npricingprobs; i++)
   {
      if ( pricerdata->pricingprobs[i] == NULL )
         continue;
      probvars = SCIPgetVars(pricerdata->pricingprobs[i]);
      nprobvars = SCIPgetNVars(pricerdata->pricingprobs[i]);
      
      for ( j = 0; j < nprobvars; j++ )
      {
         if ( pricetype == GCG_PRICETYPE_FARKAS || pricetype == GCG_PRICETYPE_INIT )
         {
            SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], 0) );
         }
         else 
         {
            vardata = SCIPvarGetData(probvars[j]);
            assert(vardata->vartype == GCG_VARTYPE_PRICING);
            assert(vardata->blocknr == i);
            assert(vardata->data.pricingvardata.origvars != NULL);
            assert(vardata->data.pricingvardata.origvars[0] != NULL);
            SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], 
                  SCIPvarGetObj(vardata->data.pricingvardata.origvars[0])) );
         }  
      }
   }

   /* compute reduced cost and update objectives in the pricing problems */
   for ( i = 0; i < nmasterconss; i++ )
   {
      /* farkas pricing */
      if ( pricetype == GCG_PRICETYPE_FARKAS )
      {
         assert(SCIPconsIsTransformed(masterconss[i]));
         pricerdata->redcost[i] = SCIPgetDualfarkasLinear(scip, masterconss[i]);
      }
      /* redcost pricing */
      if ( pricetype == GCG_PRICETYPE_REDCOST )
      {      
         pricerdata->redcost[i] = SCIPgetDualsolLinear(scip, masterconss[i]);
      }
      if ( !SCIPisFeasZero(scip, pricerdata->redcost[i]) )
      {
         /* for all variables in the constraint, modify the objective of the corresponding variable in a pricing problem */
         consvars = SCIPgetVarsLinear(origprob, origconss[i]);
         nconsvars = SCIPgetNVarsLinear(origprob, origconss[i]);
         consvals = SCIPgetValsLinear(origprob, origconss[i]);
         for ( j = 0; j < nconsvars; j++ )
         {
            vardata = SCIPvarGetData(consvars[j]);
            assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
            if ( vardata->blocknr != -1 && pricerdata->pricingprobs[vardata->blocknr] != NULL )
            {
               assert(vardata->data.origvardata.pricingvar != NULL);
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[vardata->blocknr], 
                     vardata->data.origvardata.pricingvar, -1.0 * pricerdata->redcost[i] * consvals[j]) );
            }
         }
      }
   }

   /* get the cuts of the master problem and the corresponding cuts in the original problem */
   mastercuts = GCGsepaGetMastercuts(scip);
   nmastercuts = GCGsepaGetNMastercuts(scip);
   origcuts = GCGsepaGetOrigcuts(scip);
   norigcuts = GCGsepaGetNOrigcuts(scip);
   
   assert(mastercuts != NULL);
   assert(origcuts != NULL);
   assert(norigcuts == nmastercuts);

   /* compute reduced cost and update objectives in the pricing problems */
   for ( i = 0; i < nmastercuts; i++ )
   {
      /* farkas pricing */
      if ( pricetype == GCG_PRICETYPE_FARKAS )
      {
         redcost = SCIProwGetDualfarkas(mastercuts[i]);
      }
      /* redcost pricing */
      else
      {    
         assert(pricetype == GCG_PRICETYPE_REDCOST);
         redcost = SCIProwGetDualsol(mastercuts[i]);
      }
      if ( !SCIPisFeasZero(scip, redcost) )
      {
         /* get columns and vals of the cut */
         nconsvars = SCIProwGetNNonz(origcuts[i]);
         cols = SCIProwGetCols(origcuts[i]);
         consvals = SCIProwGetVals(origcuts[i]);

         /* get the variables corresponding to the columns in the cut */
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
         for ( j = 0; j < nconsvars; j++ )
         {
            consvars[j] = SCIPcolGetVar(cols[j]);
         }

         /* for all variables in the cut, modify the objective of the corresponding variable in a pricing problem */
         for ( j = 0; j < nconsvars; j++ )
         {
            vardata = SCIPvarGetData(consvars[j]);
            assert(vardata != NULL);
            assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
            if ( vardata->blocknr != -1 && pricerdata->pricingprobs[vardata->blocknr] != NULL )
            {
               assert(vardata->data.origvardata.pricingvar != NULL);
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[vardata->blocknr], 
                     vardata->data.origvardata.pricingvar, -1.0 * redcost * consvals[j]) );
            }
         }
         SCIPfreeBufferArray(scip, &consvars);
      }
   }

   /* get dual solutions / farkas values of the convexity constraints */
   for ( i = 0; i < pricerdata->npricingprobs; i++ )
   {

      assert( GCGrelaxIsPricingprobRelevant(origprob, i) 
         == (GCGrelaxGetConvCons(origprob, i) != NULL) );
      if ( !GCGrelaxIsPricingprobRelevant(origprob, i) )
      {
         pricerdata->redcostconv[i] = -1.0 * SCIPinfinity(scip);
         continue;
      }  
      if ( pricetype == GCG_PRICETYPE_FARKAS || pricetype == GCG_PRICETYPE_INIT )
      {
         assert(SCIPconsIsTransformed(GCGrelaxGetConvCons(origprob, i)));
         pricerdata->redcostconv[i] = SCIPgetDualfarkasLinear(scip, GCGrelaxGetConvCons(origprob, i));
      }
      else
      {
         pricerdata->redcostconv[i] = SCIPgetDualsolLinear(scip, GCGrelaxGetConvCons(origprob, i));
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpconvredcost, pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &permu, pricerdata->npricingprobs) );
   for ( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      assert(GCGrelaxIsPricingprobRelevant(origprob, i) 
         == (pricerdata->redcostconv[i] != -1.0 * SCIPinfinity(scip)));

      permu[i] = i;
      tmpconvredcost[i] = pricerdata->redcostconv[i];
   }
   SCIPsortDownRealInt(tmpconvredcost, permu, pricerdata->npricingprobs);
   

   /* solve the pricing MIPs and check whether solutions corresponding to variables with negative reduced costs where found */
   for ( i = 0; i < pricerdata->npricingprobs && 
            (pricetype == GCG_PRICETYPE_REDCOST || nfoundvars < pricerdata->maxvarsroundfarkas)
            && (pricetype == GCG_PRICETYPE_FARKAS || nfoundvars < pricerdata->maxvarsroundredcost); i++)
   {
      //prob = (pricerdata->pricingprobnr + i) % pricerdata->npricingprobs;
      prob = permu[i];

      if ( pricerdata->pricingprobs[prob] == NULL )
         continue;

      /* set objective limit, such that only solutions with negative reduced costs are accepted */
      //SCIP_CALL( SCIPsetObjlimit(pricerdata->pricingprobs[prob], pricerdata->redcostconv[prob]) );

      SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );
      SCIP_CALL( SCIPstartClock(scip, pricerdata->subsolveclock) );

      /* start clock measuring the transformation effort */
      SCIP_CALL( SCIPstartClock(scip, pricerdata->transformclock) );
      SCIP_CALL( SCIPtransformProb(pricerdata->pricingprobs[prob]) );
      SCIP_CALL( SCIPstopClock(scip, pricerdata->transformclock) );

      /* start clock measuring the presolving effort */
      if ( pricetype == GCG_PRICETYPE_REDCOST )
         SCIP_CALL( SCIPstartClock(scip, pricerdata->redcostpresolveclock) );
      else
         SCIP_CALL( SCIPstartClock(scip, pricerdata->farkaspresolveclock) );

      /* presolve the pricing submip */
      SCIP_CALL( SCIPpresolve(pricerdata->pricingprobs[prob]) );

      /* stop clock measuring the presolving effort, start clock measuring the solving effort */
      if ( pricetype == GCG_PRICETYPE_REDCOST )
      {
         SCIP_CALL( SCIPstopClock(scip, pricerdata->redcostpresolveclock) );
         SCIP_CALL( SCIPstartClock(scip, pricerdata->redcostsolveclock) );
      }
      else
      {
         SCIP_CALL( SCIPstopClock(scip, pricerdata->farkaspresolveclock) );
         SCIP_CALL( SCIPstartClock(scip, pricerdata->farkassolveclock) );
      }

      /* solve the pricing submip */
      SCIP_CALL( SCIPsolve(pricerdata->pricingprobs[prob]) );

      /* stop clock measuring the solving effort */
      if ( pricetype == GCG_PRICETYPE_REDCOST )
         SCIP_CALL( SCIPstopClock(scip, pricerdata->redcostsolveclock) );
      else
         SCIP_CALL( SCIPstopClock(scip, pricerdata->farkassolveclock) );

      pricerdata->solvedsubmips++;

      SCIP_CALL( SCIPstopClock(scip, pricerdata->subsolveclock) );
      SCIP_CALL( SCIPstartClock(scip, pricerdata->owneffortclock) );

      /* so far, the pricing problem should be solved to optimality */
      assert( SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_OPTIMAL
         || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_GAPLIMIT
         || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_USERINTERRUPT 
         || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_INFEASIBLE );
      /** @todo handle userinterrupt: set result pointer (and lowerbound), handle other solution states */
   
      nsols = SCIPgetNSols(pricerdata->pricingprobs[prob]);
      sols = SCIPgetSols(pricerdata->pricingprobs[prob]);
      //printf("Pricingprob %d has found %d sols!\n", prob, nsols);

      for ( j = 0; j < nsols && 
               (pricetype == GCG_PRICETYPE_REDCOST || nfoundvars < pricerdata->maxvarsroundfarkas)
               && (pricetype == GCG_PRICETYPE_FARKAS || nfoundvars < pricerdata->maxvarsroundredcost); j++ )
      {
#ifndef NDEBUG
         SCIP_Bool feasible;
         SCIP_CALL( SCIPcheckSolOrig(pricerdata->pricingprobs[prob], sols[j], &feasible, TRUE, TRUE) );
         assert(feasible);
#endif
         /* solution value - dual value of associated convexity constraint < 0 
            --> can make the LP feasible / improve the current solution */ 
         if ( SCIPisFeasNegative(scip, SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sols[j]) - pricerdata->redcostconv[prob]) )
         {
            nfoundvars++;
            
            /* get variables of the pricing problem and their values in the current solution */
            probvars = SCIPgetOrigVars(pricerdata->pricingprobs[prob]);
            nprobvars = SCIPgetNOrigVars(pricerdata->pricingprobs[prob]);
            SCIP_CALL( SCIPgetSolVals(pricerdata->pricingprobs[prob], sols[j], nprobvars, probvars, pricerdata->solvals) );
            
            /* create data for the new variable in the master problem */
            SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
            newvardata->vartype = GCG_VARTYPE_MASTER;
            newvardata->blocknr = prob;

            /* compute objective coefficient of the variable */
            objcoeff = 0;
            for ( k = 0; k < nprobvars; k++ )
            {
               if ( !SCIPisFeasZero(scip, pricerdata->solvals[k]) )
               {
                  vardata = SCIPvarGetData(probvars[k]);
                  assert(vardata->vartype == GCG_VARTYPE_PRICING);
                  assert(vardata->data.pricingvardata.origvars != NULL);
                  assert(vardata->data.pricingvardata.origvars[0] != NULL);
                  /* add quota of original variable's objcoef to the master variable's coef */
                  objcoeff += pricerdata->solvals[k] * SCIPvarGetObj(vardata->data.pricingvardata.origvars[0]);
               }
            }

            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "p_%d_%d", prob, pricerdata->nvarsprob[prob]);
            pricerdata->nvarsprob[prob]++;

            /* create variable in the master problem */
            SCIP_CALL( SCIPcreateVar(scip, &newvar, varname, 
                  0, GCGrelaxGetNIdenticalBlocks(origprob, prob), 
                  objcoeff, pricerdata->vartype, TRUE, TRUE, NULL, NULL, gcgvardeltrans, newvardata) );

            SCIPdebugMessage("found var %s with redcost %f:\n", SCIPvarGetName(newvar), 
               SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sols[j]) - pricerdata->redcostconv[prob]);
            //printf("found var %s with redcost %f:\n", SCIPvarGetName(newvar), 
            //   SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sols[j]) - pricerdata->redcostconv[prob]);
            //SCIP_CALL( SCIPprintSol(pricerdata->pricingprobs[prob], sols[j], NULL, FALSE ) );


            /* count number of non-zeros */
            newvardata->data.mastervardata.norigvars = 0;
            for ( k = 0; k < nprobvars; k++ )
            {
               if ( !SCIPisFeasZero(scip, pricerdata->solvals[k]) )
               {
                  newvardata->data.mastervardata.norigvars++;
               }
            }

            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvars), newvardata->data.mastervardata.norigvars) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvals), newvardata->data.mastervardata.norigvars) );
               
            /* number of original variables yet saved in mastervardata */
            l = 0;

            /* update variable datas */
            for ( k = 0; k < nprobvars; k++ )
            {
               if ( !SCIPisFeasZero(scip, pricerdata->solvals[k]) )
               {
                  vardata = SCIPvarGetData(probvars[k]);
                  assert(vardata->vartype == GCG_VARTYPE_PRICING);
                  assert(vardata->data.pricingvardata.origvars != NULL);
                  assert(vardata->data.pricingvardata.origvars[0] != NULL);
                  /* save in the master problem variable's data the quota of the corresponding original variable */
                  newvardata->data.mastervardata.origvars[l] = vardata->data.pricingvardata.origvars[0];
                  newvardata->data.mastervardata.origvals[l] = pricerdata->solvals[k];
                  /* save the quota in the original variable's data */
                  SCIP_CALL( GCGpricerAddMasterVarToOrigVar(scip, vardata->data.pricingvardata.origvars[0], newvar, pricerdata->solvals[k]) );
                  l++;
               }
            }
            assert(l == newvardata->data.mastervardata.norigvars);

            /* add variable and set the lazy upper bound */
            SCIP_CALL( SCIPaddPricedVar(scip, newvar, 1.0) );
            SCIPchgVarUbLazy(scip, newvar, GCGrelaxGetNIdenticalBlocks(origprob, prob));

#if 1
            BMSclearMemoryArray(mastercoefs, nmasterconss);

            /* compute coef of the variable in the master constraints and add it to the master constraints */
            for ( l = 0; l < nprobvars; l++ )
            {
               if ( !SCIPisFeasZero(scip, pricerdata->solvals[l]) )
               {
                  SCIP_CONS* linkcons;
                  int c;
                  int idx;
                  //printf("var = %s\n", SCIPvarPrintName(probvars[l]));
                  vardata = SCIPvarGetData(probvars[l]);
                  assert(vardata != NULL);
                  assert(vardata->vartype == GCG_VARTYPE_PRICING);
                  assert(vardata->data.pricingvardata.origvars != NULL);
                  assert(vardata->data.pricingvardata.origvars[0] != NULL);
                  vardata = SCIPvarGetData(vardata->data.pricingvardata.origvars[0]);
                  assert(vardata != NULL);
                  assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
                  assert(vardata->data.origvardata.coefs != NULL);
                  
                  //printf("ncoefs = %d\n", vardata->data.origvardata.ncoefs);
                  for ( c = 0; c < vardata->data.origvardata.ncoefs; c++ )
                  {
                     assert(!SCIPisZero(scip, vardata->data.origvardata.coefs[c]));
                     SCIP_CALL( SCIPgetTransformedCons(scip, vardata->data.origvardata.linkconss[c], 
                           &linkcons) );

                     idx = (int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, linkcons);
                     assert(0 <= idx && idx < nmasterconss);
                     assert(masterconss[idx] == linkcons);
                     mastercoefs[idx] += vardata->data.origvardata.coefs[c] * pricerdata->solvals[l];
                  }
                  
               }
            }
            for ( k = 0; k < nmasterconss; k++ )
            {
               if ( !SCIPisZero(scip, mastercoefs[k]) )
               {
                  /* add variable to the correponding master constraints */
                  SCIP_CALL( SCIPaddCoefLinear(scip, masterconss[k], newvar, mastercoefs[k]) );
               }
            }
#else
            /* compute coef of the variable in the master constraints and add it to the master constraints */
            for ( k = 0; k < nmasterconss; k++ )
            {
               conscoeff = 0;
               for ( l = 0; l < nprobvars; l++ )
               {
                  if ( !SCIPisFeasZero(scip, pricerdata->solvals[l]) )
                  {
                     SCIP_CONS* linkcons;
                     int c;
                     //printf("var = %s\n", SCIPvarPrintName(probvars[l]));
                     vardata = SCIPvarGetData(probvars[l]);
                     assert(vardata != NULL);
                     assert(vardata->vartype == GCG_VARTYPE_PRICING);
                     assert(vardata->data.pricingvardata.origvars != NULL);
                     assert(vardata->data.pricingvardata.origvars[0] != NULL);
                     vardata = SCIPvarGetData(vardata->data.pricingvardata.origvars[0]);
                     assert(vardata != NULL);
                     assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
                     assert(vardata->data.origvardata.coefs != NULL);
                        
                     //printf("ncoefs = %d\n", vardata->data.origvardata.ncoefs);
                     for ( c = 0; c < vardata->data.origvardata.ncoefs; c++ )
                     {
                        assert(!SCIPisZero(scip, vardata->data.origvardata.coefs[c]));
                        //printf("linkconns[c] = %s, masterconss[k] = %s\n", SCIPconsGetName(vardata->data.origvardata.linkconss[c]), SCIPconsGetName(masterconss[k]));
                        //printf("linkconns[c] = %p, masterconss[k] = %p\n", vardata->data.origvardata.linkconss[c], masterconss[k]);
                        SCIP_CALL( SCIPgetTransformedCons(scip, vardata->data.origvardata.linkconss[c], 
                              &linkcons) );
                        if ( linkcons == masterconss[k] )
                        {
                           conscoeff += vardata->data.origvardata.coefs[c] * pricerdata->solvals[l];
                           //printf("coeff = %f\n", conscoeff);
                           break;
                        }
                           
                     }

                  }
               }

               /* add variable to the correponding master constraints */
               SCIP_CALL( SCIPaddCoefLinear(scip, masterconss[k], newvar, conscoeff) );

            }
#endif

            /* compute coef of the variable in the cuts and add it to the cuts */
            for ( k = 0; k < nmastercuts; k++ )
            {
               /* get columns of the cut and their coefficients */
               cols = SCIProwGetCols(origcuts[k]);
               consvals = SCIProwGetVals(origcuts[k]);

               conscoeff = 0;

               for ( l = 0; l < SCIProwGetNNonz(origcuts[k]); l++ )
               {
                  //printf("l = %d\n", l);
                  var = SCIPcolGetVar(cols[l]);
                  vardata = SCIPvarGetData(var);

                  assert(vardata != NULL);
                  assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
                  if ( vardata->blocknr == prob )
                  {
                     //printf("vardata->blocknr == prob\n");
                     /*printf("solval of var %s is %f\n", 
                       SCIPvarGetName(vardata->data.origvardata.pricingvar),
                       SCIPgetSolVal(pricerdata->pricingprobs[prob], sols[j],
                       vardata->data.origvardata.pricingvar));*/
                     assert(vardata->data.origvardata.pricingvar != NULL);
                     if ( !SCIPisFeasZero(scip, SCIPgetSolVal(pricerdata->pricingprobs[prob], sols[j], 
                              vardata->data.origvardata.pricingvar)) )
                     {
                        conscoeff += ( consvals[l] * SCIPgetSolVal(pricerdata->pricingprobs[prob], sols[j], 
                              vardata->data.origvardata.pricingvar));
                     }
                  }

               }
               
               if ( !SCIPisFeasZero(scip, conscoeff) )
               {
                  SCIP_CALL( SCIPaddVarToRow(scip , mastercuts[k], newvar, conscoeff) );
                  //printf("new variable has coef = %f in cut %s:\n", conscoeff, SCIProwGetName(mastercuts[k]));
                  //SCIP_CALL( SCIPprintRow(origprob, mastercuts[k], NULL) );
               }
            }

            /* add variable to convexity constraint */
            SCIP_CALL( SCIPaddCoefLinear(scip, GCGrelaxGetConvCons(origprob, prob), newvar, 1) );

            //printf("created Var %s\n", SCIPvarGetName(newvar));
#if 0            
            /* check whether the created variable already existed */
            SCIP_CALL( checkNewVar(scip, newvar, SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sols[j]) - pricerdata->redcostconv[prob], pricerdata->redcostconv[prob]) );
#endif
            SCIPreleaseVar(scip, &newvar);
         }
      }

      SCIP_CALL( SCIPstartClock(scip, pricerdata->freeclock) );
      SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[prob]) );
      SCIP_CALL( SCIPstopClock(scip, pricerdata->freeclock) );
   }

   SCIPfreeBufferArray(scip, &tmpconvredcost);
   SCIPfreeBufferArray(scip, &permu);

   SCIPfreeBufferArray(scip, &mastercoefs);

   SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );

   SCIPdebugMessage("Pricing: found %d new vars\n", nfoundvars);

   return SCIP_OKAY;
}



/*
 * Callback methods of variable pricer
 */


/** destructor of variable pricer to free user data (called when SCIP is exiting) */

static
SCIP_DECL_PRICERFREE(pricerFreeGcg)
{ 
   SCIP_PRICERDATA* pricerdata;  

   assert(scip != NULL);
  
   /* get pricerdata */
   pricerdata = SCIPpricerGetData(pricer);

   /* free memory for pricerdata*/
   if ( pricerdata != NULL )
   {
      SCIPfreeMemory(scip, &pricerdata);
   }
   
   SCIPpricerSetData(pricer, NULL);
   return SCIP_OKAY;
}



/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
static
SCIP_DECL_PRICERINITSOL(pricerInitsolGcg)
{  
   SCIP_PRICERDATA* pricerdata;
   int i;
   SCIP* origprob;
   SCIP_VAR** vars;
   int nvars;
   int v;
   SCIP_VARDATA* vardata;
   SCIP_Bool discretization;
   SCIP_CONS** masterconss;
   int nmasterconss;


   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   /* set last used pricingprobnr to 0 */
   pricerdata->pricingprobnr = 0;

   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);

   /* init array containing all pricing problems */
   pricerdata->npricingprobs = GCGrelaxGetNPricingprobs(origprob);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->pricingprobs), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nvarsprob), pricerdata->npricingprobs) );

   for ( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if ( GCGrelaxIsPricingprobRelevant(origprob, i) )
      {
         pricerdata->pricingprobs[i] = GCGrelaxGetPricingprob(origprob, i);
      }
      else
      {
         pricerdata->pricingprobs[i] = NULL;
      }
      pricerdata->nvarsprob[i] = 0;
   }
   
   /* alloc memory for arrays of reduced cost */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcost), nmasterconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcostconv), pricerdata->npricingprobs) );

   /* alloc memory for solution values of variables in pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->solvals), SCIPgetNOrigVars(origprob)) );

   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->subsolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->redcostclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->redcostsolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->redcostpresolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->farkasclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->farkassolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->farkaspresolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->owneffortclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->freeclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->transformclock)) );

   pricerdata->solvedsubmips = 0;
   pricerdata->calls = 0;
   pricerdata->redcostcalls = 0;
   pricerdata->farkascalls = 0;
   pricerdata->initcalls = 0;

   /* set variable type for master variables */
   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if ( discretization )
   {
      pricerdata->vartype = SCIP_VARTYPE_INTEGER;
   }
   else
   {
      pricerdata->vartype = SCIP_VARTYPE_CONTINUOUS;
   }

   /* for variables in the original problem that do not belong to any block, 
    * create the corresponding variable in the master problem */
   vars = SCIPgetVars(origprob);
   nvars = SCIPgetNVars(origprob);
   for ( v = 0; v < nvars; v++ )
   {
      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      if ( vardata->blocknr == -1 )
      {
         SCIP_VAR* newvar;
         SCIP_VARDATA* newvardata;

         assert(vardata->data.origvardata.pricingvar == NULL);

         SCIPdebugMessage("var %s is in no block!\n", SCIPvarGetName(vars[v]));

         /* create vardata */
         SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
         newvardata->vartype = GCG_VARTYPE_MASTER;
         newvardata->blocknr = -1;
         newvardata->data.mastervardata.norigvars = 2;

         /* save corresoponding origvar */
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, 
               &(newvardata->data.mastervardata.origvars), 2) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, 
               &(newvardata->data.mastervardata.origvals), 2) );
         newvardata->data.mastervardata.origvars[0] = vars[v];
         newvardata->data.mastervardata.origvals[0] = 1.0;
         newvardata->data.mastervardata.origvars[1] = vars[v];
         newvardata->data.mastervardata.origvals[1] = 0.0;

         /* create variable in the master problem */
         SCIP_CALL( SCIPcreateVar(scip, &newvar, SCIPvarGetName(vars[v]), 
               SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v]), SCIPvarGetObj(vars[v]), SCIPvarGetType(vars[v]), 
               TRUE, TRUE, NULL, NULL, gcgvardeltrans, newvardata) );
         SCIPaddVar(scip, newvar);

         SCIPchgVarUbLazy(scip, newvar, SCIPvarGetUbGlobal(vars[v]));
         SCIPchgVarLbLazy(scip, newvar, SCIPvarGetLbGlobal(vars[v]));

         SCIP_CALL( GCGpricerAddMasterVarToOrigVar(scip, vars[v], newvar, 1.0) );

         /* add variable in the master to the master constraints it belongs to */
         for ( i = 0; i < vardata->data.origvardata.ncoefs; i++ )
         {
            SCIP_CONS* linkcons;
            assert(!SCIPisFeasZero(scip, vardata->data.origvardata.coefs[i]));
            SCIP_CALL( SCIPgetTransformedCons(scip, vardata->data.origvardata.linkconss[i], 
                  &linkcons) );

            SCIP_CALL( SCIPaddCoefLinear(scip, linkcons, 
                  newvar, vardata->data.origvardata.coefs[i]) );
         }
         SCIPreleaseVar(scip, &newvar);

      }
   }

   SCIP_CALL( SCIPhashmapCreate(&(pricerdata->mapcons2idx), SCIPblkmem(scip), 10 * nmasterconss) );
   for ( i = 0; i < nmasterconss; i++ )
   {
      //printf("add cons %s to hashmap: pointer %p\n", SCIPconsGetName(masterconss[i]), masterconss[i]);
      SCIP_CALL( SCIPhashmapInsert(pricerdata->mapcons2idx, masterconss[i], (void*)(size_t)i) );
      assert((int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, masterconss[i]) == i);
   }

   return SCIP_OKAY;
}



/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolGcg)
{  
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   SCIPhashmapFree(&(pricerdata->mapcons2idx));
   
   SCIPfreeMemoryArray(scip, &(pricerdata->pricingprobs));
   SCIPfreeMemoryArray(scip, &(pricerdata->redcost));
   SCIPfreeMemoryArray(scip, &(pricerdata->redcostconv));
   SCIPfreeMemoryArray(scip, &(pricerdata->solvals));
   SCIPfreeMemoryArray(scip, &(pricerdata->nvarsprob));

   printf("calls = %d\n", pricerdata->calls);
   printf("time for pricing without sub-MIPs: %f\n", SCIPgetClockTime(scip, pricerdata->owneffortclock));
   printf("solved sub-MIPs = %d\n", pricerdata->solvedsubmips);
   printf("time for solving sub-MIPs for pricing: %f\n", SCIPgetClockTime(scip, pricerdata->subsolveclock));
   printf("farkas calls = %d, redcost calls = %d\n", pricerdata->farkascalls, pricerdata->redcostcalls);
   printf("time for farkas pricing (presolving): %f\n", SCIPgetClockTime(scip, pricerdata->farkaspresolveclock));
   printf("time for farkas pricing (solving): %f\n", SCIPgetClockTime(scip, pricerdata->farkassolveclock));
   printf("time for farkas pricing (total): %f\n", SCIPgetClockTime(scip, pricerdata->farkasclock));
   printf("time for redcost pricing (presolving): %f\n", SCIPgetClockTime(scip, pricerdata->redcostpresolveclock));
   printf("time for redcost pricing (solving): %f\n", SCIPgetClockTime(scip, pricerdata->redcostsolveclock));
   printf("time for redcost pricing (total): %f\n", SCIPgetClockTime(scip, pricerdata->redcostclock));
   printf("time for transformation: %f\n", SCIPgetClockTime(scip, pricerdata->transformclock));
   printf("time for freeing sub-MIPs: %f\n", SCIPgetClockTime(scip, pricerdata->freeclock));



   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->subsolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->redcostclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->redcostsolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->redcostpresolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->farkasclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->farkassolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->farkaspresolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->owneffortclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->freeclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->transformclock)) );

   
   return SCIP_OKAY;
}




/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostGcg)
{  
   SCIP_RETCODE retcode;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   assert(pricerdata != NULL);

   SCIP_CALL( SCIPstartClock(scip, pricerdata->redcostclock) );

   //printf("pricerredcost\n");
   retcode = performPricing(scip, pricer, GCG_PRICETYPE_REDCOST);

   SCIP_CALL( SCIPstopClock(scip, pricerdata->redcostclock) );

   *result = SCIP_SUCCESS;

   return retcode;
}




static
SCIP_DECL_PRICERFARKAS(pricerFarkasGcg)
{  
   SCIP_RETCODE retcode;
   SCIP_PRICERDATA* pricerdata;
   
   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   assert(pricerdata != NULL);

   SCIP_CALL( SCIPstartClock(scip, pricerdata->farkasclock) );

   //printf("pricerfarkas\n");
   retcode = performPricing(scip, pricer, GCG_PRICETYPE_FARKAS);

   SCIP_CALL( SCIPstopClock(scip, pricerdata->farkasclock) );

   return retcode;
}

/* define not used callbacks as NULL */
#define pricerInitGcg NULL
#define pricerExitGcg NULL


/*
 * variable pricer specific interface methods
 */

/** creates the gcg variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origprob            /**< SCIP data structure of the original problem */
   )
{
   SCIP_PRICERDATA* pricerdata;

   SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );
   pricerdata->origprob = origprob;

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricer(scip, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerFreeGcg, pricerInitGcg, pricerExitGcg, 
         pricerInitsolGcg, pricerExitsolGcg, pricerRedcostGcg, pricerFarkasGcg,
         pricerdata) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxvarsroundredcost",
         "maximal number of variables created in one redcost pricing round",
         &pricerdata->maxvarsroundredcost, TRUE, DEFAULT_MAXVARSROUNDREDCOST, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxvarsroundfarkas",
         "maximal number of variables created in one farkas pricing round",
         &pricerdata->maxvarsroundfarkas, TRUE, DEFAULT_MAXVARSROUNDFARKAS, 1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** returns the pointer to the scip instance representing the original problem */
SCIP* GCGpricerGetOrigprob(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   return pricerdata->origprob;
}
