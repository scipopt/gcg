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
//#define DEBUG_PRICING
//#define DEBUG_PRICING_ALL_OUTPUT
//#define CHECKNEWVAR
//#define CHECKVARBOUNDS
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
#define DEFAULT_MAXSUCCESSFULMIPSREDCOST INT_MAX
#define DEFAULT_MAXROUNDSREDCOST INT_MAX
#define DEFAULT_MAXSOLSPROB INT_MAX
#define DEFAULT_CHECKSOLS TRUE
#define DEFAULT_USEHEURPRICING FALSE
#define DEFAULT_ONLYPOSCONV FALSE


/*
 * Data structures
 */


/** variable pricer data */
struct SCIP_PricerData
{
   int npricingprobs;              /* number of pricing problems */
   SCIP** pricingprobs;            /* pointers to the pricing problems */
   SCIP_Real* dualsol;             /* array of dual solutions for the constraints of the master problem */
   SCIP_Real* dualsolconv;         /* array of dual solutions for the convexity constraints */
   SCIP_Real* objsums;             /* sum of objective coefficients for the pricing problems */
   SCIP_Real* objfactors;
   SCIP* origprob;                 /* the original program */
   SCIP_Real* solvals;             /* solution values of variables in the pricing problems */
   int* nvarsprob;                 /* number of variables created by the pricing probs */
   SCIP_Longint currnodenr;
   SCIP_HASHMAP* mapcons2idx;

   /** variables used for statistics */
   SCIP_CLOCK* subsolveclock;
   SCIP_CLOCK* heursolveclock;
   SCIP_CLOCK* redcostclock;
   SCIP_CLOCK* redcostsolveclock;
   SCIP_CLOCK* redcostpresolveclock;
   SCIP_CLOCK* farkasclock;
   SCIP_CLOCK* farkassolveclock;
   SCIP_CLOCK* farkaspresolveclock;
   SCIP_CLOCK* owneffortclock;
   SCIP_CLOCK* freeclock;
   SCIP_CLOCK* transformclock;
   int solvedsubmipsoptimal;
   int solvedsubmipsheur;
   int calls;
   int farkascalls;
   int redcostcalls;
   int nallsame;

   /** parameter values */
   SCIP_VARTYPE vartype;           /* vartype of created master variables */
   int maxvarsroundfarkas;
   int maxvarsroundredcost;
   int maxsuccessfulmipsredcost;
   int maxroundsredcost;
   int maxsolsprob;
   int nroundsredcost;
   SCIP_Bool checksols;
   SCIP_Bool useheurpricing;
   SCIP_Bool onlyposconv;
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

#ifdef CHECKNEWVAR
static
SCIP_RETCODE checkNewVar(
   SCIP*                 scip,
   SCIP_VAR*             newvar,
   SCIP_Real             redcost,
   SCIP_Real             dualsolconv
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
   for( v = 0; v < nvars; v++ )
   {
      assert(vars[v] != NULL);

      /* vars[v] is the new variable itself */
      if( vars[v] == newvar )
         continue;

      /* vars[v] has a different objective function value, may not be equal to newvar */
      if( !SCIPisEQ(scip, SCIPvarGetObj(vars[v]), SCIPvarGetObj(newvar)) )
         continue;

      vardata = SCIPvarGetData(vars[v]);
      newvardata = SCIPvarGetData(newvar);
      assert(vardata != NULL);
      assert(newvardata != NULL);

      /* vars[v] belongs to a different block, may not be equal to newvar */
      if( vardata->blocknr != newvardata->blocknr )
         continue;

      assert(vardata->vartype == GCG_VARTYPE_MASTER);
      assert(newvardata->vartype == GCG_VARTYPE_MASTER);

      /* vars[v] belongs to a different block, may not be equal to newvar */
      if( vardata->data.mastervardata.norigvars != newvardata->data.mastervardata.norigvars )
         continue;

      /* compare the parts of the original variables contained in vars[i] and newvar */
      for( i = 0; i < vardata->data.mastervardata.norigvars; i++ )
      {
         /* original variables are not equal */
         if( vardata->data.mastervardata.origvars[i] != newvardata->data.mastervardata.origvars[i] )
            break;
         /* the values of the original variable are not equal */
         if( !SCIPisEQ(scip, vardata->data.mastervardata.origvals[i], newvardata->data.mastervardata.origvals[i]) )
            break;
      }

      if( i == vardata->data.mastervardata.norigvars )
      {
         printf("var %s is equal to var %s! solval = %f, ub = %g, lazyub = %g, redcost = %f, lpsolstat = %d, dualsolconv = %f\n", 
            SCIPvarGetName(newvar), SCIPvarGetName(vars[v]), SCIPgetSolVal(scip, NULL, vars[v]), 
            SCIPvarGetUbLocal(vars[v]), SCIPvarGetUbLazy(vars[v]), redcost, SCIPgetLPSolstat(scip), dualsolconv);
      }
   }

   return SCIP_OKAY;
}
#endif

#ifndef CHECKVARBOUNDS
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
   for( v = 0; v < nvars; v++ )
   {
      assert(vars[v] != NULL);

      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->data.origvardata.pricingvar != NULL || vardata->blocknr == -1);

      if( vardata->blocknr == -1 )
         continue;

      if( !GCGrelaxIsPricingprobRelevant(origscip, vardata->blocknr) || GCGrelaxGetNIdenticalBlocks(origscip, vardata->blocknr) != 1 ) 
         continue;

      if( SCIPvarGetUbLocal(vars[v]) != SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar) )
      {
         printf("var %s: orig upper bound = %f, pricing upper bound = %f!\n", SCIPvarGetName(vars[v]),
            SCIPvarGetUbLocal(vars[v]), SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar) );
      }
      if( SCIPvarGetLbLocal(vars[v]) != SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar) )
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
   if( vardata->data.origvardata.maxmastervars == vardata->data.origvardata.nmastervars )
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

static
void updateObjFactor(
   SCIP*                 scip,
   int                   probnr,
   SCIP_Real             objvalue
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(0 <= probnr && probnr <= pricerdata->npricingprobs);
   
   if( !SCIPisZero(scip, pricerdata->objsums[probnr]) )
      pricerdata->objfactors[probnr] = 0.9 * pricerdata->objfactors[probnr] + 0.1 * (objvalue / pricerdata->objsums[probnr]);
}

static
SCIP_RETCODE checkSolNew(
   SCIP*                 scip,
   SCIP_SOL**            sols,
   int                   idx,
   int                   prob,
   SCIP_Bool*            isnew
   )
{
   SCIP* origprob;
   SCIP* pricingprob;
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_VAR** probvars;
   int nprobvars;
   SCIP_Real* newvals;

   int s;
   int i;

   assert(scip != NULL);
   assert(sols != NULL);
   assert(sols[idx] != NULL);
   assert(isnew != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(0 <= prob && prob <= pricerdata->npricingprobs);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   pricingprob = pricerdata->pricingprobs[prob];
   assert(pricingprob != NULL);

   probvars = SCIPgetVars(pricingprob);
   nprobvars = SCIPgetNVars(pricingprob);

   *isnew = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &newvals, nprobvars) );

   SCIP_CALL( SCIPgetSolVals(pricingprob, sols[idx], nprobvars, probvars, newvals) );

   for( s = 0; s < idx && *isnew == TRUE; s++ )
   {
      assert(sols[s] != NULL);
      assert(SCIPisLE(scip, SCIPgetSolOrigObj(pricingprob, sols[s]), SCIPgetSolOrigObj(pricingprob, sols[idx])));
      if( !SCIPisEQ(scip, SCIPgetSolOrigObj(pricingprob, sols[s]), SCIPgetSolOrigObj(pricingprob, sols[idx])) )
         continue;

      if( SCIPsolGetOrigin(sols[s]) != SCIP_SOLORIGIN_ORIGINAL && SCIPsolGetOrigin(sols[idx]) != SCIP_SOLORIGIN_ORIGINAL )
         continue;

      for( i = 0; i < nprobvars; i++ )
      {
         if( !SCIPisEQ(scip, SCIPgetSolVal(pricingprob, sols[s], probvars[i]), newvals[i]) )
            break;
      }
      if( i == nprobvars )
      {
         //printf("sol %d (origin %d) equals sol %d (origin %d):\n", idx, SCIPsolGetOrigin(sols[idx]), s, SCIPsolGetOrigin(sols[s]));
         //SCIP_CALL( SCIPprintSol(pricerdata->pricingprobs[prob], sols[idx], NULL, FALSE ) );
         //SCIP_CALL( SCIPprintSol(pricerdata->pricingprobs[prob], sols[s], NULL, FALSE ) );
         *isnew = FALSE;
      }
   }

   SCIPfreeBufferArray(scip, &newvals);

   return SCIP_OKAY;
}

static
SCIP_RETCODE setPricingObjs(
   SCIP*                 scip,
   GCG_PRICETYPE         pricetype
   )
{
   SCIP* origprob;
   SCIP_VARDATA* vardata;
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS** origconss;
   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_VAR** probvars;
   int nprobvars;

   SCIP_ROW** mastercuts;
   int nmastercuts;
   SCIP_ROW** origcuts;
   int norigcuts;
   SCIP_COL** cols;
   SCIP_Real* consvals;
   SCIP_Real dualsol;

   SCIP_VAR** consvars;
   int nconsvars;

   int i;
   int j;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   /* get the constraints of the master problem and the corresponding constraints in the original problem */
   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);
   origconss = GCGrelaxGetLinearOrigMasterConss(origprob);

   /* set objective value of all variables in the pricing problems to 0 (for farkas pricing) /
    * to the original objective of the variable (for redcost pricing) */
   for( i = 0; i < pricerdata->npricingprobs; i++)
   {
      if( pricerdata->pricingprobs[i] == NULL )
         continue;
      probvars = SCIPgetVars(pricerdata->pricingprobs[i]);
      nprobvars = SCIPgetNVars(pricerdata->pricingprobs[i]);

      pricerdata->objsums[i] = 0.0;
      
      for( j = 0; j < nprobvars; j++ )
      {
         if( pricetype == GCG_PRICETYPE_FARKAS )
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
            pricerdata->objsums[i] += SCIPvarGetObj(vardata->data.pricingvardata.origvars[0]);
         }  
      }
   }

   /* compute reduced cost and update objectives in the pricing problems */
   for( i = 0; i < nmasterconss; i++ )
   {
      /* farkas pricing */
      if( pricetype == GCG_PRICETYPE_REDCOST )
         pricerdata->dualsol[i] = SCIPgetDualsolLinear(scip, masterconss[i]);
      /* redcost pricing */
      else 
      {
         assert(pricetype == GCG_PRICETYPE_FARKAS);
         pricerdata->dualsol[i] = SCIPgetDualfarkasLinear(scip, masterconss[i]);
      }
      if( !SCIPisZero(scip, pricerdata->dualsol[i]) )
      {
         /* for all variables in the constraint, modify the objective of the corresponding variable in a pricing problem */
         consvars = SCIPgetVarsLinear(origprob, origconss[i]);
         consvals = SCIPgetValsLinear(origprob, origconss[i]);
         nconsvars = SCIPgetNVarsLinear(origprob, origconss[i]);
         for( j = 0; j < nconsvars; j++ )
         {
            vardata = SCIPvarGetData(consvars[j]);
            assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
            if( vardata->blocknr != -1 && pricerdata->pricingprobs[vardata->blocknr] != NULL )
            {
               assert(vardata->data.origvardata.pricingvar != NULL);
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[vardata->blocknr], 
                     vardata->data.origvardata.pricingvar, -1.0 * pricerdata->dualsol[i] * consvals[j]) );
               pricerdata->objsums[vardata->blocknr] -= pricerdata->dualsol[i] * consvals[j];
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
   for( i = 0; i < nmastercuts; i++ )
   {
      /* farkas pricing */
      if( pricetype == GCG_PRICETYPE_REDCOST )
         dualsol = SCIProwGetDualsol(mastercuts[i]);
      /* redcost pricing */
      else
      {    
         assert(pricetype == GCG_PRICETYPE_FARKAS);
         dualsol = SCIProwGetDualfarkas(mastercuts[i]);
      }
      if( !SCIPisZero(scip, dualsol) )
      {
         /* get columns and vals of the cut */
         nconsvars = SCIProwGetNNonz(origcuts[i]);
         cols = SCIProwGetCols(origcuts[i]);
         consvals = SCIProwGetVals(origcuts[i]);

         /* get the variables corresponding to the columns in the cut */
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
         for( j = 0; j < nconsvars; j++ )
            consvars[j] = SCIPcolGetVar(cols[j]);

         /* for all variables in the cut, modify the objective of the corresponding variable in a pricing problem */
         for( j = 0; j < nconsvars; j++ )
         {
            vardata = SCIPvarGetData(consvars[j]);
            assert(vardata != NULL);
            assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
            if( vardata->blocknr != -1 && pricerdata->pricingprobs[vardata->blocknr] != NULL )
            {
               assert(vardata->data.origvardata.pricingvar != NULL);
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[vardata->blocknr], 
                     vardata->data.origvardata.pricingvar, -1.0 * dualsol * consvals[j]) );
               pricerdata->objsums[vardata->blocknr] -= dualsol * consvals[j];
            }
         }
         SCIPfreeBufferArray(scip, &consvars);
      }
   }

   /* get dual solutions / farkas values of the convexity constraints */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {

      assert( GCGrelaxIsPricingprobRelevant(origprob, i) 
         == (GCGrelaxGetConvCons(origprob, i) != NULL) );
      if( !GCGrelaxIsPricingprobRelevant(origprob, i) )
      {
         pricerdata->dualsolconv[i] = -1.0 * SCIPinfinity(scip);
         continue;
      }  
      if( pricetype == GCG_PRICETYPE_REDCOST )
         pricerdata->dualsolconv[i] = SCIPgetDualsolLinear(scip, GCGrelaxGetConvCons(origprob, i));
      else
      {
         assert(pricetype == GCG_PRICETYPE_FARKAS);
         pricerdata->dualsolconv[i] = SCIPgetDualfarkasLinear(scip, GCGrelaxGetConvCons(origprob, i));
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE createNewMasterVar(
   SCIP*                 scip,
   SCIP_SOL*             sol,
   int                   prob
   )
{
   SCIP* origprob;
   SCIP_VARDATA* vardata;
   SCIP_VARDATA* newvardata;
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_Real* mastercoefs;
   SCIP_VAR** probvars;
   int nprobvars;
   char varname[SCIP_MAXSTRLEN];

   SCIP_ROW** mastercuts;
   int nmastercuts;
   SCIP_ROW** origcuts;
   int norigcuts;
   SCIP_COL** cols;
   SCIP_Real conscoef;
   SCIP_VAR* var;
   SCIP_Real* consvals;
   SCIP_Real* solvals;
   SCIP_Real objcoeff;
   SCIP_VAR* newvar;

   SCIP_CONS* linkcons;
   int c;
   int idx;

   int i;
   int j;

   assert(scip != NULL);
   assert(sol != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);
   
   /* get variables of the pricing problem and their values in the current solution */
   probvars = SCIPgetOrigVars(pricerdata->pricingprobs[prob]);
   nprobvars = SCIPgetNOrigVars(pricerdata->pricingprobs[prob]);

   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nprobvars) );
   SCIP_CALL( SCIPgetSolVals(pricerdata->pricingprobs[prob], sol, nprobvars, probvars, solvals) );
            
   /* create data for the new variable in the master problem */
   SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
   newvardata->vartype = GCG_VARTYPE_MASTER;
   newvardata->blocknr = prob;

   /* compute objective coefficient of the variable */
   objcoeff = 0;
   for( i = 0; i < nprobvars; i++ )
   {
      if( !SCIPisZero(scip, solvals[i]) )
      {
         vardata = SCIPvarGetData(probvars[i]);
         assert(vardata->vartype == GCG_VARTYPE_PRICING);
         assert(vardata->data.pricingvardata.origvars != NULL);
         assert(vardata->data.pricingvardata.origvars[0] != NULL);
         /* add quota of original variable's objcoef to the master variable's coef */
         objcoeff += solvals[i] * SCIPvarGetObj(vardata->data.pricingvardata.origvars[0]);
      }
   }

   (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "p_%d_%d", prob, pricerdata->nvarsprob[prob]);
   pricerdata->nvarsprob[prob]++;
   
   /* create variable in the master problem */
   SCIP_CALL( SCIPcreateVar(scip, &newvar, varname, 0, INT_MAX /*GCGrelaxGetNIdenticalBlocks(origprob, prob)*/, 
         objcoeff, pricerdata->vartype, TRUE, TRUE, NULL, NULL, gcgvardeltrans, newvardata) );

   SCIPdebugMessage("found var %s with redcost %f:\n", SCIPvarGetName(newvar), 
      SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sol) - pricerdata->dualsolconv[prob]);

   /* count number of non-zeros */
   newvardata->data.mastervardata.norigvars = 0;
   for( i = 0; i < nprobvars; i++ )
      if( !SCIPisZero(scip, solvals[i]) )
         newvardata->data.mastervardata.norigvars++;
   
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvars), newvardata->data.mastervardata.norigvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvals), newvardata->data.mastervardata.norigvars) );
               
   /* number of original variables yet saved in mastervardata */
   j = 0;

   /* update variable datas */
   for( i = 0; i < nprobvars; i++ )
   {
      if( !SCIPisZero(scip, solvals[i]) )
      {
         vardata = SCIPvarGetData(probvars[i]);
         assert(vardata->vartype == GCG_VARTYPE_PRICING);
         assert(vardata->data.pricingvardata.origvars != NULL);
         assert(vardata->data.pricingvardata.origvars[0] != NULL);
         /* save in the master problem variable's data the quota of the corresponding original variable */
         newvardata->data.mastervardata.origvars[j] = vardata->data.pricingvardata.origvars[0];
         newvardata->data.mastervardata.origvals[j] = solvals[i];
         /* save the quota in the original variable's data */
         SCIP_CALL( GCGpricerAddMasterVarToOrigVar(scip, vardata->data.pricingvardata.origvars[0], newvar, solvals[i]) );
         j++;
      }
   }
   assert(j == newvardata->data.mastervardata.norigvars);

   /* add variable and set the lazy upper bound */
   SCIP_CALL( SCIPaddPricedVar(scip, newvar, 1.0) );
   SCIPchgVarUbLazy(scip, newvar, GCGrelaxGetNIdenticalBlocks(origprob, prob));

   SCIP_CALL( SCIPallocBufferArray(scip, &mastercoefs, nmasterconss) );
   BMSclearMemoryArray(mastercoefs, nmasterconss);

   /* compute coef of the variable in the master constraints */
   for( i = 0; i < nprobvars; i++ )
   {
      if( !SCIPisZero(scip, solvals[i]) )
      {
         vardata = SCIPvarGetData(probvars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_PRICING);
         assert(vardata->data.pricingvardata.origvars != NULL);
         assert(vardata->data.pricingvardata.origvars[0] != NULL);
         vardata = SCIPvarGetData(vardata->data.pricingvardata.origvars[0]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         assert(vardata->data.origvardata.coefs != NULL || vardata->data.origvardata.ncoefs == 0);
                  
         /* for each coef, add coef * solval to the coef of the new variable for the corresponding constraint */
         for( c = 0; c < vardata->data.origvardata.ncoefs; c++ )
         {
            assert(!SCIPisZero(scip, vardata->data.origvardata.coefs[c]));
            SCIP_CALL( SCIPgetTransformedCons(scip, vardata->data.origvardata.linkconss[c], 
                  &linkcons) );

            idx = (int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, linkcons);
            assert(0 <= idx && idx < nmasterconss);
            assert(masterconss[idx] == linkcons);
            mastercoefs[idx] += vardata->data.origvardata.coefs[c] * solvals[i];
         }
                  
      }
   }
   /* add the variables to the master constraints */
   for( i = 0; i < nmasterconss; i++ )
   {
      if( !SCIPisZero(scip, mastercoefs[i]) )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, masterconss[i], newvar, mastercoefs[i]) );
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

   /* compute coef of the variable in the cuts and add it to the cuts */
   for( i = 0; i < nmastercuts; i++ )
   {
      /* get columns of the cut and their coefficients */
      cols = SCIProwGetCols(origcuts[i]);
      consvals = SCIProwGetVals(origcuts[i]);

      conscoef = 0;

      for( j = 0; j < SCIProwGetNNonz(origcuts[i]); j++ )
      {
         var = SCIPcolGetVar(cols[j]);
         vardata = SCIPvarGetData(var);

         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         if( vardata->blocknr == prob )
         {
            assert(vardata->data.origvardata.pricingvar != NULL);
            if( !SCIPisZero(scip, SCIPgetSolVal(pricerdata->pricingprobs[prob], sol, 
                     vardata->data.origvardata.pricingvar)) )
            {
               conscoef += ( consvals[j] * SCIPgetSolVal(pricerdata->pricingprobs[prob], sol, 
                     vardata->data.origvardata.pricingvar));
            }
         }

      }
               
      if( !SCIPisZero(scip, conscoef) )
      {
         SCIP_CALL( SCIPaddVarToRow(scip , mastercuts[i], newvar, conscoef) );
         //printf("new variable has coef = %f in cut %s:\n", conscoef, SCIProwGetName(mastercuts[i]));
         //SCIP_CALL( SCIPprintRow(origprob, mastercuts[i], NULL) );
      }
   }

   /* add variable to convexity constraint */
   SCIP_CALL( SCIPaddCoefLinear(scip, GCGrelaxGetConvCons(origprob, prob), newvar, 1) );

   SCIPfreeBufferArray(scip, &mastercoefs);
   SCIPfreeBufferArray(scip, &solvals);

#ifdef CHECKNEWVAR
   /* check whether the created variable already existed */
   SCIP_CALL( checkNewVar(scip, newvar, SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sol) - pricerdata->dualsolconv[prob], pricerdata->dualsolconv[prob]) );
#endif

   SCIPreleaseVar(scip, &newvar);

   return SCIP_OKAY;
}


/* performs the pricing routine, gets the type of pricing that should be done: farkas or redcost pricing */
static
SCIP_RETCODE performPricing(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_PRICER*          pricer,             /* the pricer */
   GCG_PRICETYPE         pricetype,          /* type of the pricing */
   SCIP_RESULT*          result,             /* result pointer */
   SCIP_Real*            lowerbound          /* lowerbound pointer */
   )
{
   SCIP_PRICERDATA* pricerdata;            /* the data of the pricer */
   SCIP* origprob;

   int i;
   int j;
   int prob;

   int nfoundvars;
   int nfoundvarsprob;
   int successfulmips;

   int nsols;
   SCIP_SOL** sols;

   SCIP_Real timelimit;

   SCIP_Real bestredcost;
   SCIP_Bool bestredcostvalid;
   SCIP_Bool newsol;

   SCIP_Real* tmpconvdualsol;
   int* permu;

   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   assert(result != NULL || pricetype == GCG_PRICETYPE_FARKAS);
   assert(lowerbound != NULL || pricetype == GCG_PRICETYPE_FARKAS);

#ifdef DEBUG_PRICING
   if( pricetype == GCG_PRICETYPE_REDCOST || SCIPgetNVars(scip) % 50 == 0 )
   {
      printf("nvars = %d, current lowerbound = %g, time = %f, node = %lld\n", SCIPgetNVars(scip), 
         SCIPgetLPObjval(scip), SCIPgetSolvingTime(scip), SCIPgetNNodes(scip));
   }
#endif

   /* check whether pricing can be aborted: if objective value is always integral
    * and the current node's current lowerbound rounded up equals the 
    * current lp objective value rounded up we don't need to continue pricing
    * since the best possible feasible solution must have at least this value
    */
   if( SCIPisObjIntegral(scip) && pricetype == GCG_PRICETYPE_REDCOST && 
      SCIPceil(scip, SCIPgetNodeDualbound(scip, SCIPgetCurrentNode(scip))) 
      == SCIPceil(scip, SCIPgetLPObjval(scip)) )
   {
#ifdef DEBUG_PRICING
      printf("pricing aborted due to integral objective: node LB = %g, LP obj = %g\n", 
         SCIPgetNodeDualbound(scip, SCIPgetCurrentNode(scip)), SCIPgetLPObjval(scip));
#endif
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( pricetype == GCG_PRICETYPE_REDCOST )
   {
      pricerdata->redcostcalls++;
      *result = SCIP_SUCCESS;
   }
   if( pricetype == GCG_PRICETYPE_FARKAS )
      pricerdata->farkascalls++;

   SCIP_CALL( SCIPstartClock(scip, pricerdata->owneffortclock) );
   pricerdata->calls++;
   nfoundvars = 0;
   successfulmips = 0;

#ifndef CHECKVARBOUNDS
   SCIP_CALL( checkVarBounds(scip) );
#endif
   /* set objectives of the variables in the pricing sub-MIPs */
   SCIP_CALL( setPricingObjs(scip, pricetype) );

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpconvdualsol, pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &permu, pricerdata->npricingprobs) );
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      assert(GCGrelaxIsPricingprobRelevant(origprob, i)
         == (pricerdata->dualsolconv[i] != -1.0 * SCIPinfinity(scip)));

      permu[i] = i;
      tmpconvdualsol[i] = pricerdata->dualsolconv[i];
#if 0
      if( pricetype == GCG_PRICETYPE_REDCOST )
         tmpconvdualsol[i] -= pricerdata->objfactors[i] * pricerdata->objsums[i];
#endif
   }
   SCIPsortDownRealInt(tmpconvdualsol, permu, pricerdata->npricingprobs);
#if 0
   if( SCIPisEQ(scip, tmpconvdualsol[0], tmpconvdualsol[pricerdata->npricingprobs-1]) )
   {
      printf("tmpconvdualsol all %g!\n", tmpconvdualsol[0]);
      pricerdata->nallsame++;
   }
#endif
   bestredcost = 0.0;
   bestredcostvalid = FALSE;

   if( pricerdata->useheurpricing )
   {
      /* solve the pricing MIPs heuristically and check whether solutions 
       * corresponding to variables with negative reduced costs where found 
       */
      for( i = 0; i < pricerdata->npricingprobs &&
              (pricetype == GCG_PRICETYPE_REDCOST || nfoundvars < pricerdata->maxvarsroundfarkas)
              && (pricetype == GCG_PRICETYPE_FARKAS || nfoundvars < pricerdata->maxvarsroundredcost); i++)
      {
         prob = permu[i];

         if( pricerdata->pricingprobs[prob] == NULL )
            continue;

         /* set objective limit, such that only solutions with negative reduced costs are accepted */
         SCIP_CALL( SCIPsetObjlimit(pricerdata->pricingprobs[prob], pricerdata->dualsolconv[prob]) );

#ifdef DEBUG_PRICING_ALL_OUTPUT
         if( pricetype == GCG_PRICETYPE_REDCOST )
         {
            SCIP_CALL( SCIPsetIntParam(pricerdata->pricingprobs[prob], "display/verblevel", SCIP_VERBLEVEL_HIGH) );
         }
#endif

         /* set time limit */
         SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
         if( !SCIPisInfinity(scip, timelimit) )
         {
            if( timelimit - SCIPgetTotalTime(scip) > 0 )
            {
               SCIP_CALL( SCIPsetRealParam(pricerdata->pricingprobs[prob], "limits/time", 
                     timelimit - SCIPgetTotalTime(scip)) );
            }
            else
            {
               *result = SCIP_DIDNOTRUN;

               /* free transformed pricing MIPs */
               for( j = 0; j < pricerdata->npricingprobs; j++ )
               {
                  if( pricerdata->pricingprobs[j] != NULL )
                  {
                     SCIP_CALL( SCIPstartClock(scip, pricerdata->freeclock) );
                     SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[j]) );
                     SCIP_CALL( SCIPstopClock(scip, pricerdata->freeclock) );
                  }
               }
               SCIPfreeBufferArray(scip, &permu);
               SCIPfreeBufferArray(scip, &tmpconvdualsol);

               SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );
            
               return SCIP_OKAY;
            }
         }

         //SCIP_CALL( SCIPsetLongintParam(pricerdata->pricingprobs[prob], "limits/stallnodes", 50) );
         //SCIP_CALL( SCIPsetLongintParam(pricerdata->pricingprobs[prob], "limits/nodes", 200) ); 
         SCIP_CALL( SCIPsetRealParam(pricerdata->pricingprobs[prob], "limits/gap", 0.2) ); 
         //SCIP_CALL( SCIPsetIntParam(pricerdata->pricingprobs[prob], "limits/bestsol", 5) ); 

         SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );
         SCIP_CALL( SCIPstartClock(scip, pricerdata->heursolveclock) );

         /* start clock measuring the transformation effort */
         SCIP_CALL( SCIPstartClock(scip, pricerdata->transformclock) );
         SCIP_CALL( SCIPtransformProb(pricerdata->pricingprobs[prob]) );
         SCIP_CALL( SCIPstopClock(scip, pricerdata->transformclock) );

         /* start clock measuring the presolving effort */
         if( pricetype == GCG_PRICETYPE_REDCOST )
            SCIP_CALL( SCIPstartClock(scip, pricerdata->redcostpresolveclock) );
         else
            SCIP_CALL( SCIPstartClock(scip, pricerdata->farkaspresolveclock) );

         /* presolve the pricing submip */
         SCIP_CALL( SCIPpresolve(pricerdata->pricingprobs[prob]) );

         /* stop clock measuring the presolving effort, start clock measuring the solving effort */
         if( pricetype == GCG_PRICETYPE_REDCOST )
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
         if( pricetype == GCG_PRICETYPE_REDCOST )
            SCIP_CALL( SCIPstopClock(scip, pricerdata->redcostsolveclock) );
         else
            SCIP_CALL( SCIPstopClock(scip, pricerdata->farkassolveclock) );

         pricerdata->solvedsubmipsheur++;

         SCIP_CALL( SCIPstopClock(scip, pricerdata->heursolveclock) );
         SCIP_CALL( SCIPstartClock(scip, pricerdata->owneffortclock) );

#ifdef DEBUG_PRICING_ALL_OUTPUT
         if( pricetype == GCG_PRICETYPE_REDCOST )
         {
            SCIP_CALL( SCIPsetIntParam(pricerdata->pricingprobs[prob], "display/verblevel", 0) );
            SCIP_CALL( SCIPprintStatistics(pricerdata->pricingprobs[prob], NULL) );
         }
#endif

         /* so far, the pricing problem should be solved to optimality */
         assert( SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_OPTIMAL
            || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_GAPLIMIT
            || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_USERINTERRUPT 
            || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_INFEASIBLE
            || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_TIMELIMIT );
         /** @todo handle userinterrupt: set result pointer (and lowerbound), handle other solution states */

         nsols = SCIPgetNSols(pricerdata->pricingprobs[prob]);
         sols = SCIPgetSols(pricerdata->pricingprobs[prob]);
         //printf("Pricingprob %d has found %d sols!\n", prob, nsols);

         nfoundvarsprob = 0;

         for( j = 0; j < nsols && nfoundvarsprob <= pricerdata->maxsolsprob &&
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
            if( SCIPisNegative(scip, SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sols[j]) - pricerdata->dualsolconv[prob]) )
            {
               if( pricerdata->checksols )
               {
                  SCIP_CALL( checkSolNew(scip, sols, j, prob, &newsol) );

                  if( !newsol )
                     continue;
               }

               nfoundvars++;
               nfoundvarsprob++;

               /* create new variable, compute objective function value and add it to the master constraints and cuts it belongs to */
               SCIP_CALL( createNewMasterVar(scip, sols[j], prob) );
            }
         }

         SCIP_CALL( SCIPsetLongintParam(pricerdata->pricingprobs[prob], "limits/stallnodes", -1) );
         SCIP_CALL( SCIPsetLongintParam(pricerdata->pricingprobs[prob], "limits/nodes", -1) ); 
         SCIP_CALL( SCIPsetRealParam(pricerdata->pricingprobs[prob], "limits/gap", 0.0) ); 
         SCIP_CALL( SCIPsetIntParam(pricerdata->pricingprobs[prob], "limits/bestsol", -1) ); 


         //SCIP_CALL( SCIPstartClock(scip, pricerdata->freeclock) );
         //SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[prob]) );
         //SCIP_CALL( SCIPstopClock(scip, pricerdata->freeclock) );
      }
   }

   /* if no variables were found so far, solve the pricing MIPs to optimality and check whether
    * solutions corresponding to variables with negative reduced costs where found
    */
   if( nfoundvars == 0 )
   {
      bestredcostvalid = ( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL ? TRUE : FALSE );

      for( i = 0; i < pricerdata->npricingprobs && (pricetype == GCG_PRICETYPE_FARKAS ||
            (nfoundvars < pricerdata->maxvarsroundredcost && successfulmips < pricerdata->maxsuccessfulmipsredcost))
              && (nfoundvars == 0 || pricerdata->dualsolconv[permu[i]] > 0 || !pricerdata->onlyposconv)
              && (pricetype == GCG_PRICETYPE_REDCOST || nfoundvars < pricerdata->maxvarsroundfarkas); i++)
      {
         prob = permu[i];

         if( pricerdata->pricingprobs[prob] == NULL )
            continue;

         /* set objective limit, such that only solutions with negative reduced costs are accepted */
         //SCIP_CALL( SCIPsetObjlimit(pricerdata->pricingprobs[prob], pricerdata->dualsolconv[prob]) );

#ifdef DEBUG_PRICING_ALL_OUTPUT
         if( pricetype == GCG_PRICETYPE_REDCOST )
         {
            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "pricingmip_%d_vars.lp", SCIPgetNVars(scip));
            SCIP_CALL( SCIPwriteOrigProblem(pricerdata->pricingprobs[prob], varname, NULL, FALSE) );

            SCIP_CALL( SCIPsetIntParam(pricerdata->pricingprobs[prob], "display/verblevel", SCIP_VERBLEVEL_HIGH) );

            SCIP_CALL( SCIPwriteParams(pricerdata->pricingprobs[prob], "pricing.set", TRUE, TRUE) );
         }
#endif

         /* set time limit */
         SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
         if( !SCIPisInfinity(scip, timelimit) )
         {
            if( timelimit - SCIPgetTotalTime(scip) > 0 )
            {
               SCIP_CALL( SCIPsetRealParam(pricerdata->pricingprobs[prob], "limits/time", 
                     timelimit - SCIPgetTotalTime(scip)) );
            }
            else
            {
               *result = SCIP_DIDNOTRUN;

               /* free transformed pricing MIPs */
               for( j = 0; j < pricerdata->npricingprobs; j++ )
               {
                  if( pricerdata->pricingprobs[j] != NULL )
                  {
                     SCIP_CALL( SCIPstartClock(scip, pricerdata->freeclock) );
                     SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[j]) );
                     SCIP_CALL( SCIPstopClock(scip, pricerdata->freeclock) );
                  }
               }
               SCIPfreeBufferArray(scip, &permu);
               SCIPfreeBufferArray(scip, &tmpconvdualsol);

               SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );
            
               return SCIP_OKAY;
            }
         }

         SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );
         SCIP_CALL( SCIPstartClock(scip, pricerdata->subsolveclock) );

         /* start clock measuring the transformation effort */
         SCIP_CALL( SCIPstartClock(scip, pricerdata->transformclock) );
         if( !pricerdata->useheurpricing )
         {
            SCIP_CALL( SCIPtransformProb(pricerdata->pricingprobs[prob]) );
         }
         SCIP_CALL( SCIPstopClock(scip, pricerdata->transformclock) );

         /* start clock measuring the presolving effort */
         if( pricetype == GCG_PRICETYPE_REDCOST )
            SCIP_CALL( SCIPstartClock(scip, pricerdata->redcostpresolveclock) );
         else
            SCIP_CALL( SCIPstartClock(scip, pricerdata->farkaspresolveclock) );

         /* presolve the pricing submip */
         if( !pricerdata->useheurpricing )
         {
            SCIP_CALL( SCIPpresolve(pricerdata->pricingprobs[prob]) );
         }

         /* stop clock measuring the presolving effort, start clock measuring the solving effort */
         if( pricetype == GCG_PRICETYPE_REDCOST )
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
         if( pricetype == GCG_PRICETYPE_REDCOST )
            SCIP_CALL( SCIPstopClock(scip, pricerdata->redcostsolveclock) );
         else
            SCIP_CALL( SCIPstopClock(scip, pricerdata->farkassolveclock) );

         pricerdata->solvedsubmipsoptimal++;

         SCIP_CALL( SCIPstopClock(scip, pricerdata->subsolveclock) );
         SCIP_CALL( SCIPstartClock(scip, pricerdata->owneffortclock) );

#ifdef DEBUG_PRICING_ALL_OUTPUT
         if( pricetype == GCG_PRICETYPE_REDCOST )
         {
            SCIP_CALL( SCIPsetIntParam(pricerdata->pricingprobs[prob], "display/verblevel", 0) );
            SCIP_CALL( SCIPprintStatistics(pricerdata->pricingprobs[prob], NULL) );
         }
#endif

         /* so far, the pricing problem should be solved to optimality */
         assert( SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_OPTIMAL
            || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_GAPLIMIT
            || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_USERINTERRUPT 
            || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_INFEASIBLE
            || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_TIMELIMIT );
         /** @todo handle userinterrupt: set result pointer (and lowerbound), handle other solution states */

         if( SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_USERINTERRUPT
            || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_TIMELIMIT )
         {
            nsols = 0;
            sols = NULL;
            bestredcostvalid = FALSE;
            if( result != NULL ) 
               *result = SCIP_DIDNOTRUN;
         } 
         else
         {
            nsols = SCIPgetNSols(pricerdata->pricingprobs[prob]);
            sols = SCIPgetSols(pricerdata->pricingprobs[prob]);
            //printf("Pricingprob %d has found %d sols!\n", prob, nsols);
         }

         if( nsols > 0 && SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_OPTIMAL )
         {
            //printf("objfactor[%d] = %g\n", prob, SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sols[0])/pricerdata->objsums[prob]);
            updateObjFactor(scip, prob, SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sols[0]));
         }
         if( nsols > 0 && SCIPisSumNegative(scip, SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sols[0]) - pricerdata->dualsolconv[prob]) )
         {
            bestredcost += GCGrelaxGetNIdenticalBlocks(origprob, prob) * 
               (SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sols[0]) - pricerdata->dualsolconv[prob]);
         }

         if( SCIPgetStatus(pricerdata->pricingprobs[prob]) != SCIP_STATUS_OPTIMAL
            && SCIPgetStatus(pricerdata->pricingprobs[prob]) != SCIP_STATUS_INFEASIBLE )
         {
            bestredcostvalid = FALSE;
         }

         nfoundvarsprob = 0;

         for( j = 0; j < nsols && nfoundvarsprob <= pricerdata->maxsolsprob &&
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
            if( SCIPisSumNegative(scip, SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sols[j]) - pricerdata->dualsolconv[prob]) )
            {
               if( pricerdata->checksols )
               {
                  SCIP_CALL( checkSolNew(scip, sols, j, prob, &newsol) );

                  if( !newsol )
                     continue;
               }

               nfoundvars++;
               nfoundvarsprob++;
               if( nfoundvarsprob == 1 )
                  successfulmips++;

               /* create new variable, compute objective function value and add it to the master constraints and cuts it belongs to */
               SCIP_CALL( createNewMasterVar(scip, sols[j], prob) );
            }
         }

         //SCIP_CALL( SCIPstartClock(scip, pricerdata->freeclock) );
         //SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[prob]) );
         //SCIP_CALL( SCIPstopClock(scip, pricerdata->freeclock) );
      }
   }

   if( bestredcostvalid && i < pricerdata->npricingprobs )
     for( j = i; j < pricerdata->npricingprobs; j++ )
       if( pricerdata->pricingprobs[permu[j]] != NULL )
         {
	   bestredcostvalid = FALSE;
         }
   
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if( pricerdata->pricingprobs[i] != NULL )
      {
         SCIP_CALL( SCIPstartClock(scip, pricerdata->freeclock) );
         SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[i]) );
         SCIP_CALL( SCIPstopClock(scip, pricerdata->freeclock) );
      }
   }

   if( pricetype == GCG_PRICETYPE_REDCOST && bestredcostvalid )
   {
      assert(lowerbound != NULL);
#ifdef DEBUG_PRICING
      printf("lower bound = %g, bestredcost = %g\n", SCIPgetLPObjval(scip) + bestredcost, bestredcost);
#endif
      *lowerbound = SCIPgetLPObjval(scip) + bestredcost;
   }

   SCIPfreeBufferArray(scip, &permu);
   SCIPfreeBufferArray(scip, &tmpconvdualsol);

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
   if( pricerdata != NULL )
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

   pricerdata->currnodenr = -1;

   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);

   /* init array containing all pricing problems */
   pricerdata->npricingprobs = GCGrelaxGetNPricingprobs(origprob);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->pricingprobs), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nvarsprob), pricerdata->npricingprobs) );

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if( GCGrelaxIsPricingprobRelevant(origprob, i) )
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
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->dualsol), nmasterconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->dualsolconv), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->objsums), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->objfactors), pricerdata->npricingprobs) );

   BMSclearMemoryArray(pricerdata->objfactors, pricerdata->npricingprobs);

   /* alloc memory for solution values of variables in pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->solvals), SCIPgetNOrigVars(origprob)) );

   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->subsolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->heursolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->redcostclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->redcostsolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->redcostpresolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->farkasclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->farkassolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->farkaspresolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->owneffortclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->freeclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->transformclock)) );

   pricerdata->solvedsubmipsoptimal = 0;
   pricerdata->solvedsubmipsheur = 0;
   pricerdata->calls = 0;
   pricerdata->redcostcalls = 0;
   pricerdata->farkascalls = 0;
   pricerdata->nallsame = 0;

   /* set variable type for master variables */
   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( discretization )
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
   for( v = 0; v < nvars; v++ )
   {
      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      if( vardata->blocknr == -1 )
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
         for( i = 0; i < vardata->data.origvardata.ncoefs; i++ )
         {
            SCIP_CONS* linkcons;
            assert(!SCIPisZero(scip, vardata->data.origvardata.coefs[i]));
            SCIP_CALL( SCIPgetTransformedCons(scip, vardata->data.origvardata.linkconss[i], 
                  &linkcons) );

            SCIP_CALL( SCIPaddCoefLinear(scip, linkcons, 
                  newvar, vardata->data.origvardata.coefs[i]) );
         }
         SCIPreleaseVar(scip, &newvar);

      }
   }

   SCIP_CALL( SCIPhashmapCreate(&(pricerdata->mapcons2idx), SCIPblkmem(scip), 10 * nmasterconss) );
   for( i = 0; i < nmasterconss; i++ )
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
   SCIPfreeMemoryArray(scip, &(pricerdata->dualsol));
   SCIPfreeMemoryArray(scip, &(pricerdata->dualsolconv));
   SCIPfreeMemoryArray(scip, &(pricerdata->objsums));
   SCIPfreeMemoryArray(scip, &(pricerdata->objfactors));
   SCIPfreeMemoryArray(scip, &(pricerdata->solvals));
   SCIPfreeMemoryArray(scip, &(pricerdata->nvarsprob));

   printf("calls = %d\n", pricerdata->calls);
   printf("time for pricing without sub-MIPs: %f\n", SCIPgetClockTime(scip, pricerdata->owneffortclock));
   printf("solved sub-MIPs heur = %d\n", pricerdata->solvedsubmipsheur);
   printf("solved sub-MIPs optimal = %d\n", pricerdata->solvedsubmipsoptimal);
   printf("time for solving sub-MIPs for pricing to optimality: %f\n", SCIPgetClockTime(scip, pricerdata->subsolveclock));
   printf("time for solving sub-MIPs for pricing heuristically: %f\n", SCIPgetClockTime(scip, pricerdata->heursolveclock));
   printf("farkas calls = %d, redcost calls = %d\n", pricerdata->farkascalls, pricerdata->redcostcalls);
   printf("time for farkas pricing (presolving): %f\n", SCIPgetClockTime(scip, pricerdata->farkaspresolveclock));
   printf("time for farkas pricing (solving): %f\n", SCIPgetClockTime(scip, pricerdata->farkassolveclock));
   printf("time for farkas pricing (total): %f\n", SCIPgetClockTime(scip, pricerdata->farkasclock));
   printf("time for redcost pricing (presolving): %f\n", SCIPgetClockTime(scip, pricerdata->redcostpresolveclock));
   printf("time for redcost pricing (solving): %f\n", SCIPgetClockTime(scip, pricerdata->redcostsolveclock));
   printf("time for redcost pricing (total): %f\n", SCIPgetClockTime(scip, pricerdata->redcostclock));
   printf("time for transformation: %f\n", SCIPgetClockTime(scip, pricerdata->transformclock));
   printf("time for freeing sub-MIPs: %f\n", SCIPgetClockTime(scip, pricerdata->freeclock));
#if 0
   printf("estimated values for all pricing MIPs same: %d times\n", pricerdata->nallsame);
#endif


   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->subsolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->heursolveclock)) );
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

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetNTotalNodes(scip) == pricerdata->currnodenr )
   {
      pricerdata->nroundsredcost++;
   }
   else
   {
      pricerdata->currnodenr = SCIPgetNTotalNodes(scip);
      pricerdata->nroundsredcost = 0;
   }

   if( pricerdata->nroundsredcost >= pricerdata->maxroundsredcost && pricerdata->currnodenr != 1)
   {
      SCIPdebugMessage("pricing aborted at node %lld\n", pricerdata->currnodenr);
      return SCIP_OKAY;
   }

   *result = SCIP_SUCCESS;

   SCIP_CALL( SCIPstartClock(scip, pricerdata->redcostclock) );

   //printf("pricerredcost\n");
   retcode = performPricing(scip, pricer, GCG_PRICETYPE_REDCOST, result, lowerbound);

   SCIP_CALL( SCIPstopClock(scip, pricerdata->redcostclock) );

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
   retcode = performPricing(scip, pricer, GCG_PRICETYPE_FARKAS, NULL, NULL);

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

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxsuccessfulmipsredcost",
         "maximal number of pricing mips leading to new variables solved solved in one redcost pricing round",
         &pricerdata->maxsuccessfulmipsredcost, FALSE, DEFAULT_MAXSUCCESSFULMIPSREDCOST, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxvarsroundredcost",
         "maximal number of variables created in one redcost pricing round",
         &pricerdata->maxvarsroundredcost, FALSE, DEFAULT_MAXVARSROUNDREDCOST, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxvarsroundfarkas",
         "maximal number of variables created in one farkas pricing round",
         &pricerdata->maxvarsroundfarkas, FALSE, DEFAULT_MAXVARSROUNDFARKAS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxroundsredcost",
         "maximal number of pricing rounds per node after the root node",
         &pricerdata->maxroundsredcost, FALSE, DEFAULT_MAXROUNDSREDCOST, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxsolsprob",
         "maximal number of variables added for each block in a pricinground",
         &pricerdata->maxsolsprob, FALSE, DEFAULT_MAXSOLSPROB, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/checksols",
         "should solutions of the pricing MIPs be checked for duplicity?",
         &pricerdata->checksols, TRUE, DEFAULT_CHECKSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/useheurpricing",
         "should pricing be performed heuristically befor solving the MIPs to optimality?",
         &pricerdata->useheurpricing, TRUE, DEFAULT_USEHEURPRICING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/onlyposconv",
         "should only pricing problems be solved with a positive dualsol of the convexity constraint, if possible?",
         &pricerdata->onlyposconv, TRUE, DEFAULT_ONLYPOSCONV, NULL, NULL) );

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
