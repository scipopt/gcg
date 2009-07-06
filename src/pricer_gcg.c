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

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "pricer_gcg.h"
#include "probdata_gcg.h"
#include "scip/cons_linear.h"
#include "scip/scip.h"
#include <stdio.h>
#include <stdlib.h>


#define PRICER_NAME            "gcg"
#define PRICER_DESC            "pricer for gcg"
#define PRICER_PRIORITY        5000000
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

#define EPS 0.0001


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
   SCIP_CLOCK* initclock;
   int solvedsubmips;
   int calls;
   int farkascalls;
   int redcostcalls;
   int initcalls;
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


/* informs an original variable, that a variable in the master problem was created, 
 * that contains a part of the original variable.
 * Saves this information in the original variable's data */
SCIP_RETCODE GCGpricerAddMasterVarToOrigVar(
   SCIP_PRICER*          pricer,
   SCIP_VAR*             origvar,
   SCIP_VAR*             var,
   SCIP_Real             val
   )
{
   SCIP_VARDATA* vardata;
   SCIP_PRICERDATA* pricerdata;

   assert(pricer != NULL);
   assert(origvar != NULL);
   assert(var != NULL);

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

   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   if ( pricetype == GCG_PRICETYPE_REDCOST )
      pricerdata->redcostcalls++;
   if ( pricetype == GCG_PRICETYPE_FARKAS )
      pricerdata->farkascalls++;

   SCIP_CALL( SCIPstartClock(scip, pricerdata->owneffortclock) );
   pricerdata->calls++;
   nfoundvars = 0;

   /* set objective value of all variables in the pricing problems to 0 (for farkas pricing) /
    * to the original objective of the variable (for redcost pricing) */
   for ( i = 0; i < pricerdata->npricingprobs; i++)
   {
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
            assert(vardata->data.pricingvardata.origvar != NULL);
            SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], 
                  SCIPvarGetObj(vardata->data.pricingvardata.origvar)) );
         }  
      }
   }

   /* get the constraints of the master problem and the corresponding constraints in the original problem */
   nmasterconss = GCGrelaxGetNMasterConss(pricerdata->origprob);
   masterconss = GCGrelaxGetMasterConss(pricerdata->origprob);
   origconss = GCGrelaxGetLinearOrigMasterConss(pricerdata->origprob);
   

   /* compute reduced cost and update objectives in the pricing problems */
   for ( i = 0; i < nmasterconss; i++ )
   {
      /* farkas pricing */
      if ( pricetype == GCG_PRICETYPE_FARKAS )
      {
         assert(SCIPconsIsTransformed(masterconss[i]));
         pricerdata->redcost[i] = SCIPgetDualfarkasLinear(scip, masterconss[i]);
         //printf("farkas value of cons %s = %f\n", SCIPconsGetName(masterconss[i]), pricerdata->redcost[i]);
      }
      /* redcost pricing */
      if ( pricetype == GCG_PRICETYPE_REDCOST )
      {      
         pricerdata->redcost[i] = SCIPgetDualsolLinear(scip, masterconss[i]);
      }
      //printf("master[%d] = %f\n", i, pricerdata->redcost[i]);
      if ( !SCIPisFeasZero(scip, pricerdata->redcost[i]) )
      {
         /* for all variables in the constraint, modify the objective of the corresponding variable in a pricing problem */
         consvars = SCIPgetVarsLinear(pricerdata->origprob, origconss[i]);
         nconsvars = SCIPgetNVarsLinear(pricerdata->origprob, origconss[i]);
         consvals = SCIPgetValsLinear(pricerdata->origprob, origconss[i]);
         for ( j = 0; j < nconsvars; j++ )
         {
            vardata = SCIPvarGetData(consvars[j]);
            assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
            if ( vardata->blocknr != -1 )
            {
               assert(vardata->data.origvardata.pricingvar != NULL);
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[vardata->blocknr], 
                     vardata->data.origvardata.pricingvar, -1.0 * pricerdata->redcost[i] * consvals[j]) );
            }
         }
      }
   }
   SCIP_CALL( SCIPwriteLP(scip, "test.lp") );

   /* get dual solutions / farkas values of the convexity constraints */
   for ( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if ( pricetype == GCG_PRICETYPE_FARKAS || pricetype == GCG_PRICETYPE_INIT )
      {
         assert(SCIPconsIsTransformed(GCGrelaxGetConvCons(pricerdata->origprob, i)));
         pricerdata->redcostconv[i] = SCIPgetDualfarkasLinear(scip, GCGrelaxGetConvCons(pricerdata->origprob, i));
         //printf("farkas value of cons %s = %f\n", SCIPconsGetName(GCGrelaxGetConvCons(pricerdata->origprob, i)), pricerdata->redcostconv[i]);
      }
      else
      {
         pricerdata->redcostconv[i] = SCIPgetDualsolLinear(scip, GCGrelaxGetConvCons(pricerdata->origprob, i));
      }
   }

   /* solve the pricing MIPs and chech whether solutions corresponding to variables with negative reduced costs where found */
   for ( i = 1; i <= pricerdata->npricingprobs && (nfoundvars == 0 || pricetype == GCG_PRICETYPE_REDCOST); i++)
   {
      prob = (pricerdata->pricingprobnr + i) % pricerdata->npricingprobs;

      SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );
      SCIP_CALL( SCIPstartClock(scip, pricerdata->subsolveclock) );

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

      /* start clock measuring the solving effort */
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
         || SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_USERINTERRUPT );
      /** @todo handle userinterrupt: set result pointer (and lowerbound), handle other solution states */

      nsols = SCIPgetNSols(pricerdata->pricingprobs[prob]);
      sols = SCIPgetSols(pricerdata->pricingprobs[prob]);
      for ( j = 0; j < nsols; j++ )
      {
         //SCIP_CALL( SCIPprintSol(pricerdata->pricingprobs[prob], NULL, NULL, FALSE) );
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
            newvardata->data.mastervardata.norigvars = nprobvars;
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvars), nprobvars) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvals), nprobvars) );

            /* compute objective coefficient of the variable */
            objcoeff = 0;
            for ( k = 0; k < nprobvars; k++ )
            {
               if ( !SCIPisFeasZero(scip, pricerdata->solvals[k]) )
               {
                  vardata = SCIPvarGetData(probvars[k]);
                  assert(vardata->vartype == GCG_VARTYPE_PRICING);
                  assert(vardata->data.pricingvardata.origvar != NULL);
                  /* add quota of original variable's objcoef to the master variable's coef */
                  objcoeff += pricerdata->solvals[k] * SCIPvarGetObj(vardata->data.pricingvardata.origvar);
               }
            }

            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "p_%d_%d", prob, pricerdata->nvarsprob[prob]);
            pricerdata->nvarsprob[prob]++;

            /* create variable in the master problem */
            SCIP_CALL( SCIPcreateVar(scip, &newvar, varname, 
                  0, 1, objcoeff, SCIP_VARTYPE_CONTINUOUS, 
                  TRUE, TRUE, NULL, NULL, gcgvardeltrans, newvardata) );

            /* update variable datas */
            for ( k = 0; k < nprobvars; k++ )
            {
               if ( !SCIPisFeasZero(scip, pricerdata->solvals[k]) )
               {
                  vardata = SCIPvarGetData(probvars[k]);
                  assert(vardata->vartype == GCG_VARTYPE_PRICING);
                  assert(vardata->data.pricingvardata.origvar != NULL);
                  /* save in the master problem variable's data the quota of the corresponding original variable */
                  newvardata->data.mastervardata.origvars[k] = vardata->data.pricingvardata.origvar;
                  newvardata->data.mastervardata.origvals[k] = pricerdata->solvals[k];
                  /* save the quota in the original variable's data */
                  SCIP_CALL( GCGpricerAddMasterVarToOrigVar(pricer, vardata->data.pricingvardata.origvar, newvar, pricerdata->solvals[k]) );
               }
               else
               {
                  vardata = SCIPvarGetData(probvars[k]);
                  assert(vardata->vartype == GCG_VARTYPE_PRICING);
                  assert(vardata->data.pricingvardata.origvar != NULL);
                  /** @todo really store variable not connected to the master variable? */
                  newvardata->data.mastervardata.origvars[k] = vardata->data.pricingvardata.origvar;
                  newvardata->data.mastervardata.origvals[k] = 0.0;
               }
            }

            SCIP_CALL( SCIPaddPricedVar(scip, newvar, 1.0) );

            SCIPchgVarUbLazy(scip, newvar, 1.0);

            for ( k = 0; k < nmasterconss; k++ )
            {
               conscoeff = 0;
               for ( l = 0; l < nprobvars; l++ )
               {
                  if ( !SCIPisFeasZero(scip, pricerdata->solvals[l]) )
                  {
                     //printf("var = %s\n", SCIPvarPrintName(probvars[l]));
                     vardata = SCIPvarGetData(probvars[l]);
                     assert(vardata != NULL);
                     assert(vardata->vartype == GCG_VARTYPE_PRICING);
                     vardata = SCIPvarGetData(vardata->data.pricingvardata.origvar);
                     assert(vardata != NULL);
                     assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
                     assert(vardata->data.origvardata.coefs != NULL);
                                          
                     conscoeff += vardata->data.origvardata.coefs[k] * pricerdata->solvals[l];
                  }
               }
               SCIP_CALL( SCIPaddCoefLinear(scip, masterconss[k], newvar, conscoeff) );
               //printf("new variable has coef = %f in constraint %d\n", conscoeff, k);
            }
            SCIP_CALL( SCIPaddCoefLinear(scip, GCGrelaxGetConvCons(pricerdata->origprob, prob), newvar, 1) );
            //printf("new variable has coef = %f in convexity constraint %d\n", 1.0, prob);

            SCIPreleaseVar(scip, &newvar);

            /* ??? */
            if ( pricetype == GCG_PRICETYPE_FARKAS && FALSE )
            {
               SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[prob]) );
               
               SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );

               //printf("Farkas-Pricing: found %d new vars\n", nfoundvars);

               return SCIP_OKAY;
            }

         }
      }

      SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[prob]) );
   }

   SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );

   //printf("Pricing: found %d new vars\n", nfoundvars);

   return SCIP_OKAY;
}

#if 0
static
SCIP_RETCODE createInitialVars(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_PRICER*          pricer              /* the pricer */
   )
{
   SCIP_PRICERDATA* pricerdata;            /* the data of the pricer */
   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_CONS** origconss;
   int norigconss;
   char varname[SCIP_MAXSTRLEN];

   int i;
   int j;
   int k;
   int l;
   int prob;

   int nfoundvars;
   
   SCIP_Real objsum;

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

   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   pricerdata->initcalls++;

   SCIP_CALL( SCIPstartClock(scip, pricerdata->initclock) );

   nfoundvars = 0;

   /* set objective value of all variables in the pricing problems to 0 */
   for ( i = 0; i < pricerdata->npricingprobs; i++)
   {

   }

   /* get the constraints of the master problem */
   GCGrelaxGetMasterConss(scip, &masterconss, &nmasterconss);

   /* get the corresponding constraints in the original problem */
   GCGrelaxGetLinearOrigMasterConss(scip, &origconss, &norigconss);
   
   assert(nmasterconss == norigconss);

   /* compute reduced cost and set objectives in the pricing problems */
      /* farkas pricing */
   for ( i = 0; i < nmasterconss; i++ )
   {
      pricerdata->redcost[i] = 2.0;
   }

   for ( prob = 0; prob < pricerdata->npricingprobs; prob++)
   {

      probvars = SCIPgetVars(pricerdata->pricingprobs[prob]);
      nprobvars = SCIPgetNVars(pricerdata->pricingprobs[prob]);
      
      /* compute new objective coeffs */
      for ( k = 0; k < nprobvars; k++ )
      {
         SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[prob], probvars[k], 0) );
      }
      for ( k = 0; k < nmasterconss; k++ )
      {
         consvars = SCIPgetVarsLinear(pricerdata->origprob, origconss[k]);
         nconsvars = SCIPgetNVarsLinear(pricerdata->origprob, origconss[k]);
         consvals = SCIPgetValsLinear(pricerdata->origprob, origconss[k]);
         for ( l = 0; l < nconsvars; l++ )
         {
            vardata = SCIPvarGetData(consvars[l]);
            assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
            assert(vardata->blocknr != -1);
            assert(vardata->data.origvardata.pricingvar != NULL);
            if ( vardata->blocknr == prob )
            {
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[prob], 
                     vardata->data.origvardata.pricingvar, -1.0 * pricerdata->redcost[k] * consvals[l]) );
            }
         }
      }
      
      objsum = 0.0;
      for ( j = 0; j < nprobvars; j++ )
      {
         objsum += SCIPvarGetObj(probvars[j]);
      }
      while ( objsum < 0 )
      {
         SCIP_CALL( SCIPpresolve(pricerdata->pricingprobs[prob]) );
         SCIP_CALL( SCIPsolve(pricerdata->pricingprobs[prob]) );
         
         assert( SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_OPTIMAL);

         nsols = SCIPgetNSols(pricerdata->pricingprobs[prob]);
         sols = SCIPgetSols(pricerdata->pricingprobs[prob]);
         for ( j = 0; j < nsols; j++ )
         {
            nfoundvars++;
            /* get variables of the pricing problem and their values in the current solution */
            SCIP_CALL( SCIPgetSolVals(pricerdata->pricingprobs[prob], sols[j], nprobvars, probvars, pricerdata->solvals) );
         
            assert(scip != NULL);
         
            /* create variable in the master problem */
            SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
            newvardata->vartype = GCG_VARTYPE_MASTER;
            newvardata->blocknr = prob;
            newvardata->data.mastervardata.norigvars = nprobvars;
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvars), nprobvars) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvals), nprobvars) );
            
            /* compute objective coefficient */
            objcoeff = 0;
            for ( k = 0; k < nprobvars; k++ )
            {
               if ( !SCIPisFeasZero(scip, pricerdata->solvals[k]) )
               {
                  vardata = SCIPvarGetData(probvars[k]);
                  assert(vardata->vartype == GCG_VARTYPE_PRICING);
                  assert(vardata->data.pricingvardata.origvar != NULL);
                  objcoeff += pricerdata->solvals[k] * SCIPvarGetObj(vardata->data.pricingvardata.origvar);
                  newvardata->data.mastervardata.origvars[k] = vardata->data.pricingvardata.origvar;
                  newvardata->data.mastervardata.origvals[k] = pricerdata->solvals[k];
               }
            }
         
         
            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "p_%d_%d", prob, pricerdata->nvarsprob[prob]);
            pricerdata->nvarsprob[prob]++;
 
            SCIP_CALL( SCIPcreateVar(scip, &newvar, varname, 
                  0, 1, objcoeff, SCIP_VARTYPE_CONTINUOUS, 
                  TRUE, TRUE, NULL, NULL, gcgvardeltrans, newvardata) );

            SCIPchgVarUbLazy(scip, newvar, 1.0);

            SCIP_CALL( SCIPaddVar(scip, newvar) );
         
            for ( k = 0; k < norigconss; k++ )
            {
               conscoeff = 0;
               for ( l = 0; l < nprobvars; l++ )
               {
                  if ( !SCIPisFeasZero(scip, pricerdata->solvals[l]) )
                  {
                     vardata = SCIPvarGetData(probvars[l]);
                     assert(vardata != NULL);
                     assert(vardata->vartype == GCG_VARTYPE_PRICING);
                     vardata = SCIPvarGetData(vardata->data.pricingvardata.origvar);
                     assert(vardata != NULL);
                     assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
                     assert(vardata->data.origvardata.coefs != NULL);
                     
                     conscoeff += vardata->data.origvardata.coefs[k] * pricerdata->solvals[l];

                     /* reduce the redcost */
                     pricerdata->redcost[k] -= vardata->data.origvardata.coefs[k] * pricerdata->solvals[l];
                  }
               }
               SCIP_CALL( SCIPaddCoefLinear(scip, masterconss[k], newvar, conscoeff) );
               //printf("new variable has coef = %f in constraint %d\n", conscoeff, k);
            }
            SCIP_CALL( SCIPaddCoefLinear(scip, GCGrelaxGetConvCons(scip, prob), newvar, 1) );
            //printf("new variable has coef = %f in convexity constraint %d\n", 1.0, prob);

            SCIPreleaseVar(scip, &newvar);
         }

         SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[prob]) );

         /* compute new objective coeffs */
         for ( k = 0; k < nprobvars; k++ )
         {
            SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[prob], probvars[k], 0) );
         }
         for ( k = 0; k < nmasterconss; k++ )
         {
            consvars = SCIPgetVarsLinear(pricerdata->origprob, origconss[k]);
            nconsvars = SCIPgetNVarsLinear(pricerdata->origprob, origconss[k]);
            consvals = SCIPgetValsLinear(pricerdata->origprob, origconss[k]);
            for ( l = 0; l < nconsvars; l++ )
            {
               vardata = SCIPvarGetData(consvars[l]);
               assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
               assert(vardata->blocknr != -1);
               assert(vardata->data.origvardata.pricingvar != NULL);
               if ( vardata->blocknr == prob )
               {
                  /* modify the objective of the corresponding variable in the pricing problem */
                  SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[vardata->blocknr], 
                        vardata->data.origvardata.pricingvar, -1.0 * pricerdata->redcost[i] * consvals[j]) );
               }
            }
         }
         objsum = 0.0;
         for ( j = 0; j < nprobvars; j++ )
         {
            objsum += SCIPvarGetObj(probvars[j]);
         }
         
      }

   }

   SCIP_CALL( SCIPstopClock(scip, pricerdata->initclock) );

   printf("createInitialVars: varscreated = %d\n", nfoundvars);

   return SCIP_OKAY;
}
#endif




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

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* set last used pricingprobnr to 0 */
   pricerdata->pricingprobnr = 0;

   /* init array containing all pricing problems */
   pricerdata->npricingprobs = GCGrelaxGetNPricingprobs(pricerdata->origprob);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->pricingprobs), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nvarsprob), pricerdata->npricingprobs) );

   for ( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      pricerdata->pricingprobs[i] = GCGrelaxGetPricingprob(pricerdata->origprob, i);
      pricerdata->nvarsprob[i] = 0;
   }
   
   /* alloc memory for arrays of reduced cost */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcost), GCGrelaxGetNMasterConss(pricerdata->origprob)) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcostconv), GCGrelaxGetNPricingprobs(pricerdata->origprob)) );

   /* alloc memory for solution values of variables in pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->solvals), SCIPgetNOrigVars(pricerdata->origprob)) );

   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->subsolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->redcostclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->redcostsolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->redcostpresolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->farkasclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->farkassolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->farkaspresolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->owneffortclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->initclock)) );

   pricerdata->solvedsubmips = 0;
   pricerdata->calls = 0;
   pricerdata->redcostcalls = 0;
   pricerdata->farkascalls = 0;
   pricerdata->initcalls = 0;

   //SCIP_CALL( createInitialVars(scip, pricer) );

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
   
   SCIPfreeMemoryArray(scip, &(pricerdata->pricingprobs));
   SCIPfreeMemoryArray(scip, &(pricerdata->redcost));
   SCIPfreeMemoryArray(scip, &(pricerdata->redcostconv));
   SCIPfreeMemoryArray(scip, &(pricerdata->solvals));
   SCIPfreeMemoryArray(scip, &(pricerdata->nvarsprob));

   printf("calls = %d\n", pricerdata->calls);
   printf("time for pricing without sub-MIPs: %f\n", SCIPgetClockTime(scip, pricerdata->owneffortclock));
   printf("solved sub-MIPs = %d\n", pricerdata->solvedsubmips);
   printf("time for solving sub-MIPs for pricing: %f\n", SCIPgetClockTime(scip, pricerdata->subsolveclock));
   printf("time for initial var creation: %f\n", SCIPgetClockTime(scip, pricerdata->initclock));
   printf("init calls = %d, farkas calls = %d, redcost calls = %d\n", pricerdata->initcalls, pricerdata->farkascalls, pricerdata->redcostcalls);
   printf("time for farkas pricing (presolving): %f\n", SCIPgetClockTime(scip, pricerdata->farkaspresolveclock));
   printf("time for farkas pricing (solving): %f\n", SCIPgetClockTime(scip, pricerdata->farkassolveclock));
   printf("time for farkas pricing (total): %f\n", SCIPgetClockTime(scip, pricerdata->farkasclock));
   printf("time for redcost pricing (presolving): %f\n", SCIPgetClockTime(scip, pricerdata->redcostpresolveclock));
   printf("time for redcost pricing (solving): %f\n", SCIPgetClockTime(scip, pricerdata->redcostsolveclock));
   printf("time for redcost pricing (total): %f\n", SCIPgetClockTime(scip, pricerdata->redcostclock));



   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->subsolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->redcostclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->redcostsolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->redcostpresolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->farkasclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->farkassolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->farkaspresolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->owneffortclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->initclock)) );

   
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

   return SCIP_OKAY;
}
