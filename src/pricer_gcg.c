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
#include <stdio.h>
#include <stdlib.h>


#define PRICER_NAME            "gcg"
#define PRICER_DESC            "pricer for gcg"
#define PRICER_PRIORITY        5000000
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

#define VARNAMELEN 20

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
   SCIP_CLOCK* owneffortclock;
   int solvedsubmips;
   int calls;
   int farkascalls;
   int redcostcalls;
};


/*
 * Vardata methods
 */

/** problem reading method of reader */
static
SCIP_DECL_VARDELORIG(gcgvardeltrans)
{  
   if ( (*vardata)->vartype == GCG_VARTYPE_MASTER )
   {
      SCIPfreeBlockMemoryArray(scip, &((*vardata)->data.mastervardata.vals), (*vardata)->data.mastervardata.norigvars);
      SCIPfreeBlockMemoryArray(scip, &((*vardata)->data.mastervardata.origvars), (*vardata)->data.mastervardata.norigvars);
   }
   
   SCIPfreeBlockMemory(scip, vardata);

   return SCIP_OKAY;
}


/*
 * Local methods
 */

static
SCIP_RETCODE performPricing(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_PRICER*          pricer,             /* the pricer */
   SCIP_Bool             redcostpricing      /* true, if redcost pricing should be performed, 
                                                false, for farkas pricing */
   )
{
   SCIP_PRICERDATA* pricerdata;            /* the data of the pricer */
   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_CONS** origconss;
   int norigconss;
   char* varname;

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

   if ( redcostpricing == TRUE )
      pricerdata->redcostcalls++;
   else
      pricerdata->farkascalls++;

   SCIP_CALL( SCIPstartClock(scip, pricerdata->owneffortclock) );
   pricerdata->calls++;
   nfoundvars = 0;

   SCIP_CALL(SCIPallocBufferArray(scip, &varname, VARNAMELEN) );

   /* set objective value of all variables in the pricing problems to 0 */
   for ( i = 0; i < pricerdata->npricingprobs; i++)
   {
      probvars = SCIPgetVars(pricerdata->pricingprobs[i]);
      nprobvars = SCIPgetNVars(pricerdata->pricingprobs[i]);

      for ( j = 0; j < nprobvars; j++ )
      {
         if ( redcostpricing == FALSE )
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

   /* get the constraints of the master problem */
   GCGprobGetMasterConss(scip, &masterconss, &nmasterconss);

   /* get the corresponding constraints in the original problem */
   GCGprobGetOrigMasterConss(scip, &origconss, &norigconss);
   

   /* compute reduced cost and set objectives in the pricing problems */
      /* farkas pricing */
   for ( i = 0; i < nmasterconss; i++ )
   {
      if ( redcostpricing == FALSE )
      {
         pricerdata->redcost[i] = SCIPgetDualfarkasLinear(scip, masterconss[i]);
      }
      else
      {      
         pricerdata->redcost[i] = SCIPgetDualsolLinear(scip, masterconss[i]);
      }
      //printf("master[%d] = %f\n", i, pricerdata->redcost[i]);
      if ( !SCIPisFeasZero(scip, pricerdata->redcost[i]) )
      {
         consvars = SCIPgetVarsLinear(pricerdata->origprob, origconss[i]);
         nconsvars = SCIPgetNVarsLinear(pricerdata->origprob, origconss[i]);
         consvals = SCIPgetValsLinear(pricerdata->origprob, origconss[i]);
         for ( j = 0; j < nconsvars; j++ )
         {
            vardata = SCIPvarGetData(consvars[j]);
            assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
            assert(vardata->blocknr != -1);
            assert(vardata->data.origvardata.pricingvar != NULL);
            /* modify the objective of the corresponding variable in the pricing problem */
            SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[vardata->blocknr], 
                  vardata->data.origvardata.pricingvar, -1.0 * pricerdata->redcost[i] * consvals[j]) );
         }
      }
   }
   //printf(".\n");
   for ( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if ( redcostpricing == FALSE )
      {
         pricerdata->redcostconv[i] = SCIPgetDualfarkasLinear(scip, GCGprobGetConvCons(scip, i));
      }
      else
      {
         pricerdata->redcostconv[i] = SCIPgetDualsolLinear(scip, GCGprobGetConvCons(scip, i));
         if ( pricerdata->redcostconv[i] > 0 )
         {
            //printf("redcostconv[%d] = %f > 0\n", i, pricerdata->redcostconv[i]);
         }
      }
      //printf("conv[%d] = %f\n", i, pricerdata->redcostconv[i]);
   }


   for ( i = 1; i <= pricerdata->npricingprobs && (nfoundvars == 0 || redcostpricing == TRUE); i++)
   {
      prob = (pricerdata->pricingprobnr + i) % pricerdata->npricingprobs;

      SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );
      SCIP_CALL( SCIPstartClock(scip, pricerdata->subsolveclock) );

      //SCIP_CALL( SCIPprintOrigProblem(pricerdata->pricingprobs[prob], NULL, NULL, TRUE) );
      SCIP_CALL( SCIPpresolve(pricerdata->pricingprobs[prob]) );
      SCIP_CALL( SCIPsolve(pricerdata->pricingprobs[prob]) );
      //printf("presolving = %fs, solving = %fs\n", SCIPgetPresolvingTime(pricerdata->pricingprobs[prob]), SCIPgetSolvingTime(pricerdata->pricingprobs[prob]));

      pricerdata->solvedsubmips++;

      SCIP_CALL( SCIPstopClock(scip, pricerdata->subsolveclock) );
      SCIP_CALL( SCIPstartClock(scip, pricerdata->owneffortclock) );

      assert( SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_OPTIMAL);
      nsols = SCIPgetNSols(pricerdata->pricingprobs[prob]);
      sols = SCIPgetSols(pricerdata->pricingprobs[prob]);
      for ( j = 0; j < nsols; j++ )
      {
         /* solution value - dual value of associated convexity constraint < 0 
            --> can make the LP feasible or improve the current solution */ 
         if ( SCIPisFeasNegative(scip, SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sols[j]) - pricerdata->redcostconv[prob]) )
         {

            nfoundvars++;
            /* get variables of the pricing problem and their values in the current solution */
            probvars = SCIPgetOrigVars(pricerdata->pricingprobs[prob]);
            nprobvars = SCIPgetNOrigVars(pricerdata->pricingprobs[prob]);
            SCIP_CALL( SCIPgetSolVals(pricerdata->pricingprobs[prob], sols[j], nprobvars, probvars, pricerdata->solvals) );
            
            assert(scip != NULL);

            /* create variable in the master problem */
            SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
            newvardata->vartype = GCG_VARTYPE_MASTER;
            newvardata->data.mastervardata.norigvars = nprobvars;
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvars), nprobvars) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.vals), nprobvars) );

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
                  newvardata->data.mastervardata.vals[k] = pricerdata->solvals[k];
               }
            }


            (void) SCIPsnprintf(varname, VARNAMELEN, "p_%d_%d", prob, pricerdata->nvarsprob[prob]);
            pricerdata->nvarsprob[prob]++;

            SCIP_CALL( SCIPcreateVar(scip, &newvar, varname, 
                  0, SCIPinfinity(scip), objcoeff, SCIP_VARTYPE_CONTINUOUS, 
                  TRUE, TRUE, NULL, NULL, gcgvardeltrans, newvardata) );
            SCIP_CALL( SCIPaddPricedVar(scip, newvar, 1.0) );

            for ( k = 0; k < norigconss; k++ )
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
            SCIP_CALL( SCIPaddCoefLinear(scip, GCGprobGetConvCons(scip, prob), newvar, 1) );
            //printf("new variable has coef = %f in convexity constraint %d\n", 1.0, prob);

            SCIPreleaseVar(scip, &newvar);

            if ( redcostpricing == FALSE && FALSE )
            {
               SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[prob]) );
               
               SCIPfreeBufferArray(scip, &varname);
               
               SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );

               return SCIP_OKAY;
            }

         }
      }

      SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[prob]) );
   }

   SCIPfreeBufferArray(scip, &varname);

   SCIP_CALL( SCIPstopClock(scip, pricerdata->owneffortclock) );

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

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* set last used pricingprobnr to 0 */
   pricerdata->pricingprobnr = 0;

   /* get rhe original program */
   pricerdata->origprob = GCGprobGetOrigprob(scip);
   
   /* init array containing all pricing problems */
   pricerdata->npricingprobs = GCGprobGetNPricingprobs(scip);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->pricingprobs), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nvarsprob), pricerdata->npricingprobs) );

   for ( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      pricerdata->pricingprobs[i] = GCGprobGetPricingprob(scip, i);
      pricerdata->nvarsprob[i] = 0;
   }
   
   /* alloc memory for arrays of reduced cost */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcost), GCGprobGetNMasterConss(scip)) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcostconv), GCGprobGetNPricingprobs(scip)) );

   /* alloc memory for solution values of variables in pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->solvals), SCIPgetNOrigVars(pricerdata->origprob)) );

   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->subsolveclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->owneffortclock)) );

   pricerdata->solvedsubmips = 0;
   pricerdata->calls = 0;
   pricerdata->redcostcalls = 0;
   pricerdata->farkascalls = 0;

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
   printf("farkas calls = %d, redcost calls = %d\n", pricerdata->farkascalls, pricerdata->redcostcalls);


   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->subsolveclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->owneffortclock)) );

   
   return SCIP_OKAY;
}




/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostGcg)
{  
   assert(scip != NULL);
   assert(pricer != NULL);

   //printf("pricerredcost\n");
   return performPricing(scip, pricer, TRUE);

   return SCIP_OKAY;
}




static
SCIP_DECL_PRICERFARKAS(pricerFarkasGcg)
{  
   assert(scip != NULL);
   assert(pricer != NULL);

   //printf("pricerfarkas\n");
   return performPricing(scip, pricer, FALSE);
}

/* define not used callbacks as NULL */
#define pricerInitGcg NULL
#define pricerExitGcg NULL


/*
 * variable pricer specific interface methods
 */

/** creates the gcg variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerGcg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;

   SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricer(scip, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerFreeGcg, pricerInitGcg, pricerExitGcg, 
         pricerInitsolGcg, pricerExitsolGcg, pricerRedcostGcg, pricerFarkasGcg,
         pricerdata) );

   return SCIP_OKAY;
}
