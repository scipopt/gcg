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

/*
 * Data structures
 */


/** variable pricer data */
struct SCIP_PricerData
{
   int pricingprobnr;
};


/*
 * Vardata methods
 */

/** problem reading method of reader */
static
SCIP_DECL_VARDELORIG(gcgdelvarorig)
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

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   
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

   return SCIP_OKAY;
}




/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostGcg)
{  
   SCIP_PRICERDATA* pricerdata;            /* the data of the pricer */


   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   printf("pricerredcost\n");

   return SCIP_OKAY;
}


static
SCIP_DECL_PRICERFARKAS(pricerFarkasGcg)
{  
   SCIP_PRICERDATA* pricerdata;            /* the data of the pricer */
   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_CONS** origconss;
   int norigconss;
   SCIP* origprob;

   SCIP_Real* farkas;
   int i;
   int j;
   int k;
   int l;

   SCIP_VAR** consvars;
   int nconsvars;
   SCIP_Real* consvals;
   SCIP_VAR** probvars;
   int nprobvars;
   SCIP_Real* solvals;

   SCIP_VARDATA* vardata;
   SCIP_VARDATA* newvardata;
   SCIP_VAR* pricingvar;

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

   /* get original problem */
   origprob = GCGprobGetOrigprob(scip);
   
   /* set objective value of all variables in the pricing problems to 0 */
   for ( i = 0; i < GCGprobGetNPricingprobs(scip); i++)
   {
      probvars = SCIPgetVars(GCGprobGetPricingprob(scip, i));
      nprobvars = SCIPgetNVars(GCGprobGetPricingprob(scip, i));
      printf("nprobvars = %d\n", nprobvars);

      for ( j = 0; j < nprobvars; j++ )
      {
         SCIP_CALL( SCIPchgVarObj(GCGprobGetPricingprob(scip, i), probvars[j], 0) );
      }
   }

   //SCIP_CALL( SCIPprintOrigProblem(GCGprobGetOrigprob(scip), NULL, NULL, TRUE) );

   /* get the constraints of the master problem */
   GCGprobGetMasterConss(scip, &masterconss, &nmasterconss);
   SCIP_CALL( SCIPallocBufferArray(scip, &farkas, nmasterconss) );

   /* get the corresponding constraints of the master problem */
   GCGprobGetOrigMasterConss(scip, &origconss, &norigconss);
   
   for ( i = 0; i < nmasterconss; i++ )
   {
      farkas[i] = SCIPgetDualfarkasLinear(scip, masterconss[i]);
      printf("%f\n", farkas[i]);
      if ( !SCIPisFeasZero(scip, farkas[i]) )
      {
         consvars = SCIPgetVarsLinear(origprob, origconss[i]);
         nconsvars = SCIPgetNVarsLinear(origprob, origconss[i]);
         consvals = SCIPgetValsLinear(origprob, origconss[i]);
         for ( j = 0; j < nconsvars; j++ )
         {
            vardata = SCIPvarGetData(consvars[j]);
            assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
            assert(vardata->blocknr != -1);
            assert(vardata->data.origvardata.pricingvar != NULL);
            /* modify the objective of the corresponding variable in the pricing problem */
            SCIP_CALL( SCIPaddVarObj(GCGprobGetPricingprob(scip, vardata->blocknr), 
                  vardata->data.origvardata.pricingvar, -1.0 * farkas[i] * consvals[j]) );
         }
      }
   }

   /* set objective value of all variables in the pricing problems to 0 */
   for ( i = 0; i < GCGprobGetNPricingprobs(scip); i++)
   {
      probvars = SCIPgetVars(GCGprobGetPricingprob(scip, i));
      nprobvars = SCIPgetNVars(GCGprobGetPricingprob(scip, i));
      printf("nprobvars = %d\n", nprobvars);

      //printf("Pricingprob no. %d: nvars = %d\n", i, nprobvars);

      for ( j = 0; j < nprobvars; j++ )
      {
         //printf( "var = %s, obj = %f\n", SCIPvarGetName(probvars[j]), SCIPvarGetObj(probvars[j]));
      }
   }

   for ( i = 0; i < GCGprobGetNPricingprobs(scip); i++)
   {
      SCIP_CALL( SCIPprintOrigProblem(GCGprobGetPricingprob(scip, i), NULL, NULL, TRUE) );
      SCIP_CALL( SCIPpresolve(GCGprobGetPricingprob(scip, i)) );

      SCIP_CALL( SCIPsolve(GCGprobGetPricingprob(scip, i)) );
      nsols = SCIPgetNSols(GCGprobGetPricingprob(scip, i));
      sols = SCIPgetSols(GCGprobGetPricingprob(scip, i));
      for ( j = 0; j < nsols; j++ )
      {
         /* solution value < 0 --> can make the LP feasible */
         if ( SCIPisFeasNegative(scip, SCIPgetSolOrigObj(GCGprobGetPricingprob(scip, i), sols[j])) )
         {
            /* get variables of the pricing problem and their values in the current solution */
            probvars = SCIPgetVars(GCGprobGetPricingprob(scip, i));
            nprobvars = SCIPgetNOrigVars(GCGprobGetPricingprob(scip, i));
            printf("a) scip = %p\n", scip);
            SCIP_CALL( SCIPgetSolVals(GCGprobGetPricingprob(scip, i), sols[j], nprobvars, probvars, solvals) );
            
            printf("b) scip = %p\n", scip);
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
               if ( !SCIPisFeasZero(scip, solvals[k]) )
               {
                  vardata = SCIPvarGetData(probvars[k]);
                  assert(vardata->vartype == GCG_VARTYPE_PRICING);
                  assert(vardata->data.pricingvardata.origvar != NULL);
                  objcoeff += solvals[k] * SCIPvarGetObj(vardata->data.pricingvardata.origvar);
                  newvardata->data.mastervardata.origvars[k] = vardata->data.pricingvardata.origvar;
                  newvardata->data.mastervardata.vals[k] = solvals[k];
               }
            }


            SCIP_CALL( SCIPcreateVar(scip, &newvar, "p", 
                  0, SCIPinfinity(scip), objcoeff, SCIP_VARTYPE_CONTINUOUS, 
                  TRUE, TRUE, gcgdelvarorig, NULL, NULL, newvardata) );
            SCIP_CALL( SCIPaddPricedVar(scip, newvar, 1.0) );
            printf("added new variable\n");

            
            for ( k = 0; k < norigconss; k++ )
            {
               conscoeff = 0;
               printf("nprobvars = %d\n", nprobvars);
               for ( l = 0; l < nprobvars; l++ )
               {
                  printf("solval = %f\n", solvals[l]);
                  if ( !SCIPisFeasZero(scip, solvals[l]) )
                  {
                     vardata = SCIPvarGetData(probvars[j]);
                     assert(vardata != NULL);
                     assert(vardata->vartype == GCG_VARTYPE_PRICING);
                     vardata = SCIPvarGetData(vardata->data.pricingvardata.origvar);
                     assert(vardata != NULL);
                     assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
                     assert(vardata->data.origvardata.coefs != NULL);
                                          
                     conscoeff += vardata->data.origvardata.coefs[k] * solvals[l];
                  }
               }
               SCIP_CALL( SCIPaddCoefLinear(scip, masterconss[k], newvar, conscoeff) );
               printf("new variable has coef = %f in constraint %d\n", conscoeff, k);
               assert(0);
            }

         }
      }


      SCIP_CALL( SCIPfreeTransform(GCGprobGetPricingprob(scip, i)) );
   }
   
   SCIPfreeBufferArray(scip, &farkas);   

   return SCIP_OKAY;
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
