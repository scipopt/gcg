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
};


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
   SCIP_CONS** conss;
   int nconss;
   SCIP_Real* farkas;
   int              i;

   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   GCGprobGetMasterConss(scip, &conss, &nconss);
   SCIP_CALL( SCIPallocBufferArray(scip, &farkas, nconss) );
   
   for ( i = 0; i < nconss; i++ )
   {
      farkas[i] = SCIPgetDualfarkasLinear(scip, conss[i]);
      printf("%f\n", farkas[i]);
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
