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

/**@file   probdata_gcg.c
 * @brief  problem data for generic column generation
 * @author Gerald Gamrath
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "probdata_gcg.h"
#include "pricer_gcg.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_setppc.h"
#include <assert.h>
#include <string.h>





struct SCIP_ProbData
{



};






/*
 * Local methods
 */



/*
 * Callback methods of probdata
 */

/** transforms the problem */
static
SCIP_DECL_PROBTRANS(probtransGcg)
{
   assert(scip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, targetdata) );

   (*targetdata)->masterprob = sourcedata->masterprob;
   (*targetdata)->pricingprobs = sourcedata->pricingprobs;
   (*targetdata)->npricingprobs = sourcedata->npricingprobs;
   (*targetdata)->origmasterconss = sourcedata->origmasterconss;
   (*targetdata)->linearmasterconss = sourcedata->linearmasterconss;
   (*targetdata)->maxorigmasterconss = sourcedata->maxorigmasterconss;
   (*targetdata)->norigmasterconss = sourcedata->norigmasterconss;

   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->masterconss), sourcedata->maxmasterconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->convconss), sourcedata->npricingprobs) );
   (*targetdata)->maxmasterconss = sourcedata->maxmasterconss;
   (*targetdata)->nmasterconss = sourcedata->nmasterconss;
   
   SCIP_CALL( SCIPtransformConss(scip, sourcedata->nmasterconss, 
         sourcedata->masterconss, (*targetdata)->masterconss) );

   SCIP_CALL( SCIPtransformConss(scip, sourcedata->npricingprobs, 
         sourcedata->convconss, (*targetdata)->convconss) );

   return SCIP_OKAY;
}


/** deletes the transformed problem */
static
SCIP_DECL_PROBDELTRANS(probdeltransGcg)
{
   int i;

   assert(scip !=  NULL);
   assert(probdata != NULL);
   assert(*probdata != NULL);

   for ( i = 0; i < (*probdata)->nmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->masterconss[i]) );
   }
   for ( i = 0; i < (*probdata)->npricingprobs; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->convconss[i]) );
   }

   SCIPfreeMemoryArray(scip, &(*probdata)->masterconss);
   SCIPfreeMemoryArray(scip, &(*probdata)->convconss);

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


#define probinitsolGcg NULL
#define probexitsolGcg NULL


static
SCIP_DECL_PROBDELORIG(probdelorigGcg)
{
   int i;
   SCIP_VAR** vars;
   int nvars;
   SCIP_VARDATA* vardata;

   assert(probdata != NULL);
   assert(*probdata != NULL);

   for ( i = 0; i < (*probdata)->norigmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->origmasterconss[i]) );
   }
   for ( i = 0; i < (*probdata)->norigmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->linearmasterconss[i]) );
   }
   for ( i = 0; i < (*probdata)->nmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons((*probdata)->masterprob, &(*probdata)->masterconss[i]) );
   }
   for ( i = 0; i < (*probdata)->npricingprobs; i++ )
   {
      SCIP_CALL( SCIPreleaseCons((*probdata)->masterprob, &(*probdata)->convconss[i]) );
   }

   /* free pricing problems */
   for ( i = (*probdata)->npricingprobs - 1; i >= 0 ; i-- )
   {
      SCIP_CALL( SCIPfreeTransform((*probdata)->pricingprobs[i]) );
      SCIP_CALL( SCIPfree(&((*probdata)->pricingprobs[i])) );
   }


   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for ( i = 0; i < nvars; i++ )
   {
      vardata = SCIPvarGetData(vars[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      //printf("free vardata->data.origvardata.coefs with length");
      SCIPfreeBlockMemoryArray(scip, &(vardata->data.origvardata.coefs), (*probdata)->norigmasterconss);
   }
   
   /* free master problem */
   SCIP_CALL( SCIPfreeTransform((*probdata)->masterprob) );
   SCIP_CALL( SCIPfree(&((*probdata)->masterprob)) );



   /* free array for constraints */
   SCIPfreeMemoryArray(scip, &((*probdata)->masterconss));
   SCIPfreeMemoryArray(scip, &((*probdata)->origmasterconss));
   SCIPfreeMemoryArray(scip, &((*probdata)->linearmasterconss));

   SCIPfreeMemoryArray(scip, &((*probdata)->pricingprobs));

   SCIPfreeMemoryArray(scip, &((*probdata)->convconss));

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}




/*
 * probdata specific interface methods
 */

/** sets up the problem data */
SCIP_RETCODE SCIPcreateProbGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< problem name */
   )
{
   char probname[SCIP_MAXSTRLEN];
   SCIP_PROBDATA* probdata;
   
   printf("Creating problem: %s\n", name);
   
   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, &probdata) );

   /* initialize data */
   probdata->maxmasterconss = 5;
   probdata->nmasterconss = 0;
   probdata->maxorigmasterconss = 5;
   probdata->norigmasterconss = 0;


   /* arrays of constraints belonging to the master problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->masterconss), probdata->maxmasterconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->origmasterconss), probdata->maxorigmasterconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->linearmasterconss), probdata->maxorigmasterconss) );

   /* initializing the scip data structure for the master problem */  
   SCIP_CALL( SCIPcreate(&(probdata->masterprob)) );
   SCIP_CALL( SCIPincludeDefaultPlugins(probdata->masterprob) );
 
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "master_%s", name);

   SCIP_CALL( SCIPcreateProb(probdata->masterprob, probname, NULL, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPactivatePricer(probdata->masterprob, SCIPfindPricer(scip, "gcg")) );


   /* create original problem in SCIP */
   SCIP_CALL( SCIPcreateProb(scip, name, probdelorigGcg, probtransGcg, probdeltransGcg, 
         probinitsolGcg, probexitsolGcg, probdata) );

   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) ); 

   return SCIP_OKAY;
}


/* ----------------------------------- external methods -------------------------- */


#if 0
/** sets up the problem data */
SCIP_RETCODE GCGprobCreateFramework(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nblocks             /**< number of blocks */
   )
{
   SCIP* origprob;
   SCIP_PROBDATA* probdata;
   char name[SCIP_MAXSTRLEN];
   SCIP_VARDATA* vardata;
   SCIP_CONSHDLR** conshdlrs;
   int nconshdlrs;
   int nactiveconss;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nvars;
   SCIP_CONS* newcons;
   SCIP_CONS* mastercons;
   SCIP_Bool success;
   int i;
   int c;
   int v;
   int b;
   SCIP_HASHMAP**   hashorig2pricingvar; /* hashmap mapping original variables to corresponding pricing variables */

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   origprob = probdata->origprob;
   assert(origprob != NULL);

   printf("Creating framework for problem: %s\n", SCIPgetProbName(origprob));
      
      
   /* initialize the pricing problems */
   probdata->npricingprobs = nblocks;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->pricingprobs), nblocks) );
   /* array for saving convexity constraints belonging to one of the pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->convconss), nblocks) );

   for ( i = 0; i < nblocks; i++ )
   {
      /* initializing the scip data structure for the original problem */  
      SCIP_CALL( SCIPcreate(&(probdata->pricingprobs[i])) );
      SCIP_CALL( SCIPincludeDefaultPlugins(probdata->pricingprobs[i]) );
 
      /* disable conflict analysis */
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/useprop", FALSE) ); 
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/useinflp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/useboundlp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/usesb", FALSE) ); 
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/usepseudo", FALSE) );
      /* disable output to console */
      SCIP_CALL( SCIPsetIntParam(probdata->pricingprobs[i], "display/verblevel", 0) );
      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "misc/catchctrlc", FALSE) );

      /* create the pricing submip */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricing_block_%d", i);
      SCIP_CALL( SCIPcreateProb(probdata->pricingprobs[i], name, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* create the corresponding convexity constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conv_block_%d", i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &(probdata->convconss[i]), name, 0, NULL, NULL, 1, 1, TRUE, 
            TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, probdata->convconss[i]) );
   }

   /* presolve original problem */
   SCIPpresolve(origprob);
   if ( SCIPisObjIntegral(origprob) )
      SCIP_CALL( SCIPsetObjIntegral(scip) );
   SCIP_CALL( SCIPsetLongintParam(origprob, "limits/nodes", 1) );
   SCIP_CALL( SCIPsetIntParam(origprob, "separating/maxroundsroot", 0) );
   SCIPsolve(origprob);
   //SCIP_CALL( SCIPprintStatistics(origprob, NULL) );
   SCIP_CALL( SCIPsetIntParam(origprob, "separating/maxroundsroot", -1) );

   /* create hashmaps for mapping from original to pricing variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &hashorig2pricingvar, nblocks+1) );
   for ( i = 0; i < nblocks+1; i++ )
   {
      SCIP_CALL( SCIPhashmapCreate(&(hashorig2pricingvar[i]), SCIPblkmem(scip), SCIPgetNVars(origprob)) );
   }
   /* create pricing variables and map them to the original variables */
   vars = SCIPgetVars(origprob);
   nvars = SCIPgetNVars(origprob);
   for ( v = 0; v < nvars; v++ )
   {
      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      if ( vardata->blocknr != -1 )
      {
         if ( vardata->data.origvardata.pricingvar == NULL )
         {
            SCIP_CALL( GCGcreatePricingVar(scip, vars[v], vardata->blocknr) );
            assert(vardata->data.origvardata.pricingvar != NULL);
            SCIP_CALL( SCIPhashmapInsert(hashorig2pricingvar[vardata->blocknr], 
                  (void*)(vars[v]), (void*)(vardata->data.origvardata.pricingvar)) );
            SCIP_CALL( SCIPhashmapInsert(hashorig2pricingvar[nblocks], 
                  (void*)(vars[v]), (void*)(vars[v])) );
         }
      }
      else
      {
         assert(vardata->data.origvardata.pricingvar == NULL);
         SCIP_CALL( SCIPhashmapInsert(hashorig2pricingvar[nblocks], 
               (void*)(vars[v]), (void*)(vars[v])) );
      }
   }

   /* copy constraints of the original problem into master/pricing problems */
   conshdlrs = SCIPgetConshdlrs(origprob);
   nconshdlrs = SCIPgetNConshdlrs(origprob);

   /* iterate over all constraint handlers */
   for ( i = 0; i < nconshdlrs; i++ )
   {
      /* if there are constraints managed by this constraint handler, iterate over these constraints */
      nactiveconss = SCIPconshdlrGetNConss(conshdlrs[i]);
      if ( nactiveconss > 0 )
      {
         conss = SCIPconshdlrGetConss(conshdlrs[i]);
         for ( c = 0; c < nactiveconss; c++ )
         {
            success = FALSE;
            for ( b = 0; b < nblocks && !success; b++ )
            {
               /* try to copy the constraint */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "p%d_%s", b, SCIPconsGetName(conss[c]));
               SCIP_CALL( SCIPcopyCons(probdata->pricingprobs[b], &newcons, name, conshdlrs[i],
                     origprob, conss[c], hashorig2pricingvar[b], 
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, &success) );

               if ( success )
               {
                  SCIP_CALL( SCIPaddCons(probdata->pricingprobs[b], newcons) );

                  SCIP_CALL( SCIPreleaseCons(probdata->pricingprobs[b], &newcons) );
               }
            } 
            if ( !success )
            {
               /* copy the constraint (dirty trick, we only need lhs and rhs, because variables are added later) */
               SCIP_CALL( SCIPcopyCons(origprob, &newcons, NULL, conshdlrs[i],
                     origprob, conss[c], hashorig2pricingvar[nblocks], 
                     FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

               assert(success);

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "m_%s", SCIPconsGetName(conss[c]));
               SCIP_CALL( SCIPcreateConsLinear(scip, &mastercons, name, 0, NULL, NULL, SCIPgetLhsLinear(scip, newcons), 
                     SCIPgetRhsLinear(scip, newcons), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

               SCIP_CALL( SCIPaddCons(scip, mastercons) );

               /* store the constraints in the arrays origmasterconss and masterconss in the problem data */
               SCIP_CALL( ensureSizeMasterConss(scip, probdata, probdata->norigmasterconss+1) );
               SCIP_CALL( SCIPcaptureCons(origprob, conss[c]) );
               probdata->origmasterconss[probdata->norigmasterconss] = conss[c];
               probdata->linearmasterconss[probdata->norigmasterconss] = newcons;
               probdata->norigmasterconss++;
               probdata->masterconss[probdata->nmasterconss] = mastercons;
               probdata->nmasterconss++;

            }
            
         }
      }
   }

   /* for original variables, save the coefficiants in the master problem in their vardata */
   vars = SCIPgetVars(origprob);
   nvars = SCIPgetNVars(origprob);
   for ( v = 0; v < nvars; v++ )
   {
      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->data.origvardata.coefs == NULL);

      /* create array for saving all the coefficiants of this variable for all the constraints */
      SCIP_CALL( SCIPallocBlockMemoryArray(origprob, &(vardata->data.origvardata.coefs), probdata->norigmasterconss) );
      vardata->data.origvardata.ncoefs = probdata->norigmasterconss;
      for ( i = 0; i < vardata->data.origvardata.ncoefs; i++ )
      {
         vardata->data.origvardata.coefs[i] = 0;
      }
   }

   /* save coefs in the vardata */   
   for ( i = 0; i < probdata->norigmasterconss; i++ )
   {
      vars = SCIPgetVarsLinear(origprob, probdata->linearmasterconss[i]);
      nvars = SCIPgetNVarsLinear(origprob, probdata->linearmasterconss[i]);
      vals = SCIPgetValsLinear(origprob, probdata->linearmasterconss[i]);
      for ( v = 0; v < nvars; v++ )
      {
         vardata = SCIPvarGetData(vars[v]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         assert(vardata->data.origvardata.coefs != NULL);
         vardata->data.origvardata.coefs[i] = vals[v];
      }

   }

   /* for variables that do not belong to any block, create the corresponding variable in the master problem */
   vars = SCIPgetVars(origprob);
   nvars = SCIPgetNVars(origprob);
   for ( v = 0; v < nvars; v++ )
   {
      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      if ( vardata->blocknr == -1 )
      {
         SCIP_VAR* newvar;
         SCIP_VARDATA* newvardata;

         assert(vardata->data.origvardata.pricingvar == NULL);

         printf("var %s is in no block!\n", SCIPvarGetName(vars[i]));
         assert(SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(vars[i]) == SCIP_VARTYPE_IMPLINT);

         /* create vardata */
         SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
         newvardata->vartype = GCG_VARTYPE_MASTER;
         newvardata->blocknr = -1;
         newvardata->data.mastervardata.norigvars = 1;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvars), 
               vardata->data.mastervardata.norigvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvals), 
               vardata->data.mastervardata.norigvars) );
         newvardata->data.mastervardata.origvars[0] = vars[v];
         newvardata->data.mastervardata.origvals[0] = 1.0;

         SCIP_CALL( SCIPcreateVar(scip, &newvar, SCIPvarGetName(vars[i]), 
               SCIPvarGetLbGlobal(vars[i]), SCIPvarGetUbGlobal(vars[i]), SCIPvarGetObj(vars[i]), SCIPvarGetType(vars[i]), 
               TRUE, TRUE, gcgvardelorig, NULL, NULL, newvardata) );
         SCIPaddVar(scip, newvar);

         for ( i = 0; i < vardata->data.origvardata.ncoefs; i++ )
         {
            if ( !SCIPisFeasZero(scip, vardata->data.origvardata.coefs[i]) )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, probdata->masterconss[i], newvar, vardata->data.origvardata.coefs[i]) );
            }
         }
         SCIPreleaseVar(scip, &newvar);

      }
   }

   /* free hashmaps for mapping from original to pricing variables */
   for ( i = 0; i < nblocks+1; i++ )
   {
      SCIPhashmapFree(&(hashorig2pricingvar[i]));
   }
   SCIPfreeMemoryArray(scip, &hashorig2pricingvar);

#if 1
   /* create initial solution */
   if ( SCIPgetNSols(origprob) > 0 )
   {
      SCIP_SOL* sol;
      
      printf("creating starting variables...\n");

      /* get feasible primal solution of the original problem */
      sol = SCIPgetBestSol(origprob);
      assert(sol != NULL);

      /* create solutions for all the pricing problems */
      for ( i = 0; i < probdata->npricingprobs; i++ )
      {
         SCIP_VAR* newvar;
         SCIP_VARDATA* newvardata;
         SCIP_VAR** probvars;
         SCIP_VAR** origvars;
         SCIP_Real* solvals;
         SCIP_Real objcoeff;
         SCIP_Real conscoeff;
         int nprobvars;
         int k;

         /* get variables of the pricing problem */
         probvars = SCIPgetOrigVars(probdata->pricingprobs[i]);
         nprobvars = SCIPgetNOrigVars(probdata->pricingprobs[i]);
         SCIP_CALL( SCIPallocBufferArray(scip, &origvars, nprobvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nprobvars) );
         for ( v = 0; v < nprobvars; v++ )
         {
            vardata = SCIPvarGetData(probvars[v]);
            assert(vardata->vartype == GCG_VARTYPE_PRICING);
            origvars[v] = vardata->data.pricingvardata.origvar;
         }

         SCIP_CALL( SCIPgetSolVals(origprob, sol, nprobvars, origvars, solvals) );

         /* create data for the new variable in the master problem */
         SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
         newvardata->vartype = GCG_VARTYPE_MASTER;
         newvardata->blocknr = i;
         newvardata->data.mastervardata.norigvars = nprobvars;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvars), nprobvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvals), nprobvars) );

         /* compute objective coefficient of the variable */
         objcoeff = 0;
         for ( k = 0; k < nprobvars; k++ )
         {
            if ( !SCIPisFeasZero(scip, solvals[k]) )
            {
               vardata = SCIPvarGetData(probvars[k]);
               assert(vardata->vartype == GCG_VARTYPE_PRICING);
               assert(vardata->data.pricingvardata.origvar != NULL);
               /* add quota of original variable's objcoef to the master variable's coef */
               objcoeff += solvals[k] * SCIPvarGetObj(vardata->data.pricingvardata.origvar);
            }
         }

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "p_%d_init", i);
         
         /* create variable in the master problem */
         SCIP_CALL( SCIPcreateVar(scip, &newvar, name, 
               0, 1, objcoeff, SCIP_VARTYPE_CONTINUOUS, 
               TRUE, TRUE, gcgvardelorig, NULL, NULL, newvardata) );
         
         /* update variable datas */
         for ( k = 0; k < nprobvars; k++ )
         {
            vardata = SCIPvarGetData(probvars[k]);
            assert(vardata->vartype == GCG_VARTYPE_PRICING);
            assert(vardata->data.pricingvardata.origvar != NULL);

            if ( !SCIPisFeasZero(scip, solvals[k]) )
            {
               /* save in the master problem variable's data the quota of the corresponding original variable */
               newvardata->data.mastervardata.origvars[k] = vardata->data.pricingvardata.origvar;
               newvardata->data.mastervardata.origvals[k] = solvals[k];
               /* save the quota in the original variable's data */
               SCIP_CALL( GCGpricerAddMasterVarToOrigVar(scip, vardata->data.pricingvardata.origvar, newvar, solvals[k]) );
            }
            else
            {
               /** @todo really store variable not connected to the master variable? */
               newvardata->data.mastervardata.origvars[k] = vardata->data.pricingvardata.origvar;
               newvardata->data.mastervardata.origvals[k] = 0.0;
            }
         }

         SCIP_CALL( SCIPaddVar(scip, newvar) );
         
         SCIPchgVarUbLazy(scip, newvar, 1.0);

         for ( c = 0; c < probdata->norigmasterconss; c++ )
         {
            conscoeff = 0;
            for ( k = 0; k < nprobvars; k++ )
            {
               if ( !SCIPisFeasZero(scip, solvals[k]) )
               {
                  //printf("var = %s\n", SCIPvarPrintName(probvars[k]));
                  vardata = SCIPvarGetData(probvars[k]);
                  assert(vardata != NULL);
                  assert(vardata->vartype == GCG_VARTYPE_PRICING);
                  vardata = SCIPvarGetData(vardata->data.pricingvardata.origvar);
                  assert(vardata != NULL);
                  assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
                  assert(vardata->data.origvardata.coefs != NULL);
                                          
                  conscoeff += vardata->data.origvardata.coefs[c] * solvals[k];
               }
            }
            SCIP_CALL( SCIPaddCoefLinear(scip, probdata->masterconss[c], newvar, conscoeff) );
            //printf("new variable has coef = %f in constraint %d\n", conscoeff, k);
         }
         SCIP_CALL( SCIPaddCoefLinear(scip, GCGprobGetConvCons(scip, i), newvar, 1) );
         //printf("new variable has coef = %f in convexity constraint %d\n", 1.0, prob);

         SCIPreleaseVar(scip, &newvar);
         SCIPfreeBufferArray(scip, &origvars);
         SCIPfreeBufferArray(scip, &solvals);

      } /* for ( i = 0; i < probdata->npricingprobs; i++ ) */

      
   } /* if ( SCIPgetNSols(origprob) > 0 ) */
#endif

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, "gcg")) );

   return SCIP_OKAY;
}
#endif











/** gets values of multiple original variables w.r.t. primal master solution */
SCIP_RETCODE GCGgetSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables to get solution value for */
   SCIP_VAR**            vars,               /**< array with variables to get value for */
   SCIP_Real*            vals                /**< array to store solution values of variables */
   )
{
   SCIP_VARDATA* vardata;
   int v;
   int i;
   
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   for ( v = 0; v < nvars; v++ )
   {
      vals[v] = 0.0;
      if ( !SCIPvarIsNegated(vars[v]) )
         vardata = SCIPvarGetData(vars[v]);
      else
         vardata = SCIPvarGetData(SCIPvarGetNegatedVar(vars[v]));
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      for ( i = 0; i < vardata->data.origvardata.nmastervars; i++ )
      {
         vals[v] += vardata->data.origvardata.mastervals[i] * 
            SCIPgetSolVal(scip, sol, vardata->data.origvardata.mastervars[i]);
      }
      if ( SCIPvarIsNegated(vars[v]) )
      {
         vals[v] = 1.0 - vals[v];
      }
   }

   return SCIP_OKAY;
}
