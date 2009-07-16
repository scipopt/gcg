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

/**@file    relax_gcg.c
 * @ingroup RELAXATORS
 * @brief   gcg relaxator
 * @author  Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "relax_gcg.h"




#define RELAX_NAME             "gcg"
#define RELAX_DESC             "relaxator for gcg project representing the master lp"
#define RELAX_PRIORITY         1
#define RELAX_FREQ             1

#define STARTMAXMASTERVARS 10


/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   SCIP*            masterprob;          /* the master problem */
   SCIP**           pricingprobs;        /* the array of pricing problems */
   int              npricingprobs;       /* the number of pricing problems */

   SCIP_CONS**      convconss;           /* array of convexity constraints, one for each block */


   SCIP_HASHMAP**   hashorig2pricingvar; /* hashmap mapping original variables to corresponding 
                                          * pricing variables */
   SCIP_HASHMAP*    hashorig2origvar;    /* hashmap mapping original variables to corresponding 
                                          * themselves */

   SCIP_CONS**      masterconss;         /* array of cons in the master problem */
   SCIP_CONS**      origmasterconss;     /* array of cons in the original problem that belong to the 
                                          * master problem */
   SCIP_CONS**      linearmasterconss;   /* array of linear constraints equivalent to the cons in 
                                          * the original problem that belong to the master problem */
   int              maxmasterconss;      /* length of the array mastercons */
   int              nmasterconss;        /* number of constraints saved in mastercons */
   
   SCIP_SOL*        currentorigsol;

   int lastnrows;
};



/*
 * Vardata methods
 */

static
SCIP_DECL_VARDELORIG(gcgvardelorig)
{
   if ( (*vardata)->vartype == GCG_VARTYPE_ORIGINAL )
   {
      SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.mastervars));
      SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.mastervals));
      if ((*vardata)->data.origvardata.coefs != NULL)
      {
         SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.coefs));
         /* ,(*vardata)->data.origvardata.ncoefs);*/
      }
      
   }
   if ( (*vardata)->vartype == GCG_VARTYPE_MASTER )
   {
      assert((*vardata)->data.mastervardata.norigvars == 1);
      SCIPfreeBlockMemoryArray(scip, &((*vardata)->data.mastervardata.origvars), 2);
      SCIPfreeBlockMemoryArray(scip, &((*vardata)->data.mastervardata.origvals), 2);
   }
   SCIPfreeBlockMemory(scip, vardata);


   return SCIP_OKAY;
}  


/*
 * Local methods
 */

static
SCIP_RETCODE ensureSizeMasterConss(
   SCIP*                 scip,
   SCIP_RELAXDATA*       relaxdata,
   int                   size
   )
{
   assert(scip != NULL);
   assert(relaxdata != NULL);
   assert(relaxdata->masterconss != NULL);

   if ( relaxdata->maxmasterconss < size )
   {
      relaxdata->maxmasterconss = MAX(relaxdata->maxmasterconss + 5, size);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(relaxdata->masterconss), relaxdata->maxmasterconss) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(relaxdata->origmasterconss), relaxdata->maxmasterconss) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(relaxdata->linearmasterconss), relaxdata->maxmasterconss) );
   }
   assert(relaxdata->maxmasterconss >= size);

   return SCIP_OKAY;
}


static
SCIP_RETCODE updateCurrentSol(
   SCIP*                 scip,
   SCIP_RELAX*           relax
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_VARDATA* vardata;
   SCIP_VAR** vars;
   int nvars;
   int v;
   int i;
   SCIP_Real val;
   SCIP_Bool stored;

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   assert(vars != NULL);

   /* free previous solution */
   if ( relaxdata->currentorigsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &(relaxdata->currentorigsol)) );
   }

   /* create new solution */
   SCIP_CALL( SCIPcreateSol(scip, &relaxdata->currentorigsol, NULL) );

   for ( v = 0; v < nvars; v++ )
   {
      val = 0.0;
      if ( !SCIPvarIsNegated(vars[v]) )
         vardata = SCIPvarGetData(vars[v]);
      else
         vardata = SCIPvarGetData(SCIPvarGetNegatedVar(vars[v]));
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      for ( i = 0; i < vardata->data.origvardata.nmastervars; i++ )
      {
         val += vardata->data.origvardata.mastervals[i] * 
            SCIPgetSolVal(relaxdata->masterprob, NULL, 
               vardata->data.origvardata.mastervars[i]);
      }
      if ( SCIPvarIsNegated(vars[v]) )
      {
         val = 1.0 - val;
      }
      SCIP_CALL( SCIPsetSolVal(scip, relaxdata->currentorigsol, vars[v], val) );
   }

   SCIP_CALL( SCIPtrySol(scip, relaxdata->currentorigsol, TRUE, TRUE, FALSE, &stored) );

   return SCIP_OKAY;
}               


static
SCIP_RETCODE flushMaster(
   SCIP*                 scip,
   SCIP*                 masterscip
   )
{
   printf("flushMaster()\n");

   return SCIP_OKAY;
}


static
SCIP_RETCODE createMaster(
   SCIP*                 scip,
   SCIP_RELAX*           relax
   )
{
   SCIP_RELAXDATA* relaxdata;
   int npricingprobs;

   char name[SCIP_MAXSTRLEN];
   SCIP_VARDATA* vardata;
   SCIP_CONSHDLR** conshdlrs;
   int nconshdlrs;
   int nactiveconss;
   SCIP_CONS** conss;
   SCIP_CONS** bufconss;
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


   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   printf("createMaster()\n");

   /* initialize relaxator data */
   relaxdata->maxmasterconss = 5;
   relaxdata->nmasterconss = 0;

   /* arrays of constraints belonging to the master problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->masterconss), relaxdata->maxmasterconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->origmasterconss), relaxdata->maxmasterconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->linearmasterconss), relaxdata->maxmasterconss) );

   /* initializing the scip data structure for the master problem */  
   SCIP_CALL( SCIPcreate(&(relaxdata->masterprob)) );
   SCIP_CALL( GCGincludeMasterPlugins(relaxdata->masterprob) );
 
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "master_%s", SCIPgetProbName(scip));

   SCIP_CALL( SCIPcreateProb(relaxdata->masterprob, name, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "display/verblevel", 0) );


   SCIP_CALL( SCIPincludePricerGcg(relaxdata->masterprob, scip) );
   SCIP_CALL( SCIPactivatePricer(relaxdata->masterprob, SCIPfindPricer(relaxdata->masterprob, "gcg")) );
   GCGnodeselMasterSetOrigscip(relaxdata->masterprob, scip);


   /* ----- initialize the pricing problems ----- */
   npricingprobs = relaxdata->npricingprobs;
   assert(npricingprobs >= 0);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->pricingprobs), npricingprobs) );
   /* array for saving convexity constraints belonging to one of the pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->convconss), npricingprobs) );

   for ( i = 0; i < npricingprobs; i++ )
   {
      /* initializing the scip data structure for the original problem */  
      SCIP_CALL( SCIPcreate(&(relaxdata->pricingprobs[i])) );
      SCIP_CALL( SCIPincludeDefaultPlugins(relaxdata->pricingprobs[i]) );
 
      /* disable conflict analysis */
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "conflict/useprop", FALSE) ); 
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "conflict/useinflp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "conflict/useboundlp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "conflict/usesb", FALSE) ); 
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "conflict/usepseudo", FALSE) );
      /* disable output to console */
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "display/verblevel", 0) );
      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "misc/catchctrlc", FALSE) );

      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "separating/maxrounds", 0) );
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "separating/maxroundsroot", 0) );

      /* create the pricing submip */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricing_block_%d", i);
      SCIP_CALL( SCIPcreateProb(relaxdata->pricingprobs[i], name, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* create the corresponding convexity constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conv_block_%d", i);
      SCIP_CALL( SCIPcreateConsLinear(relaxdata->masterprob, &(relaxdata->convconss[i]), name, 0, NULL, NULL, 
            1, 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(relaxdata->masterprob, relaxdata->convconss[i]) );
   }

   /* create hashmaps for mapping from original to pricing variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->hashorig2pricingvar), npricingprobs) );
   for ( i = 0; i < npricingprobs; i++ )
   {
      SCIP_CALL( SCIPhashmapCreate(&(relaxdata->hashorig2pricingvar[i]), 
            SCIPblkmem(scip), SCIPgetNVars(scip)) );
   }
   SCIP_CALL( SCIPhashmapCreate(&(relaxdata->hashorig2origvar), 
         SCIPblkmem(scip), 10*SCIPgetNVars(scip)) );

   /* create pricing variables and map them to the original variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for ( v = 0; v < nvars; v++ )
   {
      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      if ( vardata->blocknr != -1 )
      {
         assert(vardata->data.origvardata.pricingvar == NULL);
         if ( vardata->data.origvardata.pricingvar == NULL )
         {
            SCIP_CALL( GCGrelaxCreatePricingVar(scip, vars[v]) );
            assert(vardata->data.origvardata.pricingvar != NULL);
            SCIP_CALL( SCIPhashmapInsert(relaxdata->hashorig2pricingvar[vardata->blocknr], 
                  (void*)(vars[v]), (void*)(vardata->data.origvardata.pricingvar)) );
            SCIP_CALL( SCIPhashmapInsert(relaxdata->hashorig2origvar, 
                  (void*)(vars[v]), (void*)(vars[v])) );
         }
      }
      else
      {
         assert(vardata->data.origvardata.pricingvar == NULL);
         SCIP_CALL( SCIPhashmapInsert(relaxdata->hashorig2origvar, 
               (void*)(vars[v]), (void*)(vars[v])) );
      }
   }

   /* copy constraints of the original problem into master/pricing problems */
   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);

   /* iterate over all constraint handlers */
   for ( i = 0; i < nconshdlrs; i++ )
   {
      /* if there are constraints managed by this constraint handler, iterate over these constraints */
      nactiveconss = SCIPconshdlrGetNConss(conshdlrs[i]);
      if ( nactiveconss > 0 )
      {
         conss = SCIPconshdlrGetConss(conshdlrs[i]);

         /* copy conss array */
         SCIP_CALL( SCIPallocBufferArray(scip, &bufconss, nactiveconss) );
         for ( c = 0; c < nactiveconss; c++ )
         {
            bufconss[c] = conss[c];
         }
         for ( c = 0; c < nactiveconss; c++ )
         {
            success = FALSE;
            assert(bufconss[0] == SCIPconshdlrGetConss(conshdlrs[i])[0]);
            assert(bufconss[1] == SCIPconshdlrGetConss(conshdlrs[i])[1]);
            assert(bufconss[c] == SCIPconshdlrGetConss(conshdlrs[i])[c]);
            for ( b = 0; b < npricingprobs && !success; b++ )
            {
               /* try to copy the constraint */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "p%d_%s", b, SCIPconsGetName(bufconss[c]));
               SCIP_CALL( SCIPcopyCons(relaxdata->pricingprobs[b], &newcons, name, conshdlrs[i],
                     scip, bufconss[c], relaxdata->hashorig2pricingvar[b], 
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, &success) );

               if ( success )
               {
                  SCIP_CALL( SCIPaddCons(relaxdata->pricingprobs[b], newcons) );

                  SCIP_CALL( SCIPreleaseCons(relaxdata->pricingprobs[b], &newcons) );
               }
            } 
            if ( !success )
            {
               /* copy the constraint (dirty trick, we only need lhs and rhs, because variables are added later) */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "linear_%s", SCIPconsGetName(bufconss[c]));
               SCIP_CALL( SCIPcopyCons(scip, &newcons, name, conshdlrs[i],
                     scip, bufconss[c], relaxdata->hashorig2origvar, 
                     FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

               assert(success);

               /* create and add corresponding linear constraint in the master problem */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "m_%s", SCIPconsGetName(bufconss[c]));
               SCIP_CALL( SCIPcreateConsLinear(relaxdata->masterprob, &mastercons, name, 0, NULL, NULL, 
                     SCIPgetLhsLinear(scip, newcons), SCIPgetRhsLinear(scip, newcons), 
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

               SCIP_CALL( SCIPaddCons(relaxdata->masterprob, mastercons) );

               /* store the constraints in the arrays origmasterconss and masterconss in the problem data */
               SCIP_CALL( ensureSizeMasterConss(scip, relaxdata, relaxdata->nmasterconss+1) );
               SCIP_CALL( SCIPcaptureCons(scip, bufconss[c]) );
               relaxdata->origmasterconss[relaxdata->nmasterconss] = bufconss[c];
               relaxdata->linearmasterconss[relaxdata->nmasterconss] = newcons;
               relaxdata->masterconss[relaxdata->nmasterconss] = mastercons;
               relaxdata->nmasterconss++;
            }
            
         }
         SCIPfreeBufferArray(scip, &bufconss);
      }
   }

   /* for original variables, save the coefficiants in the master problem in their vardata */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for ( v = 0; v < nvars; v++ )
   {
      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->data.origvardata.coefs == NULL);

      /* create array for saving all the coefficiants of this variable for all the constraints */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.origvardata.coefs), 
            relaxdata->nmasterconss) );
      vardata->data.origvardata.ncoefs = relaxdata->nmasterconss;
      for ( i = 0; i < vardata->data.origvardata.ncoefs; i++ )
      {
         vardata->data.origvardata.coefs[i] = 0;
      }
   }

   /* save coefs in the vardata */   
   for ( i = 0; i < relaxdata->nmasterconss; i++ )
   {
      vars = SCIPgetVarsLinear(scip, relaxdata->linearmasterconss[i]);
      nvars = SCIPgetNVarsLinear(scip, relaxdata->linearmasterconss[i]);
      vals = SCIPgetValsLinear(scip, relaxdata->linearmasterconss[i]);
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
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for ( v = 0; v < nvars; v++ )
   {
      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      if ( vardata->blocknr == -1 )
      {
         SCIP_VAR* newvar;
         SCIP_VARDATA* newvardata;

         assert(vardata->data.origvardata.pricingvar == NULL);

         printf("var %s is in no block!\n", SCIPvarGetName(vars[v]));
         //assert(SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(vars[i]) == SCIP_VARTYPE_IMPLINT);

         /* create vardata */
         SCIP_CALL( SCIPallocBlockMemory(relaxdata->masterprob, &newvardata) );
         newvardata->vartype = GCG_VARTYPE_MASTER;
         newvardata->blocknr = -1;
         newvardata->data.mastervardata.norigvars = 1;
#if 0
         SCIP_CALL( SCIPallocBlockMemoryArray(relaxdata->masterprob, 
               &(newvardata->data.mastervardata.origvars), 2) );
         SCIP_CALL( SCIPallocBlockMemoryArray(relaxdata->masterprob, 
               &(newvardata->data.mastervardata.origvals), 2) );
         newvardata->data.mastervardata.origvars[0] = vars[v];
         newvardata->data.mastervardata.origvals[0] = 1.0;
#endif
         SCIP_CALL( SCIPcreateVar(relaxdata->masterprob, &newvar, SCIPvarGetName(vars[v]), 
               SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v]), SCIPvarGetObj(vars[v]), SCIPvarGetType(vars[v]), 
               TRUE, TRUE, gcgvardelorig, NULL, NULL, newvardata) );
         SCIPaddVar(relaxdata->masterprob, newvar);

         for ( i = 0; i < vardata->data.origvardata.ncoefs; i++ )
         {
            if ( !SCIPisFeasZero(scip, vardata->data.origvardata.coefs[i]) )
            {
               SCIP_CALL( SCIPaddCoefLinear(relaxdata->masterprob, relaxdata->masterconss[i], 
                     newvar, vardata->data.origvardata.coefs[i]) );
            }
         }
         SCIPreleaseVar(relaxdata->masterprob, &newvar);

      }
   }


   return SCIP_OKAY;
}




/*
 * Callback methods of relaxator
 */

/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreeGcg)
{  
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   SCIPfreeMemory(scip, &relaxdata);

   return SCIP_OKAY;
}



/** initialization method of relaxator (called after problem was transformed) */
static
SCIP_DECL_RELAXINIT(relaxInitGcg)
{  
   return SCIP_OKAY;
}



/** deinitialization method of relaxator (called before transformed problem is freed) */

static
SCIP_DECL_RELAXEXIT(relaxExitGcg)
{  
   SCIP_RELAXDATA* relaxdata;
   int i;

   assert(scip != NULL);
   assert(relax != NULL);
   
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* free hashmaps for mapping from original to pricing variables */
   if ( relaxdata->hashorig2pricingvar != NULL )
   {
      for ( i = 0; i < relaxdata->npricingprobs; i++ )
      {
         SCIPhashmapFree(&(relaxdata->hashorig2pricingvar[i]));
      }
      SCIPfreeMemoryArray(scip, &(relaxdata->hashorig2pricingvar));
      relaxdata->hashorig2pricingvar = NULL;
   }
   if ( relaxdata->hashorig2origvar != NULL )
   {
      SCIPhashmapFree(&(relaxdata->hashorig2origvar));
      relaxdata->hashorig2origvar = NULL;
   }

   /* free array for constraints */
   for ( i = 0; i < relaxdata->nmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &relaxdata->origmasterconss[i]) );
   }
   for ( i = 0; i < relaxdata->nmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &relaxdata->linearmasterconss[i]) );
   }
   for ( i = 0; i < relaxdata->nmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(relaxdata->masterprob, &relaxdata->masterconss[i]) );
   }
   for ( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(relaxdata->masterprob, &relaxdata->convconss[i]) );
   }
   SCIPfreeMemoryArray(scip, &(relaxdata->origmasterconss));
   SCIPfreeMemoryArray(scip, &(relaxdata->linearmasterconss));
   SCIPfreeMemoryArray(scip, &(relaxdata->masterconss));
   SCIPfreeMemoryArray(scip, &(relaxdata->convconss));

   /* free master problem */
   SCIP_CALL( SCIPfree(&(relaxdata->masterprob)) );

   /* free pricing problems */
   for ( i = relaxdata->npricingprobs - 1; i >= 0 ; i-- )
   {
      SCIP_CALL( SCIPfreeTransform(relaxdata->pricingprobs[i]) );
      SCIP_CALL( SCIPfree(&(relaxdata->pricingprobs[i])) );
   }
   SCIPfreeMemoryArray(scip, &(relaxdata->pricingprobs));

   /* free solution */
   if ( relaxdata->currentorigsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &relaxdata->currentorigsol) );
   }

   return SCIP_OKAY;
}


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
static
SCIP_DECL_RELAXINITSOL(relaxInitsolGcg)
{  
   SCIP* masterprob;
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   SCIP_CALL( createMaster(scip, relax) );

   masterprob = relaxdata->masterprob;
   assert(masterprob != NULL);

   SCIP_CALL( SCIPtransformProb(masterprob) );

   SCIP_CALL( SCIPtransformConss(masterprob, relaxdata->nmasterconss, 
         relaxdata->masterconss, relaxdata->masterconss) );
   SCIP_CALL( SCIPtransformConss(masterprob, relaxdata->npricingprobs, 
         relaxdata->convconss, relaxdata->convconss) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolGcg)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of gcg relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxExitsolGcg NULL
#endif


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecGcg)
{  
   SCIP* masterprob;
   SCIP_RELAXDATA* relaxdata;
   SCIP_Bool delayed;
   SCIP_Bool cutoff;
   SCIP_Longint oldnnodes;

   assert(scip != NULL);
   assert(relax != NULL);
   assert(result != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   masterprob = relaxdata->masterprob;
   assert(masterprob != NULL);

   *result = SCIP_DIDNOTRUN;
   
   printf("relaxexec()\n");

#if 0
   if ( relaxdata->lastnrows < SCIPgetNLPRows(scip) )
   {
      SCIP_CALL( flushMaster(scip, masterprob) );
   }
#endif

   printf("solving node %d's relaxation!\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   SCIP_CALL( SCIPgetLongintParam(masterprob, "limits/nodes", &oldnnodes) );

   SCIP_CALL( SCIPsetLongintParam(masterprob, "limits/nodes", 
         ( SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) ? 1 : oldnnodes+1)) );
   SCIP_CALL( SCIPsolve(masterprob) );
   //SCIP_CALL( SCIPprintStatistics(masterprob, NULL) );
   
   SCIP_CALL( SCIPupdateLocalLowerbound(scip, SCIPgetSolOrigObj(masterprob, NULL)) );

   SCIP_CALL( updateCurrentSol(scip, relax) );

   //SCIP_CALL( SCIPprintSol(scip, relaxdata->currentorigsol, NULL, FALSE) );

   /* separate the optimal LP-solution */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   //printf("SCIPconstructLP: cutoff = %d\n", cutoff);
   assert(!cutoff);

   SCIP_CALL( SCIPflushLP(scip) );

   SCIP_CALL( SCIPseparateSol(scip, relaxdata->currentorigsol, FALSE, FALSE, &delayed, &cutoff) );

   //printf("delayed: %d, cutoff: %d\n", delayed, cutoff);

   printf("SCIPseparateSol() found %d cuts!\n", SCIPgetNCuts(scip));

   //SCIP_CALL( SCIPprintSol(scip, relaxdata->currentorigsol, NULL, FALSE) );


   return SCIP_OKAY;
}





/*
 * relaxator specific interface methods
 */

/** creates the gcg relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxGcg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;

   /* create gcg relaxator data */
   SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );
   
   relaxdata->lastnrows = 0;
   relaxdata->npricingprobs = -1;
   relaxdata->currentorigsol = NULL;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelax(scip, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, relaxFreeGcg, relaxInitGcg, 
         relaxExitGcg, relaxInitsolGcg, relaxExitsolGcg, relaxExecGcg, relaxdata) );

   /* inform scip, that no LPs should be solved */
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) ); 

   /* add gcg relaxator parameters */


   return SCIP_OKAY;
}


/** creates a variable in a pricing problem corresponding to the given original variable */
SCIP_RETCODE GCGrelaxCreatePricingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar             /**< corresponding variable in the original program */
   )
{
   SCIP_VARDATA* vardata;
   SCIP_VARDATA* origvardata;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_VAR* var;
   char name[SCIP_MAXSTRLEN];
   int pricingprobnr;

   assert(scip != NULL);
   assert(origvar != NULL);
   
   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* get variable data of the original variable */
   origvardata = SCIPvarGetData(origvar);
   assert(origvardata != NULL);
   assert(origvardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(origvardata->data.origvardata.pricingvar == NULL);
   assert(origvardata->blocknr != -1);

   /* get the number of the pricing block, the variable belongs to */
   pricingprobnr = origvardata->blocknr;
   assert(pricingprobnr >= 0 && pricingprobnr < relaxdata->npricingprobs);

   /* create variable data */
   SCIP_CALL( SCIPallocBlockMemory(relaxdata->pricingprobs[pricingprobnr], &vardata) );
   vardata->vartype = GCG_VARTYPE_PRICING;
   vardata->blocknr = pricingprobnr;
   vardata->data.pricingvardata.origvar = origvar;

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pr%d_%s", pricingprobnr, SCIPvarGetName(origvar));
   SCIP_CALL( SCIPcreateVar(relaxdata->pricingprobs[pricingprobnr], &var, name, 
         SCIPvarGetLbGlobal(origvar), SCIPvarGetUbGlobal(origvar), 0, SCIPvarGetType(origvar), 
         TRUE, FALSE, gcgvardelorig, NULL, NULL, vardata) );

   origvardata->data.origvardata.pricingvar = var;

   SCIP_CALL( SCIPaddVar(relaxdata->pricingprobs[pricingprobnr], var) );

   /* because the variable was added to the problem, 
    * it is captured by SCIP and we can safely release it right now
    */
   SCIP_CALL( SCIPreleaseVar(relaxdata->pricingprobs[pricingprobnr], &var) );
   
   return SCIP_OKAY;
}


/** creates the data for a variable of the original program */
SCIP_RETCODE GCGrelaxCreateOrigVardata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< pointer to variable object */
   )
{
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &vardata) );
   vardata->vartype = GCG_VARTYPE_ORIGINAL;
   vardata->blocknr = -1;
   vardata->data.origvardata.pricingvar = NULL;
   vardata->data.origvardata.coefs = NULL;
   vardata->data.origvardata.ncoefs = 0;
   vardata->data.origvardata.nmastervars = 0;
   vardata->data.origvardata.maxmastervars = STARTMAXMASTERVARS;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.origvardata.mastervars), 
         vardata->data.origvardata.maxmastervars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.origvardata.mastervals), 
         vardata->data.origvardata.maxmastervars) );

   SCIPvarSetData(var, vardata);
   SCIPvarSetDelorigData(var, gcgvardelorig);

   return SCIP_OKAY;
}

/** creates the data for all variables of the original program */
SCIP_RETCODE GCGrelaxCreateOrigVarsData(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   for ( i = 0; i < nvars; i++ )
   {
      assert(vars[i] != NULL);
      SCIP_CALL( GCGrelaxCreateOrigVardata(scip, vars[i]) );
   }

   return SCIP_OKAY;
}

/* sets the number of the block, the given original variable belongs to */
SCIP_RETCODE GCGrelaxSetOriginalVarBlockNr(
   SCIP_VAR*             var,                /**< variable to set the block number for */
   int                   blocknr             /**< number of the block, the variable belongs to */
   )
{
   SCIP_VARDATA* vardata;

   assert(SCIPvarIsOriginal(var) && SCIPvarGetTransVar(var) == NULL);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   assert(vardata->blocknr == -1);

   vardata->blocknr = blocknr;

   return SCIP_OKAY;
}

/* returns the master problem */
SCIP* GCGrelaxGetMasterprob(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->masterprob;
}

/* returns the pricing problem of the given number */
SCIP* GCGrelaxGetPricingprob(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->pricingprobs[pricingprobnr];
}

/* returns the number of pricing problems */
int GCGrelaxGetNPricingprobs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->npricingprobs;
}

/* sets the number of pricing problems */
void GCGrelaxSetNPricingprobs(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   npricingprobs       /**< the number of pricing problems */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   relaxdata->npricingprobs = npricingprobs;
}

/* returns the number of constraints in the master problem */
int GCGrelaxGetNMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->nmasterconss;
}

/* returns the contraints in the master problem */
SCIP_CONS** GCGrelaxGetMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->masterconss;
}




/* returns the contraints in the original problem that correspond to the constraints in the master problem */
SCIP_CONS** GCGrelaxGetOrigMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->origmasterconss;
}

/* returns the linear counterpart of the contraints in the original problem that correspond 
 * to the constraints in the master problem */
SCIP_CONS** GCGrelaxGetLinearOrigMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->linearmasterconss;
}

/* returns the convexity constraint for the given block */
SCIP_CONS* GCGrelaxGetConvCons(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   blocknr             /**< the number of the block for which we 
                                              *   need the convexity constraint */   
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->convconss[blocknr];
}

/* returns the current solution for the original problem */
SCIP_SOL* GCGrelaxGetCurrentOrigSol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->currentorigsol;
}
