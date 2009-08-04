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
#define SCIP_DEBUG
/**@file    relax_gcg.c
 * @ingroup RELAXATORS
 * @brief   gcg relaxator
 * @author  Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "relax_gcg.h"
#include "scip/scip.h"
#include "struct_vardata.h"




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
   int*             blockrepresentative; 
   int*             nblocksidentical;

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
   SCIP_Longint     lastmasterlpiters;
   SCIP_SOL*        lastmastersol;

   SCIP_VAR**       branchcands;
   SCIP_Real*       branchcandssol;
   SCIP_Real*       branchcandsfrac;
   int              nbranchcands;

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
   if ( (*vardata)->vartype == GCG_VARTYPE_PRICING )
   {
      assert((*vardata)->data.pricingvardata.norigvars >= 1);
      SCIPfreeMemoryArray(scip, &((*vardata)->data.pricingvardata.origvars));
   }
   assert((*vardata)->vartype != GCG_VARTYPE_MASTER);
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
SCIP_Bool intArraysAreEqual(
   int*                  array1,
   int                   array1length,
   int*                  array2,
   int                   array2length
   )
{
   int i;

   assert(array1 != NULL && array2 != NULL); 
   if ( array1length != array2length )
      return FALSE;

   for ( i = 0; i < array1length; i++ )
   {
      if ( array1[i] != array2[i] )
         return FALSE;
   }

   return TRUE;

}

static 
SCIP_Bool realArraysAreEqual(
   SCIP_Real*            array1,
   int                   array1length,
   SCIP_Real*            array2,
   int                   array2length
   )
{
   int i;

   assert(array1 != NULL && array2 != NULL); 
   if ( array1length != array2length )
      return FALSE;

   for ( i = 0; i < array1length; i++ )
   {
      if ( array1[i] != array2[i] )
         return FALSE;
   }

   return TRUE;

}

static
SCIP_Bool scipEqualsScip(
   SCIP*                 scip1,              /**< SCIP data structure */
   SCIP*                 scip2,              /**< SCIP data structure */
   SCIP_HASHMAP*         varmap
   )
{
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   int nvars1;
   int nvars2;
   
   SCIP_VARDATA* vardata1;
   SCIP_VARDATA* vardata2;

   SCIP_CONS** conss1;
   SCIP_CONS** conss2;
   int nconss;

   int i;
   int j;

   if ( SCIPgetNVars(scip1) != SCIPgetNVars(scip2))
      return FALSE;
   
   if ( SCIPgetNConss(scip1) != SCIPgetNConss(scip1))
      return FALSE;
   
   /* get variables */
   SCIP_CALL( SCIPgetVarsData(scip1, &vars1, &nvars1, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(scip2, &vars2, &nvars2, NULL, NULL, NULL, NULL) );
   
   for ( i = 0; i < nvars1; i++ )
   {
      if ( SCIPvarGetObj(vars1[i]) != SCIPvarGetObj(vars2[i]) )
         return FALSE;
      if ( SCIPvarGetLbOriginal(vars1[i]) != SCIPvarGetLbOriginal(vars2[i]) )
         return FALSE;
      if ( SCIPvarGetUbOriginal(vars1[i]) != SCIPvarGetUbOriginal(vars2[i]) )
         return FALSE;
      if ( SCIPvarGetType(vars1[i]) != SCIPvarGetType(vars2[i]) )
         return FALSE;
      
      vardata1 = SCIPvarGetData(vars1[i]);
      vardata2 = SCIPvarGetData(vars2[i]);
      assert(vardata1 != NULL && vardata2 != NULL);
      assert(vardata1->vartype == GCG_VARTYPE_PRICING);
      assert(vardata2->vartype == GCG_VARTYPE_PRICING);
      //assert(vardata1->data.pricingvardata.norigvars >= 1);
      //assert(vardata2->data.pricingvardata.norigvars >= 1);
      assert(vardata1->data.pricingvardata.origvars != NULL);
      assert(vardata2->data.pricingvardata.origvars != NULL);

      if ( SCIPvarGetObj(vardata1->data.pricingvardata.origvars[0]) 
         != SCIPvarGetObj(vardata2->data.pricingvardata.origvars[0]) )
         return FALSE;

      vardata1 = SCIPvarGetData(vardata1->data.pricingvardata.origvars[0]);
      vardata2 = SCIPvarGetData(vardata2->data.pricingvardata.origvars[0]);
      assert(vardata1 != NULL && vardata2 != NULL);
      assert(vardata1->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata1->data.origvardata.coefs != NULL);
      assert(vardata2->data.origvardata.coefs != NULL);

      if ( !realArraysAreEqual(vardata1->data.origvardata.coefs, vardata1->data.origvardata.ncoefs,
            vardata2->data.origvardata.coefs, vardata2->data.origvardata.ncoefs) )
         return FALSE;

      SCIP_CALL_ABORT( SCIPhashmapInsert(varmap, (void*) vars1[i], (void*) vars2[i]) );

   }

   /* check whether the conss are the same */
   conss1 = SCIPgetConss(scip1);
   conss2 = SCIPgetConss(scip2);
   nconss = SCIPgetNConss(scip1);
   assert(nconss == SCIPgetNConss(scip2));
   for ( i = 0; i < nconss; i++ )
   {
      if ( SCIPgetNVarsLinear(scip1, conss1[i]) != SCIPgetNVarsLinear(scip2, conss2[i]) )
         return FALSE;
      if ( SCIPgetLhsLinear(scip1, conss1[i]) != SCIPgetLhsLinear(scip2, conss2[i]) )
         return FALSE;
      if ( SCIPgetRhsLinear(scip1, conss1[i]) != SCIPgetRhsLinear(scip2, conss2[i]) )
         return FALSE;
      if ( !realArraysAreEqual(SCIPgetValsLinear(scip1, conss1[i]), SCIPgetNVarsLinear(scip1, conss1[i]),
            SCIPgetValsLinear(scip2, conss2[i]), SCIPgetNVarsLinear(scip2, conss2[i])) )
         return FALSE;
      
      vars1 = SCIPgetVarsLinear(scip1, conss1[i]);
      vars2 = SCIPgetVarsLinear(scip2, conss2[i]);
      for ( j = 0; j < SCIPgetNVarsLinear(scip1, conss1[i]); j++ )
      {
         if ( (SCIP_VAR*) SCIPhashmapGetImage(varmap, (void*) vars1[j]) != vars2[j] )
            return FALSE;
      }

   }

   return TRUE;
}

/** checks whether there are identical pricing blocks */
static
SCIP_RETCODE checkIdenticalBlocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax               /**< the relaxator */
   )
{
   SCIP_RELAXDATA* relaxdata;

   SCIP_HASHMAP* varmap;
   SCIP_VARDATA* vardata;
   SCIP_VARDATA* vardata2;
   SCIP_VAR** vars;
   SCIP_VAR* origvar;
   SCIP_VAR* pricingvar;
   int nvars;

   int i;
   int j;
   int k;
   

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   for ( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      relaxdata->blockrepresentative[i] = i;
      relaxdata->nblocksidentical[i] = 1;

      for ( j = 0; j < i && relaxdata->blockrepresentative[i] == i; j++ )
      {
         if ( relaxdata->blockrepresentative[j] != j )
            continue;

         SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), 5 * SCIPgetNVars(relaxdata->pricingprobs[i])) );

         if ( scipEqualsScip(relaxdata->pricingprobs[i], relaxdata->pricingprobs[j], varmap) )
         {
            SCIPdebugMessage("Block %d is identical to block %d!\n", i, j);
            
            /* block i is now represented by block j */
            relaxdata->blockrepresentative[i] = j;
            relaxdata->nblocksidentical[i] = 0;
            relaxdata->nblocksidentical[j]++;
            /* save variables in pricing problem variable's vardata */
            vars = SCIPgetVars(relaxdata->pricingprobs[i]);
            nvars = SCIPgetNVars(relaxdata->pricingprobs[i]);
            for ( k = 0; k < nvars; k++ )
            {
               vardata = SCIPvarGetData(vars[k]);
               assert(vardata->vartype == GCG_VARTYPE_PRICING);
               assert(vardata->data.pricingvardata.origvars != NULL);
               assert(vardata->data.pricingvardata.norigvars == 1);
               assert(vardata->data.pricingvardata.origvars[0] != NULL);
               origvar = vardata->data.pricingvardata.origvars[0];
               pricingvar = (SCIP_VAR*) SCIPhashmapGetImage(varmap, (void*) vars[k]);
               vardata = SCIPvarGetData(pricingvar);
               assert(vardata->vartype == GCG_VARTYPE_PRICING);
               assert(vardata->data.pricingvardata.origvars != NULL);
               assert(vardata->data.pricingvardata.norigvars >= 1);
               vardata2 = SCIPvarGetData(origvar);
               assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);
               assert(vardata2->data.origvardata.pricingvar != NULL);
               vardata2->data.origvardata.pricingvar = pricingvar;
               if ( vardata->data.pricingvardata.norigvars >= 2 )
               {
                  SCIP_CALL( SCIPreallocMemoryArray(relaxdata->pricingprobs[vardata->blocknr], 
                        &(vardata->data.pricingvardata.origvars), vardata->data.pricingvardata.norigvars+1) );
               }
               vardata->data.pricingvardata.origvars[vardata->data.pricingvardata.norigvars] = origvar;
               vardata->data.pricingvardata.norigvars++;
            }

         }
         SCIPhashmapFree(&varmap);
         
      }
      if ( relaxdata->blockrepresentative[i] == i )
      {
         SCIPdebugMessage("Block %d is relevant!\n", i);
      }
   }

   return SCIP_OKAY;
}



/** creates the master problem and the pricing problems and copies the constraints into them */
static
SCIP_RETCODE createMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax               /**< the relaxator */
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

   SCIPdebugMessage("Creating Master Problem...\n");

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
   //SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "display/freq", 1) );
   SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "display/curdualbound/active", 2) );

   SCIP_CALL( SCIPincludePricerGcg(relaxdata->masterprob, scip) );
   SCIP_CALL( SCIPactivatePricer(relaxdata->masterprob, SCIPfindPricer(relaxdata->masterprob, "gcg")) );


   /* ----- initialize the pricing problems ----- */
   npricingprobs = relaxdata->npricingprobs;
   assert(npricingprobs >= 0);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->pricingprobs), npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->blockrepresentative), npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->nblocksidentical), npricingprobs) );
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

      /* disable expensive presolving */
      //SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "presolving/probing/maxrounds", 0) );
      //SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/linear/presolpairwise", FALSE) );
      //SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/setppc/presolpairwise", FALSE) );
      //SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/logicor/presolpairwise", FALSE) );
      //SCIP_CALL( SCIPsetRealParam(relaxdata->pricingprobs[i], "constraints/linear/maxaggrnormscale", 0.0) );

      /* disable output to console */
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "display/verblevel", 0) );
      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "misc/catchctrlc", FALSE) );

      //SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "separating/maxrounds", 0) );
      //SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "separating/maxroundsroot", 0) );

      /* create the pricing submip */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricing_block_%d", i);
      SCIP_CALL( SCIPcreateProb(relaxdata->pricingprobs[i], name, NULL, NULL, NULL, NULL, NULL, NULL) );

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
      if ( strcmp(SCIPconshdlrGetName(conshdlrs[i]), "origbranch") == 0)
      {
         continue;
      }

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

   SCIP_CALL( checkIdenticalBlocks(scip, relax) );

   for ( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      if ( relaxdata->blockrepresentative[i] != i )
         continue;
      /* create the corresponding convexity constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conv_block_%d", i);
      SCIP_CALL( SCIPcreateConsLinear(relaxdata->masterprob, &(relaxdata->convconss[i]), name, 0, NULL, NULL, 
            relaxdata->nblocksidentical[i], relaxdata->nblocksidentical[i], 
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(relaxdata->masterprob, relaxdata->convconss[i]) );
   }


   return SCIP_OKAY;
}


/* checks the consistency between original scip and master scip */
static
SCIP_RETCODE checkConsistency(
   SCIP*                 scip                /**< SCIP data structure of the original scip */
   )
{
   SCIP_CONS** origconss;
   SCIP_CONS** masterconss;
   SCIP* masterprob;
   int norigconss;
   int nmasterconss;
   int i;
   //SCIP_Bool feasible;

   assert(scip != NULL);

   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);
   assert(SCIPgetStage(masterprob) == SCIP_STAGE_SOLVING || SCIPgetStage(masterprob) == SCIP_STAGE_SOLVED); 
   //assert(SCIPgetNSols(masterprob) == 0);
   //assert(SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == SCIPnodeGetNumber(SCIPgetCurrentNode(masterprob)));

   GCGconsOrigbranchGetStack(scip, &origconss, &norigconss);
   GCGconsMasterbranchGetStack(masterprob, &masterconss, &nmasterconss);

   assert(norigconss == nmasterconss);

   for ( i = 0; i < norigconss; i++ )
   {
      assert(origconss[i] == GCGconsMasterbranchGetOrigcons(masterconss[i]));
      assert(GCGconsOrigbranchGetOrigvar(origconss[i]) == GCGconsMasterbranchGetOrigvar(masterconss[i]));
      assert(GCGconsOrigbranchGetVal(origconss[i]) == GCGconsMasterbranchGetVal(masterconss[i]));
      assert(GCGconsOrigbranchGetConssense(origconss[i]) == GCGconsMasterbranchGetConssense(masterconss[i]));
   }

   GCGconsOrigbranchCheckConsistency(scip);
   GCGconsMasterbranchCheckConsistency(masterprob);


   //SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), TRUE, FALSE, TRUE, &feasible) );
   //assert(feasible);

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

   /* free arrays for constraints */
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
      if ( relaxdata->convconss[i] != NULL )
         SCIP_CALL( SCIPreleaseCons(relaxdata->masterprob, &relaxdata->convconss[i]) );
   }

   SCIPfreeMemoryArray(scip, &(relaxdata->origmasterconss));
   SCIPfreeMemoryArray(scip, &(relaxdata->linearmasterconss));
   SCIPfreeMemoryArray(scip, &(relaxdata->masterconss));
   SCIPfreeMemoryArray(scip, &(relaxdata->convconss));

   /* free master problem */
   SCIP_CALL( SCIPprintStatistics(relaxdata->masterprob, NULL) );
   SCIP_CALL( SCIPfree(&(relaxdata->masterprob)) );

   /* free pricing problems */
   for ( i = relaxdata->npricingprobs - 1; i >= 0 ; i-- )
   {
      SCIP_CALL( SCIPfreeTransform(relaxdata->pricingprobs[i]) );
      SCIP_CALL( SCIPfree(&(relaxdata->pricingprobs[i])) );
   }
   SCIPfreeMemoryArray(scip, &(relaxdata->pricingprobs));
   SCIPfreeMemoryArray(scip, &(relaxdata->blockrepresentative));
   SCIPfreeMemoryArray(scip, &(relaxdata->nblocksidentical));

   /* free solution */
   if ( relaxdata->currentorigsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &relaxdata->currentorigsol) );
   }

   SCIPfreeMemoryArray(scip, &(relaxdata->branchcands));
   SCIPfreeMemoryArray(scip, &(relaxdata->branchcandssol));
   SCIPfreeMemoryArray(scip, &(relaxdata->branchcandsfrac));


   return SCIP_OKAY;
}


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
static
SCIP_DECL_RELAXINITSOL(relaxInitsolGcg)
{  
   SCIP* masterprob;
   SCIP_RELAXDATA* relaxdata;
   int i;

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

   for ( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      if ( relaxdata->convconss[i] != NULL) 
      {
         SCIP_CALL( SCIPtransformCons(masterprob, relaxdata->convconss[i], &(relaxdata->convconss[i])) );
      }
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->branchcands), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->branchcandssol), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->branchcandsfrac), SCIPgetNVars(scip)) );
   relaxdata->nbranchcands = 0;

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
   
   //printf("solving node %d's relaxation!\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   /* increase the node limit for the master problem by 1 */
   SCIP_CALL( SCIPgetLongintParam(masterprob, "limits/nodes", &oldnnodes) );
   SCIP_CALL( SCIPsetLongintParam(masterprob, "limits/nodes", 
         ( SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) ? 1 : oldnnodes+1)) );

   /* separate the optimal LP-solution */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   assert(!cutoff);
   SCIP_CALL( SCIPflushLP(scip) );

   SCIPdebugMessage("Solve master LP.\n");
   /* solve the next node in the master problem */
   SCIP_CALL( SCIPsolve(masterprob) );


   /* update the lower bound of the current node */
   if ( SCIPgetStage(masterprob) == SCIP_STAGE_SOLVING )
   {
      SCIP_CALL( SCIPupdateLocalLowerbound(scip, SCIPgetSolOrigObj(masterprob, NULL)) );
   }
   else
   {
      assert(SCIPgetBestSol(masterprob) != NULL);
      SCIP_CALL( SCIPupdateLocalLowerbound(scip, SCIPgetSolOrigObj(masterprob, SCIPgetBestSol(masterprob))) );
   }
   SCIPdebugMessage("Update lower bound (value = %"SCIP_REAL_FORMAT").\n", SCIPgetLocalLowerbound(scip));


   SCIPdebugMessage("Update current sol.\n");
   /* transform the current solution of the master problem to the original space and save it */
   SCIP_CALL( GCGrelaxUpdateCurrentSol(scip) );

   //SCIP_CALL( checkConsistency(scip) );

   //SCIP_CALL( SCIPprintSol(scip, relaxdata->currentorigsol, NULL, FALSE) );

   *result = SCIP_SUCCESS;

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
   
   relaxdata->npricingprobs = -1;
   relaxdata->currentorigsol = NULL;
   relaxdata->lastmastersol = NULL;
   relaxdata->lastmasterlpiters = 0;

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
   SCIP_CALL( SCIPallocMemoryArray(relaxdata->pricingprobs[pricingprobnr], 
         &(vardata->data.pricingvardata.origvars), 2) );
   vardata->data.pricingvardata.origvars[0] = origvar;
   vardata->data.pricingvardata.norigvars = 1;

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

/** returns TRUE iff the pricingproblem of the given number is relevant, that means is not identical to
 *  another and represented by it */
SCIP_Bool GCGrelaxIsPricingprobRelevant(
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

   return (relaxdata->blockrepresentative[pricingprobnr] == pricingprobnr);

}

/** returns the number of blocks in the original formulation, that are represented by 
 *  the pricingprob with the given number */
int GCGrelaxGetNIdenticalBlocks(
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

   assert((relaxdata->blockrepresentative[pricingprobnr] == pricingprobnr) 
      == (relaxdata->nblocksidentical[pricingprobnr] > 0));

   return relaxdata->nblocksidentical[pricingprobnr];

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


/** transformes the current solution of the master problem into the original problem's space 
 *  and saves this solution as currentsol in the relaxator's data */
SCIP_RETCODE GCGrelaxUpdateCurrentSol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   SCIP_VAR** origvars;
   SCIP_Real* origvals;
   int norigvars;

   SCIP_VAR** mastervars;
   SCIP_Real* mastervals;
   int nmastervars;

   SCIP_Bool stored;
   SCIP_SOL* mastersol;

   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   origvars = SCIPgetVars(scip);
   norigvars = SCIPgetNVars(scip);
   assert(origvars != NULL);
   
   mastervars = SCIPgetVars(relaxdata->masterprob);
   nmastervars = SCIPgetNVars(relaxdata->masterprob);
   assert(mastervars != NULL);
   
   SCIP_CALL( SCIPallocBufferArray(scip, &origvals, norigvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mastervals, nmastervars) );


   /* nothing has to be done, if no LP was solved after the last update */
   if ( TRUE || relaxdata->lastmasterlpiters != SCIPgetNLPIterations(relaxdata->masterprob) )
   {
      relaxdata->lastmasterlpiters = SCIPgetNLPIterations(relaxdata->masterprob);

      /* free previous solution */
      if ( relaxdata->currentorigsol != NULL )
      {
         SCIP_CALL( SCIPfreeSol(scip, &(relaxdata->currentorigsol)) );
      }

      /* create new solution */
      SCIP_CALL( SCIPcreateSol(scip, &relaxdata->currentorigsol, NULL) );

      if ( SCIPgetStage(relaxdata->masterprob) == SCIP_STAGE_SOLVING )
      {
         mastersol = NULL;
      }
      else if ( SCIPgetStage(relaxdata->masterprob) == SCIP_STAGE_SOLVED )
      {
         mastersol = SCIPgetBestSol(relaxdata->masterprob);
      }
      else 
      {
         printf("solstat in master not solving and not solved!\n");
         return SCIP_OKAY;
      }
      
      SCIP_CALL( SCIPgetSolVals(relaxdata->masterprob, mastersol, nmastervars, mastervars, mastervals) );

      SCIP_CALL( GCGrelaxTransformMastervalsToOrigvals(scip, mastervars, mastervals, nmastervars, 
            origvars, origvals, norigvars) );

      SCIP_CALL( SCIPsetSolVals(scip, relaxdata->currentorigsol, norigvars, origvars, origvals) );

      SCIP_CALL( SCIPtrySol(scip, relaxdata->currentorigsol, TRUE, TRUE, TRUE, &stored) );

      //SCIP_CALL( SCIPprintSol(scip, relaxdata->currentorigsol, NULL, FALSE) );

      SCIPdebugMessage("updated current original LP solution, %s feasible in the original problem!\n",
         (stored ? "" : "not"));

      //SCIP_CALL( SCIPprintSol(scip, relaxdata->currentorigsol, NULL, FALSE) );
   }
   if ( SCIPgetBestSol(relaxdata->masterprob) != NULL && relaxdata->lastmastersol != SCIPgetBestSol(relaxdata->masterprob) )
   {
      SCIP_SOL* newsol;

      relaxdata->lastmastersol = SCIPgetBestSol(relaxdata->masterprob);

      SCIP_CALL( SCIPcreateSol(scip, &newsol, NULL) );
      assert(newsol != NULL);

      SCIP_CALL( SCIPgetSolVals(relaxdata->masterprob, SCIPgetBestSol(relaxdata->masterprob), 
            nmastervars, mastervars, mastervals) );

      SCIP_CALL( GCGrelaxTransformMastervalsToOrigvals(scip, mastervars, mastervals, nmastervars, 
            origvars, origvals, norigvars) );

      SCIP_CALL( SCIPsetSolVals(scip, newsol, norigvars, origvars, origvals) );

      //SCIP_CALL( SCIPprintSol(scip, newsol, NULL, FALSE) );

      SCIP_CALL( SCIPtrySolFree(scip, &newsol, TRUE, TRUE, TRUE, &stored) );
      assert(stored);

      SCIPdebugMessage("updated current best primal feasible solution!\n");
   }

   SCIPfreeBufferArray(scip, &origvals);
   SCIPfreeBufferArray(scip, &mastervals);

   return SCIP_OKAY;
}        


/** returns the number of fractional variables in the relaxator's current solution */
int GCGrelaxGetNBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_VAR** vars;
   int nvars;
   int v;
   int nbranchcands;
   SCIP_Real* vals;

   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* return if no current solution exists */
   if ( relaxdata->currentorigsol == NULL )
   {
      return -1;
   }

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   assert(vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, relaxdata->currentorigsol, nvars, vars, vals) );

   nbranchcands = 0;
   for ( v = 0; v < nvars; v++ )
   {
      if ( !SCIPisFeasIntegral(scip, vals[v]) )
      {
         nbranchcands++;
      }
   }

   SCIPfreeBufferArray(scip, &vals);

   return nbranchcands;
}        

/** returns the fractional variables in the relaxator's current solution */
SCIP_RETCODE GCGrelaxGetBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           branchcands,        /**< pointer to store the array of branching candidates, or NULL */
   SCIP_Real**           branchcandssol,     /**< pointer to store the array of candidate solution values, or NULL */
   SCIP_Real**           branchcandsfrac,    /**< pointer to store the array of candidate fractionalities, or NULL */
   int*                  nbranchcands,       /**< pointer to store the number of branching candidates, or NULL */
   int*                  npriobranchcands    /**< pointer to store the number of candidates with maximal priority, or NULL */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_VAR** vars;
   int nvars;
   int v;
   int ncands;
   SCIP_Real* vals;

   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* return if no current solution exists */
   if ( relaxdata->currentorigsol == NULL )
   {
      return -1;
   }

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   assert(vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, relaxdata->currentorigsol, nvars, vars, vals) );

   ncands = 0;
   for ( v = 0; v < nvars; v++ )
   {
      if ( !SCIPisFeasIntegral(scip, vals[v]) )
      {
         relaxdata->branchcands[*nbranchcands] = vars[v];
         relaxdata->branchcandssol[*nbranchcands] = vals[v];
         relaxdata->branchcandsfrac[*nbranchcands] = SCIPfrac(scip, vals[v]);
         ncands++;
      }
   }

   SCIPfreeBufferArray(scip, &vals);

   if ( nbranchcands != NULL )
      *nbranchcands = ncands;
   if ( npriobranchcands != NULL )
      *npriobranchcands = ncands;
   if ( branchcandsfrac != NULL )
      *branchcandsfrac = relaxdata->branchcandsfrac;
   if ( branchcandssol != NULL )
      *branchcandssol = relaxdata->branchcandssol;
   if ( branchcands != NULL )
      *branchcands = relaxdata->branchcands;

   return SCIP_OKAY;
}        

/** returns solution value of the given variable in the relaxator's current solution */
SCIP_Real GCGrelaxGetVarSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* return if no current solution exists */
   if ( relaxdata->currentorigsol == NULL )
   {
      return 0.0;
   }

   return SCIPgetSolVal(scip, relaxdata->currentorigsol, var);
}        

/** returns solution value of the given variable in the relaxator's current solution */
SCIP_RETCODE GCGrelaxLinkSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_VAR** vars;
   int nvars;
   SCIP_Real* vals;

   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* return if no current solution exists */
   if ( relaxdata->currentorigsol == NULL )
   {
      SCIP_CALL( SCIPclearSol(scip, sol) );
      printf("sol cleared!\n");
   }

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   assert(vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, relaxdata->currentorigsol, nvars, vars, vals) );
   SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, vals) );

   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}

/* transforms given values of the given original variables into values of the given master variables */
void GCGrelaxTransformOrigvalsToMastervals(
   SCIP*                 scip,               /** SCIP data structure */
   SCIP_VAR**            origvars,           /** array with (subset of the) original variables */
   SCIP_Real*            origvals,           /** array with values for the given original variables */
   int                   norigvars,          /** number of given original variables */
   SCIP_VAR**            mastervars,         /** array of (all present) master variables */
   SCIP_Real*            mastervals,         /** return value: values of the master variables */
   int                   nmastervars         /** number of master variables */
   )
{
   SCIP_VARDATA* vardata;
   SCIP_VARDATA* vardata2;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(origvars != NULL);
   assert(origvals != NULL);
   assert(mastervars != NULL);
   assert(mastervals != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);


   for ( i = 0; i < nmastervars; i++ )
   {
      mastervals[i] = 0.0;
   }

   for ( i = 0; i < norigvars; i++ )
   {
      vardata = SCIPvarGetData(origvars[i]);

      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->data.origvardata.nmastervars >= 0);
      assert(vardata->data.origvardata.mastervars != NULL);
      assert(vardata->data.origvardata.mastervals != NULL);
      assert(vardata->data.origvardata.nmastervars == 1 || vardata->blocknr != -1);

      if ( vardata->blocknr == -1 )
      {
         //printf("%f, ", vardata->data.origvardata.mastervals[j]);
         for ( k = 0; k < nmastervars; k++ )
         {
            assert(!SCIPvarIsTransformedOrigvar(mastervars[k]));
            if ( mastervars[k] == vardata->data.origvardata.mastervars[0])
            {
               assert(!SCIPvarIsTransformedOrigvar(vardata->data.origvardata.mastervars[0]));
               mastervals[k] += (vardata->data.origvardata.mastervals[0] * origvals[i]);
               break;
            }
         }
         assert(k < nmastervars);
         //printf("\n");
      }
      else
      {
         vardata2 = SCIPvarGetData(vardata->data.origvardata.pricingvar);
         assert(vardata2->vartype == GCG_VARTYPE_PRICING);
         assert(vardata2->data.pricingvardata.norigvars >= 0);
         assert(vardata2->data.pricingvardata.origvars != NULL);
         assert(vardata2->data.pricingvardata.origvars[0] != NULL);

         vardata2 = SCIPvarGetData(vardata2->data.pricingvardata.origvars[0]);
         assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);
         assert(vardata2->data.origvardata.nmastervars >= 0);
         assert(vardata2->data.origvardata.mastervars != NULL);
         assert(vardata2->data.origvardata.mastervals != NULL);

         for ( j = 0; j < vardata2->data.origvardata.nmastervars; j++ )
         {
            //printf("%f, ", vardata->data.origvardata.mastervals[j]);
            for ( k = 0; k < nmastervars; k++ )
            {
               if ( mastervars[k] == vardata2->data.origvardata.mastervars[j] )
               {
                  mastervals[k] += (vardata2->data.origvardata.mastervals[j] * origvals[i]);
                  break;
               }
            }
            assert(k < nmastervars);
            //printf("\n");
         }
      }

   }
}       

/* transforms given values of the given master variables into values of the given original variables */
SCIP_RETCODE GCGrelaxTransformMastervalsToOrigvals(
   SCIP*                 scip,               /** SCIP data structure */
   SCIP_VAR**            mastervars,         /** array of (subset of the) master variables */
   SCIP_Real*            mastervals,         /** values of the given master variables */
   int                   nmastervars,        /** number of master variables */
   SCIP_VAR**            origvars,           /** array with (all present) original variables */
   SCIP_Real*            origvals,           /** return value: array with values of the original variables */
   int                   norigvars           /** number of given original variables */
   )
{
   SCIP_VARDATA* vardata;
   SCIP_VARDATA* vardata2;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_Real* blockvalue;
   SCIP_Real increaseval;
   SCIP_Real* copymastervals;
   int* blocknr;
   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(origvars != NULL);
   assert(origvals != NULL);
   assert(mastervars != NULL);
   assert(mastervals != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* set all values in the given array to 0 */
   for ( i = 0; i < norigvars; i++ )
   {
      assert(!SCIPvarIsNegated(origvars[i]));
      origvals[i] = 0.0;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &blockvalue, relaxdata->npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocknr, relaxdata->npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &copymastervals, nmastervars) );

   for ( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      blockvalue[i] = 0.0;
      blocknr[i] = 0;
   }

   for ( i = 0; i < nmastervars; i++ )
   {
      copymastervals[i] = mastervals[i];
   }


   /* loop over all given master variables */
   for ( i = 0; i < nmastervars; i++ )
   {
      /* first of all, handle the variables with integral values */
      while ( SCIPisFeasGE(scip, copymastervals[i], 1) )
      {
         vardata = SCIPvarGetData(mastervars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_MASTER);
         assert(vardata->data.mastervardata.norigvars > 0);
         assert(vardata->data.mastervardata.origvars != NULL);
         assert(vardata->data.mastervardata.origvals != NULL);

         if ( vardata->blocknr == -1 )
         {
            assert(vardata->data.mastervardata.norigvars == 2);
            assert(vardata->data.mastervardata.origvals[0] == 1.0);
            assert(vardata->data.mastervardata.origvals[1] == 0.0);

            /* search the variable in the given array and increase the corresponding value */
            for ( k = 0; k < norigvars; k++ )
            {
               if ( origvars[k] == vardata->data.mastervardata.origvars[0] )
               {
                  origvals[k] += (vardata->data.mastervardata.origvals[0] * copymastervals[i]);
                  copymastervals[i] = 0.0;
                  break;
               }
            }
            assert(k < norigvars);
         }
         else
         {
            /* loop over all original variables contained in the current master variable */
            for ( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
            {
               /* get the right original variable */
               vardata2 = SCIPvarGetData(vardata->data.mastervardata.origvars[j]);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);
               assert(vardata2->data.origvardata.pricingvar != NULL);
               vardata2 = SCIPvarGetData(vardata2->data.origvardata.pricingvar);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_PRICING);
               assert(vardata2->data.pricingvardata.norigvars > blocknr[vardata->blocknr]);


               /* search the variable in the given array and increase the corresponding value */
               for ( k = 0; k < norigvars; k++ )
               {
                  if ( origvars[k] == vardata2->data.pricingvardata.origvars[blocknr[vardata->blocknr]] )
                  {
                     origvals[k] += (vardata->data.mastervardata.origvals[j]);
                     break;
                  }
               }
               assert(k < norigvars);
            }

            copymastervals[i] = copymastervals[i] - 1.0;
            blocknr[vardata->blocknr]++;
         }
      }
   }

   /* loop over all given master variables */
   for ( i = 0; i < nmastervars; i++ )
   {
      if ( SCIPisFeasZero(scip, copymastervals[i]) || SCIPisFeasIntegral(scip, copymastervals[i]) )
      {
         continue;
      }

      assert(SCIPisFeasGE(scip, copymastervals[i], 0.0) && SCIPisFeasLT(scip, copymastervals[i], 1.0));

      while ( SCIPisFeasPositive(scip, copymastervals[i]) )
      {
         vardata = SCIPvarGetData(mastervars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_MASTER);
         assert(vardata->data.mastervardata.norigvars > 0);
         assert(vardata->data.mastervardata.origvars != NULL);
         assert(vardata->data.mastervardata.origvals != NULL);
         
         if ( vardata->blocknr == -1 )
         {
            assert(vardata->data.mastervardata.norigvars == 2);
            assert(vardata->data.mastervardata.origvals[0] == 1.0);
            assert(vardata->data.mastervardata.origvals[1] == 0.0);
            
            /* search the variable in the given array and increase the corresponding value */
            for ( k = 0; k < norigvars; k++ )
            {
               if ( origvars[k] == vardata->data.mastervardata.origvars[0] )
               {
                  origvals[k] += (vardata->data.mastervardata.origvals[0] * copymastervals[i]);
                  copymastervals[i] = 0.0;
                  break;
               }
            }
            assert(k < norigvars);
         }
         else
         {
            increaseval = MIN(copymastervals[i], 1.0 - blockvalue[vardata->blocknr]);
            /* loop over all original variables contained in the current master variable */
            for ( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
            {
               /* get the right original variable */
               vardata2 = SCIPvarGetData(vardata->data.mastervardata.origvars[j]);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);
               assert(vardata2->data.origvardata.pricingvar != NULL);
               vardata2 = SCIPvarGetData(vardata2->data.origvardata.pricingvar);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_PRICING);
               assert(vardata2->data.pricingvardata.norigvars > blocknr[vardata->blocknr]);


               /* search the variable in the given array and increase the corresponding value */
               for ( k = 0; k < norigvars; k++ )
               {
                  if ( origvars[k] == vardata2->data.pricingvardata.origvars[blocknr[vardata->blocknr]] )
                  {
                     origvals[k] += (vardata->data.mastervardata.origvals[j] * increaseval );
                     break;
                  }
               }
               assert(k < norigvars);
            }

            copymastervals[i] = copymastervals[i] - increaseval;
            if ( SCIPisFeasZero(scip, copymastervals[i]) )
            {
               copymastervals[i] = 0.0;
            }
            blockvalue[vardata->blocknr] += increaseval;
            assert(SCIPisFeasLE(scip, blockvalue[vardata->blocknr], 1.0));
            if ( SCIPisFeasEQ(scip, blockvalue[vardata->blocknr], 1.0) )
            {
               blockvalue[vardata->blocknr] = 0.0;
               blocknr[vardata->blocknr]++;
            }

         }
      }

   }

   SCIPfreeBufferArray(scip, &copymastervals);
   SCIPfreeBufferArray(scip, &blocknr);
   SCIPfreeBufferArray(scip, &blockvalue);

   for ( i = 0; i < norigvars; i++ )
   {
      if ( SCIPisFeasIntegral(scip, origvals[i]) )
      {
         origvals[i] = SCIPfeasCeil(scip, origvals[i]);
      }
   }

   return SCIP_OKAY;

}       
