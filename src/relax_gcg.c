/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"
//#define SCIP_DEBUG
//#define CHECKCONSISTENCY
/**@file    relax_gcg.c
 * @ingroup RELAXATORS
 * @brief   gcg relaxator
 * @author  Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "struct_vardata.h"
#include "struct_branchgcg.h"
#include "type_branchgcg.h"

#include "relax_gcg.h"
#include "gcgplugins.h"
#include "pricer_gcg.h"
#include "masterplugins.h"
#include "pricingplugins.h"
#include "nodesel_master.h"



#define RELAX_NAME             "gcg"
#define RELAX_DESC             "relaxator for gcg project representing the master lp"
#define RELAX_PRIORITY         -1
#define RELAX_FREQ             1

#define STARTMAXMASTERVARS 10
#define DEFAULT_DISCRETIZATION TRUE
#define DEFAULT_MERGEIDENTICALBLOCS TRUE
#define DEFAULT_DISPINFOS FALSE


/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   /* problems and convexity constraints */
   SCIP*            masterprob;          /* the master problem */
   SCIP**           pricingprobs;        /* the array of pricing problems */
   int              npricingprobs;       /* the number of pricing problems */
   int              nrelpricingprobs;    /* the number of relevantpricing problems */
   int*             blockrepresentative; /* number of the pricing problem, that represents the i-th problem */
   int*             nblocksidentical;    /* number of pricing blocks represented by the i-th pricing problem */
   SCIP_CONS**      convconss;           /* array of convexity constraints, one for each block */

   /* hashmaps for transformation */
   SCIP_HASHMAP**   hashorig2pricingvar; /* hashmap mapping original variables to corresponding 
                                          * pricing variables */
   SCIP_HASHMAP*    hashorig2origvar;    /* hashmap mapping original variables to themselves */

   /* constraint data */
   SCIP_CONS**      masterconss;         /* array of constraints in the master problem */
   SCIP_CONS**      origmasterconss;     /* array of constraints in the original problem that belong to the 
                                          * master problem */
   SCIP_CONS**      linearmasterconss;   /* array of linear constraints equivalent to the cons in 
                                          * the original problem that belong to the master problem */
   int              maxmasterconss;      /* length of the array mastercons */
   int              nmasterconss;        /* number of constraints saved in mastercons */
   
   SCIP_SOL*        currentorigsol;      /* current lp solution transformed into the original space */
   SCIP_Longint     lastmasterlpiters;   /* number of lp iterations when currentorigsol was updated the last time */
   SCIP_SOL*        lastmastersol;       /* last feasible master solution that was added to the original problem */
   SCIP_CONS**      markedmasterconss;   /* array of conss that are marked to be in the master */
   int              nmarkedmasterconss;  /* number of elements in array of conss that are marked to be in the master */
   SCIP_Longint     lastsolvednodenr;    /* node number of the node that was solved at the last call of the relaxator */

   /* branchrule data */
   GCG_BRANCHRULE** branchrules;         /* branching rules registered in the relaxator */
   int              nbranchrules;        /* number of branching rules registered in the relaxator */
   
   /* parameter data */
   SCIP_Bool        discretization;      /* TRUE: use discretization approach; FALSE: use convexification approach */
   SCIP_Bool        mergeidenticalblocks;/* should identical blocks be merged (only for discretization approach)? */
   SCIP_Bool        masterissetpart;     /* is the master a set partitioning problem? */
   SCIP_Bool        masterissetcover;    /* is the master a set covering problem? */
   SCIP_Bool        dispinfos;           /* should additional information be displayed? */

   /* data for probing */
   SCIP_Bool        masterinprobing;     /* is the master problem in probing mode? */
   SCIP_SOL*        storedorigsol;       /* orig solution that was stored from before the probing */
};



/*
 * Vardata methods
 */

static
SCIP_DECL_VARDELORIG(gcgvardelorig)
{
   if( (*vardata)->vartype == GCG_VARTYPE_ORIGINAL )
   {
      SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.mastervars));
      SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.mastervals));
      if((*vardata)->data.origvardata.ncoefs > 0)
      {
         assert((*vardata)->data.origvardata.coefs != NULL);
         assert((*vardata)->data.origvardata.linkconss != NULL);
         SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.coefs));
         SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.linkconss));
      }
      
   }
   if( (*vardata)->vartype == GCG_VARTYPE_PRICING )
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

/* ensures size of masterconss array */
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

   if( relaxdata->maxmasterconss < size )
   {
      relaxdata->maxmasterconss = MAX(relaxdata->maxmasterconss + 5, size);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(relaxdata->masterconss), relaxdata->maxmasterconss) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(relaxdata->origmasterconss), relaxdata->maxmasterconss) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(relaxdata->linearmasterconss), relaxdata->maxmasterconss) );
   }
   assert(relaxdata->maxmasterconss >= size);

   return SCIP_OKAY;
}

/* ensures size of branchrules array: enlarges the array by 1 */
static
SCIP_RETCODE ensureSizeBranchrules(
   SCIP*                 scip,
   SCIP_RELAXDATA*       relaxdata
   )
{
   assert(scip != NULL);
   assert(relaxdata != NULL);
   assert((relaxdata->branchrules == NULL) == (relaxdata->nbranchrules == 0));   

   if( relaxdata->nbranchrules == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->branchrules), 1) );
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(relaxdata->branchrules), relaxdata->nbranchrules+1) );
   }

   return SCIP_OKAY;
}

/* checks whether two arrays of SCIP_Real's are identical */
static 
SCIP_Bool realArraysAreEqual(
   SCIP_Real*            array1,
   int                   array1length,
   SCIP_Real*            array2,
   int                   array2length
   )
{
   int i;

   assert(array1 != NULL || array1length == 0);
   assert(array2 != NULL || array2length == 0); 

   if( array1length != array2length )
      return FALSE;

   for( i = 0; i < array1length; i++ )
   {
      if( array1[i] != array2[i] )
         return FALSE;
   }

   return TRUE;
}

/* checks whether two pricingproblems represent identical blocks */
static
SCIP_RETCODE pricingprobsAreIdentical(
   SCIP_RELAXDATA*       relaxdata,          /**< the relaxator's data */
   int                   probnr1,            /**< number of the first pricingproblem */
   int                   probnr2,            /**< number of the second pricingproblem */
   SCIP_HASHMAP*         varmap,             /**< hashmap mapping the variables of the second pricing problem
                                              *   to those of the first pricing problem */
   SCIP_Bool*            identical           /**< return value: are blocks identical */
   )
{
   SCIP* scip1;
   SCIP* scip2;
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

   assert(relaxdata != NULL);
   assert(0 <= probnr1 && probnr1 < relaxdata->npricingprobs);
   assert(0 <= probnr2 && probnr2 < relaxdata->npricingprobs);
   assert(varmap != NULL);
   assert(identical != NULL);

   scip1 = relaxdata->pricingprobs[probnr1];
   scip2 = relaxdata->pricingprobs[probnr2];
   assert(scip1 != NULL);
   assert(scip2 != NULL);

   *identical = FALSE;

   SCIPdebugMessage("check block %d and block %d for identity...\n", probnr1, probnr2);

   if( SCIPgetNVars(scip1) != SCIPgetNVars(scip2))
   {
      SCIPdebugMessage("--> number of variables differs!\n");
      return SCIP_OKAY;
   }
   if( SCIPgetNConss(scip1) != SCIPgetNConss(scip1))
   {
      SCIPdebugMessage("--> number of constraints differs!\n");
      return SCIP_OKAY;
   }
   /* get variables */
   SCIP_CALL( SCIPgetVarsData(scip1, &vars1, &nvars1, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(scip2, &vars2, &nvars2, NULL, NULL, NULL, NULL) );
   
   for( i = 0; i < nvars1; i++ )
   {
      if( SCIPvarGetObj(vars1[i]) != SCIPvarGetObj(vars2[i]) )
      {
         SCIPdebugMessage("--> obj differs for var %s and var %s!\n", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
         return SCIP_OKAY;
      }
      if( SCIPvarGetLbOriginal(vars1[i]) != SCIPvarGetLbOriginal(vars2[i]) )
      {
         SCIPdebugMessage("--> lb differs for var %s and var %s!\n", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
         return SCIP_OKAY;
      }
      if( SCIPvarGetUbOriginal(vars1[i]) != SCIPvarGetUbOriginal(vars2[i]) )
      {
         SCIPdebugMessage("--> ub differs for var %s and var %s!\n", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
         return SCIP_OKAY;
      }
      if( SCIPvarGetType(vars1[i]) != SCIPvarGetType(vars2[i]) )
      {
         SCIPdebugMessage("--> type differs for var %s and var %s!\n", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
         return SCIP_OKAY;
      }
      
      vardata1 = SCIPvarGetData(vars1[i]);
      vardata2 = SCIPvarGetData(vars2[i]);
      assert(vardata1 != NULL && vardata2 != NULL);
      assert(vardata1->vartype == GCG_VARTYPE_PRICING);
      assert(vardata2->vartype == GCG_VARTYPE_PRICING);
      assert(vardata1->data.pricingvardata.origvars != NULL);
      assert(vardata2->data.pricingvardata.origvars != NULL);

      if( SCIPvarGetObj(vardata1->data.pricingvardata.origvars[0]) 
         != SCIPvarGetObj(vardata2->data.pricingvardata.origvars[0]) )
      {
         SCIPdebugMessage("--> orig obj differs for var %s and var %s!\n", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
         return SCIP_OKAY;
      }

      vardata1 = SCIPvarGetData(vardata1->data.pricingvardata.origvars[0]);
      vardata2 = SCIPvarGetData(vardata2->data.pricingvardata.origvars[0]);
      assert(vardata1 != NULL && vardata2 != NULL);
      assert(vardata1->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);

      if( !realArraysAreEqual(vardata1->data.origvardata.coefs, vardata1->data.origvardata.ncoefs,
            vardata2->data.origvardata.coefs, vardata2->data.origvardata.ncoefs) )
      {
         SCIPdebugMessage("--> coefs differ for var %s and var %s!\n", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
         return SCIP_OKAY;
      }
      SCIP_CALL( SCIPhashmapInsert(varmap, (void*) vars1[i], (void*) vars2[i]) );

   }

   /* check whether the conss are the same */
   conss1 = SCIPgetConss(scip1);
   conss2 = SCIPgetConss(scip2);
   nconss = SCIPgetNConss(scip1);
   assert(nconss == SCIPgetNConss(scip2));
   for( i = 0; i < nconss; i++ )
   {
      if( SCIPgetNVarsLinear(scip1, conss1[i]) != SCIPgetNVarsLinear(scip2, conss2[i]) )
      {
         SCIPdebugMessage("--> nvars differs for cons %s and cons %s!\n", SCIPconsGetName(conss1[i]), SCIPconsGetName(conss2[i]));
         return SCIP_OKAY;
      }
      if( SCIPgetLhsLinear(scip1, conss1[i]) != SCIPgetLhsLinear(scip2, conss2[i]) )
      {
         SCIPdebugMessage("--> lhs differs for cons %s and cons %s!\n", SCIPconsGetName(conss1[i]), SCIPconsGetName(conss2[i]));
         return SCIP_OKAY;
      }
      if( SCIPgetRhsLinear(scip1, conss1[i]) != SCIPgetRhsLinear(scip2, conss2[i]) )
      {
         SCIPdebugMessage("--> rhs differs for cons %s and cons %s!\n", SCIPconsGetName(conss1[i]), SCIPconsGetName(conss2[i]));
         return SCIP_OKAY;
      }
      if( !realArraysAreEqual(SCIPgetValsLinear(scip1, conss1[i]), SCIPgetNVarsLinear(scip1, conss1[i]),
            SCIPgetValsLinear(scip2, conss2[i]), SCIPgetNVarsLinear(scip2, conss2[i])) )
      {
         SCIPdebugMessage("--> coefs differ for cons %s and cons %s!\n", SCIPconsGetName(conss1[i]), SCIPconsGetName(conss2[i]));
         return SCIP_OKAY;
      }
      vars1 = SCIPgetVarsLinear(scip1, conss1[i]);
      vars2 = SCIPgetVarsLinear(scip2, conss2[i]);
      for( j = 0; j < SCIPgetNVarsLinear(scip1, conss1[i]); j++ )
      {
         if( (SCIP_VAR*) SCIPhashmapGetImage(varmap, (void*) vars1[j]) != vars2[j] )
         {
            SCIPdebugMessage("--> vars differ for cons %s and cons %s!\n", SCIPconsGetName(conss1[i]), SCIPconsGetName(conss2[i]));
            return SCIP_OKAY;
         }
      }

   }

   SCIPdebugMessage("--> blocks are identical!\n");

   *identical = TRUE;
   return SCIP_OKAY;
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
   SCIP_Bool identical;

   int i;
   int j;
   int k;

   int nrelevant;
   

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   for( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      relaxdata->blockrepresentative[i] = i;
      relaxdata->nblocksidentical[i] = 1;
   }
   relaxdata->nrelpricingprobs = relaxdata->npricingprobs;
   nrelevant = 0;

   if( !relaxdata->discretization || !relaxdata->mergeidenticalblocks )
      return SCIP_OKAY;

   for( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      for( j = 0; j < i && relaxdata->blockrepresentative[i] == i; j++ )
      {
         if( relaxdata->blockrepresentative[j] != j )
            continue;

         SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), 5 * SCIPgetNVars(relaxdata->pricingprobs[i])+1) ); /* TODO: +1 to deal with empty subproblems */
         SCIP_CALL( pricingprobsAreIdentical(relaxdata, i, j, varmap, &identical) );
         
         if( identical )
         {
            SCIPdebugMessage("Block %d is identical to block %d!\n", i, j);
            
            /* block i will be represented by block j */
            relaxdata->blockrepresentative[i] = j;
            relaxdata->nblocksidentical[i] = 0;
            relaxdata->nblocksidentical[j]++;
            /* save variables in pricing problem variable's vardata */
            vars = SCIPgetVars(relaxdata->pricingprobs[i]);
            nvars = SCIPgetNVars(relaxdata->pricingprobs[i]);
            for( k = 0; k < nvars; k++ )
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
               if( vardata->data.pricingvardata.norigvars >= 2 )
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
      if( relaxdata->blockrepresentative[i] == i )
      {
         SCIPdebugMessage("Block %d is relevant!\n", i);
         nrelevant++;
      }
   }

   printf("Matrix has %d blocks, %d %s relevant!\n", relaxdata->npricingprobs, nrelevant, 
      (nrelevant == 1 ? "is" : "are"));

   relaxdata->nrelpricingprobs = nrelevant;

   return SCIP_OKAY;
}

/** checks whether a constrains belongs to a block */
static
SCIP_Bool consIsInBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< hashmap mapping variables of original to variables of pricing problem */
   SCIP_CONS*            cons                /**< the constraint that should be checked for correspondence to the block */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(varmap != NULL);
   assert(cons != NULL);


   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "linear") == 0 )
   {
      vars = SCIPgetVarsLinear(scip, cons);
      nvars = SCIPgetNVarsLinear(scip, cons);
   }
   else if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "setppc") == 0 )
   {
      vars = SCIPgetVarsSetppc(scip, cons);
      nvars = SCIPgetNVarsSetppc(scip, cons);
   }
   else if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "knapsack") == 0 )
   {
      vars = SCIPgetVarsKnapsack(scip, cons);
      nvars = SCIPgetNVarsKnapsack(scip, cons);
   }
   else if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "logicor") == 0 )
   {
      vars = SCIPgetVarsLogicor(scip, cons);
      nvars = SCIPgetNVarsLogicor(scip, cons);
   }   
   else if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "varbound") == 0 )
   {
      /* check whether bounded variable is contained in block */
      var = SCIPgetVarVarbound(scip, cons);
      if( SCIPvarIsNegated(var) )
         var = SCIPvarGetNegationVar(var);
      if( !SCIPhashmapExists(varmap, (void*) var) )
         return FALSE;

      /* check whether bounding variable is contained in block */
      var = SCIPgetVbdvarVarbound(scip, cons);
      if( SCIPvarIsNegated(var) )
         var = SCIPvarGetNegationVar(var);
      if( !SCIPhashmapExists(varmap, (void*) var) )
         return FALSE;

      /* both variables are in the block, return TRUE */
      return TRUE;
   }
   else
   {
      printf("constraint %s of unknown type <%s>, copy failed!\n", SCIPconsGetName(cons), SCIPconshdlrGetName(SCIPconsGetHdlr(cons)));
   }

   for( i = 0; i < nvars; i++ )
   {
      var = vars[i];
      if( SCIPvarIsNegated(var) )
      {
         var = SCIPvarGetNegationVar(var);
      }
      if( !SCIPhashmapExists(varmap, (void*) var) )
      {
         return FALSE;
      }
   }

   return TRUE;

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
   SCIP_Bool marked;
   int i;
   int c;
   int v;
   int b;

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   SCIPdebugMessage("Creating master problem...\n");

   /* initialize relaxator data */
   relaxdata->maxmasterconss = 5;
   relaxdata->nmasterconss = 0;

   /* arrays of constraints belonging to the master problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->masterconss), relaxdata->maxmasterconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->origmasterconss), relaxdata->maxmasterconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->linearmasterconss), relaxdata->maxmasterconss) );

   /* create the problem in the master scip instance */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "master_%s", SCIPgetProbName(scip));

   SCIP_CALL( SCIPcreateProb(relaxdata->masterprob, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* activate the pricer */
   SCIP_CALL( SCIPactivatePricer(relaxdata->masterprob, SCIPfindPricer(relaxdata->masterprob, "gcg")) );
   //SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "presolving/probing/maxrounds", 0) );

   SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "pricing/maxvars", INT_MAX) );
   SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "pricing/maxvarsroot", INT_MAX) );
   SCIP_CALL( SCIPsetRealParam(relaxdata->masterprob, "pricing/abortfac", 1) );

   /* ----- initialize the pricing problems ----- */
   npricingprobs = relaxdata->npricingprobs;
   assert(npricingprobs >= 0);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->pricingprobs), npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->blockrepresentative), npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->nblocksidentical), npricingprobs) );

   /* array for saving convexity constraints belonging to one of the pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->convconss), npricingprobs) );

   //SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "display/verblevel", SCIP_VERBLEVEL_FULL) );

   /* create the pricing problem */
   for( i = 0; i < npricingprobs; i++ )
   {
      relaxdata->convconss[i] = NULL;

      /* initializing the scip data structure for the pricing problem */  
      SCIP_CALL( SCIPcreate(&(relaxdata->pricingprobs[i])) );
      SCIP_CALL( SCIPincludeDefaultPlugins(relaxdata->pricingprobs[i]) );
 
      /* disable conflict analysis */
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "conflict/useprop", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "conflict/useinflp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "conflict/useboundlp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "conflict/usesb", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "conflict/usepseudo", FALSE) );

      /* reduce the effort spent for hash tables */
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "misc/usevartable", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "misc/useconstable", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "misc/usesmalltables", TRUE) );

#if 0      
      /* diable some buggy heuristics */
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "heuristics/oneopt/freq", -1) );
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "heuristics/zirounding/freq", -1) );
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "separating/rapidlearning/freq", -1) );
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "heuristics/rens/freq", -1) );
#endif

      /* disable expensive presolving */
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "presolving/probing/maxrounds", 0) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/linear/presolpairwise", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/setppc/presolpairwise", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/logicor/presolpairwise", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/linear/presolusehashing", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/setppc/presolusehashing", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/logicor/presolusehashing", FALSE) );

      /* disable output to console */
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "display/verblevel", SCIP_VERBLEVEL_NONE) );

      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "misc/catchctrlc", FALSE) );

      /* create the pricing submip */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricing_block_%d", i);
      SCIP_CALL( SCIPcreateProb(relaxdata->pricingprobs[i], name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   }

   /* create hashmaps for mapping from original to pricing variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->hashorig2pricingvar), npricingprobs) );
   for( i = 0; i < npricingprobs; i++ )
   {
      SCIP_CALL( SCIPhashmapCreate(&(relaxdata->hashorig2pricingvar[i]), 
            SCIPblkmem(scip), SCIPgetNVars(scip)) );
   }
   SCIP_CALL( SCIPhashmapCreate(&(relaxdata->hashorig2origvar), 
         SCIPblkmem(scip), 10*SCIPgetNVars(scip)) );

   /* create pricing variables and map them to the original variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for( v = 0; v < nvars; v++ )
   {
      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      if( vardata->blocknr != -1 )
      {
         assert(vardata->data.origvardata.pricingvar == NULL);
	 
	 SCIP_CALL( GCGrelaxCreatePricingVar(scip, vars[v]) );
	 assert(vardata->data.origvardata.pricingvar != NULL);

	 SCIP_CALL( SCIPhashmapInsert(relaxdata->hashorig2pricingvar[vardata->blocknr], 
	     (void*)(vars[v]), (void*)(vardata->data.origvardata.pricingvar)) );
	 SCIP_CALL( SCIPhashmapInsert(relaxdata->hashorig2origvar, 
	     (void*)(vars[v]), (void*)(vars[v])) );
      }
      else
      {
         assert(vardata->data.origvardata.pricingvar == NULL);
         SCIP_CALL( SCIPhashmapInsert(relaxdata->hashorig2origvar, 
               (void*)(vars[v]), (void*)(vars[v])) );
      }
   }

   /* ------- copy constraints of the original problem into master/pricing problems ------- */
   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);

   /* iterate over all constraint handlers */
   for( i = 0; i < nconshdlrs; i++ )
   {
      if( strcmp(SCIPconshdlrGetName(conshdlrs[i]), "origbranch") == 0)
      {
         continue;
      }

      /* if there are constraints managed by this constraint handler, iterate over these constraints */
      nactiveconss = SCIPconshdlrGetNConss(conshdlrs[i]);

      /* upgraded linear constraints that were copied before are added a second time as linear constraints in the original problem,
       * hence, we disregard the last constraints */
      if( strcmp(SCIPconshdlrGetName(conshdlrs[i]), "linear") == 0)
      {
         nactiveconss -= relaxdata->nmasterconss;
#ifndef NDEBUG
         conss = SCIPconshdlrGetConss(conshdlrs[i]);
         for( c = 0; c < relaxdata->nmasterconss; c++ )
         {
            assert(conss[nactiveconss+c] == relaxdata->linearmasterconss[c]);
         }
#endif
      }

      if( nactiveconss > 0 )
      {
         conss = SCIPconshdlrGetConss(conshdlrs[i]);

         /* copy conss array */
         SCIP_CALL( SCIPallocBufferArray(scip, &bufconss, nactiveconss) );
         for( c = 0; c < nactiveconss; c++ )
         {
            bufconss[c] = conss[c];
         }
         for( c = 0; c < nactiveconss; c++ )
         {
            marked = FALSE;
            success = FALSE;

            /* check whether the constraint is marked to be transfered to the master */
            if( relaxdata->markedmasterconss != NULL )
            {
               for( b = 0; b < relaxdata->nmarkedmasterconss; b++ )
               {
                  if( strcmp(SCIPconsGetName(relaxdata->markedmasterconss[b]), SCIPconsGetName(bufconss[c])) == 0 )
                  {
                     marked = TRUE;
                     break;
                  }
               }
            }
            /* if it is not marked, try to copy the constraints into one of the pricing blocks */
            if( !marked )
            {
               for( b = 0; b < npricingprobs && !success; b++ )
               {
                  /* check whether constraint belongs to this block */
                  if( consIsInBlock(scip, relaxdata->hashorig2pricingvar[b], bufconss[c]) )
                  {
                     /* copy the constraint */
                     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "p%d_%s", b, SCIPconsGetName(bufconss[c]));
                     SCIP_CALL( SCIPgetConsCopy(scip, relaxdata->pricingprobs[b], bufconss[c], &newcons, conshdlrs[i],
                           relaxdata->hashorig2pricingvar[b], NULL, name, 
                           TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, &success) );

                     /* constraint was successfully copied */
                     assert(success);

                     SCIP_CALL( SCIPaddCons(relaxdata->pricingprobs[b], newcons) );
                     
                     SCIP_CALL( SCIPreleaseCons(relaxdata->pricingprobs[b], &newcons) );
                  }
               } 
            }
            else
            {
               SCIPdebugMessage("cons %s forced to be in the master problem!\n", SCIPconsGetName(bufconss[c]));
            }
            /* constraint was marked to be in the master or could not be copied into one of the pricing blocks */
            if( !success )
            {
               newcons = NULL;

               assert(SCIPhashmapGetNEntries(relaxdata->hashorig2origvar) == SCIPgetNVars(scip));

               /* copy the constraint (dirty trick, we only need lhs and rhs, because variables are added later) */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "linear_%s", SCIPconsGetName(bufconss[c]));
               SCIP_CALL( SCIPgetConsCopy(scip, scip, bufconss[c], &newcons, conshdlrs[i],
                     relaxdata->hashorig2origvar, NULL, name,
                     FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, &success) );
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

   /* for original variables, save the coefficients in the master problem in their vardata */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for( v = 0; v < nvars; v++ )
   {
      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->data.origvardata.coefs == NULL);

      vardata->data.origvardata.ncoefs = 0;
   }

   /* save coefs in the vardata */   
   for( i = 0; i < relaxdata->nmasterconss; i++ )
   {
      vars = SCIPgetVarsLinear(scip, relaxdata->linearmasterconss[i]);
      nvars = SCIPgetNVarsLinear(scip, relaxdata->linearmasterconss[i]);
      vals = SCIPgetValsLinear(scip, relaxdata->linearmasterconss[i]);
      for( v = 0; v < nvars; v++ )
      {
         assert(!SCIPisZero(scip, vals[v]));
         vardata = SCIPvarGetData(vars[v]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         if( vardata->data.origvardata.ncoefs == 0 ) 
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.origvardata.coefs), 1) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.origvardata.linkconss), 1) );
         }
         else
         {
            SCIP_CALL( SCIPreallocMemoryArray(scip, &(vardata->data.origvardata.coefs), vardata->data.origvardata.ncoefs+1) );
            SCIP_CALL( SCIPreallocMemoryArray(scip, &(vardata->data.origvardata.linkconss), vardata->data.origvardata.ncoefs+1) );
         }
         assert(vardata->data.origvardata.coefs != NULL);
         assert(vardata->data.origvardata.linkconss != NULL);
         vardata->data.origvardata.coefs[vardata->data.origvardata.ncoefs] = vals[v];
         vardata->data.origvardata.linkconss[vardata->data.origvardata.ncoefs] = relaxdata->masterconss[i];
         vardata->data.origvardata.ncoefs++;
      }
   }

   /* check for identity of blocks */
   SCIP_CALL( checkIdenticalBlocks(scip, relax) );

   for( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      if( relaxdata->blockrepresentative[i] != i )
         continue;

      /* create the corresponding convexity constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conv_block_%d", i);
      SCIP_CALL( SCIPcreateConsLinear(relaxdata->masterprob, &(relaxdata->convconss[i]), name, 0, NULL, NULL, 
            relaxdata->nblocksidentical[i], relaxdata->nblocksidentical[i], 
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(relaxdata->masterprob, relaxdata->convconss[i]) );
   }

   /* set integral obj status in the extended problem, if possible */
   if( SCIPisObjIntegral(scip) )
      SCIP_CALL( SCIPsetObjIntegral(relaxdata->masterprob) );

   /* display statistics */
   for( i = 0; i < relaxdata->npricingprobs && relaxdata->dispinfos; i++ )
   {
      int nbin;
      int nint;
      int nimpl;
      int ncont;

      if( relaxdata->blockrepresentative[i] != i )
         continue;

      SCIP_CALL( SCIPgetVarsData(relaxdata->pricingprobs[i], NULL, NULL, &nbin, &nint, &nimpl, &ncont) );

      printf("pricing problem %d: %d conss, %d vars (%d bins, %d ints, %d impls and %d cont)\n", i, 
         SCIPgetNConss(relaxdata->pricingprobs[i]), SCIPgetNVars(relaxdata->pricingprobs[i]), nbin, nint, nimpl, ncont);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricingprob_%d.lp", i);

      SCIP_CALL( SCIPwriteOrigProblem(relaxdata->pricingprobs[i], name, NULL, FALSE) );
   }

   return SCIP_OKAY;
}

#ifdef CHECKCONSISTENCY
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
   int nactivemasterconss;
   int i;
   int j;
   int k;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);

   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);
   
   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);
   assert(SCIPgetStage(masterprob) == SCIP_STAGE_TRANSFORMED || SCIPgetStage(masterprob) == SCIP_STAGE_SOLVING
      || SCIPgetStage(masterprob) == SCIP_STAGE_SOLVED);

   if( SCIPgetStage(masterprob) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   GCGconsOrigbranchCheckConsistency(scip);
   GCGconsMasterbranchCheckConsistency(masterprob);

   /* check variables, constraints and coefficients */

   nactivemasterconss = 0;
   for( i = 0; i < relaxdata->nmasterconss; i++ )
   {
      if( SCIPconsIsActive(relaxdata->masterconss[i]) )
         nactivemasterconss++;
   }
   assert(SCIPgetNActiveConss(masterprob) == nactivemasterconss + relaxdata->nrelpricingprobs 
      + GCGconsMasterbranchGetNStackelements(masterprob));

   SCIPdebugMessage("consistency checked: all ok!\n");

   return SCIP_OKAY;
}
#endif



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

   /* free master problem */
   if( relaxdata->masterprob != NULL )
      SCIP_CALL( SCIPfree(&(relaxdata->masterprob)) );
   
   SCIPfreeMemory(scip, &relaxdata);

   return SCIP_OKAY;
}



/** initialization method of relaxator (called after problem was transformed) */
static
SCIP_DECL_RELAXINIT(relaxInitGcg)
{  
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

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

   /* free array for branchrules*/
   if( relaxdata->nbranchrules > 0 )
   {
      for( i = 0; i < relaxdata->nbranchrules; i++ )
      {
         SCIPfreeMemory(scip, &(relaxdata->branchrules[i]));
      }
      SCIPfreeMemoryArray(scip, &(relaxdata->branchrules));
   }

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

   relaxdata->lastsolvednodenr = -1;

   masterprob = relaxdata->masterprob;
   assert(masterprob != NULL);

   SCIP_CALL( SCIPtransformProb(masterprob) );

   SCIP_CALL( SCIPtransformConss(masterprob, relaxdata->nmasterconss, 
         relaxdata->masterconss, relaxdata->masterconss) );

   for( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      if( relaxdata->convconss[i] != NULL) 
      {
         SCIP_CALL( SCIPtransformCons(masterprob, relaxdata->convconss[i], &(relaxdata->convconss[i])) );
      }
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolGcg)
{  
   SCIP_RELAXDATA* relaxdata;
   int i;

   assert(scip != NULL);
   assert(relax != NULL);
   
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* free hashmaps for mapping from original to pricing variables */
   if( relaxdata->hashorig2pricingvar != NULL )
   {
      for( i = 0; i < relaxdata->npricingprobs; i++ )
      {
         SCIPhashmapFree(&(relaxdata->hashorig2pricingvar[i]));
      }
      SCIPfreeMemoryArray(scip, &(relaxdata->hashorig2pricingvar));
      relaxdata->hashorig2pricingvar = NULL;
   }
   if( relaxdata->hashorig2origvar != NULL )
   {
      SCIPhashmapFree(&(relaxdata->hashorig2origvar));
      relaxdata->hashorig2origvar = NULL;
   }
   if( relaxdata->markedmasterconss != NULL )
   {
      SCIPfreeMemoryArray(scip, &(relaxdata->markedmasterconss));
      relaxdata->markedmasterconss = NULL;
   }

   /* free arrays for constraints */
   for( i = 0; i < relaxdata->nmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &relaxdata->origmasterconss[i]) );
   }
   for( i = 0; i < relaxdata->nmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &relaxdata->linearmasterconss[i]) );
   }
   for( i = 0; i < relaxdata->nmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(relaxdata->masterprob, &relaxdata->masterconss[i]) );
   }
   for( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      if( relaxdata->convconss[i] != NULL )
         SCIP_CALL( SCIPreleaseCons(relaxdata->masterprob, &relaxdata->convconss[i]) );
   }

   SCIPfreeMemoryArray(scip, &(relaxdata->origmasterconss));
   SCIPfreeMemoryArray(scip, &(relaxdata->linearmasterconss));
   SCIPfreeMemoryArray(scip, &(relaxdata->masterconss));
   SCIPfreeMemoryArray(scip, &(relaxdata->convconss));

   /* free master problem */
   SCIP_CALL( SCIPfree(&(relaxdata->masterprob)) );

   /* free pricing problems */
   for( i = relaxdata->npricingprobs - 1; i >= 0 ; i-- )
   {
      SCIP_CALL( SCIPfreeTransform(relaxdata->pricingprobs[i]) );
      SCIP_CALL( SCIPfree(&(relaxdata->pricingprobs[i])) );
   }
   SCIPfreeMemoryArray(scip, &(relaxdata->pricingprobs));
   SCIPfreeMemoryArray(scip, &(relaxdata->blockrepresentative));
   SCIPfreeMemoryArray(scip, &(relaxdata->nblocksidentical));

   /* free solution */
   if( relaxdata->currentorigsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &relaxdata->currentorigsol) );
   }

   return SCIP_OKAY;
}


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecGcg)
{  
   SCIP* masterprob;
   SCIP_RELAXDATA* relaxdata;
   SCIP_Bool cutoff;
   SCIP_Longint oldnnodes;
   SCIP_Real timelimit;
   SCIP_Bool feasible;

   assert(scip != NULL);
   assert(relax != NULL);
   assert(result != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   masterprob = relaxdata->masterprob;
   assert(masterprob != NULL);

   *result = SCIP_DIDNOTRUN;
   
   SCIPdebugMessage("solving node %lld's relaxation!\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   /* construct the LP in the original problem */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   assert(!cutoff);
   SCIP_CALL( SCIPflushLP(scip) );

   /* solve the next node in the master problem */
   SCIPdebugMessage("Solve master LP.\n");
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
   {
      SCIP_CALL( SCIPsetRealParam(relaxdata->masterprob, "limits/time", 
            timelimit - SCIPgetTotalTime(scip) + SCIPgetTotalTime(masterprob)) );
   }

   /* only solve the relaxation if it was not yet solved at the current node */
   if( SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) != relaxdata->lastsolvednodenr )
   {

      /* increase the node limit for the master problem by 1 */
      SCIP_CALL( SCIPgetLongintParam(masterprob, "limits/nodes", &oldnnodes) );
      SCIP_CALL( SCIPsetLongintParam(masterprob, "limits/nodes", 
            ( SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) ? 1 : oldnnodes+1)) );

      SCIP_CALL( SCIPsolve(masterprob) );

      /* set the lower bound pointer */
      if( SCIPgetStage(masterprob) == SCIP_STAGE_SOLVING )
         *lowerbound = SCIPgetLocalLowerbound(masterprob);
      else
      {
         assert(SCIPgetBestSol(masterprob) != NULL || SCIPgetStatus(masterprob) == SCIP_STATUS_INFEASIBLE);
         if( SCIPgetStatus(masterprob) == SCIP_STATUS_OPTIMAL )
            *lowerbound = SCIPgetSolOrigObj(masterprob, SCIPgetBestSol(masterprob));
         else if( SCIPgetStatus(masterprob) == SCIP_STATUS_INFEASIBLE )
            *lowerbound = SCIPinfinity(scip);
      }
      
      SCIPdebugMessage("Update lower bound (value = %"SCIP_REAL_FORMAT").\n", *lowerbound);
   }

   /* transform the current solution of the master problem to the original space and save it */
   SCIPdebugMessage("Update current sol.\n");
   SCIP_CALL( GCGrelaxUpdateCurrentSol(scip, &feasible) );

   if( GCGconsOrigbranchGetBranchrule(GCGconsOrigbranchGetActiveCons(scip)) != NULL 
      && SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) != relaxdata->lastsolvednodenr )
   {
      SCIP_CALL( GCGrelaxBranchMasterSolved(scip, GCGconsOrigbranchGetBranchrule(GCGconsOrigbranchGetActiveCons(scip)), 
            GCGconsOrigbranchGetBranchdata(GCGconsOrigbranchGetActiveCons(scip)), *lowerbound) );
   }

#ifdef CHECKCONSISTENCY
   SCIP_CALL( checkConsistency(scip) );
#endif

   /* update the number of the last solved node */
   relaxdata->lastsolvednodenr = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   *result = SCIP_SUCCESS;

   /* if the transferred master solution is feasible, the current node is solved to optimality and can be pruned */
   if( feasible )
   {
      *result = SCIP_CUTOFF;
      SCIPdebugMessage("solution was feasible, node can be cut off!");
   }

   return SCIP_OKAY;
}

#define relaxCopyGcg NULL



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
   relaxdata->markedmasterconss = NULL;

   relaxdata->nbranchrules = 0;
   relaxdata->branchrules = NULL;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelax(scip, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, relaxCopyGcg, relaxFreeGcg, relaxInitGcg, 
         relaxExitGcg, relaxInitsolGcg, relaxExitsolGcg, relaxExecGcg, relaxdata) );

   /* inform the main scip, that no LPs should be solved */
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );

   /* initialize the scip data structure for the master problem */  
   SCIP_CALL( SCIPcreate(&(relaxdata->masterprob)) );
   SCIP_CALL( SCIPincludePricerGcg(relaxdata->masterprob, scip) );
   SCIP_CALL( GCGincludeMasterPlugins(relaxdata->masterprob) );
 
   /* include masterbranch constraint handler */
   SCIP_CALL( SCIPincludeConshdlrMasterbranch(relaxdata->masterprob) );

   /* add gcg relaxator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/gcg/discretization",
         "should discretization (TRUE) or convexification (FALSE) approach be used?",
         &(relaxdata->discretization), FALSE, DEFAULT_DISCRETIZATION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/gcg/mergeidenticalblocks",
         "should identical blocks be merged (only for discretization approach)?",
         &(relaxdata->mergeidenticalblocks), FALSE, DEFAULT_MERGEIDENTICALBLOCS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/gcg/dispinfos",
         "should additional information about the blocks be displayed?",
         &(relaxdata->dispinfos), FALSE, DEFAULT_DISPINFOS, NULL, NULL) );

   return SCIP_OKAY;
}


/*
 * relaxator specific interface methods for coordination of branching rules
 */

/** includes a branching rule into the relaxator data */
SCIP_RETCODE GCGrelaxIncludeBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule for which callback methods are saved */
   GCG_DECL_BRANCHACTIVEMASTER     ((*branchactivemaster)),      /**<  activation method for branchrule */
   GCG_DECL_BRANCHDEACTIVEMASTER   ((*branchdeactivemaster)),    /**<  deactivation method for branchrule */
   GCG_DECL_BRANCHPROPMASTER       ((*branchpropmaster)),        /**<  propagation method for branchrule */
   GCG_DECL_BRANCHMASTERSOLVED     ((*branchmastersolved)),      /**<  master solved method for branchrule */
   GCG_DECL_BRANCHDATADELETE       ((*branchdatadelete))         /**<  branchdata deletion method for branchrule */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);
   assert(branchrule != NULL);
   
   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   SCIP_CALL( ensureSizeBranchrules(scip, relaxdata) );

   /* store callback functions */
   SCIP_CALL( SCIPallocMemory(scip, &(relaxdata->branchrules[relaxdata->nbranchrules])) );
   relaxdata->branchrules[relaxdata->nbranchrules]->branchrule = branchrule;
   relaxdata->branchrules[relaxdata->nbranchrules]->branchactivemaster = branchactivemaster;
   relaxdata->branchrules[relaxdata->nbranchrules]->branchdeactivemaster = branchdeactivemaster;
   relaxdata->branchrules[relaxdata->nbranchrules]->branchpropmaster = branchpropmaster;
   relaxdata->branchrules[relaxdata->nbranchrules]->branchmastersolved = branchmastersolved;
   relaxdata->branchrules[relaxdata->nbranchrules]->branchdatadelete = branchdatadelete;
   relaxdata->nbranchrules++;

   return SCIP_OKAY;
}

/** perform activation method of the given branchrule for the given branchdata */
SCIP_RETCODE GCGrelaxBranchActiveMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata          /**< data representing the branching decision */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;
   
   assert(scip != NULL);
   assert(branchrule != NULL);
   
   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         /* call activation method of branching rule */
         if( relaxdata->branchrules[i]->branchactivemaster != NULL )
            SCIP_CALL( relaxdata->branchrules[i]->branchactivemaster(relaxdata->masterprob, branchdata) );

         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

   return SCIP_OKAY;
}

/** perform deactivation method of the given branchrule for the given branchdata */
SCIP_RETCODE GCGrelaxBranchDeactiveMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata          /**< data representing the branching decision */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;
   
   assert(scip != NULL);
   assert(branchrule != NULL);
   
   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         /* call deactivation method of branching rule */
         if( relaxdata->branchrules[i]->branchdeactivemaster != NULL )
            SCIP_CALL( relaxdata->branchrules[i]->branchdeactivemaster(relaxdata->masterprob, branchdata) );

         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

   return SCIP_OKAY;
}

/** perform popagation method of the given branchrule for the given branchdata */
SCIP_RETCODE GCGrelaxBranchPropMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation call */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;
   
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(result != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   *result = SCIP_DIDNOTRUN;
   
   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         /* call propagation method of branching rule*/
         if( relaxdata->branchrules[i]->branchpropmaster != NULL )
            SCIP_CALL( relaxdata->branchrules[i]->branchpropmaster(relaxdata->masterprob, branchdata, result) );

         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

   return SCIP_OKAY;
}

/** frees branching data created by the given branchrule */
SCIP_RETCODE GCGrelaxBranchDataDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA**      branchdata          /**< data representing the branching decision */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;
   
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(branchdata != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         /* call branchrule data deletion method of the branching rule */
         if( relaxdata->branchrules[i]->branchdatadelete != NULL )
            SCIP_CALL( relaxdata->branchrules[i]->branchdatadelete(scip, branchdata) );
         else
         {
            if( *branchdata != NULL )
            {
               SCIPfreeMemory(scip, branchdata);
            }
         }
         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

   return SCIP_OKAY;
}

/** perform method of the given branchrule that is called after the master LP is solved */
SCIP_RETCODE GCGrelaxBranchMasterSolved(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   SCIP_Real             newlowerbound       /**< the new local lowerbound */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;
   
   assert(scip != NULL);
   assert(branchrule != NULL);
   
   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         /* call master problem solved method of the branching rule */
         if( relaxdata->branchrules[i]->branchmastersolved != NULL )
            SCIP_CALL( relaxdata->branchrules[i]->branchmastersolved(scip, branchdata, newlowerbound) );

         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

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

   /* get the number of the pricing block to which the variable belongs */
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
         TRUE, FALSE, gcgvardelorig, NULL, NULL, NULL, vardata) );

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

   /* create the vardata and initialize its values */
   SCIP_CALL( SCIPallocBlockMemory(scip, &vardata) );
   vardata->vartype = GCG_VARTYPE_ORIGINAL;
   vardata->blocknr = -1;
   vardata->data.origvardata.pricingvar = NULL;
   vardata->data.origvardata.coefs = NULL;
   vardata->data.origvardata.linkconss = NULL;
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

   /* loop over the variables in the original problem */
   for( i = 0; i < nvars; i++ )
   {
      assert(vars[i] != NULL);
      SCIP_CALL( GCGrelaxCreateOrigVardata(scip, vars[i]) );
   }

   return SCIP_OKAY;
}

/** transforms a constraint of the original problem into the master variable space 
 *  and stores information about the constraints in the vardatas */
SCIP_RETCODE GCGrelaxTransOrigToMasterCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint that should be transformed */
   SCIP_CONS**           transcons           /**< pointer to store the transformed constraint */
   )
{

   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_CONS* newcons;
   SCIP_CONS* mastercons;
   char name[SCIP_MAXSTRLEN];

   SCIP_VAR** mastervars;
   int nmastervars;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;
   int v;
   int i;
   int j;
   SCIP_VARDATA* vardata;
   SCIP_Real coef;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(cons != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   newcons = NULL;

   /* copy the constraint (dirty trick, we only need lhs and rhs, because variables are added later) */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "linear_%s", SCIPconsGetName(cons));
   SCIP_CALL( SCIPgetConsCopy(scip, scip, cons, &newcons, SCIPconsGetHdlr(cons),
         relaxdata->hashorig2origvar, NULL, name, 
         FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, &success) );
   
   assert(success && newcons != NULL);

   /* create and add corresponding linear constraint in the master problem */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "m_%s", SCIPconsGetName(cons));
   SCIP_CALL( SCIPcreateConsLinear(relaxdata->masterprob, &mastercons, name, 0, NULL, NULL, 
         SCIPgetLhsLinear(scip, newcons), SCIPgetRhsLinear(scip, newcons), 
         TRUE, TRUE, TRUE, TRUE, TRUE, SCIPconsIsLocal(cons), TRUE, FALSE, FALSE, 
         SCIPconsIsStickingAtNode(cons)) );
   
   /* now compute coefficients of the master variables in the master constraint */
   mastervars = SCIPgetVars(relaxdata->masterprob);
   nmastervars = SCIPgetNVars(relaxdata->masterprob);

   consvars = SCIPgetVarsLinear(scip, cons);
   nconsvars = SCIPgetNVarsLinear(scip, cons);
   consvals = SCIPgetValsLinear(scip, cons);


   /* add coefs of the original variables in the constraint to their variable data */
   for( v = 0; v < nconsvars; v++ )
   {
         assert(!SCIPisZero(scip, consvals[v]));
         vardata = SCIPvarGetData(consvars[v]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);

         SCIP_CALL( SCIPreallocMemoryArray(scip, &(vardata->data.origvardata.coefs), vardata->data.origvardata.ncoefs+1) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(vardata->data.origvardata.linkconss), vardata->data.origvardata.ncoefs+1) );

         assert(vardata->data.origvardata.coefs != NULL);
         assert(vardata->data.origvardata.linkconss != NULL);
         vardata->data.origvardata.coefs[vardata->data.origvardata.ncoefs] = consvals[v];
         vardata->data.origvardata.linkconss[vardata->data.origvardata.ncoefs] = mastercons;
         vardata->data.origvardata.ncoefs++;
   }

   /* add master variables to the corresponding master constraint */
   for( v = 0; v < nmastervars; v++ )
   {
      coef = 0;

      vardata = SCIPvarGetData(mastervars[v]);
      assert(vardata != NULL);
      assert(vardata->data.mastervardata.norigvars >= 0);
      assert((vardata->data.mastervardata.norigvars == 0) == (vardata->data.mastervardata.origvars == NULL));
      for( i = 0; i < vardata->data.mastervardata.norigvars; i++ )
         for( j = 0; j < nconsvars; j++ )
            if( consvars[j] == vardata->data.mastervardata.origvars[i] )
               coef += consvals[j] * vardata->data.mastervardata.origvals[i];

      if( !SCIPisFeasZero(scip, coef) )
      {
         SCIP_CALL( SCIPaddCoefLinear(relaxdata->masterprob, mastercons, mastervars[v], coef) );
      }
   }

   /* store the constraints in the arrays origmasterconss and masterconss in the problem data */
   SCIP_CALL( ensureSizeMasterConss(scip, relaxdata, relaxdata->nmasterconss+1) );
   SCIP_CALL( SCIPcaptureCons(scip, cons) );
   relaxdata->origmasterconss[relaxdata->nmasterconss] = cons;
   relaxdata->linearmasterconss[relaxdata->nmasterconss] = newcons;
   relaxdata->masterconss[relaxdata->nmasterconss] = mastercons;

   SCIP_CALL( GCGpricerAddMasterconsToHashmap(relaxdata->masterprob, relaxdata->masterconss[relaxdata->nmasterconss],
         relaxdata->nmasterconss) );

   relaxdata->nmasterconss++;

   *transcons = mastercons;

   return SCIP_OKAY;
}

/* prints the given variable: name, type (original, master or pricing) block number,
 * and the list of all variables related to the given variable */
void GCGrelaxPrintVar(
   SCIP_VAR*             var                 /**< variable that shpuld be printed */
   )
{
   SCIP_VARDATA* vardata;
   int i;

   vardata = SCIPvarGetData(var);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL || vardata->vartype == GCG_VARTYPE_MASTER
      || vardata->vartype == GCG_VARTYPE_PRICING );

   if( vardata->vartype == GCG_VARTYPE_ORIGINAL )
   {
      printf("Variable %s (original): block %d\n", SCIPvarGetName(var), vardata->blocknr);
      printf("mastervars:");
      for( i = 0; i < vardata->data.origvardata.nmastervars-1; i++ )
      {
         printf("%s (%g), ", SCIPvarGetName(vardata->data.origvardata.mastervars[i]),
            vardata->data.origvardata.mastervals[i]);
      }
         printf("%s (%g)\n", SCIPvarGetName(vardata->data.origvardata.mastervars[vardata->data.origvardata.nmastervars-1]),
            vardata->data.origvardata.mastervals[vardata->data.origvardata.nmastervars-1]);
   }
   else if( vardata->vartype == GCG_VARTYPE_PRICING )
   {
      printf("Variable %s (pricing): block %d\n", SCIPvarGetName(var), vardata->blocknr);
      printf("origvars:");
      for( i = 0; i < vardata->data.pricingvardata.norigvars-1; i++ )
      {
         printf("%s, ", SCIPvarGetName(vardata->data.pricingvardata.origvars[i]));
      }
      printf("%s\n", SCIPvarGetName(vardata->data.pricingvardata.origvars[vardata->data.pricingvardata.norigvars-1]));
   }
   else if( vardata->vartype == GCG_VARTYPE_MASTER )
   {
      printf("Variable %s (master): block %d\n", SCIPvarGetName(var), vardata->blocknr);
      printf("origvars:");
      for( i = 0; i < vardata->data.mastervardata.norigvars-1; i++ )
      {
         printf("%s (%g), ", SCIPvarGetName(vardata->data.mastervardata.origvars[i]),
            vardata->data.mastervardata.origvals[i]);
      }
         printf("%s (%g)\n", SCIPvarGetName(vardata->data.mastervardata.origvars[vardata->data.mastervardata.norigvars-1]),
            vardata->data.mastervardata.origvals[vardata->data.mastervardata.norigvars-1]);
   }
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
   assert(vardata->blocknr == -1 || vardata->blocknr == blocknr);

   vardata->blocknr = blocknr;

   return SCIP_OKAY;
}

/* marks the constraint to be transferred to the master problem */
SCIP_RETCODE GCGrelaxMarkConsMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint that is forced to be in the master */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
#ifndef NDEBUG
   int i;
#endif   
   assert(scip != NULL);
   assert(cons != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* allocate array, if not yet done */
   if( relaxdata->markedmasterconss == NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->markedmasterconss), SCIPgetNConss(scip)) );
      relaxdata->nmarkedmasterconss = 0;
   }
   assert(relaxdata->nmarkedmasterconss + 1 < SCIPgetNConss(scip));

#ifndef NDEBUG
   /* check that constraints are not marked more than one time */
   for( i = 0; i < relaxdata->nmarkedmasterconss; i++ )
      assert(relaxdata->markedmasterconss[i] != cons);
#endif

   /* save constraint */
   relaxdata->markedmasterconss[relaxdata->nmarkedmasterconss] = cons;
   relaxdata->nmarkedmasterconss++;

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

/** returns TRUE iff the pricing problem of the given number is relevant, that means is not identical to
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

   assert(relaxdata->nblocksidentical[pricingprobnr] >= 0);
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
   assert(npricingprobs >= 0);

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

/* returns the linking constraints in the original problem that correspond to the constraints in the master problem */
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
   assert(blocknr >= 0);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(blocknr < relaxdata->npricingprobs);

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

/** start probing mode on master problem */
SCIP_RETCODE GCGrelaxStartProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP* masterscip;

   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(!relaxdata->masterinprobing);

   masterscip = relaxdata->masterprob;
   assert(masterscip != NULL);

   /* create probing node in master problem, propagate and solve it with pricing */
   SCIP_CALL( SCIPstartProbing(masterscip) );

   relaxdata->masterinprobing;

   return SCIP_OKAY;
}


/** for a probing node in the original problem, create a corresponding probing node in the master problem,
 *  propagate domains and solve the LP with pricing. */
SCIP_RETCODE GCGrelaxPerformProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint*         nlpiterations,      /**< pointert to store the number of used LP iterations */
   SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the probing direction is infeasible */
   SCIP_Bool*            feasible            /**< pointer to store whether the probing solution is feasible */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP* masterscip;
   SCIP_NODE* mprobingnode;
   SCIP_CONS* mprobingcons;
   SCIP_LPSOLSTAT lpsolstat;

   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   masterscip = relaxdata->masterprob;
   assert(masterscip != NULL);

   /* create probing node in master problem, propagate and solve it with pricing */
   SCIPnewProbingNode(masterscip);

   mprobingnode = SCIPgetCurrentNode(masterscip);
   assert(GCGconsMasterbranchGetActiveCons(masterscip) != NULL);
   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &mprobingcons, mprobingnode, 
         GCGconsMasterbranchGetActiveCons(masterscip)) );
   SCIP_CALL( SCIPaddConsNode(masterscip, mprobingnode, mprobingcons, NULL) );
   SCIP_CALL( SCIPreleaseCons(scip, &mprobingcons) );

   SCIP_CALL( SCIPpropagateProbing(masterscip, -1, cutoff, NULL) );
   assert(!(*cutoff));
      
   SCIP_CALL( SCIPsolveProbingLPWithPricing( masterscip, FALSE/* pretendroot */, FALSE /*displayinfo*/,
         -1 /*maxpricerounds*/, lperror ) );
   lpsolstat = SCIPgetLPSolstat(masterscip);

   *nlpiterations += SCIPgetNLPIterations(masterscip);

   if( !(*lperror) )
   {
      /* get LP solution status, objective value */
      *cutoff = *cutoff || (lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT || lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && SCIPisLPRelax(masterscip) )
      {
         SCIPdebugMessage("lpobjval = %g\n", SCIPgetLPObjval(masterscip));
         *lpobjvalue = SCIPgetLPObjval(masterscip);
         *lpsolved = TRUE;
         SCIP_CALL( GCGrelaxUpdateCurrentSol(scip, feasible) );
      }
   }
   else 
   {
      SCIPinfoMessage(scip, NULL, "something went wrong, an lp error occured\n");         
   }

   return SCIP_OKAY;
}


/** end probing mode in master problem */
SCIP_RETCODE GCGrelaxEndProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP* masterscip;

   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->masterinprobing);

   masterscip = relaxdata->masterprob;
   assert(masterscip != NULL);

   SCIP_CALL( SCIPendProbing(masterscip) );

   /* if a new primal solution was found in the master problem, transfer it to the original problem */
   if( SCIPgetBestSol(relaxdata->masterprob) != NULL && relaxdata->lastmastersol != SCIPgetBestSol(relaxdata->masterprob) )
   {
      SCIP_SOL* newsol;
      SCIP_Bool stored;

      relaxdata->lastmastersol = SCIPgetBestSol(relaxdata->masterprob);

      SCIP_CALL( GCGrelaxTransformMastersolToOrigsol(scip, relaxdata->lastmastersol, &newsol) );

      SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, TRUE, TRUE, TRUE, &stored) );
      if( !stored )
      {
         SCIP_CALL( SCIPcheckSolOrig(scip, newsol, &stored, TRUE, TRUE) );
      }
      assert(stored);
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );

      SCIPdebugMessage("probing finished in master problem\n");
   }

   /* restore old relaxation solution and branching candidates */

   /* TODO: solve master problem again */

   return SCIP_OKAY;
}


/** transforms the current solution of the master problem into the original problem's space 
 *  and saves this solution as currentsol in the relaxator's data */
SCIP_RETCODE GCGrelaxUpdateCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            feasible            /**< pointer to store whether the master problem's solution is 
                                              *   primal feasible*/
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_VAR** origvars;
   int norigvars;
   SCIP_SOL* mastersol;
   SCIP_Bool stored;
   int i;

   assert(scip != NULL);
   assert(feasible != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   origvars = SCIPgetVars(scip);
   norigvars = SCIPgetNVars(scip);
   assert(origvars != NULL);

   *feasible = FALSE;

   /* free previous solution and clear branching candidates */
   if( relaxdata->currentorigsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &(relaxdata->currentorigsol)) );
   }
   SCIPclearExternBranchCands(scip);
   
   /** @todo: remove the TRUE of the if condition and use correct abort criteria */
   /* nothing has to be done, if no LP was solved after the last update */
   /*if( TRUE || relaxdata->lastmasterlpiters != SCIPgetNLPIterations(relaxdata->masterprob) )*/
   if( SCIPgetStage(relaxdata->masterprob) == SCIP_STAGE_SOLVED || SCIPgetLPSolstat(relaxdata->masterprob) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      //printf("nlpiterations = %lld, lastlpiterations = %lld\n", SCIPgetNLPIterations(relaxdata->masterprob), relaxdata->lastmasterlpiters);
      relaxdata->lastmasterlpiters = SCIPgetNLPIterations(relaxdata->masterprob);

      /* create new solution */
      if( SCIPgetStage(relaxdata->masterprob) == SCIP_STAGE_SOLVING )
         mastersol = NULL;
      else if( SCIPgetStage(relaxdata->masterprob) == SCIP_STAGE_SOLVED )
      {
         mastersol = SCIPgetBestSol(relaxdata->masterprob);
         if( mastersol == NULL )
            return SCIP_OKAY;
      }
      else 
      {
         printf("stage in master not solving and not solved!\n");
         return SCIP_OKAY;
      }

      if( !SCIPisInfinity(scip, SCIPgetSolOrigObj(relaxdata->masterprob, mastersol)) )
      {
         /* transform the master solution to the original variable space */
         SCIP_CALL( GCGrelaxTransformMastersolToOrigsol(scip, mastersol, &(relaxdata->currentorigsol)) );

         /* store the solution as relaxation solution */
         SCIP_CALL( SCIPsetRelaxSolValsSol(scip, relaxdata->currentorigsol) );
         assert(SCIPisEQ(scip, SCIPgetRelaxSolObj(scip), SCIPgetSolTransObj(scip, relaxdata->currentorigsol)));
         
         SCIP_CALL( SCIPtrySol(scip, relaxdata->currentorigsol, FALSE, TRUE, TRUE, TRUE, &stored) );
         if( !stored )
         {
            SCIP_CALL( SCIPcheckSol(scip, relaxdata->currentorigsol, FALSE, TRUE, TRUE, TRUE, &stored) );
         }

         SCIPdebugMessage("updated current original LP solution, %s feasible in the original problem!\n",
            (stored ? "" : "not"));

         if( stored )
            *feasible = TRUE;

         /* store branching candidates */
         for( i = 0; i < norigvars; i++ )
            if( SCIPvarGetType(origvars[i]) <= SCIP_VARTYPE_INTEGER && !SCIPisFeasIntegral(scip, SCIPgetRelaxSolVal(scip, origvars[i])) )
            {
               assert(!SCIPisEQ(scip, SCIPvarGetLbLocal(origvars[i]), SCIPvarGetUbLocal(origvars[i])));
               SCIP_CALL( SCIPaddExternBranchCand(scip, origvars[i], SCIPgetRelaxSolVal(scip, 
                        origvars[i]) - SCIPfloor(scip, SCIPgetRelaxSolVal(scip, origvars[i])), 
                     SCIPgetRelaxSolVal(scip, origvars[i])) );
            }
         SCIPdebugMessage("updated relaxation branching candidates\n");
      }
   }
   /* if a new primal solution was found in the master problem, transfer it to the original problem */
   if( SCIPgetBestSol(relaxdata->masterprob) != NULL && relaxdata->lastmastersol != SCIPgetBestSol(relaxdata->masterprob) )
   {
      SCIP_SOL* newsol;

      relaxdata->lastmastersol = SCIPgetBestSol(relaxdata->masterprob);

      SCIP_CALL( GCGrelaxTransformMastersolToOrigsol(scip, relaxdata->lastmastersol, &newsol) );

      SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, &stored) );
      if( !stored )
      {
         SCIP_CALL( SCIPcheckSolOrig(scip, newsol, &stored, TRUE, TRUE) );
      }
      assert(stored);
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );

      SCIPdebugMessage("updated current best primal feasible solution!\n");
   }

   return SCIP_OKAY;
}        


/** transforms given values of the given original variables into values of the given master variables */
void GCGrelaxTransformOrigvalsToMastervals(
   SCIP*                 scip,               /** SCIP data structure */
   SCIP_VAR**            origvars,           /** array with (subset of the) original variables */
   SCIP_Real*            origvals,           /** array with values (coefs) for the given original variables */
   int                   norigvars,          /** number of given original variables */
   SCIP_VAR**            mastervars,         /** array of (all present) master variables */
   SCIP_Real*            mastervals,         /** array to store the values of the master variables */
   int                   nmastervars         /** number of master variables */
   )
{
   SCIP_VARDATA* vardata;
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

   /* set all values to 0 initially */
   for( i = 0; i < nmastervars; i++ )
      mastervals[i] = 0.0;

   /* iterate over all original variables */
   for( i = 0; i < norigvars; i++ )
   {
      vardata = SCIPvarGetData(origvars[i]);

      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->data.origvardata.nmastervars >= 0);
      assert(vardata->data.origvardata.mastervars != NULL);
      assert(vardata->data.origvardata.mastervals != NULL);
      assert(vardata->data.origvardata.nmastervars == 1 || vardata->blocknr != -1);

      /* variable belongs to no block, so it was transferred directly to the master problem, 
       * hence, we transfer the solution value directly to the corresponding master variabe
       */
      if( vardata->blocknr == -1 )
      {
         for( k = 0; k < nmastervars; k++ )
         {
            assert(!SCIPvarIsTransformedOrigvar(mastervars[k]));
            if( mastervars[k] == vardata->data.origvardata.mastervars[0])
            {
               assert(!SCIPvarIsTransformedOrigvar(vardata->data.origvardata.mastervars[0]));
               mastervals[k] += (vardata->data.origvardata.mastervals[0] * origvals[i]);
               break;
            }
         }
         assert(k < nmastervars);
      }
      /* variable belongs to a block, so we have to look at all master variables and increase their values 
       * if they contain the original variable
       */
      else
      {
         vardata = SCIPvarGetData(vardata->data.origvardata.pricingvar);
         assert(vardata->vartype == GCG_VARTYPE_PRICING);
         assert(vardata->data.pricingvardata.norigvars >= 0);
         assert(vardata->data.pricingvardata.origvars != NULL);
         assert(vardata->data.pricingvardata.origvars[0] != NULL);

         vardata = SCIPvarGetData(vardata->data.pricingvardata.origvars[0]);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         assert(vardata->data.origvardata.nmastervars >= 0);
         assert(vardata->data.origvardata.mastervars != NULL);
         assert(vardata->data.origvardata.mastervals != NULL);

         for( j = 0; j < vardata->data.origvardata.nmastervars; j++ )
         {
            for( k = 0; k < nmastervars; k++ )
               if( mastervars[k] == vardata->data.origvardata.mastervars[j] )
               {
                  mastervals[k] += (vardata->data.origvardata.mastervals[j] * origvals[i]);
                  break;
               }
            assert(k < nmastervars);
         }
      }

   }
}       

/** transforms given solution of the master problem into solution of the original problem
 *  TODO: think about types of epsilons used in this method*/
SCIP_RETCODE GCGrelaxTransformMastersolToOrigsol(
   SCIP*                 scip,               /** SCIP data structure */
   SCIP_SOL*             mastersol,          /** solution of the master problem, or NULL for current LP solution */
   SCIP_SOL**            origsol             /** pointer to store the new created original problem's solution */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_VARDATA* vardata;
   SCIP_VARDATA* vardata2;
   int* blocknr;
   SCIP_Real* blockvalue;
   SCIP_Real increaseval;
   SCIP_VAR** mastervars;
   SCIP_Real* mastervals;
   int nmastervars;
   int i;
   int j;

   assert(scip != NULL);
   assert(origsol != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   assert( !SCIPisInfinity(scip, SCIPgetSolOrigObj(relaxdata->masterprob, mastersol)) );

   SCIP_CALL( SCIPcreateSol(scip, origsol, NULL) );

   SCIP_CALL( SCIPallocBufferArray(scip, &blockvalue, relaxdata->npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocknr, relaxdata->npricingprobs) );

   /* get variables of the master problem and their solution values */
   SCIP_CALL( SCIPgetVarsData(relaxdata->masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &mastervals, nmastervars) );
   SCIP_CALL( SCIPgetSolVals(relaxdata->masterprob, mastersol, nmastervars, mastervars, mastervals) );

   /* initialize the block values for the pricing problems */
   for( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      blockvalue[i] = 0.0;
      blocknr[i] = 0;
   }

   /* loop over all given master variables */
   for( i = 0; i < nmastervars; i++ )
   {
      vardata = SCIPvarGetData(mastervars[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_MASTER);
      assert(vardata->data.mastervardata.norigvars >= 0);
      assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
      assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);

      assert(!SCIPisFeasNegative(scip, mastervals[i]));

      /* TODO: handle infinite master solution values */
      assert(!SCIPisInfinity(scip, mastervals[i]));

      /* first of all, handle variables representing rays */
      if( vardata->data.mastervardata.isray )
      {
         assert(vardata->blocknr != -1);
         /* we also want to take into account variables representing rays, that have a small value (between normal and feas eps), 
          * so we do no feas comparison here */
         if( SCIPisPositive(scip, mastervals[i]) )
         {
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
            {
               assert(!SCIPisZero(scip, vardata->data.mastervardata.origvals[j]));
               
               /* increase the corresponding value */
               SCIP_CALL( SCIPincSolVal(scip, *origsol, vardata->data.mastervardata.origvars[j], vardata->data.mastervardata.origvals[j] * mastervals[i]) );
            }
         }
         mastervals[i] = 0.0;
         continue;
      }

      /* handle the variables with integral values */
      while( SCIPisFeasGE(scip, mastervals[i], 1) )
      {

         if( vardata->blocknr == -1 )
         {
#ifdef ORIGVARS
            assert(vardata->data.mastervardata.norigvars == 2);
            assert(vardata->data.mastervardata.origvals[0] == 1.0);
            assert(vardata->data.mastervardata.origvals[1] == 0.0);
#endif
            /* increase the corresponding value */
            SCIP_CALL( SCIPincSolVal(scip, *origsol, vardata->data.mastervardata.origvars[0], vardata->data.mastervardata.origvals[0] * mastervals[i]) );
            mastervals[i] = 0.0;
         }
         else
         {
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
            {
               assert(!SCIPisZero(scip, vardata->data.mastervardata.origvals[j]));

               /* get the right original variable */
               vardata2 = SCIPvarGetData(vardata->data.mastervardata.origvars[j]);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);
               assert(vardata2->data.origvardata.pricingvar != NULL);
               vardata2 = SCIPvarGetData(vardata2->data.origvardata.pricingvar);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_PRICING);

               /* just in case a variable has a value higher than the number of blocks, it represents */
               if( vardata2->data.pricingvardata.norigvars <= blocknr[vardata->blocknr] )
               {
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, vardata2->data.pricingvardata.origvars[vardata2->data.pricingvardata.norigvars-1], mastervals[i] * vardata->data.mastervardata.origvals[j]) );
                  mastervals[i] = 1.0;
               }
               /* this should be default */
               else
               {              
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, vardata2->data.pricingvardata.origvars[blocknr[vardata->blocknr]], vardata->data.mastervardata.origvals[j]) );
               }
            }
            mastervals[i] = mastervals[i] - 1.0;
            blocknr[vardata->blocknr]++;
         }
      }
   }

   /* loop over all given master variables */
   for( i = 0; i < nmastervars; i++ )
   {
      if( SCIPisFeasZero(scip, mastervals[i]) )
      {
         continue;
      }
      assert(SCIPisFeasGE(scip, mastervals[i], 0.0) && SCIPisFeasLT(scip, mastervals[i], 1.0));

      while( SCIPisFeasPositive(scip, mastervals[i]) )
      {
         vardata = SCIPvarGetData(mastervars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_MASTER);
         assert(vardata->data.mastervardata.norigvars >= 0);
         assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
         assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);
         assert(!vardata->data.mastervardata.isray);
         
         if( vardata->blocknr == -1 )
         {
#ifdef ORIGVARS
            assert(vardata->data.mastervardata.norigvars == 2);
            assert(vardata->data.mastervardata.origvals[0] == 1.0);
            assert(vardata->data.mastervardata.origvals[1] == 0.0);
#endif            
            /* increase the corresponding value */
            SCIP_CALL( SCIPincSolVal(scip, *origsol, vardata->data.mastervardata.origvars[0], vardata->data.mastervardata.origvals[0] * mastervals[i]) );
            mastervals[i] = 0.0;
         }
         else
         {
            increaseval = MIN(mastervals[i], 1.0 - blockvalue[vardata->blocknr]);
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
            {
               if( SCIPisZero(scip, vardata->data.mastervardata.origvals[j]) )
                  continue;

               /* get the right original variable */
               vardata2 = SCIPvarGetData(vardata->data.mastervardata.origvars[j]);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);
               assert(vardata2->data.origvardata.pricingvar != NULL);
               vardata2 = SCIPvarGetData(vardata2->data.origvardata.pricingvar);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_PRICING);
               
               if( vardata2->data.pricingvardata.norigvars <= blocknr[vardata->blocknr] )
               {
                  increaseval = mastervals[i];
                  
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, vardata2->data.pricingvardata.origvars[vardata2->data.pricingvardata.norigvars-1], vardata->data.mastervardata.origvals[j] * increaseval) );
               }
               else
               {
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, vardata2->data.pricingvardata.origvars[blocknr[vardata->blocknr]], vardata->data.mastervardata.origvals[j] * increaseval) );
               }
            }

            mastervals[i] = mastervals[i] - increaseval;
            if( SCIPisFeasZero(scip, mastervals[i]) )
            {
               mastervals[i] = 0.0;
            }
            blockvalue[vardata->blocknr] += increaseval;

            /* if the value assigned to the block is equal to 1, this block is full and we take the next block */
            if( SCIPisFeasGE(scip, blockvalue[vardata->blocknr], 1.0) )
            {
               blockvalue[vardata->blocknr] = 0.0;
               blocknr[vardata->blocknr]++;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &mastervals);
   SCIPfreeBufferArray(scip, &blocknr);
   SCIPfreeBufferArray(scip, &blockvalue);

   return SCIP_OKAY;

}
