/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* #define SCIP_DEBUG */
//#define CHECKCONSISTENCY
/**@file    relax_gcg.c
 * @ingroup RELAXATORS
 * @brief   gcg relaxator
 * @author  Gerald Gamrath
 * @author  Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scipdefplugins.h"

#include "relax_gcg.h"

#include "struct_branchgcg.h"

#include "cons_origbranch.h"
#include "cons_masterbranch.h"
#include "cons_connected.h"
#include "pricer_gcg.h"
#include "masterplugins.h"
#include "nodesel_master.h"
#include "pub_gcgvar.h"

#define RELAX_NAME             "gcg"
#define RELAX_DESC             "relaxator for gcg project representing the master lp"
#define RELAX_PRIORITY         -1
#define RELAX_FREQ             1

#define DEFAULT_DISCRETIZATION FALSE
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
   int              nlinkingvars;        /* number of linking variables */
   int              nvarlinkconss;       /* number of constraints that ensure that copies of linking variables have the same value */

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

   /* solution data */
   SCIP_SOL*        origprimalsol;       /* best original primal solution */

   /* structure information */
   SCIP_Bool        hasblockdetection;   /* indicates whether the block detection code is present */
};

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

/** checks whether two arrays of SCIP_Real's are identical
 * @todo: What about using SCIPisEq()
 */

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

   SCIP_CONS** conss1;
   SCIP_CONS** conss2;
   int nconss;

   SCIP_VAR** origvars1;
   SCIP_VAR** origvars2;

   SCIP_Real* coefs1;
   int ncoefs1;
   SCIP_Real* coefs2;
   int ncoefs2;
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
      if( !SCIPisEQ(relaxdata->masterprob, SCIPvarGetObj(vars1[i]), SCIPvarGetObj(vars2[i])) )
      {
         SCIPdebugMessage("--> obj differs for var %s and var %s!\n", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
         return SCIP_OKAY;
      }
      if( !SCIPisEQ(relaxdata->masterprob, SCIPvarGetLbOriginal(vars1[i]), SCIPvarGetLbOriginal(vars2[i])) )
      {
         SCIPdebugMessage("--> lb differs for var %s and var %s!\n", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
         return SCIP_OKAY;
      }
      if( !SCIPisEQ(relaxdata->masterprob, SCIPvarGetUbOriginal(vars1[i]), SCIPvarGetUbOriginal(vars2[i])) )
      {
         SCIPdebugMessage("--> ub differs for var %s and var %s!\n", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
         return SCIP_OKAY;
      }
      if( SCIPvarGetType(vars1[i]) != SCIPvarGetType(vars2[i]) )
      {
         SCIPdebugMessage("--> type differs for var %s and var %s!\n", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
         return SCIP_OKAY;
      }

      assert(GCGvarIsPricing(vars1[i]));
      assert(GCGvarIsPricing(vars2[i]));

      origvars1 = GCGpricingVarGetOrigvars(vars1[i]);
      origvars2 = GCGpricingVarGetOrigvars(vars2[i]);

      if( !SCIPisEQ(relaxdata->masterprob, SCIPvarGetObj(origvars1[0]),SCIPvarGetObj(origvars2[0])) )
      {
         SCIPdebugMessage("--> orig obj differs for var %s and var %s!\n", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
         return SCIP_OKAY;
      }

      coefs1 = GCGoriginalVarGetCoefs(origvars1[0]);
      ncoefs1 = GCGoriginalVarGetNCoefs(origvars2[0]);
      coefs2 = GCGoriginalVarGetCoefs(origvars1[0]);
      ncoefs2 = GCGoriginalVarGetNCoefs(origvars2[0]);

      assert(GCGvarIsOriginal(origvars1[0]));
      assert(GCGvarIsOriginal(origvars2[0]));

      if( !realArraysAreEqual(coefs1, ncoefs1, coefs2, ncoefs2) )
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
      if( !SCIPisEQ(relaxdata->masterprob, SCIPgetLhsLinear(scip1, conss1[i]), SCIPgetLhsLinear(scip2, conss2[i])) )
      {
         SCIPdebugMessage("--> lhs differs for cons %s and cons %s!\n", SCIPconsGetName(conss1[i]), SCIPconsGetName(conss2[i]));
         return SCIP_OKAY;
      }
      if( !SCIPisEQ(relaxdata->masterprob, SCIPgetRhsLinear(scip1, conss1[i]), SCIPgetRhsLinear(scip2, conss2[i])) )
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
            /* save variables in pricing problem variable */
            vars = SCIPgetVars(relaxdata->pricingprobs[i]);
            nvars = SCIPgetNVars(relaxdata->pricingprobs[i]);
            for( k = 0; k < nvars; k++ )
            {
               int blocknr;
               assert(GCGvarIsPricing(vars[k]));
               origvar = GCGpricingVarGetOrigvars(vars[k])[0];

               pricingvar = (SCIP_VAR*) SCIPhashmapGetImage(varmap, (void*) vars[k]);
               blocknr = GCGvarGetBlock(pricingvar);

               assert(GCGvarIsPricing(pricingvar));
               assert(GCGvarIsOriginal(origvar));
               assert(GCGoriginalVarGetPricingVar(origvar) != NULL);
               GCGoriginalVarSetPricingVar(origvar, pricingvar);
               SCIP_CALL( GCGpricingVarAddOrigVar(relaxdata->pricingprobs[blocknr], pricingvar, origvar) );
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

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Matrix has %d blocks, %d %s relevant!\n", relaxdata->npricingprobs, nrelevant,
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

   nvars = 0; /* fix potential problems */
   vars = NULL;

   /* TODO: maybe change that to SCIPgetNVarsXXX/SCIPgetNConssXXX */
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
      var = SCIPvarGetProbvar(var);
      if( !SCIPhashmapExists(varmap, (void*) var) )
         return FALSE;

      /* check whether bounding variable is contained in block */
      var = SCIPgetVbdvarVarbound(scip, cons);
      var = SCIPvarGetProbvar(var);
      if( !SCIPhashmapExists(varmap, (void*) var) )
         return FALSE;

      /* both variables are in the block, return TRUE */
      return TRUE;
   }
   else
   {
      SCIPerrorMessage("constraint %s of unknown type <%s>!\n", SCIPconsGetName(cons), SCIPconshdlrGetName(SCIPconsGetHdlr(cons)));
   }

   assert(vars != NULL || nvars == 0);

   for( i = 0; i < nvars; i++ )
   {
      var = vars[i];
      var = SCIPvarGetProbvar(vars[i]);
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
   //SCIP_CALL( SCIPsetBoolParam(relaxdata->masterprob, "pricing/delvars", TRUE) );
   //SCIP_CALL( SCIPsetBoolParam(relaxdata->masterprob, "pricing/delvarsroot", TRUE) );
   //SCIP_CALL( SCIPsetBoolParam(relaxdata->masterprob, "lp/cleanupcols", TRUE) );
   //SCIP_CALL( SCIPsetBoolParam(relaxdata->masterprob, "lp/cleanupcolsroot", TRUE) );
   SCIP_CALL( SCIPsetRealParam(relaxdata->masterprob, "pricing/abortfac", 1.0) );
   SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "timing/clocktype", 2) );

   /* ----- initialize the pricing problems ----- */
   npricingprobs = relaxdata->npricingprobs;
   relaxdata->npricingprobs = npricingprobs > 0 ? npricingprobs : 0;
   if( npricingprobs > 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->pricingprobs), npricingprobs) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->blockrepresentative), npricingprobs) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->nblocksidentical), npricingprobs) );

      /* array for saving convexity constraints belonging to one of the pricing problems */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->convconss), npricingprobs) );

      /* create hashmaps for mapping from original to pricing variables */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->hashorig2pricingvar), npricingprobs) );
   }

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

      /* disable expensive presolving */
//      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "presolving/probing/maxrounds", 0) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/linear/presolpairwise", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/setppc/presolpairwise", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/logicor/presolpairwise", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/linear/presolusehashing", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/setppc/presolusehashing", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "constraints/logicor/presolusehashing", FALSE) );

      /* disable output to console */
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "display/verblevel", SCIP_VERBLEVEL_NONE) );
#if SCIP_VERSION > 210
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "misc/printreason", FALSE) );
#endif
      /* disable solution store */
      //SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "limits/maxorigsol", 0) );

      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(relaxdata->pricingprobs[i], "misc/catchctrlc", FALSE) );
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "timing/clocktype", 2) );

      /* disable solution store */
      //      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "limits/maxorigsol", 0) );

      /* create the pricing submip */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricing_block_%d", i);
      SCIP_CALL( SCIPcreateProb(relaxdata->pricingprobs[i], name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   }

   for( i = 0; i < npricingprobs; i++ )
   {
      SCIP_CALL( SCIPhashmapCreate(&(relaxdata->hashorig2pricingvar[i]),
            SCIPblkmem(scip), SCIPgetNVars(scip)) );
   }
   SCIP_CALL( SCIPhashmapCreate(&(relaxdata->hashorig2origvar),
         SCIPblkmem(scip), 10*SCIPgetNVars(scip)+1) );

   /* create pricing variables and map them to the original variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for( v = 0; v < nvars; v++ )
   {
      int blocknr;
      SCIP_VAR* var;
      var = SCIPvarGetProbvar(vars[v]);
      blocknr = GCGvarGetBlock(var);

      /* variable belongs to exactly one block --> create corresponding pricing variable*/
      if( blocknr >= 0 )
      {
         assert(GCGoriginalVarGetPricingVar(var) == NULL);

         SCIP_CALL( GCGrelaxCreatePricingVar(scip, var) );
         assert(GCGoriginalVarGetPricingVar(var) != NULL);

         SCIP_CALL( SCIPhashmapInsert(relaxdata->hashorig2pricingvar[blocknr], (void*)(var),
               (void*)(GCGoriginalVarGetPricingVar(var)) ));
         SCIP_CALL( SCIPhashmapInsert(relaxdata->hashorig2origvar, (void*)(var), (void*)(var)) );
      }
      /* variable is a linking variable --> create corresponding pricing variable in all linked blocks
       * and create corresponding linking constraints */
      else if( GCGvarIsLinking(var) )
      {
         SCIP_VAR** pricingvars;

         assert(GCGoriginalVarGetPricingVar(var) == NULL);

#ifndef NDEBUG
         /* checks that GCGrelaxSetOriginalVarBlockNr() worked correctly */
         {
            int count;
            int nblocks;
            SCIP_CONS** linkconss;

            pricingvars = GCGlinkingVarGetPricingVars(var);
            linkconss = GCGlinkingVarGetLinkingConss(var);
            nblocks = GCGlinkingVarGetNBlocks(var);

            count = 0;
            for( i = 0; i < npricingprobs; i++ )
            {
               if( pricingvars[i] != NULL)
               {
                  count++;
                  //assert(pricingvars[i] == vars[v]);
               }
               assert(linkconss[i] == NULL);
            }
            assert(nblocks == count);
         }
#endif
         SCIP_CALL( GCGrelaxCreateLinkingPricingVars(scip, var) );
#ifndef NDEBUG
         /* checks that GCGrelaxCreateLinkingPricingVars() worked correctly */
         {
            int count;
            int nblocks;
            SCIP_CONS** linkconss;

            pricingvars = GCGlinkingVarGetPricingVars(var);
            linkconss = GCGlinkingVarGetLinkingConss(var);
            nblocks = GCGlinkingVarGetNBlocks(var);

            count = 0;
            for( i = 0; i < npricingprobs; i++ )
            {
               if( pricingvars[i] != NULL)
               {
                  count++;
                  assert(GCGvarIsPricing(pricingvars[i]));
                  assert(linkconss[i] != NULL);
               }
               else
                  assert(linkconss[i] == NULL);
            }
            assert(nblocks == count);
         }
#endif
         assert(GCGoriginalVarGetPricingVar(var) == NULL);

         pricingvars = GCGlinkingVarGetPricingVars(var);

         for( i = 0; i < npricingprobs; i++ )
         {
            if( pricingvars[i] != NULL)
            {
               SCIP_CALL( SCIPhashmapInsert(relaxdata->hashorig2pricingvar[i], (void*)(var),
                     (void*)(pricingvars[i])) );
            }
         }
         SCIP_CALL( SCIPhashmapInsert(relaxdata->hashorig2origvar, (void*)(var), (void*)(var)) );
      }
      else
      {
         assert(GCGvarGetBlock(var) == -1);
         assert(GCGoriginalVarGetPricingVar(var) == NULL);
         SCIP_CALL( SCIPhashmapInsert(relaxdata->hashorig2origvar, (void*)(var), (void*)(var)) );
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
               SCIP_Bool copied;
               copied = FALSE;
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

                     SCIPdebugMessage("copying %s to pricing problem %d\n",  SCIPconsGetName(bufconss[c]), b);
                     copied = TRUE;
                     /* constraint was successfully copied */
                     assert(success);

                     SCIP_CALL( SCIPaddCons(relaxdata->pricingprobs[b], newcons) );

                     SCIP_CALL( SCIPreleaseCons(relaxdata->pricingprobs[b], &newcons) );
                  }
               }
               assert(copied);
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

   /* for original variables, save the coefficients in the master problem */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for( v = 0; v < nvars; v++ )
   {
      SCIP_VAR* var;
      var = SCIPvarGetProbvar(vars[v]);
      assert(GCGvarIsOriginal(var));
      assert(GCGoriginalVarGetCoefs(var) == NULL);
      GCGoriginalVarSetNCoefs(var, 0);
   }

   /* save coefs */
   for( i = 0; i < relaxdata->nmasterconss; i++ )
   {
      vars = SCIPgetVarsLinear(scip, relaxdata->linearmasterconss[i]);
      nvars = SCIPgetNVarsLinear(scip, relaxdata->linearmasterconss[i]);
      vals = SCIPgetValsLinear(scip, relaxdata->linearmasterconss[i]);
      for( v = 0; v < nvars; v++ )
      {
         SCIP_CALL( GCGoriginalVarAddCoef(scip, vars[v], vals[v], relaxdata->masterconss[i]) );
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
            relaxdata->nblocksidentical[i]*1.0, relaxdata->nblocksidentical[i]*1.0,
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

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "pricing problem %d: %d conss, %d vars (%d bins, %d ints, %d impls and %d cont)\n", i,
         SCIPgetNConss(relaxdata->pricingprobs[i]), SCIPgetNVars(relaxdata->pricingprobs[i]), nbin, nint, nimpl, ncont);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricingprob_%d.lp", i);

      SCIP_CALL( SCIPwriteOrigProblem(relaxdata->pricingprobs[i], name, NULL, FALSE) );
   }

   return SCIP_OKAY;
}

/** combines the solutions from all (disjoint) problems to one solution */
static
SCIP_RETCODE combineSolutions(
   SCIP*      scip,     /**< SCIP data structure */
   SCIP_SOL** newsol,   /**< pointer to store new solution */
   SCIP**     probs,    /**< array of (solved) subproblems */
   int        nprobs    /**< number of subproblems */
   )
{
#ifdef SCIP_DEBUG
   int i;
#endif

   int v;
   SCIP_SOL* sol;
   int nvars;

   SCIP_VAR** vars;
   assert(scip != NULL);
   assert(newsol != NULL);
   assert(probs != NULL);
   assert(nprobs > 0);

   SCIP_CALL( SCIPcreateSol(scip, newsol, NULL) );
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

#ifdef SCIP_DEBUG
   for (i = 0; i < nprobs; ++i)
   {
      if( probs[i] == NULL )
         continue;

      SCIPprintOrigProblem(probs[i], NULL, "lp", FALSE);
      SCIPprintSol(probs[i], SCIPgetBestSol(probs[i]), NULL, FALSE );
   }
#endif

   for( v = 0; v < nvars; ++v)
   {
      SCIP_VAR* pricingvar;
      int block;

      pricingvar = GCGoriginalVarGetPricingVar(vars[v]);
      block = GCGvarGetBlock(vars[v]);
      assert(block >= 0);
      assert(block < nprobs);
      assert(probs[block] != NULL);
      sol = SCIPgetBestSol(probs[block]);
      SCIP_CALL( SCIPincSolVal(scip, *newsol, vars[v], SCIPgetSolVal(probs[block], sol, pricingvar)) );
   }
   return SCIP_OKAY;
}

/** sets the pricing objective function to what is necessary */
static
SCIP_RETCODE setPricingObjsOriginal(
   SCIP*    scip,    /**< SCIP data structure */
   SCIP**   probs,   /**< array of subproblems */
   int      nprobs   /**< number of subproblems */
   )
{
   int v;
   int nvars;
   SCIP_VAR** vars;
   SCIP_VAR* origvar;

   assert(scip != NULL);
   assert(probs != NULL);
   assert(nprobs > 0);

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* pricingvar;
      SCIP_Real objvalue;
      assert(GCGvarIsOriginal(vars[v]));
      origvar = SCIPvarGetProbvar(vars[v]);
      pricingvar = GCGoriginalVarGetPricingVar(origvar);
      assert(pricingvar != NULL);

      objvalue = SCIPvarGetObj(origvar);
      /* SCIPinfoMessage(scip, NULL, "%s: %f\n", SCIPvarGetName(origvar), SCIPvarGetObj(origvar));*/
      SCIP_CALL( SCIPchgVarObj(probs[GCGvarGetBlock(origvar)], pricingvar,  SCIPvarGetObj(origvar)) );
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE solveDiagonalBlocks(
   SCIP* scip,
   SCIP_RELAXDATA* relaxdata,
   SCIP_RESULT *result,
   SCIP_Real *lowerbound
   )
{
   int i;
   SCIP_Real objvalue;
   SCIP_Real timelimit;
   SCIP_Real pricingtimelimit;
   SCIP_SOL *newsol;
   SCIP_Bool isfeasible;

   /* set objective of pricing problems to original objective */
   SCIP_CALL( setPricingObjsOriginal(scip, relaxdata->pricingprobs, relaxdata->npricingprobs) );

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   objvalue = 0.0;
   /* solve pricing problems one after the other */

   for( i = 0; i < relaxdata->npricingprobs; ++i)
   {
#ifdef SCIP_DEBUG
      char name[SCIP_MAXSTRLEN];
#endif

      if( relaxdata->pricingprobs[i] == NULL )
         continue;

      SCIPinfoMessage(scip, NULL, "Solving pricing %i.\n", i);
      SCIP_CALL( SCIPsetIntParam(relaxdata->pricingprobs[i], "display/verblevel", SCIP_VERBLEVEL_NONE) );
      /* give the pricing problem 2% more time then the original scip has left */
      if( SCIPgetStage(relaxdata->pricingprobs[i]) > SCIP_STAGE_PROBLEM )
      {
         pricingtimelimit = (timelimit - SCIPgetSolvingTime(scip)) * 1.02 + SCIPgetSolvingTime(relaxdata->pricingprobs[i]);
      }
      else
      {
         pricingtimelimit = (timelimit - SCIPgetSolvingTime(scip)) * 1.02;
      }
      SCIP_CALL( SCIPsetRealParam(relaxdata->pricingprobs[i], "limits/time", pricingtimelimit));

#ifdef SCIP_DEBUG
      SCIPsnprintf(name, SCIP_MAXSTRLEN, "block_%i.lp", i);
      SCIP_CALL( SCIPwriteOrigProblem(relaxdata->pricingprobs[i], name, "lp", FALSE) );
#endif

      SCIP_CALL( SCIPsolve(relaxdata->pricingprobs[i]) );

      switch( SCIPgetStatus(relaxdata->pricingprobs[i]) )
      {
      case SCIP_STATUS_UNBOUNDED:
      case SCIP_STATUS_INFORUNBD:
      case SCIP_STATUS_INFEASIBLE:
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
         break;
      case SCIP_STATUS_BESTSOLLIMIT:
      case SCIP_STATUS_MEMLIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_NODELIMIT:
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_TIMELIMIT:
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
         break;
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_OPTIMAL:
         objvalue += SCIPgetDualbound(relaxdata->pricingprobs[i]);
         break;
      default:
         break;
      }
   }

   /* get solution and glue it together */
   SCIP_CALL( combineSolutions(scip, &newsol, relaxdata->pricingprobs, relaxdata->npricingprobs) );

   /* update lower bound pointer and add solution such that this node will be cut off automatically */
   if(SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE)
      *lowerbound = -objvalue;
   else
      *lowerbound = objvalue;

   SCIP_CALL( SCIPcheckSol(scip, newsol, TRUE, TRUE, TRUE, TRUE, &isfeasible) );
   assert(isfeasible);

   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, &isfeasible) );

   /* maybe add a constraint to the node to indicate that it has been decomposed */

   SCIPinfoMessage(scip, NULL, "We need code for this situation here!\n");

   *result = SCIP_SUCCESS;
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

   /* free master problem */
   if( relaxdata->masterprob != NULL )
      SCIP_CALL( SCIPfree(&(relaxdata->masterprob)) );

   SCIPfreeMemory(scip, &relaxdata);

   return SCIP_OKAY;
}



/** initialization method of relaxator (called after problem was transformed) */
/*static
SCIP_DECL_RELAXINIT(relaxInitGcg)
{
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);
   assert(
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return SCIP_OKAY;
   }*/
#define relaxInitGcg NULL


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
   SCIP_VAR** vars;
   SCIP_RELAXDATA* relaxdata;
   int i;
   int nvars;

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

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "GCG                : Performing Dantzig-Wolfe with %d blocks.\n", relaxdata->npricingprobs);

   for( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      if( relaxdata->convconss[i] != NULL)
      {
         SCIP_CALL( SCIPtransformCons(masterprob, relaxdata->convconss[i], &(relaxdata->convconss[i])) );
      }
   }

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);
   /* transform the linking constraints */
   for( i = 0; i < nvars; ++i)
   {
      int j;
      assert(GCGvarIsOriginal(vars[i]));

      if( GCGvarIsLinking(vars[i]) )
      {
         SCIP_CONS** linkconss;
         linkconss = GCGlinkingVarGetLinkingConss(vars[i]);
         for( j = 0;j < relaxdata->npricingprobs; ++j )
         {
            SCIP_CONS* tempcons;
            if( linkconss[j] != NULL )
            {
               SCIP_CALL( SCIPtransformCons(masterprob, linkconss[j], &(tempcons)) );
               GCGlinkingVarSetLinkingCons(vars[i], tempcons, j);
            }
         }
      }
   }

   if( SCIPfindConshdlr(scip, "connected") != NULL )
   {
      relaxdata->hasblockdetection = TRUE;
      SCIPdebugMessage("Block detection code present.\n");
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
   SCIPfreeMemoryArrayNull(scip, &(relaxdata->convconss));

   /* free master problem */
   SCIP_CALL( SCIPfree(&(relaxdata->masterprob)) );

   /* free pricing problems */
   for( i = relaxdata->npricingprobs - 1; i >= 0 ; i-- )
   {
      SCIP_CALL( SCIPfreeTransform(relaxdata->pricingprobs[i]) );
      SCIP_CALL( SCIPfree(&(relaxdata->pricingprobs[i])) );
   }
   SCIPfreeMemoryArrayNull(scip, &(relaxdata->pricingprobs));
   SCIPfreeMemoryArrayNull(scip, &(relaxdata->blockrepresentative));
   SCIPfreeMemoryArrayNull(scip, &(relaxdata->nblocksidentical));

   /* free solutions */
   if( relaxdata->currentorigsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &relaxdata->currentorigsol) );
   }
   if( relaxdata->storedorigsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &relaxdata->storedorigsol) );
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

   /* only solve the relaxation if it was not yet solved at the current node */
   if( SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) != relaxdata->lastsolvednodenr )
   {
      if( SCIPgetBestSol(scip) != NULL && SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == 1 )
      {
         relaxdata->origprimalsol = SCIPgetBestSol(scip);
      }
      /* increase the node limit for the master problem by 1 */
      SCIP_CALL( SCIPgetLongintParam(masterprob, "limits/nodes", &oldnnodes) );
      SCIP_CALL( SCIPsetLongintParam(masterprob, "limits/nodes",
            ( SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) ? 1 : oldnnodes+1)) );


      /* loop to solve the master problem, this is a workaround and does not fix any problem */
      while( !SCIPisStopped(scip))
      {
         double mastertimelimit = SCIPinfinity(scip);
         SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
         if( !SCIPisInfinity(scip, timelimit) )
         {

            /* give the master 2% more time then the original scip has left */
            mastertimelimit = (timelimit - SCIPgetSolvingTime(scip)) * 1.02 + SCIPgetSolvingTime(masterprob);
            SCIP_CALL( SCIPsetRealParam(masterprob, "limits/time", mastertimelimit) );

            SCIPdebugMessage("Orig left: %f, limit for master %f, left %f\n",
                  timelimit - SCIPgetSolvingTime(scip),
                  mastertimelimit,
                  mastertimelimit - SCIPgetSolvingTime(masterprob));
         }

         /* if we have a blockdetection, see whether the node is block diagonal */

         if( relaxdata->hasblockdetection && SCIPisMatrixBlockDiagonal(scip) )
         {
            SCIP_CALL( solveDiagonalBlocks(scip, relaxdata, result, lowerbound) );
            if( *result == SCIP_SUCCESS)
               return SCIP_OKAY;
         }
         /* We are solving the masterproblem regularly */
         else
         {
            SCIP_CALL( SCIPsolve(masterprob) );
         }


         if(SCIPgetStatus(masterprob) != SCIP_STATUS_TIMELIMIT)
         {
            break;
         }

         if( !SCIPisInfinity(scip, timelimit) )
            SCIPinfoMessage(scip, NULL, "Masterprob was to short, extending time by %f.\n", mastertimelimit - SCIPgetSolvingTime(masterprob));
      }
      if(SCIPgetStatus(masterprob) == SCIP_STATUS_TIMELIMIT && SCIPisStopped(scip))
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }

      /* set the lower bound pointer */
      if( SCIPgetStage(masterprob) == SCIP_STAGE_SOLVING )
         *lowerbound = SCIPgetLocalDualbound(masterprob);
      else
      {
         SCIPdebugMessage("Stage: %d\n", SCIPgetStage(masterprob));
         assert(SCIPgetBestSol(masterprob) != NULL || SCIPgetStatus(masterprob) == SCIP_STATUS_INFEASIBLE);
         if( SCIPgetStatus(masterprob) == SCIP_STATUS_OPTIMAL )
            *lowerbound = SCIPgetSolOrigObj(masterprob, SCIPgetBestSol(masterprob));
         else if( SCIPgetStatus(masterprob) == SCIP_STATUS_INFEASIBLE )
         {
            double tilim;
            SCIP_CALL( SCIPgetRealParam(masterprob, "limits/time", &tilim) );
            if(tilim-SCIPgetSolvingTime(masterprob) < 0)
            {
               *result = SCIP_DIDNOTRUN;
               return SCIP_OKAY;
            }
            *lowerbound = SCIPinfinity(scip);
         }
      }

      SCIPdebugMessage("Update lower bound (value = %"SCIP_REAL_FORMAT").\n", *lowerbound);
   }

   /* transform the current solution of the master problem to the original space and save it */
   SCIPdebugMessage("Update current sol.\n");
   SCIP_CALL( GCGrelaxUpdateCurrentSol(scip, &feasible) );

   if( GCGconsOrigbranchGetBranchrule(GCGconsOrigbranchGetActiveCons(scip)) != NULL
      && SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) != relaxdata->lastsolvednodenr )
   {
      SCIP_CALL( GCGrelaxBranchMasterSolved(scip, GCGconsOrigbranchGetBranchrule(GCGconsOrigbranchGetActiveCons(scip) ),
            GCGconsOrigbranchGetBranchdata(GCGconsOrigbranchGetActiveCons(scip)), *lowerbound) );
   }

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

   relaxdata->blockrepresentative = NULL;
   relaxdata->convconss = NULL;
   relaxdata->hashorig2origvar = NULL;
   relaxdata->hashorig2pricingvar = NULL;
   relaxdata->lastsolvednodenr = 0;

   relaxdata->linearmasterconss = NULL;
   relaxdata->masterconss = NULL;

   relaxdata->npricingprobs = -1;
   relaxdata->pricingprobs = NULL;
   relaxdata->nrelpricingprobs = 0;
   relaxdata->currentorigsol = NULL;
   relaxdata->storedorigsol = NULL;
   relaxdata->origprimalsol = NULL;
   relaxdata->masterprob = NULL;
   relaxdata->nblocksidentical = NULL;

   relaxdata->lastmastersol = NULL;
   relaxdata->lastmasterlpiters = 0;
   relaxdata->markedmasterconss = NULL;
   relaxdata->masterinprobing = FALSE;

   relaxdata->nbranchrules = 0;
   relaxdata->branchrules = NULL;

   relaxdata->nlinkingvars = 0;
   relaxdata->nvarlinkconss = 0;


   /* include relaxator */
   SCIP_CALL( SCIPincludeRelax(scip, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, relaxCopyGcg, relaxFreeGcg, relaxInitGcg,
         relaxExitGcg, relaxInitsolGcg, relaxExitsolGcg, relaxExecGcg, relaxdata) );

   /* inform the main scip, that no LPs should be solved */
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );

   /* Disable restarts */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );

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


/** creates a variable in a pricing problem corresponding to the given original variable (belonging to exactly one block) */
SCIP_RETCODE GCGrelaxCreatePricingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar             /**< corresponding variable in the original program */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_VAR* var;
   int pricingprobnr;

   assert(scip != NULL);
   assert(origvar != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   pricingprobnr = GCGvarGetBlock(origvar);

   SCIP_CALL( GCGoriginalVarCreatePricingVar(relaxdata->pricingprobs[pricingprobnr], origvar, &var) );
   assert(var != NULL);

   GCGoriginalVarSetPricingVar(origvar, var);
   SCIP_CALL( SCIPaddVar(relaxdata->pricingprobs[pricingprobnr], var) );

   /* because the variable was added to the problem,
    * it is captured by SCIP and we can safely release it right now
    */
   SCIP_CALL( SCIPreleaseVar(relaxdata->pricingprobs[pricingprobnr], &var) );

   return SCIP_OKAY;
}

/** creates a variable in each of the pricing problems linked by given original variable */
SCIP_RETCODE GCGrelaxCreateLinkingPricingVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar             /**< corresponding linking variable in the original program */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_VAR* var;
   int pricingprobnr;
   SCIP_CONS* linkcons;
   SCIP_VAR** pricingvars;

   assert(scip != NULL);
   assert(origvar != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* get variable data of the original variable */
   assert(GCGvarIsOriginal(origvar));
   assert(GCGvarIsLinking(origvar));
   pricingvars = GCGlinkingVarGetPricingVars(origvar);

   for( pricingprobnr = 0; pricingprobnr < relaxdata->npricingprobs; pricingprobnr++ )
   {
      if( pricingvars[pricingprobnr] == NULL )
         continue;

      SCIP_CALL( GCGlinkingVarCreatePricingVar(relaxdata->masterprob,
            relaxdata->pricingprobs[pricingprobnr], pricingprobnr, origvar, &var, &linkcons) );
      GCGlinkingVarSetPricingVar(origvar, pricingprobnr, var);

      SCIP_CALL( SCIPaddVar(relaxdata->pricingprobs[pricingprobnr], var) );
      SCIP_CALL( SCIPaddCons(relaxdata->masterprob, linkcons) );
      relaxdata->nvarlinkconss++;

      /* because the variable was added to the problem,
       * it is captured by SCIP and we can safely release it right now
       */
      SCIP_CALL( SCIPreleaseVar(relaxdata->pricingprobs[pricingprobnr], &var) );
   }

   return SCIP_OKAY;
}

/** transforms a constraint of the original problem into the master variable space
 *  and stores information about the constraints in the variable */
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
      SCIP_CALL( GCGoriginalVarAddCoef(scip, consvars[v], consvals[v], mastercons) );
   }

   /* add master variables to the corresponding master constraint */
   for( v = 0; v < nmastervars; v++ )
   {
      SCIP_VAR** origvars;
      SCIP_Real* origvals;
      int norigvars;
      coef = 0;

      origvars = GCGmasterVarGetOrigvars(mastervars[v]);
      norigvars = GCGmasterVarGetNOrigvars(mastervars[v]);
      origvals = GCGmasterVarGetOrigvals(mastervars[v]);

      for( i = 0; i < norigvars; i++ )
         for( j = 0; j < nconsvars; j++ )
            if( consvars[j] == origvars[i] )
               coef += consvals[j] * origvals[i];

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

/**  prints the given variable: name, type (original, master or pricing) block number,
 * and the list of all variables related to the given variable
 */
void GCGrelaxPrintVar(
   SCIP_VAR*             var                 /**< variable that should be printed */
   )
{
   int i;
   int blocknr;
   assert(GCGvarIsOriginal(var) || GCGvarIsMaster(var) || GCGvarIsPricing(var));

   blocknr = GCGvarGetBlock(var);

   if( GCGvarIsOriginal(var) )
   {
      SCIP_VAR** mastervars;
      SCIP_Real* mastervals;
      int  nmastervars;

      if( GCGvarIsLinking(var) )
      {
         SCIP_VAR** pricingvars;
         int nblocks;
         int j;
         pricingvars = GCGlinkingVarGetPricingVars(var);
         nblocks = GCGlinkingVarGetNBlocks(var);
         printf("Variable %s (linking): %d block%s (", SCIPvarGetName(var), nblocks, nblocks == 1 ? "":"s" );
         for( i = 0, j = 0; j < nblocks; ++i)  /*lint --e{440}*/
         {
            if( pricingvars[i] != NULL )
            {
               printf("%d ", i);
               ++j;
            }
         }
         printf(")\n");
      }
      else
      {
         printf("Variable %s (original): block %d\n", SCIPvarGetName(var), blocknr);
      }

      mastervars = GCGoriginalVarGetMastervars(var);
      mastervals = GCGoriginalVarGetMastervals(var);
      nmastervars = GCGoriginalVarGetNMastervars(var);
      printf("mastervars:");
      for( i = 0; i < nmastervars-1; i++ )
      {
         printf("%s (%g), ", SCIPvarGetName(mastervars[i]), mastervals[i]);
      }
      printf("%s (%g)\n", SCIPvarGetName(mastervars[nmastervars-1]), mastervals[nmastervars-1]);
   }
   else if( GCGvarIsPricing(var) )
   {
      SCIP_VAR** origvars;
      int  norigvars;

      origvars = GCGpricingVarGetOrigvars(var);
      norigvars = GCGpricingVarGetNOrigvars(var);

      printf("Variable %s (pricing): block %d\n", SCIPvarGetName(var), blocknr);
      printf("origvars:");
      for( i = 0; i < norigvars-1; i++ )
      {
         printf("%s, ", SCIPvarGetName(origvars[i]));
      }
      printf("%s\n", SCIPvarGetName(origvars[norigvars-1]));
   }
   else if( GCGvarIsMaster(var) )
   {
      SCIP_VAR** origvars;
      int  norigvars;
      SCIP_Real* origvals;

      origvars = GCGmasterVarGetOrigvars(var);
      norigvars = GCGmasterVarGetNOrigvars(var);
      origvals = GCGmasterVarGetOrigvals(var);
      printf("Variable %s (master): block %d\n", SCIPvarGetName(var), blocknr);
      printf("origvars:");
      for( i = 0; i < norigvars-1; i++ )
      {
         printf("%s (%g), ", SCIPvarGetName(origvars[i]), origvals[i]);
      }
      printf("%s (%g)\n", SCIPvarGetName(origvars[norigvars-1]), origvals[norigvars-1]);
   }
}

/* sets the number of the block, the given original variable belongs to */
SCIP_RETCODE GCGrelaxSetOriginalVarBlockNr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to set the block number for */
   int                   newblock            /**< number of the block, the variable belongs to */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   int blocknr;

   assert(scip != NULL);
   assert(var != NULL);
   assert(newblock >= 0);
   assert(SCIPvarIsOriginal(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   blocknr = GCGvarGetBlock(var);
   assert(GCGvarIsOriginal(var));

   assert(GCGrelaxGetNPricingprobs(scip) > 0);
   assert(newblock < GCGrelaxGetNPricingprobs(scip));
   assert(blocknr >= -2 && blocknr < GCGrelaxGetNPricingprobs(scip));

   /* var belongs to no block so far, just set the new block number */
   if( blocknr == -1 )
      GCGvarSetBlock(var, newblock);

   /* if var already belongs to another block, it is a linking variable */
   else if ( blocknr != newblock )
   {
      if(!GCGvarIsLinking(var))
         relaxdata->nlinkingvars++;

      SCIP_CALL( GCGoriginalVarAddBlock(scip, var, newblock) );
      assert(GCGisLinkingVarInBlock(var, newblock));
   }
   blocknr = GCGvarGetBlock(var);
   assert(blocknr == -2 || blocknr == newblock);

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
      int nconss;
      nconss = SCIPgetNConss(scip);
      SCIP_CALL( SCIPallocMemoryArray(scip, &(relaxdata->markedmasterconss), nconss) );
      relaxdata->nmarkedmasterconss = 0;
   }
   assert(relaxdata->nmarkedmasterconss < SCIPgetNConss(scip));

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

   assert(relaxdata->npricingprobs >= -1);
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

/**
 *  for a given block, return the block by which it is represented
 */
int GCGrelaxGetBlockRepresentative(
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

   return relaxdata->blockrepresentative[pricingprobnr];
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

   /* start probing in the master problem */
   SCIP_CALL( SCIPstartProbing(masterscip) );

   relaxdata->masterinprobing = TRUE;

   /* remember the current original solution */
   assert(relaxdata->storedorigsol == NULL);
   if( relaxdata->currentorigsol != NULL )
      SCIP_CALL( SCIPcreateSolCopy(scip, &relaxdata->storedorigsol, relaxdata->currentorigsol) );

   return SCIP_OKAY;
}


/** for a probing node in the original problem, create a corresponding probing node in the master problem,
 *  propagate domains and solve the LP with or without pricing. */
static
SCIP_RETCODE performProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint          maxlpiterations,    /**< maximum number of lp iterations allowed */
   int                   maxpricerounds,     /**< maximum number of pricing rounds allowed */
   SCIP_Bool             usepricing,         /**< should the LP be solved with or without pricing? */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of performed LP iterations (or NULL) */
   int*                  npricerounds,       /**< pointer to store the number of performed pricing rounds (or NULL) */
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
   SCIP_Longint oldnlpiters;
   int oldpricerounds;
   SCIP_Longint nodelimit;

   assert(scip != NULL);

   /* get the relaxator */
   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   /* get the relaxator data */
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->masterinprobing);

   /* get master problem */
   masterscip = relaxdata->masterprob;
   assert(masterscip != NULL);

   /* create probing node in the master problem */
   SCIP_CALL( SCIPnewProbingNode(masterscip) );

   /* create master constraint that captures the branching decision in the original instance */
   mprobingnode = SCIPgetCurrentNode(masterscip);
   assert(GCGconsMasterbranchGetActiveCons(masterscip) != NULL);
   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &mprobingcons, mprobingnode,
         GCGconsMasterbranchGetActiveCons(masterscip)) );
   SCIP_CALL( SCIPaddConsNode(masterscip, mprobingnode, mprobingcons, NULL) );
   SCIP_CALL( SCIPreleaseCons(scip, &mprobingcons) );

   /* increase node limit for the master problem by 1 */
   SCIP_CALL( SCIPgetLongintParam(masterscip, "limits/nodes", &nodelimit) );
   SCIP_CALL( SCIPsetLongintParam(masterscip, "limits/nodes", nodelimit + 1) );

   /* propagate */
   SCIP_CALL( SCIPpropagateProbing(masterscip, -1, cutoff, NULL) );
   assert(!(*cutoff));

   /* remember LP iterations and pricing rounds before LP solving */
   oldnlpiters = SCIPgetNLPIterations(masterscip);
   oldpricerounds = SCIPgetNPriceRounds(masterscip);

   /* solve the probing LP */
   if( usepricing )
   {
      /* LP iterations are unlimited when probing LP is solved with pricing */
      assert(maxlpiterations == -1);
      SCIP_CALL( SCIPsolveProbingLPWithPricing(masterscip, FALSE/* pretendroot */, TRUE /*displayinfo*/,
            maxpricerounds, lperror) );
   }
   else
   {
      assert(maxpricerounds == 0);
      SCIP_CALL( SCIPsolveProbingLP(masterscip, maxlpiterations, lperror) );
   }
   lpsolstat = SCIPgetLPSolstat(masterscip);

   /* reset the node limit */
   SCIP_CALL( SCIPsetLongintParam(masterscip, "limits/nodes", nodelimit) );

   /* calculate number of LP iterations and pricing rounds performed */
   if( nlpiterations != NULL )
      *nlpiterations = SCIPgetNLPIterations(masterscip) - oldnlpiters;
   if( npricerounds != NULL )
      *npricerounds = SCIPgetNPriceRounds(masterscip) - oldpricerounds;

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
      SCIPdebugMessage("something went wrong, an lp error occured\n");
   }

   return SCIP_OKAY;
}


/** for a probing node in the original problem, create a corresponding probing node in the master problem,
 *  propagate domains and solve the LP without pricing. */
SCIP_RETCODE GCGrelaxPerformProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint          maxlpiterations,    /**< maximum number of lp iterations allowed */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of performed LP iterations (or NULL) */
   SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the probing direction is infeasible */
   SCIP_Bool*            feasible            /**< pointer to store whether the probing solution is feasible */
   )
{
   SCIP_CALL( performProbing(scip, maxlpiterations, 0, FALSE, nlpiterations,
         NULL, lpobjvalue, lpsolved, lperror, cutoff, feasible) );

   return SCIP_OKAY;
}


/** for a probing node in the original problem, create a corresponding probing node in the master problem,
 *  propagate domains and solve the LP with pricing. */
SCIP_RETCODE GCGrelaxPerformProbingWithPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   maxpricerounds,     /**< maximum number of pricing rounds allowed */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of performed LP iterations (or NULL) */
   int*                  npricerounds,       /**< pointer to store the number of performed pricing rounds (or NULL) */
   SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the probing direction is infeasible */
   SCIP_Bool*            feasible            /**< pointer to store whether the probing solution is feasible */
   )
{
   SCIP_CALL( performProbing(scip, -1, maxpricerounds, TRUE, nlpiterations,
         npricerounds, lpobjvalue, lpsolved, lperror, cutoff, feasible) );

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

   SCIP_VAR** vars;
   int nvars;

   int i;

   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->masterinprobing);

   masterscip = relaxdata->masterprob;
   assert(masterscip != NULL);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(vars != NULL);
   assert(nvars >= 0);

   SCIP_CALL( SCIPendProbing(masterscip) );

   relaxdata->masterinprobing = FALSE;

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
   if( relaxdata->currentorigsol != NULL )
   {
      SCIPdebugMessage("Freeing previous solution origsol\n");
      SCIP_CALL( SCIPfreeSol(scip, &(relaxdata->currentorigsol)) );
   }
   SCIPclearExternBranchCands(scip);

   if( relaxdata->storedorigsol != NULL )
   {
      SCIP_CALL( SCIPcreateSol(scip, &relaxdata->currentorigsol, NULL) );
      SCIP_CALL( SCIPsetRelaxSolValsSol(scip, relaxdata->storedorigsol) );

      for( i = 0; i < nvars; i++ )
      {
         SCIP_VAR* var;
         SCIP_Real solval;

         var = vars[i];
         solval = SCIPgetSolVal(scip, relaxdata->storedorigsol, var);

         SCIPsetSolVal(scip, relaxdata->currentorigsol, var, solval);

         if( SCIPvarGetType(var) <= SCIP_VARTYPE_INTEGER && !SCIPisFeasIntegral(scip, solval) )
         {
            assert(!SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, solval - SCIPfloor(scip, solval), solval) );
         }
      }
      assert(SCIPisFeasEQ(scip, SCIPgetRelaxSolObj(scip), SCIPgetSolTransObj(scip, relaxdata->currentorigsol)));

      SCIP_CALL( SCIPfreeSol(scip, &relaxdata->storedorigsol) );
   }

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
      SCIPdebugMessage("Freeing previous solution origsol\n");
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
      {
         SCIPdebugMessage("Masterproblem still solving, mastersol = NULL\n");
         mastersol = NULL;
      }
      else if( SCIPgetStage(relaxdata->masterprob) == SCIP_STAGE_SOLVED )
      {
         mastersol = SCIPgetBestSol(relaxdata->masterprob);
         if( mastersol == NULL )
         {
            SCIPdebugMessage("Masterproblem solved, no master sol present\n");
            return SCIP_OKAY;
         }
         SCIPdebugMessage("Masterproblem solved, mastersol = %pd\n", mastersol);
      }
      else
      {
         SCIPdebugMessage("stage in master not solving and not solved!\n");
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
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, &stored) );
#else
      SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, TRUE, TRUE, TRUE, &stored) );
#endif
      if( !stored )
      {

         SCIP_CALL( SCIPcheckSolOrig(scip, newsol, &stored, TRUE, TRUE) );
      }
      /** @todo: Martin does not see why the solution has to be accepted, numerics might bite us, so the transformation might fail.
       *  Remedy could be: Round the values or propagate changes or call a heuristic to fix it.
       */
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
      /** @todo: Martin will disable that here, because at the current stage, it does not have to be true!
       *       assert(stored);
       */
      if(stored)
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
   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(origvars != NULL);
   assert(origvals != NULL);
   assert(mastervars != NULL);
   assert(mastervals != NULL);
   assert(nmastervars >= 0);

   /* set all values to 0 initially */
   for( i = 0; i < nmastervars; i++ )
      mastervals[i] = 0.0;

   /* iterate over all original variables */
   for( i = 0; i < norigvars; i++ )
   {
      SCIP_VAR** varmastervars;
      SCIP_Real* varmastervals;
      int blocknr;

      assert(GCGvarIsOriginal(origvars[i]));
      varmastervars = GCGoriginalVarGetMastervars(origvars[i]);
      varmastervals = GCGoriginalVarGetMastervals(origvars[i]);
      blocknr = GCGvarGetBlock(origvars[i]);

      /* variable belongs to no block (or is a linking variable), so it was transferred directly to the master problem,
       * hence, we transfer the value directly to the corresponding master variabe
       */
      if( blocknr < 0 )
      {
         assert(blocknr == -1 || blocknr == -2);
         for( k = 0; k < nmastervars; k++ )
         {
            assert(!SCIPvarIsTransformedOrigvar(mastervars[k]));
            if( mastervars[k] == varmastervars[0])
            {
               assert(!SCIPvarIsTransformedOrigvar(varmastervars[0]));
               mastervals[k] += (varmastervals[0] * origvals[i]);
               break;
            }
         }
         assert(k < nmastervars);
      }
      /* variable belongs to exactly one block, so we have to look at all master variables and increase their values
       * if they contain the original variable
       */
      else
      {
         SCIP_VAR* pricingvar;
         SCIP_VAR* origvar;
         SCIP_VAR** curmastervars;
         SCIP_Real* curmastervals;
         int ncurmastervars;

         pricingvar = GCGoriginalVarGetPricingVar(origvars[i]);
         assert(GCGvarIsPricing(pricingvar));

         origvar = GCGpricingVarGetOriginalVar(pricingvar);
         assert(GCGvarIsOriginal(origvar));
         curmastervars = GCGoriginalVarGetMastervars(origvar);
         curmastervals = GCGoriginalVarGetMastervals(origvar);
         ncurmastervars = GCGoriginalVarGetNMastervars(origvar);

         for( j = 0; j < ncurmastervars; j++ )
         {
            for( k = 0; k < nmastervars; k++ )
               if( mastervars[k] == curmastervars[j] )
               {
                  mastervals[k] += (curmastervals[j] * origvals[i]);
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
   int* blocknrs;
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
   SCIP_CALL( SCIPallocBufferArray(scip, &blocknrs, relaxdata->npricingprobs) );

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
      blocknrs[i] = 0;
   }

   /* loop over all given master variables */
   for( i = 0; i < nmastervars; i++ )
   {
      SCIP_VAR** origvars;
      int norigvars;
      SCIP_Real* origvals;
      SCIP_Bool isray;
      int blocknr;

      origvars = GCGmasterVarGetOrigvars(mastervars[i]);
      norigvars = GCGmasterVarGetNOrigvars(mastervars[i]);
      origvals = GCGmasterVarGetOrigvals(mastervars[i]);
      blocknr = GCGvarGetBlock(mastervars[i]);
      isray = GCGmasterVarIsRay(mastervars[i]);

      assert(GCGvarIsMaster(mastervars[i]));
      assert(!SCIPisFeasNegative(scip, mastervals[i]));

      /* TODO: handle infinite master solution values */
      assert(!SCIPisInfinity(scip, mastervals[i]));

      /* first of all, handle variables representing rays */
      if( isray )
      {
         assert(blocknr >= 0);
         /* we also want to take into account variables representing rays, that have a small value (between normal and feas eps),
          * so we do no feas comparison here */
         if( SCIPisPositive(scip, mastervals[i]) )
         {
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < norigvars; j++ )
            {
               if(SCIPisZero(scip, origvals[j]))
                  break;

               assert(!SCIPisZero(scip, origvals[j]));

               /* the original variable is a linking variable: just transfer the solution value of the direct copy (this is done later) */
               if( GCGvarIsLinking(origvars[j]) )
                  continue;

//               SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origvars[j]), origvals[j] * mastervals[i], SCIPvarGetName(mastervars[i]));
               /* increase the corresponding value */
               SCIP_CALL( SCIPincSolVal(scip, *origsol, origvars[j], origvals[j] * mastervals[i]) );
            }
         }
         mastervals[i] = 0.0;
         continue;
      }

      /* handle the variables with value >= 1 to get integral values in original solution */
      while( SCIPisFeasGE(scip, mastervals[i], 1.0) )
      {
         /* variable was directly transferred to the master problem (only in linking conss or linking variable) */
         /* TODO: this may be the wrong place for this case, handle it before the while loop
          * and remove the similar case in the next while loop */
         if( blocknr == -1 )
         {
            assert(norigvars == 1);
            assert(origvals[0] == 1.0);

            /* increase the corresponding value */
//            SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origvars[0]), origvals[0] * mastervals[i],  SCIPvarGetName(mastervars[i]));
            SCIP_CALL( SCIPincSolVal(scip, *origsol, origvars[0], origvals[0] * mastervals[i]) );
            mastervals[i] = 0.0;
         }
         else
         {
            assert(blocknr >= 0);
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < norigvars; j++ )
            {
               SCIP_VAR* pricingvar;
               int norigpricingvars;
               SCIP_VAR** origpricingvars;
               if(SCIPisZero(scip, origvals[j]))
                  break;
               assert(!SCIPisZero(scip, origvals[j]));

               /* the original variable is a linking variable: just transfer the solution value of the direct copy (this is done above) */
               if( GCGvarIsLinking(origvars[j]) )
                  continue;

               pricingvar = GCGoriginalVarGetPricingVar(origvars[j]);
               assert(GCGvarIsPricing(pricingvar));

               norigpricingvars = GCGpricingVarGetNOrigvars(pricingvar);
               origpricingvars = GCGpricingVarGetOrigvars(pricingvar);

               /* just in case a variable has a value higher than the number of blocks, it represents */
               if( norigpricingvars <= blocknrs[blocknr] )
               {
//                  SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[norigpricingvars-1]), mastervals[i] * origvals[j], SCIPvarGetName(mastervars[i]));
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, origpricingvars[norigpricingvars-1], mastervals[i] * origvals[j]) );
                  mastervals[i] = 1.0;
               }
               /* this should be default */
               else
               {
//                  SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[blocknrs[blocknr]]), origvals[j], SCIPvarGetName(mastervars[i]) );
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, origpricingvars[blocknrs[blocknr]], origvals[j]) );
               }
            }
            mastervals[i] = mastervals[i] - 1.0;
            blocknrs[blocknr]++;
         }
      }
   }

   /* loop over all given master variables */
   for( i = 0; i < nmastervars; i++ )
   {
      SCIP_VAR** origvars;
      int norigvars;
      SCIP_Real* origvals;
      int blocknr;

      origvars = GCGmasterVarGetOrigvars(mastervars[i]);
      norigvars = GCGmasterVarGetNOrigvars(mastervars[i]);
      origvals = GCGmasterVarGetOrigvals(mastervars[i]);
      blocknr = GCGvarGetBlock(mastervars[i]);

      if( SCIPisFeasZero(scip, mastervals[i]) )
      {
         continue;
      }
      assert(SCIPisFeasGE(scip, mastervals[i], 0.0) && SCIPisFeasLT(scip, mastervals[i], 1.0));

      while( SCIPisFeasPositive(scip, mastervals[i]) )
      {
         assert(GCGvarIsMaster(mastervars[i]));
         assert(!GCGmasterVarIsRay(mastervars[i]));

         if( blocknr == -1 )
         {
            assert(norigvars == 1);
            assert(origvals[0] == 1.0);

//            SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origvars[0]), origvals[0] * mastervals[i], SCIPvarGetName(mastervars[i]) );
            /* increase the corresponding value */
            SCIP_CALL( SCIPincSolVal(scip, *origsol, origvars[0], origvals[0] * mastervals[i]) );
            mastervals[i] = 0.0;
         }
         else
         {
            increaseval = MIN(mastervals[i], 1.0 - blockvalue[blocknr]);
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < norigvars; j++ )
            {
               SCIP_VAR* pricingvar;
               int norigpricingvars;
               SCIP_VAR** origpricingvars;

               if( SCIPisZero(scip, origvals[j]) )
                  continue;

               /* the original variable is a linking variable: just transfer the solution value of the direct copy (this is done above) */
               if( GCGvarIsLinking(origvars[j]) )
                  continue;

               pricingvar = GCGoriginalVarGetPricingVar(origvars[j]);
               assert(GCGvarIsPricing(pricingvar));

               norigpricingvars = GCGpricingVarGetNOrigvars(pricingvar);
               origpricingvars = GCGpricingVarGetOrigvars(pricingvar);

               if( norigpricingvars <= blocknrs[blocknr] )
               {
                  increaseval = mastervals[i];

//                  SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[norigpricingvars-1]), origvals[j] * increaseval, SCIPvarGetName(mastervars[i]) );
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, origpricingvars[norigpricingvars-1], origvals[j] * increaseval) );
               }
               else
               {
                  /* increase the corresponding value */
//                  SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[blocknrs[blocknr]]), origvals[j] * increaseval, SCIPvarGetName(mastervars[i]) );
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, origpricingvars[blocknrs[blocknr]], origvals[j] * increaseval) );
               }
            }

            mastervals[i] = mastervals[i] - increaseval;
            if( SCIPisFeasZero(scip, mastervals[i]) )
            {
               mastervals[i] = 0.0;
            }
            blockvalue[blocknr] += increaseval;

            /* if the value assigned to the block is equal to 1, this block is full and we take the next block */
            if( SCIPisFeasGE(scip, blockvalue[blocknr], 1.0) )
            {
               blockvalue[blocknr] = 0.0;
               blocknrs[blocknr]++;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &mastervals);
   SCIPfreeBufferArray(scip, &blocknrs);
   SCIPfreeBufferArray(scip, &blockvalue);

   return SCIP_OKAY;

}

/* returns the stored primal solution of the original problem  */
SCIP_SOL* GCGrelaxGetOrigPrimalSol(
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

   return relaxdata->origprimalsol;
}

/* sets the stored primal solution of the original problem  */
void GCGrelaxSetOrigPrimalSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< solution */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);

   relax = SCIPfindRelax(scip, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   relaxdata->origprimalsol = sol;
}

