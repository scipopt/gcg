/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_connected.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for connected constraints
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */

#include <assert.h>
#include <string.h>

#include "cons_connected.h"
#include "cons_decomp.h"
#include "scip_misc.h"
#include "scip/clock.h"
#include "pub_decomp.h"
#include "cons_decomp.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "connected"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                          *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA         TRUE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP         TRUE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP


/*
 * Data structures
 */

/** constraint data for connected constraints */
struct SCIP_ConsData
{
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_HASHMAP* constoblock;
   SCIP_HASHMAP* vartoblock;
   SCIP_Bool blockdiagonal;

   DECDECOMP* decdecomp;
   SCIP_CLOCK* clock;
   int nblocks;
   SCIP_Bool enable;
   SCIP_Bool *consismaster;
   SCIP_Bool setppcinmaster;
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/* returns whether the constraint belongs to GCG or not */
static
SCIP_Bool isConsGCGCons(
   SCIP_CONS* cons   /**< constraint to check */
   )
{
   SCIP_CONSHDLR* conshdlr;
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   if( strcmp("origbranch", SCIPconshdlrGetName(conshdlr)) == 0 )
      return TRUE;
   else if( strcmp("masterbranch", SCIPconshdlrGetName(conshdlr)) == 0 )
      return TRUE;

   return FALSE;
}


/* returns true if the constraint should be a master constraint and false otherwise */
static
SCIP_Bool isConsMaster(
   SCIP*      scip,  /**< SCIP data structure */
   SCIP_CONS* cons   /**< constraint to check */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int i;
   int nvars;
   SCIP_Bool relevant = TRUE;
   assert(scip != NULL);
   assert(cons != NULL);

   SCIPdebugMessage("cons %s is ", SCIPconsGetName(cons));

   if(SCIPconsGetType(cons) == setcovering || SCIPconsGetType(cons) == setpartitioning || SCIPconsGetType(cons) == logicor)
   {
      SCIPdebugPrintf("setcov, part or logicor.\n");
      return TRUE;
   }
   nvars = SCIPgetNVarsXXX(scip, cons);
   vars = NULL;
   vals = NULL;
   if( nvars > 0)
   {
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vars, nvars) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vals, nvars) );
      SCIP_CALL_ABORT( SCIPgetVarsXXX(scip, cons, vars, nvars) );
      SCIP_CALL_ABORT( SCIPgetValsXXX(scip, cons, vals, nvars) );
   }

   /* check vars and vals for integrality */
   for(i = 0; i < nvars && relevant; ++i)
   {
      assert(vars != NULL);
      assert(vals != NULL);

      if( !SCIPvarIsIntegral(vars[i]) && !SCIPvarIsBinary(vars[i]) )
      {
         SCIPdebugPrintf("(%s is not integral) ", SCIPvarGetName(vars[i]) );
         relevant = FALSE;
      }
      if( !SCIPisEQ(scip, vals[i], 1.0) )
      {
         SCIPdebugPrintf("(coeff for var %s is %.2f != 1.0) ", SCIPvarGetName(vars[i]), vals[i] );
         relevant = FALSE;
      }
   }

   /* free temporary data  */
   SCIPfreeBufferArrayNull(scip, &vals);
   SCIPfreeBufferArrayNull(scip, &vars);

   SCIPdebugPrintf("%s master\n", relevant ? "in" : "not in");
   return relevant;
}


/** looks for connected components in the constraints in conshdlrdata */
static
SCIP_RETCODE findConnectedComponents(
   SCIP*              scip,         /**< SCIP data structure */
   SCIP_CONSHDLRDATA* conshdlrdata, /**< constraint handler data structure */
   SCIP_RESULT*       result        /**< result pointer to indicate success oder failuer */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int nconss;
   SCIP_VAR** curvars;
   int ncurvars;
   SCIP_CONS* cons;
   SCIP_CONS** conss;

   int i;
   int j;
   int k;
   int tempblock;

   int* blockrepresentative;
   int nextblock;
   int *vartoblock;
   SCIP_HASHMAP *constoblock;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Trying to detect block diagonal matrix.\n");

   /* initialize data structures */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   nextblock = 1; /* start at 1 in order to see whether the hashmap has a key*/

   SCIP_CALL( SCIPallocBufferArray(scip, &vartoblock, nvars+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blockrepresentative, nconss+1) );
   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss+1) );
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->constoblock, SCIPblkmem(scip), nconss) );

   for( i = 0; i < nvars; ++i )
   {
      vartoblock[i] = -1;
   }

   for( i = 0; i < nconss+1; ++i )
   {
      blockrepresentative[i] = -1;
   }

   blockrepresentative[0] = 0;
   blockrepresentative[1] = 1;
   assert(nconss >= 1);

   /* in a first preprocessing step, indicate which constraints should go in the master */

   for ( i = 0; i < nconss; ++i )
   {
      if(conshdlrdata->setppcinmaster)
         conshdlrdata->consismaster[i] = isConsMaster(scip, conss[i]);
      else
         conshdlrdata->consismaster[i] = FALSE;

   }

   /* go through the all constraints */
   for( i = 0; i < nconss; ++i )
   {
      int consblock;
      cons = conss[i];
      assert(cons != NULL);
      if( isConsGCGCons(cons) )
         continue;

      if( conshdlrdata->consismaster[i] )
         continue;

      ncurvars = SCIPgetNVarsXXX(scip, cons);
      curvars = NULL;
      if( ncurvars > 0)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, cons, curvars, ncurvars) );
      }
      assert(ncurvars >= 0);
      assert(ncurvars <= nvars);
      assert(curvars != NULL || ncurvars == 0);

      assert(SCIPhashmapGetImage(constoblock, cons) == NULL);
      if(ncurvars == 0)
         consblock = 1;
      else
         consblock = nextblock;

      /* go through all variables */
      for( j = 0; j < ncurvars; ++j)
      {
         SCIP_VAR* probvar;
         int varindex;
         int varblock;

         assert(curvars != NULL);
         probvar = SCIPvarGetProbvar(curvars[j]);
         assert(probvar != NULL);

         varindex = SCIPvarGetProbindex(probvar);
         assert(varindex >= 0);
         assert(varindex < nvars);

         /** @todo what about deleted variables? */
         varblock = vartoblock[varindex];

         /* if variable is assigned to a block, assign constraint to that block */
         if( varblock != -1 )
         {
            if(consblock == nextblock)
               consblock = varblock;

            /* if variable is assigned to a different block, merge the blocks */
            if( varblock != consblock  )
            {
               /* always take the lower one of both as the representative*/
               if(varblock < consblock)
                  blockrepresentative[consblock] = varblock;
               else
                  blockrepresentative[varblock] = consblock;

               assert(blockrepresentative[consblock] >= 1);
               assert(blockrepresentative[consblock] <= nextblock);
               assert(blockrepresentative[varblock] >= 1);
               assert(blockrepresentative[varblock] <= nextblock);
            }

            /* assign all previous variables of this constraint to this block */
            for (k = j; k >= 0; --k)
            {
               /** @todo what about deleted variables? */
               assert(consblock >= 1);
               assert(consblock <= nextblock);
               vartoblock[SCIPvarGetProbindex(SCIPvarGetProbvar(curvars[k]))] = consblock;
            }
         }
         else
         {
            /* if variable is free, assign it to the new block for this constraint */
            varblock = consblock;
            assert(varblock > 0);
            assert(varblock <= nextblock);
            vartoblock[varindex] = varblock;
         }
      }

      /* if the constraint belongs to a new block, mark it as such */
      if( consblock == nextblock )
      {
         blockrepresentative[consblock] = consblock;
         assert(blockrepresentative[consblock] >= 0);
         assert(blockrepresentative[consblock] <= nextblock);
         ++nextblock;
      }

      SCIPfreeBufferArrayNull(scip, &curvars);
      assert(consblock >= 1);
      assert(consblock <= nextblock);

      /* store the constraint block */
      SCIP_CALL( SCIPhashmapInsert(constoblock, cons, (void*)(size_t)consblock) );
   }

   tempblock = 1;

   SCIPdebugPrintf("Blocks: ");
   /* postprocess blockrepresentatives */
   for( i = 1; i < nextblock; ++i )
   {
      /* forward replace the representatives */
      assert(blockrepresentative[i] >= 0);
      assert(blockrepresentative[i] < nextblock);
      if(blockrepresentative[i] != i)
         blockrepresentative[i] = blockrepresentative[blockrepresentative[i]];
      else
      {
         blockrepresentative[i] = tempblock;
         ++tempblock;
      }
      /* It is crucial that this condition holds */
      assert(blockrepresentative[i] <= i);
      SCIPdebugPrintf("%d ", blockrepresentative[i]);
   }
   SCIPdebugPrintf("\n");

   /* convert temporary data to conshdlrdata */
   for(i = 0; i < nconss; ++i)
   {
      int consblock;

      cons = conss[i];
      if( isConsGCGCons(cons) )
         continue;

      if( conshdlrdata->consismaster[i] )
         continue;

      consblock = (size_t)SCIPhashmapGetImage(constoblock, cons);
      assert(consblock > 0);
      consblock = blockrepresentative[consblock];
      assert(consblock < tempblock);
      SCIP_CALL( SCIPhashmapInsert(conshdlrdata->constoblock, cons, (void*)(size_t)consblock) );
      SCIPdebugMessage("%d %s\n", consblock, SCIPconsGetName(cons));
   }

   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->vartoblock, SCIPblkmem(scip), nvars+1) );

   for( i = 0; i < nvars; ++i)
   {
      int varindex;
      int varblock;
      varindex = SCIPvarGetProbindex(SCIPvarGetProbvar(vars[i]));
      assert(varindex >= 0);
      assert(varindex < nvars);

      assert(vartoblock[varindex] < nextblock);
      if(vartoblock[varindex] < 0)
         continue;

      varblock = blockrepresentative[vartoblock[varindex]];
      assert(varblock == -1 || varblock >= 0);
      if(varblock > 0)
      {
         assert(varblock < tempblock);
         SCIP_CALL( SCIPhashmapInsert(conshdlrdata->vartoblock, SCIPvarGetProbvar(vars[i]),
               (void*)(size_t)(varblock)) );
      }
   }

   /* free method data */
   SCIPfreeBufferArray(scip, &vartoblock);
   SCIPfreeBufferArray(scip, &blockrepresentative);
   SCIPhashmapFree(&constoblock);
   conshdlrdata->nblocks = tempblock-1;

   if(conshdlrdata->nblocks > 1)
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/* copy conshdldata data to decdecomp */
static
SCIP_RETCODE copyToDecdecomp(
   SCIP*              scip,         /**< SCIP data structure */
   SCIP_CONSHDLRDATA* conshdlrdata, /**< constraint handler data structure */
   DECDECOMP*         decdecomp     /**< decdecomp data structure */
   )
{
   SCIP_CONS** conss;
   int nconss;
   SCIP_VAR** vars;
   int nvars;
   int i;

   SCIP_CONS*** subscipconss;
   int* nsubscipconss;
   SCIP_CONS** linkingconss;
   int nlinkingconss;
   SCIP_VAR*** subscipvars;
   int* nsubscipvars;
   int nblocks;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(decdecomp != NULL);

   assert(DECdecdecompGetType(decdecomp) == DEC_DECTYPE_UNKNOWN);

   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);
   nlinkingconss = 0;
   nblocks = conshdlrdata->nblocks;

   SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipvars, nblocks) ) ;
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linkingconss, nconss) );

   for( i = 0; i < nblocks; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars[i], nvars) );
      nsubscipvars[i] = 0;
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss[i], nconss) );
      nsubscipconss[i] = 0;
   }

   DECdecdecompSetNBlocks(decdecomp, nblocks);
   DECdecdecompSetConstoblock(decdecomp, conshdlrdata->constoblock);
   DECdecdecompSetVartoblock(decdecomp, conshdlrdata->vartoblock);

   for( i = 0; i < nconss; ++i)
   {
      size_t consblock;
      if( isConsGCGCons(conss[i]) )
         continue;
      if( conshdlrdata->consismaster[i] )
      {
         linkingconss[nlinkingconss] = conss[i];
         ++(nlinkingconss);
         continue;
      }

      consblock = (size_t) SCIPhashmapGetImage(conshdlrdata->constoblock, conss[i]);
      assert(consblock > 0);
      assert(nblocks >= 0);
      assert(consblock <= (size_t)nblocks);

      subscipconss[consblock-1][nsubscipconss[consblock-1]] = conss[i];
      ++(nsubscipconss[consblock-1]);
   }

   for( i = 0; i < nvars; ++i)
   {
      size_t varblock;
      varblock = (size_t) SCIPhashmapGetImage(conshdlrdata->vartoblock, vars[i]);
      assert(varblock > 0);
      assert(nblocks >= 0);
      assert(varblock <= (size_t) nblocks);

      subscipvars[varblock-1][nsubscipvars[varblock-1]] = vars[i];
      ++(nsubscipvars[varblock-1]);
   }

   if(nlinkingconss > 0)
   {
      SCIP_CALL( DECdecdecompSetLinkingconss(scip, decdecomp, linkingconss, nlinkingconss) );
      DECdecdecompSetType(decdecomp, DEC_DECTYPE_BORDERED);
   }
   else
   {
      DECdecdecompSetType(decdecomp, DEC_DECTYPE_DIAGONAL);
   }
   SCIP_CALL( DECdecdecompSetSubscipconss(scip, decdecomp, subscipconss, nsubscipconss) );
   SCIP_CALL( DECdecdecompSetSubscipvars(scip, decdecomp, subscipvars, nsubscipvars) );


   for( i = 0; i < conshdlrdata->nblocks; ++i)
   {
      SCIPfreeBufferArray(scip, &subscipvars[i]);
      SCIPfreeBufferArray(scip, &subscipconss[i]);
   }

   SCIPfreeBufferArray(scip, &subscipvars);
   SCIPfreeBufferArray(scip, &subscipconss);
   SCIPfreeBufferArray(scip, &nsubscipvars);
   SCIPfreeBufferArray(scip, &nsubscipconss);
   SCIPfreeBufferArray(scip, &linkingconss);

   conshdlrdata->vartoblock = NULL;
   conshdlrdata->constoblock = NULL;

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/* unneeded methods */

#define conshdlrCopyConnected NULL
#define consDeleteConnected NULL
#define consTransConnected NULL
#define consInitlpConnected NULL
#define consSepalpConnected NULL
#define consSepasolConnected NULL
#define consPropConnected NULL
#define consPresolConnected NULL
#define consRespropConnected NULL
#define consActiveConnected NULL
#define consDeactiveConnected NULL
#define consEnableConnected NULL
#define consDisableConnected NULL
#define consDelvarConnected NULL
#define consPrintConnected NULL
#define consCopyConnected NULL
#define consParseConnected NULL
#define consInitConnected NULL
#define consExitConnected NULL
#define consInitpreConnected NULL
#define consExitpreConnected NULL


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeConnected)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA *conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, &conshdlrdata);
   return SCIP_OKAY;

}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolConnected)
{  /*lint --e{715}*/

   SCIP_CONSHDLRDATA *conshdlrdata;
   SCIP_RESULT result;
   int runs;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( !conshdlrdata->enable )
      return SCIP_OKAY;

   /* apparently, there is a structure, which means we don't try to detect one */
   if( SCIPconshdlrDecompGetNDecdecomps(scip) > 0 )
   {
      return SCIP_OKAY;
   }

   result = SCIP_DIDNOTFIND;
   nconss = SCIPgetNConss(scip);

   runs = conshdlrdata->setppcinmaster ? 1:2;
   SCIP_CALL( SCIPallocBufferArray(scip, &conshdlrdata->consismaster, nconss) );

   for (i = 0; i < runs && result != SCIP_SUCCESS; ++i)
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting %s structure:", conshdlrdata->setppcinmaster ? "extended":"blockdiagonal" );
      SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->clock) );
      SCIP_CALL( SCIPstartClock(scip, conshdlrdata->clock) );

      SCIP_CALL( findConnectedComponents(scip, conshdlrdata, &result) );

      SCIP_CALL( SCIPstopClock(scip, conshdlrdata->clock) );

      SCIPdebugMessage("Detection took %fs.\n", SCIPclockGetTime(conshdlrdata->clock));
      if(result == SCIP_SUCCESS)
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " found %d blocks.\n", conshdlrdata->nblocks);
         conshdlrdata->blockdiagonal = !conshdlrdata->setppcinmaster;
         SCIP_CALL( DECdecdecompCreate(scip, &conshdlrdata->decdecomp) );
         SCIP_CALL( copyToDecdecomp(scip, conshdlrdata, conshdlrdata->decdecomp) );
         SCIP_CALL( SCIPconshdlrDecompAddDecdecomp(scip, conshdlrdata->decdecomp) );
      }
      else
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " not found.\n");
      }
      if( conshdlrdata->setppcinmaster == FALSE && result != SCIP_SUCCESS )
      {
         conshdlrdata->setppcinmaster = TRUE;
      }
   }
   SCIPfreeBufferArray(scip, &conshdlrdata->consismaster);

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolConnected)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA *conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->clock != NULL )
      SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->clock) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpConnected)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsConnected)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckConnected)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockConnected)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for connected constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrConnected(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create connected constraint handler data */
   conshdlrdata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   assert(conshdlrdata != NULL);
   conshdlrdata->clock = NULL;
   conshdlrdata->constoblock = NULL;
   conshdlrdata->vartoblock = NULL;
   conshdlrdata->blockdiagonal = FALSE;

   conshdlrdata->nblocks = 0;
   conshdlrdata->enable = TRUE;
   conshdlrdata->decdecomp = NULL;
   conshdlrdata->consismaster = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyConnected,
         consFreeConnected, consInitConnected, consExitConnected,
         consInitpreConnected, consExitpreConnected, consInitsolConnected, consExitsolConnected,
         consDeleteConnected, consTransConnected, consInitlpConnected,
         consSepalpConnected, consSepasolConnected, consEnfolpConnected, consEnfopsConnected, consCheckConnected,
         consPropConnected, consPresolConnected, consRespropConnected, consLockConnected,
         consActiveConnected, consDeactiveConnected,
         consEnableConnected, consDisableConnected, consDelvarConnected,
         consPrintConnected, consCopyConnected, consParseConnected,
         conshdlrdata) );

   /* add connected constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/connected/enable", "Controls whether block diagonal detection is enabled", &conshdlrdata->enable, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/connected/setppcinmaster", "Controls whether SETPPC constraints chould be ignored while detecting", &conshdlrdata->setppcinmaster, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a connected constraint */
SCIP_RETCODE SCIPcreateConsConnected(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name                /**< name of constraint */
   )
{

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the connected constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("connected constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   /** @todo create and store constraint specific data here */

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         TRUE, TRUE, FALSE, TRUE, TRUE) );

   return SCIP_OKAY;
}


/** returns whether a block diagonal structure was found */
SCIP_Bool SCIPisMatrixBlockDiagonal(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA *conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   return conshdlrdata->blockdiagonal;
}
