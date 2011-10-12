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
#define SCIP_DEBUG

#include <assert.h>
#include <string.h>

#include "cons_connected.h"
#include "cons_decomp.h"
#include "scip_misc.h"
#include "struct_decomp.h"
#include "scip/clock.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "connected"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
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

/* TODO: fill in the necessary constraint data */

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
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static
SCIP_Bool isConsGCGCons(
   SCIP_CONS* cons
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

static
SCIP_RETCODE findConnectedComponents(
   SCIP* scip,
   SCIP_CONSHDLRDATA* conshdlrdata,
   SCIP_RESULT* result
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

   SCIP_CALL( SCIPallocMemoryArray(scip, &vartoblock, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &blockrepresentative, nconss+1) );
   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->constoblock, SCIPblkmem(scip), nconss) );

   for( i = 0; i < nvars; ++i )
   {
      vartoblock[i] = -1;
   }

   for( i = 0; i < nconss+1; ++i )
   {
      blockrepresentative[i] = -1;
   }

   assert(nconss >= 1);

   /* process the first constraint */
   cons = conss[0];
   ncurvars = SCIPgetNVarsXXX(scip, cons);
   curvars = SCIPgetVarsXXX(scip, cons);
   assert(ncurvars >= 0);
   assert(ncurvars <= nvars);
   assert(curvars != NULL || ncurvars == 0);

   SCIP_CALL( SCIPhashmapInsert(constoblock, cons, (void*)(size_t)nextblock) );
   for( j = 0; j < ncurvars; ++j)
   {
      SCIP_VAR* probvar;
      int varindex;
      assert(curvars != NULL);

      probvar = SCIPvarGetProbvar(curvars[j]);
      assert(probvar != NULL);

      varindex = SCIPvarGetProbindex(probvar);
      assert(varindex >= 0);
      assert(varindex < nvars);
      vartoblock[varindex] = nextblock;

   }
   blockrepresentative[0] = nextblock;

   /* prepare consblock for the next costraint */
   ++nextblock;

   SCIPfreeMemoryArrayNull(scip, &curvars);
   /* go through the remaining constraints */
   for( i = 1; i < nconss; ++i )
   {
      int consblock;
      cons = conss[i];
      assert(cons != NULL);
      if( isConsGCGCons(cons) )
         continue;

      ncurvars = SCIPgetNVarsXXX(scip, cons);
      curvars = SCIPgetVarsXXX(scip, cons);
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

         /** @todo: what about deleted variables? */
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
            }
            /* assign all previous variables of this constraint to this block */
            for (k = j; k >= 0; --k)
            {
               /** @todo: what about deleted variables? */
               vartoblock[SCIPvarGetProbindex(SCIPvarGetProbvar(curvars[k]))] = consblock;
            }
         }
         else
         {
            /* if variable is free, assign it to the new block for this constraint */
            varblock = consblock;
            vartoblock[varindex] = varblock;
         }
         if( consblock == nextblock )
         {
            blockrepresentative[consblock] = consblock;
            ++nextblock;
         }
      }

      SCIPfreeMemoryArrayNull(scip, &curvars);
      assert(consblock >= 1);
      SCIP_CALL( SCIPhashmapInsert(constoblock, cons, (void*)(size_t)consblock) );

   }

   tempblock = 1;
   /* postprocess blockrepresentatives */
   for( i = 1; i < nextblock; ++i )
   {
      /* forward replace the representatives */
      if(blockrepresentative[i] != i)
         blockrepresentative[i] = blockrepresentative[blockrepresentative[i]];
      else
      {
         blockrepresentative[i] = tempblock;
         ++tempblock;
      }
      /* It is crucial that this condition holds */
      assert(blockrepresentative[i] <= i);
   }



   /* convert temporary data to conshdlrdata */
   for(i = 0; i < nconss; ++i)
   {
      int consblock;

      cons = conss[i];
      if( isConsGCGCons(cons) )
         continue;

      consblock = (size_t)SCIPhashmapGetImage(constoblock, cons);
      assert(consblock > 0);
      consblock = blockrepresentative[consblock];
      assert(consblock < nextblock);
      SCIP_CALL( SCIPhashmapInsert(conshdlrdata->constoblock, cons, (void*)(size_t)consblock) );
      SCIPdebugMessage("%d %s\n", consblock, SCIPconsGetName(cons));
   }

   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->vartoblock, SCIPblkmem(scip), nvars) );

   for( i = 0; i < nvars; ++i)
   {
      int varindex;
      int varblock;
      varindex = SCIPvarGetProbindex(SCIPvarGetProbvar(vars[i]));
      assert(varindex >= 0);
      assert(varindex < nvars);

      varblock = vartoblock[varindex];
      assert(varblock == -1 || varblock > 0);
      if(varblock > 0)
      {
         assert(varblock > 0);
         SCIP_CALL( SCIPhashmapInsert(conshdlrdata->vartoblock, SCIPvarGetProbvar(vars[i]), (void*)(size_t)(varblock)) );
      }

   }

   /* free method data */
   SCIPfreeMemoryArray(scip, &vartoblock);
   SCIPfreeMemoryArray(scip, &blockrepresentative);
   conshdlrdata->constoblock = constoblock;
   conshdlrdata->nblocks = tempblock;

   if(conshdlrdata->nblocks > 1)
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

static
SCIP_RETCODE copyToDecdecomp(
   SCIP* scip,
   SCIP_CONSHDLRDATA* conshdlrdata,
   DECDECOMP* decdecomp
   )
{
   SCIP_CONS** conss;
   int nconss;
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(decdecomp != NULL);

   if( !conshdlrdata->blockdiagonal )
      return SCIP_OKAY;

   assert(decdecomp->type == DEC_UNKNOWN);

   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   SCIPallocMemoryArray(scip, &decdecomp->subscipvars, conshdlrdata->nblocks);
   SCIPallocMemoryArray(scip, &decdecomp->nsubscipvars, conshdlrdata->nblocks);
   SCIPallocMemoryArray(scip, &decdecomp->subscipconss, conshdlrdata->nblocks);
   SCIPallocMemoryArray(scip, &decdecomp->nsubscipconss, conshdlrdata->nblocks);

   for( i = 0; i < conshdlrdata->nblocks; ++i)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->subscipvars[i], nvars) );
      decdecomp->nsubscipvars[i] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->subscipconss[i], nconss) );
      decdecomp->nsubscipconss[i] = 0;
   }

   decdecomp->nlinkingcuts = 0;
   decdecomp->nlinkingconss = 0;
   decdecomp->nblocks = conshdlrdata->nblocks;
   decdecomp->type = DEC_DIAGONAL;

   assert(decdecomp->constoblock == NULL);
   assert(decdecomp->vartoblock == NULL);

   decdecomp->constoblock = conshdlrdata->constoblock;
   decdecomp->vartoblock = conshdlrdata->vartoblock;

   for( i = 0; i < nconss; ++i)
   {
      size_t consblock;
      if( isConsGCGCons(conss[i]) )
         continue;

      consblock = (size_t) SCIPhashmapGetImage(decdecomp->constoblock, conss[i]);
      assert(consblock > 0);
      assert(decdecomp->nblocks >= 0);
      assert(consblock <= (size_t)decdecomp->nblocks);

      decdecomp->subscipconss[consblock-1][decdecomp->nsubscipconss[consblock-1]] = conss[i];
      ++(decdecomp->nsubscipconss[consblock-1]);
   }

   for( i = 0; i < nvars; ++i)
   {
      size_t varblock;
      varblock = (size_t) SCIPhashmapGetImage(decdecomp->vartoblock, vars[i]);
      assert(varblock > 0);
      assert(decdecomp->nblocks >= 0);
      assert(varblock <= (size_t) decdecomp->nblocks);

      decdecomp->subscipvars[varblock-1][decdecomp->nsubscipvars[varblock-1]] = vars[i];
      ++(decdecomp->nsubscipvars[varblock-1]);
   }

   for( i = 0; i < conshdlrdata->nblocks; ++i)
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &decdecomp->subscipvars[i], decdecomp->nsubscipvars[i]) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &decdecomp->subscipconss[i], decdecomp->nsubscipconss[i]) );
   }

   conshdlrdata->vartoblock = NULL;
   conshdlrdata->constoblock = NULL;

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyConnected NULL
#endif

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


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitConnected NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitConnected NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreConnected NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreConnected NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolConnected)
{  /*lint --e{715}*/

   SCIP_CONSHDLRDATA *conshdlrdata;
   SCIP_RESULT result;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( !conshdlrdata->enable )
      return SCIP_OKAY;

   if( conshdlrdata->decdecomp == NULL )
   {
      conshdlrdata->decdecomp = SCIPconshdlrDecompGetDecdecomp(scip);
   }
   /* apparently, there is a structure, which means we don't try to detect one */
   if( conshdlrdata->decdecomp->type != DEC_UNKNOWN )
   {
      return SCIP_OKAY;
   }

   assert(conshdlrdata->decdecomp != NULL);

   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->clock) );
   SCIP_CALL( SCIPstartClock(scip, conshdlrdata->clock) );

   SCIP_CALL( findConnectedComponents(scip, conshdlrdata, &result) );

   SCIP_CALL( SCIPstopClock(scip, conshdlrdata->clock) );

   SCIPdebugMessage("Detection took %f s.\n", SCIPclockGetTime(conshdlrdata->clock));

   if(result == SCIP_SUCCESS)
   {
      SCIPdebugMessage("Found block diagonal structure with %d blocks.\n", conshdlrdata->nblocks);
      conshdlrdata->blockdiagonal = TRUE;
   }
   else
   {
      SCIPdebugMessage("No block diagonal structure found.\n");
   }

   SCIP_CALL( copyToDecdecomp(scip, conshdlrdata, conshdlrdata->decdecomp) );
   SCIP_CALL( DECOMPconvertStructToGCG(scip, conshdlrdata->decdecomp) );
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


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteConnected NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransConnected NULL
#endif


/** LP initialization method of constraint handler */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpConnected NULL
#endif


/** separation method of constraint handler for LP solutions */
#define consSepalpConnected NULL

/** separation method of constraint handler for arbitrary primal solutions */
#define consSepasolConnected NULL


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


/** domain propagation method of constraint handler */
#define consPropConnected NULL

/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolConnected NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropConnected NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockConnected)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveConnected NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveConnected NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableConnected NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableConnected NULL
#endif

#define consDelvarConnected NULL

/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintConnected NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyConnected NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseConnected)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of connected constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseConnected NULL
#endif




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

   conshdlrdata->constoblock = NULL;
   conshdlrdata->vartoblock = NULL;
   conshdlrdata->blockdiagonal = FALSE;

   conshdlrdata->nblocks = 0;
   conshdlrdata->enable = TRUE;
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
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/connected/enable", "Controls whether block diagonal detection is enabled", &conshdlrdata->enable, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a connected constraint */
SCIP_RETCODE SCIPcreateConsConnected(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name                /**< name of constraint */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsConnected() call, if you don't need all the information */

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
   /* TODO: create and store constraint specific data here */

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         TRUE, TRUE, FALSE, TRUE, TRUE) );

   return SCIP_OKAY;
}

void SCIPconsConnectedSetDecomp(
   SCIP* scip,
   DECDECOMP* decdecomp
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA *conshdlrdata;

   assert(scip != NULL);
   assert(decdecomp != NULL);

   conshdlr = SCIPfindConshdlr(scip, "CONSHDLR_NAME");
   assert(conshdlr != NULL);

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   conshdlrdata->decdecomp = decdecomp;
}

/** returns whether a block diagonal structure was found */
extern
SCIP_Bool SCIPisMatrixBlockDiagonal(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA *conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, "CONSHDLR_NAME");
   assert(conshdlr != NULL);

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   return conshdlrdata->blockdiagonal;
}
