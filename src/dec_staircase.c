/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_staircase.c
 * @ingroup DETECTORS
 * @brief  detector for staircase matrices
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */

#include <assert.h>
#include <string.h>

#include "dec_staircase.h"
#include "cons_decomp.h"
#include "scip_misc.h"
#include "pub_decomp.h"
#include "dijkstra/dijkstra.h"

/* constraint handler properties */
#define DEC_DETECTORNAME         "staircase"    /**< name of detector */
#define DEC_PRIORITY             0              /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              'S'            /**< display character of detector */
#define DEC_ENABLED              TRUE           /**< should the detection be enabled */

/*
 * Data structures
 */

/** constraint handler data */
struct DEC_DetectorData
{
   SCIP_HASHMAP* constoblock;
   SCIP_HASHMAP* vartoblock;
   Dijkstra_Graph graph;

   SCIP_CLOCK* clock;
   int nblocks;
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

/** looks for staircase components in the constraints in detectordata */
static
SCIP_RETCODE findStaircaseComponents(
   SCIP*              scip,         /**< SCIP data structure */
   DEC_DETECTORDATA* detectordata, /**< constraint handler data structure */
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
   assert(detectordata != NULL);
   assert(result != NULL);

   /* initialize data structures */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   nextblock = 1; /* start at 1 in order to see whether the hashmap has a key*/

   SCIP_CALL( SCIPallocBufferArray(scip, &vartoblock, nvars+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blockrepresentative, nconss+1) );
   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss+1) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip), nconss) );

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

   /* go through the all constraints */
   for( i = 0; i < nconss; ++i )
   {
      int consblock;

      cons = conss[i];
      assert(cons != NULL);
      if( isConsGCGCons(cons) )
         continue;

      /* get variables of constraint */
      ncurvars = SCIPgetNVarsXXX(scip, cons);
      curvars = NULL;
      if( ncurvars > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, cons, curvars, ncurvars) );
      }
      assert(ncurvars >= 0);
      assert(ncurvars <= nvars);
      assert(curvars != NULL || ncurvars == 0);

      assert(SCIPhashmapGetImage(constoblock, cons) == NULL);

      /* if there are no variables, put it in the first block, otherwise put it in the next block */
      if( ncurvars == 0 )
         consblock = -1;
      else
         consblock = nextblock;

      /* go through all variables */
      for( j = 0; j < ncurvars; ++j )
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
         /* get block of variable */
         varblock = vartoblock[varindex];
         SCIPdebugMessage("\tVar %s (%d): ", SCIPvarGetName(probvar), varblock);
         /* if variable is assigned to a block, assign constraint to that block */
         if( varblock > -1 && varblock != consblock )
         {
            consblock = MIN(consblock, varblock);
            SCIPdebugPrintf("still in block %d.\n",  varblock);
         }
         else if( varblock == -1 )
         {
            /* if variable is free, assign it to the new block for this constraint */
            varblock = consblock;
            assert(varblock > 0);
            assert(varblock <= nextblock);
            vartoblock[varindex] = varblock;
            SCIPdebugPrintf("new in block %d.\n",  varblock);
         }
         else
         {
            assert((varblock > 0) && (consblock == varblock));
            SCIPdebugPrintf("no change.\n");
         }
      }

      /* if the constraint belongs to a new block, mark it as such */
      if( consblock == nextblock )
      {
         assert(consblock > 0);
         blockrepresentative[consblock] = consblock;
         assert(blockrepresentative[consblock] > 0);
         assert(blockrepresentative[consblock] <= nextblock);
         ++nextblock;
      }

      SCIPdebugMessage("Cons %s will be in block %d (next %d)\n", SCIPconsGetName(cons), consblock, nextblock);
      for( k = 0; k < ncurvars; ++k )
      {
         int curvarindex;
         SCIP_VAR* curprobvar;
         int oldblock;
         assert(curvars != NULL);

         curprobvar = SCIPvarGetProbvar(curvars[k]);
         curvarindex = SCIPvarGetProbindex(curprobvar);
         oldblock = vartoblock[curvarindex];
         assert((oldblock > 0) && (oldblock <= nextblock));
         SCIPdebugMessage("\tVar %s ", SCIPvarGetName(curprobvar));
         if( oldblock != consblock )
         {
            SCIPdebugPrintf("reset from %d to block %d.\n", oldblock, consblock);
            vartoblock[curvarindex] = consblock;

            if( (blockrepresentative[oldblock] != -1) && (blockrepresentative[oldblock] > consblock))
            {
               int oldrepr;
               oldrepr = blockrepresentative[oldblock];
               SCIPdebugMessage("\t\tBlock representative from block %d changed from %d to %d.\n", oldblock, blockrepresentative[oldblock], consblock );
               assert(consblock > 0);
               blockrepresentative[oldblock] = consblock;
               if( (oldrepr != consblock) && (oldrepr != oldblock) )
               {
                  blockrepresentative[oldrepr] = consblock;
                  SCIPdebugMessage("\t\tBlock representative from block %d changed from %d to %d.\n", oldrepr, blockrepresentative[oldrepr], consblock );
               }
            }
         }
         else
         {
            SCIPdebugPrintf("will not be changed from %d to %d.\n", oldblock, consblock);
         }
      }

      SCIPfreeBufferArrayNull(scip, &curvars);
      assert(consblock >= 1 || consblock == -1);
      assert(consblock <= nextblock);

      /* store the constraint block */
      if(consblock != -1)
      {
         SCIPdebugMessage("cons %s in block %d\n", SCIPconsGetName(cons), consblock);
         SCIP_CALL( SCIPhashmapInsert(constoblock, cons, (void*)(size_t)consblock) );
      }
      else
      {
         SCIPdebugMessage("ignoring %s\n", SCIPconsGetName(cons));
      }
   }

   tempblock = 1;

   SCIPdebugPrintf("Blocks: ");
   /* postprocess blockrepresentatives */
   for( i = 1; i < nextblock; ++i )
   {
      /* forward replace the representatives */
      assert(blockrepresentative[i] >= 0);
      assert(blockrepresentative[i] < nextblock);
      if( blockrepresentative[i] != i )
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

   /* convert temporary data to detectordata */
   for( i = 0; i < nconss; ++i )
   {
      int consblock;

      cons = conss[i];
      if( isConsGCGCons(cons) )
         continue;

      if(!SCIPhashmapExists(constoblock, cons))
         continue;

      consblock = (int)(size_t) SCIPhashmapGetImage(constoblock, cons); /*lint !e507*/
      assert(consblock > 0);
      consblock = blockrepresentative[consblock];
      assert(consblock < tempblock);
      SCIP_CALL( SCIPhashmapInsert(detectordata->constoblock, cons, (void*)(size_t)consblock) );
      SCIPdebugMessage("%d %s\n", consblock, SCIPconsGetName(cons));
   }

   SCIP_CALL( SCIPhashmapCreate(&detectordata->vartoblock, SCIPblkmem(scip), nvars+1) );

   for( i = 0; i < nvars; ++i )
   {
      int varindex;
      int varblock;
      varindex = SCIPvarGetProbindex(SCIPvarGetProbvar(vars[i]));
      assert(varindex >= 0);
      assert(varindex < nvars);

      assert(vartoblock[varindex] < nextblock);
      if( vartoblock[varindex] < 0 )
         continue;

      varblock = blockrepresentative[vartoblock[varindex]];
      assert(varblock == -1 || varblock > 0);
      if( varblock > 0 )
      {
         assert(varblock < tempblock);
         SCIPdebugMessage("Var %s in block %d\n", SCIPvarGetName(SCIPvarGetProbvar(vars[i])), varblock-1);
         SCIP_CALL( SCIPhashmapInsert(detectordata->vartoblock, SCIPvarGetProbvar(vars[i]),
               (void*)(size_t)(varblock)) );
      }
   }

   /* free method data */
   SCIPfreeBufferArray(scip, &vartoblock);
   SCIPfreeBufferArray(scip, &blockrepresentative);
   SCIPhashmapFree(&constoblock);
   detectordata->nblocks = tempblock-1;

   if( detectordata->nblocks > 1 )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/* copy conshdldata data to decdecomp */
static
SCIP_RETCODE copyToDecdecomp(
   SCIP*              scip,         /**< SCIP data structure */
   DEC_DETECTORDATA* detectordata, /**< constraint handler data structure */
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
   assert(detectordata != NULL);
   assert(decdecomp != NULL);

   assert(DECdecdecompGetType(decdecomp) == DEC_DECTYPE_UNKNOWN);

   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);
   nlinkingconss = 0;
   nblocks = detectordata->nblocks;

   SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipvars, nblocks) ) ;
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linkingconss, nconss) );

   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars[i], nvars) ); /*lint !e866*/
      nsubscipvars[i] = 0;
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss[i], nconss) ); /*lint !e866*/
      nsubscipconss[i] = 0;
   }

   DECdecdecompSetNBlocks(decdecomp, nblocks);
   DECdecdecompSetConstoblock(decdecomp, detectordata->constoblock);
   DECdecdecompSetVartoblock(decdecomp, detectordata->vartoblock);

   for( i = 0; i < nconss; ++i )
   {
      size_t consblock;
      if( isConsGCGCons(conss[i]) )
         continue;

      consblock = (size_t) SCIPhashmapGetImage(detectordata->constoblock, conss[i]); /*lint !e507*/
      assert(consblock > 0);
      assert(nblocks >= 0);
      assert(consblock <= (size_t)nblocks);

      subscipconss[consblock-1][nsubscipconss[consblock-1]] = conss[i];
      ++(nsubscipconss[consblock-1]);
   }

   for( i = 0; i < nvars; ++i )
   {
      size_t varblock;
      SCIP_VAR* var;
      var = SCIPvarGetProbvar(vars[i]);
      if(var == NULL)
         continue;
      varblock = (size_t) SCIPhashmapGetImage(detectordata->vartoblock, SCIPvarGetProbvar(vars[i])); /*lint !e507*/

      assert(varblock > 0);
      assert(nblocks >= 0);
      assert(varblock <= (size_t) nblocks);

      subscipvars[varblock-1][nsubscipvars[varblock-1]] = SCIPvarGetProbvar(vars[i]);
      ++(nsubscipvars[varblock-1]);
   }

   if( nlinkingconss > 0 )
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


   for( i = 0; i < detectordata->nblocks; ++i )
   {
      SCIPfreeBufferArray(scip, &subscipvars[i]);
      SCIPfreeBufferArray(scip, &subscipconss[i]);
   }

   SCIPfreeBufferArray(scip, &subscipvars);
   SCIPfreeBufferArray(scip, &subscipconss);
   SCIPfreeBufferArray(scip, &nsubscipvars);
   SCIPfreeBufferArray(scip, &nsubscipconss);
   SCIPfreeBufferArray(scip, &linkingconss);

   detectordata->vartoblock = NULL;
   detectordata->constoblock = NULL;

   return SCIP_OKAY;
}

static
SCIP_RETCODE initDijkstraGraph(
   SCIP* scip,
   DEC_DETECTORDATA* detectordata
   )
{
   int i;
   int j;
   int nconss;
   int nvars;
   int nedges;
   SCIP_CONS** conss;

   Dijkstra_Graph* g = &(detectordata->graph);
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);
   nedges = 0;
   g->nodes = nconss;
   g->arcs = 2*nconss*nvars; /**@todo make better */

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->outbeg), (int) g->nodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->outcnt), (int) g->nodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->head), (int) g->arcs) );

   for( i = 0; i < nconss; ++i )
   {
      SCIP_VAR** curvars;
      int ncurvars;
      ncurvars = SCIPgetNVarsXXX(scip, conss[i]);
      g->outcnt[i] = ncurvars;
      SCIP_CALL( SCIPallocMemoryArray(scip, &curvars, ncurvars) );
      SCIP_CALL( SCIPgetVarsXXX(scip, conss[i], curvars, ncurvars) );
      nedges += ncurvars;
      for( j = 0; j < ncurvars; ++j )
      {
         int pindex = SCIPvarGetProbindex(SCIPvarGetProbvar(curvars[j]));
         vartocons[pindex][nvartoconss[pindex]] = i;
      }

      SCIPfreeMemoryArray(scip, &curvars);
   }

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
DEC_DECL_INITDETECTOR(initStaircase)
{  /*lint --e{715}*/

   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   detectordata->clock = NULL;
   detectordata->constoblock = NULL;
   detectordata->vartoblock = NULL;

   detectordata->nblocks = 0;
   SCIP_CALL( initDijkstraGraph(scip, detectordata) );
   SCIP_CALL( SCIPcreateClock(scip, &detectordata->clock) );

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
DEC_DECL_EXITDETECTOR(exitStaircase)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   if( detectordata->clock != NULL )
      SCIP_CALL( SCIPfreeClock(scip, &detectordata->clock) );

   SCIPfreeMemory(scip, &detectordata);
   return SCIP_OKAY;

}

static
DEC_DECL_DETECTSTRUCTURE(detectStaircase)
{
   *result = SCIP_DIDNOTFIND;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting staircase structure:");

   SCIP_CALL( SCIPstartClock(scip, detectordata->clock) );

   SCIP_CALL( findStaircaseComponents(scip, detectordata, result) );

   SCIP_CALL( SCIPstopClock(scip, detectordata->clock) );

   SCIPdebugMessage("Detection took %fs.\n", SCIPgetClockTime(scip, detectordata->clock));
   if( *result == SCIP_SUCCESS )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " found %d blocks.\n", detectordata->nblocks);
      SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, 1) ); /*lint !e506*/
      SCIP_CALL( DECdecdecompCreate(scip, &((*decdecomps)[0])) );
      SCIP_CALL( copyToDecdecomp(scip, detectordata, (*decdecomps)[0]) );
      *ndecdecomps = 1;
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " not found.\n");
      SCIPhashmapFree(&detectordata->constoblock);
      SCIPhashmapFree(&detectordata->vartoblock);
   }

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for staircase constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionStaircase(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /* create staircase constraint handler data */
   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_PRIORITY, DEC_ENABLED, detectordata, detectStaircase, initStaircase, exitStaircase) );

   return SCIP_OKAY;
}
