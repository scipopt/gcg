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
#define DEC_DESC                 "Staircase detection via shortest paths" /**< description of detector */
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
   DIJKSTRA_GRAPH graph;

   int** varconss;
   int* nvarconss;

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

/** perform BFS on the graph, storing distance information in the user supplied array */
static
SCIP_RETCODE doBFS(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< constraint handler data structure */
   unsigned int          startnode,          /**< starting node */
   int**                 distances           /**< triangular matrix to store the distance when starting from node i */
   )
{
   unsigned int *queue;
   SCIP_Bool* marked;
   int squeue;
   int equeue;
   unsigned int i;
   unsigned int j;

   DIJKSTRA_GRAPH* graph;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(distances != NULL);

   graph = &detectordata->graph;
   assert(i > 0 && i < graph->nodes);

   squeue = 0;
   equeue = 0;

   SCIP_CALL( SCIPallocMemoryArray(scip, &queue, graph->nodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &marked, graph->nodes) );

   for( i = 0; i < graph->nodes; ++i)
   {
      marked[i] = FALSE;
   }

   queue[equeue] = startnode;
   ++equeue;

   distances[startnode][startnode] = 0;
   marked[startnode] = TRUE;

   while(equeue > squeue)
   {
      unsigned int currentnode;

      /* dequeue new node */
      currentnode = queue[squeue];
      SCIPdebugMessage("Dequeueing %ud\n", currentnode);

      assert(currentnode < graph->nodes);
      ++squeue;

      /* go through all neighbours */
      for( j = 0; j < graph->outcnt[currentnode]; ++j )
      {
         int eindex;
         unsigned int targetnode;

         eindex = graph->outbeg[currentnode+j];
         targetnode = graph->head[eindex];

         if( !marked[targetnode] )
         {
            int curdistance;

            /* little magic for triangular distance matrix */
            if( startnode < currentnode )
            {
               curdistance = distances[startnode][currentnode];
            }
            else
            {
               curdistance = distances[currentnode][startnode];
            }

            marked[targetnode] = TRUE;
            queue[equeue] = targetnode;
            if(targetnode > startnode)
               distances[startnode][targetnode] = curdistance+1;
            else if(startnode > targetnode)
               distances[targetnode][startnode] = curdistance+1;

            ++equeue;
         }
      }
   }

   SCIPfreeMemoryArray(scip, &queue);
   SCIPfreeMemoryArray(scip, &marked);
   return SCIP_OKAY;
}

/** finds the maximal shortest path by inspecting the distance array and returns the path in start and end*/
SCIP_RETCODE findMaximalPath(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< constraint handler data structure */
   int**                 distance,           /**< distance matrix of the maximal s-t path starting from s to any node t*/
   unsigned int*         start,              /**< start vertex */
   unsigned int*         end                 /**< end vertex */
   )
{
   unsigned int i;
   unsigned int j;
   int max;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(distance != NULL);
   assert(start != NULL);
   assert(end != NULL);

   max = 0;

   for(i = 0; i < detectordata->graph.nodes; ++i)
   {
      for(j = 0; j < i; ++j)
      {
         if(distance[i][j] > max)
         {
            max = distance[i][j];
            *start = i;
            *end = j;
         }
      }
   }

   return SCIP_OKAY;
}

/** This method will construct the cuts based on the longest shortest path and the distance matrix */
SCIP_RETCODE constructCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< constraint handler data structure */
   unsigned int          start,              /**< start vertex */
   unsigned int          end,                /**< end vertex */
   int**                 distance,           /**< distance matrix giving the distance from any constraint to any constraint */
   SCIP_VAR****          cuts                /**< which variables should be in the cuts */
   )
{
   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(detectordata->graph != NULL);
   assert(start < detectordata->graph.nodes);
   assert(end < detectordata->graph.nodes);
   assert(distance != NULL);
   assert(cuts != NULL);

   return SCIP_OKAY;
}


/** converts the cuts to a structure that GCG can understand */
SCIP_RETCODE convertCutsToDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< constraint handler data structure */
   SCIP_VAR***           cuts                /**< which variables should be in the cuts */
   )
{
   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(cuts != NULL);

   return SCIP_OKAY;
}



/** looks for staircase components in the constraints in detectordata */
static
SCIP_RETCODE findStaircaseComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< constraint handler data structure */
   SCIP_RESULT*          result              /**< result pointer to indicate success oder failuer */
   )
{
   int nconss;
   int nvars;
   int** distance;
   int i;
   unsigned int start;
   unsigned int end;
   SCIP_VAR*** cuts;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(result != NULL);

   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);

   /* allocate triangular distance matrix */
   SCIP_CALL( SCIPallocMemoryArray(scip, &distance, nconss) );
   for( i = 0; i < nconss; ++i)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &distance[i], nconss-i) );
   }

   for( i = 0; i < nconss; ++i)
   {
      SCIP_CALL( doBFS(scip, detectordata, i, distance) );
   }

   SCIP_CALL( findMaximalPath(scip, detectordata, distance, &start, &end) );

   SCIP_CALL( constructCuts(scip, detectordata, start, end, distance, &cuts) );

   SCIP_CALL( convertCutsToDecomp(scip, detectordata, cuts) );

   if( detectordata->nblocks > 1 )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/* copy conshdldata data to decdecomp */
static
SCIP_RETCODE copyToDecdecomp(
   SCIP*              scip,                  /**< SCIP data structure */
   DEC_DETECTORDATA*  detectordata,          /**< constraint handler data structure */
   DEC_DECOMP*        decdecomp              /**< decdecomp data structure */
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

   SCIP_Bool valid;

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

   DECdecompSetNBlocks(decdecomp, nblocks);
   DECdecompSetConstoblock(decdecomp, detectordata->constoblock, &valid);
   assert(valid);
   DECdecompSetVartoblock(decdecomp, detectordata->vartoblock, &valid);
   assert(valid);

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
      SCIP_CALL( DECdecompSetLinkingconss(scip, decdecomp, linkingconss, nlinkingconss, &valid) );
      assert(valid);
      DECdecompSetType(decdecomp, DEC_DECTYPE_BORDERED, &valid);
      assert(valid);

   }
   else
   {
      DECdecompSetType(decdecomp, DEC_DECTYPE_DIAGONAL, &valid);
      assert(valid);
   }

   SCIP_CALL( DECdecompSetSubscipconss(scip, decdecomp, subscipconss, nsubscipconss, &valid) );
   assert(valid);

   SCIP_CALL( DECdecompSetSubscipvars(scip, decdecomp, subscipvars, nsubscipvars, &valid) );
   assert(valid);


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

   DIJKSTRA_GRAPH* g = &(detectordata->graph);
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
         detectordata->varconss[pindex][detectordata->nvarconss[pindex]] = i;
         ++(detectordata->nvarconss[pindex]);
      }
      SCIPfreeMemoryArray(scip, &curvars);
   }

   /* reallocate to necessary size */
   for( i = 0; i < nvars; ++i )
   {
      if( detectordata->nvarconss[i] > 0)
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, &detectordata->varconss[i], detectordata->nvarconss[i]) );
      }
      else
      {
         SCIPfreeMemoryArray(scip, &detectordata->varconss[i] );
         assert(detectordata->nvarconss[i] == 0);
      }

   }
   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
DEC_DECL_INITDETECTOR(initStaircase)
{  /*lint --e{715}*/

   DEC_DETECTORDATA *detectordata;
   int i;
   int nvars;
   int nconss;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);

   detectordata->clock = NULL;
   detectordata->constoblock = NULL;
   detectordata->vartoblock = NULL;

   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->nvarconss, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varconss, nvars) );

   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varconss[i], nconss) );
      detectordata->nvarconss[i] = 0;
   }

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
   int i;
   int nvars;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   nvars = SCIPgetNVars(scip);

   if( detectordata->clock != NULL )
      SCIP_CALL( SCIPfreeClock(scip, &detectordata->clock) );

   for( i = 0; i < nvars; ++i )
   {
      SCIPfreeMemoryArray(scip, &detectordata->varconss[i]);
   }

   SCIPfreeMemoryArray(scip, &detectordata->varconss);
   SCIPfreeMemoryArray(scip, &detectordata->nvarconss);

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
      SCIP_CALL( DECdecompCreate(scip, &((*decdecomps)[0])) );
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

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, detectordata, detectStaircase, initStaircase, exitStaircase) );

   return SCIP_OKAY;
}
