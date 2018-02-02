/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_staircase.c
 * @ingroup DETECTORS
 * @brief  detector for staircase matrices
 * @author Martin Bergner
 *
 * This detector detects staircase structures in the constraint matrix by searching for the longest shortest path
 * in the row graph of the matrix.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "dec_staircase.h"
#include "cons_decomp.h"
#include "scip_misc.h"
#include "pub_decomp.h"
#include "tclique/tclique.h"

/* constraint handler properties */
#define DEC_DETECTORNAME         "staircase"    /**< name of detector */
#define DEC_DESC                 "Staircase detection via shortest paths" /**< description of detector */
#define DEC_PRIORITY             200            /**< priority of the detector */
#define DEC_DECCHAR              'S'            /**< display character of detector */
#define DEC_ENABLED              TRUE           /**< should the detection be enabled */
#define DEC_SKIP                 FALSE          /**< should detector be skipped if others found detections */


#define TCLIQUE_CALL(x) do                                                                                    \
                       {                                                                                      \
                          if((x) != TRUE )                                                      \
                          {                                                                                   \
                             SCIPerrorMessage("Error in function call\n");                                    \
                             return SCIP_ERROR;                                                               \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

/*
 * Data structures
 */

/** constraint handler data */
struct DEC_DetectorData
{
   SCIP_HASHMAP* constoblock;
   SCIP_HASHMAP* vartoblock;
   TCLIQUE_GRAPH* graph;
   int nblocks;
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static SCIP_DECL_SORTPTRCOMP(cmp)
{
   if( elem1 == elem2 )
      return 0;
   else if( elem1 < elem2 )
      return -1;
   else {
      assert(elem1 > elem2);
      return 1;
   }
}

/** creates the graph from the constraint matrix */
static
SCIP_RETCODE createGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH**       graph               /**< Graph data structure */
   )
{
   int i;
   int j;
   int v;
   int nconss;
   SCIP_CONS** conss;

   assert(scip != NULL);
   assert(graph != NULL);

   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);

   TCLIQUE_CALL( tcliqueCreate(graph) );
   assert(*graph != NULL);

   for( i = 0; i < nconss; ++i )
   {
      TCLIQUE_CALL( tcliqueAddNode(*graph, i, 0) );
   }

   /* Be aware: the following has n*n*m*log(m) complexity but doesn't need any additional memory
    * With additional memory, we can get it down to probably n*m + m*m*n
    */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_VAR** curvars1;
      int ncurvars1;

      ncurvars1 = GCGconsGetNVars(scip, conss[i]);
      if( ncurvars1 == 0 )
         continue;

      SCIP_CALL( SCIPallocBufferArray(scip, &curvars1, ncurvars1) );

      SCIP_CALL( GCGconsGetVars(scip, conss[i], curvars1, ncurvars1) );

      SCIPsortPtr((void**)curvars1, cmp, ncurvars1);
      for( j = i+1; j < nconss; ++j )
      {
         SCIP_VAR** curvars2;
         int ncurvars2;

         ncurvars2 = GCGconsGetNVars(scip, conss[j]);
         if( ncurvars2 == 0 )
            continue;

         SCIP_CALL( SCIPallocBufferArray(scip, &curvars2, ncurvars2) );

         SCIP_CALL( GCGconsGetVars(scip, conss[j], curvars2, ncurvars2) );

         SCIPsortPtr((void**)curvars1, cmp, ncurvars1);
         for( v = 0; v < ncurvars2; ++v )
         {
            int pos;

            if( SCIPsortedvecFindPtr((void*)curvars1, cmp, curvars2[v], ncurvars1, &pos) )
            {
               assert(curvars1[pos] == curvars2[v]);
               TCLIQUE_CALL( tcliqueAddEdge(*graph, i, j) );
               break;
            }
         }
         SCIPfreeBufferArray(scip, &curvars2);
      }

      SCIPfreeBufferArray(scip, &curvars1);
   }

   TCLIQUE_CALL( tcliqueFlush(*graph) );
   SCIPdebug(tcliquePrintGraph(*graph));
   return SCIP_OKAY;
}


/** returns the distance between vertex i and j based on the distance matrix */
static
int getDistance(
   int                   i,                  /**< vertex i */
   int                   j,                  /**< vertex j */
   int**                 distance            /**< triangular distance matrix */
   )
{
   assert(distance != NULL);

   if( i >= j )
      return distance[i][j];
   else if (i < j)
      return distance[j][i];
   else
      return 0;
}


/** perform BFS on the graph, storing distance information in the user supplied array */
static
SCIP_RETCODE doBFS(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< constraint handler data structure */
   int                   startnode,          /**< starting node */
   int**                 distances           /**< triangular matrix to store the distance when starting from node i */
   )
{
   int *queue;
   int nnodes;
   SCIP_Bool* marked;
   int squeue;
   int equeue;
   int i;
   int* node;

   TCLIQUE_GRAPH* graph;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(distances != NULL);
   assert(detectordata->graph != NULL);
   graph = detectordata->graph;
   nnodes = tcliqueGetNNodes(graph);

   assert(startnode < tcliqueGetNNodes(graph));

   squeue = 0;
   equeue = 0;

   SCIP_CALL( SCIPallocMemoryArray(scip, &queue, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &marked, nnodes) );

   for( i = 0; i < nnodes; ++i )
   {
      marked[i] = FALSE;
   }

   queue[equeue] = startnode;
   ++equeue;

   distances[startnode][startnode] = 0;
   marked[startnode] = TRUE;

   while( equeue > squeue )
   {
      int currentnode;
      int* lastneighbour;

      /* dequeue new node */
      currentnode = queue[squeue];

      assert(currentnode < nnodes);
      ++squeue;

      lastneighbour = tcliqueGetLastAdjedge(graph, currentnode);
      /* go through all neighbours */
      for( node = tcliqueGetFirstAdjedge(graph, currentnode); node <= lastneighbour; ++node )
      {
         if( !marked[*node] )
         {
            int curdistance;

            curdistance = getDistance(startnode, currentnode, distances);

            marked[*node] = TRUE;
            queue[equeue] = *node;
            if( *node < startnode )
               distances[startnode][*node] = curdistance+1;
            else if( *node > startnode )
               distances[*node][startnode] = curdistance+1;

            ++equeue;
         }
      }
   }

   SCIPfreeMemoryArray(scip, &queue);
   SCIPfreeMemoryArray(scip, &marked);

   return SCIP_OKAY;
}

/** finds the maximal shortest path by inspecting the distance array and returns the path in start and end*/
static
SCIP_RETCODE findMaximalPath(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< constraint handler data structure */
   int**                 distance,           /**< distance matrix of the maximal s-t path starting from s to any node t*/
   int*                  start,              /**< start vertex */
   int*                  end                 /**< end vertex */
   )
{
   int i;
   int j;
   int max;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(distance != NULL);
   assert(start != NULL);
   assert(end != NULL);

   max = -1;
   *start = -1;
   *end = -1;

   for( i = 0; i < tcliqueGetNNodes(detectordata->graph); ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( distance[i][j] > max )
         {
            max = distance[i][j];
            *start = i;
            *end = j;
         }
      }
   }
   assert(*start >= 0);
   assert(*end >= 0);

   SCIPdebugMessage("Path from %d to %d is longest %d.\n", *start, *end, max);
   detectordata->nblocks = max+1;

   return SCIP_OKAY;
}

/** this method will construct the cuts based on the longest shortest path and the distance matrix */
static
SCIP_RETCODE constructCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< constraint handler data structure */
   int                   start,              /**< start vertex */
   int                   end,                /**< end vertex */
   int**                 distance,           /**< distance matrix giving the distance from any constraint to any constraint */
   SCIP_VAR****          cuts                /**< which variables should be in the cuts */
   )
{

   int nnodes;
   SCIP_CONS** conss;
   int i;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(start >= 0);
   assert(end >= 0);
   assert(distance != NULL);
   assert(cuts != NULL);

   assert(detectordata->graph != NULL);
   nnodes = tcliqueGetNNodes(detectordata->graph);
   conss = SCIPgetConss(scip);
   assert(start < nnodes);
   assert(end < nnodes);

   /* The cuts will be generated on a trivial basis:
    * The vertices  of distance i will be in block i
    */

   for( i = 0; i < nnodes; ++i )
   {
      int dist;
      dist = getDistance(start, i, distance);
      SCIPdebugPrintf("from %d to %d = %d (%s = %d)\n", start, i, dist, SCIPconsGetName(conss[i]), dist+1 );
      SCIP_CALL( SCIPhashmapInsert(detectordata->constoblock, conss[i], (void*) (size_t) (dist+1)) );
   }

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
   int** distance;
   int i;
   int start;
   int end;
   SCIP_VAR*** cuts;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(result != NULL);

   nconss = SCIPgetNConss(scip);

   /* allocate triangular distance matrix */
   SCIP_CALL( SCIPallocMemoryArray(scip, &distance, nconss) );
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(distance[i]), (size_t)i+1) ); /*lint !e866*/
      BMSclearMemoryArray(distance[i], (size_t)i+1); /*lint !e866*/
   }

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( doBFS(scip, detectordata, i, distance) );
   }

   SCIP_CALL( findMaximalPath(scip, detectordata, distance, &start, &end) );
   SCIP_CALL( constructCuts(scip, detectordata, start, end, distance, &cuts) );

   for( i = 0; i < nconss; ++i )
   {
      SCIPfreeMemoryArray(scip, &distance[i]);
   }
   SCIPfreeMemoryArray(scip, &distance);

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

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(decdecomp != NULL);

   SCIP_CALL( DECfilloutDecompFromConstoblock(scip, decdecomp, detectordata->constoblock, detectordata->nblocks, TRUE) );

   return SCIP_OKAY;
}

/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(detectorFreeStaircase)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** detector initialization method (called after the problem has been transformed) */
static
DEC_DECL_INITDETECTOR(detectorInitStaircase)
{  /*lint --e{715}*/

   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   detectordata->graph = NULL;
   detectordata->constoblock = NULL;
   detectordata->vartoblock = NULL;
   detectordata->nblocks = 0;

   return SCIP_OKAY;
}

/** detector deinitialization method (called before the transformed problem is freed) */
static
DEC_DECL_EXITDETECTOR(detectorExitStaircase)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   if( detectordata->graph != NULL )
   {
      tcliqueFree(&detectordata->graph);
   }

   return SCIP_OKAY;
}

/** detector structure detection method, tries to detect a structure in the problem */
static
DEC_DECL_DETECTSTRUCTURE(detectorDetectStaircase)
{
   *result = SCIP_DIDNOTFIND;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting staircase structure:");

   SCIP_CALL( createGraph(scip, &(detectordata->graph)) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip), SCIPgetNConss(scip)) );

   SCIP_CALL( findStaircaseComponents(scip, detectordata, result) );

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
      if( detectordata->constoblock != NULL )
         SCIPhashmapFree(&detectordata->constoblock);
      if( detectordata->vartoblock != NULL )
         SCIPhashmapFree(&detectordata->vartoblock);
   }

   tcliqueFree(&detectordata->graph);

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for staircase constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorStaircase(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /* create staircase constraint handler data */
   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP,
      detectordata, detectorDetectStaircase, detectorFreeStaircase, detectorInitStaircase, detectorExitStaircase) );

   return SCIP_OKAY;
}
