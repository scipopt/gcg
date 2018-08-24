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

/**@file   dec_cutpacking.c
 * @ingroup DETECTORS
 * @brief  staircase detector via recursive partitioning (uses hmetis)
 * @author Friederike Menge
 * @author Martin Bergner
 * @author Christian Puchert
 *
 * This detector tries to detect staircase structures by recursively partitioning the
 * rowgraph of the matrix by using hmetis.
 *
 * This detector needs hmetis and works only under Linux/MacOS. It further needs the Z-shell (zsh)
 * to enforce memory and time limits on hmetis as this is the only shell reliably doing that.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "dec_cutpacking.h"

#if !defined(_WIN32) && !defined(_WIN64)

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>


#include "cons_decomp.h"
#include "struct_decomp.h"
#include "pub_decomp.h"
#include "scip_misc.h"

#ifdef HMETIS_HEADER
#include "hmetis.h"
#else
#define HMETIS_EXECUTABLE "hmetis"
#endif

#define DEC_DETECTORNAME         "cutpacking"   /**< name of the detector */
#define DEC_DESC                 "detects staircase matrices via graph partioning and cutpacking" /**< detector description */
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              1100           /**< priority of the detector */
#define DEC_DECCHAR               'c'            /**< display character of detector */
#define DEC_ENABLED               FALSE          /**< should detector be called by default */
#define DEC_ENABLEDORIGINAL       FALSE  /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING      FALSE          /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE          /**< should detector be skipped if others found detections */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */
#define DEC_LEGACYMODE            FALSE       /**< should (old) DETECTSTRUCTURE method also be used for detection */

/* Default parameter settings */
#define DEFAULT_RANDSEED              1      /**< random seed for the hmetis call */
#define DEFAULT_TIDY               TRUE      /**< whether to clean up afterwards */
#define DEFAULT_FIXEDBLOCKS       FALSE      /**< whether the blocks should consist of a given number of constraints */
#define DEFAULT_BLOCKSIZE           200      /**< number of constraints per block */
#define DEFAULT_USEMETIS           TRUE      /**< should metis be used (TRUE) or the Stoer-Wagner algorithm */

#define DEFAULT_METIS_UBFACTOR      5.0      /**< default unbalance factor given to metis on the commandline */
#define DEFAULT_METIS_VERBOSE     FALSE      /**< should metis be verbose */
#define DEFAULT_METISUSEPTYPE_RB   TRUE      /**< Should metis use the rb or kway partitioning algorithm */

/*
 * Data structures
 */

/** adjacency list of a single vertex (representing a constraint) in a graph */
struct AdjList
{
   SCIP_CONS**           conss;
   int*                  weights;
   int                   nconss;
   int                   maxconss;
};
typedef struct AdjList ADJLIST;

/** graph structure, where
 *  - each vertex represents a constraint;
 *  - an edge between two vertices exists iff the constraints corresponding to the vertices
 *    have at least one variable in common;
 *  - the weight of an edge is equal to the number of these variables
 */
struct Graph
{
   ADJLIST**             adjlists;           /**< adjacencylists of the graph */
   SCIP_CONS**           conss;              /**< constraints (each constraint represents a vertex of the graph)*/
   int                   nconss;             /**< number of vertices */
   SCIP_HASHMAP*         constopos;          /**< assigns constraints to their position in conss */

   int                   nedges;             /**< number of edges */

   SCIP_CONS*            cons1;              /**< first linking constraint */
   SCIP_CONS*            cons2;              /**< second linking constraint */
};
typedef struct Graph GRAPH;

/** detector data */
struct DEC_DetectorData
{
   /* general parameters */
   SCIP_Bool             tidy;
   SCIP_Bool             usemetis;
   int                   blocksize;
   SCIP_Bool             fixedblocks;

   /* stuff for the algorithm */
   GRAPH**               graphs;
   int                   ngraphs;
   int                   maxgraphs;
   SCIP_CONS***          subscipconss;
   int*                  nsubscipconss;
   SCIP_HASHMAP*         constoblock;

   int                   nblocks;
   int                   position;
   int                   startblock;
   int*                  partition;

   SCIP_CONS***          mergedconss;
   int*                  nmergedconss;
   SCIP_CONS**           representatives;
   int                   nrepresentatives;

   /* general stuff */
   SCIP_HASHMAP*         vartopos;
   int*                  nvarinconss;
   SCIP_CONS***          varinconss;
   SCIP_VAR**            relvars;
   int                   nrelvars;
   int                   nrelconss;

   /* graph stuff for hmetis */
   int                   randomseed;
   SCIP_Real             metisubfactor;
   SCIP_Bool             metisverbose;
   SCIP_Bool             metisuseptyperb;

};


/*
 * Local methods
 */


/** initialize data and allocate memory, called just before detection */
static
SCIP_RETCODE initData(
   SCIP*                 scip,
   DEC_DETECTORDATA*     detectordata
   )
{
   int i;
   int j;
   int k;
   int nallvars;
   int nconss;
   SCIP_Bool ishandled;
   SCIP_CONS** conss;
   SCIP_CONS** newconss;
   SCIP_VAR** curvars;
   SCIP_VAR** allvars;
   SCIP_VAR** relvars;
   int ncurvars;
   SCIP_HASHMAP* vartopos;
   SCIP_CONS*** varinconss;
   int* nvarinconss;

   assert(scip != NULL);
   assert(detectordata != NULL);

   nallvars = SCIPgetNVars(scip);
   allvars = SCIPgetVars(scip);
   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);

   detectordata->nblocks = 0;
   detectordata->ngraphs = 0;
   detectordata->position = -1;
   detectordata->nrepresentatives = 0;
   detectordata->nrelvars = 0;
   detectordata->startblock = -1;

   /* get number of relevant variables */
   /* vartopos */
   SCIP_CALL( SCIPhashmapCreate(&detectordata->vartopos, SCIPblkmem(scip),nallvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->relvars, nallvars) ); /*lint !e522*/
   vartopos = detectordata->vartopos;
   relvars = detectordata->relvars;
   j = 0;

   for( i = 0; i < nallvars; ++i )
   {
      if( GCGisVarRelevant(allvars[i]) )
      {
         relvars[j] = SCIPvarGetProbvar(allvars[i]);
         SCIP_CALL( SCIPhashmapInsert(vartopos, SCIPvarGetProbvar(allvars[i]), (void*) (size_t) j) );
         j++;
      }
   }
   detectordata->nrelvars = j;
   SCIP_CALL( SCIPreallocMemoryArray(scip, &(detectordata->relvars), j) ); /*lint !e522*/

   /* get number of relevant conss */
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->graphs, (size_t)nconss+1) );
   SCIP_CALL( SCIPallocMemory(scip, &detectordata->graphs[0]) );
   for( i = 1; i < nconss + 1; ++i )
      detectordata->graphs[i] = NULL;
   detectordata->maxgraphs = nconss + 1;
   SCIP_CALL( SCIPallocMemoryArray(scip, &newconss, nconss) ); /*lint !e522*/
   k = 0;
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);
      if( !SCIPconsIsActive(conss[i]) )
      {
         continue;
      }

      ncurvars = GCGconsGetNVars(scip, conss[i]);
      if( ncurvars > 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &curvars, ncurvars) );
         SCIP_CALL( GCGconsGetVars(scip, conss[i], curvars, ncurvars) );

         ishandled = FALSE;

         for( j = 0; (j < ncurvars) && (ishandled == FALSE); ++j )
            ishandled = GCGisVarRelevant(curvars[j]);

         if( ishandled )
         {
            newconss[k] = conss[i];
            k++;
         }

         SCIPfreeMemoryArrayNull(scip, &curvars);
      }
   }

   SCIP_CALL( SCIPreallocMemoryArray(scip, &newconss, k) ); /*lint !e522*/
   detectordata->nrelconss = k;
   detectordata->graphs[0]->nconss = k;
   detectordata->graphs[0]->conss = newconss;

   /* alloc */
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->partition, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->representatives, k) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip),k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->subscipconss, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->mergedconss, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->nmergedconss, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->nsubscipconss, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->nvarinconss, detectordata->nrelvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varinconss, detectordata->nrelvars) );

   varinconss = detectordata->varinconss;
   nvarinconss = detectordata->nvarinconss;

   for( i = 0; i < k; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->subscipconss[i], k) ); /*lint !e866*/
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->mergedconss[i], k) ); /*lint !e866*/
   }

   /* varinconss */
   for( i = 0; i < detectordata->nrelvars; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varinconss[i], k) ); /*lint !e866*/
      detectordata->nvarinconss[i] = 0;
   }

   for( i = 0; i < k; ++i )
   {
      curvars = NULL;
      ncurvars = GCGconsGetNVars(scip, detectordata->graphs[0]->conss[i]);
      if( ncurvars > 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &curvars, ncurvars) );
         SCIP_CALL( GCGconsGetVars(scip,detectordata->graphs[0]->conss[i], curvars, ncurvars) );

         for( j = 0; j < ncurvars; ++j )
         {
            if( GCGisVarRelevant(curvars[j]) && GCGisVarRelevant(curvars[j]) )
            {
               int varpos;

               varpos = (int) (size_t) SCIPhashmapGetImage(vartopos, SCIPvarGetProbvar(curvars[j])); /*lint !e507*/

               (varinconss[varpos])[nvarinconss[varpos]] = detectordata->graphs[0]->conss[i];
               ++nvarinconss[varpos];
            }
         }

         SCIPfreeMemoryArrayNull(scip, &curvars);
      }
   }

   return SCIP_OKAY;
}

/** free data after detection */
static
SCIP_RETCODE freeData(
   SCIP*                 scip,
   DEC_DETECTORDATA*     detectordata
   )
{
   int i;

   assert(scip != NULL);
   assert(detectordata != NULL);

   /* free detector data */
   for( i = 0; i < detectordata->nrelconss; ++i )
   {
      SCIPfreeMemoryArray(scip, &detectordata->subscipconss[i]);
   }
   for( i = 0; i < detectordata->nrelvars; ++i )
   {
      SCIPfreeMemoryArray(scip, &detectordata->varinconss[i]);
   }

   SCIPfreeMemoryArray(scip, &detectordata->nsubscipconss);
   SCIPfreeMemoryArray(scip, &detectordata->subscipconss);
   SCIPfreeMemoryArray(scip, &detectordata->representatives);
   SCIPfreeMemoryArray(scip, &detectordata->partition);
   for( i = 0; i < detectordata->maxgraphs; ++i)   /* @todo: This should not be necessary anymore */
   {
      SCIPfreeMemoryNull(scip, &detectordata->graphs[i]);
   }
   SCIPfreeMemoryArray(scip, &detectordata->graphs);
   SCIPfreeMemoryArray(scip, &detectordata->varinconss);
   SCIPfreeMemoryArray(scip, &detectordata->nvarinconss);
   SCIPfreeMemoryArray(scip, &detectordata->relvars);
   SCIPhashmapFree(&detectordata->vartopos);

   SCIPfreeMemoryArray(scip, &detectordata->nmergedconss);
   for( i = 0; i < detectordata->nrelconss; i++ )
   {
      SCIPfreeMemoryArray(scip, &detectordata->mergedconss[i]);
   }
   SCIPfreeMemoryArray(scip, &detectordata->mergedconss);

   return SCIP_OKAY;
}

/** creates a new adjacency list for a vertex */
static
SCIP_RETCODE createAdjlist(
   SCIP*                 scip,
   ADJLIST**             adjlist
   )
{
   SCIP_CONS** conss;
   int* weights;
   int size;

   SCIP_CALL( SCIPallocMemory(scip, adjlist) );

   size = SCIPcalcMemGrowSize(scip, 1);
   SCIP_CALL( SCIPallocMemoryArray(scip, &conss, size) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &weights, size) );

   BMSclearMemoryArray(conss, size);
   BMSclearMemoryArray(weights, size);

   (*adjlist)->conss = conss;
   (*adjlist)->weights = weights;
   (*adjlist)->nconss = 0;
   (*adjlist)->maxconss = size;

   return SCIP_OKAY;
}

/** frees an adjacency list for a vertex */
static
SCIP_RETCODE freeAdjlist(
   SCIP*                 scip,
   ADJLIST**             adjlist
   )
{
   SCIPfreeMemoryArray(scip, &((*adjlist)->weights) );
   SCIPfreeMemoryArray(scip, &((*adjlist)->conss) );
   SCIPfreeMemory(scip, adjlist);
   *adjlist = NULL;

   return SCIP_OKAY;
}

/** increases an entry in an adjacency list */
static
SCIP_RETCODE adjlistIncreaseEntry(
   SCIP*                 scip,
   ADJLIST*              adjlist,
   SCIP_CONS*            cons,
   int                   incval
   )
{
   int i;

   for( i = 0; i < adjlist->nconss; ++i )
   {
      if( adjlist->conss[i] == cons )
      {
         adjlist->weights[i] += incval;
         return SCIP_OKAY;
      }
   }

   if( adjlist->nconss == adjlist->maxconss )
   {
      int newsize = SCIPcalcMemGrowSize(scip, adjlist->maxconss+1);
      assert(newsize > adjlist->maxconss);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(adjlist->conss), newsize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(adjlist->weights), newsize) );
      adjlist->maxconss = newsize;
   }

   assert(adjlist->nconss < adjlist->maxconss);
   adjlist->conss[i] = cons;
   adjlist->weights[i] = incval;
   ++adjlist->nconss;

   return SCIP_OKAY;
}

/** removes an entry from an adjacency list */
static
void adjlistRemoveEntry(
   ADJLIST*              adjlist,
   SCIP_CONS*            cons
   )
{
   int pos;
   int i;

   for( pos = 0; pos < adjlist->nconss; ++pos )
      if( adjlist->conss[pos] == cons )
         break;

   if( pos == adjlist->nconss )
      return;

   for( i = pos; i < adjlist->nconss-1; ++i )
   {
      adjlist->conss[i] = adjlist->conss[i+1]; /*lint !e679*/
      adjlist->weights[i] = adjlist->weights[i+1]; /*lint !e679*/
   }
   --adjlist->nconss;

   return;
}

/** for a given vertex, get an entry in its adjacency list */
static
int adjlistGetEntry(
   ADJLIST*              adjlist,
   SCIP_CONS*            cons
   )
{
   int i;

   for( i = 0; i < adjlist->nconss; ++i )
   {
      if( adjlist->conss[i] == cons )
         return adjlist->weights[i];
   }

   return 0;
}

/** copies an adjacency list into another */
static
SCIP_RETCODE copyAdjlist(
   SCIP*                 scip,
   ADJLIST*              source,
   ADJLIST*              target,
   ADJLIST*              linkadjlist,
   SCIP_Bool*            consslink,
   SCIP_HASHMAP*         constopos,
   SCIP_CONS*            sourcecons
   )
{
   int cost;
   int i;

   assert(target->nconss == 0);

   if( source->nconss > target->maxconss )
   {
      int newsize = SCIPcalcMemGrowSize(scip, source->nconss);
      assert(newsize > target->maxconss);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(target->conss), newsize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(target->weights), newsize) );
      target->maxconss = newsize;
   }

   cost = 0;
   for( i = 0; i < source->nconss; ++i )
   {
      int idx = (int) (size_t) SCIPhashmapGetImage(constopos, source->conss[i]); /*lint !e507*/

      if( consslink == NULL || !consslink[idx] )
      {
         target->conss[target->nconss] = source->conss[i];
         target->weights[target->nconss] = source->weights[i];
         ++target->nconss;
      }
      else if( consslink[idx] )
         cost += source->weights[i];
   }

   if( linkadjlist != NULL && cost > 0 )
   {
      SCIP_CALL( adjlistIncreaseEntry(scip, linkadjlist, sourcecons, cost) );
   }

   return SCIP_OKAY;
}

/** builds the graph from the given SCIP instance */
static
SCIP_RETCODE buildGraphStructure(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< detectordata data structure */
   )
{
   int i;
   int j;
   int k;
   GRAPH* graph;
   SCIP_CONS*** varinconss;

   graph = detectordata->graphs[0];

   SCIP_CALL( SCIPhashmapCreate(&(graph->constopos), SCIPblkmem(scip),detectordata->nrelconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &graph->adjlists, detectordata->nrelconss) );
   for( i = 0; i < detectordata->nrelconss; ++i )
   {
      SCIP_CALL( createAdjlist(scip, &graph->adjlists[i]) );
   }

   /* initialize constopos */
   assert( graph->nconss > 0 );
   for( i = 0; i < graph->nconss; ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(graph->constopos, graph->conss[i], (void*) (size_t) i) );
   }

   /* initialize adjacency list */

   varinconss = detectordata->varinconss;

   for( i = 0; i < detectordata->nrelvars; ++i )
   {
      for( j = 0; j < detectordata->nvarinconss[i]; ++j )
      {
         int idx;
         ADJLIST* adjlist;

         idx = (int) (size_t) SCIPhashmapGetImage(graph->constopos, varinconss[i][j]); /*lint !e507*/
         assert(idx < detectordata->nrelconss);
         adjlist = graph->adjlists[idx];

         for( k = j + 1; k < detectordata->nvarinconss[i]; ++k )
         {
            SCIP_CALL( adjlistIncreaseEntry(scip, adjlist, varinconss[i][k], 1) );
         }
      }
   }

   graph->cons1 = NULL;
   graph->cons2 = NULL;

   /* compute number of edges */
   graph->nedges = 0;
   for( i = 0; i < detectordata->nrelconss; ++i )
   {
      graph->nedges += graph->adjlists[i]->nconss;
   }

   detectordata->ngraphs = 1;

   return SCIP_OKAY;
}

/** inserts element into hashmap if it doesn't already exist */
static
SCIP_RETCODE hashmapInsert(
   SCIP_HASHMAP*         hm,                 /**< pointer to hashmap */
   void*                 origin,             /**< key to store */
   void*                 image               /**< image to store */
   )
{
   if( !SCIPhashmapExists(hm, origin) )
   {
      SCIP_CALL( SCIPhashmapInsert(hm, origin, image) );
   }

   return SCIP_OKAY;
}


/** builds a new adjacency list for a subgraph
 *  by copying the adjacency list from the current graph
 *  and merging the linking constraints into a single one
 */
static
SCIP_RETCODE buildNewAdjacencyList(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data data structure */
   GRAPH*                graph,              /**< current graph */
   GRAPH*                subgraph,           /**< subgraph whose adjacency list is to be built */
   int*                  partition,          /**< partition of the graph */
   int                   partidx,            /**< side of the partition the subgraph belongs to (0 or 1) */
   SCIP_Bool*            consslink,          /**< array of linking constraints */
   int                   nconsslink          /**< number of linking constraints */
   )
{
   int i;
   int j;
   SCIP_CONS* representative;                /* constraint representing the merged linking constraints */

   assert(consslink != NULL);

   representative = NULL;

   /* If there are linking constraints, they need to be merged into a single one */
   if( nconsslink > 0 )
   {
      int k;   /* counter for non-linking constraints */
      int l;   /* counter for linking constraints */

      k = 0;
      l = 0;
      subgraph->nedges = 0;

      for( i = 0; i < graph->nconss; ++i )
      {
         ADJLIST* adjlist;

         if( partition[i] != partidx )
            continue;

         adjlist = graph->adjlists[i];

         /* if the constraint is not a linking constraint, just copy it together with its adjacency list */
         if( !consslink[i] )
         {
            ++k;

            subgraph->conss[k] = graph->conss[i];
            SCIP_CALL( SCIPhashmapInsert(subgraph->constopos, graph->conss[i], (void*) (size_t) k) );

            SCIP_CALL( copyAdjlist(scip, adjlist, subgraph->adjlists[k], subgraph->adjlists[0], consslink, graph->constopos, graph->conss[i]) );
            subgraph->nedges += subgraph->adjlists[k]->nconss;
         }
         /* otherwise, the constraint will be merged with the other linking constraints;
          * therefore, add its incident edges to the merged adjacency list
          */
         else
         {
            representative = graph->conss[i];
            detectordata->mergedconss[detectordata->nrepresentatives][l] = graph->conss[i];
            ++l;

            for( j = 0; j < adjlist->nconss; ++j )
            {
               int idx = (int) (size_t) SCIPhashmapGetImage(graph->constopos, adjlist->conss[j]); /*lint !e507*/

               if( partition[idx] == partidx && !consslink[idx] )
               {
                  SCIP_CALL( adjlistIncreaseEntry(scip, subgraph->adjlists[0], adjlist->conss[j], adjlist->weights[j]) );
               }
            }
         }
      }
      assert(l == nconsslink);
      subgraph->conss[0] = representative;
      SCIP_CALL( SCIPhashmapInsert(subgraph->constopos, representative, (void*) (size_t) 0) );
      subgraph->nedges += subgraph->adjlists[0]->nconss;
      subgraph->nconss = k + 1;

      SCIP_CALL( SCIPreallocMemoryArray(scip, &detectordata->representatives, detectordata->nrepresentatives+5) );
      detectordata->nmergedconss[detectordata->nrepresentatives] = nconsslink;
      detectordata->representatives[detectordata->nrepresentatives] = representative;
      ++detectordata->nrepresentatives;

      /* free unnecessary adjacency lists */
      for( i = subgraph->nconss; i < subgraph->nconss + nconsslink - 1; ++i )
      {
         SCIP_CALL( freeAdjlist(scip, &subgraph->adjlists[i]) );
      }
   }
   /* If there are no linking constraints, just copy the vertices (constraints) and edges */
   else
   {
      int k = 0;

      subgraph->nedges = 0;
      for( i = 0; i < graph->nconss; ++i )
      {
         if( partition[i] != partidx )
            continue;

         subgraph->conss[k] = graph->conss[i];
         SCIP_CALL( SCIPhashmapSetImage(subgraph->constopos, graph->conss[i], (void*) (size_t) k) );

         SCIP_CALL( copyAdjlist(scip, graph->adjlists[i], subgraph->adjlists[k], NULL, NULL, NULL, NULL) );
         subgraph->nedges += subgraph->adjlists[k]->nconss;
      }
   }

   return SCIP_OKAY;
}

/** frees graph at position pos */
static
SCIP_RETCODE freeGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detectordata data structure */
   int                   pos,                /**< position of graph */
   int                   nconss              /**< number of vertices */
   )
{
   int i;

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( freeAdjlist(scip, &detectordata->graphs[pos]->adjlists[i]) );
   }
   SCIPfreeMemoryArray(scip, &detectordata->graphs[pos]->adjlists);
   SCIPhashmapFree(&detectordata->graphs[pos]->constopos);
   SCIPfreeMemoryArray(scip, &detectordata->graphs[pos]->conss);
   SCIPfreeMemory(scip, &detectordata->graphs[pos]);
   detectordata->graphs[pos] = NULL;

   return SCIP_OKAY;
}

/** allocates memory at position pos to facilitate saving a graph with nconss vertices */
static
SCIP_RETCODE allocateGraph(
   SCIP*                scip,             /**< SCIP data structure */
   DEC_DETECTORDATA*    detectordata,     /**< detectordata data structure */
   int                  pos,              /**< position in graph array */
   int                  nconss            /**< number of vertices */
   )
{
   int i;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata->graphs[pos]) ); /*lint !e866*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->graphs[pos]->conss, nconss) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->graphs[pos]->constopos, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->graphs[pos]->adjlists, nconss) );

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( createAdjlist(scip, &detectordata->graphs[pos]->adjlists[i]) );
   }
   detectordata->graphs[pos]->nconss = nconss;

   return SCIP_OKAY;
}

/** assigns the right linking constraint to the subgraph */
static
SCIP_RETCODE setLinkingCons(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                subgraph,           /**< subgraph */
   int                   cas,                /**< first parameter for case distiction */
   int                   cas2,               /**< second case distinction */
   SCIP_CONS*            cons1,              /**< first constraint */
   SCIP_CONS*            cons2               /**< second constraint */
   )
{
   assert(cas2 == 1 || cas2 == 2);

   switch( cas2 )
   {
   case 1:
      if( cas )
      {
         subgraph->cons1 = subgraph->conss[0];
         subgraph->cons2 = cons2;
      }
      else
      {
         subgraph->cons1 = cons1;
         subgraph->cons2 = subgraph->conss[0];
      }
      break;
   case 2:
      if( !cas )
      {
         subgraph->cons1 = subgraph->conss[0];
         subgraph->cons2 = cons1;
      }
      else
      {
         subgraph->cons1 = cons2;
         subgraph->cons2 = subgraph->conss[0];
      }
      break;
   default:
      break;
   }

   return SCIP_OKAY;
}

/** sets the startblock */
static
SCIP_RETCODE setStartBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detectordata data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   if( cons == NULL )
   {
      detectordata->startblock = detectordata->nblocks;
   }
   return SCIP_OKAY;
}

/** copies constraints of a the graph at position pos to subscipconss */
static
SCIP_RETCODE copyConss(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detectordata data structure */
   GRAPH*                graph,              /**< graph whose constraints should be copied */
   int*                  partition,          /**< partition of the graph, or NULL if all constraints are copied */
   int                   partidx             /**< side of the partition that is copied */
   )
{
   int nconss;
   int i;

   nconss = 0;
   for( i = 0; i < graph->nconss; ++i )
   {
      if( partition == NULL || partition[i] == partidx )
      {
         detectordata->subscipconss[detectordata->nblocks][nconss] = graph->conss[i];
         ++nconss;
      }
   }
   detectordata->nsubscipconss[detectordata->nblocks] = nconss;
   detectordata->nblocks++;

   return SCIP_OKAY;
}


/** builds the new graphs which result from the last found cut */
static
SCIP_RETCODE buildNewGraphs(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< detectordata data structure */
   )
{
   int i;
   int pos1;
   int pos2;
   int cas;
   int cons1idx;
   int cons2idx;
   SCIP_Bool stop1;
   SCIP_Bool stop2;

   int* partition;
   GRAPH* graph;

   int nconss1;
   int nconss2;
   SCIP_Bool* consslink;
   int nconsslink1;
   int nconsslink2;

   cas = -1;
   nconss1 = 0;
   nconss2 = 0;
   stop1 = FALSE;
   stop2 = FALSE;

   /* build partitions */
   partition = detectordata->partition;
   graph = detectordata->graphs[detectordata->position];

   /* obtain indices for storing the two new graphs */
   pos1 = -1;
   pos2 = -1;
   for( i = 0; i < detectordata->maxgraphs; ++i )
   {
      if( detectordata->graphs[i] == NULL )
      {
         if( pos1 == -1 )
            pos1 = i;
         else
         {
            pos2 = i;
            break;
         }
      }
   }
   assert(i < detectordata->maxgraphs);
   assert((pos1 != -1) && (pos2 != -1));

   /* for each of the two new subgraphs according to the partition,
    * count the number of constraints and linking constraints;
    * also, find the indices of cons1 and cons2, respectively in the graph
    */
   SCIP_CALL( SCIPallocMemoryArray(scip, &consslink, graph->nconss) );
   BMSclearMemoryArray(consslink, graph->nconss);
   nconsslink1 = 0;
   nconsslink2 = 0;
   cons1idx = -1;
   cons2idx = -1;
   for( i = 0; i < graph->nconss; ++i )
   {
      int j;

      assert(partition[i] == 0 || partition[i] == 1);

      if( partition[i] == 0 )
         nconss1++;
      else
         nconss2++;

      /* check if the constraint is a linking constraint */
      for( j = 0; j < graph->adjlists[i]->nconss; ++j )
      {
         int idx = (int) (size_t) SCIPhashmapGetImage(graph->constopos, graph->adjlists[i]->conss[j]); /*lint !e507*/
         if( partition[i] != partition[idx] )
         {
            if( !consslink[i] )
            {
               consslink[i] = TRUE;
               if( partition[i] == 0 )
                  ++nconsslink1;
               else
                  ++nconsslink2;
            }
            if( !consslink[idx] )
            {
               consslink[idx] = TRUE;
               if( partition[idx] == 0 )
                  ++nconsslink1;
               else
                  ++nconsslink2;
            }
         }
      }

      if( graph->conss[i] == graph->cons1 )
         cons1idx = i;
      if( graph->conss[i] == graph->cons2 )
         cons2idx = i;
   }
   assert(nconss1 + nconss2 == graph->nconss);
   assert((nconss1 != 0) && (nconss2 != 0));

   /* if subgraphs consist of linking constraints only, do not partition them further */
   if( nconsslink1 == nconss1 )
      stop1 = TRUE;
   if( nconsslink2 == nconss2 )
      stop2 = TRUE;

   /* test whether the cut is feasible */
   if( (graph->cons1 != NULL) && (graph->cons2 != NULL) )
   {
      assert(cons1idx != -1);
      assert(cons2idx != -1);

      if( partition[cons1idx] == partition[cons2idx] )
      {
         SCIP_CALL( copyConss(scip, detectordata, graph, NULL, -1) );
         SCIP_CALL( freeGraph(scip, detectordata, detectordata->position, graph->nconss) );
         detectordata->ngraphs--;
         SCIPfreeMemoryArray(scip, &consslink);
         return SCIP_OKAY;
      }

      if( partition[cons1idx] == 0 )
      {
         cas = 0;
         stop1 = partition[cons1idx] == 0 && consslink[cons1idx];
         stop2 = partition[cons2idx] == 1 && consslink[cons2idx];
      }
      else
      {
         cas = 1;
         stop1 = partition[cons2idx] == 0 && consslink[cons2idx];
         stop2 = partition[cons1idx] == 1 && consslink[cons1idx];
      }
   }

   /* test right or left */
   if( (graph->cons1 != NULL) && (graph->cons2 == NULL) )
   {
      assert(cons1idx != -1);
      assert(cons2idx == -1);

      if( partition[cons1idx] == 0 )
      {
         cas = 0;
         stop1 = consslink[cons1idx];
      }
      else
      {
         cas = 1;
         stop2 = consslink[cons1idx];
      }
   }
   else if( (graph->cons1 == NULL) && (graph->cons2 != NULL) )
   {
      assert(cons1idx == -1);
      assert(cons2idx != -1);

      if( partition[cons2idx] == 1 )
      {
         cas = 0;
         stop2 = consslink[cons2idx];
      }
      else
      {
         cas = 1;
         stop1 = consslink[cons2idx];
      }
   }
   else if( (graph->cons1 == NULL) && (graph->cons2 == NULL) )
      cas = 1;

   /* if the first subgraph is large enough, create it and build its adjacency list */
   if( (nconss1 > 1) && !stop1 )
   {
      SCIP_CALL( allocateGraph(scip, detectordata, pos1, nconss1) );
      SCIP_CALL( buildNewAdjacencyList(scip, detectordata, graph, detectordata->graphs[pos1], partition, 0, consslink, nconsslink1) );
      SCIP_CALL( setLinkingCons(scip, detectordata->graphs[pos1], cas, 1, graph->cons1, graph->cons2) );
      detectordata->ngraphs++;
   }
   else
   {
      SCIP_CALL( copyConss(scip, detectordata, graph, partition, 0) );
      SCIP_CALL( setStartBlock(scip, detectordata, cas ? graph->cons2 : graph->cons1) );
   }

   /* if the second subgraph is large enough, create it and build its adjacency list */
   if( (nconss2 > 1) && !stop2 )
   {
      SCIP_CALL( allocateGraph(scip, detectordata, pos2, nconss2) );
      SCIP_CALL( buildNewAdjacencyList(scip, detectordata, graph, detectordata->graphs[pos2], partition, 1, consslink, nconsslink2) );
      SCIP_CALL( setLinkingCons(scip, detectordata->graphs[pos2], cas, 2, graph->cons2, graph->cons1) );
      detectordata->ngraphs++;
   }
   else
   {
      SCIP_CALL( copyConss(scip, detectordata, graph, partition, 1) );
      SCIP_CALL( setStartBlock(scip, detectordata, cas ? graph->cons1 : graph->cons2) );
   }

   /* free the old graph */
   SCIP_CALL( freeGraph(scip, detectordata, detectordata->position, graph->nconss) );
   detectordata->ngraphs--;

   SCIPfreeMemoryArray(scip, &consslink);

   return SCIP_OKAY;
}

/** adds the merged constraints to the right blocks */
static
SCIP_RETCODE getMergedConss(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< detectordata data structure */
   )
{
   int i;
   int j;
   SCIP_CONS*** mergedconss;
   int* nmergedconss;
   SCIP_CONS** representatives;
   SCIP_CONS*** subscipconss;
   SCIP_HASHMAP* constoblock;

   mergedconss = detectordata->mergedconss;
   nmergedconss = detectordata->nmergedconss;
   representatives = detectordata->representatives;
   subscipconss = detectordata->subscipconss;
   constoblock = detectordata->constoblock;

   /* constoblock */
   for( i = 0; i < detectordata->nblocks; ++i )
   {
      for( j = 0; j < detectordata->nsubscipconss[i]; ++j )
      {
         SCIP_CALL( hashmapInsert(constoblock, subscipconss[i][j], (void *) (size_t) (i+1)) );
      }
   }

   for( i = detectordata->nrepresentatives - 1; i >= 0; --i )
   {
      SCIP_CONS* cons;
      int block;

      cons = representatives[i];
      block = (int) (size_t) SCIPhashmapGetImage(constoblock, cons); /*lint !e507*/

      for( j = 0; j < nmergedconss[i]; ++j )
      {
         SCIP_CONS* cons2 = mergedconss[i][j];

         if( cons != cons2 )
         {
            subscipconss[block - 1][detectordata->nsubscipconss[block - 1]] = cons2;
            detectordata->nsubscipconss[block - 1]++;
            SCIP_CALL( hashmapInsert(constoblock, cons2, (void*) (size_t) block) );
         }

      }
   }
   assert(SCIPhashmapGetNEntries(constoblock) == detectordata->nrelconss);

#ifndef NDEBUG
   j = 0;
   for( i = 0; i < detectordata->nblocks; ++i )
   {
      j += detectordata->nsubscipconss[i];
   }
   assert( j == detectordata->nrelconss );
#endif

   return SCIP_OKAY;
}

/** arranges the constraints as prescribed by the cuts */
static
SCIP_RETCODE getConsIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data structure */
   DEC_DECOMP*           decdecomp           /**< decdecomp pointer */
)
{
   int i;
   int j;
   int k;
   int block;
   int no;

   int newblock;
   int oldblock;
   int actblock;
   int counter;

   SCIP_CONS*** newsubscipconss;
   int* nnewsubscipconss;
   SCIP_VAR*** stairlinkingvars;
   int* nstairlinkingvars;
   SCIP_VAR** linkingvars;


   newblock = 0;
   oldblock = 0;

   SCIP_CALL( SCIPallocMemoryArray(scip, &stairlinkingvars, detectordata->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &linkingvars, detectordata->nrelvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nstairlinkingvars, detectordata->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &newsubscipconss, detectordata->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nnewsubscipconss, detectordata->nblocks) );

   for( i = 0; i < detectordata->nblocks; ++i )
   {
      SCIP_VAR** tmpvars;

      SCIP_CALL( SCIPallocMemoryArray(scip, &tmpvars, detectordata->nrelvars) );
      stairlinkingvars[i] = tmpvars;
      nstairlinkingvars[i] = 0;
      nnewsubscipconss[i] = 0;
   }

   block = detectordata->startblock;
   for( counter = detectordata->nblocks; counter > 0; --counter )
   {
      SCIP_CONS** conss;
      int nconss;
      SCIP_CONS** newconss;

      conss = detectordata->subscipconss[block-1];
      nconss = detectordata->nsubscipconss[block-1];

      for( i = 0; i < nconss; ++i )
      {
         SCIP_VAR** consvars;
         int nconsvars;

         nconsvars = GCGconsGetNVars(scip, conss[i]);
         SCIP_CALL( SCIPallocMemoryArray(scip, &consvars, nconsvars) );
         SCIP_CALL( GCGconsGetVars(scip, conss[i], consvars, nconsvars) );

         for( j = 0; j < nconsvars; ++j )
         {
            if( (!GCGisVarRelevant(consvars[j])) && (j < nconsvars) )
               continue;

            no = (int) (size_t) SCIPhashmapGetImage(detectordata->vartopos, SCIPvarGetProbvar(consvars[j])); /*lint !e507*/

            for( k = 0; (k < detectordata->nvarinconss[no]); ++k )
            {
               actblock = (int) (size_t) SCIPhashmapGetImage(detectordata->constoblock, detectordata->varinconss[no][k]); /*lint !e507*/
               if( actblock != block && actblock != oldblock )
                  newblock = actblock;
            }
         }

         SCIPfreeMemoryArrayNull(scip, &consvars);
      }
      assert((newblock != block) || (counter < 2));

      SCIP_CALL( SCIPallocMemoryArray(scip, &newconss, nconss) );
      for( j = 0; j < nconss; ++j )
         newconss[j] = conss[j];

      newsubscipconss[detectordata->nblocks-counter] = newconss; /*lint !e679*/
      nnewsubscipconss[detectordata->nblocks-counter] = nconss; /*lint !e679*/

      oldblock = block;
      block = newblock;
      assert(0 < block);
      assert(block <= detectordata->nblocks);
   }

   for( i = 0; i < detectordata->nblocks; ++i )
   {
      for( j = 0; j < nnewsubscipconss[i]; ++j )
      {
         SCIP_CALL( SCIPhashmapSetImage(detectordata->constoblock, newsubscipconss[i][j], (void*) (size_t) (i+1)) );
      }
   }

   if( !detectordata->fixedblocks )
   {

      SCIP_CALL( DECfilloutDecompFromConstoblock(scip, decdecomp, detectordata->constoblock, detectordata->nblocks, TRUE) );

      for( i = 0; i < detectordata->nblocks; ++i )
      {
         SCIPfreeMemoryArray(scip, &newsubscipconss[i]);
      }
      SCIPfreeMemoryArray(scip, &newsubscipconss);
      SCIPfreeMemoryArray(scip, &nnewsubscipconss);
      SCIPfreeMemoryArray(scip, &linkingvars);
   }
   for( i = 0; i < detectordata->nblocks; ++i )
   {
      SCIPfreeMemoryArray(scip, &stairlinkingvars[i]);
   }
   SCIPfreeMemoryArray(scip, &stairlinkingvars);
   SCIPfreeMemoryArray(scip, &nstairlinkingvars);

   return SCIP_OKAY;
}


/** gets the variables which are linking */
static
SCIP_RETCODE getLinkingVars(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data structure */
   DEC_DECOMP*           decdecomp           /**< decdecomp pointer */
   )
{
   int i;
   int j;
   int newblock;

   SCIP_VAR** vars;
   int nvars;
   SCIP_VAR** linkingvars;
   int nlinkingvars;
   SCIP_CONS*** varinconss;
   int* nvarinconss;
   SCIP_VAR*** subscipvars;
   int* nsubscipvars;
   SCIP_HASHMAP* vartoblock;

   varinconss = detectordata->varinconss;
   nvarinconss = detectordata->nvarinconss;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &nsubscipvars, detectordata->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &subscipvars, detectordata->nblocks) );
   SCIP_CALL( SCIPhashmapCreate(&vartoblock, SCIPblkmem(scip),detectordata->nrelvars) );
   for( i = 0; i < detectordata->nblocks; ++i )
   {
      SCIP_VAR** subvars;

      SCIP_CALL( SCIPallocMemoryArray(scip, &subvars, detectordata->nrelvars) );
      subscipvars[i] = subvars;
      nsubscipvars[i] = 0;
   }

   /* bulid linkingvars, subscipvars and vartoblock */

   nlinkingvars = 0;
   SCIP_CALL( SCIPallocMemoryArray(scip, &linkingvars, detectordata->nrelvars) );

   for( i = 0; i < nvars; ++i )
   {
      int varpos;
      int oldblock;
      SCIP_Bool stop;

      if( !GCGisVarRelevant(vars[i]) )
         continue;

      varpos = (int) (size_t) SCIPhashmapGetImage(detectordata->vartopos, vars[i]); /*lint !e507*/
      oldblock = (int) (size_t) SCIPhashmapGetImage(detectordata->constoblock, varinconss[varpos][0]); /*lint !e507*/
      stop = FALSE;

      for( j = 1; !stop && (j < nvarinconss[varpos]); ++j )
      {
         newblock = (int) (size_t) SCIPhashmapGetImage(detectordata->constoblock, varinconss[varpos][j]); /*lint !e507*/

         if( newblock != oldblock )
            stop = TRUE;
      }

      if( stop )
      {
         linkingvars[nlinkingvars] = vars[i];
         ++nlinkingvars;
         SCIP_CALL( SCIPhashmapInsert(vartoblock, vars[i], (void*) (size_t) (detectordata->nblocks+1)) );
      }
      else
      {
         SCIP_CALL( SCIPhashmapInsert(vartoblock, vars[i], (void*) (size_t) oldblock) );
         subscipvars[oldblock-1][nsubscipvars[oldblock-1]] = vars[i];
         ++nsubscipvars[oldblock-1];
      }
   }

#ifdef SCIP_DEBUG
   j = 0;
   for( i = 0; i < detectordata->nblocks; ++i )
   {
      j += nsubscipvars[i];
   }
   SCIPdebugMessage("var %d %d %d \n",nlinkingvars, j, detectordata->nrelvars);
   assert(nlinkingvars +j == detectordata->nrelvars);
   assert(detectordata->nrelvars == SCIPhashmapGetNEntries(vartoblock));
#endif

   SCIP_CALL( DECdecompSetSubscipvars(scip, decdecomp, subscipvars, nsubscipvars) );
   SCIP_CALL( DECdecompSetLinkingvars(scip, decdecomp, linkingvars, nlinkingvars, 0, 0) );
   DECdecompSetVartoblock(decdecomp, vartoblock);

   for( i = 0; i < detectordata->nblocks; ++i )
   {
      SCIPfreeMemoryArrayNull(scip, &subscipvars[i]);
   }
   SCIPfreeMemoryArray(scip, &subscipvars);
   SCIPfreeMemoryArray(scip, &nsubscipvars);
   SCIPfreeMemoryArray(scip, &linkingvars);

   return SCIP_OKAY;
}

/** assigns the constraints and variables to blocks if a fixed blocksize is given */
static
SCIP_RETCODE fixedBlocks(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data structure */
   DEC_DECOMP*           decdecomp           /**< decdecomp pointer */
   )
{
   int block;
   int blocksize;
   int indexcons;
   int newblock;

   SCIP_CONS** conss;
   int nconss;
   SCIP_HASHMAP* consindex;

   int i;

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   SCIP_CALL( SCIPhashmapCreate(&consindex, SCIPblkmem(scip), detectordata->nrelconss) );

   /* get nblocks */
   blocksize = detectordata->blocksize;
   block = detectordata->nrelconss/blocksize + 1;
   if( detectordata->nrelconss/blocksize * blocksize == detectordata->nrelconss )
      block = detectordata->nrelconss/blocksize;
   assert(block > 0);

   detectordata->nblocks = block;

   SCIPdebugMessage("nblocks %d; nconss %d; blocksize %d; \n", block, detectordata->nrelconss, blocksize);

   SCIP_CALL( SCIPhashmapRemoveAll(detectordata->constoblock) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &detectordata->nsubscipconss, block) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &detectordata->subscipconss, block) );
   for( i = 0; i < block; ++i )
   {
      SCIP_CONS** subconss;

      subconss = detectordata->subscipconss[i];
      SCIP_CALL( SCIPreallocMemoryArray(scip, &subconss, blocksize) ); /*lint !e522*/
      detectordata->subscipconss[i] = subconss;
      detectordata->nsubscipconss[i] = 0;
   }

   /* build subscipcons and constoblock */
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPhashmapExists(consindex, conss[i]) )
      {
         assert( SCIPhashmapExists(consindex, conss[i]) );

         indexcons = (int) (size_t) SCIPhashmapGetImage(consindex, conss[i]); /*lint !e507*/
         newblock = indexcons/blocksize + 1;
         if( indexcons/blocksize * blocksize == indexcons )
            newblock = indexcons/blocksize;

         assert(newblock > 0);
         assert(newblock < block+1);

         SCIP_CALL( SCIPhashmapInsert(detectordata->constoblock, conss[i], (void*) (size_t) newblock) );

         assert(detectordata->subscipconss[newblock-1] != NULL);
         assert(detectordata->nsubscipconss != NULL);

         detectordata->subscipconss[newblock-1][detectordata->nsubscipconss[newblock-1]] = conss[i];
         ++detectordata->nsubscipconss[newblock-1];
      }
   }

   assert( SCIPhashmapGetNEntries(detectordata->constoblock) == detectordata->nrelconss);

   SCIP_CALL( DECfilloutDecompFromConstoblock(scip, decdecomp, detectordata->constoblock, detectordata->nblocks, TRUE) );

   return SCIP_OKAY;
}

/** will find a minimum cut via the Stoer-Wagner algorithm */
static
SCIP_RETCODE applyStoerWagner(
   SCIP*                 scip,               /**< SCIP data struture */
   DEC_DETECTORDATA*     detectordata        /**< detector data data structure */
   )
{
   int i;
   int j;
   int lastpos;
   int nextpos;
   int tight;
   SCIP_Real cutval;
   int* tightness;
   int ntight;
   SCIP_Bool* merged;
   int nmerged;
   SCIP_HASHMAP* represconss;
   ADJLIST** adjlists;
   SCIP_CONS** mincut;
   int nmincut;
   int nrepresconss;
   SCIP_CONS*** mergedconss;
   int* nmergedconss;
   GRAPH* graph;
   SCIP_CONS* s;
   SCIP_CONS* t;
   SCIP_CONS* cut;
   SCIP_CONS* last;
   SCIP_CONS* nexttolast;
   SCIP_CONS* representative;

   graph = detectordata->graphs[detectordata->position];
   nrepresconss = 1;

   cut = NULL;

   SCIP_CALL( SCIPallocMemoryArray(scip, &merged, graph->nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &tightness, graph->nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mincut, graph->nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nmergedconss, graph->nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mergedconss, graph->nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &adjlists, graph->nconss) );
   for( i = 0; i < graph->nconss; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(mergedconss[i]), graph->nconss) ); /*lint !e866*/
      SCIP_CALL( createAdjlist(scip, &adjlists[i]) );
   }
   SCIP_CALL( SCIPhashmapCreate(&represconss, SCIPblkmem(scip), graph->nconss) );

   /* copy adjacency lists */
   for( i = 0; i < graph->nconss; ++i )
   {
      SCIP_CALL( copyAdjlist(scip, graph->adjlists[i], adjlists[i], NULL, NULL, graph->constopos, NULL) );
   }

   SCIPdebugMessage("apply Stoer-Wagner...\n");

   t = NULL;
   cutval = SCIPinfinity(scip);
   if( graph->cons1 != NULL )
   {
      s = graph->cons1;
      if( graph->cons2 != NULL )
         t = graph->cons2;
   }
   else
   {
      if( graph->cons2 != NULL )
         s = graph->cons2;
      else
         s = graph->conss[0];
   }

   assert(s != NULL);
   representative = t;
   BMSclearMemoryArray(merged, graph->nconss);
   nmerged = 0;

   while( nmerged < graph->nconss - 1 )
   {
      SCIP_Real curcutval;

      last = s;
      nexttolast = s;
      lastpos = (int) (size_t) SCIPhashmapGetImage(graph->constopos, last); /*lint !e507*/
      nextpos = lastpos;

      BMSclearMemoryArray(tightness, graph->nconss);
      tightness[lastpos] = -1;
      ntight = 1;

      while( ntight < graph->nconss - nmerged )
      {
         nexttolast = last;
         nextpos = lastpos;

         for( j = 0; j < adjlists[nextpos]->nconss; ++j )
         {
            int idx;

            idx = (int) (size_t) SCIPhashmapGetImage(graph->constopos, adjlists[nextpos]->conss[j]); /*lint !e507*/
            assert(!merged[idx]);

            if( tightness[idx] != -1 )
               tightness[idx] += adjlists[nextpos]->weights[j];
         }
         for( j = 0; j < graph->nconss; ++j )
         {
            if( !merged[j] && tightness[j] != -1 )
               tightness[j] += adjlistGetEntry(adjlists[j], nexttolast);
         }

         /* choose the most tight */
         tight = 0;
         for( j = 0; j < graph->nconss; ++j )
         {
            if( !merged[j] && tightness[j] >= tight )
            {
               tight = tightness[j];
               last = graph->conss[j];
               lastpos = j;
            }
         }

         assert(nexttolast != last);
         assert(tightness[lastpos] != -1);
         tightness[lastpos] = -1;
         ++ntight;
      }

      /* calculate the value of the current cut */
      curcutval = 0;
      for( j = 0; j < adjlists[lastpos]->nconss; ++j )
      {
         int idx = (int) (size_t) SCIPhashmapGetImage(graph->constopos, adjlists[lastpos]->conss[j]); /*lint !e507*/
         assert(adjlists[lastpos]->conss[j] != last);
         if( !merged[idx] )
            curcutval += adjlists[lastpos]->weights[j];
      }
      for( j = 0; j < graph->nconss; ++j )
         if( !merged[j] )
            curcutval += adjlistGetEntry(adjlists[j], last);

      if( curcutval < cutval )
      {
         if( (t != NULL) && (representative == last) )
         {
            /* test if act_cut ist a s-t-cut */
            representative = nexttolast;
            cutval = curcutval;
            cut = last;
            assert(cutval != 0);
         }
         else if( t == NULL )
         {
            cutval = curcutval;
            cut = last;
            assert(cutval != 0);
         }
      }

      /* merging */
      if( SCIPhashmapGetNEntries(represconss) == 0 )
      {
         SCIP_CALL( SCIPhashmapInsert(represconss, nexttolast, (void*) (size_t) nrepresconss) );
         nrepresconss++; /* starts with 1 */
         nmergedconss[nrepresconss - 2] = 1;
         mergedconss[nrepresconss - 2][0] = last;
      }
      else if( !SCIPhashmapExists(represconss, nexttolast) )
      {
         SCIP_CALL( SCIPhashmapInsert(represconss, nexttolast, (void*) (size_t) nrepresconss) );
         nrepresconss++; /* starts with 1 */
         nmergedconss[nrepresconss - 2] = 1;
         mergedconss[nrepresconss - 2][0] = last;
      }
      else
      {
         int idx = (int) (size_t) SCIPhashmapGetImage(represconss, nexttolast); /*lint !e507*/
         mergedconss[idx-1][nmergedconss[idx-1]] = last;
         nmergedconss[idx-1]++;
      }

      /* connect last and nexttolast */
      for( j = 0; j < adjlists[lastpos]->nconss; ++j )
      {
         int idx = (int) (size_t) SCIPhashmapGetImage(graph->constopos, adjlists[lastpos]->conss[j]); /*lint !e507*/
         assert(!merged[idx]);

         if( adjlists[lastpos]->conss[j] != nexttolast )
         {
            if( adjlistGetEntry(adjlists[idx], nexttolast) > 0 )
            {
               SCIP_CALL( adjlistIncreaseEntry(scip, adjlists[idx], nexttolast, adjlists[lastpos]->weights[j]) );
            }
            else
            {
               SCIP_CALL( adjlistIncreaseEntry(scip, adjlists[nextpos], adjlists[lastpos]->conss[j], adjlists[lastpos]->weights[j]) );
            }
         }
      }
      for( j = 0; j < graph->nconss; ++j )
      {
         int lastweight;

         if( merged[j] )
            continue;

         lastweight = adjlistGetEntry(adjlists[j], last);
         if( lastweight > 0 && graph->conss[j] != nexttolast )
         {
            if( adjlistGetEntry(adjlists[j], nexttolast) > 0 )
            {
               SCIP_CALL( adjlistIncreaseEntry(scip, adjlists[j], nexttolast, lastweight) );
            }
            else
            {
               SCIP_CALL( adjlistIncreaseEntry(scip, adjlists[nextpos], graph->conss[j], lastweight) );
            }
         }
      }

      merged[lastpos] = TRUE;
      ++nmerged;

      /* delete 'last' in adjacency list */
      for( i = 0; i < graph->nconss; ++i )
         if( !merged[i] )
            adjlistRemoveEntry(adjlists[i], last);
   }

   /* complete cut */
   nmincut = 1;
   mincut[0] = cut;
   for( i = 0; i < nmincut; ++i )
   {
      if( SCIPhashmapExists(represconss, mincut[i]) )
      {
         int idx;

         idx = (int) (size_t) SCIPhashmapGetImage(represconss, mincut[i]); /*lint !e507*/

         for( j = 0; j < nmergedconss[idx-1]; ++j )
         {
            mincut[nmincut] = mergedconss[idx-1][j];
            nmincut++;
         }
      }
   }

   /* create partition */
   for( i = 0; i < graph->nconss; ++i )
      detectordata->partition[i] = 0;
   for( i = 0; i < nmincut; ++i )
      detectordata->partition[(int) (size_t) SCIPhashmapGetImage(graph->constopos, mincut[i])] = 1; /*lint !e507*/

   SCIPfreeMemoryArray(scip, &tightness);
   SCIPfreeMemoryArray(scip, &merged);
   for( i = 0; i < graph->nconss; ++i )
   {
      SCIP_CALL( freeAdjlist(scip, &adjlists[i]) );
      SCIPfreeMemoryArray(scip, &(mergedconss[i]));
   }
   SCIPfreeMemoryArrayNull(scip, &adjlists);
   SCIPhashmapFree(&represconss);
   SCIPfreeMemoryArray(scip, &mincut);
   SCIPfreeMemoryArray(scip, &nmergedconss);
   SCIPfreeMemoryArray(scip, &mergedconss);

   return SCIP_OKAY;
}


/** will call hmetis via a system call */
static
SCIP_RETCODE callMetis(
   SCIP*                 scip,               /**< SCIP data struture */
   DEC_DETECTORDATA*     detectordata,       /**< detector data data structure */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   char metiscall[SCIP_MAXSTRLEN];
   char metisout[SCIP_MAXSTRLEN];
   char line[SCIP_MAXSTRLEN];
   char tempfile[SCIP_MAXSTRLEN];

   int i;
   int j;
   int status;
   int nvertices;
   int nedges;
   int entry;
   int cost;
   int* partition;
   ADJLIST** adjlists;
   SCIP_HASHMAP* constopos;
   SCIP_FILE *zfile;
   FILE* file;
   int temp_filedes;

   assert(scip != NULL);
   assert(detectordata != NULL);

   *result = SCIP_DIDNOTRUN;

   adjlists = detectordata->graphs[detectordata->position]->adjlists;
   constopos = detectordata->graphs[detectordata->position]->constopos;
   nvertices = detectordata->graphs[detectordata->position]->nconss;
   nedges = detectordata->graphs[detectordata->position]->nedges;

   assert(adjlists != NULL);
   assert(constopos != NULL);

   (void) SCIPsnprintf(tempfile, SCIP_MAXSTRLEN, "gcg.metis.XXXXXX");
   temp_filedes = mkstemp(tempfile);
   if( temp_filedes < 0 )
   {
      SCIPerrorMessage("Error creating temporary file: %s\n", strerror(errno));
      return SCIP_FILECREATEERROR;
   }

   SCIPdebugMessage("temporary filename: %s\n", tempfile);

   file = fdopen(temp_filedes, "w");
   if( file == NULL )
   {
      SCIPerrorMessage("Could not open temporary metis file!\n");
      return SCIP_FILECREATEERROR;
   }

   SCIPinfoMessage(scip, file, "%d %d 1\n", nedges, nvertices);

   for( i = 0; i < nvertices - 1; ++i )
   {
      assert(adjlists[i] != NULL);

      for( j = 0; j < adjlists[i]->nconss; ++j )
      {
         entry = (int) (size_t) SCIPhashmapGetImage(constopos, adjlists[i]->conss[j]); /*lint !e507*/
         cost = adjlists[i]->weights[j];

         SCIPinfoMessage(scip, file, "%d ", cost);
         if( entry > i )
         {
            SCIPinfoMessage(scip, file, "%d ", entry + 1);
            SCIPinfoMessage(scip, file, "%d \n", i + 1);
         }
         else
         {
            SCIPinfoMessage(scip, file, "%d ", i + 1);
            SCIPinfoMessage(scip, file, "%d \n", entry + 1);
         }

      }
   }
   status = fclose(file);

   if( status == -1 )
   {
      SCIPerrorMessage("Could not close '%s'\n", tempfile);
      return SCIP_WRITEERROR;
   }

   (void) SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"" HMETIS_EXECUTABLE " %s %d -seed %d -ptype %s -ufactor %f %s\"", tempfile, 2,
      detectordata->randomseed, detectordata->metisuseptyperb ? "rb" : "kway", detectordata->metisubfactor,
      detectordata->metisverbose ? "" : "> /dev/null");

   /* SCIP_CALL( SCIPresetClock(scip, detectordata->metisclock) ); */
   /* SCIP_CALL( SCIPstartClock(scip, detectordata->metisclock) ); */
   SCIPdebugMessage("Calling metis with: %s\n", metiscall);

   status = system(metiscall);

   /* SCIP_CALL( SCIPstopClock(scip, detectordata->metisclock) ); */
   /* SCIPdebugMessage("time left before metis started: %f, time metis spend %f, remainingtime: %f\n", remainingtime, SCIPgetClockTime(scip, detectordata->metisclock),  remainingtime-SCIPgetClockTime(scip, detectordata->metisclock) ); */

   /* check error codes */
   if( status == -1 )
   {
      SCIPerrorMessage("System call did not succed: %s\n", strerror(errno));
      SCIPerrorMessage
      ("Call was %s\n", metiscall);
   }
   else if( status != 0 )
   {

      SCIPerrorMessage("Calling hmetis unsuccessful! See the above error message for more details.\n");
      SCIPerrorMessage
      ("Call was %s\n", metiscall);
   }

   /* exit gracefully in case of errors */
   if( status != 0 )
   {
      if( detectordata->tidy )
      {
         status = unlink(tempfile);
         if( status == -1 )
         {
            SCIPerrorMessage("Could not remove metis input file: ", strerror(errno));
            return SCIP_WRITEERROR;
         }
      }
      return SCIP_ERROR;
   }

   /*
    * parse the output into the vector
    * alloc the memory
    */
   if( detectordata->partition == NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->partition, nvertices) );
   }

   assert(detectordata->partition != NULL);
   partition = detectordata->partition;

   (void) SCIPsnprintf(metisout, SCIP_MAXSTRLEN, "%s.part.%d", tempfile, 2);

   zfile = SCIPfopen(metisout, "r");
   i = 0;
   while( !SCIPfeof(zfile) && i < nvertices )
   {
      int temp;
      if( SCIPfgets(line, SCIP_MAXSTRLEN, zfile) == NULL )
      {
         SCIPerrorMessage("Line could not be read\n");
         return SCIP_READERROR;
      }

      temp = atoi(line);
      assert( temp >= 0 && temp <= 1 );
      partition[i] = temp;
      ++i;
   }
   SCIPfclose(zfile);

   /* if desired delete the temoprary metis file */
   if( detectordata->tidy )
   {
      status = unlink(tempfile);
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis input file: %s\n", strerror(errno));
         return SCIP_WRITEERROR;
      }
      status = unlink(metisout);
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis output file: %s\n", strerror(errno));
         return SCIP_WRITEERROR;
      }
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "Temporary file is in: %s\n", tempfile);
   }
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * detector callback methods
 */

/** detector deinitialization method (called before the transformed problem is freed) */
static
DEC_DECL_FREEDETECTOR(detectorFreeCutpacking)
{
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);
   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** detector initialization method (called after problem was transformed) */
static
DEC_DECL_INITDETECTOR(detectorInitCutpacking)
{
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);
   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   detectordata->partition = NULL;
   detectordata->nblocks = -1;

   return SCIP_OKAY;
}

/** detector structure detection method, tries to detect a structure in the problem */
static
DEC_DECL_DETECTSTRUCTURE(detectorDetectCutpacking)
{
   int i;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(decdecomps != NULL);
   assert(ndecdecomps != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting cutpacking structure: ");

   SCIP_CALL( initData(scip, detectordata) );
   SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, 1) ); /*lint !e506*/
   *ndecdecomps = 1;

   /* build the hypergraph structure from the original problem */
   SCIP_CALL( buildGraphStructure(scip, detectordata) );
   SCIPdebugMessage("Building graph structure successful.\n");

   SCIP_CALL( DECdecompCreate(scip, &(*decdecomps)[0]) );

   /* get the partitions for the new variables from metis */
   while( detectordata->ngraphs > 0 )
   {
      for( i = 0; i < detectordata->maxgraphs; i++ )
      {
         if( detectordata->graphs[i] != NULL )
         {
            detectordata->position = i;

            if( detectordata->usemetis )
            {
               SCIP_CALL( callMetis(scip, detectordata,result) );
               SCIPdebugMessage("  -> metis successful.\n");
            }
            else
            {
               SCIP_CALL( applyStoerWagner(scip, detectordata) );
               SCIPdebugMessage("  -> Stoer-Wagner successful.\n");
            }

            SCIPdebugMessage("Creating two new graphs according to partition...\n");
            SCIP_CALL( buildNewGraphs(scip, detectordata) );
            SCIPdebugMessage("  -> buildNewGraphs successful.\n");
         }
      }
   }
   /* add merged conss */
   SCIP_CALL( getMergedConss(scip, detectordata) );
   SCIPdebugMessage("getMergedConss successful.\n");

   /* get subscipvars, copy data to decdecomp */
   SCIP_CALL( getConsIndex(scip, detectordata, (*decdecomps)[0]) );

   if( detectordata->fixedblocks )
   {
        SCIP_CALL( fixedBlocks(scip, detectordata, (*decdecomps)[0]) );
        SCIP_CALL( getLinkingVars(scip, detectordata, (*decdecomps)[0]) );
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "found %d blocks.\n", DECdecompGetNBlocks((*decdecomps)[0]));

   SCIP_CALL( freeData(scip, detectordata) );

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}
#endif

#define detectorPropagateSeeedCutpacking NULL
#define detectorFinishSeeedCutpacking NULL
#define detectorExitCutpacking NULL

#define detectorPostprocessSeeedCutpacking NULL

#define setParamAggressiveCutpacking NULL
#define setParamDefaultCutpacking NULL
#define setParamFastCutpacking NULL



/** creates the cutpacking detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorCutpacking(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
#if !defined(_WIN32) && !defined(_WIN64)
   DEC_DETECTORDATA *detectordata;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);


   /* include structure detector */
   SCIP_CALL( DECincludeDetector(scip,
      DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDORIGINAL, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, DEC_LEGACYMODE,
      detectordata, detectorDetectCutpacking, detectorFreeCutpacking, detectorInitCutpacking, detectorExitCutpacking, detectorPropagateSeeedCutpacking, NULL, NULL, detectorFinishSeeedCutpacking, detectorPostprocessSeeedCutpacking, setParamAggressiveCutpacking, setParamDefaultCutpacking, setParamFastCutpacking) );


   /* add cutpacking detector parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/cutpacking/algorithm",
      "should the Stoer-Wagner algorithm or metis be used for finding a minimal cut",
      &detectordata->usemetis, FALSE, DEFAULT_USEMETIS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/cutpacking/fixedblocks",
      "Should the blocks consist of a certain number of constraints",
      &detectordata->fixedblocks, FALSE, DEFAULT_FIXEDBLOCKS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/cutpacking/blocksize",
      "number of constraints per block",
      &detectordata->blocksize, FALSE, DEFAULT_BLOCKSIZE, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/cutpacking/tidy",
      "Whether to clean up temporary files",
      &detectordata->tidy, FALSE, DEFAULT_TIDY, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/cutpacking/randomseed",
      "random seed for hmetis",
      &detectordata->randomseed, FALSE, DEFAULT_RANDSEED, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detection/detectors/cutpacking/ubfactor",
      "Unbalance factor for metis",
      &detectordata->metisubfactor, FALSE, DEFAULT_METIS_UBFACTOR, 0.0, 1E20, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/cutpacking/metisverbose",
      "Should the metis output be displayed",
      &detectordata->metisverbose, FALSE, DEFAULT_METIS_VERBOSE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/cutpacking/metisuseptyperb",
      "Should the rb or kway method be used for partitioning by metis",
      &detectordata->metisuseptyperb, FALSE, DEFAULT_METISUSEPTYPE_RB, NULL, NULL) );
#endif

   return SCIP_OKAY;
}
