/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2014 Operations Research, RWTH Aachen University       */
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
 *
 * This detector tries to detect staircase structures by recursively partitioning the
 * rowgraph of the matrix by using hmetis.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#include "dec_cutpacking.h"

#include "cons_decomp.h"
#include "struct_decomp.h"
#include "pub_decomp.h"
#include "scip_misc.h"

#define DEC_DETECTORNAME      "cutpacking"   /**< name of the detector */
#define DEC_DESC              "detects staircase matrices via graph partioning and cutpacking" /**< detector description */
#define DEC_PRIORITY          1100           /**< priority of the detector */
#define DEC_DECCHAR           'c'            /**< display character of detector */
#define DEC_ENABLED           FALSE          /**< should detector be called by default */
#define DEC_SKIP              FALSE          /**< should detector be skipped if others found detections */

/* Default parameter settings */
#define DEFAULT_RANDSEED                  1     /**< random seed for the hmetis call */
#define DEFAULT_TIDY                     TRUE   /**< whether to clean up afterwards */
#define DEFAULT_FIXEDBLOCKS              FALSE   /**< whether the blocks should consist of a given number of constraints */
#define DEFAULT_BLOCKSIZE                 200    /**< number of constraints per block */
#define DEFAULT_ALGORITHM_METIS           1     /**< Should be used metis or the stoer-wagner algorithm */

#define DEFAULT_METIS_UBFACTOR            5.0   /**< default unbalance factor given to metis on the commandline */
#define DEFAULT_METIS_VERBOSE             FALSE /**< should metis be verbose */
#define DEFAULT_METISUSEPTYPE_RB          TRUE  /**< Should metis use the rb or kway partitioning algorithm */

/*
 * Data structures
 */

/** graphstructure **/
struct graphstructure
{
   SCIP_HASHMAP** adjacencylist; /**< adjacencylists of the graph */
   SCIP_CONS** conss;            /**< constraints (each constraint represents a vertice of the graph)*/
   int nconss;                   /**< number of vertices */
   SCIP_HASHMAP* constopos;      /**< assigns constraints to their poistion in conss */

   int nedges;                   /**< number of edges */

   SCIP_CONS* cons1;             /**< */
   SCIP_CONS* cons2;             /**< */
};
typedef struct graphstructure Graph;

/** detector data */
struct DEC_DetectorData
{
   int iter;
   int algorithm;
   int blocksize;
   SCIP_Bool fixedblocks;

   /* Stuff for the algorithm */
   Graph* graphs;
   int ngraphs;
   SCIP_CONS*** subscipconss;
   int* nsubscipconss;
   SCIP_HASHMAP* constoblock;

   int nblocks;
   SCIP_HASHMAP* occupied;
   int position;
   int startblock;
   int* partition;

   SCIP_HASHMAP** mergedconss;
   SCIP_HASHMAP* representatives;
   int nrepresentatives;

   /* general stuff */
   SCIP_HASHMAP* vartopos;
   int* nvarinconss;
   SCIP_CONS*** varinconss;
   SCIP_VAR** relvars;
   int nrelvars;
   int nrelconss;

   /* general parameters */
   SCIP_Bool tidy;

   /* Graph stuff for hmetis */
   int randomseed;
   SCIP_Real metisubfactor;
   SCIP_Bool metisverbose;
   SCIP_Bool metisuseptyperb;
   SCIP_Bool found;

};

/*
 * Local methods
 */
/* put your local methods here, and declare them static */

static DEC_DECL_INITDETECTOR(initCutpacking)
{
   int i;
   int j;
   int k;
   int nallvars;
   int nconss;
   SCIP_Bool ishandled;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_VAR** allvars;
   SCIP_VAR** relvars;
   int nvars;
   SCIP_HASHMAP* vartopos;
   SCIP_CONS*** varinconss;
   int* nvarinconss;

   DEC_DETECTORDATA* detectordata;
   assert( scip != NULL );

   detectordata = DECdetectorGetData(detector);
   assert( detectordata != NULL );
   assert( strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0 );

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
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->relvars, nallvars) );
   vartopos = detectordata->vartopos;
   relvars = detectordata->relvars;
   j = 0;

   for( i = 0; i < nallvars; ++i )
   {
      if( SCIPisVarRelevant(allvars[i]) )
      {
         relvars[j] = SCIPvarGetProbvar(allvars[i]);
         SCIP_CALL( SCIPhashmapInsert(vartopos, SCIPvarGetProbvar(allvars[i]), (void*) (size_t) j) );
         j++;
      }
   }
   detectordata->nrelvars = j;
   SCIPreallocMemoryArray(scip, &(detectordata->relvars), j);

   /* get number of relevant conss */
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->graphs, nconss+1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(detectordata->graphs[0].conss), nconss) );
   k = 0;
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);
      if( !SCIPconsIsActive(conss[i]) )
      {
         continue;
      }

      nvars = SCIPgetNVarsXXX(scip, conss[i]);
      if( nvars > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, nvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, conss[i], vars, nvars) );
      }
      ishandled = FALSE;

      for( j = 0; (j < nvars) && (ishandled == FALSE); ++j )
      {
         ishandled = SCIPisVarRelevant(vars[j]);
      }

      if( ishandled )
      {
         detectordata->graphs[0].conss[k] = conss[i];
         k++;
      }
      /* SCIPfreeMemoryArrayNull(scip, &vars); */
   }
   detectordata->nrelconss = k;
   detectordata->graphs[0].nconss = k;
   SCIPreallocMemoryArray(scip, &(detectordata->graphs[0].conss), k);

   /* alloc */
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->partition, k) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip),k) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->representatives, SCIPblkmem(scip),k) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->occupied, SCIPblkmem(scip),k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->subscipconss, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->mergedconss, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->nsubscipconss, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->nvarinconss, detectordata->nrelvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varinconss, detectordata->nrelvars) );

   varinconss = detectordata->varinconss;
   nvarinconss = detectordata->nvarinconss;

   for( i = 0; i < k; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->subscipconss[i], k) );
      SCIP_CALL( SCIPhashmapCreate(&detectordata->mergedconss[i], SCIPblkmem(scip),k) );
   }

   /* varinconss */
   for( i = 0; i < detectordata->nrelvars; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varinconss[i], k) );
      detectordata->nvarinconss[i] = 0;
   }

   for( i = 0; i < k; ++i )
   {
      nvars = SCIPgetNVarsXXX(scip, detectordata->graphs[0].conss[i]);
      SCIP_CALL( SCIPallocMemoryArray(scip,&vars,nvars) );
      SCIP_CALL( SCIPgetVarsXXX(scip,detectordata->graphs[0].conss[i], vars, nvars) );
      for( j = 0; j < nvars; ++j )
      {
         if( SCIPisVarRelevant(vars[j]) )
         {
            (varinconss[(long int)SCIPhashmapGetImage(vartopos, SCIPvarGetProbvar(vars[j]))])[nvarinconss[(long int)SCIPhashmapGetImage(vartopos, SCIPvarGetProbvar(vars[j]))]] = detectordata->graphs[0].conss[i];
            ++nvarinconss[(long int)SCIPhashmapGetImage(vartopos, SCIPvarGetProbvar(vars[j]))];
         }
      }
      SCIPfreeMemoryArrayNull(scip, &vars);
   }

   return SCIP_OKAY;
}

/** presolving deinitialization method of presolver (called after presolving has been finished) */
static DEC_DECL_EXITDETECTOR(exitCutpacking)
{
   int i;
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);
   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   /* copy data to decomp structure */
   if( !detectordata->found )
   {
      SCIPfreeMemory(scip, &detectordata);
      return SCIP_OKAY;
   }

   /* free presolver data */

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
   SCIPfreeMemoryArray(scip, &detectordata->partition);
   SCIPfreeMemoryArray(scip, &detectordata->graphs);
   SCIPfreeMemoryArray(scip, &detectordata->varinconss);
   SCIPfreeMemoryArray(scip, &detectordata->nvarinconss);
   SCIPfreeMemoryArray(scip, &detectordata->relvars);
   SCIPhashmapFree(&detectordata->vartopos);
   SCIPhashmapFree(&detectordata->representatives);
   SCIPhashmapFree(&detectordata->occupied);

   for( i = 0; i < detectordata->nrelconss; i++ )
   {
      SCIPhashmapFree(&detectordata->mergedconss[i]);
   }
   SCIPfreeMemoryArray(scip, &detectordata->mergedconss);
   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** builds the graph from the given scip instance */
static SCIP_RETCODE buildGraphStructure(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< detectordata data structure */
   )
{

   int i;
   int j;
   int k;
   int nedges;
   long int cost;
   Graph graph;
   SCIP_CONS*** varinconss;

   graph = detectordata->graphs[0];

   SCIP_CALL( SCIPhashmapCreate(&(graph.constopos), SCIPblkmem(scip),detectordata->nrelconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &graph.adjacencylist, detectordata->nrelconss) );
   for( i = 0; i < detectordata->nrelconss; ++i )
   {
      SCIP_CALL( SCIPhashmapCreate(&graph.adjacencylist[i], SCIPblkmem(scip),detectordata->nrelconss) );
   }

   nedges = 0;
   /* initialize constopos */
   assert( graph.nconss > 0 );

   for( i = 0; i < graph.nconss; ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(graph.constopos, graph.conss[i], (void*) (size_t) i) );
   }

   /* initialize adjacencylist */

   varinconss = detectordata->varinconss;

   for( i = 0; i < detectordata->nrelvars; ++i )
   {
      for( j = 0; j < detectordata->nvarinconss[i]; ++j )
      {
         for( k = j + 1; k < detectordata->nvarinconss[i]; ++k )
         {
            if( SCIPhashmapExists( graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,varinconss[i][j])], varinconss[i][k]) )
            {
               cost = (long int)SCIPhashmapGetImage( graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, varinconss[i][j])], varinconss[i][k]) + 1;
               SCIP_CALL( SCIPhashmapSetImage(graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, varinconss[i][j])], varinconss[i][k], (void*) (size_t) cost) );
               cost = (long int)SCIPhashmapGetImage(graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, varinconss[i][k])], varinconss[i][j]) + 1;
               SCIP_CALL( SCIPhashmapSetImage(graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, varinconss[i][k])], varinconss[i][j], (void*) (size_t) cost) );
            }
            else
            {
               SCIP_CALL( SCIPhashmapInsert(graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, varinconss[i][j])], varinconss[i][k],(void*) (size_t) 1) );
               SCIP_CALL( SCIPhashmapInsert(graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, varinconss[i][k])], varinconss[i][j],(void*) (size_t) 1) );
               nedges++;
            }
         }
      }
   }

   graph.cons1 = NULL;
   graph.cons2 = NULL;
   graph.nedges = nedges;

   detectordata->graphs[0] = graph;
   detectordata->ngraphs = 1;
   SCIP_CALL( SCIPhashmapInsert(detectordata->occupied,(void*)(size_t)1,NULL) );

   return SCIP_OKAY;
}

/** copies hashmap hm1 to hashmap hm2 */
static SCIP_RETCODE copyhashmap(
   SCIP_HASHMAP*         hm1,                /**< pointer to first hashmap */
   SCIP_HASHMAP*         hm2                 /**< pointer to second hashmap */
)
{
   int i;

   assert( hm1 != NULL );
   assert( hm2 != NULL );

   for( i = 0; i < SCIPhashmapGetNLists(hm1); ++i )
   {
      SCIP_HASHMAPLIST* list = SCIPhashmapGetList(hm1, i);
      if( SCIPhashmapListGetNEntries(list) > 0 )
      {
         while( list != NULL )
         {
            SCIP_CALL( SCIPhashmapInsert(hm2, SCIPhashmapListGetOrigin(list),SCIPhashmapListGetImage(list)) );
            list = SCIPhashmapListGetNext(list);
         }
      }
   }

   return SCIP_OKAY;
}

/** returns the next element of the hashmap hm */
static SCIP_HASHMAPLIST* hashmapiteration(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detectordata data structure */
   SCIP_HASHMAP*         hm,                 /**< the hashmap */
   SCIP_HASHMAPLIST*     list                /**< current iteration list */
   )
{
      assert(detectordata->iter < SCIPhashmapGetNLists(hm)+1);
      assert((detectordata->iter == 0)||(list != 0));

      if( list != NULL )
      {
         list = SCIPhashmapListGetNext(list);
         if( list != NULL )
            return list;
      }

      if( list == NULL )
      {
         int j;
         for( j = detectordata->iter; j < SCIPhashmapGetNLists(hm); ++j )
         {
            list = SCIPhashmapGetList(hm,j);
            ++detectordata->iter;
            if( SCIPhashmapListGetNEntries(list) > 0 )
            {
               assert(list != NULL);
               return list;
            }
         }
      }
   return NULL;
}

/** inserts element into hashmap if it doesn't already exist */
static SCIP_RETCODE hashmapinsert(
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


/** builds a new adjacencylist */
static SCIP_RETCODE buildnewadjacencylist(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< presolver data data structure */
   int                   pos,                /**< position */
   int                   nconss,             /**< number of constraints */
   Graph                 graph,              /**< current graph */
   SCIP_HASHMAP*         consslink,          /**< */
   SCIP_HASHMAP*         consslink2          /**< */
   )
{
   int i;
   int cost;
   int nedges;
   SCIP_HASHMAPLIST* list;
   SCIP_CONS* representative;
   Graph newgraph;

   newgraph = detectordata->graphs[pos];
   representative = NULL;

   if( SCIPhashmapGetNEntries(consslink) > 0 )
   {
      int j = 1;
      SCIP_HASHMAP* newadja = detectordata->graphs[pos].adjacencylist[0];

      for( i = 0; i < nconss; ++i )
      {
         assert( consslink != NULL );
         if( !SCIPhashmapExists(consslink, newgraph.conss[i]) )
         {
            SCIP_CALL( copyhashmap(graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, newgraph.conss[i])],newgraph.adjacencylist[j]) );
            SCIP_CALL( SCIPhashmapSetImage(newgraph.constopos, newgraph.conss[i], (void*) (size_t) j) );
            ++j;
         }
         else
         {
            representative = newgraph.conss[i];
            SCIP_CALL( SCIPhashmapRemove(newgraph.constopos, newgraph.conss[i]) );
            detectordata->iter = 0;
            list = NULL;
            do
            {
               list = hashmapiteration(scip, detectordata, graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, newgraph.conss[i])], list);
               if( list == NULL )
                  break;
               if( !SCIPhashmapExists(consslink, SCIPhashmapListGetOrigin(list)) )
               {
                 if( SCIPhashmapExists(newadja, SCIPhashmapListGetOrigin(list)) )
                 {
                    cost = (long int)SCIPhashmapGetImage(newadja, SCIPhashmapListGetOrigin(list));
                    cost += (long int)SCIPhashmapListGetImage(list);
                    SCIP_CALL( SCIPhashmapSetImage(newadja, SCIPhashmapListGetOrigin(list), (void*) (size_t) cost) );
                 }
                 else if( !SCIPhashmapExists(consslink2, SCIPhashmapListGetOrigin(list)) )
                 {
                    SCIP_CALL( SCIPhashmapInsert(newadja, SCIPhashmapListGetOrigin(list), SCIPhashmapListGetImage(list)) );
                 }
              }
            } while (list != NULL);
         }
      }
      SCIP_CALL( SCIPhashmapInsert(newgraph.constopos, representative, (void*) (size_t) 0) );
      nconss = j;

      /* insert representative */
      SCIP_CALL( SCIPhashmapInsert(detectordata->representatives, (void*) (size_t) (detectordata->nrepresentatives + 1), representative ) );
      SCIP_CALL( copyhashmap(consslink,detectordata->mergedconss[detectordata->nrepresentatives]) );
      detectordata->nrepresentatives++;

   }

   if( SCIPhashmapGetNEntries(consslink) == 0 )
   {
      for( i = 0; i < nconss; ++i )
      {
         SCIP_CALL( copyhashmap(graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,newgraph.conss[i])],newgraph.adjacencylist[i]) );
         SCIP_CALL( SCIPhashmapSetImage(newgraph.constopos, newgraph.conss[i], (void*) (size_t) i) );
      }
   }

   nedges = 0;
   /* delete merged conss */
   for( i = 1; i < nconss; ++i )
   {
      cost = 0;
      detectordata->iter = 0;
      list = NULL;
      do
      {
         list = hashmapiteration(scip, detectordata, newgraph.adjacencylist[i], list);
         if( list == NULL )
            break;
         if( SCIPhashmapExists(consslink, SCIPhashmapListGetOrigin(list)) )
         {
            cost += (long int)SCIPhashmapListGetImage(list);
            SCIP_CALL( SCIPhashmapRemove(newgraph.adjacencylist[i], SCIPhashmapListGetOrigin(list)) );
         }
         else
            ++nedges;
      } while (list != NULL);
      if( cost > 0 )
      {
         SCIP_CALL( SCIPhashmapInsert(newgraph.adjacencylist[i], representative, (void*) (size_t) cost) );
         nedges += 2;
      }
   }

   /* arranges conss */
   detectordata->iter = 0;
   list = NULL;
   do
   {
      list = hashmapiteration(scip, detectordata, newgraph.constopos, list);
      if( list == NULL )
         break;
      newgraph.conss[(long int)SCIPhashmapListGetImage(list)] = (SCIP_CONS*)SCIPhashmapListGetOrigin(list);
   } while (list != NULL);



   for( i = nconss; i < nconss + SCIPhashmapGetNEntries(consslink) - 1; ++i )
   {
      SCIPhashmapFree(&newgraph.adjacencylist[i]);
   }

   /* SCIPreallocMemoryArray(scip, &detectordata->graphs[pos].adjacencylist, nconss); */

   newgraph.nconss = nconss;
   newgraph.nedges = nedges / 2;
   assert( 2*newgraph.nedges == nedges );

   detectordata->graphs[pos] = newgraph;

   return SCIP_OKAY;
}

/** frees graph at position pos */
static SCIP_RETCODE FreeGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detectordata data structure */
   int                   pos,                /**< */
   int                   nconss              /**< */
   )
{
   int i;

   for( i = 0; i < nconss; ++i )
   {
      SCIPhashmapFree(&detectordata->graphs[pos].adjacencylist[i]);
   }
   SCIPhashmapFree(&detectordata->graphs[pos].constopos);
   SCIPfreeMemoryArray(scip, &detectordata->graphs[pos].adjacencylist);
   SCIPfreeMemoryArray(scip, &detectordata->graphs[pos].conss);

   return SCIP_OKAY;
}

/** allocates memory at position pos to facilitate saving a graph with nconss vertices */
static SCIP_RETCODE AllocateMemoryGraph(
   SCIP*                scip,             /**< SCIP data structure */
   DEC_DETECTORDATA*    detectordata,     /**< detectordata data structure */
   int                  pos,              /**< */
   int                  nconss            /**< */
   )
{
   int i;

   SCIP_CALL( SCIPhashmapCreate(&detectordata->graphs[pos].constopos, SCIPblkmem(scip),nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->graphs[pos].adjacencylist, nconss) );

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(detectordata->graphs[pos].constopos, detectordata->graphs[pos].conss[i], NULL) );
      SCIP_CALL( SCIPhashmapCreate(&detectordata->graphs[pos].adjacencylist[i], SCIPblkmem(scip),nconss) );
   }

   return SCIP_OKAY;
}

/** assigns the right linking constraint to the graph at position pos */
static SCIP_RETCODE SetLinkingCons(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< presolver data data structure */
   int                   cas,                /**< */
   int                   cas2,               /**< */
   int                   pos,                /**< */
   SCIP_CONS*            cons1,              /**< */
   SCIP_CONS*            cons2               /**< */
   )
{
   switch( cas2 )
   {
   case 0:
      if( cas )
      {
         detectordata->graphs[pos].cons1 = cons2;
         detectordata->graphs[pos].cons2 = cons2;
      }
      else
      {
         detectordata->graphs[pos].cons1 = cons1;
         detectordata->graphs[pos].cons2 = cons1;
      }
      break;
   case 1:
      if( cas )
      {
         detectordata->graphs[pos].cons1 = detectordata->graphs[pos].conss[0];
         detectordata->graphs[pos].cons2 = cons2;
      }
      else
      {
         detectordata->graphs[pos].cons1 = cons1;
         detectordata->graphs[pos].cons2 = detectordata->graphs[pos].conss[0];
      }
      break;
   case 2:
      if( !cas )
      {
         detectordata->graphs[pos].cons1 = detectordata->graphs[pos].conss[0];
         detectordata->graphs[pos].cons2 = cons1;
      }
      else
      {
         detectordata->graphs[pos].cons1 = cons2;
         detectordata->graphs[pos].cons2 = detectordata->graphs[pos].conss[0];
      }
      break;
   default:
      break;
   }

   return SCIP_OKAY;
}

/** sets the startblock */
static SCIP_RETCODE SetStartblock(
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
static SCIP_RETCODE CopyConss(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detectordata data structure */
   int                   pos,                /**< position of constraint */
   int                   nconss              /**< number of constraints */
   )
{
   int i;
   for( i = 0; i < nconss; ++i )
   {
      detectordata->subscipconss[detectordata->nblocks][i] = detectordata->graphs[pos].conss[i];
   }
   detectordata->nsubscipconss[detectordata->nblocks] = nconss;
   detectordata->nblocks++;

   return SCIP_OKAY;
}


/** builds the new graphs which result from the last found cut */
static SCIP_RETCODE buildnewgraphs(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< detectordata data structure */
   )
{
   int i;
   int pos1;
   int pos2;
   int cas;
   int stop1;
   int stop2;

   int* partition;
   Graph graph;

   int nconss1;
   SCIP_HASHMAP* consslink1;
   int nconss2;
   SCIP_HASHMAP* consslink2;

   cas = -1;
   nconss1 = 0;
   nconss2 = 0;
   stop1 = 0;
   stop2 = 0;

   /* build partitions */
   partition = detectordata->partition;
   graph = detectordata->graphs[detectordata->position];

   pos1 = -1;
   pos2 = -1;
   for( i = 0; i < detectordata->nrelconss + 1; ++i )
   {
      if( (!SCIPhashmapExists(detectordata->occupied, (void*)(size_t)(i + 1))) )
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
   assert( i != detectordata->nrelconss+1 );
   assert( (pos1 != -1) && (pos2 != -1) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->graphs[pos1].conss, graph.nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->graphs[pos2].conss, graph.nconss) );
   SCIP_CALL( SCIPhashmapRemove(detectordata->occupied, (void*)(size_t)(detectordata->position + 1)) );

   for( i = 0; i < graph.nconss; ++i )
   {
      assert( (-1 < partition[i]) && (partition[i] < 2) );
      if( partition[i] == 0 )
      {
         detectordata->graphs[pos1].conss[nconss1] = graph.conss[i];
         nconss1++;
      }
      else
      {
         detectordata->graphs[pos2].conss[nconss2] = graph.conss[i];
         nconss2++;
      }
   }

   SCIP_CALL( AllocateMemoryGraph(scip, detectordata, pos1, nconss1) );
   SCIP_CALL( AllocateMemoryGraph(scip, detectordata, pos2, nconss2) );

   assert( nconss1 + nconss2 == graph.nconss );
   assert( (nconss1 != 0) && (nconss2 != 0) );

   /* get linking conss */
   SCIP_CALL( SCIPhashmapCreate(&consslink1, SCIPblkmem(scip), nconss1) );
   SCIP_CALL( SCIPhashmapCreate(&consslink2, SCIPblkmem(scip), nconss2) );
   assert( SCIPhashmapIsEmpty(consslink1) );
   assert( SCIPhashmapIsEmpty(consslink2) );

   for( i = 0; i < nconss1; ++i )
   {
      SCIP_HASHMAPLIST* list = NULL;
      SCIP_HASHMAP* adja = graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, detectordata->graphs[pos1].conss[i])];
      detectordata->iter = 0;

      do
      {
         list = hashmapiteration(scip, detectordata, adja, list);
         if( list == NULL )
            break;
         if( SCIPhashmapExists(detectordata->graphs[pos2].constopos, SCIPhashmapListGetOrigin(list)) )
         {
            SCIP_CALL( hashmapinsert(consslink1, detectordata->graphs[pos1].conss[i], NULL) );
            SCIP_CALL( hashmapinsert(consslink2, SCIPhashmapListGetOrigin(list), NULL) );
         }
      }
      while (list != NULL);
   }

   if( SCIPhashmapGetNEntries(consslink1) == nconss1 )
      stop1 = 1;
   if( SCIPhashmapGetNEntries(consslink2) == nconss2 )
      stop2 = 1;

   /* test whether the cut is feasible */

   if( (graph.cons1 != NULL) && (graph.cons2 != NULL) )
   {
      if( (SCIPhashmapExists(detectordata->graphs[pos1].constopos, graph.cons1 )
         && SCIPhashmapExists(detectordata->graphs[pos1].constopos, graph.cons2))
         || (SCIPhashmapExists(detectordata->graphs[pos2].constopos, graph.cons1)
         && SCIPhashmapExists(detectordata->graphs[pos2].constopos, graph.cons2)) )
      {
         SCIP_CALL( CopyConss(scip, detectordata, detectordata->position, graph.nconss) );
         detectordata->ngraphs--;
         SCIP_CALL( FreeGraph(scip, detectordata, pos1, nconss1) );
         SCIP_CALL( FreeGraph(scip, detectordata, pos2, nconss2) );
         SCIPhashmapFree(&consslink1);
         SCIPhashmapFree(&consslink2);
         return SCIP_OKAY;
      }

      if( SCIPhashmapExists(detectordata->graphs[pos1].constopos, graph.cons1) )
      {
        cas = 0;
        stop1 = (int)SCIPhashmapExists(consslink1, graph.cons1);
        stop2 = (int)SCIPhashmapExists(consslink2, graph.cons2);
     }
     else
     {
        cas = 1;
        stop1 = (int)SCIPhashmapExists(consslink1, graph.cons2);
        stop2 = (int)SCIPhashmapExists(consslink2, graph.cons1);
     }
   }

   /* test right or left*/

   if( (graph.cons1 != NULL) && (graph.cons2 == NULL) )
   {
      if( SCIPhashmapExists(detectordata->graphs[pos1].constopos, graph.cons1) )
      {
         cas = 0;
         stop1 = (int) SCIPhashmapExists(consslink1, graph.cons1);
      }
      else
      {
         cas = 1;
         stop2 = (int)SCIPhashmapExists(consslink2, graph.cons1);
      }
   }
   else if( (graph.cons1 == NULL) && (graph.cons2 != NULL) )
   {
      if( SCIPhashmapExists(detectordata->graphs[pos2].constopos, graph.cons2) )
      {
         cas = 0;
         stop2 = SCIPhashmapExists(consslink2, graph.cons2);
      }
      else
      {
         cas = 1;
         stop1 = SCIPhashmapExists(consslink1, graph.cons2);
      }
   }
   else if( (graph.cons1 == NULL) && (graph.cons2 == NULL) )
      cas = 1;


   if( (nconss1 > 1) && !stop1 )
   {
      SCIP_CALL( buildnewadjacencylist(scip, detectordata, pos1, nconss1, graph, consslink1, consslink2) );
      SCIP_CALL( SetLinkingCons(scip, detectordata, cas, 1, pos1, graph.cons1, graph.cons2) );
   }
   else if( stop1 )
   {
      SCIP_CALL( SetLinkingCons(scip, detectordata, cas, 0, pos1, graph.cons1, graph.cons2) );
   }

   if( (nconss2 > 1) && !stop2 )
   {
      SCIP_CALL( buildnewadjacencylist(scip, detectordata, pos2, nconss2, graph, consslink2, consslink1) );
      SCIP_CALL( SetLinkingCons(scip, detectordata, cas, 2, pos2, graph.cons2, graph.cons1) );
   }
   else if( stop2 )
   {
      SCIP_CALL( SetLinkingCons(scip, detectordata, cas, 0, pos2, graph.cons2, graph.cons1) );
   }

   if( ((nconss1 < 2) && (nconss2 < 2)) || (stop1 && stop2 ) )
   {
      SCIP_CALL( CopyConss(scip, detectordata, pos1, nconss1) );
      SCIP_CALL( SetStartblock(scip, detectordata,detectordata->graphs[pos1].cons1) );
      SCIP_CALL( CopyConss(scip, detectordata, pos2, nconss2) );
      detectordata->ngraphs--;
      SCIP_CALL( SetStartblock(scip, detectordata,detectordata->graphs[pos2].cons1) );
      SCIP_CALL( FreeGraph(scip, detectordata, pos1, nconss1) );
      SCIP_CALL( FreeGraph(scip, detectordata, pos2, nconss2) );
   }
   else if( (nconss1 < 2) || (stop1 && (stop2 == 0)) )
   {
      SCIP_CALL( CopyConss(scip, detectordata, pos1, nconss1) );
      SCIP_CALL( SetStartblock(scip, detectordata,detectordata->graphs[pos1].cons1) );
      SCIP_CALL( SCIPhashmapInsert(detectordata->occupied,(void*) (size_t) (pos2 + 1), NULL) );
      SCIP_CALL( FreeGraph(scip, detectordata, pos1, nconss1) );
   }
   else if( (nconss2 < 2) || ((stop1 == 0) && stop2) )
   {
      SCIP_CALL( CopyConss(scip, detectordata, pos2, nconss2) );
      SCIP_CALL( SetStartblock(scip, detectordata,detectordata->graphs[pos2].cons1) );
      SCIP_CALL( SCIPhashmapInsert(detectordata->occupied,(void*) (size_t) (pos1 + 1), NULL) );
      SCIP_CALL( FreeGraph(scip, detectordata, pos2, nconss2) );
   }
   else
   {
      SCIP_CALL( SCIPhashmapInsert(detectordata->occupied,(void*) (size_t) (pos1 + 1), NULL) );
      SCIP_CALL( SCIPhashmapInsert(detectordata->occupied,(void*) (size_t) (pos2 + 1), NULL) );
      detectordata->ngraphs++;
   }

   SCIP_CALL( FreeGraph(scip, detectordata, detectordata->position, graph.nconss) );
   SCIPhashmapFree(&consslink1);
   SCIPhashmapFree(&consslink2);

   return SCIP_OKAY;
}

/** adds the merged constraints to the right blocks */
static
SCIP_RETCODE getmergedconss(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< detectordata data structure */
   )
{
   int i;
   int j;
   SCIP_HASHMAP** mergedconss;
   SCIP_HASHMAP* representatives;
   SCIP_CONS*** subscipconss;
   SCIP_HASHMAP* constoblock;

   mergedconss = detectordata->mergedconss;
   representatives = detectordata->representatives;
   subscipconss = detectordata->subscipconss;
   constoblock = detectordata->constoblock;

   /* constoblock */
   for( i = 0; i < detectordata->nblocks; ++i )
   {
      for( j = 0; j < detectordata->nsubscipconss[i]; ++j )
      {
         SCIP_CALL( hashmapinsert(constoblock, subscipconss[i][j], (void *) (size_t) (i+1)) );
      }
   }

   for( i = detectordata->nrepresentatives; i > 0; --i )
   {
      int block = (long int)SCIPhashmapGetImage(constoblock, SCIPhashmapGetImage(representatives, (void*)(size_t)i));
      SCIP_HASHMAPLIST* list = NULL;
      detectordata->iter = 0;

      do
      {
         list = hashmapiteration(scip, detectordata, mergedconss[i - 1], list);
         if( list == NULL )
            break;
         if( (SCIP_CONS*)SCIPhashmapGetImage(representatives, (void*)(size_t)i) != (SCIP_CONS*)SCIPhashmapListGetOrigin(list) )
         {
            subscipconss[block - 1][detectordata->nsubscipconss[block - 1]] = (SCIP_CONS*)SCIPhashmapListGetOrigin(list);
            detectordata->nsubscipconss[block - 1]++;
            SCIP_CALL( hashmapinsert(constoblock,(SCIP_CONS*)SCIPhashmapListGetOrigin(list),(void*) (size_t) block) );
         }
      } while (list != NULL);
   }
   assert( SCIPhashmapGetNEntries(constoblock) == detectordata->nrelconss );

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
SCIP_RETCODE GetConsindex(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< presolver data structure */
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
   int linkingblock;
   int counter;

   int nvars;
   SCIP_VAR** vars;

   SCIP_CONS*** subscipconss;
   SCIP_CONS*** subscipconss2;
   int* nsubscipconss2;
   SCIP_VAR*** stairlinkingvars;
   int* nstairlinkingvars;
   SCIP_VAR** linkingvars;


   newblock = 0;
   oldblock = 0;
   linkingblock = 0;

   subscipconss = detectordata->subscipconss;

   SCIP_CALL( SCIPallocMemoryArray(scip, &stairlinkingvars, detectordata->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &linkingvars, detectordata->nrelvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nstairlinkingvars, detectordata->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &subscipconss2, detectordata->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nsubscipconss2, detectordata->nblocks) );

   for( i = 0; i < detectordata->nblocks; ++i )
   {
     SCIP_CALL( SCIPallocMemoryArray(scip, &stairlinkingvars[i], detectordata->nrelvars) );
     nstairlinkingvars[i] = 0;
     nsubscipconss2[i] = 0;
   }

   counter = detectordata->nblocks;
   block = detectordata->startblock;

   while( counter > 0 )
   {
      for( i = 0; i < detectordata->nsubscipconss[block - 1]; ++i )
      {
         nvars = SCIPgetNVarsXXX(scip, detectordata->subscipconss[block - 1][i]);
         SCIP_CALL( SCIPallocMemoryArray(scip, &vars, nvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, detectordata->subscipconss[block - 1][i], vars, nvars) );

         for( j = 0; j < nvars; ++j )
         {
            while( (!SCIPisVarRelevant(vars[j])) && (j < nvars) )
            {
               ++j;
            }
            if( j >= nvars )
               break;
            no = (long int)SCIPhashmapGetImage(detectordata->vartopos, SCIPvarGetProbvar(vars[j]));
            for( k = 0; (k < detectordata->nvarinconss[no]); ++k )
            {
               actblock = (long int)SCIPhashmapGetImage(detectordata->constoblock, detectordata->varinconss[no][k]);
               if( actblock != block )
               {
                  if( actblock != oldblock )
                  {
                     newblock = actblock;
                  }
               }
            }
         }
         SCIPfreeMemoryArrayNull(scip, &vars);
      }
      assert( (newblock != block)||(counter<2) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &subscipconss2[detectordata->nblocks-counter], detectordata->nsubscipconss[block - 1]) );
      for( j=0; j < detectordata->nsubscipconss[block - 1]; ++j )
      {
         subscipconss2[detectordata->nblocks-counter][nsubscipconss2[detectordata->nblocks-counter]] = subscipconss[block - 1][j];
         ++nsubscipconss2[detectordata->nblocks-counter];
      }
      oldblock = block;
      block = newblock;
      assert( 0 < block );
      assert( block <= detectordata->nblocks );

      detectordata->iter = 0;

      ++linkingblock;
      --counter;
   }

   for( i = 0; i < detectordata->nblocks; ++i )
   {
      for( j = 0; j < nsubscipconss2[i]; ++j )
      {
         SCIP_CALL( SCIPhashmapSetImage(detectordata->constoblock, subscipconss2[i][j],(void*)(size_t)(i+1)) );
      }
   }

   if( !detectordata->fixedblocks )
   {

      SCIP_CALL( DECfilloutDecompFromConstoblock(scip, decdecomp, detectordata->constoblock, detectordata->nblocks, TRUE) );

      for( i = 0; i < detectordata->nblocks; ++i )
      {
         SCIPfreeMemoryArray(scip, &subscipconss2[i]);
      }
      SCIPfreeMemoryArray(scip, &subscipconss2);
      SCIPfreeMemoryArray(scip, &nsubscipconss2);
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
static SCIP_RETCODE GetLinkingVars(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< presolver data structure */
   DEC_DECOMP*           decdecomp           /**< decdecomp pointer */
   )
{
   int i;
   int j;
   int newblock;

   SCIP_VAR** linkingvars;
   int nlinkingvars;
   SCIP_CONS*** varinconss;
   int* nvarinconss;
   SCIP_VAR*** subscipvars;
   int* nsubscipvars;
   SCIP_HASHMAP* vartoblock;

   varinconss = detectordata->varinconss;
   nvarinconss = detectordata->nvarinconss;


   SCIP_CALL( SCIPallocMemoryArray(scip, &nsubscipvars, detectordata->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &subscipvars, detectordata->nblocks) );
   SCIP_CALL( SCIPhashmapCreate(&vartoblock, SCIPblkmem(scip),detectordata->nrelvars) );
   for( i = 0; i < detectordata->nblocks; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &subscipvars[i], detectordata->nrelvars) );
   }
   for( i = 0; i < detectordata->nblocks; ++i )
   {
      nsubscipvars[i] = 0;
   }

   /* bulid linkingvars, subscipvars and vartoblock */

   nlinkingvars = 0;
   SCIP_CALL( SCIPallocMemoryArray(scip, &linkingvars, detectordata->nrelvars) );

   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      int oldblock;
      SCIP_Bool stop;

      while( !SCIPisVarRelevant(SCIPgetVars(scip)[i]) )
      {
         ++i;
      }
      if( i >= SCIPgetNVars(scip) )
         break;
      oldblock = (int) (size_t) SCIPhashmapGetImage(detectordata->constoblock, varinconss[(long int) SCIPhashmapGetImage(detectordata->vartopos, SCIPgetVars(scip)[i])][0]);
      stop = FALSE;
      for( j = 1; !stop && (j < nvarinconss[(size_t) SCIPhashmapGetImage(detectordata->vartopos, SCIPgetVars(scip)[i])]); ++j )
      {
         newblock = (int) (size_t) SCIPhashmapGetImage(detectordata->constoblock, varinconss[(long int) SCIPhashmapGetImage(detectordata->vartopos, SCIPgetVars(scip)[i])][j]);
         if( newblock != oldblock )
         {
            stop = TRUE;
         }
      }

      if( stop )
      {
         linkingvars[nlinkingvars] = SCIPgetVars(scip)[i];
         ++nlinkingvars;
         SCIP_CALL( SCIPhashmapInsert(vartoblock,SCIPgetVars(scip)[i],(void*) (size_t) (detectordata->nblocks+1)) );
      }
      else
      {
         SCIP_CALL( SCIPhashmapInsert(vartoblock,SCIPgetVars(scip)[i],(void*) (size_t) oldblock) );
         subscipvars[oldblock-1][nsubscipvars[oldblock-1]] = SCIPgetVars(scip)[i];
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
   SCIP_CALL( DECdecompSetLinkingvars(scip, decdecomp, linkingvars, nlinkingvars) );
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
static SCIP_RETCODE FixedBlocks(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< presolver data structure */
   DEC_DECOMP*           decdecomp           /**< decdecomp pointer */
   )
{
   int i;
   int block;
   int blocksize;
   int indexcons;
   int newblock;
   SCIP_HASHMAP* consindex;

   SCIP_CALL( SCIPhashmapCreate(&consindex, SCIPblkmem(scip),detectordata->nrelconss) );

   /* get nblocks */
   blocksize = detectordata->blocksize;
   block = detectordata->nrelconss/blocksize + 1;
   if( detectordata->nrelconss/blocksize * blocksize == detectordata->nrelconss )
      block = detectordata->nrelconss/blocksize;
   assert( block > 0 );

   detectordata->nblocks = block;

   SCIPdebugMessage("nblocks %d; nconss %d; blocksize %d; \n", block, detectordata->nrelconss, blocksize);

   SCIP_CALL( SCIPhashmapRemoveAll(detectordata->constoblock) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &detectordata->nsubscipconss, block) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &detectordata->subscipconss, block) );
   for( i = 0; i < block; ++i )
   {
      SCIPreallocMemoryArray(scip, &detectordata->subscipconss[i], blocksize);
      detectordata->nsubscipconss[i] = 0;
   }

   /* build subscipcons and constoblock */
   for( i = 0; i < SCIPgetNConss(scip); ++i )
   {
      if( SCIPhashmapExists(consindex,SCIPgetConss(scip)[i]) )
      {
         assert( SCIPhashmapExists(consindex, SCIPgetConss(scip)[i]) );

         indexcons = (long int) SCIPhashmapGetImage(consindex, SCIPgetConss(scip)[i]);
         newblock = indexcons/blocksize + 1;
         if( indexcons/blocksize * blocksize == indexcons )
            newblock = indexcons/blocksize;

         assert( newblock > 0 );
         assert( newblock < block+1 );

         SCIP_CALL( SCIPhashmapInsert(detectordata->constoblock, SCIPgetConss(scip)[i], (void*) (size_t) newblock ) );

         assert(detectordata->subscipconss[newblock-1] != NULL);
         assert(detectordata->nsubscipconss != NULL);

         detectordata->subscipconss[newblock-1][detectordata->nsubscipconss[newblock-1]] = SCIPgetConss(scip)[i];
         ++detectordata->nsubscipconss[newblock-1];
      }
   }

   assert( SCIPhashmapGetNEntries(detectordata->constoblock) == detectordata->nrelconss);

   SCIP_CALL( DECfilloutDecompFromConstoblock(scip, decdecomp, detectordata->constoblock, detectordata->nblocks, TRUE) );
   return SCIP_OKAY;
}

/** will find a minimum cut via the Stoer-Wagner algorithm */
static SCIP_RETCODE StoerWagner(
   SCIP*                 scip,               /**< SCIP data struture */
   DEC_DETECTORDATA*     detectordata        /**< presolver data data structure */
   )
{
   int i;
   int j;
   int tight;
   double value_cut;
   SCIP_HASHMAP* tightness;
   SCIP_HASHMAP* repres_conss;
   SCIP_HASHMAP* constopos;
   SCIP_HASHMAP** adja;
   SCIP_CONS** mincut;
   int nmincut;
   int nrepres_conss;
   SCIP_CONS*** merged_conss;
   int* nmerged_conss;
   Graph graph;
   SCIP_CONS* s;
   SCIP_CONS* t;
   SCIP_CONS* cut;
   SCIP_CONS* last;
   SCIP_CONS* next_to_last;
   SCIP_CONS* represent_t;

   graph = detectordata->graphs[detectordata->position];
   nrepres_conss = 1;

   cut = NULL;

   SCIP_CALL( SCIPhashmapCreate(&tightness, SCIPblkmem(scip), graph.nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mincut, graph.nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nmerged_conss, graph.nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &merged_conss, graph.nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &adja, graph.nconss) );
   for( i = 0; i < graph.nconss; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(merged_conss[i]), graph.nconss) );
      SCIP_CALL( SCIPhashmapCreate(&(adja[i]), SCIPblkmem(scip), graph.nconss) );
   }
   SCIP_CALL( SCIPhashmapCreate(&constopos, SCIPblkmem(scip), graph.nconss) );
   SCIP_CALL( SCIPhashmapCreate(&repres_conss, SCIPblkmem(scip), graph.nconss) );

   /* copy constopos */
   SCIP_CALL( copyhashmap(graph.constopos, constopos) );

   /* copy adja */
   for( i = 0; i < graph.nconss; ++i )
   {
      SCIP_CALL( copyhashmap(graph.adjacencylist[i], adja[i]) );
   }

   SCIPdebugMessage("stoerwagner \n");

   t = NULL;
   value_cut = INFINITY;
   if( graph.cons1 != NULL )
   {
      s = graph.cons1;
      if( graph.cons2 != NULL )
      {
         t = graph.cons2;
      }
   }
   else
   {
      if( graph.cons2 != NULL )
      {
         s = graph.cons2;
      }
      else
      {
         s = graph.conss[0];
      }
   }

   assert(s != NULL);
   represent_t = t;
   while( SCIPhashmapGetNEntries(constopos) > 1 )
   {
      SCIP_HASHMAPLIST* list = NULL;
      double value_act_cut;

      SCIP_CALL( SCIPhashmapRemoveAll(tightness) );
      detectordata->iter = 0;

      do
      {
         list = hashmapiteration(scip, detectordata, constopos, list);
         if( list == NULL )
            break;
         if( SCIPhashmapListGetOrigin(list) != s )
         {
           SCIP_CALL( SCIPhashmapInsert(tightness, SCIPhashmapListGetOrigin(list), (void*)0) );
         }
      } while (list != NULL);

      assert( SCIPhashmapGetNEntries(tightness) + 1 == SCIPhashmapGetNEntries(constopos) );
      last = s;
      next_to_last = s;

      while( (SCIPhashmapGetNEntries(tightness) > 0) )
      {
         next_to_last = last;
         /* actualize tightness */
         detectordata->iter = 0;
         list = NULL;
         do
         {
            list = hashmapiteration(scip, detectordata, adja[(long int)SCIPhashmapGetImage(constopos, last)], list);
            if( list == NULL )
               break;
            assert( SCIPhashmapExists(constopos, SCIPhashmapListGetOrigin(list)) );
            if( SCIPhashmapExists(tightness, SCIPhashmapListGetOrigin(list)) )
            {
               j = ((long int)SCIPhashmapGetImage(tightness, SCIPhashmapListGetOrigin(list))) + ((long int)SCIPhashmapListGetImage(list));
               SCIP_CALL( SCIPhashmapSetImage(tightness,SCIPhashmapListGetOrigin(list),(void*) (size_t) j ) );
            }
         } while (list != NULL);
         /* choose the most tight */
         tight = 0;
         detectordata->iter = 0;
         list = NULL;
         do
         {
            list = hashmapiteration(scip, detectordata, tightness, list);
            if( list == NULL )
               break;
            if( (long int)SCIPhashmapListGetImage(list) >= tight )
            {
               last = (SCIP_CONS*)SCIPhashmapListGetOrigin(list);
               tight = (long int)SCIPhashmapListGetImage(list);
            }
         } while (list != NULL);
         assert(next_to_last != last);
         assert(SCIPhashmapExists(tightness, last));
         SCIP_CALL( SCIPhashmapRemove(tightness, last) );
      }
      /* calculate the value of the current cut */
      value_act_cut = 0;
      detectordata->iter = 0;
      list = NULL;
      do
      {
         list = hashmapiteration(scip, detectordata, adja[(long int)SCIPhashmapGetImage(constopos, last)], list);
         if( list == NULL )
            break;
         assert( SCIPhashmapListGetOrigin(list) != last );
         if( (SCIPhashmapExists(constopos, SCIPhashmapListGetOrigin(list))) )
         {
            value_act_cut += (long int)SCIPhashmapListGetImage(list);
         }
      } while (list != NULL);

      if( (value_act_cut < value_cut) )
      {
         if( (t != NULL) && (represent_t == last) )
         {
            /* test if act_cut ist a s-t-cut */
            represent_t = next_to_last;
            value_cut = value_act_cut;
            cut = last;
            assert(value_cut != 0);
         }
         else if( t == NULL )
         {
            value_cut = value_act_cut;
            cut = last;
            assert(value_cut != 0);
         }
      }

      /* merging */
      if( SCIPhashmapGetNEntries(repres_conss) == 0 )
      {
         SCIP_CALL( SCIPhashmapInsert(repres_conss, next_to_last, (void*) (size_t) nrepres_conss) );
         nrepres_conss++; /* starts with 1 */
         nmerged_conss[nrepres_conss - 2] = 1;
         merged_conss[nrepres_conss - 2][0] = last;
      }
      else if( !SCIPhashmapExists(repres_conss, next_to_last) )
      {
         SCIP_CALL( SCIPhashmapInsert(repres_conss, next_to_last, (void*) (size_t) nrepres_conss) );
         nrepres_conss++; /* starts with 1 */
         nmerged_conss[nrepres_conss - 2] = 1;
         merged_conss[nrepres_conss - 2][0] = last;
      }
      else
      {
         merged_conss[((long int)SCIPhashmapGetImage(repres_conss, next_to_last)) - 1][nmerged_conss[((long int)SCIPhashmapGetImage(
            repres_conss, next_to_last)) - 1]] = last;
         nmerged_conss[((long int)SCIPhashmapGetImage(repres_conss, next_to_last)) - 1]++;
      }

      /* in adja: connect last and next_to_last */
      detectordata->iter = 0;
      list = NULL;
      do
      {
         list = hashmapiteration(scip, detectordata, adja[(long int)SCIPhashmapGetImage(constopos, last)], list);
         if( list == NULL )
            break;
         if( SCIPhashmapExists(adja[(long int)SCIPhashmapGetImage(constopos, next_to_last)], SCIPhashmapListGetOrigin(list)) )
         {
            j = (long int)SCIPhashmapGetImage(adja[(long int)SCIPhashmapGetImage(constopos, next_to_last)],
               SCIPhashmapListGetOrigin(list)) + (long int)SCIPhashmapListGetImage(list);
            SCIP_CALL( SCIPhashmapSetImage(adja[(long int) SCIPhashmapGetImage(constopos, next_to_last)], SCIPhashmapListGetOrigin(list),(void*) (size_t) j) );
            SCIP_CALL( SCIPhashmapSetImage(adja[(long int) SCIPhashmapGetImage(constopos, SCIPhashmapListGetOrigin(list))], next_to_last, (void*)(size_t) j) );
         }
         else
         {
            if( SCIPhashmapListGetOrigin(list) != next_to_last )
            {
               SCIP_CALL( SCIPhashmapInsert(adja[(long int) SCIPhashmapGetImage(constopos, next_to_last)],SCIPhashmapListGetOrigin(list), SCIPhashmapListGetImage(list)) );
               SCIP_CALL( SCIPhashmapInsert(adja[(long int) SCIPhashmapGetImage(constopos, SCIPhashmapListGetOrigin(list))],next_to_last, SCIPhashmapListGetImage(list)) );
            }
         }
      } while (list != NULL);

      SCIP_CALL( SCIPhashmapRemove(constopos, last) );

      /* delete last in adja */
      for( i = 0; i < graph.nconss; ++i )
      {
         if( SCIPhashmapExists(constopos, graph.conss[i]) )
         {
            if( SCIPhashmapExists(adja[i], last) )
               SCIP_CALL( SCIPhashmapRemove(adja[i], last) );
         }
      }
   }

   /* complete cut */
   nmincut = 1;
   mincut[0] = cut;
   for( i = 0; i < nmincut; ++i )
   {
      if( SCIPhashmapExists(repres_conss, mincut[i]) )
      {
         for( j = 0; j < nmerged_conss[((long int)SCIPhashmapGetImage(repres_conss, mincut[i])) - 1]; ++j )
         {
            mincut[nmincut] = merged_conss[((long int)SCIPhashmapGetImage(repres_conss, mincut[i])) - 1][j];
            nmincut++;
         }
      }
   }

   /* create partition */
   for( i = 0; i < graph.nconss; ++i )
   {
      detectordata->partition[i] = 0;
   }
   for( i = 0; i < nmincut; ++i )
   {
      detectordata->partition[(long int)SCIPhashmapGetImage(graph.constopos, mincut[i])] = 1;
   }

   SCIPhashmapFree(&tightness);
   SCIPhashmapFree(&constopos);
   for( i = 0; i < graph.nconss; ++i )
   {
      SCIPhashmapFree(&(adja[i]));
      SCIPfreeMemoryArray(scip, &(merged_conss[i]));
   }
   SCIPfreeMemoryNull(scip, &adja);
   SCIPhashmapFree(&repres_conss);
   SCIPfreeMemoryArray(scip, &mincut);
   SCIPfreeMemoryArray(scip, &nmerged_conss);
   SCIPfreeMemoryArray(scip, &merged_conss);

   return SCIP_OKAY;
}


/** will call hmetis via a system call */
static SCIP_RETCODE callMetis(
   SCIP*                 scip,               /**< SCIP data struture */
   DEC_DETECTORDATA*     detectordata,       /**< presolver data data structure */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   char metiscall[SCIP_MAXSTRLEN];
   char metisout[SCIP_MAXSTRLEN];
   char line[SCIP_MAXSTRLEN];
   char tempfile[SCIP_MAXSTRLEN];

   int i;
   int status;
   int nvertices;
   int nedges;
   void* entry;
   void* cost;
   int* partition;
   SCIP_HASHMAP** adja;
   SCIP_HASHMAP* constopos;
   SCIP_FILE *zfile;
   FILE* file;
   int temp_filedes = -1;

   assert(scip != NULL);
   assert(detectordata != NULL);

   *result = SCIP_DIDNOTRUN;

   adja = detectordata->graphs[detectordata->position].adjacencylist;
   constopos = detectordata->graphs[detectordata->position].constopos;
   nvertices = detectordata->graphs[detectordata->position].nconss;
   nedges = detectordata->graphs[detectordata->position].nedges;

   assert( adja != NULL );
   assert( constopos != NULL );

   SCIPsnprintf(tempfile, SCIP_MAXSTRLEN, "gcg-metis-XXXXXX");
   if( (temp_filedes = mkstemp(tempfile)) < 0 )
   {
      SCIPerrorMessage("Error creating temporary file: %s\n", strerror(errno));
      return SCIP_FILECREATEERROR;
   }

   SCIPdebugMessage("Temporary filename: %s\n", tempfile);

   file = fdopen(temp_filedes, "w");
   if( file == NULL )
   {
      SCIPerrorMessage("Could not open temporary metis file!\n");
      return SCIP_FILECREATEERROR;
   }

   SCIPinfoMessage(scip, file, "%d %d 1\n", nedges, nvertices);

   for( i = 0; i < nvertices - 1; ++i )
   {
      SCIP_HASHMAPLIST* list = NULL;
      assert(adja[i] != NULL);
      detectordata->iter = 0;

      do
      {
         list = hashmapiteration(scip, detectordata, adja[i], list);
         if( list == NULL )
            break;
         entry = SCIPhashmapGetImage(constopos, SCIPhashmapListGetOrigin(list));
         if( (long int)entry > i )
         {
            cost = SCIPhashmapListGetImage(list);
            SCIPinfoMessage(scip, file, "%d ", (long int)cost);
            SCIPinfoMessage(scip, file, "%d ", (long int)entry + 1);
            SCIPinfoMessage(scip, file, "%d \n", i + 1);
         }
      } while (list != NULL);
   }
   status = fclose(file);

   if( status == -1 )
   {
      SCIPerrorMessage("Could not close '%s'\n", tempfile);
      return SCIP_WRITEERROR;
   }

   SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"", tempfile, 2,
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

   SCIPsnprintf(metisout, SCIP_MAXSTRLEN, "%s.part.%d", tempfile, 2);

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


static
DEC_DECL_DETECTSTRUCTURE(detectAndBuildCutpacking)
{
   int i;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(decdecomps != NULL);
   assert(ndecdecomps != NULL);

   *ndecdecomps = 1;
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting cutpacking structure: ");

   SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, 1) );

   /* build the hypergraph structure from the original problem */
   SCIP_CALL( buildGraphStructure(scip, detectordata) );

   SCIP_CALL( DECdecompCreate(scip, &(*decdecomps)[0]) );

   /* get the partitions for the new variables from metis */
   while( detectordata->ngraphs > 0 )
   {
      for( i = 0; i < detectordata->nrelconss + 1; i++ )
      {
         if( SCIPhashmapExists(detectordata->occupied, (void*)(size_t)(i + 1)) )
         {
            detectordata->position = i;

            if( detectordata->algorithm )
            {
               SCIP_CALL( callMetis(scip, detectordata,result) );
               SCIPdebugMessage("Metis successful \n");
            }
            else
            {
               SCIP_CALL( StoerWagner(scip, detectordata) );
               SCIPdebugMessage("StoerWagner successful \n");
            }

            SCIP_CALL( buildnewgraphs(scip, detectordata) );
            SCIPdebugMessage("buildnewgraphs successful \n");
         }
      }
   }
   /** add merged conss */
   SCIP_CALL( getmergedconss(scip, detectordata) );
   SCIPdebugMessage("getmergedconss successful \n");

   /** get subscipvars, copy data to decdecomp */
   SCIP_CALL( GetConsindex(scip, detectordata, (*decdecomps)[0]) );

   if( detectordata->fixedblocks )
   {
        SCIP_CALL( FixedBlocks(scip, detectordata, (*decdecomps)[0]) );
        SCIP_CALL( GetLinkingVars(scip, detectordata, (*decdecomps)[0]) );
   }

   detectordata->found = TRUE;
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "found %d blocks.\n", DECdecompGetNBlocks((*decdecomps)[0]));

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** creates the cutpacking presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionCutpacking(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA *detectordata;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );

   assert(detectordata != NULL);
   detectordata->found = FALSE;
   detectordata->partition = NULL;
   detectordata->nblocks = -1;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP, detectordata, detectAndBuildCutpacking, initCutpacking, exitCutpacking) );

   /* add staircase presolver parameters */
      SCIP_CALL( SCIPaddIntParam(scip, "staircase/algorithm", "should the stoer-wagner algorithm or metis be used for finding a minimal cut", &detectordata->algorithm, FALSE, DEFAULT_ALGORITHM_METIS,  INT_MIN, INT_MAX, NULL, NULL) );
      SCIP_CALL( SCIPaddBoolParam(scip, "staircase/fixedblocks", "Should the blocks consist of a certain number of constraints", &detectordata->fixedblocks, FALSE, DEFAULT_FIXEDBLOCKS, NULL, NULL) );
      SCIP_CALL( SCIPaddIntParam(scip, "staircase/blocksize", "number of constraints per block", &detectordata->blocksize, FALSE, DEFAULT_BLOCKSIZE,  1.0, INT_MAX, NULL, NULL) );
      SCIP_CALL( SCIPaddBoolParam(scip, "staircase/tidy", "Whether to clean up temporary files", &detectordata->tidy, FALSE, DEFAULT_TIDY, NULL, NULL) );
      SCIP_CALL( SCIPaddIntParam(scip, "staircase/randomseed", "random seed for hmetis", &detectordata->randomseed, FALSE, DEFAULT_RANDSEED, -1, INT_MAX, NULL, NULL) );
      SCIP_CALL( SCIPaddRealParam(scip, "staircase/ubfactor", "Unbalance factor for metis", &detectordata->metisubfactor, FALSE, DEFAULT_METIS_UBFACTOR, 0.0, 1E20, NULL, NULL ) );
      SCIP_CALL( SCIPaddBoolParam(scip, "staircase/metisverbose", "Should the metis output be displayed", &detectordata->metisverbose, FALSE, DEFAULT_METIS_VERBOSE, NULL, NULL ) );
      SCIP_CALL( SCIPaddBoolParam(scip, "staircase/metisuseptyperb", "Should the rb or kway method be used for partitioning by metis", &detectordata->metisuseptyperb, FALSE, DEFAULT_METISUSEPTYPE_RB, NULL, NULL) );

   return SCIP_OKAY;
}
