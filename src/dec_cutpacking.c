/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_cutpacking.c
 * @ingroup DETECTORS
 * @brief  cutpacking presolver
 * @author
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>

#include "dec_cutpacking.h"

#include "cons_decomp.h"
#include "struct_decomp.h"
#include "scip_misc.h"

#define DEC_DETECTORNAME      "cutpacking"   /**< name of the detector */
#define DEC_PRIORITY          -500             /**< priority of the detector */

/* Default parameter settings */
#define DEFAULT_PRIORITY                  DEC_PRIORITY

/*
 * Data structures
 */

struct graphstructure
{
   SCIP_HASHMAP** adjacencylist;
   SCIP_CONS** conss;
   int nconss;
   SCIP_HASHMAP* constopos;

   int nedges;

   SCIP_CONS* cons1;
   SCIP_CONS* cons2;
};
typedef struct graphstructure Graph;

/** detector data */
struct DEC_DetectorData
{
   DECDECOMP* decdecomp;

   int nrelconss;
   int nblocks;
   SCIP_CONS*** subscipconss;
   int* nsubscipconss;
   SCIP_VAR*** subscipvars;
   int* nsubscipvars;
   SCIP_VAR** linkingvars;
   int nlinkingvars;

   SCIP_HASHMAP* constoblock;
   SCIP_HASHMAP* varstoblock;
   SCIP_HASHMAP* consindex;
   SCIP_HASHMAP* varindex;

   Graph* graphs;
   int ngraphs;

   SCIP_HASHMAP* occupied;
   int position;
   int* partition;

   SCIP_HASHMAP** mergedconss;
   SCIP_HASHMAP* representatives;
   int nrepresentatives;

   SCIP_HASHMAP* vartopos;
   int* nvarinconss;
   SCIP_CONS*** varinconss;
   SCIP_VAR** relvars;
   int nrelvars;

   int priority;

   int found; //noch default parameter setzten; ende;

   int startblock;
};

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
   detectordata->nlinkingvars = 0;
   detectordata->position = -1;
   detectordata->nrepresentatives = 0;
   detectordata->nrelvars = 0;
   detectordata->startblock = -1;

   /*** get number of relevant variables ***/
   /** vartopos */
   SCIP_CALL( SCIPhashmapCreate(&detectordata->vartopos, SCIPblkmem(scip),nallvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->relvars, nallvars) );
   vartopos = detectordata->vartopos;
   relvars = detectordata->relvars;
   j = 0;
   for( i = 0; i < nallvars; i++ )
   {
      if( isVarRelevant(allvars[i]) )
      {
         relvars[j] = SCIPvarGetProbvar(allvars[i]);
         SCIP_CALL( SCIPhashmapInsert(vartopos, SCIPvarGetProbvar(allvars[i]), (void*) (size_t) j) );
         j++;
      }
   }
   detectordata->nrelvars = j;
   SCIPreallocMemoryArray(scip, &(detectordata->relvars), j);

   /*** get number of relevant conss **/
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->graphs, nconss+1));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(detectordata->graphs[0].conss), nconss));
   k = 0;
   for( i = 0; i < nconss; i++ )
   {
      vars = SCIPgetVarsXXX(scip, conss[i]);
      nvars = SCIPgetNVarsXXX(scip, conss[i]);
      ishandled = FALSE;

      for( j = 0; (j < nvars) && (ishandled == FALSE); j++ )
      {
         ishandled = isVarRelevant(vars[j]);
      }

      if( ishandled )
      {
         detectordata->graphs[0].conss[k] = conss[i];
         k++;
      }SCIPfreeMemoryArrayNull(scip, &vars);
   }
   detectordata->nrelconss = k;
   detectordata->graphs[0].nconss = k;
   SCIPreallocMemoryArray(scip, &(detectordata->graphs[0].conss), k);

   /** alloc **/
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->partition, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->subscipvars, k) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip),k) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->representatives, SCIPblkmem(scip),k) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->occupied, SCIPblkmem(scip),k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->subscipconss, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->mergedconss, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->nsubscipconss, k) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->nsubscipvars, k) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->varstoblock, SCIPblkmem(scip),detectordata->nrelvars) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->varindex, SCIPblkmem(scip),detectordata->nrelvars) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->consindex, SCIPblkmem(scip),detectordata->nrelconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->nvarinconss, detectordata->nrelvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varinconss, detectordata->nrelvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->linkingvars, detectordata->nrelvars) );

   for( i = 0; i < k; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->subscipconss[i], k) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->subscipvars[i], detectordata->nrelvars) );
      SCIP_CALL( SCIPhashmapCreate(&detectordata->mergedconss[i], SCIPblkmem(scip),k) );
      detectordata->nsubscipvars[i] = 0;
   }

   /** varinconss */
   for( i = 0; i < detectordata->nrelvars; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varinconss[i], k) );
      detectordata->nvarinconss[i] = 0;
   }

   for( i = 0; i < k; i++ )
   {
      vars = SCIPgetVarsXXX(scip, detectordata->graphs[0].conss[i]);
      nvars = SCIPgetNVarsXXX(scip, detectordata->graphs[0].conss[i]);
      for( j = 0; j < nvars; j++ )
      {
         if( isVarRelevant(vars[j]) )
         {
            (detectordata->varinconss[(long int)SCIPhashmapGetImage(vartopos, SCIPvarGetProbvar(vars[j]))])[detectordata->nvarinconss[(long int)SCIPhashmapGetImage(
               vartopos, SCIPvarGetProbvar(vars[j]))]] = detectordata->graphs[0].conss[i];
            detectordata->nvarinconss[(long int)SCIPhashmapGetImage(vartopos, SCIPvarGetProbvar(vars[j]))]++;
         }
      }
      SCIPfreeMemoryArrayNull(scip, &vars);
   }

   //printf("nrelvars: %d    nrelconss: %d \n",detectordata->nrelvars, detectordata->nrelconss);

   return SCIP_OKAY;
}

/** copies the variable and block information to the decomp structure */
static SCIP_RETCODE copyDetectorDataToDecomp(SCIP* scip, /**< SCIP data structure */
DEC_DETECTORDATA* detectordata, /**< presolver data data structure */
DECDECOMP* decomp /**< DECOMP data structure */
)
{
   int i;
   SCIP_HASHMAPLIST* list;
   assert(scip != 0);
   assert(detectordata != 0);
   assert(decomp != 0);

   SCIP_CALL( SCIPallocMemoryArray(scip, &decomp->subscipvars, detectordata->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decomp->subscipconss, detectordata->nblocks) );

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &decomp->linkingvars, detectordata->linkingvars, detectordata->nlinkingvars) );
   decomp->nlinkingvars = detectordata->nlinkingvars;

   for( i = 0; i < detectordata->nblocks; i++ )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decomp->subscipconss[i], detectordata->subscipconss[i], detectordata->nsubscipconss[i]) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decomp->subscipvars[i], detectordata->subscipvars[i], detectordata->nsubscipvars[i]) );
   }

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &decomp->nsubscipconss, detectordata->nsubscipconss, detectordata->nblocks) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &decomp->nsubscipvars, detectordata->nsubscipvars, detectordata->nblocks) );

   assert(!SCIPhashmapIsEmpty(detectordata->constoblock));
   assert(!SCIPhashmapIsEmpty(detectordata->varstoblock));

   SCIP_CALL( SCIPhashmapCreate(&decomp->constoblock, SCIPblkmem(scip),detectordata->nrelconss) );
   SCIP_CALL( SCIPhashmapCreate(&decomp->vartoblock, SCIPblkmem(scip),detectordata->nrelvars) );
   SCIP_CALL( SCIPhashmapCreate(&decomp->consindex, SCIPblkmem(scip),detectordata->nrelconss) );
   SCIP_CALL( SCIPhashmapCreate(&decomp->varindex, SCIPblkmem(scip),detectordata->nrelvars) );

   for( i = 0; i < SCIPhashmapGetNLists(detectordata->constoblock); ++i )
   {
      list = SCIPhashmapGetList(detectordata->constoblock, i);
      if( SCIPhashmapListGetNEntries(list) != 0 )
      {
         while( list != NULL )
         {
            SCIP_CALL( SCIPhashmapInsert(decomp->constoblock, SCIPhashmapListGetOrigin(list),SCIPhashmapListGetImage(list)) );
            list = SCIPhashmapListGetNext(list);
         }
      }
   }assert(!SCIPhashmapIsEmpty(decomp->constoblock));
   for( i = 0; i < SCIPhashmapGetNLists(detectordata->varstoblock); ++i )
   {
      list = SCIPhashmapGetList(detectordata->varstoblock, i);

      if( SCIPhashmapListGetNEntries(list) != 0 )
      {
         while( list != NULL )
         {
            SCIP_CALL( SCIPhashmapInsert(decomp->vartoblock, SCIPhashmapListGetOrigin(list), SCIPhashmapListGetImage(list)) );
            list = SCIPhashmapListGetNext(list);
         }
      }
   }assert(!SCIPhashmapIsEmpty(decomp->vartoblock));

   for( i = 0; i < SCIPhashmapGetNLists(detectordata->consindex); ++i )
   {
      list = SCIPhashmapGetList(detectordata->consindex, i);
      if( SCIPhashmapListGetNEntries(list) != 0 )
      {
         while( list != NULL )
         {
            SCIP_CALL( SCIPhashmapInsert(decomp->consindex, SCIPhashmapListGetOrigin(list),SCIPhashmapListGetImage(list)) );
            list = SCIPhashmapListGetNext(list);
         }
      }
   }assert(!SCIPhashmapIsEmpty(decomp->consindex));

   for( i = 0; i < SCIPhashmapGetNLists(detectordata->varindex); ++i )
   {
      list = SCIPhashmapGetList(detectordata->varindex, i);
      if( SCIPhashmapListGetNEntries(list) != 0 )
      {
         while( list != NULL )
         {
            SCIP_CALL( SCIPhashmapInsert(decomp->varindex, SCIPhashmapListGetOrigin(list), SCIPhashmapListGetImage(list)) );
            list = SCIPhashmapListGetNext(list);
         }
      }
   }
   assert(!SCIPhashmapIsEmpty(decomp->varindex));

   decomp->nblocks = detectordata->nblocks;
   decomp->type = DEC_STAIRCASE;

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

   for( i = 0; i < detectordata->nrelconss; i++ )
   {
      SCIPfreeMemoryArray(scip, &detectordata->subscipconss[i]);
      SCIPfreeMemoryArray(scip, &detectordata->subscipvars[i]);
   }

   for( i = 0; i < detectordata->nrelvars; ++i )
   {
      SCIPfreeMemoryArray(scip, &detectordata->varinconss[i]);
   }

   SCIPfreeMemoryArray(scip, &detectordata->nsubscipconss);
   SCIPfreeMemoryArray(scip, &detectordata->subscipconss);
   SCIPfreeMemoryArray(scip, &detectordata->nsubscipvars);
   SCIPfreeMemoryArray(scip, &detectordata->subscipvars);
   SCIPfreeMemoryArray(scip, &detectordata->linkingvars);
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

   SCIPhashmapFree(&detectordata->varstoblock);
   SCIPhashmapFree(&detectordata->constoblock);
   SCIPhashmapFree(&detectordata->consindex);
   SCIPhashmapFree(&detectordata->varindex);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/**
 * Builds a graph structure out of the matrix.
 *
 * The function will create a vertice for every constraint and an edge between two constraints if they have a common variable....
 */
static SCIP_RETCODE buildGraphStructure(SCIP* scip, /**< SCIP data structure */
DEC_DETECTORDATA* detectordata /**< presolver data data structure */
)
{
   int i;
   int j;
   int k;
   int nedges;
   long int cost;

   SCIP_CONS*** varinconss;

   SCIP_CALL( SCIPhashmapCreate(&(detectordata->graphs[0].constopos), SCIPblkmem(scip),detectordata->nrelconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->graphs[0].adjacencylist, detectordata->nrelconss) );
   for( i = 0; i < detectordata->nrelconss; ++i )
   {
      SCIP_CALL( SCIPhashmapCreate(&detectordata->graphs[0].adjacencylist[i], SCIPblkmem(scip),detectordata->nrelconss) );
   }

   nedges = 0;

   /** constopos */
   assert(detectordata->graphs[0].nconss > 0);

   for( i = 0; i < detectordata->graphs[0].nconss; i++ )
   {
      SCIP_CALL(SCIPhashmapInsert(detectordata->graphs[0].constopos, detectordata->graphs[0].conss[i], (void*) (size_t) i));
   }

   /** adjacencylist */

   varinconss = detectordata->varinconss;

   for( i = 0; i < detectordata->nrelvars; i++ )
   {
      for( j = 0; j < detectordata->nvarinconss[i]; j++ )
      {
         for( k = j + 1; k < detectordata->nvarinconss[i]; k++ )
         {
            if( SCIPhashmapExists(
               detectordata->graphs[0].adjacencylist[(long int)SCIPhashmapGetImage(detectordata->graphs[0].constopos,
                  varinconss[i][j])], varinconss[i][k]) )
            {
               cost = (long int)SCIPhashmapGetImage(
                  detectordata->graphs[0].adjacencylist[(long int)SCIPhashmapGetImage(detectordata->graphs[0].constopos,
                     varinconss[i][j])], varinconss[i][k]) + 1;
               SCIP_CALL( SCIPhashmapSetImage(detectordata->graphs[0].adjacencylist[(long int)SCIPhashmapGetImage(detectordata->graphs[0].constopos, varinconss[i][j])], varinconss[i][k], (void*) (size_t) cost) );
               cost = (long int)SCIPhashmapGetImage(
                  detectordata->graphs[0].adjacencylist[(long int)SCIPhashmapGetImage(detectordata->graphs[0].constopos,
                     varinconss[i][k])], varinconss[i][j]) + 1;
               SCIP_CALL( SCIPhashmapSetImage(detectordata->graphs[0].adjacencylist[(long int)SCIPhashmapGetImage(detectordata->graphs[0].constopos, varinconss[i][k])], varinconss[i][j], (void*) (size_t) cost) );
            }
            else
            {
               SCIP_CALL( SCIPhashmapInsert(detectordata->graphs[0].adjacencylist[(long int)SCIPhashmapGetImage(detectordata->graphs[0].constopos, varinconss[i][j])], varinconss[i][k],(void*) (size_t) 1) );
               SCIP_CALL( SCIPhashmapInsert(detectordata->graphs[0].adjacencylist[(long int)SCIPhashmapGetImage(detectordata->graphs[0].constopos, varinconss[i][k])], varinconss[i][j],(void*) (size_t) 1) );
               nedges++;
            }
         }
      }
   }

   detectordata->graphs[0].cons1 = NULL;
   detectordata->graphs[0].cons2 = NULL;
   detectordata->graphs[0].nedges = nedges;

   detectordata->ngraphs = 1;
   SCIP_CALL( SCIPhashmapInsert(detectordata->occupied,(void*)(size_t)1,NULL) );

   return SCIP_OKAY;
}

static SCIP_RETCODE StoerWagner(SCIP* scip, /**< SCIP data struture */
DEC_DETECTORDATA* detectordata /**< presolver data data structure */
)
{
   int i;
   int j;
   int tight;
   double value_act_cut;
   double value_cut;
   SCIP_HASHMAP* tightness;
   SCIP_HASHMAPLIST* list;
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

   SCIPhashmapCreate(&tightness, SCIPblkmem(scip), graph.nconss);
   SCIPallocMemoryArray(scip, &mincut, graph.nconss);
   SCIPallocMemoryArray(scip, &nmerged_conss, graph.nconss);
   SCIPallocMemoryArray(scip, &merged_conss, graph.nconss);
   SCIPallocMemoryArray(scip, &adja, graph.nconss);
   for( i = 0; i < graph.nconss; ++i )
   {
      SCIPallocMemoryArray(scip, &(merged_conss[i]), graph.nconss);
      SCIPhashmapCreate(&(adja[i]), SCIPblkmem(scip), graph.nconss);
   }
   SCIPhashmapCreate(&constopos, SCIPblkmem(scip), graph.nconss);
   SCIPhashmapCreate(&repres_conss, SCIPblkmem(scip), graph.nconss);

   //copy constopos

   for( i = 0; i < SCIPhashmapGetNLists(graph.constopos); ++i )
   {
      list = SCIPhashmapGetList(graph.constopos, i);

      if( SCIPhashmapListGetNEntries(list) > 0 )
      {
         while( list != NULL )
         {
            SCIP_CALL( SCIPhashmapInsert(constopos, SCIPhashmapListGetOrigin(list), SCIPhashmapListGetImage(list)) );
            list = SCIPhashmapListGetNext(list);
         }
      }
   }

   //copy adja

   for( i = 0; i < graph.nconss; ++i )
   {
      for( j = 0; j < SCIPhashmapGetNLists(graph.adjacencylist[i]); ++j )
      {
         list = SCIPhashmapGetList(graph.adjacencylist[i], j);

         if( SCIPhashmapListGetNEntries(list) > 0 )
         {
            while( list != NULL )
            {
               SCIP_CALL( SCIPhashmapInsert(adja[i], SCIPhashmapListGetOrigin(list), SCIPhashmapListGetImage(list)) );
               list = SCIPhashmapListGetNext(list);
            }
         }
      }

   }

   SCIPdebugMessage("stoerwagner \n");

   t = NULL;
   value_cut = INFINITY;
   value_act_cut = INFINITY;
   if( graph.cons1 != NULL)
   {
      s = graph.cons1;
      if( graph.cons2 != NULL)
      {
         t = graph.cons2;
      }
   }
   else
   {
      if( graph.cons2 != NULL)
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
      SCIPhashmapRemoveAll(tightness);
      for( i = 0; i < SCIPhashmapGetNLists(constopos); ++i )
      {
         list = SCIPhashmapGetList(constopos, i);
         while( list != NULL )
         {
            if( SCIPhashmapListGetNEntries(list) > 0 )
            {
               if( SCIPhashmapListGetOrigin(list) != s )
               {
                  SCIP_CALL( SCIPhashmapInsert(tightness, SCIPhashmapListGetOrigin(list), (void*)0) );
               }
            }
            list = SCIPhashmapListGetNext(list);
         }
      }

      assert(SCIPhashmapGetNEntries(tightness) + 1 == SCIPhashmapGetNEntries(constopos));
      last = s;
      next_to_last = s;

      while( (SCIPhashmapGetNEntries(tightness) > 0) )
      {
         next_to_last = last;
         //actualize tightness
         assert( SCIPhashmapExists(constopos, last));
         for( i = 0; i < SCIPhashmapGetNLists(adja[(long int)SCIPhashmapGetImage(constopos, last)]); ++i )
         {
            list = SCIPhashmapGetList(adja[(long int)SCIPhashmapGetImage(constopos, last)], i);

            while( list != NULL )
            {
               if( SCIPhashmapListGetNEntries(list) > 0 )
               {
                  assert( SCIPhashmapExists(constopos, SCIPhashmapListGetOrigin(list)));
                  if( SCIPhashmapExists(tightness, SCIPhashmapListGetOrigin(list)) )
                  {
                     j = ((long int)SCIPhashmapGetImage(tightness, SCIPhashmapListGetOrigin(list)))
                        + ((long int)SCIPhashmapListGetImage(list));
                     SCIP_CALL( SCIPhashmapSetImage(tightness,SCIPhashmapListGetOrigin(list),(void*) (size_t) j ));
                  }
               }
               list = SCIPhashmapListGetNext(list);
            }
         }
         // choose the most tight
         tight = 0;
         for( i = 0; i < SCIPhashmapGetNLists(tightness); ++i )
         {
            list = SCIPhashmapGetList(tightness, i);
            while( list != NULL )
            {
               if( SCIPhashmapListGetNEntries(list) != 0 )
               {
                  if( (long int)SCIPhashmapListGetImage(list) >= tight )
                  {
                     last = (SCIP_CONS*)SCIPhashmapListGetOrigin(list);
                     tight = (long int)SCIPhashmapListGetImage(list);
                  }
               }
               list = SCIPhashmapListGetNext(list);
            }
         }

         assert(next_to_last != last);
         assert(SCIPhashmapExists(tightness, last));
         SCIPhashmapRemove(tightness, last);
      }
      //value of cut
      value_act_cut = 0;

      for( i = 0; i < SCIPhashmapGetNLists(adja[(long int)SCIPhashmapGetImage(constopos, last)]); ++i )
      {
         list = SCIPhashmapGetList(adja[(long int)SCIPhashmapGetImage(constopos, last)], i);
         while( list != NULL )
         {
            if( SCIPhashmapListGetNEntries(list) != 0 )
            {
               assert(SCIPhashmapListGetOrigin(list) != last);
               if( (SCIPhashmapExists(constopos, SCIPhashmapListGetOrigin(list))) )
               {
                  value_act_cut += (long int)SCIPhashmapListGetImage(list);
               }
            }
            list = SCIPhashmapListGetNext(list);
         }
      }
      if( (value_act_cut < value_cut) )
      {
         /*if( (t != NULL) && (represent_t == last) )
         {
            //test if act_cut ist a s-t-cut
            represent_t = next_to_last;
            value_cut = value_act_cut;
            cut = last;
         }
         else if( t == NULL)
         {*/
            value_cut = value_act_cut;
            cut = last;
        // }
      }

      //mergen
      if( SCIPhashmapGetNEntries(repres_conss) == 0 )
      {
         SCIP_CALL( SCIPhashmapInsert(repres_conss, next_to_last, (void*) (size_t) nrepres_conss));
         nrepres_conss++; //starts with 1
         nmerged_conss[nrepres_conss - 2] = 1;
         merged_conss[nrepres_conss - 2][0] = last;
      }
      else if( !SCIPhashmapExists(repres_conss, next_to_last) )
      {
         SCIP_CALL( SCIPhashmapInsert(repres_conss, next_to_last, (void*) (size_t) nrepres_conss));
         nrepres_conss++; //starts with 1
         nmerged_conss[nrepres_conss - 2] = 1;
         merged_conss[nrepres_conss - 2][0] = last;
      }
      else
      {
         merged_conss[((long int)SCIPhashmapGetImage(repres_conss, next_to_last)) - 1][nmerged_conss[((long int)SCIPhashmapGetImage(
            repres_conss, next_to_last)) - 1]] = last;
         nmerged_conss[((long int)SCIPhashmapGetImage(repres_conss, next_to_last)) - 1]++;
      }

      //in adja: connect last and next_to_last

      for( i = 0; i < SCIPhashmapGetNLists(adja[(long int)SCIPhashmapGetImage(constopos, last)]); ++i )
      {
         list = SCIPhashmapGetList(adja[(long int)SCIPhashmapGetImage(constopos, last)], i);
         while( list != NULL )
         {
            if( SCIPhashmapListGetNEntries(list) > 0 )
            {
               if( SCIPhashmapExists(adja[(long int)SCIPhashmapGetImage(constopos, next_to_last)],
                  SCIPhashmapListGetOrigin(list)) )
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
            }
            list = SCIPhashmapListGetNext(list);
         }
      }

      SCIPhashmapRemove(constopos, last);

      // delete last in adja
      for( i = 0; i < graph.nconss; ++i )
      {
         if( SCIPhashmapExists(constopos, graph.conss[i]) )
         {
            if( SCIPhashmapExists(adja[i], last) )
            {
               SCIPhashmapRemove(adja[i], last);
            }
         }
      }
   }

   //complete cut
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

   //create partition
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

static SCIP_RETCODE buildnewgraphs(SCIP* scip, DEC_DETECTORDATA* detectordata)
{
   int i;
   int j;
   int k;
   int pos1;
   int pos2;
   int cas;
   int cost;
   int nedges;
   int nsconss1;
   int nsconss2;
   int stop1;
   int stop2;

   int* partition;
   Graph graph;
   SCIP_CONS* representative;

   SCIP_HASHMAP* adja;
   SCIP_HASHMAP* newadja1;
   SCIP_HASHMAP* newadja2;
   SCIP_HASHMAPLIST* list;

   int nconss1;
   SCIP_HASHMAP* consslink1;
   int nconss2;
   SCIP_HASHMAP* consslink2;

   cas = -1;
   nconss1 = 0;
   nconss2 = 0;
   stop1 = 0;
   stop2 = 0;

   /** build partitions */
   partition = detectordata->partition;

   graph = detectordata->graphs[detectordata->position];

   j = 0;
   pos1 = -1;
   pos2 = -1;
   for( i = 0; (j < 2) && (i < detectordata->nrelconss + 1); i++ )
   {
      if( (!SCIPhashmapExists(detectordata->occupied, (void*)(size_t)(i + 1))) )
      {
         if( pos1 == -1 )
         {
            pos1 = i;
            j++;
         }
         else
         {
            pos2 = i;
            j++;
         }
      }
   }assert(i != detectordata->nrelconss+1);
   assert((pos1 != -1) && (pos2 != -1));

   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->graphs[pos1].conss, graph.nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->graphs[pos2].conss, graph.nconss));
   SCIP_CALL( SCIPhashmapRemove(detectordata->occupied, (void*)(size_t)(detectordata->position + 1)));

   for( i = 0; i < graph.nconss; i++ )
   {
      assert((-1 < partition[i])&&(partition[i]<2));
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

   nsconss1 = nconss1;
   nsconss2 = nconss2;
   SCIP_CALL( SCIPhashmapCreate(&detectordata->graphs[pos1].constopos, SCIPblkmem(scip),nconss1) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->graphs[pos2].constopos, SCIPblkmem(scip),nconss2) );
   SCIPallocMemoryArray(scip, &detectordata->graphs[pos1].adjacencylist, nconss1);
   SCIPallocMemoryArray(scip, &detectordata->graphs[pos2].adjacencylist, nconss2);

   for( i = 0; i < nconss1; ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(detectordata->graphs[pos1].constopos, detectordata->graphs[pos1].conss[i], NULL) );
      SCIP_CALL( SCIPhashmapCreate(&detectordata->graphs[pos1].adjacencylist[i], SCIPblkmem(scip),nconss1) );
   }
   for( i = 0; i < nconss2; ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(detectordata->graphs[pos2].constopos, detectordata->graphs[pos2].conss[i], NULL) );
      SCIP_CALL( SCIPhashmapCreate(&detectordata->graphs[pos2].adjacencylist[i], SCIPblkmem(scip),nconss2) );
   }

   assert(nconss1 + nconss2 == graph.nconss);
   assert((nconss1 != 0) && (nconss2 != 0));

   /** get linking conss **/
   SCIP_CALL( SCIPhashmapCreate(&consslink1, SCIPblkmem(scip), graph.nconss) );
   SCIP_CALL( SCIPhashmapCreate(&consslink2, SCIPblkmem(scip), graph.nconss) );
   assert(SCIPhashmapIsEmpty(consslink1));
   assert(SCIPhashmapIsEmpty(consslink2));

   for( i = 0; i < nconss1; ++i )
   {
      adja = graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, detectordata->graphs[pos1].conss[i])];

      for( j = 0; j < SCIPhashmapGetNLists(adja); ++j )
      {
         list = SCIPhashmapGetList(adja, j);
         if( SCIPhashmapListGetNEntries(list) != 0 )
         {
            while( list != NULL )
            {
               if( SCIPhashmapExists(detectordata->graphs[pos2].constopos, SCIPhashmapListGetOrigin(list)) )
               {
                  if( !SCIPhashmapExists(consslink1, (void*)detectordata->graphs[pos1].conss[i]) )
                  {
                     SCIP_CALL( SCIPhashmapInsert(consslink1, detectordata->graphs[pos1].conss[i], NULL) );
                  }
                  if( !SCIPhashmapExists(consslink2, SCIPhashmapListGetOrigin(list)) )
                  {
                     SCIP_CALL( SCIPhashmapInsert(consslink2, SCIPhashmapListGetOrigin(list), NULL) );
                  }
               }
               list = SCIPhashmapListGetNext(list);
            }
         }
      }
   }
   //printf("nconsslink1: %d   nconsslink2: %d \n", SCIPhashmapGetNEntries(consslink1), SCIPhashmapGetNEntries(consslink2));

   /** tests whether the cut is feasible? **/
   //SCIPdebugMessage("cut feasible? \n");
   if( (graph.cons1 != NULL) && (graph.cons2 != NULL) )
   {
      if( (SCIPhashmapExists(detectordata->graphs[pos1].constopos, graph.cons1)
         && SCIPhashmapExists(detectordata->graphs[pos1].constopos, graph.cons2))
         || (SCIPhashmapExists(detectordata->graphs[pos2].constopos, graph.cons1)
            && SCIPhashmapExists(detectordata->graphs[pos2].constopos, graph.cons2)) )
      {
         for( i = 0; i < graph.nconss; ++i )
         {
            detectordata->subscipconss[detectordata->nblocks][i] = graph.conss[i];
         }
         detectordata->nsubscipconss[detectordata->nblocks] = graph.nconss;
         detectordata->nblocks++;
         detectordata->ngraphs--;
         for( i = 0; i < nconss1; ++i )
         {
            SCIPhashmapFree(&detectordata->graphs[pos1].adjacencylist[i]);
         }
         for( i = 0; i < nconss2; ++i )
         {
            SCIPhashmapFree(&detectordata->graphs[pos2].adjacencylist[i]);
         }
         SCIPhashmapFree(&detectordata->graphs[pos1].constopos);
         SCIPhashmapFree(&detectordata->graphs[pos2].constopos);
         SCIPfreeMemoryArray(scip, &detectordata->graphs[pos1].adjacencylist);
         SCIPfreeMemoryArray(scip, &detectordata->graphs[pos2].adjacencylist);
         SCIPfreeMemoryArray(scip, &detectordata->graphs[pos1].conss);
         SCIPfreeMemoryArray(scip, &detectordata->graphs[pos2].conss);
      }
      else
      {
         if( SCIPhashmapExists(detectordata->graphs[pos1].constopos, graph.cons1) )
         {
            cas = 0;
            if( SCIPhashmapExists(consslink1, graph.cons1) )
            {
               stop1 = 1;
            }
            if( SCIPhashmapExists(consslink2, graph.cons2) )
            {
               stop2 = 1;
            }
         }
         else
         {
            cas = 1;
            if( SCIPhashmapExists(consslink1, graph.cons2) )
            {
               stop1 = 1;
            }
            if( SCIPhashmapExists(consslink2, graph.cons1) )
            {
               stop2 = 1;
            }
         }
      }
   }
   else if( (graph.cons1 != NULL) && (graph.cons2 == NULL) )
   {
      if( SCIPhashmapExists(detectordata->graphs[pos1].constopos, graph.cons1) )
      {
         cas = 0;
         if( SCIPhashmapExists(consslink1, graph.cons1) )
         {
            stop1 = 1;
         }
      }
      else
      {
         cas = 1;
         if( SCIPhashmapExists(consslink2, graph.cons1) )
         {
            stop2 = 1;
         }
      }
   }
   else if( (graph.cons1 == NULL) && (graph.cons2 != NULL) )
   {
      if( SCIPhashmapExists(detectordata->graphs[pos2].constopos, graph.cons2) )
      {
         cas = 0;
         if( SCIPhashmapExists(consslink2, graph.cons2) )
         {
            stop2 = 1;
         }
      }
      else
      {
         cas = 1;
         if( SCIPhashmapExists(consslink1, graph.cons2) )
         {
            stop1 = 1;
         }
      }
   }
   else
   {
      cas = 1;
   }

   /** get contraints to merge */
   //SCIPdebugMessage("contraints to merge \n");
   if( cas != -1 )
   {

      /*********************/
      /** build new graphs */
      /*********************/

      //SCIPdebugMessage("build new graphs \n");
      /** graph1 */
      /** adjacencylist */

      if( (nconss1 > 1) && (nsconss1 != SCIPhashmapGetNEntries(consslink1)) && (stop1 == 0) )
      {
         if( SCIPhashmapGetNEntries(consslink1) > 0 )
         {
            newadja1 = detectordata->graphs[pos1].adjacencylist[0];
            representative = NULL;
            j = 1;
            for( i = 0; i < nconss1; i++ )
            {
               assert(consslink1 != NULL);
               if( !SCIPhashmapExists(consslink1, detectordata->graphs[pos1].conss[i]) )
               {
                  for( k = 0;
                     k
                        < SCIPhashmapGetNLists(
                           graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,
                              detectordata->graphs[pos1].conss[i])]); ++k )
                     {
                     list = SCIPhashmapGetList(
                        graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,
                           detectordata->graphs[pos1].conss[i])], k);
                     if( SCIPhashmapListGetNEntries(list) != 0 )
                     {
                        while( list != NULL )
                        {
                           SCIP_CALL(
                              SCIPhashmapInsert(detectordata->graphs[pos1].adjacencylist[j],SCIPhashmapListGetOrigin(list),SCIPhashmapListGetImage(list)));
                           list = SCIPhashmapListGetNext(list);
                        }
                     }
                  }

                  SCIP_CALL(
                     SCIPhashmapSetImage(detectordata->graphs[pos1].constopos, detectordata->graphs[pos1].conss[i], (void*) (size_t) j));
                  j++;
               }
               else
               {
                  representative = detectordata->graphs[pos1].conss[i];
                  SCIP_CALL(SCIPhashmapRemove(detectordata->graphs[pos1].constopos, detectordata->graphs[pos1].conss[i]));
                  for( k = 0;
                     k
                        < SCIPhashmapGetNLists(
                           graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,
                              detectordata->graphs[pos1].conss[i])]); ++k )
                     {
                     list = SCIPhashmapGetList(
                        graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,
                           detectordata->graphs[pos1].conss[i])], k);
                     if( SCIPhashmapListGetNEntries(list) != 0 )
                     {
                        while( list != NULL )
                        {
                           if( !SCIPhashmapExists(consslink1, SCIPhashmapListGetOrigin(list)) )
                           {
                              if( SCIPhashmapExists(newadja1, SCIPhashmapListGetOrigin(list)) )
                              {
                                 cost = (long int)SCIPhashmapGetImage(newadja1, SCIPhashmapListGetOrigin(list));
                                 cost += (long int)SCIPhashmapListGetImage(list);
                                 SCIP_CALL( SCIPhashmapSetImage(newadja1, SCIPhashmapListGetOrigin(list), (void*) (size_t) cost) );
                              }
                              else if( !SCIPhashmapExists(consslink2, SCIPhashmapListGetOrigin(list)) )
                              {
                                 SCIP_CALL( SCIPhashmapInsert(newadja1, SCIPhashmapListGetOrigin(list), SCIPhashmapListGetImage(list)) );
                              }
                           }
                           list = SCIPhashmapListGetNext(list);
                        }
                     }
                  }
               }
            }
            SCIP_CALL( SCIPhashmapInsert(detectordata->graphs[pos1].constopos, representative, (void*) (size_t) 0) );
            nconss1 = j;
            if( nsconss1 != SCIPhashmapGetNEntries(consslink1) )
            {
               SCIP_CALL( SCIPhashmapInsert(detectordata->representatives, (void*) (size_t) (detectordata->nrepresentatives + 1), representative) );
               for( i = 0; i < SCIPhashmapGetNLists(consslink1); ++i )
               {
                  list = SCIPhashmapGetList(consslink1, i);
                  if( SCIPhashmapListGetNEntries(list) != 0 )
                  {
                     while( list != NULL )
                     {
                        SCIP_CALL(
                           SCIPhashmapInsert(detectordata->mergedconss[detectordata->nrepresentatives], SCIPhashmapListGetOrigin(list),SCIPhashmapListGetImage(list)));
                        list = SCIPhashmapListGetNext(list);
                     }
                  }
               }
               detectordata->nrepresentatives++;
            }
         }
         else
         {
            for( i = 0; i < nconss1; i++ )
            {
               for( i = 0; i < SCIPhashmapGetNLists( graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, detectordata->graphs[pos1].conss[i])]); ++i )
               {
                  list = SCIPhashmapGetList(graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos, detectordata->graphs[pos1].conss[i])], i);
                  if( SCIPhashmapListGetNEntries(list) != 0 )
                  {
                     while( list != NULL )
                     {
                        SCIP_CALL( SCIPhashmapInsert(detectordata->graphs[pos1].adjacencylist[i], SCIPhashmapListGetOrigin(list),SCIPhashmapListGetImage(list)));
                        list = SCIPhashmapListGetNext(list);
                     }
                  }
               }
               SCIP_CALL(
                  SCIPhashmapSetImage(detectordata->graphs[pos1].constopos, detectordata->graphs[pos1].conss[i], (void*) (size_t) i));
            }
         }

         nedges = 0;
         //SCIPdebugMessage("delete merged conss\n");
         /** delete merged conss */
         for( i = 1; i < nconss1; i++ )
         {
            cost = 0;
            for( j = 0; j < SCIPhashmapGetNLists(detectordata->graphs[pos1].adjacencylist[i]); j++ )
            {
               list = SCIPhashmapGetList(detectordata->graphs[pos1].adjacencylist[i], j);
               if( SCIPhashmapListGetNEntries(list) != 0 )
               {
                  while( list != NULL )
                  {
                     if( SCIPhashmapExists(consslink1, SCIPhashmapListGetOrigin(list)) )
                     {
                        cost += (long int)SCIPhashmapListGetImage(list);
                        SCIP_CALL(
                           SCIPhashmapRemove(detectordata->graphs[pos1].adjacencylist[i],SCIPhashmapListGetOrigin(list)));
                     }
                     else
                     {
                        nedges++;
                     }
                     list = SCIPhashmapListGetNext(list);
                  }

               }
            }
            if( cost > 0 )
            {
               SCIP_CALL(
                  SCIPhashmapInsert(detectordata->graphs[pos1].adjacencylist[i],representative, (void*) (size_t) cost));
               nedges += 2;
            }
         }

         /** arranges conss1 */
         for( i = 0; i < SCIPhashmapGetNLists(detectordata->graphs[pos1].constopos); ++i )
         {
            list = SCIPhashmapGetList(detectordata->graphs[pos1].constopos, i);
            if( SCIPhashmapListGetNEntries(list) != 0 )
            {
               while( list != NULL )
               {
                  detectordata->graphs[pos1].conss[(long int)SCIPhashmapListGetImage(list)] =
                     (SCIP_CONS*)SCIPhashmapListGetOrigin(list);
                  list = SCIPhashmapListGetNext(list);
               }
            }
         }

         assert(nconss1 + SCIPhashmapGetNEntries(consslink1) -1 == graph.nconss - nconss2);
         for( i = nconss1; i < nconss1 + SCIPhashmapGetNEntries(consslink1) - 1; ++i )
         {
            SCIPhashmapFree(&detectordata->graphs[pos1].adjacencylist[i]);
         }

         SCIPreallocMemoryArray(scip, &detectordata->graphs[pos1].adjacencylist, nconss1);

         /** builds graph1 */
         detectordata->graphs[pos1].nconss = nconss1;
         //printf("nconss1: %d \n", nconss1);
         detectordata->graphs[pos1].nedges = nedges / 2;
         assert( 2*detectordata->graphs[pos1].nedges == nedges);

         switch( cas )
         {
         case 0:
            detectordata->graphs[pos1].cons1 = graph.cons1;
            detectordata->graphs[pos1].cons2 = representative;
            break;

         case 1:
            detectordata->graphs[pos1].cons1 = representative;
            detectordata->graphs[pos1].cons2 = graph.cons2;
            break;

         default:
            break;

            assert(representative != NULL);
            assert(SCIPhashmapExists(detectordata->graphs[pos1].constopos, representative));
         }

      }
      else if( (nsconss1 == SCIPhashmapGetNEntries(consslink1)) )
      {
         switch( cas )
         {
         case 0:
            detectordata->graphs[pos1].cons1 = graph.cons1;
            detectordata->graphs[pos1].cons2 = detectordata->graphs[pos1].conss[0];
            break;

         case 1:
            detectordata->graphs[pos1].cons1 = detectordata->graphs[pos1].conss[0];
            ;
            detectordata->graphs[pos1].cons2 = graph.cons2;
            break;

         default:
            break;
         }
      }
      else if( (nsconss1 != SCIPhashmapGetNEntries(consslink1)) && (stop1 == 1) )
      {
         switch( cas )
         {
         case 0:
            detectordata->graphs[pos1].cons1 = graph.cons1;
            detectordata->graphs[pos1].cons2 = graph.cons1;
            break;

         case 1:
            detectordata->graphs[pos1].cons1 = graph.cons2;
            detectordata->graphs[pos1].cons2 = graph.cons2;
            break;

         default:
            break;
         }
      }

      /** graph 2 */
      /** adjacencylist */

      if( (nconss2 > 1) && (nsconss2 != SCIPhashmapGetNEntries(consslink2)) && (stop2 == 0) )
      {
         if( SCIPhashmapGetNEntries(consslink2) > 0 )
         {
            representative = NULL;
            j = 1;
            newadja2 = detectordata->graphs[pos2].adjacencylist[0];
            for( i = 0; i < nconss2; i++ )
            {
               assert(consslink2 != NULL);
               if( !SCIPhashmapExists(consslink2, detectordata->graphs[pos2].conss[i]) )
               {
                  for( k = 0;
                     k
                        < SCIPhashmapGetNLists(
                           graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,
                              detectordata->graphs[pos2].conss[i])]); ++k )
                     {
                     list = SCIPhashmapGetList(
                        graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,
                           detectordata->graphs[pos2].conss[i])], k);
                     if( SCIPhashmapListGetNEntries(list) != 0 )
                     {
                        while( list != NULL )
                        {
                           SCIP_CALL(
                              SCIPhashmapInsert(detectordata->graphs[pos2].adjacencylist[j],SCIPhashmapListGetOrigin(list),SCIPhashmapListGetImage(list)));
                           list = SCIPhashmapListGetNext(list);
                        }
                     }
                  }SCIP_CALL(
                     SCIPhashmapSetImage(detectordata->graphs[pos2].constopos, detectordata->graphs[pos2].conss[i], (void*) (size_t) j));
                  j++;
               }
               else
               {
                  representative = detectordata->graphs[pos2].conss[i];
                  SCIP_CALL(SCIPhashmapRemove(detectordata->graphs[pos2].constopos, detectordata->graphs[pos2].conss[i]));
                  for( k = 0;
                     k
                        < SCIPhashmapGetNLists(
                           graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,
                              detectordata->graphs[pos2].conss[i])]); ++k )
                     {
                     list = SCIPhashmapGetList(
                        graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,
                           detectordata->graphs[pos2].conss[i])], k);
                     if( SCIPhashmapListGetNEntries(list) != 0 )
                     {
                        while( list != NULL )
                        {
                           if( !SCIPhashmapExists(consslink2, SCIPhashmapListGetOrigin(list)) )
                           {
                              if( SCIPhashmapExists(newadja2, SCIPhashmapListGetOrigin(list)) )
                              {
                                 cost = (long int)SCIPhashmapGetImage(newadja2, SCIPhashmapListGetOrigin(list));
                                 cost += (long int)SCIPhashmapListGetImage(list);
                                 SCIP_CALL(
                                    SCIPhashmapSetImage(newadja2, SCIPhashmapListGetOrigin(list), (void*) (size_t) cost));
                              }
                              else if( !SCIPhashmapExists(consslink1, SCIPhashmapListGetOrigin(list)) )
                              {
                                 SCIP_CALL(
                                    SCIPhashmapInsert(newadja2, SCIPhashmapListGetOrigin(list), SCIPhashmapListGetImage(list)));
                              }
                           }
                           list = SCIPhashmapListGetNext(list);
                        }
                     }
                  }
               }
            }
            nconss2 = j;
            SCIP_CALL(SCIPhashmapInsert(detectordata->graphs[pos2].constopos, representative,(void*) (size_t) 0));
            if( nsconss2 != SCIPhashmapGetNEntries(consslink2) )
            {
               SCIP_CALL(
                  SCIPhashmapInsert(detectordata->representatives, (void*) (size_t) (detectordata->nrepresentatives +1), representative ));
               for( i = 0; i < SCIPhashmapGetNLists(consslink2); ++i )
               {
                  list = SCIPhashmapGetList(consslink2, i);
                  if( SCIPhashmapListGetNEntries(list) != 0 )
                  {
                     while( list != NULL )
                     {
                        SCIP_CALL(
                           SCIPhashmapInsert(detectordata->mergedconss[detectordata->nrepresentatives], SCIPhashmapListGetOrigin(list),SCIPhashmapListGetImage(list)));
                        list = SCIPhashmapListGetNext(list);
                     }
                  }
               }
               detectordata->nrepresentatives++;
            }
         }
         else
         {
            for( i = 0; i < nconss2; i++ )
            {
               for( i = 0;
                  i
                     < SCIPhashmapGetNLists(
                        graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,
                           detectordata->graphs[pos2].conss[i])]); ++i )
                  {
                  list =
                     SCIPhashmapGetList(
                        graph.adjacencylist[(long int)SCIPhashmapGetImage(graph.constopos,
                           detectordata->graphs[pos2].conss[i])], i);
                  if( SCIPhashmapListGetNEntries(list) != 0 )
                  {
                     while( list != NULL )
                     {
                        SCIP_CALL(
                           SCIPhashmapInsert(detectordata->graphs[pos2].adjacencylist[i], SCIPhashmapListGetOrigin(list),SCIPhashmapListGetImage(list)));
                        list = SCIPhashmapListGetNext(list);
                     }
                  }
               }
               SCIP_CALL(
                  SCIPhashmapSetImage(detectordata->graphs[pos2].constopos, detectordata->graphs[pos2].conss[i], (void*) (size_t) i));
            }
         }

         /** delete merged conss */
         nedges = 0;
         for( i = 1; i < nconss2; i++ )
         {
            cost = 0;
            for( j = 0; j < SCIPhashmapGetNLists(detectordata->graphs[pos2].adjacencylist[i]); j++ )
            {
               list = SCIPhashmapGetList(detectordata->graphs[pos2].adjacencylist[i], j);
               if( SCIPhashmapListGetNEntries(list) != 0 )
               {
                  while( list != NULL )
                  {
                     if( SCIPhashmapExists(consslink2, SCIPhashmapListGetOrigin(list)) )
                     {
                        cost += (long int)SCIPhashmapListGetImage(list);
                        SCIP_CALL(
                           SCIPhashmapRemove(detectordata->graphs[pos2].adjacencylist[i],SCIPhashmapListGetOrigin(list)));
                     }
                     else
                     {
                        nedges++;
                     }
                     list = SCIPhashmapListGetNext(list);
                  }

               }
            }
            if( cost > 0 )
            {
               SCIP_CALL(
                  SCIPhashmapInsert(detectordata->graphs[pos2].adjacencylist[i],representative, (void*) (size_t) cost));
               nedges += 2;
            }
         }

         /** arranges conss2 */
         for( i = 0; i < SCIPhashmapGetNLists(detectordata->graphs[pos2].constopos); ++i )
         {
            list = SCIPhashmapGetList(detectordata->graphs[pos2].constopos, i);
            if( SCIPhashmapListGetNEntries(list) != 0 )
            {
               while( list != NULL )
               {
                  detectordata->graphs[pos2].conss[(long int)SCIPhashmapListGetImage(list)] =
                     (SCIP_CONS*)SCIPhashmapListGetOrigin(list);
                  list = SCIPhashmapListGetNext(list);
               }
            }
         }

         for( i = nconss2; i < nconss2 + SCIPhashmapGetNEntries(consslink2) - 1; ++i )
         {
            SCIPhashmapFree(&detectordata->graphs[pos2].adjacencylist[i]);
         }

         SCIPreallocMemoryArray(scip, &detectordata->graphs[pos2].adjacencylist, nconss2);
         /** builds graph 2 */
         detectordata->graphs[pos2].nconss = nconss2;
         //printf("nconss2: %d \n", nconss2);
         detectordata->graphs[pos2].nedges = nedges / 2;
         assert( 2 * detectordata->graphs[pos2].nedges == nedges);

         if( (nconss1 != SCIPhashmapGetNEntries(consslink1)) && (stop1 == 0) )
         {
            assert(
               nconss1 + SCIPhashmapGetNEntries(consslink1)+ nconss2 + SCIPhashmapGetNEntries(consslink2) -2 == graph.nconss);
         }
         else
         {
            assert(nconss1 + nconss2 + SCIPhashmapGetNEntries(consslink2)-1 == graph.nconss);
         }

         switch( cas )
         {
         case 0:
            detectordata->graphs[pos2].cons1 = detectordata->graphs[pos2].conss[0];
            detectordata->graphs[pos2].cons2 = graph.cons2;
            break;

         case 1:
            detectordata->graphs[pos2].cons1 = graph.cons1;
            detectordata->graphs[pos2].cons2 = detectordata->graphs[pos2].conss[0];
            break;

         default:
            break;
         }assert(SCIPhashmapExists(detectordata->graphs[pos2].constopos, representative));
         assert(representative != NULL);
      }
      else if( nsconss2 == SCIPhashmapGetNEntries(consslink2) )
      {
         switch( cas )
         {
         case 0:
            detectordata->graphs[pos2].cons1 = detectordata->graphs[pos2].conss[0];
            detectordata->graphs[pos2].cons2 = graph.cons2;
            break;

         case 1:
            detectordata->graphs[pos2].cons1 = graph.cons1;
            detectordata->graphs[pos2].cons2 = detectordata->graphs[pos2].conss[0];
            break;

         default:
            break;
         }
      }
      else if( (nsconss2 != SCIPhashmapGetNEntries(consslink2)) && (stop2 == 1) )
      {
         switch( cas )
         {
         case 0:
            detectordata->graphs[pos2].cons1 = graph.cons2;
            detectordata->graphs[pos2].cons2 = graph.cons2;
            break;

         case 1:
            detectordata->graphs[pos2].cons1 = graph.cons1;
            detectordata->graphs[pos2].cons2 = graph.cons1;
            break;

         default:
            break;
         }
      }

      if( ((nconss1 < 2) && (nconss2 < 2))
         || ((nsconss1 == SCIPhashmapGetNEntries(consslink1)) && (nsconss2 == SCIPhashmapGetNEntries(consslink2)))
         || ((stop1 == 1) && (stop2 == 1)) || ((nsconss1 == SCIPhashmapGetNEntries(consslink1)) && (stop2 == 1))
         || ((nsconss2 == SCIPhashmapGetNEntries(consslink2)) && (stop1 == 1)) )
      {
         for( i = 0; i < nconss1; i++ )
         {
            detectordata->subscipconss[detectordata->nblocks][i] = detectordata->graphs[pos1].conss[i];
         }
         detectordata->nsubscipconss[detectordata->nblocks] = nconss1;
         detectordata->nblocks++;
         if( (detectordata->graphs[pos1].cons1 == NULL) )
         {
            detectordata->startblock = detectordata->nblocks;
         }
         for( i = 0; i < nconss2; i++ )
         {
            detectordata->subscipconss[detectordata->nblocks][i] = detectordata->graphs[pos2].conss[i];
         }
         detectordata->nsubscipconss[detectordata->nblocks] = nconss2;
         detectordata->nblocks++;
         detectordata->ngraphs--;
         if( (detectordata->graphs[pos2].cons1 == NULL) )
         {
            detectordata->startblock = detectordata->nblocks;
         }

         SCIPfreeMemoryArray(scip, &detectordata->graphs[pos1].conss);
         SCIPfreeMemoryArray(scip, &detectordata->graphs[pos2].conss);
         SCIPhashmapFree(&detectordata->graphs[pos1].constopos);
         SCIPhashmapFree(&detectordata->graphs[pos2].constopos);
         for( i = 0; i < nconss1; ++i )
         {
            SCIPhashmapFree(&detectordata->graphs[pos1].adjacencylist[i]);
         }
         for( i = 0; i < nconss2; ++i )
         {
            SCIPhashmapFree(&detectordata->graphs[pos2].adjacencylist[i]);
         }SCIPfreeMemoryArray(scip, &detectordata->graphs[pos1].adjacencylist);
         SCIPfreeMemoryArray(scip, &detectordata->graphs[pos2].adjacencylist);
      }
      else if( (nconss1 < 2) || ((stop1 == 1) && (stop2 == 0))
         || ((nsconss1 == SCIPhashmapGetNEntries(consslink1)) && (nsconss2 != SCIPhashmapGetNEntries(consslink2)))
         || ((nsconss1 == SCIPhashmapGetNEntries(consslink1)) && (stop2 == 0))
         || ((stop1 == 1) && (nsconss2 != SCIPhashmapGetNEntries(consslink2))) )
      {
         for( i = 0; i < nconss1; i++ )
         {
            detectordata->subscipconss[detectordata->nblocks][i] = detectordata->graphs[pos1].conss[i];
         }
         detectordata->nsubscipconss[detectordata->nblocks] = nconss1;
         detectordata->nblocks++;
         if( detectordata->graphs[pos1].cons1 == NULL)
         {
            detectordata->startblock = detectordata->nblocks;
         }

         SCIP_CALL( SCIPhashmapInsert(detectordata->occupied,(void*) (size_t) (pos2 + 1), NULL));

         SCIPfreeMemoryArray(scip, &detectordata->graphs[pos1].conss);
         SCIPhashmapFree(&detectordata->graphs[pos1].constopos);
         for( i = 0; i < nconss1; ++i )
         {
            SCIPhashmapFree(&detectordata->graphs[pos1].adjacencylist[i]);
         }SCIPfreeMemoryArray(scip, &detectordata->graphs[pos1].adjacencylist);
      }
      else if( (nconss2 < 2) || ((stop1 == 0) && (stop2 == 1))
         || ((nsconss2 == SCIPhashmapGetNEntries(consslink2)) && (nsconss1 != SCIPhashmapGetNEntries(consslink1)))
         || ((nsconss2 == SCIPhashmapGetNEntries(consslink2)) && (stop1 == 0))
         || ((stop2 == 1) && (nsconss1 != SCIPhashmapGetNEntries(consslink1))) )
      {
         for( i = 0; i < nconss2; i++ )
         {
            detectordata->subscipconss[detectordata->nblocks][i] = detectordata->graphs[pos2].conss[i];
         }
         detectordata->nsubscipconss[detectordata->nblocks] = nconss2;
         detectordata->nblocks++;
         if( detectordata->graphs[pos2].cons1 == NULL)
         {
            detectordata->startblock = detectordata->nblocks;
         }

         SCIP_CALL( SCIPhashmapInsert(detectordata->occupied,(void*) (size_t) (pos1 + 1), NULL));
         SCIPfreeMemoryArray(scip, &detectordata->graphs[pos2].conss);
         SCIPhashmapFree(&detectordata->graphs[pos2].constopos);
         for( i = 0; i < nconss2; ++i )
         {
            SCIPhashmapFree(&detectordata->graphs[pos2].adjacencylist[i]);
         }SCIPfreeMemoryArray(scip, &detectordata->graphs[pos2].adjacencylist);
      }
      else
      {
         SCIP_CALL( SCIPhashmapInsert(detectordata->occupied,(void*) (size_t) (pos1 + 1), NULL));
         SCIP_CALL( SCIPhashmapInsert(detectordata->occupied,(void*) (size_t) (pos2 + 1), NULL));
         detectordata->ngraphs++;
      }
   }

   SCIPhashmapFree(&detectordata->graphs[detectordata->position].constopos);
   for( i = 0; i < graph.nconss; ++i )
   {
      SCIPhashmapFree(&detectordata->graphs[detectordata->position].adjacencylist[i]);
   }SCIPfreeMemoryArray(scip, &detectordata->graphs[detectordata->position].adjacencylist);
   SCIPfreeMemoryArray(scip, &detectordata->graphs[detectordata->position].conss);

   SCIPhashmapFree(&consslink1);
   SCIPhashmapFree(&consslink2);

   return SCIP_OKAY;
}

static SCIP_RETCODE getmergedconss(SCIP* scip, DEC_DETECTORDATA* detectordata)
{
   int i;
   int j;
   int block;
   SCIP_HASHMAP** mergedconss;
   SCIP_HASHMAP* representatives;
   SCIP_CONS*** subscipconss;
   SCIP_HASHMAPLIST* list;
   SCIP_HASHMAP* constoblock;

   mergedconss = detectordata->mergedconss;
   representatives = detectordata->representatives;
   subscipconss = detectordata->subscipconss;
   constoblock = detectordata->constoblock;

   /*** constoblock ***/
   for( i = 0; i < detectordata->nblocks; i++ )
   {
      for( j = 0; j < detectordata->nsubscipconss[i]; j++ )
      {
         if( !SCIPhashmapExists(constoblock, subscipconss[i][j]) )
         {
            SCIP_CALL( SCIPhashmapInsert(constoblock, subscipconss[i][j], (void *) (size_t) (i+1) ));
         }
      }
   }

   for( i = detectordata->nrepresentatives; i > 0; --i )
   {

      block = (long int)SCIPhashmapGetImage(constoblock, SCIPhashmapGetImage(representatives, (void*)(size_t)i));
      for( j = 0; j < SCIPhashmapGetNLists(mergedconss[i - 1]); ++j )
      {
         list = SCIPhashmapGetList(mergedconss[i - 1], j);
         if( SCIPhashmapListGetNEntries(list) != 0 )
         {
            while( list != NULL )
            {
               if( (SCIP_CONS*)SCIPhashmapGetImage(representatives, (void*)(size_t)i)
                  != (SCIP_CONS*)SCIPhashmapListGetOrigin(list) )
               {
                  subscipconss[block - 1][detectordata->nsubscipconss[block - 1]] = (SCIP_CONS*)SCIPhashmapListGetOrigin(
                     list);
                  detectordata->nsubscipconss[block - 1]++;
                  if( !SCIPhashmapExists(constoblock, (SCIP_CONS*)SCIPhashmapListGetOrigin(list)) )
                  {
                     SCIP_CALL(
                        SCIPhashmapInsert(constoblock,(SCIP_CONS*)SCIPhashmapListGetOrigin(list),(void*) (size_t) block));
                  }
               }

               list = SCIPhashmapListGetNext(list);
            }
         }
      }
   }
   //printf("nconstoblock: %d und nrelconss: %d \n",SCIPhashmapGetNEntries(constoblock),detectordata->nrelconss);
   assert(SCIPhashmapGetNEntries(constoblock) == detectordata->nrelconss);

   j = 0;
   for( i = 0; i < detectordata->nblocks; ++i )
   {
      j += detectordata->nsubscipconss[i];
   }assert(j == detectordata->nrelconss);

   return SCIP_OKAY;
}

/** Builds the transformed problem in the new scip instance */
static SCIP_RETCODE buildTransformedProblem(SCIP* scip, /**< SCIP data structure */
DEC_DETECTORDATA* detectordata /**< presolver data data structure */
)
{
   int i;
   int j;
   int k;
   int block;
   int stop;
   int indexcons;
   int indexvar;
   SCIP_HASHMAP* vartoblock;
   SCIP_HASHMAP* constoblock;
   SCIP_CONS*** subscipconss;
   int* nsubscipconss;
   SCIP_CONS*** varinconss;
   int* nvarinconss;
   SCIP_VAR*** subscipvars;
   int* nsubscipvars;
   SCIP_VAR** linkingvars;
   int nlinkingvars;
   SCIP_VAR** vars;
   int nvars;
   int no;
   //SCIP_VAR** linking;
   //int nlinking;
   SCIP_HASHMAP* linking;
   SCIP_HASHMAPLIST* list;
   int newblock;
   int oldblock;
   int zahl;

   //SCIPallocMemoryArray(scip, &linking, 2*detectordata->nrelvars);
   SCIPhashmapCreate(&linking, SCIPblkmem(scip), detectordata->nrelvars);

   subscipconss = detectordata->subscipconss;
   nsubscipconss = detectordata->nsubscipconss;
   varinconss = detectordata->varinconss;
   nvarinconss = detectordata->nvarinconss;
   nlinkingvars = 0;

   constoblock = detectordata->constoblock;
   vartoblock = detectordata->varstoblock;
   subscipvars = detectordata->subscipvars;
   nsubscipvars = detectordata->nsubscipvars;
   linkingvars = detectordata->linkingvars;

   indexcons = 1;
   indexvar = 1;
   //nlinking = 0;
   /*
    for( i = 0; i < detectordata->nblocks; i++ )
    {
    printf("nsubscipconss: %d \n",nsubscipconss[i] );
    }*/

   newblock = 0;
   oldblock = 0;
   stop = 0;

   zahl = detectordata->nblocks;
   block = detectordata->startblock;
   while( zahl > 0 )
   {
      for( i = 0; i < detectordata->nsubscipconss[block - 1]; ++i )
      {
         SCIP_CALL(
            SCIPhashmapInsert(detectordata->consindex, detectordata->subscipconss[block-1][i], (void*) (size_t)indexcons));
         indexcons++;

         vars = SCIPgetVarsXXX(scip, detectordata->subscipconss[block - 1][i]);
         nvars = SCIPgetNVarsXXX(scip, detectordata->subscipconss[block - 1][i]);

         for( j = 0; j < nvars; j++ )
         {
            if( isVarRelevant(vars[j]) )
            {
               stop = 0;
               no = (long int)SCIPhashmapGetImage(detectordata->vartopos, SCIPvarGetProbvar(vars[j]));
               for( k = 0; (k < detectordata->nvarinconss[no]); ++k )
               {
                  if( (long int)SCIPhashmapGetImage(constoblock, detectordata->varinconss[no][k]) != block )
                  {
                     if( (long int)SCIPhashmapGetImage(constoblock, detectordata->varinconss[no][k]) != oldblock )
                     {
                        newblock = (long int)SCIPhashmapGetImage(constoblock, detectordata->varinconss[no][k]);
                     }
                     stop = 1;
                  }
               }
               if( stop > 0 )
               {
                  //linking[nlinking] = SCIPvarGetProbvar(vars[j]);
                  //nlinking++;
                  SCIPhashmapInsert(linking, SCIPvarGetProbvar(vars[j]), NULL);
               }
               else
               {
                  if( !SCIPhashmapExists(detectordata->varindex, SCIPvarGetProbvar(vars[j])) )
                  {
                     SCIP_CALL(
                        SCIPhashmapInsert(detectordata->varindex, SCIPvarGetProbvar(vars[j]), (void*) (size_t)indexvar));
                     indexvar++;
                  }
               }
            }
         }

         SCIPfreeMemoryArrayNull(scip, &vars);
      }

      oldblock = block;
      block = newblock;
      assert( 0 < block);
      assert( block <= detectordata->nblocks);
      for( k = 0; k < SCIPhashmapGetNLists(linking); k++ )
      {
         list = SCIPhashmapGetList(linking, k);
         while( list != NULL )
         {
            if( !SCIPhashmapExists(detectordata->varindex, SCIPhashmapListGetOrigin(list)) )
            {
               SCIP_CALL(
                  SCIPhashmapInsert(detectordata->varindex, SCIPhashmapListGetOrigin(list), (void*) (size_t)indexvar));
               indexvar++;
            }
            list = SCIPhashmapListGetNext(list);
         }
      }
      //nlinking = 0;
      SCIPhashmapRemoveAll(linking);
      zahl--;

   }

   /*
    for( k = 0; k < nlinking; k++ )
    {
    if( !SCIPhashmapExists(detectordata->varindex, linking[k]) )
    {
    SCIP_CALL( SCIPhashmapInsert(detectordata->varindex, linking[k], (void*) (size_t)indexvar));
    indexvar++;
    }
    }
    */

   /** linkingvars **/

   for( i = 0; i < detectordata->nrelvars; i++ )
   {
      stop = 0;
      block = (long int)SCIPhashmapGetImage(constoblock, varinconss[i][0]);
      for( j = 1; (j < nvarinconss[i]) && (stop == 0); j++ )
      {
         if( block != (long int)SCIPhashmapGetImage(constoblock, varinconss[i][j]) )
         {
            stop = 1;
         }
      }

      if( stop )
      {
         linkingvars[nlinkingvars] = detectordata->relvars[i];
         nlinkingvars++;
      }
      else
      {
         SCIP_CALL(SCIPhashmapInsert(vartoblock, detectordata->relvars[i], (void*) (size_t) block));
         subscipvars[block - 1][nsubscipvars[block - 1]] = detectordata->relvars[i];
         nsubscipvars[block - 1]++;
      }
   }

   if( nlinkingvars == detectordata->nrelvars ) ///hier noch not found
   {
      for( i = 0; i < detectordata->nrelvars; i++ )
      {
         SCIP_CALL(SCIPhashmapInsert(vartoblock, detectordata->relvars[i], (void*) (size_t) 1));
      }
   }

   j = 0;
   for( i = 0; i < detectordata->nblocks; i++ )
   {
      j += nsubscipvars[i];
   }assert(detectordata->nrelvars == (j + nlinkingvars));
   assert(detectordata->nrelvars == SCIPhashmapGetNEntries(detectordata->varindex));

   detectordata->nlinkingvars = nlinkingvars;

   /*
    for( i = 0; i < detectordata->nblocks; i++ )
    {
    printf("nsubscipvars: %d \n",detectordata->nsubscipvars[i] );
    }
    */

   //SCIPfreeMemoryArrayNull(scip, &linking);
   SCIPhashmapFree(&linking);

   return SCIP_OKAY;
}

static DEC_DECL_DETECTSTRUCTURE(detectAndBuildBordered)
{
   int i;
   DEC_DETECTOR* cutpacking;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);

   cutpacking = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(cutpacking);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(cutpacking), DEC_DETECTORNAME) == 0);
   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);

   /* build the hypergraph structure from the original problem */
   SCIP_CALL(buildGraphStructure(scip, detectordata));
   SCIPdebugMessage("buildGraphstructure successful \n");

   /* get the partitions for the new variables from metis */
   while( detectordata->ngraphs > 0 )
   {
      for( i = 0; i < detectordata->nrelconss + 1; i++ )
      {
         if( SCIPhashmapExists(detectordata->occupied, (void*)(size_t)(i + 1)) )
         {
            detectordata->position = i;

            SCIP_CALL(StoerWagner(scip, detectordata));
            SCIPdebugMessage("StoerWagner successful \n");

            SCIP_CALL( buildnewgraphs(scip, detectordata));
            SCIPdebugMessage("buildnewgraphs successful \n");
         }
      }
   }
   /** add merged conss */
   SCIP_CALL(getmergedconss(scip, detectordata));
   SCIPdebugMessage("getmergedconss successful \n");

   /** get subscipvars */
   SCIP_CALL(buildTransformedProblem(scip, detectordata));
   SCIPdebugMessage("buildTransformedProblem successful \n");

   detectordata->found = TRUE;

   /** copy data to decdecomp */
   SCIP_CALL(copyDetectorDataToDecomp(scip, detectordata, detectordata->decdecomp));
   SCIPdebugMessage("copyDetectorDataToDecomp successful \n");

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** set the decomp structure */
static
DEC_DECL_SETSTRUCTDECOMP(CutpackingSetDecomp)
{
   DEC_DETECTOR* cutpacking;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   cutpacking = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(cutpacking);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(cutpacking), DEC_DETECTORNAME) == 0);
   SCIPdebugMessage("Setting decdecomp\n");
   detectordata->decdecomp = decdecomp;

}

static
DEC_DECL_GETPRIORITY(getPriority)
{
   DEC_DETECTOR* arrowheur;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);
   arrowheur = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(arrowheur);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(arrowheur), DEC_DETECTORNAME) == 0);
   return detectordata->priority;
}

/** creates the cutpacking presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionCutpacking(SCIP* scip /**< SCIP data structure */

)
{
   DEC_DETECTORDATA *detectordata;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata));

   assert(detectordata != NULL);
   detectordata->found = FALSE;
   detectordata->partition = NULL;
   detectordata->nblocks = -1;

   SCIP_CALL(
      DECincludeDetector(scip, DEC_DETECTORNAME, detectordata, detectAndBuildBordered, CutpackingSetDecomp, initCutpacking, exitCutpacking, getPriority));

   /* add cutpacking presolver parameters */

   SCIP_CALL(
      SCIPaddIntParam(scip, "cutpacking/priority", "random seed for hmetis", &detectordata->priority, FALSE, DEFAULT_PRIORITY, INT_MIN, INT_MAX, NULL, NULL));

   return SCIP_OKAY;
}
