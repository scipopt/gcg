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

#include "dec_cutpacking.h"

#include "cons_decomp.h"
#include "struct_decomp.h"
#include "scip_misc.h"

#define DEC_DETECTORNAME      "cutpacking"   /**< name of the detector */
#define DEC_PRIORITY          -50              /**< priority of the detector */

/* Default parameter settings */
#define DEFAULT_RANDSEED                  1     /**< random seed for the hmetis call */
#define DEFAULT_TIDY                      TRUE  /**< whether to clean up afterwards */

#define DEFAULT_METIS_UBFACTOR            5.0   /**< default unbalance factor given to metis on the commandline */
#define DEFAULT_METIS_VERBOSE             FALSE /**< should metis be verbose */
#define DEFAULT_METISUSEPTYPE_RB          TRUE  /**< Should metis use the rb or kway partitioning algorithm */
#define DEFAULT_PRIORITY                  DEC_PRIORITY

//#define DWSOLVER_REFNAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_ref.txt", (name), (blocks), (cons), (dummy)

//#define GP_NAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_%d.gp", (name), (blocks), (cons), (dummy)

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

   int nblocks;
   SCIP_CONS*** subscipconss;
   int* nsubscipconss;
   SCIP_VAR*** subscipvars;
   int* nsubscipvars;
   SCIP_VAR** linkingvars;
   int nlinkingvars;

   SCIP_HASHMAP* constoblock;
   SCIP_HASHMAP* varstoblock;

   Graph** graphs;
   int ngraphs;

   SCIP_Bool delete;
   int position;
   int* partition;

   SCIP_HASHMAP** mergedconss;
   SCIP_HASHMAP* representatives;
   int nrepresentatives;

   SCIP_HASHMAP* vartopos;
   int* nvarinconss;
   SCIP_CONS*** varinconss;

   /* Graph stuff for hmetis */
   int randomseed;
   SCIP_Real metisubfactor;
   SCIP_Bool metisverbose;
   SCIP_Bool metisuseptyperb;
   SCIP_Bool found;
   SCIP_Bool tidy;

   int priority;
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static DEC_DECL_INITDETECTOR(initCutpacking)
{

   int i;
   int j;
   int nallvars;
   int nconss;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_VAR** allvars;
   int nvars;
   SCIP_HASHMAP* vartopos;

   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   nallvars = SCIPgetNVars(scip);
   allvars = SCIPgetVars(scip);
   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);

   detectordata->nblocks = 0;
   detectordata->ngraphs = 0;
   detectordata->nlinkingvars = 0;
   detectordata->position = -1;
   detectordata->nrepresentatives = 0;

   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->partition, nconss));
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->graphs, nconss));
   //reallocate
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->graphs[i]->conss, nconss));
      SCIP_CALL(SCIPhashmapCreate(&detectordata->graphs[i]->constopos, SCIPblkmem(scip),nconss));
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->graphs[i]->adjacencylist, nconss));
      for( j = 0; j < nconss; j++ )
      {
         SCIP_CALL(SCIPhashmapCreate(&detectordata->graphs[i]->adjacencylist[j], SCIPblkmem(scip),nconss));
      }
   }

   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->nvarinconss, nallvars));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->vartopos, SCIPblkmem(scip),nallvars));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varinconss, nallvars));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->subscipvars, nallvars));

   for( i = 0; i < nallvars; i++ )
   {
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varinconss[i], nconss));
   }

   SCIP_CALL(SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip),nconss));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->varstoblock, SCIPblkmem(scip),nallvars));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->representatives, SCIPblkmem(scip),nconss));

   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->subscipconss, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->mergedconss, nconss));
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->subscipconss[i], nconss));
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->subscipvars[i], nallvars));
      SCIP_CALL(SCIPhashmapCreate(&detectordata->mergedconss[i], SCIPblkmem(scip),nconss));
   }

   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->nsubscipconss, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->subscipvars, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->linkingvars, nallvars));

   /** vartopos */

   vartopos = detectordata->vartopos;

   for( i = 0; i < nallvars; i++ )
   {
      SCIP_CALL( SCIPhashmapInsert(vartopos, allvars[i], (void*) (size_t) i) );
   }

   /** varinconss */

   for( i = 0; i < nconss; i++ )
   {
      vars = SCIPgetVarsXXX(scip, conss[i]);
      nvars = SCIPgetNVarsXXX(scip, conss[i]);

      for( j = 0; j < nvars; j++ )
      {
         (detectordata->varinconss[(long int)SCIPhashmapGetImage(vartopos, vars[i])])[detectordata->nvarinconss[(long int)SCIPhashmapGetImage(vartopos, vars[i])]] = conss[i];
         detectordata->nvarinconss[(long int)SCIPhashmapGetImage(vartopos, vars[i])]++;
      }
   }

   return SCIP_OKAY;
}

/** copies the variable and block information to the decomp structure */
static SCIP_RETCODE copyDetectorDataToDecomp(
   SCIP* scip,                         /**< SCIP data structure */
   DEC_DETECTORDATA* detectordata,     /**< presolver data data structure */
   DECDECOMP* decomp                   /**< DECOMP data structure */
)
{
   int i;
   assert(scip != 0);
   assert(detectordata != 0);
   assert(decomp != 0);

   SCIP_CALL(SCIPallocMemoryArray(scip, &decomp->subscipvars, detectordata->nblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &decomp->subscipconss, detectordata->nblocks));

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->linkingvars, detectordata->linkingvars, detectordata->nlinkingvars));
   decomp->nlinkingvars = detectordata->nlinkingvars;

   for( i = 0; i < detectordata->nblocks; i++ )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decomp->subscipconss[i], detectordata->subscipconss[i], detectordata->nsubscipconss[i]) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decomp->subscipvars[i], detectordata->subscipvars[i], detectordata->nsubscipvars[i]) );
   }

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipconss, detectordata->nsubscipconss, detectordata->nblocks));
   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipvars, detectordata->nsubscipvars, detectordata->nblocks));

   decomp->constoblock = detectordata->constoblock;
   decomp->vartoblock = detectordata->varstoblock;
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

   for( i = 0; i < SCIPgetNConss(scip); i++ )
   {
      SCIPfreeMemoryArray(scip, &detectordata->subscipconss[i]);
      SCIPfreeMemoryArray(scip, &detectordata->subscipvars[i]);
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

   SCIPhashmapFree(&detectordata->vartopos);
   SCIPhashmapFree(&detectordata->representatives);

   for( i = 0; i < detectordata->nrepresentatives; i++ )
   {
      SCIPhashmapFree(&detectordata->mergedconss[i]);
   }SCIPfreeMemoryArray(scip, &detectordata->mergedconss);

   //SCIPhashmapFree(&detectordata->varstoblock);
   //SCIPhashmapFree(&detectordata->consstoblock);

   return SCIP_OKAY;
}

/**
 * Builds a graph structure out of the matrix.
 *
 * The function will create a vertice for every constraint and an edge between two constraints if they have/own a common variable....
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

   SCIP_HASHMAP* constopos;
   SCIP_HASHMAP** adjacencylist;
   Graph* graph;

   graph = detectordata->graphs[0];
   constopos = detectordata->graphs[0]->constopos;
   adjacencylist = detectordata->graphs[0]->adjacencylist;

   nedges = 0;

   /** constopos */

   for( i = 0; i < SCIPgetNConss(scip); i++ )
   {
      SCIP_CALL(SCIPhashmapInsert(constopos, SCIPgetConss(scip)[i], (void*) (size_t) i));
   }

   /** adjacencylist */

   varinconss = detectordata->varinconss;

   for( i = 0; i < SCIPgetNVars(scip); i++ )
   {
      for( j = 0; j < detectordata->nvarinconss[i]; j++ )
      {
         for( k = j; k < detectordata->nvarinconss[i]; k++ )
         {
            if( SCIPhashmapExists(adjacencylist[(long int)SCIPhashmapGetImage(constopos, varinconss[i][j])],
               varinconss[i][k]) )
            {
               cost = (long int)SCIPhashmapGetImage(
                  adjacencylist[(long int)SCIPhashmapGetImage(constopos, varinconss[i][j])], varinconss[i][k]);
               cost++;
               SCIP_CALL(
                  SCIPhashmapSetImage(adjacencylist[(long int)SCIPhashmapGetImage(constopos, varinconss[i][j])], varinconss[i][k], (void*) (size_t) cost));
               cost = (long int)SCIPhashmapGetImage(
                  adjacencylist[(long int)SCIPhashmapGetImage(constopos, varinconss[i][k])], varinconss[i][j]);
               cost++;
               SCIP_CALL(
                  SCIPhashmapSetImage(adjacencylist[(long int)SCIPhashmapGetImage(constopos, varinconss[i][k])], varinconss[i][j], (void*) (size_t) cost));
            }
            else
            {
               SCIP_CALL(
                  SCIPhashmapInsert(adjacencylist[(long int)SCIPhashmapGetImage(constopos, varinconss[i][j])], varinconss[i][k],(void*) (size_t) 1));
               SCIP_CALL(
                  SCIPhashmapInsert(adjacencylist[(long int)SCIPhashmapGetImage(constopos, varinconss[i][k])], varinconss[i][j],(void*) (size_t) 1));
               nedges++;
            }
         }
      }
   }

   graph->conss = SCIPgetConss(scip);
   graph->nconss = SCIPgetNConss(scip);
   graph->cons1 = NULL;
   graph->cons2 = NULL;
   graph->constopos = constopos;
   graph->adjacencylist = adjacencylist;
   graph->nedges = nedges;

   //detectordata->graphs[0] = graph;
   detectordata->ngraphs = 1;

   return SCIP_OKAY;
}

/** Will call hmetis via a system call */
static SCIP_RETCODE callMetis(SCIP* scip, /**< SCIP data struture */
DEC_DETECTORDATA* detectordata, /**< presolver data data structure */
SCIP_RESULT *result)
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
   void* entry;
   int* partition;
   SCIP_CONS** conss;
   SCIP_HASHMAP** adja;
   SCIP_HASHMAP* constopos;
   SCIP_HASHMAPLIST* list;

   nvertices = detectordata->graphs[detectordata->position]->nconss;
   nedges = detectordata->graphs[detectordata->position]->nedges;

   SCIP_FILE *zfile;
   FILE* file;
   int temp_filedes = -1;
   //SCIP_Real remainingtime;

   assert(scip != NULL);
   assert(detectordata != NULL);

   *result = SCIP_DIDNOTRUN;

   adja = detectordata->graphs[detectordata->position]->adjacencylist;
   conss = detectordata->graphs[detectordata->position]->conss;
   constopos = detectordata->graphs[detectordata->position]->constopos;

   //remainingtime = DECgetRemainingTime(scip);

   //if( remainingtime <= 0 )
   //{
   // return SCIP_OKAY;
   //}

   SCIPsnprintf(tempfile, SCIP_MAXSTRLEN, "gcg-metis-XXXXXX");
   if( (temp_filedes = mkstemp(tempfile)) < 0 )
   {
      SCIPerrorMessage("Error creating temporary file: %s\n", strerror(errno));
      return SCIP_FILECREATEERROR;
   }

   SCIPdebugMessage("Temporary filename: %s\n", tempfile);

   file = fdopen(temp_filedes, "w");
   if( file == NULL)
   {
      SCIPerrorMessage("Could not open temporary metis file!\n");
      return SCIP_FILECREATEERROR;
   }

   SCIPinfoMessage(scip, file, "%d %d 001\n", nvertices, nedges);
   i = 0;
   for( i = 0; i < nvertices; i++ )
   {
      for( j = 0; j < SCIPhashmapGetNLists(adja[i]); ++j )
      {
         list = SCIPhashmapGetList(adja[i], j);

         if( SCIPhashmapListGetNEntries(list) != 0 )
         {
            while( list != NULL )
            {
               entry = SCIPhashmapGetImage(constopos, SCIPhashmapListGetOrigin(list));
               SCIPinfoMessage(scip, file, "%d ", (long int)entry);
               entry = SCIPhashmapListGetImage(list);
               SCIPinfoMessage(scip, file, "%d ", (long int)entry);
               list = SCIPhashmapListGetNext(list);
            }
         }
      }

      SCIPinfoMessage(scip, file, "\n");

   }
   status = fclose(file);

   if( status == -1 )
   {
      SCIPerrorMessage("Could not close '%s'\n", tempfile);
      return SCIP_WRITEERROR;
   }

   /* call metis via syscall as there is no library usable ... */
   /*if(!SCIPisInfinity(scip, remainingtime))
    {
    SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "timeout %.0fs ./hmetis %s %d -seed %d -ptype %s -ufactor %f %s",
    remainingtime,
    tempfile,
    2,
    detectordata->randomseed,
    detectordata->metisuseptyperb ? "rb" : "kway",
    detectordata->metisubfactor,
    detectordata->metisverbose ? "" : "> /dev/null" );
    }
    else
    {*/
   SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "./hmetis %s %d -seed %d -ptype %s -ufactor %f %s", tempfile, 2,
      detectordata->randomseed, detectordata->metisuseptyperb ? "rb" : "kway", detectordata->metisubfactor,
      detectordata->metisverbose ? "" : "> /dev/null");
   //}

   //SCIP_CALL(SCIPresetClock(scip, detectordata->metisclock));
   //SCIP_CALL(SCIPstartClock(scip, detectordata->metisclock));
   //SCIPdebugMessage("Calling metis with: %s\n", metiscall);

   status = system(metiscall);

   //SCIP_CALL(SCIPstopClock(scip, detectordata->metisclock));
   //SCIPdebugMessage("time left before metis started: %f, time metis spend %f, remainingtime: %f\n", remainingtime, SCIPgetClockTime(scip, detectordata->metisclock),  remainingtime-SCIPgetClockTime(scip, detectordata->metisclock) );

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
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->partition, nvertices));
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
      assert(temp >= 0 && temp <= 1);
      partition[i] = temp;
      i++;
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

static SCIP_RETCODE buildnewgraphs(SCIP* scip, DEC_DETECTORDATA* detectordata)
{
   int i;
   int j;
   int cas;
   int cost;
   int nedges;

   int* partition;
   Graph* graph;
   SCIP_CONS* representative;

   SCIP_HASHMAP* adja;
   SCIP_HASHMAP* newadja1; // reicht eins?
   SCIP_HASHMAP* newadja2;
   SCIP_HASHMAPLIST* list;

   Graph* graph1;
   SCIP_CONS** conss1;
   int nconss1;
   SCIP_HASHMAP* constopos1;
   SCIP_HASHMAP** adjacencylist1;
   SCIP_HASHMAP* consslink1;

   Graph* graph2;
   SCIP_CONS** conss2;
   int nconss2;
   SCIP_HASHMAP* constopos2;
   SCIP_HASHMAP** adjacencylist2;
   SCIP_HASHMAP* consslink2;

   cas = -1;
   nconss1 = 0;
   nconss2 = 0;

   /** build partitions */
   partition = detectordata->partition;
   graph = detectordata->graphs[detectordata->position];

   graph1 = detectordata->graphs[detectordata->position + 1];
   conss1 = graph1->conss;
   constopos1 = graph1->constopos;
   adjacencylist1 = graph1->adjacencylist;

   graph2 = detectordata->graphs[detectordata->position + 2]; //koennte ein graph zu viel werden, falls jede cons ein eigener block ist
   conss2 = graph2->conss;
   constopos2 = graph2->constopos;
   adjacencylist2 = graph2->adjacencylist;

   for( i = 0; i < graph->nconss; i++ )
   {
      assert((-1 < partition[i])&&(partition[i]<2));
      if( partition[i] == 0 )
      {
         SCIP_CALL(SCIPhashmapInsert(constopos1, graph->conss[i], NULL));
         conss1[nconss1] = graph->conss[i];
         nconss1++;
      }
      else
      {
         SCIP_CALL(SCIPhashmapInsert(constopos2, graph->conss[i], NULL));
         conss2[nconss2] = graph->conss[i];
         nconss2++;
      }
   }

   /** cut feasible? */

   if( (graph->cons1 != NULL) && (graph->cons2 != NULL) )
   {

      if( (SCIPhashmapExists(constopos1, graph->cons1) && SCIPhashmapExists(constopos1, graph->cons2))
         || (SCIPhashmapExists(constopos2, graph->cons1) && SCIPhashmapExists(constopos2, graph->cons2)) )
      {
         detectordata->subscipconss[detectordata->nblocks] = graph->conss;
         detectordata->nsubscipconss[detectordata->nblocks] = graph->nconss;
         detectordata->nblocks++;
         detectordata->ngraphs--;
         for( i = detectordata->position; i < detectordata->ngraphs; i++ )
         {
            detectordata->graphs[i] = detectordata->graphs[i + 1];
         }
         detectordata->delete = TRUE;
      }
      else
      {
         if( SCIPhashmapExists(constopos1, graph->cons1) )
         {
            cas = 0;
         }
         else
         {
            cas = 1;
         }
      }
   }
   else if( graph->cons1 != NULL)
   {
      if( SCIPhashmapExists(constopos1, graph->cons1) )
      {
         cas = 0;
      }
      else
      {
         cas = 1;
      }
   }
   else
   {
      if( SCIPhashmapExists(constopos2, graph->cons2) )
      {
         cas = 0;
      }
      else
      {
         cas = 1;
      }
   }

   /** get contraints to merge */

   if( cas != -1 )
   {

      SCIP_CALL(SCIPhashmapCreate(&consslink1, SCIPblkmem(scip), graph->nconss));
      SCIP_CALL(SCIPhashmapCreate(&consslink2, SCIPblkmem(scip), graph->nconss));
      SCIP_CALL(SCIPhashmapCreate(&adja, SCIPblkmem(scip), graph->nconss));
      SCIP_CALL(SCIPhashmapCreate(&newadja1, SCIPblkmem(scip), graph->nconss));
      SCIP_CALL(SCIPhashmapCreate(&newadja2, SCIPblkmem(scip), graph->nconss));

      for( i = 0; i < nconss1; ++i )
      {
         adja = graph->adjacencylist[(long int)SCIPhashmapGetImage(graph->constopos, conss1[i])];

         for( i = 0; i < SCIPhashmapGetNLists(adja); ++i )
         {
            list = SCIPhashmapGetList(adja, i);
            if( SCIPhashmapListGetNEntries(list) != 0 )
            {
               while( list != NULL )
               {
                  if( SCIPhashmapExists(constopos2, SCIPhashmapListGetOrigin(list)) )
                  {
                     if( !SCIPhashmapExists(consslink1, (void*)conss1[i]) )
                     {
                        SCIP_CALL( SCIPhashmapInsert(consslink1, conss1[i], NULL));
                     }
                     if( !SCIPhashmapExists(consslink2, SCIPhashmapListGetOrigin(list)) )
                     {
                        SCIP_CALL( SCIPhashmapInsert(consslink2, SCIPhashmapListGetOrigin(list), NULL));
                     }
                  }
                  list = SCIPhashmapListGetNext(list);
               }
            }
         }
      }
      /*********************/
      /** build new graphs */
      /*********************/

      /** graph1 */
      /** adjacencylist */
      if( nconss1 > 0 )
      {
         representative = NULL;
         j = 0;
         nedges = 0;
         for( i = 0; i < nconss1; i++ )
         {
            if( !SCIPhashmapExists(consslink1, conss1[i]) )
            {
               adjacencylist1[i] = graph->adjacencylist[(long int)SCIPhashmapGetImage(graph->constopos, conss1[i])];
               SCIP_CALL(SCIPhashmapSetImage(constopos1, conss1[i], (void*) (size_t) j));
               j++;
            }
            else
            {
               //SCIP_CALL(SCIPhashmapRemoveAll(adja));
               representative = conss1[i];
               SCIP_CALL(SCIPhashmapRemove(constopos1, conss1[i]));
               for( i = 0;
                  i < SCIPhashmapGetNLists(graph->adjacencylist[(long int)SCIPhashmapGetImage(graph->constopos, conss1[i])]);
                  ++i )
                  {
                  list = SCIPhashmapGetList(
                     graph->adjacencylist[(long int)SCIPhashmapGetImage(graph->constopos, conss1[i])], i);
                  if( SCIPhashmapListGetNEntries(list) != 0 )
                  {
                     while( list != NULL )
                     {
                        if( !SCIPhashmapExists(consslink1, SCIPhashmapListGetOrigin(list)) )
                        {
                           if( SCIPhashmapExists(adja, SCIPhashmapListGetOrigin(list)) )
                           {
                              cost = (long int)SCIPhashmapGetImage(newadja1, SCIPhashmapListGetOrigin(list));
                              cost += (long int)SCIPhashmapListGetImage(list);
                              SCIP_CALL(
                                 SCIPhashmapSetImage(newadja1, SCIPhashmapListGetOrigin(list), (void*) (size_t) cost));
                           }
                           else
                           {
                              SCIP_CALL(
                                 SCIPhashmapInsert(newadja1, SCIPhashmapListGetOrigin(list), SCIPhashmapListGetImage(list)));
                           }
                        }
                        list = SCIPhashmapListGetNext(list);
                     }
                  }
               }
            }
         }

         adjacencylist1[j] = newadja1;
         SCIP_CALL(SCIPhashmapInsert(constopos1, representative,(void*) (size_t) j));
         nconss1 = j + 1;
         /** delete merged conss */
         for( i = 0; i < nconss1 - 1; i++ )
         {
            for( j = 0; j < SCIPhashmapGetNLists(adjacencylist1[i]); j++ )
            {
               list = SCIPhashmapGetList(adjacencylist1[i], j);
               cost = 0;
               if( SCIPhashmapListGetNEntries(list) != 0 )
               {
                  while( list != NULL )
                  {
                     if( SCIPhashmapExists(consslink1, SCIPhashmapListGetOrigin(list)) )
                     {
                        representative = (SCIP_CONS*)SCIPhashmapListGetOrigin(list);
                        cost += (long int)SCIPhashmapListGetImage(list);
                        SCIP_CALL(SCIPhashmapRemove(adjacencylist1[i],SCIPhashmapListGetOrigin(list)));
                     }
                     else
                     {
                        nedges++;
                     }
                     list = SCIPhashmapListGetNext(list);
                  }

               }
               if( cost != 0 )
               {
                  SCIP_CALL(SCIPhashmapInsert(adjacencylist1[i],representative, (void*) (size_t) cost));
                  nedges += 2;
               }
            }
         }

         if( representative != NULL)
         {
            SCIP_CALL(
               SCIPhashmapInsert(detectordata->representatives, representative, (void*) (size_t) detectordata->nrepresentatives));
            detectordata->mergedconss[detectordata->nrepresentatives] = consslink1;
            ;
            detectordata->nrepresentatives++;
         }

         /** arranges conss1 */
         for( i = 0; i < SCIPhashmapGetNLists(constopos1); ++i )
         {
            list = SCIPhashmapGetList(constopos1, i);
            if( SCIPhashmapListGetNEntries(list) != 0 )
            {
               while( list != NULL )
               {
                  conss1[(long int)SCIPhashmapListGetImage(list)] = (SCIP_CONS*)SCIPhashmapListGetOrigin(list);
                  list = SCIPhashmapListGetNext(list);
               }
            }
         }

         /** builds graph1 */
         //graph1->adjacencylist = adjacencylist1;
         //graph1->constopos = constopos1;
         //graph1->conss = conss1;
         graph1->nconss = nconss1;
         graph1->nedges = nedges / 2;
         assert( 2*graph1->nedges == nedges);

         switch( cas )
         {
         case 0:
            graph1->cons1 = graph->cons1;
            graph1->cons2 = representative;
            break;

         case 1:
            graph1->cons1 = representative;
            graph1->cons2 = graph->cons2;
            break;

         default:
            break;
         }

      }
      /** graph 2 */
      /** adjacencylist */
      if( nconss2 > 1 )
      {
         representative = NULL;
         j = 0;
         nedges = 0;

         for( i = 0; i < nconss2; i++ )
         {
            if( !SCIPhashmapExists(consslink2, conss2[i]) )
            {
               adjacencylist1[i] = graph->adjacencylist[(long int)SCIPhashmapGetImage(graph->constopos, conss2[i])];
               SCIP_CALL(SCIPhashmapSetImage(constopos2, conss2[i], (void*) (size_t) j));
               j++;
            }
            else
            {
               representative = conss2[i];
               SCIP_CALL(SCIPhashmapRemove(constopos2, conss2[i]));
               for( i = 0;
                  i < SCIPhashmapGetNLists(graph->adjacencylist[(long int)SCIPhashmapGetImage(graph->constopos, conss2[i])]);
                  ++i )
                  {
                  list = SCIPhashmapGetList(
                     graph->adjacencylist[(long int)SCIPhashmapGetImage(graph->constopos, conss2[i])], i);
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
                           else
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

         adjacencylist1[j] = newadja2;
         SCIP_CALL(SCIPhashmapInsert(constopos2, representative,(void*) (size_t) j));
         nconss1 = j + 1;

         /** delete merged conss */
         for( i = 0; i < nconss2 - 1; i++ )
         {
            for( j = 0; j < SCIPhashmapGetNLists(adjacencylist2[i]); j++ )
            {
               list = SCIPhashmapGetList(adjacencylist2[i], j);
               cost = 0;
               if( SCIPhashmapListGetNEntries(list) != 0 )
               {
                  while( list != NULL )
                  {
                     if( SCIPhashmapExists(consslink2, SCIPhashmapListGetOrigin(list)) )
                     {
                        cost += (long int)SCIPhashmapListGetImage(list);
                        SCIP_CALL(SCIPhashmapRemove(adjacencylist2[i],SCIPhashmapListGetOrigin(list)));
                     }
                     else
                     {
                        nedges++;
                     }
                     list = SCIPhashmapListGetNext(list);
                  }

               }
               if( cost > 0 )
               {
                  SCIP_CALL(SCIPhashmapInsert(adjacencylist2[i],representative, (void*) (size_t) cost));
                  nedges += 2;
               }
            }
         }

         SCIP_CALL(
            SCIPhashmapInsert(detectordata->representatives, representative, (void*) (size_t) detectordata->nrepresentatives));
         detectordata->mergedconss[detectordata->nrepresentatives] = consslink1;
         ;
         detectordata->nrepresentatives++;

         /** arranges conss2 */
         for( i = 0; i < SCIPhashmapGetNLists(constopos2); ++i )
         {
            list = SCIPhashmapGetList(constopos2, i);
            if( SCIPhashmapListGetNEntries(list) != 0 )
            {
               while( list != NULL )
               {
                  conss2[(long int)SCIPhashmapListGetImage(list)] = (SCIP_CONS*)SCIPhashmapListGetOrigin(list);
                  list = SCIPhashmapListGetNext(list);
               }
            }
         }

         /** builds graph 2 */
         //graph2->adjacencylist = adjacencylist2;
         //graph2->constopos = constopos2;
         //graph2->conss = conss2;
         graph2->nconss = nconss2;
         graph2->nedges = nedges / 2;
         assert( 2 * graph2->nedges == nedges);

         switch( cas )
         {
         case 0:
            graph2->cons1 = representative;
            graph2->cons2 = graph->cons2;
            break;

         case 1:
            graph2->cons1 = graph->cons1;
            graph2->cons2 = representative;
            break;

         default:
            break;
         }

      }

      if( (nconss1 < 2) && (nconss2 < 2) )
      {
         detectordata->subscipconss[detectordata->nblocks] = conss1;
         detectordata->nsubscipconss[detectordata->nblocks] = nconss1;
         detectordata->nblocks++;
         detectordata->subscipconss[detectordata->nblocks] = conss2;
         detectordata->nsubscipconss[detectordata->nblocks] = nconss2;
         detectordata->nblocks++;
         detectordata->ngraphs--;

         for( i = detectordata->position; i < detectordata->ngraphs; i++ )
         {
            detectordata->graphs[i] = detectordata->graphs[i + 1];
         }

         detectordata->delete = TRUE;
      }
      else if( nconss1 < 2 )
      {
         detectordata->subscipconss[detectordata->nblocks] = conss1;
         detectordata->nsubscipconss[detectordata->nblocks] = nconss1;
         detectordata->nblocks++;

         detectordata->graphs[detectordata->position] = graph2;
      }
      else if( nconss2 < 2 )
      {
         detectordata->subscipconss[detectordata->nblocks] = conss2;
         detectordata->nsubscipconss[detectordata->nblocks] = nconss2;
         detectordata->nblocks++;

         detectordata->graphs[detectordata->position] = graph1;
      }
      else
      {
         detectordata->ngraphs++;
         detectordata->graphs[detectordata->position] = detectordata->graphs[detectordata->ngraphs];
      }

      SCIPhashmapFree(&consslink1);
      SCIPhashmapFree(&consslink2);
      SCIPhashmapFree(&adja);
      SCIPhashmapFree(&newadja1);
      SCIPhashmapFree(&newadja2);
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE getmergedconss(SCIP* scip, DEC_DETECTORDATA* detectordata)
{
   int i;
   int j;
   int k;
   SCIP_HASHMAP** mergedconss;
   SCIP_HASHMAP* representatives;
   SCIP_CONS*** subscipconss;
   SCIP_HASHMAPLIST* list;
   int no;

   mergedconss = detectordata->mergedconss;
   representatives = detectordata->representatives;
   subscipconss = detectordata->subscipconss;

   // �berglegen, wie viele representatives pro block m�glich

   for( i = 0; i < detectordata->nblocks; i++ )
   {
      for( j = 0; j < detectordata->nsubscipconss[i]; j++ )
      {
         if( SCIPhashmapExists(representatives, subscipconss[i][j]) )
         {
            no = (long int)SCIPhashmapGetImage(representatives, subscipconss[i][j]);

            for( k = 0; k < SCIPhashmapGetNLists(mergedconss[no]); ++k )
            {
               list = SCIPhashmapGetList(mergedconss[no], k);
               if( SCIPhashmapListGetNEntries(list) != 0 )
               {
                  while( list != NULL )
                  {
                     if( subscipconss[i][j] != (SCIP_CONS*)SCIPhashmapListGetOrigin(list) )
                     {
                        subscipconss[i][detectordata->nsubscipconss[i]] = (SCIP_CONS*)SCIPhashmapListGetOrigin(list);
                        detectordata->nsubscipconss[i]++;
                     }

                     list = SCIPhashmapListGetNext(list);
                  }
               }
            }
         }
      }
   }

   //detectordata->subscipconss = subscipconss;

   return SCIP_OKAY;
}

/** Builds the transformed problem in the new scip instance */
static SCIP_RETCODE buildTransformedProblem(SCIP* scip, /**< SCIP data structure */
DEC_DETECTORDATA* detectordata /**< presolver data data structure */
)
{
   int i;
   int j;
   long int block;
   int stop;
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

   subscipconss = detectordata->subscipconss;
   nsubscipconss = detectordata->nsubscipconss;
   varinconss = detectordata->varinconss;
   nvarinconss = detectordata->nvarinconss;

   constoblock = detectordata->constoblock;
   vartoblock = detectordata->varstoblock;
   subscipvars = detectordata->subscipvars;
   nsubscipvars = detectordata->nsubscipvars;
   linkingvars = detectordata->linkingvars;

   for( i = 0; i < detectordata->nblocks; i++ )
   {
      nsubscipvars[i] = 0;
      for( j = 0; j < nsubscipconss[i]; j++ )
      {
         SCIP_CALL(SCIPhashmapInsert(constoblock, subscipconss[i][j], (void*) (size_t) i));
      }
   }

   for( i = 0; i < SCIPgetNVars(scip); i++ )
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
         linkingvars[nlinkingvars] = SCIPgetVars(scip)[i];
         nlinkingvars++;
      }
      else
      {
         SCIP_CALL(SCIPhashmapInsert(vartoblock, SCIPgetVars(scip)[i], (void*) (size_t) block));
         subscipvars[block][nsubscipvars[block]] = SCIPgetVars(scip)[i];
         nsubscipvars[block]++;
      }
   }

   //detectordata->constoblock = constoblock;
   //detectordata->vartoblock = vartoblock;
   //detectordata->subscipvars = subscipvars;
   //detectordata->nscubscipvars = nsubscipvars;
   //detectordata->linkingvars = linkingvars;
   detectordata->nlinkingvars = nlinkingvars;

   return SCIP_OKAY;
}

static DEC_DECL_DETECTSTRUCTURE(detectAndBuildBordered)
{
   int i;
   int no;
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

   /* get the partitions for the new variables from metis */
   while( detectordata->ngraphs > 0 )
   {
      no = detectordata->ngraphs;
      for( i = 0; i < no; i++ )
      {
         detectordata->delete = FALSE;
         detectordata->position = i;

         SCIP_CALL(callMetis(scip, detectordata, result));

         if( *result != SCIP_SUCCESS )
         {
            *result = SCIP_DIDNOTFIND;
            return SCIP_OKAY;
         }

         SCIP_CALL( buildnewgraphs(scip, detectordata));

         if(detectordata->delete)
         {
            i--;
            no--;
         }
      }
   }
   /** add merged conss */
   SCIP_CALL(getmergedconss(scip, detectordata));

   /** get subscipvars */
   SCIP_CALL(buildTransformedProblem(scip, detectordata));
   // SCIP_CALL(evaluateDecomposition(scip, detectordata, &score));

   detectordata->found = TRUE;

   /** copy data to decdecomp */
   SCIP_CALL(copyDetectorDataToDecomp(scip, detectordata, detectordata->decdecomp));

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

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );

   assert(detectordata != NULL);
   detectordata->found = FALSE;
   detectordata->partition = NULL;
   detectordata->nblocks = -1;

   SCIP_CALL(
      DECincludeDetector(scip, DEC_DETECTORNAME, detectordata, detectAndBuildBordered, CutpackingSetDecomp, initCutpacking, exitCutpacking, getPriority));

   /* add cutpacking presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "cutpacking/tidy", "Whether to clean up temporary files", &detectordata->tidy, FALSE, DEFAULT_TIDY, NULL, NULL));
   SCIP_CALL( SCIPaddIntParam(scip, "cutpacking/randomseed", "random seed for hmetis", &detectordata->randomseed, FALSE, DEFAULT_RANDSEED, -1, INT_MAX, NULL, NULL));
   SCIP_CALL( SCIPaddRealParam(scip, "cutpacking/ubfactor", "Unbalance factor for metis", &detectordata->metisubfactor, FALSE, DEFAULT_METIS_UBFACTOR, 0.0, 1E20, NULL, NULL ));
   SCIP_CALL( SCIPaddBoolParam(scip, "cutpacking/metisverbose", "Should the metis output be displayed", &detectordata->metisverbose, FALSE, DEFAULT_METIS_VERBOSE, NULL, NULL ));
   SCIP_CALL( SCIPaddBoolParam(scip, "cutpacking/metisuseptyperb", "Should the rb or kway method be used for partitioning by metis", &detectordata->metisuseptyperb, FALSE, DEFAULT_METISUSEPTYPE_RB, NULL, NULL));
   SCIP_CALL( SCIPaddIntParam(scip, "cutpacking/priority", "random seed for hmetis", &detectordata->priority, FALSE, DEFAULT_PRIORITY, INT_MIN, INT_MAX, NULL, NULL));

   return SCIP_OKAY;
}
