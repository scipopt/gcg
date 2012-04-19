/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_borderheur.c
 * @brief  borderheur detector
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#include "dec_borderheur.h"

#include "cons_decomp.h"
#include "struct_decomp.h"
#include "pub_decomp.h"
#include "scip_misc.h"
#include "scip/pub_misc.h"

#define DEC_DETECTORNAME         "borderheur"   /**< name of the detector */
#define DEC_PRIORITY             0              /**< priority of the detector */
#define DEC_DECCHAR              'b'            /**< display character of detector */
#define DEC_ENABLED              TRUE           /**< should detector be called by default */

/* Default parameter settings */
#define DEFAULT_CONSWEIGHT       5              /**< weight for constraint hyperedges */
#define DEFAULT_RANDSEED         1              /**< random seed for the hmetis call */
#define DEFAULT_TIDY             TRUE           /**< whether to clean up afterwards */
#define DEFAULT_DUMMYNODES       0.2            /**< percentage of dummy vertices*/

#define DEFAULT_MAXBLOCKS        20             /**< value for the maximum number of blocks to be considered */
#define DEFAULT_MINBLOCKS        2              /**< value for the minimum number of blocks to be considered */

#define DEFAULT_METIS_UBFACTOR   5.0            /**< default unbalance factor given to metis on the commandline */
#define DEFAULT_METIS_VERBOSE    FALSE          /**< should metis be verbose */
#define DEFAULT_METISUSEPTYPE_RB TRUE           /**< should metis use the rb or kway partitioning algorithm */

/*
 * Data structures
 */

/** hypeedge data structure for hmetis */
struct hyperedge
{
   SCIP_CONS* cons;           /**< the original pointer to the constraint */
   int        cost;           /**< cost of the hyperedge */

};
typedef struct hyperedge HyperEdge;

/** detector data */
struct DEC_DetectorData
{
   /* Graph stuff for hmetis */
   HyperEdge* hedges;         /**< array of hyperedges */
   int*       partition;      /**< array storing vertex partitions */
   int        nvertices;      /**< number of vertices */
   int        nhyperedges;    /**< number of hyperedges */
   int*       varpart;        /**< array storing variable partition */

   /* general parameters */
   SCIP_Bool tidy;            /**< whether temporary metis files should be cleaned up */
   int       maxblocks;       /**< maximal number of blocks to test */
   int       minblocks;       /**< minimal number of blocks to test */
   int       consWeight;      /**< weight of a constraint hyperedge */

   SCIP_Real dummynodes;      /**< percentage of dummy nodes */

   /* metis parameters */
   int       randomseed;      /**< metis random seed */
   SCIP_Real metisubfactor;   /**< metis unbalance factor*/
   SCIP_Bool metisverbose;    /**< shoud the metis out be displayed */
   SCIP_Bool metisuseptyperb; /**< flag to indicate whether metis uses kway or rb partitioning */

   /* various data */
   SCIP_CLOCK* metisclock;    /**< clock to measure metis time */
   int         blocks;        /**< indicates the current block */
   SCIP_Bool   found;         /**< indicates whethere a decomposition has been found */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** detector initialization method */
static
DEC_DECL_INITDETECTOR(initBorderheur)
{

   int i;
   int nvars;
   int nconss;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);
   /* initialize variables and constraints per block structures*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varpart, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      detectordata->varpart[i] = -1;
   }

   detectordata->nhyperedges = 0;
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->hedges, nconss) );

   /* initialize hyperedge array */
   for( i = 0; i < nconss; ++i )
   {
      detectordata->hedges[i].cost = 0;
      detectordata->hedges[i].cons = NULL;
   }

   SCIP_CALL( SCIPcreateWallClock(scip, &detectordata->metisclock) );

   return SCIP_OKAY;
}

/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
DEC_DECL_EXITDETECTOR(exitBorderheur)
{
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

   SCIPfreeMemoryArray(scip, &detectordata->partition);
   SCIPfreeMemoryArray(scip, &detectordata->varpart);

   SCIP_CALL( SCIPfreeClock(scip, &detectordata->metisclock) );
   SCIPfreeMemoryArray(scip, &detectordata->hedges);
   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** compute weight of a hyperedge */
static
SCIP_RETCODE computeHyperedgeWeight(
   SCIP*             scip,          /**< SCIP data structure */
   DEC_DETECTORDATA* detectordata,  /**< detector data data structure */
   SCIP_CONS*        cons,          /**< constraint belonging to the hyperegde */
   int*              cost           /**< pointer storing the hyperedge cost */
   )
{ /*lint --e{715}*/
   *cost = detectordata->consWeight;

   return SCIP_OKAY;
}

/**
 * builds a graph structure out of the matrix.
 *
 * The function will create an HyperEdge for every constraint or an hyperedge for every variable depending on the type of bordered searched.
 * The weight of the hyperedges can be specified.
 *
 * @todo The nonzeroness is not checked, all variables in the variable array are considered
 */
static SCIP_RETCODE buildGraphStructure(
      SCIP*             scip,          /**< SCIP data structure */
      DEC_DETECTORDATA* detectordata   /**< presolver data data structure */
      )
{
   SCIP_Bool valid;
   SCIP_CONS **conss;
   int nconss = SCIPgetNConss( scip );
   int nvars = SCIPgetNVars( scip );
   int i;

   assert(scip != NULL);
   assert(detectordata != NULL);

   conss = SCIPgetConss(scip);
   detectordata->nhyperedges = nconss;
   /* go through all constraints */
   valid = FALSE;
   for( i = 0; i < nconss; ++i )
   {
      int ncurvars;
      SCIP_VAR** curvars;
      int j;
      assert(detectordata->hedges[i].cons == NULL);
      assert(detectordata->hedges[i].cost == 0);


      ncurvars = SCIPgetNVarsXXX(scip, conss[i] );
      curvars = NULL;
      if( ncurvars > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, conss[i], curvars, ncurvars) );
      }

      for( j = 0; j < ncurvars; ++j )
      {
         assert(curvars != NULL);
         if( SCIPisVarRelevant(curvars[j]) )
         {
            valid = TRUE;
            break;
         }
      }
      if( !valid )
      {
         --(detectordata->nhyperedges);
         continue;
      }

      /* compute its weight*/
      SCIP_CALL( computeHyperedgeWeight(scip, detectordata, conss[i], &(detectordata->hedges[i].cost)) );

      detectordata->hedges[i].cons = conss[i];
      SCIPfreeBufferArrayNull(scip, &curvars);
      valid = FALSE;

   }
   detectordata->nvertices = nvars;
   return SCIP_OKAY;
}

/** will call hmetis via a system call */
static
SCIP_RETCODE callMetis(
   SCIP*             scip,          /**< SCIP data struture */
   DEC_DETECTORDATA* detectordata,  /**< presolver data data structure */
   SCIP_RESULT*      result         /**< result indicating whether the detection was successful */
   )
{
   char metiscall[SCIP_MAXSTRLEN];
   char metisout[SCIP_MAXSTRLEN];
   char line[SCIP_MAXSTRLEN];
   char tempfile[SCIP_MAXSTRLEN];

   int status;
   int i;
   int j;
   int nvertices;
   int nhyperedges;
   int ndummyvertices;
   int* partition;

   SCIP_FILE *zfile;
   FILE* file;
   int temp_filedes;
   SCIP_Real remainingtime;

   assert(scip != NULL);
   assert(detectordata != NULL);

   *result = SCIP_DIDNOTRUN;

   remainingtime = DECgetRemainingTime(scip);

   if( remainingtime <= 0 )
   {
      return SCIP_OKAY;
   }

   nvertices = detectordata->nvertices;
   nhyperedges = detectordata->nhyperedges;
   /*lint --e{524}*/
   ndummyvertices = SCIPceil(scip, detectordata->dummynodes*nvertices);

   (void) SCIPsnprintf(tempfile, SCIP_MAXSTRLEN, "gcg-metis-XXXXXX");
   if( (temp_filedes = mkstemp(tempfile)) < 0 )
   {
      SCIPerrorMessage("Error creating temporary file: %s\n", strerror( errno ));
      return SCIP_FILECREATEERROR;
   }

   SCIPdebugMessage("Temporary filename: %s\n", tempfile);

   file = fdopen(temp_filedes, "w");
   if( file == NULL )
   {
      SCIPerrorMessage("Could not open temporary metis file!\n");
      return SCIP_FILECREATEERROR;
   }

   SCIPinfoMessage(scip, file, "%d %d 1\n", nhyperedges, nvertices+ndummyvertices);
   for( i = 0; i < SCIPgetNConss(scip); i++ )
   {
      int ncurvars;
      SCIP_VAR** curvars;
      HyperEdge hedge;
      hedge = detectordata->hedges[i];
      if( hedge.cons == NULL )
         continue;
      SCIPinfoMessage(scip, file, "%d ", hedge.cost);
      ncurvars = SCIPgetNVarsXXX(scip, hedge.cons);
      curvars = NULL;
      if( ncurvars > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, hedge.cons, curvars, ncurvars) );
      }

      for( j = 0; j < ncurvars; ++j )
      {
         int ind;
         assert(curvars != NULL);
         ind = SCIPvarGetProbindex(SCIPvarGetProbvar(curvars[j]));

         assert(ind < SCIPgetNVars(scip));
         if( ind >= 0 )
            SCIPinfoMessage(scip, file, "%d ", ind + 1 );
      }
      SCIPfreeBufferArrayNull(scip, &curvars);
      SCIPinfoMessage(scip, file, "\n");
   }
   status = fclose(file);

   if( status == -1 )
   {
      SCIPerrorMessage("Could not close '%s'\n", tempfile);
      return SCIP_WRITEERROR;
   }

   /* call metis via syscall as there is no library usable ... */
   if( !SCIPisInfinity(scip, remainingtime) )
   {
      (void) SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"ulimit -t %.0f; hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"",
               remainingtime,
               tempfile,
               detectordata->blocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }
   else
   {
      (void) SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"",
               tempfile,
               detectordata->blocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }

   SCIP_CALL( SCIPresetClock(scip, detectordata->metisclock) );
   SCIP_CALL( SCIPstartClock(scip, detectordata->metisclock) );
   SCIPdebugMessage("Calling metis with: %s\n", metiscall);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " %d", detectordata->blocks );
/*   extern char **environ;
   char **env;
   for( env = environ; *env; ++env )
        printf("%s\n", *env);
*/
   status = system( metiscall );

   SCIP_CALL( SCIPstopClock(scip, detectordata->metisclock) );
   SCIPdebugMessage("time left before metis started: %f, time metis spend %f, remainingtime: %f\n", remainingtime, SCIPgetClockTime(scip, detectordata->metisclock),  remainingtime-SCIPgetClockTime(scip, detectordata->metisclock) );

   /* check error codes */
   if( status == -1 )
   {
      SCIPerrorMessage("System call did not succed: %s\n", strerror( errno ));
      SCIPerrorMessage("Call was %s\n", metiscall);
   }
   else if( status != 0 )
   {

      SCIPerrorMessage("Calling hmetis unsuccessful! See the above error message for more details.\n");
      SCIPerrorMessage("Call was %s\n", metiscall);
   }

   /* exit gracefully in case of errors */
   if( status != 0 )
   {
      if( detectordata->tidy )
      {
         status = unlink( tempfile );
         if( status == -1 )
         {
            SCIPerrorMessage("Could not remove metis input file: ", strerror( errno ));
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

   (void) SCIPsnprintf(metisout, SCIP_MAXSTRLEN, "%s.part.%d",tempfile, detectordata->blocks);

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
      assert(temp >= 0 && temp <= detectordata->blocks-1);
      partition[i] = temp;
      i++;
   }
   SCIPfclose(zfile);

   /* if desired delete the temoprary metis file */
   if( detectordata->tidy )
   {
      status = unlink( tempfile );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis input file: %s\n", strerror( errno ));
         return SCIP_WRITEERROR;
      }
      status = unlink( metisout );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis output file: %s\n", strerror( errno ));
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

/** maps the partitions for the disaggregated vertices to the original vertices */
static
SCIP_RETCODE assignBlocksToOriginalVariables(
   SCIP*             scip,          /**< SCIP data structure */
   DEC_DETECTORDATA* detectordata   /**< presolver data data structure */
   )
{

   int i;
   int *partition;
   int *origpart;
   int nvertices;
   assert(scip != NULL);
   assert(detectordata != NULL);

   nvertices = detectordata->nvertices;
   partition = detectordata->partition;
   origpart = detectordata->varpart;
#ifndef NDEBUG
   {
      int nvars;
      nvars = SCIPgetNVars( scip );
      assert(nvertices == nvars);
   }
#endif
    /* go through the new vertices */
   for( i = 0; i < nvertices ; ++i )
   {
      int originalId;
      /* find out the original id (== index of the var in the vars array) */
      originalId = i;
      origpart[originalId] = partition[i];

      assert(origpart[originalId] == -2 || origpart[originalId] >= 0);
      assert(origpart[originalId] <= detectordata->blocks);
   }

   return SCIP_OKAY;
}

/** builds the transformed problem in the new scip instance */
static SCIP_RETCODE buildTransformedProblem(
   SCIP*                    scip,           /**< SCIP data structure */
   DEC_DETECTORDATA*        detectordata,   /**< presolver data data structure */
   DECDECOMP*               decdecomp,      /**< decdecomp data structure */
   int                      nblocks,        /**< number of blocks for this decomposition */
   SCIP_RESULT*             result          /**< indicates whether a structure was found*/
   )
{
   SCIP_Bool *isVarHandled;

   SCIP_VAR*** subscipvars;
   SCIP_CONS*** subscipconss;
   SCIP_CONS** linkingconss;
   SCIP_HASHMAP* vartoblock;
   SCIP_HASHMAP* constoblock;

   int *nsubscipvars;
   int *nsubscipconss;
   int nlinkingconss;
   int i;
   int j;
   SCIP_CONS **conss;
   int nconss;
   SCIP_VAR **vars;
   int nvars;
   SCIP_Bool emptyblocks = FALSE;

   assert(scip != NULL);
   assert(detectordata != NULL);

   nconss = SCIPgetNConss( scip );
   nvars = SCIPgetNVars( scip );

   conss = SCIPgetConss( scip );
   vars = SCIPgetVars( scip );

   SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars, nblocks) );

   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipvars, nblocks) );

   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss[i], nconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars[i], nvars) );

      nsubscipconss[i] = 0;
      nsubscipvars[i] = 0;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &linkingconss, nconss) );
   nlinkingconss = 0;

   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPhashmapCreate(&vartoblock, SCIPblkmem(scip), nconss) );

   SCIP_CALL( SCIPallocBufferArray(scip, &isVarHandled, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      isVarHandled[i] = FALSE;
   }

   /* go through all of the constraints */
   for( i = 0; i < nconss; i++ )
   {
      int consblock = -1;
      int ncurvars;
      SCIP_VAR **curvars;

      if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])), "origbranch") == 0 )
         continue;

      /* sort the variables into corresponding buckets */
      ncurvars = SCIPgetNVarsXXX( scip, conss[i] );
      curvars = NULL;
      if( ncurvars > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, conss[i], curvars, ncurvars) );
      }

      for( j = 0; j < ncurvars; j++ )
      {
         SCIP_VAR* var;
         int varblock;
         assert(curvars != NULL);

         if( !SCIPisVarRelevant(curvars[j]) )
            continue;

         var = SCIPvarGetProbvar(curvars[j]);

         assert(var != NULL);
         assert(SCIPvarIsActive(var));
         assert(!SCIPvarIsDeleted(var));
         /*
          * if the variable has already been handled, we do not need to look
          * at it again and only need to set the constraint
          */
         if( !isVarHandled[SCIPvarGetProbindex(var)] )
         {
            isVarHandled[SCIPvarGetProbindex( var )] = TRUE;
            /*
             * if the following assertion fails, the picture and the mapping is
             * certainly wrong!
             */
            assert(vars[SCIPvarGetProbindex(var)] == var);
            assert(SCIPvarGetProbindex(var) < nvars);
            /* get the number of blocks the current variable is in */

            /* the variable is in exactly one block */
            /* then the partition is given */
            varblock = detectordata->varpart[SCIPvarGetProbindex(var)];
            assert(varblock < detectordata->blocks);
            subscipvars[varblock][nsubscipvars[varblock]] = var;
            ++(nsubscipvars[varblock]);

            /* finally set the hashmap image */
            assert(!SCIPhashmapExists(vartoblock, var));
            SCIP_CALL( SCIPhashmapInsert(vartoblock, var, (void*) (size_t) varblock) );
         }
         else
         {
            varblock = (int)(size_t)SCIPhashmapGetImage(vartoblock, var);
            assert(varblock == detectordata->varpart[SCIPvarGetProbindex(var)] ||  detectordata->varpart[SCIPvarGetProbindex(var)] == -2);
         }

         assert(varblock < detectordata->blocks);
         assert(varblock >= 0 || varblock == -2);
         /*
          * if the variable is not a linking variable, add it to the correct
          * block and update the block of the constraint, if necessary
          */
         if( varblock <= detectordata->blocks && varblock >= 0 )
         {
            /*
             * if the block of the constraint has not been set yet, set it to
             * the block of the current variable
             */
            if( consblock == -1 )
            {
               consblock = varblock;
            }
            /*
             * if the block has been set but the current variable belongs to a
             * different block, we have a linking constraint
             */
            else if( consblock >= 0 && consblock != varblock )
            {
               consblock = -2;
            }

            /*
             * At this point the constraint is either in the border or it
             * belongs to the same block as the current variable
             */
            assert(consblock == -2 ||consblock == varblock);

         }
      }
      SCIPfreeBufferArrayNull(scip, &curvars);

      /*
       *  sort the constraints into the corresponding bucket
       *
       *  if the constraint is linking, put it there
       */
      if( consblock < 0 )
      {
         size_t block;
         assert(detectordata->blocks >= 0);
         /*lint --e{732} */
         block = detectordata->blocks +1;
         linkingconss[nlinkingconss] = conss[i];
         ++nlinkingconss;
         assert(!SCIPhashmapExists(constoblock, conss[i]));
         SCIP_CALL( SCIPhashmapInsert(constoblock, conss[i], (void*)(block)) );

      }
      /* otherwise put it in its block */
      else
      {
         subscipconss[consblock][nsubscipconss[consblock]] = conss[i];
         assert(!SCIPhashmapExists(constoblock, conss[i]));
         SCIP_CALL( SCIPhashmapInsert(constoblock, conss[i], (void*) (size_t) consblock) );
         ++(nsubscipconss[consblock]);
      }

   }

   /*
    * go through all variables and look at the not handled variables and add
    * them to the correct partition
    */
   for( i = 0; i < nvars; i++ )
   {
      int partitionOfVar;
      if( isVarHandled[i] )
         continue;

      partitionOfVar = -1;
      if( detectordata->varpart[i] >= 0 )
      {
         partitionOfVar = detectordata->varpart[i];
      }

      subscipvars[partitionOfVar][nsubscipvars[partitionOfVar]] = vars[i];
      ++nsubscipvars[partitionOfVar];

   }

   SCIPfreeBufferArray(scip, &isVarHandled);

   /* do some elimentary checks and report errors */
   /* first, make sure that there are constraints in every block, otherwise the hole thing is useless */
   for( i = 0; i < detectordata->blocks; ++i )
   {
      if( nsubscipconss[i] == 0 )
      {
         SCIPdebugMessage("Block %d does not have any constraints!\n", i);
         emptyblocks = TRUE;
      }
   }

   if( !emptyblocks )
   {
      /* copy the local data to the decomp structure */
      DECdecdecompSetNBlocks(decdecomp, nblocks);
      DECdecdecompSetType(decdecomp, DEC_DECTYPE_DIAGONAL);
      SCIP_CALL( DECdecdecompSetSubscipvars(scip, decdecomp, subscipvars, nsubscipvars) );
      SCIP_CALL( DECdecdecompSetSubscipconss(scip, decdecomp, subscipconss, nsubscipconss) );
      if( nlinkingconss > 0 )
      {
         SCIP_CALL( DECdecdecompSetLinkingconss(scip, decdecomp, linkingconss, nlinkingconss) );
         DECdecdecompSetType(decdecomp, DEC_DECTYPE_BORDERED);
      }
      DECdecdecompSetVartoblock(decdecomp, vartoblock);
      DECdecdecompSetConstoblock(decdecomp, constoblock);
   }
   else {
      SCIPhashmapFree(&constoblock);
      SCIPhashmapFree(&vartoblock);
   }

   /* free all local data */
   for( i = 0; i < nblocks; ++i )
   {
      SCIPfreeBufferArray(scip, &subscipconss[i]);
      SCIPfreeBufferArray(scip, &subscipvars[i]);
   }

   SCIPfreeBufferArray(scip, &subscipconss);
   SCIPfreeBufferArray(scip, &subscipvars);
   SCIPfreeBufferArray(scip, &nsubscipconss);
   SCIPfreeBufferArray(scip, &nsubscipvars);
   SCIPfreeBufferArray(scip, &linkingconss);

   *result = emptyblocks? SCIP_DIDNOTFIND:SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** detection call back method */
static
DEC_DECL_DETECTSTRUCTURE(detectAndBuildBordered)
{
   int i;
   int j;
   int ndecs;
   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(decdecomps != NULL);
   assert(ndecdecomps != NULL);

   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);
   ndecs = detectordata->maxblocks-detectordata->minblocks+1;
   *ndecdecomps = 0;
   /* allocate space for output data */
   assert(detectordata->maxblocks >= detectordata->minblocks);
   SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, ndecs) );

   /* build the hypergraph structure from the original problem */
   SCIP_CALL( buildGraphStructure(scip, detectordata) );

   for( i = 0; i < ndecs; ++i )
   {
      SCIP_CALL( DECdecdecompCreate(scip, &(*decdecomps)[i]) );
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting bordered structure:");
   for( j = 0, i = detectordata->minblocks; i <= detectordata->maxblocks; ++i )
   {
      detectordata->blocks = i;
      /* get the partitions for the new variables from metis */
      SCIP_CALL( callMetis(scip, detectordata, result) );

      if( *result != SCIP_SUCCESS )
      {
         *result = SCIP_DIDNOTFIND;
         return SCIP_OKAY;
      }
      else
      {
         detectordata->found = TRUE;
      }
      /* deduce the partitions for the original variables */
      SCIP_CALL( assignBlocksToOriginalVariables( scip, detectordata) );

      SCIP_CALL( buildTransformedProblem(scip, detectordata, (*decdecomps)[j], i, result) );
      if( *result == SCIP_SUCCESS )
      {
         *ndecdecomps += 1;
         ++j;
      }
   }
   for( i = *ndecdecomps; i < ndecs; ++i )
   {
      DECdecdecompFree(scip,  &(*decdecomps)[i]);
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d decompositions found.\n", *ndecdecomps );
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/** creates the borderheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionBorderheur(
   SCIP*                 scip              /**< SCIP data structure */

   )
{
   DEC_DETECTORDATA *detectordata;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );

   assert(detectordata != NULL);
   detectordata->found = FALSE;
   detectordata->partition = NULL;
   detectordata->blocks = -1;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_PRIORITY, DEC_ENABLED, detectordata, detectAndBuildBordered, initBorderheur, exitBorderheur) );

   /* add borderheur presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/borderheur/maxblocks", "The maximal number of blocks", &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/borderheur/minblocks", "The minimal number of blocks", &detectordata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/borderheur/consWeight", "Weight of a constraint hyperedge", &detectordata->consWeight, FALSE, DEFAULT_CONSWEIGHT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/borderheur/tidy", "Whether to clean up temporary files", &detectordata->tidy, FALSE, DEFAULT_TIDY, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/borderheur/randomseed", "random seed for hmetis", &detectordata->randomseed, FALSE, DEFAULT_RANDSEED, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detectors/borderheur/dummynodes", "percentage of dummy nodes for metis", &detectordata->dummynodes, FALSE, DEFAULT_DUMMYNODES, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detectors/borderheur/ubfactor", "Unbalance factor for metis", &detectordata->metisubfactor, FALSE, DEFAULT_METIS_UBFACTOR, 0.0, 1E20, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/borderheur/metisverbose", "Should the metis output be displayed", &detectordata->metisverbose, FALSE, DEFAULT_METIS_VERBOSE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/borderheur/metisuseptyperb", "Should the rb or kway method be used for partitioning by metis", &detectordata->metisuseptyperb, FALSE, DEFAULT_METISUSEPTYPE_RB, NULL, NULL) );

   return SCIP_OKAY;
}
