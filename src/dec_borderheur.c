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
 * @ingroup DETECTORS
 * @brief  borderheur presolver
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG

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

#define DEC_DETECTORNAME      "borderheur"   /**< name of the detector */
#define DEC_PRIORITY          0              /**< priority of the detector */

/* Default parameter settings */
#define DEFAULT_CONSWEIGHT                5     /**< weight for constraint hyperedges */
#define DEFAULT_RANDSEED                  1     /**< random seed for the hmetis call */
#define DEFAULT_TIDY                      TRUE  /**< whether to clean up afterwards */
#define DEFAULT_DUMMYNODES                0.2   /**< percentage of dummy vertices*/

#define DEFAULT_MAXBLOCKS                 20    /**< value for the maximum number of blocks to be considered */
#define DEFAULT_MINBLOCKS                 2     /**< value for the minimum number of blocks to be considered */

#define DEFAULT_METIS_UBFACTOR            5.0   /**< default unbalance factor given to metis on the commandline */
#define DEFAULT_METIS_VERBOSE             FALSE /**< should metis be verbose */
#define DEFAULT_METISUSEPTYPE_RB          TRUE  /**< Should metis use the rb or kway partitioning algorithm */
#define DEFAULT_PRIORITY                  DEC_PRIORITY

#define DWSOLVER_REFNAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_ref.txt", (name), (blocks), (cons), (dummy)

/*
 * Data structures
 */

struct hyperedge
{
   SCIP_CONS* cons;   ///< the original pointer to the constraint
   int cost;

};
typedef struct hyperedge HyperEdge;

/** score data structure **/
struct Dec_BorderheurScores
{
   SCIP_Real borderscore;
   SCIP_Real minkequicutscore;
   SCIP_Real equicutscorenormalized;
   SCIP_Real densityscore;
   SCIP_Real linkingscore;
};
typedef struct Dec_BorderheurScores DEC_BORDERHEURSCORES;

/** detector data */
struct DEC_DetectorData
{
   /* Graph stuff for hmetis */
   HyperEdge *hedges;
   int *partition;
   int nvertices;
   int nhyperedges;
   int *varpart;

   /* Stuff to get the dw-solver to work*/
   SCIP_HASHMAP *constolpid;

   SCIP_Bool tidy;
   int blocks;
   int maxblocks;
   int minblocks;
   int consWeight;
   int randomseed;
   SCIP_Bool found;
   SCIP_Real dummynodes;

   SCIP_Real metisubfactor;
   SCIP_Bool metisverbose;
   SCIP_Bool metisuseptyperb;
   SCIP_CLOCK *metisclock;
   int priority;

};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** Prints the score of the decomposition */
static
SCIP_RETCODE printBorderheurScores(
      SCIP*                 scip,           /**< SCIP data structure */
      DEC_DETECTORDATA*     detectordata,   /**< detectordata data structure */
      DEC_BORDERHEURSCORES* scores          /**< score data structure */
      )
{
   char name[SCIP_MAXSTRLEN];
   SCIPsnprintf(name, SCIP_MAXSTRLEN, DWSOLVER_REFNAME(SCIPgetProbName(scip),
         detectordata->blocks,
         detectordata->consWeight,
         detectordata->dummynodes));

   return SCIP_OKAY;
}

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
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varpart, nvars));
   for( i = 0; i < nvars; ++i )
   {
      detectordata->varpart[i] = -1;
   }

   detectordata->nhyperedges = 0;
   SCIP_CALL(SCIPhashmapCreate(&detectordata->constolpid, SCIPblkmem(scip), nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->hedges, nconss));

   /* initialise consttolpid hashmap and hyper edges array */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL(SCIPhashmapInsert(detectordata->constolpid, SCIPgetConss(scip)[i], (void*)(size_t)i));
      detectordata->hedges[i].cost = 0;
      detectordata->hedges[i].cons = NULL;
   }

   SCIP_CALL(SCIPcreateWallClock(scip, &detectordata->metisclock));

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
   if( !detectordata->found)
   {
      SCIPfreeMemory(scip, &detectordata);
      return SCIP_OKAY;
   }

   SCIPfreeMemoryArray(scip, &detectordata->partition);
   SCIPfreeMemoryArray(scip, &detectordata->varpart);

   /* free hash map */
   SCIPhashmapFree(&detectordata->constolpid);

   SCIP_CALL(SCIPfreeClock(scip, &detectordata->metisclock));
   SCIPfreeMemoryArray(scip, &detectordata->hedges);
   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

static
int computeHyperedgeWeight(
   SCIP*             scip,
   DEC_DETECTORDATA* detectordata,
   SCIP_CONS* cons,
   int *cost
   )
{
   *cost = detectordata->consWeight;

   return SCIP_OKAY;
}

/**
 * Builds a graph structure out of the matrix.
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
      curvars = SCIPgetVarsXXX(scip, conss[i] );
      for(j = 0; j < ncurvars; ++j)
      {
         if( isVarRelevant(curvars[j]) )
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
      SCIP_CALL(computeHyperedgeWeight(scip, detectordata, conss[i], &(detectordata->hedges[i].cost)));

      detectordata->hedges[i].cons = conss[i];
      SCIPfreeMemoryArray(scip, &curvars);
      valid = FALSE;

   }
   detectordata->nvertices = nvars;
   return SCIP_OKAY;
}

/** Will call hmetis via a system call */
static
SCIP_RETCODE callMetis(
      SCIP*              scip,          /**< SCIP data struture */
      DEC_DETECTORDATA*  detectordata,  /**< presolver data data structure */
      SCIP_RESULT*       result         /**< result indicating whether the detection was successful */
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
   ndummyvertices = detectordata->dummynodes*nvertices;

   SCIPsnprintf(tempfile, SCIP_MAXSTRLEN, "gcg-metis-XXXXXX");
   if( (temp_filedes = mkstemp(tempfile)) < 0 )
   {
      SCIPerrorMessage("Error creating temporary file: %s\n", strerror( errno ));
      return SCIP_FILECREATEERROR;
   }

   SCIPdebugMessage("Temporary filename: %s\n", tempfile);

   file = fdopen(temp_filedes, "w");
   if(file == NULL)
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
      if(hedge.cons == NULL)
         continue;
      SCIPinfoMessage(scip, file, "%d ", hedge.cost);
      ncurvars = SCIPgetNVarsXXX(scip, hedge.cons);
      curvars = SCIPgetVarsXXX(scip, hedge.cons);
      for( j = 0; j < ncurvars; ++j )
      {
         int ind = SCIPvarGetProbindex(SCIPvarGetProbvar(curvars[j]));
         assert(ind < SCIPgetNVars(scip));
         if( ind >= 0)
            SCIPinfoMessage(scip, file, "%d ", ind + 1 );
      }
      SCIPfreeMemoryArray(scip, &curvars);
      SCIPinfoMessage(scip, file, "\n");
   }
   status = fclose(file);

   if(status == -1)
   {
      SCIPerrorMessage("Could not close '%s'\n", tempfile);
      return SCIP_WRITEERROR;
   }

   /* call metis via syscall as there is no library usable ... */
   if(!SCIPisInfinity(scip, remainingtime))
   {
      SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"ulimit -t %.0f; hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"",
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
      SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"",
               tempfile,
               detectordata->blocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }

   SCIP_CALL(SCIPresetClock(scip, detectordata->metisclock));
   SCIP_CALL(SCIPstartClock(scip, detectordata->metisclock));
   SCIPdebugMessage("Calling metis with: %s\n", metiscall);
/*   extern char **environ;
   char **env;
   for (env = environ; *env; ++env)
        printf("%s\n", *env);
*/
   status = system( metiscall );

   SCIP_CALL(SCIPstopClock(scip, detectordata->metisclock));
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
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->partition, nvertices));
   }

   assert(detectordata->partition != NULL);
   partition = detectordata->partition;

   SCIPsnprintf(metisout, SCIP_MAXSTRLEN, "%s.part.%d",tempfile, detectordata->blocks);

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

/** Maps the partitions for the disaggregated vertices to the original vertices */
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

/** Builds the transformed problem in the new scip instance */
static SCIP_RETCODE buildTransformedProblem(
   SCIP*                    scip,           /**< SCIP data structure */
   DEC_DETECTORDATA*        detectordata,   /**< presolver data data structure */
   DECDECOMP*               decdecomp,      /**< decdecomp data structure */
   int                      nblocks,        /**< number of blocks for this decomposition */
   DEC_BORDERHEURSCORES*    score,          /**< scores */
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
   assert(score != NULL);

   nconss = SCIPgetNConss( scip );
   nvars = SCIPgetNVars( scip );

   conss = SCIPgetConss( scip );
   vars = SCIPgetVars( scip );

   SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars, nblocks) );

   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipvars, nblocks) );

   for (i = 0; i < nblocks; ++i )
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

   score->minkequicutscore = 0;
   score->equicutscorenormalized = 0;

   SCIP_CALL(SCIPallocMemoryArray(scip, &isVarHandled, nvars));
   for( i = 0; i < nvars; ++i )
   {
      isVarHandled[i] = FALSE;
   }

   /* go through all of the constraints */
   for( i = 0; i < nconss; i++ )
   {
      long int consblock = -1;

      /* sort the variables into corresponding buckets */
      int ncurvars = SCIPgetNVarsXXX( scip, conss[i] );
      SCIP_VAR **curvars = SCIPgetVarsXXX( scip, conss[i] );
      for( j = 0; j < ncurvars; j++ )
      {
         SCIP_VAR* var;
         long int varblock;
         if( !isVarRelevant(curvars[j]) )
         {
//            SCIPprintVar(scip, curvars[j], NULL);
            continue;
         }
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
            SCIP_CALL(SCIPhashmapInsert(vartoblock, var, (void*)varblock));
         }
         else
         {
            varblock = (long int)SCIPhashmapGetImage(vartoblock, var);
            assert(varblock == detectordata->varpart[SCIPvarGetProbindex(var)] ||  detectordata->varpart[SCIPvarGetProbindex(var)] == -2);
         }

         assert(varblock < detectordata->blocks);
         assert(varblock >= 0 || varblock == -2);
         /*
          * if the variable is not a linking variable, add it to the correct
          * block and update the block of the constraint, if necessary
          */
         if( varblock <= detectordata->blocks && varblock >= 0)
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
      SCIPfreeMemoryArrayNull(scip, &curvars);

      /*
       *  sort the constraints into the corresponding bucket
       *
       *  if the constraint is linking, put it there
       */
      if( consblock < 0 )
      {
         size_t block;
         assert(detectordata->blocks >= 0);
         block = detectordata->blocks +1;
         linkingconss[nlinkingconss] = conss[i];
         ++nlinkingconss;
         assert(!SCIPhashmapExists(constoblock, conss[i]));
         SCIP_CALL(SCIPhashmapInsert(constoblock, conss[i], (void*)(block)));

      }
      /* otherwise put it in its block */
      else
      {
         subscipconss[consblock][nsubscipconss[consblock]] = conss[i];
         assert(!SCIPhashmapExists(constoblock, conss[i]));
         SCIP_CALL(SCIPhashmapInsert(constoblock, conss[i], (void*)consblock));
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
      {
         continue;

      }

      partitionOfVar = -1;
      if( detectordata->varpart[i] >= 0 )
      {
         partitionOfVar = detectordata->varpart[i];
      }

      subscipvars[partitionOfVar][nsubscipvars[partitionOfVar]] = vars[i];
      ++nsubscipvars[partitionOfVar];

   }

   SCIPfreeMemoryArray(scip, &isVarHandled);

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
      DECdecdecompSetType(decdecomp, DEC_DECTYPE_BORDERED);
      SCIP_CALL( DECdecdecompSetSubscipvars(scip, decdecomp, subscipvars, nsubscipvars) );
      SCIP_CALL( DECdecdecompSetSubscipconss(scip, decdecomp, subscipconss, nsubscipconss) );
      SCIP_CALL( DECdecdecompSetLinkingconss(scip, decdecomp, linkingconss, nlinkingconss) );
      DECdecdecompSetVartoblock(decdecomp, vartoblock);
      DECdecdecompSetConstoblock(decdecomp, constoblock);
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

static
SCIP_RETCODE evaluateDecomposition(
      SCIP*                 scip,           /**< SCIP data structure */
      DEC_DETECTORDATA*     detectordata,   /**< presolver data data structure */
      DECDECOMP*            decdecomp,      /**< decomposition data structure */
      DEC_BORDERHEURSCORES* score           /**< returns the score of the decomposition */
      )
{
   char name[SCIP_MAXSTRLEN];
   long int matrixarea;
   long int borderarea;
   int nvars;
   int nconss;
   int i;
   int j;
   int k;
   int blockarea;
   SCIP_Real varratio;
   int* nzblocks;
   int nblocks;
   int* nlinkvarsblocks;
   int* nvarsblocks;
   SCIP_Real* blockdensities;
   int* blocksizes;
   SCIP_Real density;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(score != NULL);

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);

   nblocks = DECdecdecompGetNBlocks(decdecomp);

   /* get the right name */
   SCIPsnprintf(name, SCIP_MAXSTRLEN,
         DWSOLVER_REFNAME(SCIPgetProbName(scip),
            nblocks,
            detectordata->consWeight,
            detectordata->dummynodes)
      );

   SCIP_CALL(SCIPallocMemoryArray(scip, &nzblocks, nblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &nlinkvarsblocks, nblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &blockdensities, nblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &blocksizes, nblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &nvarsblocks, nblocks));
   /*
    * 3 Scores
    *
    * - Area percentage (min)
    * - block density (max)
    * - \pi_b {v_b|v_b is linking}/#vb (min)
    */

   /* calculate matrix area */
   matrixarea = nvars*nconss;

   /* calculate slave sizes, nonzeros and linkingvars */
   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CONS** curconss;
      int ncurconss;
      int nvarsblock;
      SCIP_Bool *ishandled;

      SCIP_CALL(SCIPallocMemoryArray(scip, &ishandled, nvars));
      nvarsblock = 0;
      nzblocks[i] = 0;
      nlinkvarsblocks[i] = 0;
      for( j = 0; j < nvars; ++j )
      {
         ishandled[j] = FALSE;
      }
      curconss = DECdecdecompGetSubscipconss(decdecomp)[i];
      ncurconss = DECdecdecompGetNSubscipconss(decdecomp)[i];

      for( j = 0; j < ncurconss; ++j )
      {
         SCIP_VAR** curvars;
         SCIP_VAR* var;
         int ncurvars;

         curvars = SCIPgetVarsXXX(scip, curconss[j]);
         ncurvars = SCIPgetNVarsXXX(scip, curconss[j]);
         for( k = 0; k < ncurvars; ++k )
         {
            long int block;
            if( !isVarRelevant(curvars[k]) )
               continue;

            var = SCIPvarGetProbvar(curvars[k]);
            assert(var != NULL);
            assert(SCIPvarIsActive(var));
            assert(!SCIPvarIsDeleted(var));
            ++(nzblocks[i]);
            assert(SCIPhashmapExists(DECdecdecompGetVartoblock(decdecomp), var));
            block = (long int) SCIPhashmapGetImage(DECdecdecompGetVartoblock(decdecomp), var);

            if( block == nblocks + 1 && ishandled[SCIPvarGetProbindex(var)] == FALSE )
            {
               ++(nlinkvarsblocks[i]);
            }
            ishandled[SCIPvarGetProbindex(var)] = TRUE;
         }

         SCIPfreeMemoryArray(scip, &curvars);
      }

      for( j = 0; j < nvars; ++j )
      {
         if( ishandled[j] )
         {
            ++nvarsblock;
         }
      }

      blocksizes[i] = nvarsblock*ncurconss;
      nvarsblocks[i] = nvarsblock;
      if(blocksizes[i] > 0)
      {
         blockdensities[i] = 1.0*nzblocks[i]/blocksizes[i];
      }
      else
      {
         blockdensities[i] = 0.0;
      }

      assert(blockdensities[i] >= 0 && blockdensities[i] <= 1.0);
      SCIPfreeMemoryArray(scip, &ishandled);
   }

   /* calculate border area */
   borderarea = DECdecdecompGetNLinkingconss(decdecomp)*nvars;

   blockarea = 0;
   density = 1E20;
   varratio = 1.0;
   for( i = 0; i < nblocks; ++i )
   {
      /* calculate block area */
      blockarea += blocksizes[i];


      /* calculate density */
      density = MIN(density, blockdensities[i]);
   }

   score->linkingscore = (0.5+0.5*varratio);
   score->borderscore = (1.0*(borderarea)/matrixarea);
   score->densityscore = (1-density);

   SCIPfreeMemoryArray(scip, &nzblocks);
   SCIPfreeMemoryArray(scip, &nlinkvarsblocks);
   SCIPfreeMemoryArray(scip, &blockdensities);
   SCIPfreeMemoryArray(scip, &blocksizes);
   SCIPfreeMemoryArray(scip, &nvarsblocks);
   return SCIP_OKAY;

}


static
DEC_DECL_DETECTSTRUCTURE(detectAndBuildBordered)
{

   DEC_BORDERHEURSCORES* scores;
   SCIP_Real* cumscores;
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
   SCIP_CALL( SCIPallocMemoryArray(scip, &scores, ndecs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &cumscores, ndecs) );

   /* build the hypergraph structure from the original problem */
   SCIP_CALL(buildGraphStructure(scip, detectordata));

   for(i = 0; i < ndecs; ++i)
   {
      SCIP_CALL_ABORT( DECdecdecompCreate(scip, &(*decdecomps)[i]) );
   }


   for(j = 0, i = detectordata->minblocks; i <= detectordata->maxblocks; ++i )
   {
      detectordata->blocks = i;
      /* get the partitions for the new variables from metis */
      SCIP_CALL(callMetis(scip, detectordata, result));

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

      SCIP_CALL( buildTransformedProblem(scip, detectordata, (*decdecomps)[j], i, &scores[j], result) );
      if( *result == SCIP_SUCCESS )
      {
         SCIP_CALL( evaluateDecomposition(scip, detectordata, (*decdecomps)[j], &scores[j]) );
         SCIP_CALL( printBorderheurScores(scip, detectordata, &scores[j]) );

         cumscores[j] = scores[j].borderscore*scores[j].linkingscore*scores[j].densityscore;
         *ndecdecomps += 1;
         ++j;
      }
   }

   SCIPsortRealPtr(cumscores,  (void**) *decdecomps, *ndecdecomps);

   SCIPfreeMemoryArray(scip, &cumscores);
   SCIPfreeMemoryArray(scip, &scores);

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/** gets the priority of the detector */
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

/** creates the borderheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionBorderheur(
   SCIP*                 scip              /**< SCIP data structure */

   )
{
   DEC_DETECTORDATA *detectordata;
   assert(scip != NULL);

   SCIP_CALL(SCIPallocMemory(scip, &detectordata));

   assert(detectordata != NULL);
   detectordata->found = FALSE;
   detectordata->partition = NULL;
   detectordata->blocks = -1;

   SCIP_CALL(DECincludeDetector(scip, DEC_DETECTORNAME, detectordata, detectAndBuildBordered, initBorderheur, exitBorderheur, getPriority));

   /* add borderheur presolver parameters */
   SCIP_CALL(SCIPaddIntParam(scip, "borderheur/maxblocks", "The maximal number of blocks", &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "borderheur/minblocks", "The minimal number of blocks", &detectordata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "borderheur/consWeight", "Weight of a constraint hyperedge", &detectordata->consWeight, FALSE, DEFAULT_CONSWEIGHT, 0, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddBoolParam(scip, "borderheur/tidy", "Whether to clean up temporary files", &detectordata->tidy, FALSE, DEFAULT_TIDY, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "borderheur/randomseed", "random seed for hmetis", &detectordata->randomseed, FALSE, DEFAULT_RANDSEED, -1, INT_MAX, NULL, NULL));
   SCIP_CALL(SCIPaddRealParam(scip, "borderheur/dummynodes", "percentage of dummy nodes for metis", &detectordata->dummynodes, FALSE, DEFAULT_DUMMYNODES, 0.0, 1.0, NULL, NULL));
   SCIP_CALL(SCIPaddRealParam(scip, "borderheur/ubfactor", "Unbalance factor for metis", &detectordata->metisubfactor, FALSE, DEFAULT_METIS_UBFACTOR, 0.0, 1E20, NULL, NULL ));
   SCIP_CALL(SCIPaddBoolParam(scip, "borderheur/metisverbose", "Should the metis output be displayed", &detectordata->metisverbose, FALSE, DEFAULT_METIS_VERBOSE, NULL, NULL ));
   SCIP_CALL(SCIPaddBoolParam(scip, "borderheur/metisuseptyperb", "Should the rb or kway method be used for partitioning by metis", &detectordata->metisuseptyperb, FALSE, DEFAULT_METISUSEPTYPE_RB, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "borderheur/priority", "random seed for hmetis", &detectordata->priority, FALSE, DEFAULT_PRIORITY, INT_MIN, INT_MAX, NULL, NULL));
   return SCIP_OKAY;
}
