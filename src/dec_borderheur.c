/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: presol_borderheur.c,v 1.24 2010/01/04 20:35:45 bzfheinz Exp $"

/**@file   presol_borderheur.c
 * @ingroup PRESOLVERS
 * @brief  borderheur presolver
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <assert.h>

#include "scip/scipdefplugins.h"

#include "dec_borderheur.h"
#include "cons_decomp.h"
#include "struct_decomp.h"
#include "scip_misc.h"

//#include <cstdio>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#define DEC_DETECTORNAME      "borderheur"   /**< name of the detector */
#define DEC_PRIORITY          0              /**< priority of the detector */

/* Default parameter settings*/
#define DEFAULT_BLOCKS                    2     /**< number of blocks */
#define DEFAULT_CONSWEIGHT                5     /**< weight for constraint hyperedges */
#define DEFAULT_RANDSEED                  1     /**< random seed for the hmetis call */
#define DEFAULT_TIDY                      TRUE  /**< whether to clean up afterwards */
#define DEFAULT_DUMMYNODES	              0.2   /**< percentage of dummy vertices*/

#define DEFAULT_MAXBLOCKS                 10    /**< value for the maximum number of blocks to be considered */
#define DEFAULT_MINBLOCKS                 2     /**< value for the minimum number of blocks to be considered */

#define DEFAULT_METIS_UBFACTOR            5.0   /**< default unbalance factor given to metis on the commandline */
#define DEFAULT_METIS_VERBOSE             FALSE /**< should metis be verbose */
#define DEFAULT_METISUSEPTYPE_RB          TRUE  /**< Should metis use the rb or kway partitioning algorithm */

#define DWSOLVER_REFNAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_ref.txt", (name), (blocks), (cons), (dummy)

#define GP_NAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_%d.gp", (name), (blocks), (cons), (dummy)

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
struct SCIP_BorderheurScores
{
   SCIP_Real borderscore;
   SCIP_Real minkequicutscore;
   SCIP_Real equicutscorenormalized;
   SCIP_Real densityscore;
   SCIP_Real linkingscore;
};
typedef struct SCIP_BorderheurScores SCIP_BORDERHEURSCORES;
/** presolver data */
struct DEC_DetectorData
{
   DECDECOMP* decdecomp;
   SCIP_VAR*** varsperblock;
   int* nvarsperblock;
   SCIP_CONS*** consperblock;
   int *nconsperblock;
   SCIP_CONS** linkingconss;
   int nlinkingconss;

   SCIP_HASHMAP* constoblock;
   SCIP_HASHMAP* varstoblock;

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

};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static
SCIP_RETCODE printBorderheurScores(
      SCIP*                 scip,
      DEC_DETECTORDATA*   detectordata,
      SCIP_BORDERHEURSCORES* scores
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
   DEC_DETECTOR* borderheur;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);

   borderheur = DECfindDetector(scip, DEC_DETECTORNAME);
   assert(borderheur != NULL);

   detectordata = DECdetectorGetData(borderheur);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(borderheur), DEC_DETECTORNAME) == 0);

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);
   /* initialize variables and constraints per block structures*/
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->consperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->nconsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->nvarsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->linkingconss, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varpart, nvars));
   for(i = 0; i < nvars; ++i)
   {
      detectordata->varpart[i] = -1;
   }
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->hedges, nconss));
   for( i = 0; i < detectordata->maxblocks; ++i)
   {
      detectordata->nvarsperblock[i] = 0;
      detectordata->nconsperblock[i] = 0;
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->consperblock[i], nconss));
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varsperblock[i], nvars));
   }

   detectordata->nlinkingconss = 0;
   detectordata->nhyperedges = 0;
   /* create variable and constraint hash tables */
   SCIP_CALL(SCIPhashmapCreate(&detectordata->varstoblock, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip), nconss));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->constolpid, SCIPblkmem(scip), nconss));

   /* initialise consttolpid hashmap and hyper edges array */
   for( i = 0; i < nconss; ++i)
   {
      SCIP_CALL(SCIPhashmapInsert(detectordata->constolpid, SCIPgetConss(scip)[i], (void*)(size_t)i));
      detectordata->hedges[i].cost = 0;
      detectordata->hedges[i].cons = NULL;
   }

   return SCIP_OKAY;
}

/** copies the variable and block information to the decomp structure */
static
SCIP_RETCODE copyDetectorDataToDecomp(
      SCIP*             scip,       /**< SCIP data structure */
      DEC_DETECTORDATA*  detectordata, /**< presolver data data structure */
      DECDECOMP*        decomp      /**< DECOMP data structure */
      )
{
   int i;
   assert(scip != 0);
   assert(detectordata != 0);
   assert(decomp != 0);

   SCIP_CALL(SCIPallocMemoryArray(scip, &decomp->subscipvars, detectordata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &decomp->subscipconss, detectordata->blocks));

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->linkingconss, detectordata->linkingconss, detectordata->nlinkingconss));
   decomp->nlinkingconss = detectordata->nlinkingconss;

   for( i = 0; i < detectordata->blocks; ++i)
   {
      SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->subscipconss[i], detectordata->consperblock[i], detectordata->nconsperblock[i]));
      SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->subscipvars[i], detectordata->varsperblock[i], detectordata->nvarsperblock[i]));
   }

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipconss, detectordata->nconsperblock, detectordata->blocks));
   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipvars, detectordata->nvarsperblock, detectordata->blocks));

   decomp->constoblock = detectordata->constoblock;
   decomp->vartoblock = detectordata->varstoblock;
   decomp->nblocks = detectordata->blocks;
   decomp->type = DEC_BORDERED;
   return SCIP_OKAY;
}

/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
DEC_DECL_EXITDETECTOR(exitBorderheur)
{

   int i;
   DEC_DETECTOR* borderheur;
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);

   borderheur = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(borderheur);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(borderheur), DEC_DETECTORNAME) == 0);

   /* copy data to decomp structure */
   if( !detectordata->found)
   {
      SCIPfreeMemory(scip, &detectordata);
      return SCIP_OKAY;
   }

   /* free presolver data */
   for( i = 0; i < detectordata->maxblocks; ++i)
   {
      SCIPfreeMemoryArray(scip, &detectordata->consperblock[i]);
      SCIPfreeMemoryArray(scip, &detectordata->varsperblock[i]);
   }
   SCIPfreeMemoryArray(scip, &detectordata->varsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->nvarsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->consperblock);
   SCIPfreeMemoryArray(scip, &detectordata->nconsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->linkingconss);
   SCIPfreeMemoryArray(scip, &detectordata->partition);
   SCIPfreeMemoryArray(scip, &detectordata->varpart);

   /* free hash map */
   SCIPhashmapFree(&detectordata->constolpid);
   /* TODO: Hashmap is not copied! but shallow copied, so do not free here!
   SCIPhashmapFree(&detectordata->varstoblock);
   SCIPhashmapFree(&detectordata->constoblock);
   */

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
         if(SCIPvarIsActive(curvars[j]))
         {
            valid = TRUE;
            break;
         }
      }
      if( !valid)
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
      SCIP*             scip,       /**< SCIP data struture */
      DEC_DETECTORDATA*  detectordata  /**< presolver data data structure */
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
   int temp_filedes = -1;

   assert(scip != NULL);
   assert(detectordata != NULL);

   nvertices = detectordata->nvertices;
   nhyperedges = detectordata->nhyperedges;
   ndummyvertices = detectordata->dummynodes*nvertices;
//   SCIPinfoMessage(scip, NULL, "DUMMY: %.2f", detectordata->dummynodes);
   SCIPsnprintf(tempfile, SCIP_MAXSTRLEN, "gcg-metis-XXXXXX");
   if ( (temp_filedes=mkstemp(tempfile)) <0 )
   {
      SCIPerrorMessage("Error creating temporary file: %s", strerror( errno ));
      return SCIP_FILECREATEERROR;
   }

   SCIPdebugMessage("Temporary filename: %s", tempfile);

   file = fdopen(temp_filedes, "w");
   if(file == NULL)
   {
      SCIPerrorMessage("Could not open temporary metis file!");
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
         int ind = SCIPvarGetProbindex(curvars[j]);
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
      SCIPerrorMessage("Could not close '%s'", tempfile);
      return SCIP_WRITEERROR;
   }

   /* call metis via syscall as there is no library usable ... */
   if(detectordata->metisverbose)
   {
      if(detectordata->metisuseptyperb)
      {
         SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "./hmetis %s %d -seed %d -ptype rb -ufactor %f", tempfile, detectordata->blocks,  detectordata->randomseed, detectordata->metisubfactor );
      }
      else
      {
         SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "./hmetis %s %d -seed %d -ptype kway -ufactor %f", tempfile, detectordata->blocks,  detectordata->randomseed, detectordata->metisubfactor );
      }
   }
   else
   {
      if(detectordata->metisuseptyperb)
      {
         SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "./hmetis %s %d -seed %d -ptype rb -ufactor %f >/dev/null", tempfile, detectordata->blocks,  detectordata->randomseed, detectordata->metisubfactor );
      }
      else
      {
         SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "./hmetis %s %d -seed %d -ptype kway -ufactor %f >/dev/null", tempfile, detectordata->blocks,  detectordata->randomseed, detectordata->metisubfactor );
      }
   }
//   SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "./hmetis.sh metis.temp %d -seed %d",  detectordata->blocks,  detectordata->randomseed );
   //SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\nCalling metis with '%s'.\n", metiscall);

   status = system( metiscall );

   /* check error codes */
   if( status == -1 )
   {
      SCIPerrorMessage("System call did not succed: %s", strerror( errno ));
   }
   else if( status != 0 )
   {
      SCIPerrorMessage("Calling hmetis unsuccessful! See the above error message for more details.");
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
         }
      }
      return SCIP_ERROR;
   }

   /*
    * parse the output into the vector
    * alloc the memory
    */
   if(detectordata->partition == NULL)
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
         SCIPerrorMessage("Line could not be read");
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
         SCIPerrorMessage("Could not remove metis input file: %s", strerror( errno ));
      }
      status = unlink( metisout );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis output file: %s", strerror( errno ));
      }
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "Temporary file is in: %s\n", tempfile);
   }

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
   int nvars;
   int *partition;
   int *origpart;
   int nvertices;
   assert(scip != NULL);
   assert(detectordata != NULL);

   nvertices = detectordata->nvertices;
   partition = detectordata->partition;
   origpart = detectordata->varpart;
   nvars = SCIPgetNVars( scip );

   assert(nvertices == nvars);
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
   DEC_DETECTORDATA*      detectordata,  /**< presolver data data structure */
   SCIP_BORDERHEURSCORES*    score           /**< scores */
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
   int ncons;
   SCIP_VAR **vars;
   int nvars;
   assert(scip != NULL);

   subscipconss = detectordata->consperblock;
   subscipvars = detectordata->varsperblock;

   nsubscipconss = detectordata->nconsperblock;
   nsubscipvars = detectordata->nvarsperblock;

   linkingconss = detectordata->linkingconss;
   nlinkingconss = detectordata->nlinkingconss;

   constoblock = detectordata->constoblock;
   vartoblock = detectordata->varstoblock;

   ncons = SCIPgetNConss( scip );
   nvars = SCIPgetNVars( scip );

   conss = SCIPgetConss( scip );
   vars = SCIPgetVars( scip );

   score->minkequicutscore = 0;
   score->equicutscorenormalized = 0;

   SCIP_CALL(SCIPallocMemoryArray(scip, &isVarHandled, nvars));
   for(i = 0; i < nvars; ++i)
   {
      isVarHandled[i] = FALSE;
   }

   /* go through all of the constraints */
   for( i = 0; i < ncons; i++ )
   {
      long int consblock = -1;

      /* sort the variables into corresponding buckets */
      int ncurvars = SCIPgetNVarsXXX( scip, conss[i] );
      SCIP_VAR **curvars = SCIPgetVarsXXX( scip, conss[i] );
      for( j = 0; j < ncurvars; j++ )
      {
         SCIP_VAR* var;
         long int varblock;
         if(!SCIPvarIsActive(curvars[j]))
         {
            continue;
         }
         assert(SCIPvarIsActive(curvars[j]));
         assert(!SCIPvarIsDeleted(curvars[j]));

         var = curvars[j];
         /*
          * if the variable has already been handled, we do not need to look
          * at it again and only need to set the constraint
          */
         if(!isVarHandled[SCIPvarGetProbindex( var )])
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
   detectordata->nlinkingconss = nlinkingconss;
   /* do some elimentary checks and report errors */

   /* first, make sure that there are constraints in every block, otherwise the hole thing is useless */
   for( i = 0; i < detectordata->blocks; ++i)
   {
      if(nsubscipconss[i] == 0)
      {
         SCIPerrorMessage("Block %d does not have any constraints!\n", i);
      }
   }
   return SCIP_OKAY;
}

//static
//SCIP_RETCODE writeDWsolverOutput(
//      SCIP*             scip,       /**< SCIP data structure */
//      DEC_DETECTORDATA*  detectordata  /**< presolver data data structure */
//   )
//{
//   FILE* file;
//   int i;
//   int j;
//   char name[SCIP_MAXSTRLEN];
//
//   SCIPsnprintf(name, SCIP_MAXSTRLEN,
//         DWSOLVER_REFNAME(SCIPgetProbName(scip),
//                          detectordata->blocks,
//                          detectordata->consWeight,
//                          detectordata->dummynodes)
//                          );
//
//   file = fopen(name, "w");
//   if(file == NULL)
//      return SCIP_FILECREATEERROR;
//
//   SCIPinfoMessage(scip, file, "%d ", detectordata->blocks );
//   for(i = 0; i < detectordata->blocks; ++i)
//   {
//      SCIPinfoMessage(scip, file, "%d ", detectordata->nconsperblock[i]);
//   }
//   SCIPinfoMessage(scip, file, "\n");
//
//   for(i = 0; i < detectordata->blocks; ++i)
//   {
//      for(j = 0; j < detectordata->nconsperblock[i]; ++j)
//      {
//         long int consindex = (long int) SCIPhashmapGetImage(detectordata->constolpid, detectordata->consperblock[i][j]);
//         SCIPinfoMessage(scip, file, "%d ", consindex);
//      }
//
//      SCIPinfoMessage(scip, file, "\n");
//   }
//
//   fclose(file);
//
//   return SCIP_OKAY;
//}

static
SCIP_RETCODE evaluateDecomposition(
      SCIP*                 scip,           /**< SCIP data structure */
      DEC_DETECTORDATA*   detectordata,  /**< presolver data data structure */
      SCIP_BORDERHEURSCORES* score           /**< returns the score of the decomposition */
      )
{
   char name[SCIP_MAXSTRLEN];
   unsigned int matrixarea;
   unsigned int borderarea;
   int nvars;
   int nconss;
   int i;
   int j;
   int k;
   int blockarea;
   SCIP_Real varratio;
   int* nzblocks;
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

   /* get the right name */
   SCIPsnprintf(name, SCIP_MAXSTRLEN,
         DWSOLVER_REFNAME(SCIPgetProbName(scip),
                          detectordata->blocks,
                          detectordata->consWeight,
                          detectordata->dummynodes)
                                  );

   SCIP_CALL(SCIPallocMemoryArray(scip, &nzblocks, detectordata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &nlinkvarsblocks, detectordata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &blockdensities, detectordata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &blocksizes, detectordata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &nvarsblocks, detectordata->blocks));
   /*
    * 3 Scores
    *
    * - Area percentage (min)
    * - block density (max)
    * - \pi_b {v_b|v_b is linking}/#vb (min)
    */

   /* calculate matrix area */
   matrixarea = nvars*nconss;

   //SCIPinfoMessage(scip, NULL, "Sizes: %d x %d (%d, %d)\n", nvars, nconss, detectordata->nlinkingvars, detectordata->nlinkingconss);

   /* calculate slave sizes, nonzeros and linkingvars */
   for (i = 0; i < detectordata->blocks; ++i)
   {
      SCIP_CONS** curconss;
      int ncurconss;
      int nvarsblock;
      SCIP_Bool *ishandled;

      SCIP_CALL(SCIPallocMemoryArray(scip, &ishandled, nvars));
      nvarsblock = 0;
      nzblocks[i] = 0;
      nlinkvarsblocks[i] = 0;
      for( j = 0; j < nvars; ++j)
      {
         ishandled[j] = FALSE;
      }
      curconss = detectordata->consperblock[i];
      ncurconss = detectordata->nconsperblock[i];

      for( j = 0; j < ncurconss; ++j)
      {
         SCIP_VAR** curvars;
         int ncurvars;

         curvars = SCIPgetVarsXXX(scip, curconss[j]);
         ncurvars = SCIPgetNVarsXXX(scip, curconss[j]);
         for( k = 0; k < ncurvars; ++k)
         {
            long int block;
            if(!SCIPvarIsActive(curvars[k]))
               continue;


            ++(nzblocks[i]);
            assert(SCIPhashmapExists(detectordata->varstoblock, curvars[k]));
            block = (long int) SCIPhashmapGetImage(detectordata->varstoblock, curvars[k]);
            //SCIPinfoMessage(scip, NULL, "b: %d", block);
            if(block == detectordata->blocks+1 && ishandled[SCIPvarGetProbindex(curvars[k])] == FALSE)
            {

               ++(nlinkvarsblocks[i]);
            }
            ishandled[SCIPvarGetProbindex(curvars[k])] = TRUE;
         }

         SCIPfreeMemoryArray(scip, &curvars);
      }

      for( j = 0; j < nvars; ++j)
      {
         if(ishandled[j])
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

   borderarea = detectordata->nlinkingconss*nvars;
   //SCIPinfoMessage(scip, NULL, "c: %d, v: %d", detectordata->nlinkingconss*nvars, detectordata->nlinkingvars*(nconss-detectordata->nlinkingconss));

   blockarea = 0;
   density = 1E20;
   varratio = 1.0;
   for( i = 0; i < detectordata->blocks; ++i )
   {
      /* calculate block area */
      blockarea += blocksizes[i];


      /* calculate density */
      density = MIN(density, blockdensities[i]);
   }

   score->linkingscore = (0.5+0.5*varratio);
   // score->borderscore = (1.0*(blockarea+borderarea)/matrixarea);
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

   SCIP_BORDERHEURSCORES score;
   int i;
   //char filename[SCIP_MAXSTRLEN];
   DEC_DETECTOR* borderheur;
   DEC_DETECTORDATA* detectordata;
   //SCIPinfoMessage(scip, NULL, "detectandbuild arrowhead:\n");
   assert(scip != NULL);

   borderheur = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(borderheur);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(borderheur), DEC_DETECTORNAME) == 0);
   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);

   /* build the hypergraph structure from the original problem */
   SCIP_CALL(buildGraphStructure(scip, detectordata));

   if( detectordata->minblocks == detectordata->maxblocks)
   {
      detectordata->blocks = detectordata->minblocks;
     /* get the partitions for the new variables from metis */
     SCIP_CALL(callMetis(scip, detectordata));

     /* deduce the partitions for the original variables */
     SCIP_CALL(assignBlocksToOriginalVariables( scip, detectordata));

     SCIP_CALL(buildTransformedProblem(scip, detectordata, &score));
     SCIP_CALL(evaluateDecomposition(scip, detectordata, &score));
//     SCIP_CALL(writeDWsolverOutput(scip, detectordata));

     detectordata->found = TRUE;
     SCIP_CALL(printBorderheurScores(scip, detectordata, &score));

//     SCIPsnprintf(filename, SCIP_MAXSTRLEN,
//           GP_NAME(SCIPgetProbName(scip),
//              detectordata->blocks,
//              detectordata->consWeight,
//              detectordata->dummynodes)
//              );
//
//     SCIP_CALL(SCIPwriteOrigProblem(scip, filename, "gp", FALSE));

   }
   else
   {
      SCIP_Real bestscore = 1E20;
      int bestsetting = -1;
      for( i = detectordata->minblocks; i <= detectordata->maxblocks; ++i)
      {
         SCIP_Real cumscore;
         detectordata->blocks = i;

         /* get the partitions for the new variables from metis */
         SCIP_CALL(callMetis(scip, detectordata));

         /* deduce the partitions for the original variables */
         SCIP_CALL(assignBlocksToOriginalVariables( scip, detectordata));

         SCIP_CALL(buildTransformedProblem(scip, detectordata,  &score));
         SCIP_CALL(evaluateDecomposition(scip, detectordata, &score));

         cumscore = score.borderscore*score.linkingscore*score.densityscore;
         if (cumscore < bestscore)
         {
            bestscore = cumscore;
            bestsetting = i;
         }

         SCIPhashmapFree(&detectordata->varstoblock);
         SCIPhashmapFree(&detectordata->constoblock);
         SCIP_CALL(SCIPhashmapCreate(&detectordata->varstoblock, SCIPblkmem(scip), SCIPgetNVars(scip)));
         SCIP_CALL(SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip), SCIPgetNConss(scip)));
         for(i = 0; i <  SCIPgetNVars(scip); ++i)
         {
            detectordata->varpart[i] = -1;
         }

         for( i = 0; i < detectordata->blocks; ++i)
         {
            detectordata->nvarsperblock[i] = 0;
            detectordata->nconsperblock[i] = 0;
         }

         detectordata->nlinkingconss = 0;
      }

      detectordata->found = TRUE;
      assert(bestsetting >= 0);
      detectordata->blocks = bestsetting;

      /* get the partitions for the new variables from metis */
      SCIP_CALL(callMetis(scip, detectordata));

      /* deduce the partitions for the original variables */
      SCIP_CALL(assignBlocksToOriginalVariables( scip, detectordata));

      SCIP_CALL(buildTransformedProblem(scip, detectordata, &score));
      SCIP_CALL(evaluateDecomposition(scip, detectordata, &score));
//      SCIP_CALL(writeDWsolverOutput(scip, detectordata));

      detectordata->found = TRUE;
      SCIP_CALL(printBorderheurScores(scip, detectordata, &score));

//      SCIPsnprintf(filename, SCIP_MAXSTRLEN,
//           GP_NAME(SCIPgetProbName(scip),
//              detectordata->blocks,
//              detectordata->consWeight,
//              detectordata->dummynodes)
//              );
//
//      SCIP_CALL(SCIPwriteOrigProblem(scip, filename, "gp", FALSE));
   }
   SCIP_CALL(copyDetectorDataToDecomp(scip, detectordata, detectordata->decdecomp));
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** set the decomp structure */
static
DEC_DECL_SETSTRUCTDECOMP(BorderheurSetDecomp)
{
   DEC_DETECTOR* borderheur;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   borderheur = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(borderheur);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(borderheur), DEC_DETECTORNAME) == 0);
   SCIPdebugMessage("Setting decdecomp\n");
   detectordata->decdecomp = decdecomp;

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

   SCIP_CALL(DECincludeDetector(scip, DEC_DETECTORNAME, DEC_PRIORITY, detectordata, detectAndBuildBordered, BorderheurSetDecomp, initBorderheur, exitBorderheur));

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

   return SCIP_OKAY;
}
