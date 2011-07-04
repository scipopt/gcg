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

#include "dec_borderheur.h"
#include "scip_misc.h"
#include "scip/scipdefplugins.h"

//#include <cstdio>
#include <math.h>
#include <string.h>
#include <errno.h>

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
struct SCIP_BorderheurData
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
      SCIP_BORDERHEURDATA*   borderheurdata,
      SCIP_BORDERHEURSCORES* scores
      )
{
   char name[SCIP_MAXSTRLEN];
   SCIPsnprintf(name, SCIP_MAXSTRLEN, DWSOLVER_REFNAME(SCIPgetProbName(scip),
         borderheurdata->blocks,
         borderheurdata->consWeight,
         borderheurdata->dummynodes));

   return SCIP_OKAY;
}

static
SCIP_RETCODE initBorderheurData(
      SCIP* scip,
      SCIP_BORDERHEURDATA* borderheurdata
      )

{

   int i;
   int nvars;
   int nconss;

   assert(scip != NULL);
   assert(borderheurdata != NULL);

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   borderheurdata->maxblocks = MIN(nconss, borderheurdata->maxblocks);
   /* initialize variables and constraints per block structures*/
   SCIP_CALL(SCIPallocMemoryArray(scip, &borderheurdata->consperblock, borderheurdata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &borderheurdata->varsperblock, borderheurdata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &borderheurdata->nconsperblock, borderheurdata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &borderheurdata->nvarsperblock, borderheurdata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &borderheurdata->linkingconss, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &borderheurdata->varpart, nvars));
   for(i = 0; i < nvars; ++i)
   {
      borderheurdata->varpart[i] = -1;
   }
   SCIP_CALL(SCIPallocMemoryArray(scip, &borderheurdata->hedges, nconss));
   for( i = 0; i < borderheurdata->maxblocks; ++i)
   {
      borderheurdata->nvarsperblock[i] = 0;
      borderheurdata->nconsperblock[i] = 0;
      SCIP_CALL(SCIPallocMemoryArray(scip, &borderheurdata->consperblock[i], nconss));
      SCIP_CALL(SCIPallocMemoryArray(scip, &borderheurdata->varsperblock[i], nvars));
   }

   borderheurdata->nlinkingconss = 0;
   borderheurdata->nhyperedges = 0;
   /* create variable and constraint hash tables */
   SCIP_CALL(SCIPhashmapCreate(&borderheurdata->varstoblock, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&borderheurdata->constoblock, SCIPblkmem(scip), nconss));
   SCIP_CALL(SCIPhashmapCreate(&borderheurdata->constolpid, SCIPblkmem(scip), nconss));

   /* initialise consttolpid hashmap and hyper edges array */
   for( i = 0; i < nconss; ++i)
   {
      SCIP_CALL(SCIPhashmapInsert(borderheurdata->constolpid, SCIPgetConss(scip)[i], (void*)(size_t)i));
      borderheurdata->hedges[i].cost = 0;
      borderheurdata->hedges[i].cons = NULL;
   }

   return SCIP_OKAY;
}

/** copies the variable and block information to the decomp structure */
static
SCIP_RETCODE copyBorderheurDataToDecomp(
      SCIP*             scip,       /**< SCIP data structure */
      SCIP_BORDERHEURDATA*  borderheurdata, /**< presolver data data structure */
      DECDECOMP*        decomp      /**< DECOMP data structure */
      )
{
   int i;
   assert(scip != 0);
   assert(borderheurdata != 0);
   assert(decomp != 0);

   SCIP_CALL(SCIPallocMemoryArray(scip, &decomp->subscipvars, borderheurdata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &decomp->subscipconss, borderheurdata->blocks));

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->linkingconss, borderheurdata->linkingconss, borderheurdata->nlinkingconss));
   decomp->nlinkingconss = borderheurdata->nlinkingconss;

   for( i = 0; i < borderheurdata->blocks; ++i)
   {
      SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->subscipconss[i], borderheurdata->consperblock[i], borderheurdata->nconsperblock[i]));
      SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->subscipvars[i], borderheurdata->varsperblock[i], borderheurdata->nvarsperblock[i]));
   }

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipconss, borderheurdata->nconsperblock, borderheurdata->blocks));
   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipvars, borderheurdata->nvarsperblock, borderheurdata->blocks));

   decomp->constoblock = borderheurdata->constoblock;
   decomp->vartoblock = borderheurdata->varstoblock;
   decomp->nblocks = borderheurdata->blocks;
   decomp->type = DEC_BORDERED;
   return SCIP_OKAY;
}

/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
SCIP_RETCODE freeBorderheurDataData
(
      SCIP* scip,
      SCIP_BORDERHEURDATA* borderheurdata
)
{

   int i;

   assert(scip != NULL);
   assert(borderheurdata != NULL);

   /* quick fix to correct the pictures */
//   SCIP_CALL(detectAndBuildArrowHead(scip, borderheurdata));

   /* copy data to decomp structure */
   if( !borderheurdata->found)
   {
      return SCIP_OKAY;
   }

   SCIP_CALL(copyBorderheurDataToDecomp(scip, borderheurdata, borderheurdata->decdecomp));
   /* free presolver data */
   for( i = 0; i < borderheurdata->maxblocks; ++i)
   {
      SCIPfreeMemoryArray(scip, &borderheurdata->consperblock[i]);
      SCIPfreeMemoryArray(scip, &borderheurdata->varsperblock[i]);
   }
   SCIPfreeMemoryArray(scip, &borderheurdata->varsperblock);
   SCIPfreeMemoryArray(scip, &borderheurdata->nvarsperblock);
   SCIPfreeMemoryArray(scip, &borderheurdata->consperblock);
   SCIPfreeMemoryArray(scip, &borderheurdata->nconsperblock);
   SCIPfreeMemoryArray(scip, &borderheurdata->linkingconss);
   SCIPfreeMemoryArray(scip, &borderheurdata->partition);
   SCIPfreeMemoryArray(scip, &borderheurdata->varpart);

   /* free hash map */
   SCIPhashmapFree(&borderheurdata->constolpid);
   /* TODO: Hashmap is not copied! but shallow copied, so do not free here!
   SCIPhashmapFree(&borderheurdata->varstoblock);
   SCIPhashmapFree(&borderheurdata->constoblock);
   */

   SCIPfreeMemoryArray(scip, &borderheurdata->hedges);

   return SCIP_OKAY;
}

static
int computeHyperedgeWeight(SCIP *scip, SCIP_BORDERHEURDATA *borderheurdata, SCIP_CONS* cons, int *cost)
{
   *cost = borderheurdata->consWeight;

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
      SCIP*             scip,       /**< SCIP data structure */
      SCIP_BORDERHEURDATA*  borderheurdata  /**< presolver data data structure */
      )
{
   SCIP_Bool valid;
   SCIP_CONS **conss;
   int nconss = SCIPgetNConss( scip );
   int nvars = SCIPgetNVars( scip );
   int i;

   assert(scip != NULL);
   assert(borderheurdata != NULL);

   conss = SCIPgetConss(scip);
   borderheurdata->nhyperedges = nconss;
   /* go through all constraints */
   valid = FALSE;
   for( i = 0; i < nconss; ++i )
   {
      int ncurvars;
      SCIP_VAR** curvars;
      int j;
      assert(borderheurdata->hedges[i].cons == NULL);
      assert(borderheurdata->hedges[i].cost == 0);


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
         --(borderheurdata->nhyperedges);
         continue;
      }
      /* compute its weight*/
      SCIP_CALL(computeHyperedgeWeight(scip, borderheurdata, conss[i], &(borderheurdata->hedges[i].cost)));

      borderheurdata->hedges[i].cons = conss[i];
      SCIPfreeMemoryArray(scip, &curvars);
      valid = FALSE;

   }
   borderheurdata->nvertices = nvars;
   return SCIP_OKAY;
}

/** Will call hmetis via a system call */
static
SCIP_RETCODE callMetis(
      SCIP*             scip,       /**< SCIP data struture */
      SCIP_BORDERHEURDATA*  borderheurdata  /**< presolver data data structure */
      )
{
   char metiscall[SCIP_MAXSTRLEN];
   char metisout[SCIP_MAXSTRLEN];
   char line[SCIP_MAXSTRLEN];
   const char *tempfile = "metis.temp";

   int status;
   int i;
   int j;
   int nvertices;
   int nhyperedges;
   int ndummyvertices;
   int* partition;

   SCIP_FILE *zfile;
   FILE* file;

   assert(scip != NULL);
   assert(borderheurdata != NULL);
   assert(!SCIPfileExists(tempfile));

   nvertices = borderheurdata->nvertices;
   nhyperedges = borderheurdata->nhyperedges;
   ndummyvertices = borderheurdata->dummynodes*nvertices;
//   SCIPinfoMessage(scip, NULL, "DUMMY: %.2f", borderheurdata->dummynodes);

   file = fopen(tempfile, "w");
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
      hedge = borderheurdata->hedges[i];
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
   if(borderheurdata->metisverbose)
   {
      if(borderheurdata->metisuseptyperb)
      {
         SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "./hmetis metis.temp %d -seed %d -ptype rb -ufactor %f",  borderheurdata->blocks,  borderheurdata->randomseed, borderheurdata->metisubfactor );
      }
      else
      {
         SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "./hmetis metis.temp %d -seed %d -ptype kway -ufactor %f",  borderheurdata->blocks,  borderheurdata->randomseed, borderheurdata->metisubfactor );
      }
   }
   else
   {
      if(borderheurdata->metisuseptyperb)
      {
         SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "./hmetis metis.temp %d -seed %d -ptype rb -ufactor %f >/dev/null",  borderheurdata->blocks,  borderheurdata->randomseed, borderheurdata->metisubfactor );
      }
      else
      {
         SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "./hmetis metis.temp %d -seed %d -ptype kway -ufactor %f >/dev/null",  borderheurdata->blocks,  borderheurdata->randomseed, borderheurdata->metisubfactor );
      }
   }
//   SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "./hmetis.sh metis.temp %d -seed %d",  borderheurdata->blocks,  borderheurdata->randomseed );
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
      if( borderheurdata->tidy )
      {
         status = remove( tempfile );
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
   if(borderheurdata->partition == NULL)
   {
      SCIP_CALL(SCIPallocMemoryArray(scip, &borderheurdata->partition, nvertices));
   }

   assert(borderheurdata->partition != NULL);
   partition = borderheurdata->partition;

   SCIPsnprintf(metisout, SCIP_MAXSTRLEN, "metis.temp.part.%d",borderheurdata->blocks);

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
      assert(temp >= 0 && temp <= borderheurdata->blocks-1);
      partition[i] = temp;
      i++;
   }
   SCIPfclose(zfile);

   /* if desired delete the temoprary metis file */
   if( borderheurdata->tidy )
   {
      status = remove( tempfile );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis input file: ", strerror( errno ));
      }
      status = remove( metisout );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis output file: ", strerror( errno ));
      }
   }

   return SCIP_OKAY;
}
/** Maps the partitions for the disaggregated vertices to the original vertices */
static
SCIP_RETCODE assignBlocksToOriginalVariables(
      SCIP*             scip,       /**< SCIP data structure */
      SCIP_BORDERHEURDATA*  borderheurdata  /**< presolver data data structure */
      )
{

   int i;
   int nvars;
   int *partition;
   int *origpart;
   int nvertices;
   assert(scip != NULL);
   assert(borderheurdata != NULL);

   nvertices = borderheurdata->nvertices;
   partition = borderheurdata->partition;
   origpart = borderheurdata->varpart;
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
      assert(origpart[originalId] <= borderheurdata->blocks);
   }

   return SCIP_OKAY;
}

/** Builds the transformed problem in the new scip instance */
static SCIP_RETCODE buildTransformedProblem(
   SCIP*                    scip,           /**< SCIP data structure */
   SCIP_BORDERHEURDATA*      borderheurdata,  /**< presolver data data structure */
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

   subscipconss = borderheurdata->consperblock;
   subscipvars = borderheurdata->varsperblock;

   nsubscipconss = borderheurdata->nconsperblock;
   nsubscipvars = borderheurdata->nvarsperblock;

   linkingconss = borderheurdata->linkingconss;
   nlinkingconss = borderheurdata->nlinkingconss;

   constoblock = borderheurdata->constoblock;
   vartoblock = borderheurdata->varstoblock;

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
         varblock = -1;
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
            varblock = borderheurdata->varpart[SCIPvarGetProbindex(var)];
            assert(varblock < borderheurdata->blocks);
            subscipvars[varblock][nsubscipvars[varblock]] = var;
            ++(nsubscipvars[varblock]);

            /* finally set the hashmap image */
            assert(!SCIPhashmapExists(vartoblock, var));
            SCIP_CALL(SCIPhashmapInsert(vartoblock, var, (void*)varblock));
         }
         else
         {
            varblock = (long int)SCIPhashmapGetImage(vartoblock, var);
            assert(varblock == borderheurdata->varpart[SCIPvarGetProbindex(var)] ||  borderheurdata->varpart[SCIPvarGetProbindex(var)] == -2);
         }


         /*
          * if the variable is not a linking variable, add it to the correct
          * block and update the block of the constraint, if necessary
          */
         if( varblock <= borderheurdata->blocks )
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
         assert(borderheurdata->blocks >= 0);
         block = borderheurdata->blocks +1;
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
      if( borderheurdata->varpart[i] >= 0 )
      {
         partitionOfVar = borderheurdata->varpart[i];
      }

      subscipvars[partitionOfVar][nsubscipvars[partitionOfVar]] = vars[i];
      ++nsubscipvars[partitionOfVar];


   }
   SCIPfreeMemoryArray(scip, &isVarHandled);
   borderheurdata->nlinkingconss = nlinkingconss;
   /* do some elimentary checks and report errors */

   /* first, make sure that there are constraints in every block, otherwise the hole thing is useless */
   for( i = 0; i < borderheurdata->blocks; ++i)
   {
      if(nsubscipconss[i] == 0)
      {
         SCIPerrorMessage("Block %d does not have any constraints!\n", i);
      }
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE writeDWsolverOutput(
      SCIP*             scip,       /**< SCIP data structure */
      SCIP_BORDERHEURDATA*  borderheurdata  /**< presolver data data structure */
   )
{
   FILE* file;
   int i;
   int j;
   char name[SCIP_MAXSTRLEN];

   SCIPsnprintf(name, SCIP_MAXSTRLEN,
         DWSOLVER_REFNAME(SCIPgetProbName(scip),
                          borderheurdata->blocks,
                          borderheurdata->consWeight,
                          borderheurdata->dummynodes)
                          );

   file = fopen(name, "w");
   if(file == NULL)
      return SCIP_FILECREATEERROR;

   SCIPinfoMessage(scip, file, "%d ", borderheurdata->blocks );
   for(i = 0; i < borderheurdata->blocks; ++i)
   {
      SCIPinfoMessage(scip, file, "%d ", borderheurdata->nconsperblock[i]);
   }
   SCIPinfoMessage(scip, file, "\n");

   for(i = 0; i < borderheurdata->blocks; ++i)
   {
      for(j = 0; j < borderheurdata->nconsperblock[i]; ++j)
      {
         long int consindex = (long int) SCIPhashmapGetImage(borderheurdata->constolpid, borderheurdata->consperblock[i][j]);
         SCIPinfoMessage(scip, file, "%d ", consindex);
      }

      SCIPinfoMessage(scip, file, "\n");
   }

   fclose(file);

   return SCIP_OKAY;
}

static
SCIP_RETCODE evaluateDecomposition(
      SCIP*                 scip,           /**< SCIP data structure */
      SCIP_BORDERHEURDATA*   borderheurdata,  /**< presolver data data structure */
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
   assert(borderheurdata != NULL);
   assert(score != NULL);

   matrixarea = 0;
   borderarea = 0;
   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);

   /* get the right name */
   SCIPsnprintf(name, SCIP_MAXSTRLEN,
         DWSOLVER_REFNAME(SCIPgetProbName(scip),
                          borderheurdata->blocks,
                          borderheurdata->consWeight,
                          borderheurdata->dummynodes)
                                  );

   SCIP_CALL(SCIPallocMemoryArray(scip, &nzblocks, borderheurdata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &nlinkvarsblocks, borderheurdata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &blockdensities, borderheurdata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &blocksizes, borderheurdata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &nvarsblocks, borderheurdata->blocks));
   /*
    * 3 Scores
    *
    * - Area percentage (min)
    * - block density (max)
    * - \pi_b {v_b|v_b is linking}/#vb (min)
    */

   /* calculate matrix area */
   matrixarea = nvars*nconss;

   //SCIPinfoMessage(scip, NULL, "Sizes: %d x %d (%d, %d)\n", nvars, nconss, borderheurdata->nlinkingvars, borderheurdata->nlinkingconss);

   /* calculate slave sizes, nonzeros and linkingvars */
   for (i = 0; i < borderheurdata->blocks; ++i)
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
      curconss = borderheurdata->consperblock[i];
      ncurconss = borderheurdata->nconsperblock[i];

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
            assert(SCIPhashmapExists(borderheurdata->varstoblock, curvars[k]));
            block = (long int) SCIPhashmapGetImage(borderheurdata->varstoblock, curvars[k]);
            //SCIPinfoMessage(scip, NULL, "b: %d", block);
            if(block == borderheurdata->blocks+1 && ishandled[SCIPvarGetProbindex(curvars[k])] == FALSE)
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

   borderarea = borderheurdata->nlinkingconss*nvars;
   //SCIPinfoMessage(scip, NULL, "c: %d, v: %d", borderheurdata->nlinkingconss*nvars, borderheurdata->nlinkingvars*(nconss-borderheurdata->nlinkingconss));

   blockarea = 0;
   density = 1E20;
   varratio = 1.0;
   for( i = 0; i < borderheurdata->blocks; ++i )
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

extern
SCIP_RETCODE createBorderheurData(
      SCIP*         scip,
      SCIP_BORDERHEURDATA** borderheurdata
      )
{
   assert(scip != NULL);
   assert(borderheurdata != NULL);
   SCIP_CALL(SCIPallocBlockMemory(scip, borderheurdata));
   return SCIP_OKAY;
}

extern
void freeBorderheurData(
      SCIP* scip,
      SCIP_BORDERHEURDATA** borderheurdata
   )
{
   assert(scip != NULL);
   assert(borderheurdata != NULL);

   SCIPfreeBlockMemoryNull(scip, borderheurdata);
}

extern
SCIP_RETCODE detectAndBuildBordered(
      SCIP*                scip,          /**< SCIP data structure */
      SCIP_BORDERHEURDATA*  borderheurdata, /**< presolver data data structure */
      SCIP_RESULT*         result
      )
{

   SCIP_BORDERHEURSCORES score;
   int i;
   char filename[SCIP_MAXSTRLEN];

   //SCIPinfoMessage(scip, NULL, "detectandbuild arrowhead:\n");
   assert(scip != NULL);

   SCIP_CALL(initBorderheurData(scip, borderheurdata));
   /* build the hypergraph structure from the original problem */
   SCIP_CALL(buildGraphStructure(scip, borderheurdata));

   if( borderheurdata->minblocks == borderheurdata->maxblocks)
   {
      borderheurdata->blocks = borderheurdata->minblocks;
     /* get the partitions for the new variables from metis */
     SCIP_CALL(callMetis(scip, borderheurdata));

     /* deduce the partitions for the original variables */
     SCIP_CALL(assignBlocksToOriginalVariables( scip, borderheurdata));

     SCIP_CALL(buildTransformedProblem(scip, borderheurdata, &score));
     SCIP_CALL(evaluateDecomposition(scip, borderheurdata, &score));
     SCIP_CALL(writeDWsolverOutput(scip, borderheurdata));

     borderheurdata->found = TRUE;
     SCIP_CALL(printBorderheurScores(scip, borderheurdata, &score));
     SCIP_CALL(freeBorderheurDataData(scip, borderheurdata));

     SCIPsnprintf(filename, SCIP_MAXSTRLEN,
           GP_NAME(SCIPgetProbName(scip),
              borderheurdata->blocks,
              borderheurdata->consWeight,
              borderheurdata->dummynodes)
              );

     SCIP_CALL(SCIPwriteOrigProblem(scip, filename, "gp", FALSE));

   }
   else
   {
      SCIP_Real bestscore = 1E20;
      int bestsetting;
      for( i = borderheurdata->minblocks; i <= borderheurdata->maxblocks; ++i)
      {
         SCIP_Real cumscore;
         borderheurdata->blocks = i;

         /* get the partitions for the new variables from metis */
         SCIP_CALL(callMetis(scip, borderheurdata));

         /* deduce the partitions for the original variables */
         SCIP_CALL(assignBlocksToOriginalVariables( scip, borderheurdata));

         SCIP_CALL(buildTransformedProblem(scip, borderheurdata,  &score));
         SCIP_CALL(evaluateDecomposition(scip, borderheurdata, &score));

         cumscore = score.borderscore*score.linkingscore*score.densityscore;
         if (cumscore < bestscore)
         {
            bestscore = cumscore;
            bestsetting = i;
         }

         SCIPhashmapFree(&borderheurdata->varstoblock);
         SCIPhashmapFree(&borderheurdata->constoblock);
         SCIP_CALL(SCIPhashmapCreate(&borderheurdata->varstoblock, SCIPblkmem(scip), SCIPgetNVars(scip)));
         SCIP_CALL(SCIPhashmapCreate(&borderheurdata->constoblock, SCIPblkmem(scip), SCIPgetNConss(scip)));
         for(i = 0; i <  SCIPgetNVars(scip); ++i)
         {
            borderheurdata->varpart[i] = -1;
         }

         for( i = 0; i < borderheurdata->blocks; ++i)
         {
            borderheurdata->nvarsperblock[i] = 0;
            borderheurdata->nconsperblock[i] = 0;
         }

         borderheurdata->nlinkingconss = 0;
      }

      borderheurdata->found = TRUE;
      borderheurdata->blocks = bestsetting;

      /* get the partitions for the new variables from metis */
      SCIP_CALL(callMetis(scip, borderheurdata));

      /* deduce the partitions for the original variables */
      SCIP_CALL(assignBlocksToOriginalVariables( scip, borderheurdata));

      SCIP_CALL(buildTransformedProblem(scip, borderheurdata, &score));
      SCIP_CALL(evaluateDecomposition(scip, borderheurdata, &score));
      SCIP_CALL(writeDWsolverOutput(scip, borderheurdata));

      borderheurdata->found = TRUE;
      SCIP_CALL(printBorderheurScores(scip, borderheurdata, &score));
      SCIP_CALL(freeBorderheurDataData(scip, borderheurdata));

      SCIPsnprintf(filename, SCIP_MAXSTRLEN,
           GP_NAME(SCIPgetProbName(scip),
              borderheurdata->blocks,
              borderheurdata->consWeight,
              borderheurdata->dummynodes)
              );

      SCIP_CALL(SCIPwriteOrigProblem(scip, filename, "gp", FALSE));
   }
   *result = SCIP_OKAY;
   return SCIP_OKAY;
}

/** set the decomp structure */
extern
SCIP_RETCODE SCIPBorderheurSetDecomp(
   SCIP*       scip,       /**< SCIP data structure */
   SCIP_BORDERHEURDATA* borderheurdata,
   DECDECOMP*  decdecomp   /**< DECOMP data structure */

   )
{

   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(borderheurdata != NULL);

   borderheurdata->decdecomp = decdecomp;

   return SCIP_OKAY;
}

/** creates the borderheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionBorderheur(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_BORDERHEURDATA*   borderheurdata
   )
{
   assert(scip != NULL);
   assert(borderheurdata != NULL);

   borderheurdata->found = FALSE;
   borderheurdata->partition = NULL;
   borderheurdata->blocks = -1;
   /* add borderheur presolver parameters */
   SCIP_CALL(SCIPaddIntParam(scip, "borderheur/maxblocks", "The maximal number of blocks", &borderheurdata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "borderheur/minblocks", "The minimal number of blocks", &borderheurdata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "borderheur/consWeight", "Weight of a constraint hyperedge", &borderheurdata->consWeight, FALSE, DEFAULT_CONSWEIGHT, 0, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddBoolParam(scip, "borderheur/tidy", "Whether to clean up temporary files", &borderheurdata->tidy, FALSE, DEFAULT_TIDY, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "borderheur/randomseed", "random seed for hmetis", &borderheurdata->randomseed, FALSE, DEFAULT_RANDSEED, -1, INT_MAX, NULL, NULL));
   SCIP_CALL(SCIPaddRealParam(scip, "borderheur/dummynodes", "percentage of dummy nodes for metis", &borderheurdata->dummynodes, FALSE, DEFAULT_DUMMYNODES, 0.0, 1.0, NULL, NULL));
   SCIP_CALL(SCIPaddRealParam(scip, "borderheur/ubfactor", "Unbalance factor for metis", &borderheurdata->metisubfactor, FALSE, DEFAULT_METIS_UBFACTOR, 0.0, 1E20, NULL, NULL ));
   SCIP_CALL(SCIPaddBoolParam(scip, "borderheur/metisverbose", "Should the metis output be displayed", &borderheurdata->metisverbose, FALSE, DEFAULT_METIS_VERBOSE, NULL, NULL ));
   SCIP_CALL(SCIPaddBoolParam(scip, "borderheur/metisuseptyperb", "Should the rb or kway method be used for partitioning by metis", &borderheurdata->metisuseptyperb, FALSE, DEFAULT_METISUSEPTYPE_RB, NULL, NULL));
   return SCIP_OKAY;
}
