/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
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

/**@file   reader_cluster.c
 * @brief  CLUSTER file reader
 *
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_DEBUG */

/* @todo really needed? #include <stdlib.h> */
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/ /* needed for strcasecmp() */
#endif
#include <ctype.h>

#include "reader_cluster.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"

#include "cons_decomp.h"
#include "pub_decomp.h"

#define READER_NAME             "clusterreader"
#define READER_DESC             "file reader for blocks in cluster format"
#define READER_EXTENSION        "cluster"


/*
 * Data structures
 */

/** data for dec reader */
struct SCIP_ReaderData
{
   DEC_DECOMP* decomp;      /**< decomposition data structure*/


   /* Graph stuff for hmetis */
   SCIP_PTRARRAY* hedges;                    /**< variable array of hyperedges */
   SCIP_INTARRAY* copytooriginal;            /**< array mapping copied to original variables */
   int*           partition;                 /**< array storing vertex partitions */
   int            nvertices;                 /**< number of vertices */
   int*           varpart;                   /**< array storing variable partition */

   int            blocks;                    /**< number of blocks */
};

enum htype
{
   VARIABLE, CONSTRAINT
};
typedef enum htype hType;


/** hyper edge data structure */
struct hyperedge
{

   hType type;                /**< the type of the hyperegde (is it a split variable or a real constraint) */
   int *variableIds;          /**< the associated variable IDs that appear in the hyperedge */
   int nvariableIds;          /**< number of variable ids */
   int originalId;            /**< the original SCIP ID of this constraint or variable */
};

typedef struct hyperedge HyperEdge;

/** initializes reader data structure */
static SCIP_RETCODE initReaderdata(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_READERDATA*      readerdata          /**< reader data structure */
   )
{
   int i;

   assert(scip != NULL);
   assert(readerdata != NULL);

   SCIP_CALL( SCIPcreatePtrarray(scip, &readerdata->hedges) );
   SCIP_CALL( SCIPcreateIntarray(scip, &readerdata->copytooriginal) );
   SCIP_CALL( DECdecompCreate(scip, &readerdata->decomp) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &readerdata->varpart, SCIPgetNVars(scip)) );

   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      readerdata->varpart[i] = -1;
   }

   /*
    * parse the output into the vector
    * alloc the memory
    */
   readerdata->partition = NULL;
   readerdata->blocks = 0;
   readerdata->nvertices = 0;

   return SCIP_OKAY;
}

/**
 * builds a graph structure out of the matrix.
 *
 * The function will create an HyperEdge for every constraint and every variable.
 * It will additionally create vertices for every variable and in particular
 * a copy of this variable for every constraint in which the variable has a
 * nonzero coefficient. The copies will be connected by the hyperedge for
 * the particular constraint and all copies of a variable will be connected by
 * the hyperedge belonging to that variable. The weight of these variable
 * hyperedges can be specified.
 *
 * @todo The nonzeroness is not checked, all variables in the variable array are considered
 */
static SCIP_RETCODE buildGraphStructure(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_READERDATA*      readerdata          /**< reader data structure */
   )
{
   SCIP_CONS **conss;
   SCIP_PTRARRAY *hedges;
   SCIP_INTARRAY *copytoorig;
   SCIP_INTARRAY** maporigtocopies;
   int *nmaporigtocopies;
   int nconss = SCIPgetNConss( scip );
   int nvars = SCIPgetNVars( scip );
   int id = 0;
   int nhyperedges = 0;
   int nvertices;
   int i;
   int j;
   int *copies;

   HyperEdge *hedge;
   assert(scip != NULL);
   assert(readerdata != NULL);

   conss = SCIPgetConss(scip);
   hedges = readerdata->hedges;
   copytoorig = readerdata->copytooriginal;

   SCIP_CALL( SCIPallocMemoryArray(scip, &copies, nconss) );

   /* we need at least nconss + nvars hyperedges */
   SCIP_CALL( SCIPextendPtrarray(scip, hedges,  0, nconss+nvars) );

   /* we have at least nvars may copy vertices */
   SCIP_CALL( SCIPextendIntarray(scip, copytoorig, 0, nvars) );

   /* map the original variable to all of its copies */
   SCIP_CALL( SCIPallocMemoryArray(scip, &maporigtocopies, nvars) );

   /* these are the number of copies for the given variable */
   SCIP_CALL( SCIPallocMemoryArray(scip, &nmaporigtocopies, nvars) );

   /* initialize for every variable the list of copies */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPcreateIntarray(scip, &maporigtocopies[i]) );
      nmaporigtocopies[i] = 0;
   }

   /* go through all constraints */
   for( i = 0; i < nconss; ++i )
   {
      int *varids;
      SCIP_VAR **vars;
      /* get the number of nonzeros in this constraint  */

      int ncurvars = SCIPgetNVarsXXX( scip, conss[i] );

      /* if there are no variables, skip the constraint */
      if( ncurvars == 0 )
         continue;

      /*
       * may work as is, as we are copying the constraint later regardless
       * if there are variables in it or not
       */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, ncurvars) );
      SCIP_CALL( SCIPgetVarsXXX(scip, conss[i], vars, ncurvars) );

      /** @todo skip all variables that have a zero coeffient or where all coefficients add to zero */
      /** @todo Do more then one entry per variable actually work? */

      /* allocate a hyperedge for the constraint */
      SCIP_CALL( SCIPallocMemory(scip, &hedge) );
      hedge->type = CONSTRAINT;
      hedge->originalId = i;

      /* lets collect the variable ids of the variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &varids, ncurvars) );
      hedge->nvariableIds = 0;

      for( j = 0; j < ncurvars; ++j )
      {
         SCIP_VAR* var;
         int varIndex;

         /* if the variable is inactive, skip it */
         if( !SCIPisVarRelevant(vars[j]) )
            continue;

         var = SCIPvarGetProbvar(vars[j]);
         assert(var != NULL);
         varIndex = SCIPvarGetProbindex(var);
         /* assert that the variable is active and not multiaggregated, otherwise, the mapping will be wrong */
         /* the multiaggregation is useless, if we don't presolve, it might be interesting otherwise */
         assert(SCIPvarIsActive(var));
         assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

         /* add the variable id to the end of the variable id array and update the size */
         varids[hedge->nvariableIds] = id;
         ++(hedge->nvariableIds);

         /* put the copied id (started from 0, should be the index of the nonzero entry) to the end of the map of original to copy ids  */
         SCIP_CALL( SCIPsetIntarrayVal(scip, maporigtocopies[varIndex], nmaporigtocopies[varIndex], id) );
         ++(nmaporigtocopies[varIndex]);
         SCIPdebugMessage("Adding %d at %d to copytoorig.\n", varIndex, id);
         SCIP_CALL( SCIPsetIntarrayVal(scip, copytoorig, id, varIndex) );
         ++id;
         /* Check the mapping here */
#ifdef SCIP_DEBUG
         {
            int k;
            SCIPdebugMessage("Cons %s (%d): ", SCIPconsGetName(conss[i]), i);
            SCIPdebugPrintf("Var %s (%d): ", SCIPvarGetName(var), varIndex);
            for( k = 0; k < nmaporigtocopies[varIndex]; ++k )
            {
               int copy;
               int orig;
               copy = SCIPgetIntarrayVal(scip, maporigtocopies[varIndex], k);
               orig = SCIPgetIntarrayVal(scip, copytoorig, copy);
               SCIPdebugPrintf("%d (%d), ", copy+1, orig);
               assert(varIndex == orig);
            }
            SCIPdebugPrintf("\n");
         }
#endif
      }

      if( hedge->nvariableIds > 1 )
      {
         /* if the hyperedge contains more then 0 variables, add it to the end */
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &hedge->variableIds, varids, hedge->nvariableIds) );
         SCIP_CALL( SCIPsetPtrarrayVal(scip, hedges, nhyperedges, hedge) );
         ++nhyperedges;
      }
      else
      {
         SCIPfreeMemory(scip, &hedge);
      }
      SCIPfreeBufferArray(scip, &varids);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* build variable hyperedges */
   for( i = 0; i < nvars; i++ )
   {
      /*
       * handle variables only appearing in the objective
       * this is in a sense useless, as it increases the work for metis, but
       * it is consistent this way
       */
      int size;
      size = nmaporigtocopies[i];
      assert( size >= 0);

      /* if the hyperedge contains only 1 variable or less, skip it */
      if( size <= 1 )
         continue;


      SCIP_CALL( SCIPallocMemory(scip, &hedge) );
      hedge->type = VARIABLE;
      hedge->originalId = i;
      hedge->nvariableIds = 0;
      hedge->variableIds = NULL;

      SCIP_CALL( SCIPsetPtrarrayVal(scip, hedges, nhyperedges, hedge) );
      ++nhyperedges;

      /* link the copies together */
      SCIP_CALL( SCIPallocMemoryArray(scip, &hedge->variableIds, size) );
      hedge->nvariableIds = size;

      SCIPdebugMessage("nvars hedge: ");

      for( j = 0; j < size; ++j )
      {
         hedge->variableIds[j] = SCIPgetIntarrayVal(scip, maporigtocopies[i], j);
         SCIPdebugPrintf("%d, ", hedge->variableIds[j]+1);
      }
      SCIPdebugPrintf("\n");
      SCIP_CALL( SCIPfreeIntarray(scip, &maporigtocopies[i]) );
   }
   SCIPfreeMemoryArray(scip, &maporigtocopies);
   SCIPfreeMemoryArray(scip, &nmaporigtocopies);
   SCIPfreeMemoryArray(scip, &copies);
   nvertices = id;
   readerdata->nvertices = nvertices;

   assert(SCIPgetPtrarrayMaxIdx(scip, hedges)+1 == nhyperedges);
   //   assert(nvertices == SCIPgetIntarrayMaxIdx(scip, copytoorig)+1);


   return SCIP_OKAY;
}

/** reads cluster from file */
static
SCIP_RETCODE readClusterFromFile(
   SCIP*                 scip,               /**< SCIP data struture */
   SCIP_READERDATA*      readerdata,         /**< presolver data data structure */
   const char*           inputfile,          /**< input file */
   SCIP_RESULT*          result              /**< result indicating whether the detection was successful */
   )
{
   char line[SCIP_MAXSTRLEN];

   int i;
   int nvertices;
   int* partition;


   SCIP_FILE *zfile;
   SCIP_Real remainingtime;

   assert(scip != NULL);
   assert(readerdata != NULL);
   assert(inputfile != NULL);

   *result = SCIP_DIDNOTRUN;

   remainingtime = DECgetRemainingTime(scip);
   nvertices = readerdata->nvertices;

   if( remainingtime <= 0 )
   {
      return SCIP_OKAY;
   }

   if( readerdata->partition == NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &readerdata->partition, readerdata->nvertices) );
   }


   assert(readerdata->partition != NULL);
   partition = readerdata->partition;

   zfile = SCIPfopen(inputfile, "r");
   i = 0;
   while( !SCIPfeof(zfile) && i < nvertices )
   {
      int temp;
      int t1;
      if( SCIPfgets(line, SCIP_MAXSTRLEN, zfile) == NULL )
      {
         SCIPerrorMessage("Line could not be read\n");
         return SCIP_READERROR;
      }

      sscanf(line, "%d %d", &t1, &temp);
      assert(t1 == i+1 || t1 == i);
      /*      temp = atoi(line); */
      assert(temp >= 0);
      readerdata->blocks = MAX(temp+1, readerdata->blocks);
      partition[i] = temp;
      SCIPdebugMessage("%d: %d\n", i, temp);
      i++;
   }

   if( i != nvertices )
   {
      SCIPerrorMessage("Couldn't read partition for all vertices.\n");
      return SCIP_READERROR;
   }

   SCIPfclose(zfile);

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** maps the partitions for the disaggregated vertices to the original vertices */
static
SCIP_RETCODE assignBlocksToOriginalVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*     readerdata        /**< presolver data data structure */
   )
{

   int i;
#ifndef NDEBUG
   int nvars;
#endif
   int *partition;
   int *origpart;
   int nvertices;
   assert(scip != NULL);
   assert(readerdata != NULL);

   nvertices = readerdata->nvertices;
   partition = readerdata->partition;
   origpart = readerdata->varpart;

#ifndef NDEBUG
   nvars = SCIPgetNVars( scip );
#endif

    /* go through the new vertices */
   for( i = 0; i < nvertices ; ++i )
   {
      int originalId;
      /* find out the original id (== index of the var in the vars array) */
      originalId = SCIPgetIntarrayVal(scip, readerdata->copytooriginal, i);

      /* add the id to the set of ids for the original vertex */
      assert(originalId >= 0 && originalId < nvars);
      assert(partition[i] >= 0);
      if( origpart[originalId] == -1 )
      {
         origpart[originalId] = partition[i];
      }
      else if( origpart[originalId] >= 0 )
      {
         if( origpart[originalId] != partition[i] )
         {
            origpart[originalId] = -2;
         }
      }
      assert(origpart[originalId] == -2 || origpart[originalId] >= 0);
      assert(origpart[originalId] <= readerdata->blocks);
   }

   return SCIP_OKAY;
}

/** builds the transformed problem in the new scip instance */
static SCIP_RETCODE buildTransformedProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< presolver data data structure */
   DEC_DECOMP*           decomp,             /**< decomp data structure */
   int                   nblocks,            /**< number of blocks for this decomposition */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_Bool *isVarHandled;
   SCIP_VAR*** subscipvars;
   SCIP_VAR** linkingvars;
   SCIP_CONS*** subscipconss;
   SCIP_CONS** linkingconss;
   SCIP_HASHMAP* vartoblock;
   SCIP_HASHMAP* constoblock;

   int *nsubscipvars;
   int *nsubscipconss;
   int nlinkingvars;
   int nlinkingconss;
   int i;
   int j;
   SCIP_CONS **conss;
   int nconss;
   SCIP_VAR **vars;
   int nvars;
   SCIP_Bool emptyblocks = FALSE;

   assert(scip != NULL);
   assert(readerdata != NULL);

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

   SCIP_CALL( SCIPallocBufferArray(scip, &linkingvars, nvars) );
   nlinkingvars = 0;

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
         SCIP_CALL( SCIPgetVarsXXX( scip, conss[i], curvars, ncurvars) );
      }

      for( j = 0; j < ncurvars; j++ )
      {
         SCIP_VAR* var;
         int varblock = -1;
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
            assert(readerdata->varpart[SCIPvarGetProbindex(var)] < readerdata->blocks);
            /* get the number of blocks the current variable is in */
            assert(readerdata->varpart[SCIPvarGetProbindex(var)] == -2 ||
                   readerdata->varpart[SCIPvarGetProbindex(var)] >= 0);
            /* if the variable is in exactly one block */
            if( readerdata->varpart[SCIPvarGetProbindex(var)] != -2 )
            {
               /* then the partition is given */
               varblock = readerdata->varpart[SCIPvarGetProbindex(var)];
               assert(varblock < readerdata->blocks);
               subscipvars[varblock][nsubscipvars[varblock]] = var;
               //               SCIPdebugMessage("v: %s\n", SCIPvarGetName(var));
               ++(nsubscipvars[varblock]);
            }
            /*
             * if the variable is a linking variable, don't update the constraint
             * block and add the variable to the linking variables
             */
            else
            {
               varblock = readerdata->blocks+1;
               //               SCIPdebugMessage("v: %s\n", SCIPvarGetName(var));
               linkingvars[nlinkingvars] = var;
               ++nlinkingvars;
            }

            /* finally set the hashmap image */
            assert(!SCIPhashmapExists(vartoblock, var));
            SCIP_CALL( SCIPhashmapInsert(vartoblock, var, (void*) (size_t) varblock) );
         }
         else
         {
            varblock = (int)(size_t)SCIPhashmapGetImage(vartoblock, var);
            assert(varblock == readerdata->varpart[SCIPvarGetProbindex(var)] ||  readerdata->varpart[SCIPvarGetProbindex(var)] == -2);
         }


         /*
          * if the variable is not a linking variable, add it to the correct
          * block and update the block of the constraint, if necessary
          */
         if( varblock <= readerdata->blocks )
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

         block = readerdata->blocks +1;
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
      if( readerdata->varpart[i] >= 0 )
      {
         partitionOfVar = readerdata->varpart[i];
      }

      if( partitionOfVar != -1 )
      {
         subscipvars[partitionOfVar][nsubscipvars[partitionOfVar]] = SCIPvarGetProbvar(vars[i]);
         //         SCIPdebugMessage("v: %s\n", SCIPvarGetName(SCIPvarGetProbvar(vars[i])));
         ++nsubscipvars[partitionOfVar];
      }
      else
      {
         //         SCIPdebugMessage("v: %s\n", SCIPvarGetName(SCIPvarGetProbvar(vars[i])));
         linkingvars[nlinkingvars] = SCIPvarGetProbvar(SCIPvarGetProbvar(vars[i]));
         ++nlinkingvars;
      }

   }

   SCIPfreeBufferArray(scip, &isVarHandled);

   /* first, make sure that there are constraints in every block, otherwise the hole thing is useless */
   for( i = 0; i < readerdata->blocks; ++i )
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
      DECdecompSetNBlocks(decomp, nblocks);
      DECdecompSetType(decomp, DEC_DECTYPE_BORDERED);
      SCIP_CALL( DECdecompSetSubscipvars(scip, decomp, subscipvars, nsubscipvars) );
      SCIP_CALL( DECdecompSetSubscipconss(scip, decomp, subscipconss, nsubscipconss) );
      if( nlinkingconss > 0 )
      {
         SCIP_CALL( DECdecompSetLinkingconss(scip, decomp, linkingconss, nlinkingconss) );
         DECdecompSetType(decomp, DEC_DECTYPE_BORDERED);
      }
      if( nlinkingvars > 0 )
      {
         DECdecompSetType(decomp, DEC_DECTYPE_ARROWHEAD);
         SCIP_CALL( DECdecompSetLinkingvars(scip, decomp, linkingvars, nlinkingvars) );
      }
      DECdecompSetVartoblock(decomp, vartoblock);
      DECdecompSetConstoblock(decomp, constoblock);
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
   SCIPfreeBufferArray(scip, &linkingvars);

   *result = emptyblocks? SCIP_DIDNOTFIND:SCIP_SUCCESS;
   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeCluster)
{
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   /* free decomp structure and readerdata */
   if( DECdecompGetType(readerdata->decomp) == DEC_DECTYPE_UNKNOWN )
      DECdecompFree(scip, &readerdata->decomp);
   SCIPfreeMemory(scip, &readerdata);

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCluster)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPreadCluster(scip, filename, result) );

   return SCIP_OKAY;
}

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteCluster)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);

   /* SCIP_CALL( SCIPwriteCluster(scip, file, TRUE) ); */
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the cluster file reader in SCIP */
SCIP_RETCODE
SCIPincludeReaderCluster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create cluster reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );

   /* include cluster reader */
   SCIP_CALL(SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION, NULL,
           readerFreeCluster, readerReadCluster, readerWriteCluster, readerdata));

   return SCIP_OKAY;
}

/* reads problem from file */
SCIP_RETCODE SCIPreadCluster(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   assert(scip != NULL);
   assert(filename != NULL);
   assert(result != NULL);

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);

   SCIP_CALL( initReaderdata(scip, readerdata) );

   SCIP_CALL( buildGraphStructure(scip, readerdata) );

   SCIP_CALL( readClusterFromFile(scip, readerdata, filename, result) );

   SCIP_CALL( assignBlocksToOriginalVariables(scip, readerdata) );

   SCIP_CALL( buildTransformedProblem(scip, readerdata, readerdata->decomp, readerdata->blocks, result) );

   SCIP_CALL( SCIPconshdlrDecompAddDecdecomp(scip, readerdata->decomp) );

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}
