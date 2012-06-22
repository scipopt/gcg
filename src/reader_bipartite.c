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
/**@file   reader_bipartite.c
 * @brief  BIPARTITE file reader
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

#include "reader_bipartite.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"
#include "tclique/tclique.h"

#include "cons_decomp.h"
#include "pub_decomp.h"

#define READER_NAME             "bipartitereader"
#define READER_DESC             "file reader for blocks in bipartite format"
#define READER_EXTENSION        "bip"

#define DEFAULT_VARWEIGHT         1          /**< weight for variable nodes */
#define DEFAULT_VARWEIGHTBIN      2          /**< weight for binary variable nodes */
#define DEFAULT_VARWEIGHTINT      2          /**< weight for integer variable nodes */
#define DEFAULT_VARWEIGHTCONT     1          /**< weight for continous variable nodes */
#define DEFAULT_VARWEIGHTIMPL     2          /**< weight for implicit integer variable nodes */
#define DEFAULT_CONSWEIGHT        5          /**< weight for constraint hyperedges */


#define TCLIQUE_CALL(x)   do                                                                                  \
                       {                                                                                      \
                          SCIP_Bool _restat_;                                                                 \
                          if( (_restat_ = (x)) != TRUE )                                                      \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             return SCIP_ERROR;                                                               \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

/*
 * Data structures
 */

/** data for dec reader */
struct SCIP_ReaderData
{
   DEC_DECOMP*    decomp;         /**< decomposition data structure*/
   TCLIQUE_GRAPH* graph;          /**< graph */
   int            nblocks;        /**< number of blocks */
   int*           partition;      /**< the partition of the graph */

   /* weight parameters */
   int       varWeight;             /**< weight of a variable vertex */
   int       varWeightBinary;       /**< weight of a binary variable vertex */
   int       varWeightContinous;    /**< weight of a continuous variable vertex */
   int       varWeightInteger;      /**< weight of an integer variable vertex */
   int       varWeightImplint;      /**< weight of an implicit integer variable vertex */
   int       consWeight;            /**< weight of a constraint vertex */
};

/** calculates weights for variables */
static
TCLIQUE_WEIGHT calculateVarWeight(
   SCIP_READERDATA*      readerdata,         /**< reader data structure*/
   SCIP_VAR*             var                 /**< variable for which the weight is computed */
   )
{
   int varweight;

   switch ( SCIPvarGetType(var) ) {
   case SCIP_VARTYPE_CONTINUOUS:
      varweight = readerdata->varWeightContinous;
      break;
   case SCIP_VARTYPE_INTEGER:
      varweight = readerdata->varWeightInteger;
      break;
   case SCIP_VARTYPE_IMPLINT:
      varweight = readerdata->varWeightImplint;
      break;
   case SCIP_VARTYPE_BINARY:
      varweight = readerdata->varWeightBinary;
      break;
   default:
      varweight = readerdata->varWeight;
      break;
   }

   return varweight;
}

/** calculates weights for constraints */
static
TCLIQUE_WEIGHT calculateConsWeight(
   SCIP_READERDATA*      readerdata,         /**< reader data structure*/
   SCIP_CONS*            cons                /**< constraint for which the weight is computed */
   )
{
   /*lint --e{715} */
   return  readerdata->consWeight;
}



/** initializes reader data structure */
static SCIP_RETCODE initReaderdata(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_READERDATA*      readerdata          /**< reader data structure */
   )
{
   assert(scip != NULL);
   assert(readerdata != NULL);

   SCIP_CALL( DECdecompCreate(scip, &readerdata->decomp) );
   TCLIQUE_CALL( tcliqueCreate(&readerdata->graph) );
   readerdata->nblocks = 0;

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
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   int nconss = SCIPgetNConss( scip );
   int nvars = SCIPgetNVars( scip );
   int i;
   int j;

   assert(scip != NULL);
   assert(readerdata != NULL);

   conss = SCIPgetConss(scip);
   vars = SCIPgetVars(scip);

   for( i = 0; i < nvars + nconss; ++i )
   {
      TCLIQUE_WEIGHT weight;

      /* note that the first nvars nodes correspond to variables */
      if( i < nvars )
         weight = calculateVarWeight(readerdata, vars[i]);
      else
         weight = calculateConsWeight(readerdata, conss[i]);

      TCLIQUE_CALL( tcliqueAddNode(readerdata->graph, i, weight) );
   }

   /* go through all constraints */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_VAR **curvars;
      /* get the number of nonzeros in this constraint  */

      int ncurvars = SCIPgetNVarsXXX( scip, conss[i] );

      /* if there are no variables, skip the constraint */
      if( ncurvars == 0 )
         continue;

      /*
       * may work as is, as we are copying the constraint later regardless
       * if there are variables in it or not
       */
      SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
      SCIP_CALL( SCIPgetVarsXXX(scip, conss[i], curvars, ncurvars) );

      /** @todo skip all variables that have a zero coeffient or where all coefficients add to zero */
      /** @todo Do more then one entry per variable actually work? */

      /* allocate a hyperedge for the constraint */

      for( j = 0; j < ncurvars; ++j )
      {
         SCIP_VAR* var;
         int varIndex;

         /* if the variable is inactive, skip it */
         if( !SCIPisVarRelevant(curvars[j]) )
            continue;

         var = SCIPvarGetProbvar(curvars[j]);

         assert(var != NULL);
         varIndex = SCIPvarGetProbindex(var);

         TCLIQUE_CALL( tcliqueAddEdge(readerdata->graph, varIndex, nvars+i) );
      }
      SCIPfreeBufferArray(scip, &curvars);
   }

   TCLIQUE_CALL( tcliqueFlush(readerdata->graph) );

   return SCIP_OKAY;
}

/** reads bipartite from file */
static
SCIP_RETCODE readBipartiteFromFile(
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
   nvertices = tcliqueGetNNodes(readerdata->graph);

   if( remainingtime <= 0 )
   {
      return SCIP_OKAY;
   }

   if( readerdata->partition == NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &readerdata->partition, nvertices) );
   }


   assert(readerdata->partition != NULL);
   partition = readerdata->partition;

   zfile = SCIPfopen(inputfile, "r");
   i = 0;
   while( !SCIPfeof(zfile) && i < nvertices )
   {
      int temp;
      if( SCIPfgets(line, SCIP_MAXSTRLEN, zfile) == NULL )
      {
         SCIPerrorMessage("Line could not be read\n");
         return SCIP_READERROR;
      }

      sscanf(line, "%d", &temp);
      assert(temp >= 0);

      readerdata->nblocks = MAX(temp+1, readerdata->nblocks);
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
   readerdata->nblocks -= 1;
   *result = SCIP_SUCCESS;
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
      SCIP_CALL( SCIPallocBufferArray(scip, &(subscipconss[i]), nconss) ); /*lint !e866 */
      SCIP_CALL( SCIPallocBufferArray(scip, &(subscipvars[i]), nvars) );   /*lint !e866 */

      nsubscipconss[i] = 0;
      nsubscipvars[i] = 0;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &linkingconss, nconss) );
   nlinkingconss = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &linkingvars, nvars) );
   nlinkingvars = 0;

   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPhashmapCreate(&vartoblock, SCIPblkmem(scip), nconss) );

   /* go through all of the constraints */
   for( i = 0; i < nconss; i++ )
   {
      int conspart;
      int consindex;

      consindex = nvars+i;

      assert(consindex > nvars && consindex < nconss+nvars);
      conspart = readerdata->partition[consindex];

      if( conspart > -1 && conspart < readerdata->nblocks)
      {
         subscipconss[conspart][nsubscipconss[conspart]] = conss[i];
         ++nsubscipconss[conspart];
      }
      else
      {
         linkingconss[nlinkingconss] = conss[i];
         ++nlinkingconss;
      }
   }

   /* go through all variables */

   for( i = 0; i < nvars; i++ )
   {
      int varpart;
      varpart = readerdata->partition[i];

      if( !SCIPisVarRelevant(vars[i]) )
         continue;

      if( varpart > -1 && varpart < readerdata->nblocks )
      {
         subscipvars[varpart][nsubscipvars[varpart]] = vars[i];
         ++nsubscipvars[varpart];
      }
      else
      {
         linkingvars[nlinkingvars] = vars[i];
         ++nlinkingvars;
      }

   }

   /* first, make sure that there are constraints in every block, otherwise the hole thing is useless */
   for( i = 0; i < readerdata->nblocks; ++i )
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
      DECdecompSetType(decomp, DEC_DECTYPE_DIAGONAL);
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
   for( i = nblocks-1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &(subscipvars[i]));
      SCIPfreeBufferArray(scip, &(subscipconss[i]));
   }

   SCIPfreeBufferArray(scip, &nsubscipvars);
   SCIPfreeBufferArray(scip, &nsubscipconss);
   SCIPfreeBufferArray(scip, &subscipvars);
   SCIPfreeBufferArray(scip, &subscipconss);

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
SCIP_DECL_READERFREE(readerFreeBipartite)
{
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   /* free decomp structure and readerdata */
   if( DECdecompGetType(readerdata->decomp) == DEC_DECTYPE_UNKNOWN )
      DECdecompFree(scip, &readerdata->decomp);
   tcliqueFree(&readerdata->graph);
   SCIPfreeMemory(scip, &readerdata);

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadBipartite)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPreadBipartite(scip, filename, result) );

   return SCIP_OKAY;
}

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteBipartite)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   SCIP_CALL( SCIPtransformProb(scip) );
   SCIP_CALL( SCIPwriteBipartite(scip, file) );
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the bipartite file reader in SCIP */
SCIP_RETCODE
SCIPincludeReaderBipartite(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create bipartite reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );

   /* include bipartite reader */
   SCIP_CALL(SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION, NULL,
           readerFreeBipartite, readerReadBipartite, readerWriteBipartite, readerdata));

   SCIP_CALL( SCIPaddIntParam(scip, "reader/bipartite/varWeight", "Weight of a variable hyperedge", &readerdata->varWeight, FALSE, DEFAULT_VARWEIGHT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "reader/bipartite/varWeightBinary", "Weight of a binary variable hyperedge", &readerdata->varWeightBinary, FALSE, DEFAULT_VARWEIGHTBIN, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "reader/bipartite/varWeightContinous", "Weight of a continuos variable hyperedge", &readerdata->varWeightContinous, FALSE, DEFAULT_VARWEIGHTCONT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "reader/bipartite/varWeightImplint", "Weight of a implicit integer variable hyperedge", &readerdata->varWeightImplint, FALSE, DEFAULT_VARWEIGHTIMPL, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "reader/bipartite/varWeightInteger", "Weight of a integer variable hyperedge", &readerdata->varWeightInteger, FALSE, DEFAULT_VARWEIGHTINT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "reader/bipartite/consWeight", "Weight of a constraint hyperedge", &readerdata->consWeight, FALSE, DEFAULT_CONSWEIGHT, 0, 1000000, NULL, NULL) );

   return SCIP_OKAY;
}

/* reads problem from file */
SCIP_RETCODE SCIPreadBipartite(
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

   SCIP_CALL( readBipartiteFromFile(scip, readerdata, filename, result) );

   SCIP_CALL( buildTransformedProblem(scip, readerdata, readerdata->decomp, readerdata->nblocks, result) );

   SCIP_CALL( SCIPconshdlrDecompAddDecdecomp(scip, readerdata->decomp) );

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/** writes problem to file */
SCIP_RETCODE SCIPwriteBipartite(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< file pointer to the file where the problem is written to */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   int nnodes;
   int nedges;
   int nvars;
   int i;
   const TCLIQUE_WEIGHT* weights;

   assert(scip != NULL);
   assert(file != NULL);

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);

   SCIP_CALL( initReaderdata(scip, readerdata) );

   SCIP_CALL( buildGraphStructure(scip, readerdata) );

   nnodes = tcliqueGetNNodes(readerdata->graph);
   nedges = tcliqueGetNEdges(readerdata->graph);
   weights = tcliqueGetWeights(readerdata->graph);
   nvars = SCIPgetNVars(scip);

   assert(nedges % 2 == 0);

   /* write out graph */
   SCIPinfoMessage(scip, file, "%d %d 10 2\n", nnodes, nedges/2);

   for( i = 0; i < nnodes; ++i )
   {
      int* firstedge;
      int* lastedge;
      firstedge = tcliqueGetFirstAdjedge(readerdata->graph, i);
      lastedge = tcliqueGetLastAdjedge(readerdata->graph, i);

      SCIPinfoMessage(scip, file, "%d %d", weights[i], i < nvars? 0:1);
      while( firstedge <= lastedge )
      {
         SCIPinfoMessage(scip, file, " %d", *firstedge+1);
         ++firstedge;
      }
      SCIPinfoMessage(scip, file, "\n");
   }

   return SCIP_OKAY;
}
