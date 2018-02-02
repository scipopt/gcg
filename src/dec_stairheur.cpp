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

/**@file   dec_stairheur.cpp
 * @brief  detector for staircase structures via ROC algorithms
 * @author Martin Bergner
 * @author Mathias Luers
 * @ingroup DETECTORS
 *
 * This detector is based on Jayakumar, Maliyakal D., and Ranga V. Ramasesh.
 * A clustering heuristic to detect staircase structures in large scale
 * linear programming models. European journal of operational research 76.1 (1994): 229-239.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <cstring>
#include <algorithm>
#include <vector>

#include "dec_stairheur.h"
#include "cons_decomp.h"
#include "struct_decomp.h"
#include "pub_decomp.h"
#include "scip_misc.h"
#include "scip/pub_misc.h"
#include "scip/struct_var.h"

#define DEC_DETECTORNAME      "stairheur"    /**< name of the detector */
#define DEC_DESC              "detects staircase matrices via matrix reordering" /**< detector description */
#define DEC_PRIORITY          1200           /**< priority of the detector */
#define DEC_DECCHAR           's'            /**< display character of detector */
#define DEC_ENABLED           FALSE          /**< should detector be called by default */
#define DEC_SKIP              FALSE          /**< should detector be skipped if others found detections */

/* Default parameter settings*/
#define DEFAULT_NCONSSPERBLOCK               32       /**< number of constraints per block (static blocking only) */
#define DEFAULT_MAXBLOCKS                    20       /**< value for the maximum number of blocks to be considered */
#define DEFAULT_MINBLOCKS                    2        /**< value for the minimum number of blocks to be considered */
#define DEFAULT_DESIREDBLOCKS                0        /**< value for the desired number of blocks (for all blocking types). Set to zero for self determination of block number */
#define DEFAULT_DYNAMICBLOCKING              FALSE    /**< Enable blocking type 'dynamic' */
#define DEFAULT_STATICBLOCKING               TRUE     /**< Enable blocking type 'static' */
#define DEFAULT_BLOCKINGASSOONASPOSSIBLE     FALSE    /**< Enable blocking type 'as soon as possible' */
#define DEFAULT_MULTIPLEDECOMPS              TRUE     /**< Enables multiple decompositions for all enabled blocking types. Ranging from minblocks to maxblocks' */
#define DEFAULT_MAXITERATIONSROC             1000000  /**< The maximum of iterations of the ROC-algorithm. -1 for no iteration limit */

using std::find;
using std::vector;
using std::swap;

#ifdef WRITEALLOUTPUT
#define DWSOLVER_REFNAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_ref.txt", (name), (blocks), (cons), (dummy)
#define GP_NAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_%d.gp", (name), (blocks), (cons), (dummy)
#endif
/*
 * Data structures
 */

/** A struct that contains 4 hashmaps, which maps variables and constraints to their position in the constraint matrix (Ax<=b) and vice versa */
struct IndexMap
{
   SCIP_HASHMAP*         indexcons;          /**< index in problem -> constraint */
   SCIP_HASHMAP*         consindex;          /**< constraint -> index in problem */
   SCIP_HASHMAP*         indexvar;           /**< index in problem -> variable */
   SCIP_HASHMAP*         varindex;           /** variable -> index in problem */
};
typedef struct IndexMap INDEXMAP;

/** detector data */
struct DEC_DetectorData
{
   SCIP_HASHMAP*         constoblock;        /**< hashmap mapping constraints to blocks */
   int                   blocks;             /**< number of blocks */
   int                   nconssperblock;     /**< number of constraints per block (static blocking only) */
   int                   maxblocks;          /**< maximum number of constraints per block */
   int                   minblocks;          /**< minimum number of constraints per block */
   INDEXMAP*             indexmap;           /**< index map (contains 4 hashmaps) */
   int*                  ibegin;             /**< array, ibegin[i]: index of first nonzero entry in row i */
   int*                  iend;               /**< array, iend[i]: index of last nonzero entry in row i */
   int*                  jbegin;             /**< array, jbegin[j]: index of first nonzero entry in column j */
   int*                  jend;               /**< array, jend[j]: index of last nonzero entry in column j */
   int*                  jmin;               /**< array, jmin[i]: index of first nonzero column of the i-th row */
   int*                  jmax;               /**< array, jmax[i]: the last nonzero entry among all rows prior to and including the i-th row */
   int*                  minV;               /**< array, minV[i]: number of linking variables corresponding to a partitioning after the i-th row */
   int*                  width;              /**< array, width[i]: width of the band (of nonzero entries after ROC) at row i */
   int*                  hashmapindices;     /**< array with integers running from 0 to maximum(nvars, ncons)+1 (for usage of hash maps) */
   vector<int>*          rowsWithConstrictions;
   vector<int>*          blockedAfterrow;
   int                   desiredblocks;
   SCIP_Bool             dynamicblocking;    /**< Enable blocking type 'dynamic' */
   SCIP_Bool             staticblocking;     /**< Enable blocking type 'static' */
   SCIP_Bool             blockingassoonaspossible; /**< Enable blocking type 'as soon as possible' */
   SCIP_Bool             multipledecomps;    /**< Enables multiple decompositions for all enabled blocking types. Ranging from minblocks to maxblocks */
   int                   maxiterationsROC;
};

/*
 * Local methods
 */


/* debugging methods */
void printnested(
   vector<vector<int> >  l
   )
{
   vector<int>::iterator inner;
   vector<vector<int> >::iterator outer;
   SCIPdebugPrintf("S:");
   for( outer = l.begin(); outer != l.end(); ++outer)
   {
      SCIPdebugPrintf("\t");
      for( inner = outer->begin(); inner != outer->end(); ++inner)
      {
         SCIPdebugPrintf(" %d", *inner);
      }
      SCIPdebugPrintf(".\n");
   }
   SCIPdebugPrintf("Done\n");
}

void printvector(
   vector<int >          l
   )
{
   vector<int>::iterator inner;
   for( inner = l.begin(); inner != l.end(); ++inner)
   {
      SCIPdebugPrintf(" %d", *inner);
   }
}


/** TODO:
 * currently, all vars from the first column where a linking var appears until the end of the block are considered as linking vars, although there might be empty columns. This could be changed so that these empty columns are considered as subscipvars and not linking vars.
 *
 * In some cases a block can consist of linking vars exclusively. This makes no real sense.
 *
 * For some instances the assertion regarding the consistency of the arrays ibegin and jbegin fails
 * */


/** creates a list with integers running from 'from' to 'to'. */
static
vector<int> SCIPvectorCreateInt(
   int                   from,               /**< Start index */
   int                   to                  /**< End index */
   )
{
   vector<int> v;

   for( int i = from; i <= to; ++i )
   {
      v.push_back(i);
   }
   return v;
}


/** rearranges elements of vector according to the ordering of order.
 *
 * example: vector = (a b c d); order = (3 2 4 1)
 * after calling SCIPvectorRearrange(vector, order): vector = (c b d a)
 * both vectors must have the same size
 * order must have elements from 1 to vector->size */
static
void SCIPvectorRearrange(
   vector<vector<int> >  &l,
   vector<int>           order)
{

   vector<vector<int> > new_vector;
   vector<int>::iterator it1;

   for( it1 = order.begin(); it1 != order.end(); ++it1 )
   {
      new_vector.push_back(l[(size_t)*it1-1]);

   }
   l.swap(new_vector);

}

/** allocates memory for an indexmap. */
static
SCIP_RETCODE indexmapCreate(
   SCIP*                 scip,               /**< SCIP data structure  */
   INDEXMAP**            indexmap,           /**< address of the pointer which shall store the index map*/
   int                   nconss,             /**< number of constraints */
   int                   nvars               /**< number of variables */
)
{
   INDEXMAP* imap = NULL;
   assert(scip != NULL);
   assert(nconss > 0);
   assert(nvars > 0);
   SCIP_CALL( SCIPallocMemory(scip, &imap) );
   assert(imap != NULL);

   SCIP_CALL( SCIPhashmapCreate(&imap->indexvar, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&imap->varindex, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&imap->indexcons, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPhashmapCreate(&imap->consindex, SCIPblkmem(scip), nconss) );

   *indexmap = imap;
   return SCIP_OKAY;
}

/** deallocates memory of indexmap. */
static
void indexmapFree(
   SCIP*                 scip,               /**< SCIP data structure */
   INDEXMAP**            indexmap            /**< index map */
   )
{
   SCIPhashmapFree(&(*indexmap)->indexvar);
   SCIPhashmapFree(&(*indexmap)->varindex);
   SCIPhashmapFree(&(*indexmap)->indexcons);
   SCIPhashmapFree(&(*indexmap)->consindex);
   SCIPfreeMemory(scip, indexmap);
}

/** initialization method for the indexmap */
static
SCIP_RETCODE indexmapInit(
   INDEXMAP*             indexmap,           /**< index map */
   SCIP_VAR**            vars,               /**< array of variables */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   int*                  hashmapindices      /**< indices of variables and constraints */
   )
{
   int i;
   int* hashmapindex;

   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* var = vars[i];
      /* careful: hashmapindex+1, because '0' is treated as an empty hashmap entry, which causes an error */
      hashmapindex = hashmapindices + i+1;
      assert( ! SCIPhashmapExists(indexmap->indexvar, (void*) hashmapindex));
      SCIP_CALL( SCIPhashmapInsert(indexmap->indexvar, (void*) hashmapindex, (void*) var) );
      assert( ! SCIPhashmapExists(indexmap->varindex, (void*) var));
      SCIP_CALL( SCIPhashmapInsert(indexmap->varindex, (void*) var, (void*) hashmapindex) );
   }

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CONS* cons = conss[i];
      /* careful: hashmapindex+1, because '0' is treated as an empty hashmap entry, which causes an error */
      hashmapindex = hashmapindices + i+1;
      assert( ! SCIPhashmapExists(indexmap->indexcons, (void*) hashmapindex));
      SCIP_CALL( SCIPhashmapInsert(indexmap->indexcons, (void*) hashmapindex, (void*) cons) );
      assert( ! SCIPhashmapExists(indexmap->consindex, (void*) cons));
      SCIP_CALL( SCIPhashmapInsert(indexmap->consindex, (void*) cons, (void*) hashmapindex) );
   }

   return SCIP_OKAY;
}

/* debug ? */
#ifdef WRITEALLOUTPUT
/** returns the problem name without the path */
static const char* getProbNameWithoutPath(
   SCIP*                 scip
)
{
   const char* pname;
   /* remove '/' from problem name */
   pname = strrchr(SCIPgetProbName(scip), '/');
   if( pname == NULL )
   {
      pname = SCIPgetProbName(scip);
   }
   else
   {
      pname = pname+1;
   }
   return pname;
}


static void checkConsistencyOfIndexarrays(DEC_DETECTORDATA* detectordata, int nvars, int nconss)
{
   int i;
   for( i = 0; i < nconss - 1; ++i )
   {
      assert(detectordata->ibegin[i] <= detectordata->ibegin[i+1]);
   }
   for( i = 0; i < nvars - 1; ++i )
   {
      assert(detectordata->jbegin[i] <= detectordata->jbegin[i+1]);
   }
}


/* debug ? */
/** creates a data and a gnuplot file for the initial problem.
 * @param scip < SCIP data structure
 * @param detectordata < detector data data structure
 * @param filename name of the output files (without any filename extension) */
static
SCIP_RETCODE plotInitialProblem(
   SCIP*                 scip,
   DEC_DETECTORDATA*     detectordata,
   char*                 filename
)
{
   FILE* output;
   char datafile[256];
   char gpfile[256];
   char pdffile[256];

   int nconss;
   nconss = SCIPgetNConss(scip);

   /* filenames */
   sprintf(datafile, "%s.dat", filename);
   sprintf(gpfile, "%s.gp", filename);
   sprintf(pdffile, "%s.pdf", filename);
   output = fopen(datafile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for output in plotProblem!\n");
   }
   else
   {
      int i;

      for( i = 0; i < nconss; ++i )
      {
         int j;
         SCIP_Bool success;
         SCIP_VAR** curvars;
         int ncurvars;
         int consindex;
         SCIP_CONS* cons = SCIPgetConss(scip)[i];
         consindex = *((int*)  SCIPhashmapGetImage(detectordata->indexmap->consindex, cons));
         assert(SCIPhashmapExists(detectordata->indexmap->consindex, cons));

         /* Get array of variables from constraint */
         SCIP_CALL( SCIPgetConsNVars(scip, cons, &ncurvars, &success) );
         assert(success);
         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetConsVars(scip, cons, curvars, ncurvars, &success) );
         assert(success);
         for( j = 0; j < ncurvars; ++j )
         {
            SCIP_VAR* var = SCIPvarGetProbvar(curvars[j]);
            assert(SCIPhashmapExists(detectordata->indexmap->varindex, var));
            int varindex = *((int*) SCIPhashmapGetImage(detectordata->indexmap->varindex, var));
            assert(varindex <= SCIPgetNVars(scip));
            assert(varindex > 0);
            assert(consindex <= SCIPgetNConss(scip));
            assert(consindex > 0);
            fprintf(output, "%d %d\n", varindex, consindex);
         }
         SCIPfreeBufferArray(scip, &curvars);
      }
   }
   fclose(output);

   /* write Gnuplot file */
   output = fopen(gpfile, "w");
   fprintf(output, "set terminal pdf\nset output \"%s\"\nunset xtics\nunset ytics\nunset border\nset pointsize 0.05\nset xrange [0:%i]\nset yrange[%i:0]\nplot '%s' lt 0 pt 5 notitle", pdffile, SCIPgetNVars(scip), nconss, datafile);
   fclose(output);
   return SCIP_OKAY;
}

/** creates a data and a gnuplot file for the graph representing the array minV (number of linking variables).
 * @param detectordata < detector data data structure
 * @param filename name of the output files (without any filename extension) */
static
void plotMinV(
   SCIP*                 scip,
   DEC_DETECTORDATA*     detectordata,
   char*                 filename
)
{
   FILE* output;
   char datafile[256];
   char blockingfile[256];
   char gpfile[256];
   char pdffile[256];
   int nconss;
   vector<int>::iterator it1;

   nconss = SCIPgetNConss(scip);

   /* filenames */
   sprintf(datafile, "%s.dat", filename);
   sprintf(blockingfile, "%s_blocked_at.dat", filename);
   sprintf(gpfile, "%s.gp", filename);
   sprintf(pdffile, "%s.pdf", filename);

   /* datafile */
   output = fopen(datafile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for output in plotMinV!\n");
   }
   else
   {
      int i;

      /* write data to datafile */
      for( i = 0; i < nconss -1; ++i )
      {
         fprintf(output, "%d\n", detectordata->minV[i]);
      }
   }
   fclose(output);

   /* blocking points */
   output = fopen(blockingfile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for blocking output in plotMinV!\n");
   }
   else
   {
      /* write data to blockingfile */
      for( it1 = detectordata->blockedAfterrow->begin(); it1 != detectordata->blockedAfterrow->end(); ++it1 )
      {
         fprintf(output, "%d %d\n", *it1-1, detectordata->minV[*it1-1]);
      }
   }
   fclose(output);
   /* write Gnuplot file */
   output = fopen(gpfile, "w");
   fprintf(output, "set terminal pdf\nset output \"%s\"\nset style line 1 lt 1 lc rgb \"black\"\nplot '%s' title '# verb. Variablen' ls 1 with lines, \\\n '%s' lt 0 pt 4 with points title \"Blockgrenze\"", pdffile, datafile, blockingfile);
   fclose(output);
}

#endif



/** initialization method of detector data */
static
SCIP_RETCODE initData(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< detector data structure */
   )
{
   int i;
   int nvars;
   int nconss;

   assert(scip != NULL);
   assert(detectordata != NULL);

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);

   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->ibegin, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->iend, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->jbegin, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->jend, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->jmin, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->jmax, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->minV, nconss-1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->width, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->hashmapindices, (size_t)MAX(nvars, nconss) + 1) );
   for( i = 0; i < MAX(nvars, nconss)+1; ++i )
   {
      detectordata->hashmapindices[i] = i;
   }
   detectordata->rowsWithConstrictions = new vector<int>();
   detectordata->blockedAfterrow = new vector<int>();

   /* create hash tables */
   SCIP_CALL( indexmapCreate(scip, &detectordata->indexmap, nconss, nvars) );

   return SCIP_OKAY;
}

/** deinitialization method of detector data (called after detection has been finished) */
static
SCIP_RETCODE freeData(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< detector data structure */
   )
{
   assert(scip != NULL);
   assert(detectordata != NULL);

   indexmapFree(scip, &detectordata->indexmap);

   /* delete vectors */
   delete detectordata->rowsWithConstrictions;
   detectordata->rowsWithConstrictions = NULL;
   delete detectordata->blockedAfterrow;
   detectordata->blockedAfterrow = NULL;

   SCIPfreeMemoryArray(scip, &detectordata->ibegin);
   SCIPfreeMemoryArray(scip, &detectordata->iend);
   SCIPfreeMemoryArray(scip, &detectordata->jbegin);
   SCIPfreeMemoryArray(scip, &detectordata->jend);
   SCIPfreeMemoryArray(scip, &detectordata->jmin);
   SCIPfreeMemoryArray(scip, &detectordata->jmax);
   SCIPfreeMemoryArray(scip, &detectordata->minV);
   SCIPfreeMemoryArray(scip, &detectordata->width);
   SCIPfreeMemoryArray(scip, &detectordata->hashmapindices);

   return SCIP_OKAY;
}

/** creates a nested vector with the indices of the nonzero entries of each row.
 *
 * example:
 * constraint matrix:
 *
 *  1 1 0 1 0
 *
 *  0 1 1 0 0
 *
 *  0 0 0 0 1
 *
 *  resulting vector:
 *  ( (1 2 4)
 *    (2 3)
 *    (5)    )
 */
static
SCIP_RETCODE createRowindexList(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data data structure */
   SCIP_HASHMAP*         indexcons,          /**< hashmap index -> constraint */
   SCIP_HASHMAP*         varindex,           /**< hashmap variable -> index*/
   vector<vector<int> >  &rowindices         /**< vector to store the row indices vector*/
      )
{
   int i;
   int nconss = SCIPgetNConss(scip);

   for( i = 0; i < nconss; ++i )
   {
      int j;
      SCIP_CONS* cons;
      SCIP_VAR** vars;
      int nvars;
      int* probindices = NULL;
      vector<int> rowindices_row;
      int* hashmapindex = &detectordata->hashmapindices[(size_t)i+1];

      cons = (SCIP_CONS*) SCIPhashmapGetImage(indexcons, (void*) hashmapindex);
      nvars = GCGconsGetNVars(scip, cons);

      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      SCIP_CALL( GCGconsGetVars(scip, cons, vars, nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &probindices, nvars) );

      /* fill the array with the indices of the variables of the current constraint */
      for( j = 0; j < nvars; ++j )
      {
         probindices[j] = *(int*) SCIPhashmapGetImage(varindex, SCIPvarGetProbvar(vars[j]));
      }

      std::sort(probindices, probindices+nvars);

      /* store a copy of the elements of probindices in the vector rowindices_row */
      for( j = 0; j < nvars; ++j )
      {
         rowindices_row.push_back(probindices[j]);
      }

      SCIPfreeMemoryArray(scip, &probindices);
      SCIPfreeBufferArray(scip, &vars);

      /* add rowindices_row to the vector rowindices */
      rowindices.push_back(rowindices_row);
      rowindices_row.clear();
   }

   return SCIP_OKAY;
}

/** creates a nested vector with the indices of the nonzero entries of each column.
 *
 * example:
 *
 * constraint matrix:
 *
 *  1 1 0 1 0
 *
 *  0 1 1 0 0
 *
 *  0 0 0 0 1
 *
 *  resulting vector:
 *  ( (1)
 *    (1 2)
 *    (2)
 *    (1)
 *    (3)    )
 */
static
SCIP_RETCODE createColumnindexList(
   SCIP*                 scip,               /**< SCIP data structure */
   vector<vector<int> >  &rowindices,        /**< A vector with the row indices (achieved from calling rowindices_vector() ) */
   vector<vector<int> >  &columnindices      /**< vector to store the column indices vector*/
)
{
   int position;
   int nvars;
   int i;
   vector<vector<int> >::iterator outer;
   vector<int>::iterator inner;
   nvars = SCIPgetNVars(scip);

   vector<vector<int> > columnindices_array((size_t)nvars);

   for( outer = rowindices.begin(), i = 1; outer != rowindices.end(); ++outer, ++i )
   {
      for( inner = outer->begin(); inner != outer->end(); ++inner )
      {
         position = (*inner)-1;
         columnindices_array[(size_t) position].push_back(i);
      }
   }

   /* create a columnindices vector instead of an array */
   for( i = 0; i < nvars; ++i )
   {
      columnindices.push_back(columnindices_array[(size_t)i]); /** @todo broken */
   }

   return SCIP_OKAY;
}

/** does the row ordering of the ROC2 algorithm.
 *
 * It also works for the column ordering. In this case the terms row<->column have to be exchanged.
 *
 * @return A vector with the new row order. E.g. (2 3 1) means the second row comes first now, and so on. */
static
vector<int> rowOrdering(
   SCIP*                 scip,               /**< SCIP data structure */
   vector<vector<int> >  &columnindices,     /**< A vector of the nonzero entries in each column */
   int                   nrows               /**< The number of rows of the constraint matrix (=number of relevant constraints) */
   )
{
   vector<int> roworder;
   vector<int> new_roworder;
   vector<vector<int> >::reverse_iterator it1;
   vector<int>::reverse_iterator it2;

   /* create a vector for the order of the rows ( 1 2 3 ... nrows ) */
   roworder = SCIPvectorCreateInt(1, nrows);
   new_roworder = SCIPvectorCreateInt(1, nrows);

   /* first from back to front */
   for( it1 = columnindices.rbegin(); it1 != columnindices.rend(); ++it1 )
   {

      /* second from back to front */
      for( it2 = it1->rbegin(); it2 != it1->rend(); ++it2 )
      {
         vector<int>::iterator tmp;

         tmp = std::find(new_roworder.begin(), new_roworder.end(), *it2); /*lint !e864*/
         std::rotate(new_roworder.begin(), tmp, tmp+1); /*lint !e1702 !e747*/
      }
      roworder = new_roworder;
   }

   return roworder;
}

/** stores the first and last entry of the i-th column(row) in begin[i] and end[i] respectively.
 *
 * @param begin Array to store the first nonzero entry of the i-th column (row)
 * @param end Array to store the last nonzero entry of the i-th column (row)
 * @param indices columnindices vector (rowindices vector) */
static
SCIP_RETCODE formIndexArray(
   int*                  begin,
   int*                  end,
   vector<vector<int> >  &indices
   )
{
   vector<vector<int> >::iterator it1;
   int i;
   assert(begin != NULL && end != NULL);
   for( it1 = indices.begin(), i = 0; it1 != indices.end(); ++it1, ++i )
   {
      /* case: vector not empty */
      if( !it1->empty() )
      {
         begin[i] = it1->front();
         end[i] = it1->back();
      }
      /* case: vector empty */
      else
      {
         begin[i] = 0;
         end[i] = 0;
      }
   }
   return SCIP_OKAY;
}


/**returns FALSE if at least one entry of new_array and old_array are different.*/
static
SCIP_Bool arraysAreEqual(
   int*                  new_array,          /**< new array */
   int*                  old_array,          /**< old array */
   int                   num_elements        /**< length of arrays */
   )
{
   int i;
   for( i = 0; i < num_elements; ++i )
   {
      if( new_array[i] != old_array[i] )
      {
         return FALSE;
      }
   }
   /* case: all entries of old and new are equal */
   return TRUE;
}

/**permutes the order of rows and columns in inputmap and stores the result in outputmap.
 *
 *  One call of this function is equivalent to one iteration of the ROC2-algortihm. */
static
SCIP_RETCODE rankOrderClusteringIteration(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data data structure */
   INDEXMAP*             inputmap,           /**< indexmap for input */
   INDEXMAP*             outputmap           /**< indexmap for output */
   )
{
   vector<int> roworder;
   vector<int> columnorder;
   vector<vector<int> > rowindices;
   vector<vector<int> > columnindices;
   vector<int>::iterator it1;
   int nvars;
   int ncons;
   size_t i;
   int position;
   int* hashmapindex;

   SCIPdebugMessage("Entering rankOrderClusteringIteration\n");

   assert(scip != NULL);
   assert(detectordata != NULL);
   nvars = SCIPgetNVars(scip);
   ncons = SCIPgetNConss(scip);

   /* create the vectors containing the positions of nonzero entries; row and column ordering */
   SCIP_CALL( createRowindexList(scip, detectordata, inputmap->indexcons, inputmap->varindex, rowindices) );
   SCIP_CALL( createColumnindexList(scip, rowindices, columnindices) );

   roworder = rowOrdering(scip, columnindices, ncons);
   SCIPvectorRearrange(rowindices, roworder);

   columnorder = rowOrdering(scip, rowindices, nvars);
   SCIPvectorRearrange(columnindices, columnorder);

   /* consindex and indexcons */
   for( it1 = roworder.begin(), i = 0; it1 != roworder.end() && i < (size_t) ncons; ++i,++it1 )
   {
      SCIP_CONS* cons;
      position = *it1;
      hashmapindex = &detectordata->hashmapindices[position];

      cons = (SCIP_CONS*) SCIPhashmapGetImage(inputmap->indexcons, (void*) hashmapindex);
      assert ( cons != NULL);

      /* consindex */
      hashmapindex = &detectordata->hashmapindices[i+1];
      assert( SCIPhashmapExists(outputmap->consindex, (void*) cons));
      assert(*hashmapindex <= ncons);
      SCIP_CALL( SCIPhashmapSetImage(outputmap->consindex, (void*) cons, (void*) hashmapindex) );

      /* indexcons */
      assert( SCIPhashmapExists(outputmap->indexcons, (void*) hashmapindex ));
      SCIP_CALL( SCIPhashmapSetImage(outputmap->indexcons, (void*) hashmapindex, cons) );
   }

   /* varindex and indexvar */
   for( it1 = columnorder.begin(), i = 0; it1 != columnorder.end() &&i < (size_t) nvars; ++i, ++it1 )
   {
      SCIP_VAR* var;
      position = *it1;
      hashmapindex = &detectordata->hashmapindices[position];
      var = (SCIP_VAR*) SCIPhashmapGetImage(inputmap->indexvar, (void*) hashmapindex);
      assert ( var != NULL);
      assert(*hashmapindex <= nvars);
      /* varindex */
      hashmapindex = &detectordata->hashmapindices[i+1];
      assert( SCIPhashmapExists(outputmap->varindex, (void*) var) );
      SCIP_CALL( SCIPhashmapSetImage(outputmap->varindex, (void*) var, (void*) hashmapindex) );

      /* indexvar */
      assert( SCIPhashmapExists(outputmap->indexvar, (void*) hashmapindex ));
      SCIP_CALL( SCIPhashmapSetImage(outputmap->indexvar, (void*) hashmapindex, var) );
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE rankOrderClustering(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data structure */
   int                   max_iterations,     /**< number of maximal iterations */
   int*                  iterations          /**< number of performed iterations */
   )
{
   int i;
   INDEXMAP* indexmap_permuted;
   vector<vector<int> > rowindices;
   vector<vector<int> > columnindices;
   int* ibegin_permuted = NULL;
   int* iend_permuted = NULL;
   int* jbegin_permuted = NULL;
   int* jend_permuted = NULL;
   assert(scip != NULL);
   assert(detectordata != NULL);

   int nvars = SCIPgetNVars(scip);
   int nconss = SCIPgetNConss(scip);

   if( iterations != NULL )
      *iterations = -1;

   if( max_iterations <= 0 )
   {
      return SCIP_OKAY;
   }



   SCIP_CALL( indexmapCreate(scip, &indexmap_permuted, nconss, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &ibegin_permuted, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &iend_permuted, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &jbegin_permuted, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &jend_permuted, nvars) );
   SCIP_CALL( indexmapInit(indexmap_permuted, SCIPgetVars(scip), nvars, SCIPgetConss(scip), nconss, detectordata->hashmapindices) );
   i = 0;

   do
   {
      ++i;
      /* not more than max_iterations loops. no iteration limit for max_iterations == -1 */
      if(i > max_iterations && max_iterations != -1)
         break;

      SCIPdebugMessage("Iteration # %i of ROC2\n", i);
      SCIP_CALL( rankOrderClusteringIteration(scip, detectordata, detectordata->indexmap, indexmap_permuted) );

      /* form the new index arrays after the permutation */
      SCIP_CALL( createRowindexList(scip, detectordata, indexmap_permuted->indexcons, indexmap_permuted->varindex, rowindices) );
      SCIP_CALL( createColumnindexList(scip, rowindices, columnindices) );
      SCIP_CALL( formIndexArray(ibegin_permuted, iend_permuted, rowindices) );
      SCIP_CALL( formIndexArray(jbegin_permuted, jend_permuted, columnindices) );
      rowindices.clear();
      columnindices.clear();

      /* switch between index arrays containing new and old indices */
      swap( detectordata->ibegin, ibegin_permuted);
      swap( detectordata->iend, iend_permuted);
      swap( detectordata->jbegin, jbegin_permuted);
      swap( detectordata->jend, jend_permuted);

      /* switch between hash maps containing new and old indices */
      swap(detectordata->indexmap, indexmap_permuted);
   }
   /* while Index Arrays change */
   while( ! (arraysAreEqual(detectordata->ibegin, ibegin_permuted, nconss )
         && arraysAreEqual(detectordata->iend, iend_permuted, nconss)
         && arraysAreEqual(detectordata->jbegin, jbegin_permuted, nvars)
         && arraysAreEqual(detectordata->jend, jend_permuted, nvars)));

   indexmapFree(scip, &indexmap_permuted);
   SCIPfreeMemoryArray(scip, &ibegin_permuted);
   SCIPfreeMemoryArray(scip, &iend_permuted);
   SCIPfreeMemoryArray(scip, &jbegin_permuted);
   SCIPfreeMemoryArray(scip, &jend_permuted);


   if( iterations != NULL )
      *iterations = i-1;

   return SCIP_OKAY;
}

/** finds rows with local minima regarding the number of linking variables and stores them in detectordata->rowsWithConstrictions */
static
SCIP_RETCODE rowsWithConstriction(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< detector data structure */
   )
{
   /* if blocking is performed after row i+1; local minima */
   size_t i;
   size_t nconss = (size_t) SCIPgetNConss(scip);
   for( i = 1; i < nconss - 2; ++i )
   {
      /* is minV[i] a local minimum?    < or <=   ? What does make more sense? */
      if( detectordata->minV[i] < detectordata->minV[i-1] && detectordata->minV[i] < detectordata->minV[i+1] )
      {
         detectordata->rowsWithConstrictions->push_back((int)(i+1));
      }
   }
   return SCIP_OKAY;
}

/** assigns variables to a block, divided into linking variables and nonlinking variables.*/


/** assigns constraints in the interval [first_cons, last_cons] to 'block'. */
static
SCIP_RETCODE assignConsToBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data structure */
   int                   block,              /**< id of block */
   int                   first_cons,         /**< index of first constraint to assign */
   int                   last_cons           /**< index of last constraint to assign */
   )
{
   int i;

   /* assign the constraints to the current block */
   for( i = first_cons; i <= last_cons; ++i )
   {
      int* hashmapindex = &detectordata->hashmapindices[i];
      SCIP_CONS* cons = (SCIP_CONS*) SCIPhashmapGetImage(detectordata->indexmap->indexcons, (void*) hashmapindex);
      assert(cons != NULL);
      /* insert cons into hash map vartoblock */
      assert(!SCIPhashmapExists(detectordata->constoblock, cons));
      SCIP_CALL( SCIPhashmapInsert(detectordata->constoblock, cons, (void*) (size_t) (detectordata->hashmapindices[block])) );
   }
   detectordata->blockedAfterrow->push_back(detectordata->hashmapindices[last_cons]);
   return SCIP_OKAY;
}

/** returns the largest column index of a nonzero entry between rows [from_row, to_row] */
static
int getMaxColIndex(
   DEC_DETECTORDATA*     detectordata,       /**< detector data structure */
   int                   from_row,           /**< index of starting row to check */
   int                   to_row              /**< index of last row to check */
   )
{
   /* some pointer arithmetic */
   return std::max_element(detectordata->iend + (from_row), detectordata->iend+((size_t)to_row + 1))-detectordata->iend; /*lint !e712*/
}

/** returns the column index of the first nonzero entry in 'row'. Rows start counting at 1, not 0. */
static
int getMinColIndex(
   DEC_DETECTORDATA*     detectordata,       /**< detector data structure */
   int                   row                 /**< index of row */
   )
{
   return detectordata->ibegin[row-1];
}

/** determines if a blocking at 'block_at_row' is a valid blocking
 *
 * @return TRUE if blocking is valid, else FALSE
 */
static
SCIP_Bool isValidBlocking(
   DEC_DETECTORDATA*     detectordata,       /**< detector data structure */
   int                   prev_block_first_row, /**< first row of the previous block */
   int                   prev_block_last_row, /**< last row of the previous block */
   int                   block_at_row        /**< the row for which you want to determine if the blocking is valid */
   )
{
   int last_column_prev_block;
   int first_column_current_block;

   /* if the function is called for the first block, the blocking is always valid */
   if( prev_block_last_row == 0 )
   {
      return TRUE;
   }
   last_column_prev_block = getMaxColIndex(detectordata, prev_block_first_row, prev_block_last_row);
   first_column_current_block = getMinColIndex(detectordata, block_at_row);
   return ( first_column_current_block > last_column_prev_block ? TRUE : FALSE);
}

/** this functions looks for rows to block at, which creates block of size min_block_size or bigger
 *
 * @return Iterator pointing to a node which contains a suitable row for blocking; If the iterator points after the last element, no candidate was found
 */
static
vector<int>::iterator findBlockingCandidate(
   vector<int>::iterator it_constrictions,   /**< Iterator pointing to a vector of constraints (detectordata->rowsWithConstrictions) */
   vector<int>*          it_vector,          /**< minimum number of rows to be in a block */
   int                   min_block_size,     /**< the last row of the preceding block */
   int                   prev_block_last_row /**< the last row of the preceding block */
   )
{
   for( ;; )
   {
      /* end of the vector? */
      if( it_constrictions == it_vector->end() )
      {
         return it_constrictions;
      }
      /* does a blocking to the next row forming a constriction comprise more rows than min_block_size */
      if( (*it_constrictions - prev_block_last_row) >= min_block_size )
      {
         return it_constrictions;
      }
      /* advance iterator to next element */
      ++it_constrictions;
   }
}

/** this functions determines the next row to block at
 *
 * @return Iterator pointing to a node which contains a suitable row for blocking; If the iterator points after the last element, no row was found
 */
static
vector<int>::iterator nextRowToBlockAt(
   DEC_DETECTORDATA*     detectordata,       /**< detector data structure */
   vector<int>::iterator it_constrictions,   /**< Iterator pointing to a vector of constraints (detectordata->rowsWithConstrictions) */
   vector<int>*          it_vector,          /**< vector of constrictions */
   int                   min_block_size,     /**< minimum number of rows to be in a block */
   int                   prev_block_first_row, /**< the first row of the preceding block */
   int                   prev_block_last_row /**< the last row of the preceding block*/
   )
{

   /* end of the constriction vector? */
   if( it_constrictions == it_vector->end() )
   {
      return it_constrictions;
   }

   for( ;; )
   {
      /* find a blocking candidate */
      it_constrictions = findBlockingCandidate(it_constrictions, it_vector, min_block_size, prev_block_last_row);
      /* case: no candidate found */
      if( it_constrictions == it_vector->end() )
      {
         break;
      }
      /* case: candidate found */
      else
      {
         /* valid blocking */
         if( isValidBlocking(detectordata, prev_block_first_row, prev_block_last_row, *it_constrictions) )
         {
            break;
         }
         /* invalid blocking */
         else
         {
            ++it_constrictions;
         }
      }
   }
   return it_constrictions;
}

/** calculate the number of decompositions in order to allocate decomps array */
static
int calculateNdecompositions(
   DEC_DETECTORDATA*     detectordata        /**< detector data structure */
   )
{
   int nblockingtypes;
   int nblockingspertype;

   nblockingtypes = 0;
   /* get the number of enabled blocking types */
   if( detectordata->dynamicblocking )
   {
      ++nblockingtypes;
   }
   if( detectordata->staticblocking )
   {
      ++nblockingtypes;
   }
   if( detectordata->blockingassoonaspossible )
   {
      ++nblockingtypes;
   }

   /* get the number of blockings per blocking type */
   if( detectordata->multipledecomps )
   {
      nblockingspertype = detectordata->maxblocks - detectordata->minblocks + 1;
   }
   else
   {
      nblockingspertype = 1;
   }

   return nblockingtypes * nblockingspertype;
}

/** check the consistency of the parameters */
static
void checkParameterConsistency(
   DEC_DETECTORDATA*     detectordata,       /**< detector data structure */
   SCIP_RESULT*          result              /**< pointer to store result */
   )
{
   /* maxblocks < nRelevantsCons? */

   /* desired blocks <= maxblocks? */

   /* is  minblocks <= maxblocks? */
   if( detectordata->multipledecomps )
   {
      if( detectordata->minblocks > detectordata->maxblocks )
      {
         SCIPerrorMessage("minblocks > maxblocks. Setting minblocks = maxblocks.\n");
         detectordata->minblocks = detectordata->maxblocks;
      }
   }

   /* is at least one blocking type enabled? */
   if( ! detectordata->blockingassoonaspossible && ! detectordata->staticblocking &&! detectordata->dynamicblocking )
   {
      SCIPerrorMessage("No blocking type enabled, cannot perform blocking.\n");
      *result = SCIP_DIDNOTRUN;
   }
}

/** tries to dynamically divide the problem into subproblems (blocks)*/
static
SCIP_RETCODE blockingDynamic(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data data structure */
   int                   tau,                /**< desired number of blocks */
   int                   nvars               /**< number of variables in the problem*/
   )
{ /*lint -e715*/
   int block;
   int prev_block_first_row;
   int prev_block_last_row;
   int min_block_size;
   /* notation: i=current block; im1=i-1=previous block; ip1=i+1=next block */
#ifdef SCIP_DEBUG
   int max_col_index_im1 = 0;
#endif
   vector<int>::iterator it1;
   /* debug */
   SCIPdebugMessage("Starting Blocking...\n");
   SCIPdebugMessage("Max blocks: %i\n", detectordata->maxblocks);
   block = 1;
   prev_block_first_row = 0;
   prev_block_last_row = 0;

   assert(tau > 0);
   min_block_size = (int) SCIPround(scip, ((SCIP_Real)SCIPgetNConss(scip)) / (2.0 * tau ));
   it1 = detectordata->rowsWithConstrictions->begin();

   for( it1 = nextRowToBlockAt(detectordata, it1, detectordata->rowsWithConstrictions, min_block_size, prev_block_first_row, prev_block_last_row);
         it1 != detectordata->rowsWithConstrictions->end() && block < detectordata->maxblocks;
         it1 = nextRowToBlockAt(detectordata, it1, detectordata->rowsWithConstrictions, min_block_size, prev_block_first_row, prev_block_last_row) )
   {
      int current_row = * it1;
#ifdef SCIP_DEBUG
      int max_col_index_i = getMaxColIndex(detectordata, prev_block_last_row + 1, current_row);
      int min_col_index_ip1 = getMinColIndex(detectordata, current_row + 1);
      SCIPdebugMessage("vars in block: %i - %i, linking vars: %i - %i\n", max_col_index_im1+1, max_col_index_i, min_col_index_ip1, max_col_index_i);
#endif
      /* assign the variables and constraints to block */
      SCIP_CALL( assignConsToBlock(scip, detectordata, block, prev_block_last_row + 1, current_row) );
      /* update variables in the while loop */
      prev_block_first_row = prev_block_last_row + 1;
      prev_block_last_row = current_row;
      ++block;

#ifdef SCIP_DEBUG
      max_col_index_im1 = max_col_index_i;
#endif
   }
   /* assign the remaining (< M/2tau) cons and vars to the last block; no new linking vars are added */
#ifdef SCIP_DEBUG
   SCIPdebugMessage("last time: vars in block: %i - %i, linking vars: %i - %i\n", max_col_index_im1+1, nvars, nvars+1, nvars);
#endif
   SCIP_CALL( assignConsToBlock(scip, detectordata, block, prev_block_last_row + 1, SCIPgetNConss(scip)) );
   detectordata->blockedAfterrow->pop_back();

   detectordata->blocks = block;

   /* debug plot the blocking  plot for [i=1:2:1] 'test.dat' every :::i::i lt i pt 5 */
#ifdef WRITEALLOUTPUT
   {
      char filename[256];
      char paramfile[256];

      sprintf(filename, "%s_dynamic_minV", getProbNameWithoutPath(scip));
      sprintf(paramfile, "%s_dynamic.params", getProbNameWithoutPath(scip));
      plotMinV(scip, detectordata, filename);
   }
#endif

   return SCIP_OKAY;
}

/** creates blocks with the same number of rows*/
static
SCIP_RETCODE blockingStatic(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< detector data data structure */
   )
{
   int nblocks;
   int block;
   int prev_block_last_row;
   int current_row;
   int nconss;
   /* notation: i=current block; im1=i-1=previous block; ip1=i+1=next block */

   assert(scip != NULL);
   assert(detectordata != NULL);
   nconss = SCIPgetNConss(scip);
   nblocks = nconss/detectordata->nconssperblock;
   prev_block_last_row = 0;
   current_row = 0;

   /* blocks 1 to (desired_blocks-1) */
   for( block = 1; block <= nblocks; ++block )
   {
      current_row += detectordata->nconssperblock;
      SCIPdebugMessage("Blocking from %d to %d in block %d",prev_block_last_row + 1,  current_row, block);
      SCIP_CALL( assignConsToBlock(scip, detectordata, block, prev_block_last_row + 1, current_row) );


      prev_block_last_row = current_row;
   }
   /* last block */
   /* assign the remaining cons and vars to the last block; no new linking vars are added */
   if( prev_block_last_row + 1 <= nconss)
   {
      SCIPdebugMessage("last time: assignVarsToBlock: block, from_row, to_row: %i, %i, %i\n", block, prev_block_last_row + 1, SCIPgetNConss(scip));

      SCIP_CALL( assignConsToBlock(scip, detectordata, block, prev_block_last_row + 1, nconss) );
      detectordata->blockedAfterrow->pop_back();
      ++block;
   }
   detectordata->blocks = block-1;

#ifdef WRITEALLOUTPUT
   {
      char filename[256];
      char paramfile[256];

      /* debug */
      sprintf(filename, "%s_static_minV_%i", getProbNameWithoutPath(scip), detectordata->blocks);
      sprintf(paramfile, "%s_static.params", getProbNameWithoutPath(scip));
      plotMinV(scip, detectordata, filename);
   }
#endif

   return SCIP_OKAY;
}

static
SCIP_RETCODE blockingAsSoonAsPossible(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data data structure */
   int                   desired_blocks,     /**< desired number of blocks */
   int                   nvars               /**< number of variables in the problem*/
)
{
   int block;
   block = 0;
   detectordata->blocks = block;

   return SCIP_OKAY;
}

/** resets detectordata such that it can be used for the next decomposition */
static
SCIP_RETCODE resetDetectordata(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< detector data structure */
   )
{
   if(detectordata->constoblock != NULL)
   {
      SCIP_CALL( SCIPhashmapRemoveAll(detectordata->constoblock) );
   }
   else
   {
      SCIP_CALL( SCIPhashmapCreate(&(detectordata->constoblock), SCIPblkmem(scip), SCIPgetNConss(scip)) );
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE blocking(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data structure */
   DEC_DECOMP***         decdecomps,         /**< array of decompositions */
   int*                  ndecdecomps,        /**< size of decompositions array */
   int                   nvars,              /**< number of variables */
   int                   ncons,              /**< number of constraints */
   SCIP_RESULT*          result              /**< pointer to store result */
   )
{
   int tau = 1; /*  desired number of blocks */

   assert(*ndecdecomps == 0);

   SCIPdebugMessage("Entering Blocking\n");

   /* if multiple decompositions disabled */
   if( detectordata->multipledecomps == FALSE )
   {
      /* if desiredblocks == 0 let the algorithm determine the desired number of blocks */
      if( detectordata->desiredblocks == 0 )
      {
         int n = *std::max_element(detectordata->width, detectordata->width+ncons);
         int v = *std::min_element(detectordata->width, detectordata->width+ncons);
         tau = (int)SCIPround(scip, ((SCIP_Real)nvars - v)/((SCIP_Real)n - v));
         SCIPdebugMessage("<n><v><tau>: <%i><%i><%i>\n", n, v, tau);
         if( tau > detectordata->maxblocks )
         {
            tau = detectordata->maxblocks;
         }

         SCIPdebugMessage("detectordata->enablemultipledecomps = 0. detectordata->desiredblocks == 0. Calculating tau = %i\n", tau);
         /* continue only if tau >= 2 */
         if( tau < 2 )
         {
            *result = SCIP_DIDNOTFIND;
            return SCIP_OKAY;
         }
      }
      else
      {
         tau = detectordata->desiredblocks;
      }
   }

   /* variant 2 */
   /* dynamic blocking */
   if( detectordata->dynamicblocking )
   {
      SCIPdebugMessage("detectordata->enableblockingdynamic = 1.\n");
      SCIP_CALL( rowsWithConstriction(scip, detectordata) );

      SCIPdebugMessage("detectordata->enablemultipledecomps = %ud.\n", detectordata->multipledecomps);

      if( detectordata->multipledecomps )
      {
         for( tau = detectordata->minblocks; tau <= detectordata->maxblocks; ++tau )
         {
            SCIPdebugMessage("tau = %i, dec = %i\n", tau, *ndecdecomps);
            SCIP_CALL( resetDetectordata(scip, detectordata) );

            SCIP_CALL( blockingDynamic(scip, detectordata, tau, nvars) );
            if( detectordata->blocks <= 1 )
               continue;

            SCIP_CALL( DECdecompCreate(scip, &((*decdecomps)[*ndecdecomps])) );
            SCIP_CALL( DECfilloutDecompFromConstoblock(scip, (*decdecomps)[*ndecdecomps], detectordata->constoblock, detectordata->blocks, TRUE) );
            detectordata->constoblock = NULL;

            (*ndecdecomps) += 1;
         }
      }
      else
      {
         SCIPdebugMessage("tau = %i, dec = %i\n", tau, *ndecdecomps);
         SCIP_CALL( resetDetectordata(scip, detectordata) );

         SCIP_CALL( blockingDynamic(scip, detectordata, tau, nvars) );
         if( detectordata->blocks > 1 )
         {
            SCIP_CALL( DECdecompCreate(scip, &((*decdecomps)[*ndecdecomps])) );
            SCIP_CALL( DECfilloutDecompFromConstoblock(scip, (*decdecomps)[*ndecdecomps], detectordata->constoblock, detectordata->blocks, TRUE) );
            detectordata->constoblock = NULL;

            (*ndecdecomps) += 1;
         }
      }
   }

   /* static blocking */
   SCIPdebugMessage("detectordata->staticblocking = %ud. \n", detectordata->staticblocking);

   if( detectordata->staticblocking )
   {
      SCIPdebugMessage("nconssperblock = %i, dec = %i\n", detectordata->nconssperblock, *ndecdecomps);

      SCIP_CALL( resetDetectordata(scip, detectordata) );
      SCIP_CALL( blockingStatic(scip, detectordata) );

      if( detectordata->blocks > 1 )
      {
         SCIP_CALL( DECdecompCreate(scip, &((*decdecomps)[*ndecdecomps])) );
         SCIP_CALL( DECfilloutDecompFromConstoblock(scip, (*decdecomps)[*ndecdecomps], detectordata->constoblock, detectordata->blocks, TRUE) );
         detectordata->constoblock = NULL;

         (*ndecdecomps) += 1;
      }
   }

   /* blocking ASAP */
   SCIPdebugMessage("detectordata->blockingassoonaspossible = %ud. \n", detectordata->blockingassoonaspossible);

   if( detectordata->blockingassoonaspossible )
   {
      if( detectordata->multipledecomps )
      {
         for( tau = detectordata->minblocks; tau <= detectordata->maxblocks; ++tau )
         {
            SCIPdebugMessage("tau = %i, dec = %i\n", tau, *ndecdecomps);
            SCIP_CALL( resetDetectordata(scip, detectordata) );

            SCIP_CALL( blockingAsSoonAsPossible(scip, detectordata, tau, nvars) );
            if( detectordata->blocks <= 1)
               continue;

            SCIP_CALL( DECdecompCreate(scip, &((*decdecomps)[*ndecdecomps])) );
            SCIP_CALL( DECfilloutDecompFromConstoblock(scip, (*decdecomps)[*ndecdecomps], detectordata->constoblock, detectordata->blocks, TRUE) );
            detectordata->constoblock = NULL;

            *ndecdecomps += 1;
         }
      }
      else
      {
         SCIPdebugMessage("tau = %i, dec = %i\n", tau, *ndecdecomps);
         SCIP_CALL( resetDetectordata(scip, detectordata) );
         SCIP_CALL( blockingAsSoonAsPossible(scip, detectordata, tau, nvars) );
         if( detectordata->blocks > 1)
         {
            SCIP_CALL( DECdecompCreate(scip, &((*decdecomps)[*ndecdecomps])) );
            SCIP_CALL( DECfilloutDecompFromConstoblock(scip, (*decdecomps)[*ndecdecomps], detectordata->constoblock, detectordata->blocks, TRUE) );
            detectordata->constoblock = NULL;

            *ndecdecomps += 1;
         }
      }
   }
   return SCIP_OKAY;
}

/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(detectorFreeStairheur)
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
DEC_DECL_INITDETECTOR(detectorInitStairheur)
{
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   detectordata->constoblock = NULL;
   detectordata->ibegin = NULL;
   detectordata->iend = NULL;
   detectordata->jbegin = NULL;
   detectordata->jend = NULL;
   detectordata->jmin = NULL;
   detectordata->jmax = NULL;
   detectordata->minV = NULL;
   detectordata->width = NULL;
   detectordata->hashmapindices = NULL;
   detectordata->indexmap = NULL;
   detectordata->rowsWithConstrictions = NULL;
   detectordata->blockedAfterrow = NULL;

   return SCIP_OKAY;
}

/** detection method of stairheur detector */
static
DEC_DECL_DETECTSTRUCTURE(detectorDetectStairheur)
{
   int i;
   int nconss; /* number of constraints in the problem */
   int nvars; /* number of variables in the problem */
   int ndecs;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   vector<vector<int> > rowindices;
   vector<vector<int> > columnindices;
#ifdef WRITEALLOUTPUT
   int ROC_iterations;
#endif

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(decdecomps != NULL);
   assert(ndecdecomps != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting stairheur structure:");

   SCIP_CALL( initData(scip, detectordata) );

   checkParameterConsistency(detectordata, result);
   ndecs = calculateNdecompositions(detectordata);
   SCIPdebugMessage("%i decompositions will be created\n", ndecs);
   *ndecdecomps = 0;

   /* allocate space for output data */
   SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, ndecs) );

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);
   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);

   /* initialize hash maps for keeping track of variables and constraints and their corresponding indices after being permuted by the ROC2-algorithm */
   SCIP_CALL( indexmapInit(detectordata->indexmap, vars, nvars, conss, nconss, detectordata->hashmapindices) );

#ifdef WRITEALLOUTPUT
   {
      char filename[256];
      sprintf(filename, "%s_initial_problem", getProbNameWithoutPath(scip));
      //plotInitialProblem(scip, detectordata, filename);
   }
#endif

   /* initialize index arrays ibegin, iend, jbegin, jend */
   SCIP_CALL( createRowindexList(scip, detectordata, detectordata->indexmap->indexcons, detectordata->indexmap->varindex, rowindices) );

   SCIP_CALL( createColumnindexList(scip, rowindices, columnindices) );

   SCIP_CALL( formIndexArray(detectordata->ibegin, detectordata->iend, rowindices) );
   SCIP_CALL( formIndexArray(detectordata->jbegin, detectordata->jend, columnindices) );

   /* ==================== */
   /* ===ROC2 algorithm=== */
   /* ==================== */
   SCIPdebugMessage("starting ROC2 algorithm\n");


#ifdef WRITEALLOUTPUT
   SCIP_CALL( rankOrderClustering(scip, detectordata, detectordata->maxiterationsROC, &ROC_iterations) );

   /* check conditions for arrays ibegin and jbegin: ibegin[i]<=ibegin[i+k] for all positive k */
   if( ROC_iterations < detectordata->maxiterationsROC || detectordata->maxiterationsROC  == -1 )
   {
      checkConsistencyOfIndexarrays(detectordata, nvars, nconss);
   }
#else
   SCIP_CALL( rankOrderClustering(scip, detectordata, detectordata->maxiterationsROC, NULL) );
#endif
   /* arrays jmin, jmax and minV */
   SCIPdebugMessage("calculating index arrays\n");
   detectordata->jmin[0] = detectordata->ibegin[0];
   detectordata->jmax[0] = detectordata->iend[0];
   detectordata->width[0] = detectordata->iend[0] - detectordata->ibegin[0];
   for( i = 1; i < nconss; ++i )
   {
      detectordata->width[i] = detectordata->iend[i] - detectordata->ibegin[i];
      detectordata->jmin[i] = detectordata->ibegin[i];
      detectordata->jmax[i] = MAX(detectordata->iend[i], detectordata->jmax[i-1]);
      detectordata->minV[i-1]=1 + (detectordata->jmax[i-1] - detectordata->jmin[i]);
   }
   /* ==================== */
   /* =====BLOCKING======= */
   /* ==================== */

   SCIP_CALL( blocking(scip, detectordata, decdecomps, ndecdecomps, nvars, nconss, result) );
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " found %d decompositions.\n", *ndecdecomps);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " \tBlocks:", *ndecdecomps);
#ifdef WRITEALLOUTPUT
   {
      char filename[256];
      sprintf(filename, "%s_ROC", getProbNameWithoutPath(scip));
      plotInitialProblem(scip, detectordata, filename);
   }
#endif
   for( i = 0; i < *ndecdecomps; ++i )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " %i", DECdecompGetNBlocks( (*decdecomps)[i] ));
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");

   SCIP_CALL( SCIPreallocMemoryArray(scip, decdecomps, *ndecdecomps) );

   SCIP_CALL( freeData(scip, detectordata) );

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** creates the stairheur detector and includes it in SCIP */
extern "C"
SCIP_RETCODE SCIPincludeDetectorStairheur(
   SCIP*                 scip              /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA *detectordata = NULL;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP,
      detectordata, detectorDetectStairheur, detectorFreeStairheur, detectorInitStairheur, NULL) );

   /* add stairheur detector parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/stairheur/nconssperblock",
      "The number of constraints per block (static blocking only)",
      &detectordata->nconssperblock, FALSE, DEFAULT_NCONSSPERBLOCK, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/stairheur/maxblocks",
      "The maximal number of blocks",
      &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/stairheur/minblocks", "The minimal number of blocks",
      &detectordata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/stairheur/desiredblocks",
      "The desired number of blocks. 0 means automatic determination of the number of blocks.",
      &detectordata->desiredblocks, FALSE, DEFAULT_DESIREDBLOCKS, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/stairheur/dynamicblocking",
      "Enable blocking type 'dynamic'",
      &detectordata->dynamicblocking, FALSE, DEFAULT_DYNAMICBLOCKING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/stairheur/staticblocking",
      "Enable blocking type 'static'",
      &detectordata->staticblocking, FALSE, DEFAULT_STATICBLOCKING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/stairheur/blockingassoonaspossible",
      "Enable blocking type 'as soon as possible", &detectordata->blockingassoonaspossible,
      FALSE, DEFAULT_BLOCKINGASSOONASPOSSIBLE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/stairheur/multipledecomps",
      "Enables multiple decompositions for all enabled blocking types. Ranging from minblocks to maxblocks",
      &detectordata->multipledecomps, FALSE, DEFAULT_MULTIPLEDECOMPS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/stairheur/maxiterationsROC",
      "The maximum number of iterations of the ROC-algorithm. -1 for no limit",
      &detectordata->maxiterationsROC, FALSE, DEFAULT_MAXITERATIONSROC, -1, 1000000, NULL, NULL) );

   return SCIP_OKAY;
}
