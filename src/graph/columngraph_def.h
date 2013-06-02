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

/**@file   columngraph.cpp
 * @brief  A row graph where each column is a node and columns are adjacent if they appear in one row
 * @author Martin Bergner
 * @author Annika Thome
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_COLUMNGRAPH_DEF_H_
#define GCG_COLUMNGRAPH_DEF_H_

#include "columngraph.h"
#include <algorithm>

namespace gcg {

template <class T>
ColumnGraph<T>::ColumnGraph(
   SCIP*                 scip,              /**< SCIP data structure */
   Weights               w                  /**< weights for the given graph */
   ) : MatrixGraph<T>(scip, w), graph(scip),nconss(0),nvars(0),nnonzeroes(0)
{
   this->graphiface = &graph;
   this->name = std::string("columngraph");
}

template <class T>
ColumnGraph<T>::~ColumnGraph()
{
   // TODO Auto-generated destructor stub
}


/** writes column graph to file */
template <class T>
SCIP_RETCODE ColumnGraph<T>::writeToFile(
   const char*        filename,           /**< filename where the graph should be written to */
   SCIP_Bool          writeweights         /**< whether to write weights */
   )
{
   int nedges;
   int* nrealneighbors;
   int** realneighbors;

   SCIP_Bool* handled;
   FILE* file;
   assert(filename != NULL);
   file = fopen(filename, "w");
   if( file == NULL )
      return SCIP_FILECREATEERROR;

   nrealneighbors = 0;
   nedges = 0;

   SCIP_CALL( SCIPallocMemoryArray(this->scip_, &handled, this->nvars) );
   SCIP_CALL( SCIPallocMemoryArray(this->scip_, &realneighbors, this->nvars) );
   SCIP_CALL( SCIPallocMemoryArray(this->scip_, &nrealneighbors, this->nvars) );

   for( int i = 0; i < this->nvars; ++i )
   {
      BMSclearMemoryArray(handled, this->nvars);
      handled[i] = TRUE;
      nrealneighbors[i] = 0;

      SCIP_CALL( SCIPallocMemoryArray(this->scip_, &realneighbors[i], this->nvars) );
      int nneighbors = graph.getNNeighbors(i);

      SCIPdebugMessage("%d has %d neighbors\n", i, nneighbors);

      std::vector<int> neighbors = graph.getNeighbors(i);
      for( int j = 0; j < nneighbors; ++j )
      {
         int neighbor = neighbors[j];
         int nneighborneighbors = graph.getNNeighbors(neighbor);
         std::vector<int> neighborneighbors = graph.getNeighbors(neighbor);
         SCIPdebugMessage("\tneighbor %d has %d neighbors\n", neighbor, nneighborneighbors);
         for( int k = 0; k < nneighborneighbors; ++k )
         {
            int neighborneighbor = neighborneighbors[k];

            SCIPdebugMessage("\t\t%d->%d->%d (", i, neighbor, neighborneighbor);
            if( !handled[neighborneighbor] )
            {
               SCIPdebugPrintf("x)\n");
               realneighbors[i][nrealneighbors[i]] = neighborneighbor;
               ++(nrealneighbors[i]);

               handled[neighborneighbor] = TRUE;
               ++nedges;
            }
            else
            {
               SCIPdebugPrintf("-)\n");
            }
         }
      }
   }

   SCIPinfoMessage(this->scip_, file, "%d %d\n", this->nvars, nedges);

   for( int i = 0; i < this->nvars; ++i)
   {
      for( int j = 0; j < nrealneighbors[i]; ++j )
      {
         SCIPinfoMessage(this->scip_, file, "%d ", realneighbors[i][j]+1);
      }
      SCIPinfoMessage(this->scip_, file, "\n");
      SCIPfreeMemoryArray(this->scip_, &realneighbors[i]);
   }

   for( int i = 0; i < this->dummynodes; ++i )
   {
      SCIPinfoMessage(this->scip_, file, "\n");
   }

   SCIPfreeMemoryArray(this->scip_, &handled);
   SCIPfreeMemoryArray(this->scip_, &realneighbors);
   SCIPfreeMemoryArray(this->scip_, &nrealneighbors);

   return SCIP_OKAY;
}

template <class T>
SCIP_RETCODE ColumnGraph<T>::createDecompFromPartition(
   DEC_DECOMP**         decomp
)
{
   int nblocks;
   SCIP_HASHMAP* constoblock;

   int *nsubscipconss;
   int i;
   SCIP_CONS **conss;
   SCIP_VAR **vars;
   SCIP_Bool emptyblocks = FALSE;
   std::vector<int> partition = graph.getPartition();
   conss = SCIPgetConss(this->scip_);
   vars = SCIPgetVars(this->scip_);
   nblocks = *(std::max_element(partition.begin(), partition.end()))+1;

   SCIP_CALL( SCIPallocBufferArray(this->scip_, &nsubscipconss, nblocks) );
   BMSclearMemoryArray(nsubscipconss, nblocks);

   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(this->scip_), this->nconss) );

   /* assign constraints to partition */
   for( i = 0; i < this->nconss; i++ )
   {
      int block = partition[i];
      SCIP_CALL( SCIPhashmapInsert(constoblock, conss[i], (void*) (size_t) (block +1)) );
      ++(nsubscipconss[block]);
   }

   /* first, make sure that there are constraints in every block, otherwise the hole thing is useless */
   for( i = 0; i < nblocks; ++i )
   {
      if( nsubscipconss[i] == 0 )
      {
         SCIPdebugMessage("Block %d does not have any constraints!\n", i);
         emptyblocks = TRUE;
      }
   }

   if( !emptyblocks )
   {
      SCIP_CALL( DECdecompCreate(this->scip_, decomp) );
      SCIP_CALL( DECfilloutDecdecompFromConstoblock(this->scip_, *decomp, constoblock, nblocks, vars, this->nvars, conss, this->nconss, FALSE) );
   }
   else {
      SCIPhashmapFree(&constoblock);
      *decomp = NULL;
   }

   SCIPfreeBufferArray(this->scip_, &nsubscipconss);
   return SCIP_OKAY;
}

template <class T>
SCIP_RETCODE ColumnGraph<T>::createFromMatrix(
   SCIP_CONS**           conss,              /**< constraints for which graph should be created */
   SCIP_VAR**            vars,               /**< variables for which graph should be created */
   int                   nconss_,             /**< number of constraints */
   int                   nvars_               /**< number of variables */
   )
{
   this->nvars = nvars_;
   this->nconss = nconss_;
   SCIP_CALL( graph.createFromMatrix(conss, vars, nconss_, nvars_) );
   return SCIP_OKAY;
}

} /* namespace gcg */

#endif
