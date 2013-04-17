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

/**@file   hyperrowgraph.cpp
 * @brief  Column hypergraph
 * @author Martin Bergner
 * @author Annika Thome
 *
 * Hypergraph with a node for every variable and a hyperedge for every constraint
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include "hyperrowgraph.h"
#include "scip_misc.h"
#include <set>
#include <fstream>
#include <algorithm>

using std::ifstream;

namespace gcg
{
template <class T>
HyperrowGraph<T>::HyperrowGraph(
   SCIP*                 scip,              /**< SCIP data structure */
   Weights               w                  /**< weights for the given graph */
): BipartiteGraph<T>(scip, w)
{
   this->name = std::string("hyperrow");
}

template <class T>
HyperrowGraph<T>::~HyperrowGraph()
{
   // TODO Auto-generated destructor stub
}


/** writes the graph to the given file.
 *  The format is graph dependent
 */
template <class T>
SCIP_RETCODE HyperrowGraph<T>::writeToFile(
   const char*        filename,           /**< filename where the graph should be written to */
   SCIP_Bool          edgeweights = FALSE /**< whether to write edgeweights */
 )
{
   FILE* file;
   assert(filename != NULL);
   file = fopen(filename, "w");
   if( file == NULL )
      return SCIP_FILECREATEERROR;

   SCIPinfoMessage(this->scip_, file, "%d %d %d\n", getNEdges(), getNNodes()+this->dummynodes, edgeweights ? 1 :0);

   for( int i = 0; i < getNEdges(); ++i )
   {
      std::vector<int> neighbors = getHyperedgeNodes(i);
      if( edgeweights )
      {
         SCIPinfoMessage(this->scip_, file, "%d ", Graph<T>::getWeight(i+this->nvars));
      }
      for( size_t j = 0; j < neighbors.size(); ++j )
      {
         SCIPinfoMessage(this->scip_, file, "%d ",neighbors[j]+1);
      }
      SCIPinfoMessage(this->scip_, file, "\n");
   }

   if( !fclose(file) )
      return SCIP_OKAY;
   else
      return SCIP_WRITEERROR;
}

template <class T>
int HyperrowGraph<T>::getNEdges()
{
   return this->nconss;
}

template <class T>
int HyperrowGraph<T>::getNNodes()
{
   return this->nvars;
}

template <class T>
int HyperrowGraph<T>::getNNeighbors(
   int i
)
{
   assert(i >= 0);
   assert(i < getNNodes());

   return Graph<T>::getNNeighbors(i);
}

template <class T>
std::vector<int> HyperrowGraph<T>::getNeighbors(
   int i
)
{
   assert(i >= 0);
   assert(i < getNNodes());

   std::vector<int>::iterator it;
   std::set<int> neighbors;

   std::vector<int> immediateneighbors = Graph<T>::getNeighbors(i);
   for( size_t j = 0; j < immediateneighbors.size(); ++j)
   {
      std::vector<int> alternateneighbor = Graph<T>::getNeighbors(immediateneighbors[j]);
      neighbors.insert(alternateneighbor.begin(), alternateneighbor.end() );
   }

   std::vector<int> r(neighbors.begin(), neighbors.end());
   it = std::remove(r.begin(), r.end(), i);

   return std::vector<int>(r.begin(), it);
}

template <class T>
std::vector<int> HyperrowGraph<T>::getHyperedgeNodes(
   int i
)
{
   assert(i >= 0);
   assert(i < getNEdges());

   std::vector<int> neighbors = Graph<T>::getNeighbors(i+this->nvars);
   return neighbors;
}

template <class T>
SCIP_RETCODE HyperrowGraph<T>::createDecompFromPartition(
   DEC_DECOMP**       decomp              /**< decomposition structure to generate */
)
{
   int nblocks;
      SCIP_HASHMAP* constoblock;

      int *nsubscipconss;
      int i;
      SCIP_CONS **conss;
      SCIP_VAR **vars;
      SCIP_Bool emptyblocks = FALSE;

      conss = SCIPgetConss(this->scip_);
      vars = SCIPgetVars(this->scip_);
      nblocks = *(std::max_element(this->partition.begin(), this->partition.end()))+1;

      SCIP_CALL( SCIPallocBufferArray(this->scip_, &nsubscipconss, nblocks) );
      BMSclearMemoryArray(nsubscipconss, nblocks);

      SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(this->scip_), this->nconss) );

      /* assign constraints to partition */
      for( i = 0; i < this->nconss; i++ )
      {

         std::set<int> blocks;
         std::vector<int> neighbors = getHyperedgeNodes(i);
         for( size_t k = 0; k < neighbors.size(); ++k )
         {
            if( this->partition[neighbors[k]] >= 0 )
               blocks.insert(this->partition[neighbors[k]]);
         }
         if( blocks.size() > 1 )
         {
            SCIP_CALL( SCIPhashmapInsert(constoblock, conss[i], (void*) (size_t) (nblocks+1)) );
         }
         else
         {
            int block = *(blocks.begin());
            SCIP_CALL( SCIPhashmapInsert(constoblock, conss[i], (void*) (size_t) (block +1)) );
            ++(nsubscipconss[block]);
         }
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
} /* namespace gcg */
