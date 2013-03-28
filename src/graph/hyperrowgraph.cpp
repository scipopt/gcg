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
 *
 * Hypergraph with a node for every variable and a hyperedge for every constraint
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "hyperrowgraph.h"

#include <set>
#include <fstream>
#include <algorithm>

using std::ifstream;

namespace gcg
{

HyperrowGraph::HyperrowGraph(
   SCIP*                 scip,              /**< SCIP data structure */
   Weights               w                  /**< weights for the given graph */
): BipartiteGraph(scip, w)
{
   name = std::string("hyperrow");
}

HyperrowGraph::~HyperrowGraph()
{
   // TODO Auto-generated destructor stub
}


/** writes the graph to the given file.
 *  The format is graph dependent
 */
SCIP_RETCODE HyperrowGraph::writeToFile(
   const char*        filename,           /**< filename where the graph should be written to */
   SCIP_Bool          edgeweights = FALSE /**< whether to write edgeweights */
 )
{
   FILE* file;
   assert(filename != NULL);
   file = fopen(filename, "w");
   if( file == NULL )
      return SCIP_FILECREATEERROR;

   SCIPinfoMessage(scip_, file, "%d %d %d\n", getNEdges(), getNNodes()+dummynodes, edgeweights ? 1 :0);

   for( int i = 0; i < getNEdges(); ++i )
   {
      std::vector<int> neighbors = getHyperedgeNodes(i);
      if( edgeweights )
      {
         SCIPinfoMessage(scip_, file, "%d ", Graph::getWeight(i+nvars));
      }
      for( size_t j = 0; j < neighbors.size(); ++j )
      {
         SCIPinfoMessage(scip_, file, "%d ",neighbors[j]+1);
      }
      SCIPinfoMessage(scip_, file, "\n");
   }

   if( !fclose(file) )
      return SCIP_OKAY;
   else
      return SCIP_WRITEERROR;
}

int HyperrowGraph::getNEdges()
{
   return nconss;
}


int HyperrowGraph::getNNodes()
{
   return nvars;
}

std::vector<int> HyperrowGraph::getNeighbors(
   int i
)
{
   assert(i >= 0);
   assert(i < getNNodes());

   std::vector<int>::iterator it;
   std::set<int> neighbors;

   std::vector<int> immediateneighbors = Graph::getNeighbors(i);
   for( size_t j = 0; j < immediateneighbors.size(); ++j)
   {
      std::vector<int> alternateneighbor = Graph::getNeighbors(immediateneighbors[j]);
      neighbors.insert(alternateneighbor.begin(), alternateneighbor.end() );
   }

   std::vector<int> r(neighbors.begin(), neighbors.end());
   it = std::remove(r.begin(), r.end(), i);

   return std::vector<int>(r.begin(), it);
}

std::vector<int> HyperrowGraph::getHyperedgeNodes(
   int i
)
{
   assert(i >= 0);
   assert(i < getNEdges());

   std::vector<int> neighbors = Graph::getNeighbors(i+nvars);
   return neighbors;
}

SCIP_RETCODE HyperrowGraph::createDecompFromPartition(
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

      conss = SCIPgetConss(scip_);
      vars = SCIPgetVars(scip_);
      nblocks = *(std::max_element(partition.begin(), partition.end()))+1;

      SCIP_CALL( SCIPallocBufferArray(scip_, &nsubscipconss, nblocks) );
      BMSclearMemoryArray(nsubscipconss, nblocks);

      SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip_), nconss) );

      /* assign constraints to partition */
      for( i = 0; i < nconss; i++ )
      {
         std::set<int> blocks;
         std::vector<int> neighbors = getHyperedgeNodes(i);
         for( size_t k = 0; k < neighbors.size(); ++k )
         {
            blocks.insert(partition[neighbors[k]]);
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
         SCIP_CALL( DECdecompCreate(scip_, decomp) );
         SCIP_CALL( DECfilloutDecdecompFromConstoblock(scip_, *decomp, constoblock, nblocks, vars, nvars, conss, nconss, FALSE) );
      }
      else {
         SCIPhashmapFree(&constoblock);
         *decomp = NULL;
      }

      SCIPfreeBufferArray(scip_, &nsubscipconss);
      return SCIP_OKAY;
}
} /* namespace gcg */
