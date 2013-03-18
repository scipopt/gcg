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
   Weights               &w                 /**< weights for the given graph */
): BipartiteGraph(scip, w)
{

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
   file = fopen(filename, "wx");
   if( file == NULL )
      return SCIP_FILECREATEERROR;

   SCIPinfoMessage(scip_, file, "%d %d %d\n", getNEdges(), getNNodes(), edgeweights ? 1 :0);

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

} /* namespace gcg */
