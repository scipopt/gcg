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

/**@file   hypercolgraph.cpp
 * @brief  Column hypergraph
 * @author Martin Bergner
 * @author Annika Thome
 *
 * A hypergraph structure with a node for every constraint and a hyperedge for every variable.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_HYPERCOLGRAPH_DEF_H_
#define GCG_HYPERCOLGRAPH_DEF_H_

#include "hypercolgraph.h"
#include <set>
#include <fstream>
#include <algorithm>

using std::ifstream;
namespace gcg
{
template <class T>
HypercolGraph<T>::HypercolGraph(
   SCIP*                 scip,              /**< SCIP data structure */
   Weights               w                  /**< weights for the given graph */
):  BipartiteGraph<T>(scip, w)
{
   this->name = std::string("hypercol");
}

template <class T>
HypercolGraph<T>::~HypercolGraph()
{
   // TODO Auto-generated destructor stub
}


/** writes the graph to the given file.
 *  The format is graph dependent
 */
template <class T>
SCIP_RETCODE HypercolGraph<T>::writeToFile(
   const char*        filename,           /**< filename where the graph should be written to */
   SCIP_Bool          edgeweights = FALSE /**< whether to write edgeweights */
 )
{
   function f(this->nvars);
   FILE* file;
   assert(filename != NULL);
   file = fopen(filename, "w");
   if( file == NULL )
      return SCIP_FILECREATEERROR;

   SCIPinfoMessage(this->scip_, file, "%d %d %d\n", getNEdges(), getNNodes()+this->dummynodes, edgeweights ? 1 :0);

   for( int i = 0; i < getNEdges(); ++i )
   {
      std::vector<int> neighbors = getHyperedgeNodes(i);
      int nneighbors = Graph<T>::getNNeighbors(i);
      if( edgeweights )
      {
         SCIPinfoMessage(this->scip_, file, "%d ", Graph<T>::getWeight(i));
      }
      for( int j = 0; j < nneighbors; ++j )
      {
         SCIPinfoMessage(this->scip_, file, "%d ", neighbors[j]+1);
      }
      SCIPinfoMessage(this->scip_, file, "\n");
   }

   if( !fclose(file) )
      return SCIP_OKAY;
   else
      return SCIP_WRITEERROR;
}

template <class T>
int HypercolGraph<T>::getNEdges()
{
   return this->nvars;
}

template <class T>
int HypercolGraph<T>::getNNodes()
{
   return this->nconss;
}

template <class T>
std::vector<int> HypercolGraph<T>::getNeighbors(
   int i
)
{
   assert(i >= 0);
   assert(i < this->nconss);
   function f(this->nvars);
   std::vector<int>::iterator it;
   std::set<int> neighbors;
   std::vector<int> immediateneighbors = Graph<T>::getNeighbors(i+this->nvars);
   for( size_t j = 0; j < immediateneighbors.size(); ++j)
   {
      std::vector<int> alternateneighbor = Graph<T>::getNeighbors(immediateneighbors[j]);
      neighbors.insert(alternateneighbor.begin(), alternateneighbor.end() );
   }
   std::vector<int> r(neighbors.size(), 0);
   std::transform(neighbors.begin(), neighbors.end(), r.begin(), f);
   it = std::remove(r.begin(), r.end(), i);

   return std::vector<int>(r.begin(), it);
}

template <class T>
std::vector<int> HypercolGraph<T>::getHyperedgeNodes(
   int i
)
{
   function f(this->nvars);
   assert(i >= 0);
   assert(i < this->nvars);

   std::vector<int> neighbors = Graph<T>::getNeighbors(i);
   std::transform(neighbors.begin(), neighbors.end(), neighbors.begin(), f);
   return neighbors;
}

} /* namespace gcg */

#endif
