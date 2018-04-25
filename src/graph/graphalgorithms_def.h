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

/**@file    graphalgorithms_def.h
 * @brief   several metrics for graphs
 * @author  Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_GRAPHALGORITHMS_DEF_H_
#define GCG_GRAPHALGORITHMS_DEF_H_

#include "graphalgorithms.h"
#include <vector>
#include <algorithm>
#include "graph/graph_tclique.h"
using std::vector;

namespace gcg {

/** compute weighted sum of external degrees */
template<class T>
double GraphAlgorithms<T>::computeSoed(
      Hypergraph<T>&     graph               /**< the hypergraph */
)
{
   SCIP_Real soed = 0.0;
   size_t nedges = graph.getNHyperedges();
   vector<int> partition = vector<int>(graph.getPartition());

   for( size_t i = 0; i < nedges; ++i )
   {
      vector<int> nodes = graph.getHyperedgeNodes(i);
      for( auto &it : nodes)
      {
         it = partition[it];
      }
      auto end = std::unique(nodes.begin(), nodes.end());
      if( end - nodes.begin() > 1)
         soed += ( end - nodes.begin())*graph.getHyperedgeWeight(i);
   }
   return soed;
}

/** compute minimum hyperedge cut */
template<class T>
double GraphAlgorithms<T>::computeMincut(
      Hypergraph<T>&     graph               /**< the hypergraph */
)
{
   SCIP_Real mincut = 0.0;
   size_t nedges = graph.getNHyperedges();
   vector<int> partition = vector<int>(graph.getPartition());

   for( size_t i = 0; i < nedges; ++i )
   {
      vector<int> nodes = graph.getHyperedgeNodes(i);
      for( auto &it : nodes)
      {
         it = partition[it];
      }
      auto end = std::unique(nodes.begin(), nodes.end());

      if( end - nodes.begin() > 1)
         mincut += graph.getHyperedgeWeight(i);
   }

   return mincut;
}

/** compute k-1 metric */
template<class T>
double GraphAlgorithms<T>::computekMetric(
      Hypergraph<T>&     graph               /**< the hypergraph */
)
{
   SCIP_Real kmetric = 0.0;
   size_t nedges = graph.getNHyperedges();
   vector<int> partition = vector<int>(graph.getPartition());

   for( size_t i = 0; i < nedges; ++i )
   {
      vector<int> nodes = graph.getHyperedgeNodes(i);
      for( auto &it : nodes)
      {
         it = partition[it];
      }
      auto end = std::unique(nodes.begin(), nodes.end());

      kmetric += ( end - nodes.begin() -1)*graph.getHyperedgeWeight(i);
   }

   return kmetric;
}

} /* namespace gcg */
#endif
