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

/**@file   hypercolgraph.h
 * @brief  Column hypergraph
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_HYPERCOLGRAPH_H_
#define GCG_HYPERCOLGRAPH_H_

#include "bipartitegraph.h"

namespace gcg
{
template <class T>
class HypercolGraph: public gcg::BipartiteGraph<T>
{
protected:
   class function {
      int diff;
   public:
      function(int i):diff(i) {}
      int operator()(int i) { return i-diff;}
   };

public:
   HypercolGraph(
      SCIP*                 scip,              /**< SCIP data structure */
      Weights               w                  /**< weights for the given graph */
   );

   virtual ~HypercolGraph();

   /** writes the graph to the given file.
    *  The format is graph dependent
    */
   SCIP_RETCODE writeToFile(
      const char*        filename,           /**< filename where the graph should be written to */
      SCIP_Bool          edgeweights         /**< whether to write edgeweights */
    );

   /** return the number of nodes */
   virtual int getNNodes();

   /** return the number of edges (or hyperedges) */
   virtual int getNEdges();

   virtual std::vector<int> getNeighbors(
         int i
      );

   virtual std::vector<int> getHyperedgeNodes(
         int i
      );
};

} /* namespace gcg */
#endif /* GCG_HYPERCOLGRAPH_H_ */
