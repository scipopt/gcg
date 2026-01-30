/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   hypergraph.h
 * @brief  miscellaneous hypergraph methods for structure detection
 * @author Martin Bergner
 * @author Annika Thome
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



#ifndef GCG_HYPERGRAPH_H_
#define GCG_HYPERGRAPH_H_
#include "objscip/objscip.h"
#include "tclique/tclique.h"
#include "graph/weights.h"
#include "gcg/pub_decomp.h"
#include "graph/graph.h"
#include "graph/graph_interface.h"

#include <exception>
#include <vector>
#include <string>

namespace gcg {

template <class T>
class Hypergraph : public GraphInterface {
public:
   std::string name;
protected:
   GCG* gcg;
   Graph<T>* graph;
   std::vector<int> nodes;
   std::vector<int> hedges;
   std::vector<int> mapping;
   int lastnode;
   int dummynodes;

public:
   /** Constructor */
   Hypergraph(
      GCG*                  gcgstruct
   );

   void swap(Hypergraph & other) // the swap member function (should never fail!)
   {
      // swap all the members (and base subobject, if applicable) with other
      std::swap(partition, other.partition);
      std::swap(gcg , other.gcg);
      std::swap(graph , other.graph);
      std::swap(hedges , other.hedges);
      std::swap(nodes , other.nodes);
      std::swap(lastnode, other.lastnode);
      std::swap(dummynodes, other.dummynodes);
   }

   Hypergraph& operator=(Hypergraph other) // note: argument passed by value!
   {
      // swap this with other
      swap(other);

      return *this;
   }

   /** Destruktor */
   ~Hypergraph();

   /** adds the node with the given weight to the graph */
   SCIP_RETCODE addNode(int i,int weight);

   /** adds the edge to the graph */
   SCIP_RETCODE addHyperedge(std::vector<int> &edge, int weight);

   /** adds the edge to the graph */
   SCIP_RETCODE addNodeToHyperedge(int node, int hedge);

   /** return the number of nodes */
   int getNNodes();

   /** return the number of edges (or hyperedges) */
   int getNHyperedges();

   /** return the number of neighbor nodes of given node */
   int getNNeighbors(
      int                i                   /**< the given node */
      );

   /** return the neighboring nodes of a given node */
   std::vector<int> getNeighbors(
      int                i                   /**< the given node */
      );

   /** return the nodes spanned by hyperedge */
   std::vector<int> getHyperedgeNodes(
      int i
      );

   /** return the number of nodes spanned by hyperedge */
   int getNHyperedgeNodes(
      int i
      );

   /** assigns partition to a given node*/
   void setPartition(int i, int ID);

   /** writes the hypergraph to the given file.
    *  The format is hypergraph dependent
    */
   SCIP_RETCODE writeToFile(
      int                fd,                  /**< filename where the graph should be written to */
      SCIP_Bool          writeweights        /**< whether to write weights */
    );

   /**
    * reads the partition from the given file.
    * The format is hypergraph dependent. The default is a file with one line for each node a
    */
   SCIP_RETCODE readPartition(
      const char*        filename            /**< filename where the partition is stored */
   );

   /** return the weight of given node */
   int getWeight(
      int                i                   /**< the given node */
      );

   /** return the weight of given hyperedge */
   int getHyperedgeWeight(
      int                i                   /**< the given hyperedge */
      );

   /** set the number of dummy nodes */
   void setDummynodes(int dummynodes_)
   {
      dummynodes = dummynodes_;
   }


   int getDummynodes() const
   {
      return dummynodes;
   }

   SCIP_RETCODE flush();
private:
   int computeNodeId(int i);
};

}

#endif
