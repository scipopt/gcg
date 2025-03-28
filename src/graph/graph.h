/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   graph.h
 * @brief  miscellaneous graph methods for structure detection
 * @author Martin Bergner
 * @author Annika Thome
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



#ifndef GCG_GRAPH_H_
#define GCG_GRAPH_H_
#include "objscip/objscip.h"
#include "tclique/tclique.h"
#include "graph/weights.h"
#include "gcg/pub_decomp.h"
#include "graph/bridge.h"
#include "graph/graph_interface.h"
#include <exception>
#include <vector>
#include <string>

namespace gcg {

template <class T>
class Graph : public GraphInterface {
public:
   std::string name;
protected:
   GCG* gcg;
   Bridge* graph;
   int nconss;
   int nvars;
   int nnonzeroes;
   int dummynodes;

public:
   /** Constructor */
   Graph(
      GCG*                  gcgstruct                /**< GCG data structure */
   );

   void swap(Graph & other) // the swap member function (should never fail!)
   {
      // swap all the members (and base subobject, if applicable) with other
      std::swap(partition, other.partition);
      std::swap(gcg , other.gcg);
      std::swap(graph , other.graph);
      std::swap(nconss , other.nconss);
      std::swap(nvars , other.nvars);
      std::swap(nnonzeroes , other.nnonzeroes);
      std::swap(dummynodes, other.dummynodes);
   }

   Graph& operator=(Graph other) // note: argument passed by value!
   {
      // swap this with other
      swap(other);

      return *this;
   }

   /** Destruktor */
   virtual ~Graph();

   /** adds n nodes in the graph at the same time. it is much faster than to call addNode() many times */
   SCIP_RETCODE addNNodes(int _n_nodes);

   /** adds n nodes in the graph at the same time. it is much faster than to call addNode() many times. weights represent node weights */
   SCIP_RETCODE addNNodes(int _n_nodes, std::vector<int> weights);

   /** adds the node with the given weight to the graph */
   SCIP_RETCODE addNode(int i,int weight);

   /** adds the node with the 0 weight to the graph */
   SCIP_RETCODE addNode();

   /** adds the edge to the graph */
   SCIP_RETCODE addEdge(int i, int j);

   /** adds the weighted edge to the graph */
   SCIP_RETCODE addEdge(int i, int j, double weight);

   /** sets the weight of the edge in the graph */
   SCIP_RETCODE setEdge(int i, int j, double weight);

   /** returns the weight of the edge in the graph */
   double getEdgeWeight(int i, int j);

   std::vector<std::pair<int, double> > getNeighborWeights(int i);

   /** return the number of nodes */
   int getNNodes();

   /** return the number of edges (or hyperedges) */
   int getNEdges();

   /** get list of edges in the graph (not defined how edges are implemented) */
   SCIP_RETCODE getEdges(std::vector<void*>& edges);

   /** returns whether there is an edge between nodes i and j */
   virtual int edge(int i, int j);

   /** return the number of neighbor nodes of given node */
   virtual int getNNeighbors(
      int                i                   /**< the given node */
      );

   /** return the neighboring nodes of a given node */
   virtual std::vector<int> getNeighbors(
      int                i                   /**< the given node */
      );

   /** assigns partition to a given node*/
   virtual void setPartition(int i, int ID);

   /** create graph from the matrix, to be overriden by the implementation*/
   virtual SCIP_RETCODE createFromMatrix(
      SCIP_CONS**        conss,              /**< constraints for which graph should be created */
      SCIP_VAR**         vars,               /**< variables for which graph should be created */
      int                nconss_,            /**< number of constraints */
      int                nvars_              /**< number of variables */
   ) { return SCIP_ERROR; }

   /** writes the graph to the given file.
    *  The format is graph dependent
    */
   virtual SCIP_RETCODE writeToFile(
      int                fd,                  /**< filename where the graph should be written to */
      SCIP_Bool          writeweights        /**< whether to write weights */
    );


   /**
    * reads the partition from the given file.
    * The format is graph dependent. The default is a file with one line for each node a
    */
   virtual SCIP_RETCODE readPartition(
      const char*        filename            /**< filename where the partition is stored */
   );

   int getNNonzeroes() const
   {
      return nnonzeroes;
   }

   /** return the weight of given node */
   virtual int getWeight(
      int                i                   /**< the given node */
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

   SCIP_RETCODE normalize();

   virtual double getEdgeWeightPercentile(double q);

#ifdef WITH_GSL
   /** function needed for MST clustering */
   virtual void expand(int factor);

   /** function needed for MST clustering */
   virtual void inflate(double factor);

   /** function needed for MST clustering */
   virtual void colL1Norm();

   /** function needed for MST clustering */
   virtual void prune();

   /** function needed for MST clustering */
   virtual bool stopMCL(int iter);

   /** function needed for MST clustering */
   virtual std::vector<int> getClustersMCL();

   /** function needed for MST clustering */
   virtual void initMCL();

   /** function needed for MST clustering */
   virtual void clearMCL();
#endif

};

}

#endif
