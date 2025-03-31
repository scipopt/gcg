/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
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

/**@file   graph_tclique.h
 * @brief  interface to the SCIP tclique graph library
 * @author Annika Thome
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_GRAPH_TCLIQUE_H_
#define GCG_GRAPH_TCLIQUE_H_

#include "graph/bridge.h"
#include "tclique/tclique.h"

namespace gcg {

class GraphTclique: public gcg::Bridge
{
private:
   TCLIQUE_GRAPH* graph;

public:

   GraphTclique();

   virtual ~GraphTclique();
   virtual SCIP_RETCODE addNNodes(int _n_nodes);
   virtual SCIP_RETCODE addNNodes(int _n_nodes, std::vector<int> weights);
   virtual int getNNodes();
   virtual int getNEdges();
   virtual SCIP_RETCODE getEdges(std::vector<void*>& edges);
   virtual SCIP_Bool isEdge(int i, int j);
   virtual int getNNeighbors(int i);
   virtual std::vector<int> getNeighbors(int i);
   virtual SCIP_RETCODE addNode(int i, int weight);
   virtual SCIP_RETCODE addNode();
   virtual SCIP_RETCODE deleteNode(int i);
   virtual SCIP_RETCODE addEdge(int i, int j);
   virtual SCIP_RETCODE addEdge(int i, int j, double weight);
   virtual SCIP_RETCODE setEdge(int i, int j, double weight);
   virtual double getEdgeWeight(int i, int j);
   virtual std::vector<std::pair<int, double> > getNeighborWeights(int i);
   virtual SCIP_RETCODE deleteEdge(int i, int j);
   virtual int graphGetWeights(int i);

   virtual SCIP_RETCODE flush();
   virtual SCIP_RETCODE normalize();
   virtual double getEdgeWeightPercentile(double q);

#ifdef WITH_GSL
   void expand(int factor) {return;}
   void inflate(double factor) {return;}
   void colL1Norm() {return;}
   void prune() {return;}
   bool stopMCL(int iter) {return true;}
   std::vector<int> getClustersMCL() {return std::vector<int>();}
   virtual void initMCL() {return;}
   virtual void clearMCL() {return;}
#endif
};

} /* namespace gcg */
#endif /* GCG_GRAPH_TCLIQUE_H_ */
