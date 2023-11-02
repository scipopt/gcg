/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       */
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

/**@file    graphalgorithm_test.cpp
 * @brief   description
 * @author  Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "test.h"
#include "scip/scip.h"
#include "graph/hypergraph.h"
#include "graph/graphalgorithms.h"
#include "graph/graphalgorithms_def.h"
#include "graph/graph_tclique.h"
#include "graph/rowgraph_weighted.h"
#include "graph/graph_gcg.h"
#include "graph/graph.h"

using gcg::GraphAlgorithms;
using gcg::GraphTclique;
using gcg::Hypergraph;

using namespace std;

class GraphAlgorithmEmptyTest : public ::testing::Test {

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      graph = new Hypergraph<GraphTclique>(scip);
   }

   virtual void TearDown() {
      SCIPfree(&scip);
      delete graph;
   }

protected:
   gcg::Hypergraph<GraphTclique> *graph;
   SCIP* scip;
};

TEST_F(GraphAlgorithmEmptyTest, EmptySoed) {
   ASSERT_NEAR(0.0, GraphAlgorithms<GraphTclique>::computeSoed(*graph), 1e-6);
}
TEST_F(GraphAlgorithmEmptyTest, EmptyMincut) {
   ASSERT_NEAR(0.0, GraphAlgorithms<GraphTclique>::computeMincut(*graph), 1e-6);
}
TEST_F(GraphAlgorithmEmptyTest, EmptyKmetric) {
   ASSERT_NEAR(0.0, GraphAlgorithms<GraphTclique>::computekMetric(*graph), 1e-6);
}

class GraphAlgorithmSmallTest : public ::testing::Test {

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      graph = new Hypergraph<GraphTclique>(scip);

      graph->addNode(0,1);
      graph->addNode(1,1);
      graph->flush();

      int edge[2] = {0,1};
      std::vector<int> hedge(edge, edge+2);

      graph->addHyperedge(hedge, 1);
      graph->setPartition(0, 1);
      graph->setPartition(1, 1);
      graph->flush();

   }

   virtual void TearDown() {
      SCIPfree(&scip);
      delete graph;
   }

protected:
   gcg::Hypergraph<GraphTclique> *graph;
   SCIP* scip;
};

TEST_F(GraphAlgorithmSmallTest, SmallSoed) {
   ASSERT_NEAR(0.0, GraphAlgorithms<GraphTclique>::computeSoed(*graph), 1e-6);
}
TEST_F(GraphAlgorithmSmallTest, SmallMincut) {
   ASSERT_NEAR(0.0, GraphAlgorithms<GraphTclique>::computeMincut(*graph), 1e-6);
}
TEST_F(GraphAlgorithmSmallTest, SmallKmetric) {
   ASSERT_NEAR(0.0, GraphAlgorithms<GraphTclique>::computekMetric(*graph), 1e-6);
}

class GraphAlgorithmSmallCutTest : public ::testing::Test {

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      graph = new Hypergraph<GraphTclique>(scip);

      graph->addNode(0,1);
      graph->addNode(1,1);
      graph->flush();
      int edge[2] = {0,1};
      std::vector<int> hedge(edge, edge+2);

      graph->addHyperedge(hedge, 1);
      graph->setPartition(0, 1);
      graph->setPartition(1, 2);
      graph->flush();
   }

   virtual void TearDown() {
      SCIPfree(&scip);
      delete graph;
   }

protected:
   gcg::Hypergraph<GraphTclique> *graph;
   SCIP* scip;
};

TEST_F(GraphAlgorithmSmallCutTest, SmallCutSoed) {
   ASSERT_NEAR(2.0, GraphAlgorithms<GraphTclique>::computeSoed(*graph), 1e-6);
}
TEST_F(GraphAlgorithmSmallCutTest, SmallCutMincut) {
   ASSERT_NEAR(1.0, GraphAlgorithms<GraphTclique>::computeMincut(*graph), 1e-6);
}
TEST_F(GraphAlgorithmSmallCutTest, SmallCutKmetric) {
   ASSERT_NEAR(1.0, GraphAlgorithms<GraphTclique>::computekMetric(*graph), 1e-6);
}


class GraphAlgorithmMediumCutTest : public ::testing::Test {

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      graph = new Hypergraph<GraphTclique>(scip);

      graph->addNode(0,1);
      graph->addNode(1,1);
      graph->addNode(2,1);
      graph->addNode(3,1);
      graph->flush();

      int edge[] = {0,1};
      std::vector<int> hedge(edge, edge+2);
      graph->addHyperedge(hedge, 1);

      int edge2[] = {1,2};
      hedge = std::vector<int>(edge2, edge2+2);
      graph->addHyperedge(hedge, 1);

      int edge3[] = {2,3};
      hedge = std::vector<int>(edge3, edge3+2);
      graph->addHyperedge(hedge, 1);

      graph->setPartition(0, 1);
      graph->setPartition(1, 1);
      graph->setPartition(2, 2);
      graph->setPartition(3, 2);
      graph->flush();
   }

   virtual void TearDown() {
      SCIPfree(&scip);
      delete graph;
   }

protected:
   gcg::Hypergraph<GraphTclique> *graph;
   SCIP* scip;
};

TEST_F(GraphAlgorithmMediumCutTest, MediumCutSoed) {
   ASSERT_NEAR(2.0, GraphAlgorithms<GraphTclique>::computeSoed(*graph), 1e-6);
}
TEST_F(GraphAlgorithmMediumCutTest, MediumCutMincut) {
   ASSERT_NEAR(1.0, GraphAlgorithms<GraphTclique>::computeMincut(*graph), 1e-6);
}
TEST_F(GraphAlgorithmMediumCutTest, MediumCutKmetric) {
   ASSERT_NEAR(1.0, GraphAlgorithms<GraphTclique>::computekMetric(*graph), 1e-6);
}

class GraphAlgorithmMediumMultiCutTest : public ::testing::Test {

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      graph = new Hypergraph<GraphTclique>(scip);

      graph->addNode(0,1);
      graph->addNode(1,1);
      graph->addNode(2,1);
      graph->flush();

      int edge[2] = {0,1};
      std::vector<int> hedge(edge, edge+2);
      graph->addHyperedge(hedge, 1);

      int edge2[] = {0,1,2};
      hedge = std::vector<int>(edge2, edge2+3);
      graph->addHyperedge(hedge, 1);

      graph->setPartition(0, 1);
      graph->setPartition(1, 1);
      graph->setPartition(2, 2);
      graph->flush();
   }

   virtual void TearDown() {
      SCIPfree(&scip);
      delete graph;
   }

protected:
   gcg::Hypergraph<GraphTclique> *graph;
   SCIP* scip;
};

TEST_F(GraphAlgorithmMediumMultiCutTest, MediumMultiCutSoed) {
   ASSERT_NEAR(2.0, GraphAlgorithms<GraphTclique>::computeSoed(*graph), 1e-6);
}
TEST_F(GraphAlgorithmMediumMultiCutTest, MediumMultiCutMincut) {
   ASSERT_NEAR(1.0, GraphAlgorithms<GraphTclique>::computeMincut(*graph), 1e-6);
}
TEST_F(GraphAlgorithmMediumMultiCutTest, MediumMultiCutKmetric) {
   ASSERT_NEAR(1.0, GraphAlgorithms<GraphTclique>::computekMetric(*graph), 1e-6);
}

class GraphAlgorithmMediumMultiEdgeCutTest : public ::testing::Test {

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      graph = new Hypergraph<GraphTclique>(scip);

      graph->addNode(0,1);
      graph->addNode(1,1);
      graph->addNode(2,1);
      graph->addNode(3,1);
      graph->flush();

      int edge[2] = {0,1};
      std::vector<int> hedge(edge, edge+2);
      graph->addHyperedge(hedge, 1);

      int edge2[] = {0,1,2,3};
      hedge = std::vector<int>(edge2, edge2+4);
      graph->addHyperedge(hedge, 1);

      graph->setPartition(0, 1);
      graph->setPartition(1, 1);
      graph->setPartition(2, 2);
      graph->setPartition(2, 3);
      graph->flush();
   }

   virtual void TearDown() {
      SCIPfree(&scip);
      delete graph;
   }

protected:
   gcg::Hypergraph<GraphTclique> *graph;
   SCIP* scip;
};

TEST_F(GraphAlgorithmMediumMultiEdgeCutTest, MediumMultiEdgeCutSoed) {
   ASSERT_NEAR(3.0, GraphAlgorithms<GraphTclique>::computeSoed(*graph), 1e-6);
}
TEST_F(GraphAlgorithmMediumMultiEdgeCutTest, MediumMultiEdgeCutMincut) {
   ASSERT_NEAR(1.0, GraphAlgorithms<GraphTclique>::computeMincut(*graph), 1e-6);
}
TEST_F(GraphAlgorithmMediumMultiEdgeCutTest, MediumMultiEdgeCutKmetric) {
   ASSERT_NEAR(2.0, GraphAlgorithms<GraphTclique>::computekMetric(*graph), 1e-6);
}

class GraphAlgorithmWeigthedMulticutTest : public ::testing::Test {

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      graph = new Hypergraph<GraphTclique>(scip);

      graph->addNode(0,1);
      graph->addNode(1,1);
      graph->addNode(2,1);
      graph->addNode(3,1);
      graph->flush();

      int edge[2] = {0,1};
      std::vector<int> hedge(edge, edge+2);
      graph->addHyperedge(hedge, 1);

      int edge2[] = {0,1,2,3};
      hedge = std::vector<int>(edge2, edge2+4);
      graph->addHyperedge(hedge, 2);

      int edge3[] = {1,2};
      hedge = std::vector<int>(edge3, edge3+2);
      graph->addHyperedge(hedge, 4);

      graph->setPartition(0, 1);
      graph->setPartition(1, 1);
      graph->setPartition(2, 2);
      graph->setPartition(2, 3);
      graph->flush();
   }

   virtual void TearDown() {
      SCIPfree(&scip);
      delete graph;
   }

protected:
   gcg::Hypergraph<GraphTclique> *graph;
   SCIP* scip;
};

TEST_F(GraphAlgorithmWeigthedMulticutTest, WeigthedMulticutSoed) {
   ASSERT_NEAR(14.0, GraphAlgorithms<GraphTclique>::computeSoed(*graph), 1e-6);
}
TEST_F(GraphAlgorithmWeigthedMulticutTest, WeigthedMulticutMincut) {
   ASSERT_NEAR(6.0, GraphAlgorithms<GraphTclique>::computeMincut(*graph), 1e-6);
}
TEST_F(GraphAlgorithmWeigthedMulticutTest, WeigthedMulticutKmetric) {
   ASSERT_NEAR(8.0, GraphAlgorithms<GraphTclique>::computekMetric(*graph), 1e-6);
}

class GraphAlgorithmMSTTest : public ::testing::Test {

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      graph = new gcg::Graph<gcg::GraphGCG>(scip);

      graph->addNNodes(20);

      double w1 = 0.3;
      double w2 = 0.6;
      eps = 0.5;
      for(int i = 0; i< graph->getNNodes() -1; i++)
      {
         graph->addEdge(i, i+1, w1);

         if(i == 6 || i == 12)
         {
            graph->setEdge(i, i+1, w2);
         }
      }
      graph->addEdge(2, 17, w2);
      graph->addEdge(2, 8, w2);
      graph->addEdge(17, 8, w2);

   }

   virtual void TearDown() {
      SCIPfree(&scip);
      delete graph;
   }

protected:
   gcg::Graph<gcg::GraphGCG>* graph;
   SCIP* scip;
   double eps;
};

TEST_F(GraphAlgorithmMSTTest, MSTmainTest) {
   std::cout << "This is MST test..." << std::endl;
   std::vector<int> labels = GraphAlgorithms<gcg::GraphGCG>::mst(*graph, eps);
   for(auto label: labels)
   {
      std::cout << "Label = " << label << std::endl;
   }

   std::cout << "Total nodes: " << graph->getNNodes() << std::endl;
   for(int i = 0; i< graph->getNNodes(); i++)
   {
      auto ns = graph->getNeighborWeights(i);
      std::cout << "Node " << i << ": ";
      for(auto n: ns)
      {
         std::cout << n.first << ", ";
      }
      std::cout << "" << std::endl;
   }
   std::cout << "Now we print all the edges that are saved in the list...." << std::endl;
   std::vector<void*> edges;
   graph->getEdges(edges);
   for(int i = 0 ; i < (int)edges.size(); i++)
   {
      gcg::EdgeGCG next_edge = *(gcg::EdgeGCG *)(edges[i]);
      std::cout << "Edge: " << next_edge.src << ", " << next_edge.dest << std::endl;
   }

   std::cout << "Edges total: " << graph->getNEdges() << std::endl;
}


class GraphAlgorithmMCLTest : public ::testing::Test {

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );


      graph = new gcg::Graph<gcg::GraphGCG>(scip);

      graph->addNNodes(12);

      double w = 1.0;
      graph->addEdge(0, 1, w);
      graph->addEdge(0, 5, w);
      graph->addEdge(0, 6, w);
      graph->addEdge(0, 9, w);
      graph->addEdge(1, 2, w);
      graph->addEdge(1, 4, w);
      graph->addEdge(2, 3, w);
      graph->addEdge(2, 4, w);
      graph->addEdge(3, 7, w);
      graph->addEdge(3, 8, w);
      graph->addEdge(3, 10, w);
      graph->addEdge(4, 6, w);
      graph->addEdge(4, 7, w);
      graph->addEdge(5, 9, w);
      graph->addEdge(6, 9, w);
      graph->addEdge(7, 8, w);
      graph->addEdge(7, 10, w);
      graph->addEdge(8, 10, w);
      graph->addEdge(8, 11, w);
      graph->addEdge(10, 11, w);

      graph->flush();

   }

   virtual void TearDown() {
      SCIPfree(&scip);
      delete graph;
   }

protected:
   gcg::Graph<gcg::GraphGCG>* graph;
   SCIP* scip;
   double eps;
};


#ifdef WITH_GSL
TEST_F(GraphAlgorithmMCLTest, MCLmainTest) {
   std::cout << "This is MST test..." << std::endl;
   int inflate_fac = 2;
   std::vector<int> labels = GraphAlgorithms<gcg::GraphGCG>::mcl(*graph, inflate_fac);
   for(auto label: labels)
   {
      cout << "Label = " << label << endl;
   }

   /*std::cout << "Total nodes: " << graph->getNNodes() << std::endl;
   for(int i = 0; i< graph->getNNodes(); i++)
   {
      auto ns = graph->getNeighborWeights(i);
      std::cout << "Node " << i << ": ";
      for(auto n: ns)
      {
         std::cout << n.first << ", ";
      }
      std::cout << "" << std::endl;
   }
   std::cout << "Now we print all the edges that are saved in the list...." << std::endl;
   std::vector<void*> edges;
   graph->getEdges(edges);
   for(int i = 0 ; i < (int)edges.size(); i++)
   {
      gcg::EdgeGCG next_edge = *(gcg::EdgeGCG *)(edges[i]);
      std::cout << "Edge: " << next_edge.src << ", " << next_edge.dest << std::endl;
   }

   std::cout << "Edges total: " << graph->getNEdges() << std::endl;*/
}
#endif
