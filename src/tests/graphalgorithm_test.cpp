/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2014 Operations Research, RWTH Aachen University       */
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

using gcg::GraphAlgorithms;
using gcg::GraphTclique;
using gcg::Hypergraph;

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
