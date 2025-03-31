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

/**@file   hyperrowgraph_test.cpp
 * @brief  Unit tests for Column hypergraph
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "graph/hyperrowgraph.h"
#include "test.h"
#include "graphtest.h"
#include <fstream>
#include <algorithm>

class HyperrowTest : public GraphTest {

};


TEST_F(HyperrowTest, CreateTest) {
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[implicit] <x3>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continous] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] +1<x4>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x1>[I] +2<x2>[I] +3<x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x1>[I] +1<x3>[I] == 1") );
   gcg::Weights weights(1.0, 2, 3, 4, 5, 6);
   gcg::HyperrowGraph<gcg::GraphTclique> graph(scip, weights );

   ASSERT_EQ(SCIP_OKAY, graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );

}

TEST_F(HyperrowTest, WriteFileTest) {
   FILE* file = fopen("hypergraph.g", "wx");
   ASSERT_TRUE(file != NULL);
   int fd= fileno(file);
   ASSERT_NE(fd, -1);
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[implicit] <x3>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continous] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] +1<x4>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x1>[I] +2<x2>[I] +3<x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x1>[I] +1<x3>[I] == 1") );
   gcg::Weights weights(1, 2, 3, 4, 5, 6);
   gcg::HyperrowGraph<gcg::GraphTclique> graph(scip, weights );

   ASSERT_EQ(SCIP_OKAY, graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );
   ASSERT_EQ( SCIP_OKAY, graph.writeToFile(fd, 0) );
   fclose(file);
   ASSERT_TRUE( SCIPfileExists("hypergraph.g") );

   int tmp[] = {3,4,0,1,2,4,1,2,3,1,3};

   std::vector<int> array(&tmp[0], &tmp[0]+11);

   if( SCIPfileExists("hypergraph.g") )
   {
      parseFile("hypergraph.g", array);
      remove("hypergraph.g");
   }

}

TEST_F(HyperrowTest, WriteFileWeightsTest) {
   FILE* file = fopen("hypergraph.g", "wx");
   ASSERT_TRUE(file != NULL);
   int fd= fileno(file);
   ASSERT_NE(fd, -1);
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[implicit] <x3>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continous] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] +1<x4>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x1>[I] +2<x2>[I] +3<x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x1>[I] +1<x3>[I] == 1") );
   gcg::Weights weights(1, 2, 3, 4, 5, 6);
   gcg::HyperrowGraph<gcg::GraphTclique> graph(scip, weights );

   ASSERT_EQ(SCIP_OKAY, graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );
   ASSERT_EQ( SCIP_OKAY, graph.writeToFile(fd, 1) );
   fclose(file);
   ASSERT_TRUE( SCIPfileExists("hypergraph.g") );

   int tmp[] = {3,4,1,6,1,2,4,6,1,2,3,6,1,3};

   std::vector<int> array(&tmp[0], &tmp[0]+14);

   if( SCIPfileExists("hypergraph.g") )
   {
      parseFile("hypergraph.g", array);
      remove("hypergraph.g");
   }

}

TEST_F(HyperrowTest, ReadPartitionTest) {
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[implicit] <x3>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continous] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] +1<x4>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x1>[I] +2<x2>[I] +3<x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x1>[I] +1<x3>[I] == 1") );
   gcg::Weights weights(1, 2, 3, 4, 5, 6);
   gcg::HyperrowGraph<gcg::GraphTclique> graph(scip, weights );
   SCIP_CALL_EXPECT( graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );

   std::ofstream out;
   out.open("partition.part");
   for( int i = 0; i < graph.getNNodes(); ++i )
   {
      out << i << std::endl;
   }
   out.close();
   SCIP_CALL_EXPECT( graph.readPartition("partition.part") );

   std::vector<int> partition = graph.getPartition();

   for( int i = 0; i < graph.getNNodes(); ++i )
   {
      ASSERT_EQ(i, partition[i]);
   }

   remove("partition.part");
}

TEST_F(HyperrowTest, GetHyperedgeNodesTest) {
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[implicit] <x3>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continous] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] +1<x4>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x1>[I] +2<x2>[I] +3<x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x1>[I] +1<x3>[I] == 1") );
   gcg::Weights weights(1, 2, 3, 4, 5, 6);
   gcg::HyperrowGraph<gcg::GraphTclique> graph(scip, weights );
   SCIP_CALL_EXPECT( graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );

   int array[3][4] = {
      {0,1,3,0},
      {0,1,2,0},
      {0,2,0,0}
   };

   ASSERT_EQ(3, graph.getNEdges());
   ASSERT_EQ(4, graph.getNNodes());
   for( int i = 0; i < graph.getNEdges(); ++i )
   {
      std::vector<int> nodes = graph.getHyperedgeNodes(i);
      std::sort(nodes.begin(), nodes.end());
      for( size_t j = 0; j < nodes.size(); ++j)
      {
         ASSERT_EQ(array[i][j], nodes[j]);
      }
   }
}

TEST_F(HyperrowTest, GetNeighborTest) {
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[implicit] <x3>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continous] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] +1<x4>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x1>[I] +2<x2>[I] +3<x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x1>[I] +1<x3>[I] == 1") );
   gcg::Weights weights(1, 2, 3, 4, 5, 6);
   gcg::HyperrowGraph<gcg::GraphTclique> graph(scip, weights );
   SCIP_CALL_EXPECT( graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );

   int array[4][4] = {
      {1,2,3,0},
      {0,2,3,0},
      {0,1,0,0},
      {0,1,0,0}
   };

   ASSERT_EQ(3, graph.getNEdges());
   ASSERT_EQ(4, graph.getNNodes());
   for( int i = 0; i < graph.getNNodes(); ++i )
   {
      std::vector<int> nodes = graph.getNeighbors(i);
      std::sort(nodes.begin(), nodes.end());
      for( size_t j = 0; j < nodes.size(); ++j)
      {
         ASSERT_EQ(array[i][j], nodes[j]);
      }
   }
}
