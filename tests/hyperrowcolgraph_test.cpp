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

/**@file   hyperrowcolgraph_test.cpp
 * @brief  Unit tests for row-column hypergraph
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "graph/hyperrowcolgraph.h"
#include "test.h"
#include "graphtest.h"
#include <fstream>
#include <algorithm>
class HyperrowcolTest : public GraphTest {

};


TEST_F(HyperrowcolTest, CreateTest) {
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[implicit] <x3>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continous] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] +1<x4>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x1>[I] +2<x2>[I] +3<x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x1>[I] +1<x3>[I] == 1") );
   gcg::Weights weights(1.0, 2, 3, 4, 5, 6);
   gcg::HyperrowcolGraph<gcg::GraphTclique> graph(scip, weights );

   ASSERT_EQ(SCIP_OKAY, graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );

}

TEST_F(HyperrowcolTest, WriteFileTest) {
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
   gcg::HyperrowcolGraph<gcg::GraphTclique> graph(scip, weights );

   ASSERT_EQ(SCIP_OKAY, graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );
   ASSERT_EQ( SCIP_OKAY, graph.writeToFile(fd, 0) );
   fclose(file);
   ASSERT_TRUE( SCIPfileExists("hypergraph.g") );

   int tmp[] = {7, 8, 0, 1, 4, 7, 2, 5, 6, 8, 3, 1, 2, 3, 4, 5, 6, 7, 8};

   std::vector<int> array(&tmp[0], &tmp[0]+19);

   if( SCIPfileExists("hypergraph.g") )
   {
      parseFile("hypergraph.g", array);
      remove("hypergraph.g");
   }

}

TEST_F(HyperrowcolTest, WriteFileWeightsTest) {
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
   gcg::HyperrowcolGraph<gcg::GraphTclique> graph(scip, weights );

   ASSERT_EQ(SCIP_OKAY, graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );
   ASSERT_EQ( SCIP_OKAY, graph.writeToFile(fd, 1) );
   fclose(file);
   ASSERT_TRUE( SCIPfileExists("hypergraph.g") );

   int tmp[] = {7, 8, 1, 2, 1, 4, 7, 4, 2, 5, 5, 6, 8, 3, 3, 6, 1, 2, 3, 6, 4, 5, 6, 6, 7, 8};

   std::vector<int> array(&tmp[0], &tmp[0]+26);

   if( SCIPfileExists("hypergraph.g") )
   {
      parseFile("hypergraph.g", array);
      remove("hypergraph.g");
   }

}

TEST_F(HyperrowcolTest, ReadPartitionTest) {
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[implicit] <x3>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continous] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] +1<x4>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x1>[I] +2<x2>[I] +3<x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x1>[I] +1<x3>[I] == 1") );
   gcg::Weights weights(1, 2, 3, 4, 5, 6);
   gcg::HyperrowcolGraph<gcg::GraphTclique> graph(scip, weights );
   SCIP_CALL_EXPECT( graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );

   std::ofstream out;
   out.open("partition.part");
   for( int i = 0; i < 8; ++i )
   {
      out << i << std::endl;
   }
   out.close();
   SCIP_CALL_EXPECT( graph.readPartition("partition.part") );

   std::vector<int> partition = graph.getPartition();

   for( int i = 0; i < 8; ++i )
   {
      ASSERT_EQ(i, partition[i]);
   }

   remove("partition.part");
}

TEST_F(HyperrowcolTest, GetHyperedgeNodesTest) {
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[implicit] <x3>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continous] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] +1<x4>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x1>[I] +2<x2>[I] +3<x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x1>[I] +1<x3>[I] == 1") );
   gcg::Weights weights(1, 2, 3, 4, 5, 6);
   gcg::HyperrowcolGraph<gcg::GraphTclique> graph(scip, weights );
   SCIP_CALL_EXPECT( graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );

   int array[7][7] = {
      {0,3,6,0,0,0,0},
      {1,4,0,0,0,0,0},
      {5,7,0,0,0,0,0},
      {2,0,0,0,0,0,0},
      {0,1,2,0,0,0,0},
      {3,4,5,0,0,0,0},
      {6,7,0,0,0,0,0} };

   for( int i = 0; i < SCIPgetNVars(scip)+SCIPgetNConss(scip); ++i )
   {
      std::vector<int> nodes = graph.getHyperedgeNodes(i);
      std::sort(nodes.begin(), nodes.end());
      for( size_t j = 0; j < nodes.size(); ++j)
      {
         ASSERT_EQ(array[i][j], nodes[j]);
      }
   }
}

TEST_F(HyperrowcolTest, GetNeighborTest) {
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[implicit] <x3>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continous] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] +1<x4>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x1>[I] +2<x2>[I] +3<x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x1>[I] +1<x3>[I] == 1") );
   gcg::Weights weights(1, 2, 3, 4, 5, 6);
   gcg::HyperrowcolGraph<gcg::GraphTclique> graph(scip, weights );
   SCIP_CALL_EXPECT( graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );

   int array[8][8] = {
      {1,2,3,6,0,0,0,0},
      {0,2,4,0,0,0,0,0},
      {0,1,0,0,0,0,0,0},
      {0,4,5,6,0,0,0,0},
      {1,3,5,0,0,0,0,0},
      {3,4,7,0,0,0,0,0},
      {0,3,7,0,0,0,0,0},
      {5,6,0,0,0,0,0,0},
   };

   for( int i = 0; i < graph.getNNonzeroes(); ++i )
   {
      std::vector<int> nodes = graph.getNeighbors(i);
      std::sort(nodes.begin(), nodes.end());
      for( size_t j = 0; j < nodes.size(); ++j)
      {
         ASSERT_EQ(array[i][j], nodes[j]);
      }
   }
}
