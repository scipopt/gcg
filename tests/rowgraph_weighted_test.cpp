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

/**@file   rowgraph_test.cpp
 * @brief  unit tests for row graph
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "graph/rowgraph_weighted.h"
#include "test.h"
#include "graphtest.h"

using namespace std;

class RowWeightedTest : public GraphTest
{

};


// Test the implementation of the GraphGCG and RowGraphWeighted (incl. similarity measures)
TEST_F(RowWeightedTest, test_createFromMatrix) {

   if( SCIPfileExists("rowWeightedGraph.g") )
   {
      remove("rowWeightedGraph.g");
   }

   //FILE* file = fopen("rowWeightedGraph.g", "wx");
   //ASSERT_TRUE(file != NULL);
   //int fd= fileno(file);
   //ASSERT_NE(fd, -1);
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x5>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x6>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x7>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x8>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x9>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x10>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x11>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x12>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x13>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x14>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c0>: 3<x1>[I] +5<x2>[I] +1<x3>[I] +6<x5>[I] +1<x6>[I] +1<x7>[I] +1<x8>[I] +1<x9>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 2<x1>[I] +4<x2>[I] +1<x6>[I] +2<x7>[I] +1<x14>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 1<x2>[I] +4<x3>[I] +1<x4>[I] +2<x8>[I] +1<x12>[I] +1<x13>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 3<x2>[I] +5<x4>[I] +1<x7>[I] +6<x8>[I] +1<x10>[I] +1<x11>[I] +1<x12>[I] +1<x14>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: 3<x1>[I] +5<x2>[I] +1<x4>[I] +6<x6>[I] +1<x7>[I] +1<x8>[I] +1<x10>[I] +1<x11>[I] +1<x12>[I] +1<x14>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c5>: 1<x12>[I] +1<x14>[I] <= 2") );

   gcg::Weights weights(1, 1, 1, 1, 1, 1);
   gcg::RowGraphWeighted<gcg::GraphGCG> graph(scip, weights );

   std::cout << "make a graph..." << std::endl;
   SCIP_CALL_EXPECT( graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip),
      gcg::DISTANCE_MEASURE::INTERSECTION, gcg::WEIGHT_TYPE::SIM));

   //ASSERT_EQ( SCIP_OKAY, graph.writeToFile(fd, FALSE) );
   //fclose(file);
   //ASSERT_TRUE( SCIPfileExists("rowWeightedGraph.g") );

}


#ifdef WITH_GSL
TEST_F(RowWeightedTest, GraphGCG_test) {

   gcg::GraphGCG graph;

   std::cout << "GraphGCG_test...." << std::endl;

   SCIP_CALL_EXPECT( graph.addNNodes(4) );


   graph.addEdge(0, 1, 1);
   graph.addEdge(0, 2, 1);
   graph.addEdge(0, 3, 1);
   graph.addEdge(3, 1, 1);
   graph.flush();
   graph.colL1Norm();
   graph.expand(2);
   cout << "expanded:" << endl;
   for(int i = 0; i < graph.getAdjMatrix()->size1; i++)
   {
      for(int j = 0; j < graph.getAdjMatrix()->size2; j++)
      {
         std::cout << "data["<< i <<","<< j<<"] = " << gsl_spmatrix_get(graph.getAdjMatrix(), i, j) << ",   ";
      }
      cout << endl;
      //std::cout << "data["<< i <<"] = " << graph.getAdjMatrix()->data[i] << std::endl;
   }
   graph.inflate(2);
   cout << "inflated:" << endl;
   for(int i = 0; i < graph.getAdjMatrix()->size1; i++)
   {
      for(int j = 0; j < graph.getAdjMatrix()->size2; j++)
      {
         std::cout << "data["<< i <<","<< j<<"] = " << gsl_spmatrix_get(graph.getAdjMatrix(), i, j) << ",   ";
      }
      cout << endl;
      //std::cout << "data["<< i <<"] = " << graph.getAdjMatrix()->data[i] << std::endl;
   }
   cout << "checking the neighbors....:" << endl;
   for (int col = 0; col < graph.getAdjMatrix()->size1; col++){
      size_t begin_ = graph.getAdjMatrix()->p[col];
      const size_t end_ = graph.getAdjMatrix()->p[col+1];
      vector<double> col_vals(end_ - begin_);
      vector<int> row_inds(end_ - begin_);
      size_t curr = 0;
      cout << "For column " << col << ", n neighbors: " << (int)(end_ - begin_) << endl;
      while(begin_ < end_){
         col_vals[curr] = graph.getAdjMatrix()->data[begin_];
         row_inds[curr++] = graph.getAdjMatrix()->i[begin_++];
      }
      assert(curr == row_inds.size());


      for(size_t it = 0; it < col_vals.size() ; it++ ){
         cout << " row: " << row_inds[it] << ", value = " << col_vals[it] << endl;
      }
   }
}


/*TEST_F(RowWeightedTest, GraphGCG_test_no_GSL) {

   gcg::GraphGCG graph;

   std::cout << "GraphGCG_test...." << std::endl;
   // was tested with 10k and it uses appr. 1.5GB
   int n_nodes = 100;

   SCIP_CALL_EXPECT( graph.addNNodes(n_nodes) );


   for(int i = 0; i < n_nodes; i++)
   {
      for(int j = 0; j < i; j++)
      {
         graph.addEdge(i, j, 1.0/3.0) ;
      }
   }
   graph.flush();

   //std::vector<std::vector<double>> mat = graph.getAdjMatrix();
   std::cout << "N edges: " << graph.getNEdges() << std::endl;
   for(int i = 0; i < graph.getNNodes(); i++)
   {
      int nneighbors = graph.getNNeighbors(i);
      std::cout << "i, nneighbors :" << i << ", " << nneighbors << std::endl;
      ASSERT_TRUE(nneighbors != 0);
      std::vector<int> neighbors = graph.getNeighbors(i);
      ASSERT_EQ(nneighbors, (int)neighbors.size());

      for( int j :  neighbors)
      {
         double weight = graph.getEdgeWeight(i,j);
         //std::cout << "(" << j << ", w=" << weight << ") ";
         ASSERT_NE(weight, 0.0);
      }
   }

}

TEST_F(RowWeightedTest, GraphGCG_test_GSL) {

   gcg::GraphGCG graph;

   std::cout << "New GraphGCG test...." << std::endl;
   // was tested with 10k and it uses appr. 1.5GB
   int n_nodes = 30;

   SCIP_CALL_EXPECT(graph.addNNodes(n_nodes) );

   int count = 0;
   for(int i = 0; i < n_nodes; i++)
   {
      for(int j = 0; j < n_nodes; j++)
      {
         if(count % 20 == 0 && i != j)
         {
            graph.addEdge(i, j, 1.0/3.0);
         }

         count++;
      }
   }
   graph.flush();

   //std::vector<std::vector<double>> mat = graph.getAdjMatrix();
   std::cout << "N edges: " << graph.getNEdges() << std::endl;
   for(int i = 0; i < graph.getNNodes(); i++)
   {
      std::cout << "node " << i << ": " << std::endl;
      int nneighbors = graph.getNNeighbors(i);
      //ASSERT_NE(nneighbors, 0);
      std::vector<std::pair<int,double> > neighbors = graph.getNeighborWeights(i);
      ASSERT_EQ(nneighbors, (int)neighbors.size());

      for( auto j :  neighbors)
      {
         //double weight = graph.getEdgeWeight(i,j);
         std::cout << "    (" << j.first << ", w=" << j.second << ") " <<std:: endl;
         //ASSERT_NE(weight, 0.0);
      }


      //std::cout << "OK node " << i << std::endl;
      std::cout << std::endl;
   }

}*/
#endif
