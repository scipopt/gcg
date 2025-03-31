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

#include "graph/rowgraph.h"
#include "test.h"
#include "graphtest.h"

class RowTest : public GraphTest
{

};

TEST_F(RowTest, WriteFileTest) {
   FILE* file = fopen("rowgraph.g", "wx");
   ASSERT_TRUE(file != NULL);
   int fd= fileno(file);
   ASSERT_NE(fd, -1);
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=1.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=1.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] +1<x3>[I]<= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x1>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x3>[I] == 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: 1<x1>[I] +1<x2>[I] == 1") );
   gcg::Weights weights(1.0, 2, 3, 4, 5, 6);
   gcg::RowGraph<gcg::GraphTclique> graph(scip, weights );

   SCIP_CALL_EXPECT( graph.createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );

   ASSERT_EQ( SCIP_OKAY, graph.writeToFile(fd, FALSE) );
   fclose(file);
   ASSERT_TRUE( SCIPfileExists("rowgraph.g") );

   int tmp[] = {4, 4, 2, 3, 4, 1, 4, 1, 1, 2};

   std::vector<int> array(&tmp[0], &tmp[0]+10);

   if( SCIPfileExists("rowgraph.g") )
   {
      parseFile("rowgraph.g", array);
      remove("rowgraph.g");
   }
}
