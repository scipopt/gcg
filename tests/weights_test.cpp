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

/**@file   weights_test.cpp
 * @brief  Unit tests for Weights class
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "graph/weights.h"
#include "test.h"
class WeightTest : public ::testing::Test {
 protected:
  static SCIP *scip;

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(scip) );
      SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "name") );
   }

   virtual void TearDown() {
     SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
};

SCIP* WeightTest::scip = NULL;


TEST_F(WeightTest, BinaryTest) {
   SCIP_VAR* var;
   gcg::Weights Weights(1, 2, 1, 1, 1, 1);
   SCIP_CALL_EXPECT( SCIPcreateVarBasic(scip, &var, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );

   ASSERT_EQ(2, Weights.calculate(var));

   SCIP_CALL_EXPECT( SCIPreleaseVar(scip, &var) );
}

TEST_F(WeightTest, IntegerTest) {
   SCIP_VAR* var;
   gcg::Weights Weights(1, 1, 1, 2, 1, 1);
   SCIP_CALL_EXPECT( SCIPcreateVarBasic(scip, &var, "x1", 0.0, 3.0, 1.0, SCIP_VARTYPE_INTEGER) );

   ASSERT_EQ(2, Weights.calculate(var));

   SCIP_CALL_EXPECT( SCIPreleaseVar(scip, &var) );
}

TEST_F(WeightTest, ImplintTest) {
   SCIP_VAR* var;
   gcg::Weights Weights(1, 1, 1, 1, 2, 1);
   SCIP_CALL_EXPECT( SCIPcreateVarBasic(scip, &var, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_IMPLINT) );

   ASSERT_EQ(2, Weights.calculate(var));

   SCIP_CALL_EXPECT( SCIPreleaseVar(scip, &var) );
}

TEST_F(WeightTest, ContinousTest) {
   SCIP_VAR* var;
   gcg::Weights Weights(1, 1, 2, 1, 1, 1);
   SCIP_CALL_EXPECT( SCIPcreateVarBasic(scip, &var, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );

   ASSERT_EQ(2, Weights.calculate(var));

   SCIP_CALL_EXPECT( SCIPreleaseVar(scip, &var) );
}

TEST_F(WeightTest, ConsTest) {
   SCIP_CONS* cons;
   gcg::Weights Weights(1, 1, 1, 1, 1, 2);
   SCIP_CALL_EXPECT( SCIPcreateConsBasicLinear(scip, &cons, "c1", 0, NULL, NULL, 1.0, 1.0) );

   ASSERT_EQ(2, Weights.calculate(cons));

   SCIP_CALL_EXPECT( SCIPreleaseCons(scip, &cons) );
}
