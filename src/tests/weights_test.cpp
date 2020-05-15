/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2020 Operations Research, RWTH Aachen University       */
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
