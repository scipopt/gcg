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

/**@file   aggregation_tests.cpp
 * @brief  Aggregation unit tests
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "test.h"
#include "gcgplugins.h"
#include "relax_gcg.h"
#include "gcg.h"
#include "pub_decomp.h"
#include "cons_decomp.h"

class GcgAggregationTest : public ::testing::Test {
 protected:
  static SCIP *scip;

   virtual void SetUp() {
     SCIP_CALL_ABORT( SCIPcreate(&scip) );
     SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
     SCIP_CALL_EXPECT( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );
     SCIP_CALL_EXPECT( SCIPsetIntParam(scip, "propagating/maxrounds", 0) );
     SCIP_CALL_ABORT( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrgpartition/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrcgpartition/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hcgpartition/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/random/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/staircase/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
     SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );
   }

   virtual void TearDown() {
     SCIP_CALL_ABORT( SCIPfree(&scip) );
   }

   SCIP_RETCODE createVar(const char * str) {
      SCIP_VAR* var;
      char* endptr;
      SCIP_Bool success;
      SCIP_CALL( SCIPparseVar(scip, &var, str, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL, &endptr, &success) );
      assert(success);
      SCIP_CALL( SCIPaddVar(scip, var) );
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
      return SCIP_OKAY;
   }

   SCIP_RETCODE createCons(const char * str) {
      SCIP_CONS* cons;
      SCIP_Bool success;
      SCIP_CALL( SCIPparseCons(scip, &cons, str, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      assert(success);
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      return SCIP_OKAY;
   }
};

SCIP* GcgAggregationTest::scip = NULL;

TEST_F(GcgAggregationTest, AggregateTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons;

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 2<x1>[I] +2<x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x3>[I] +2<x4>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] +<x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons = SCIPfindCons(scip, "c3");
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, &(mastercons), 1) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(2, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(2, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(0, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(FALSE, GCGisPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
}

TEST_F(GcgAggregationTest, WrongObjTest) {
   SCIP_RESULT result;
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 2<x1>[I] +2<x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x3>[I] +2<x4>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] +<x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( GCGdetectStructure(scip, &result) );
   ASSERT_EQ(SCIP_SUCCESS, result);

   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_EQ(2, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
}



TEST_F(GcgAggregationTest, WrongTypeTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons;
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continuous] <x4>: obj=2.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 2<x1>[I] +2<x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x3>[I] +2<x4>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] +<x3>[C] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons = SCIPfindCons(scip, "c3");
   assert(mastercons != NULL);
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, &(mastercons), 1) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_EQ(2, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
}

TEST_F(GcgAggregationTest, WrongBoundTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons;
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,2]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 1<x3>[I] +1<x4>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] +<x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons = SCIPfindCons(scip, "c3");
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, &(mastercons), 1) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(2, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
}

TEST_F(GcgAggregationTest, WrongCoeffSubproblemTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons;
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 2<x1>[I] +2<x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x3>[I] +4<x4>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] +<x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons = SCIPfindCons(scip, "c3");
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, &(mastercons), 1) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(2, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
}

TEST_F(GcgAggregationTest, WrongCoeffMasterTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons[2];
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,4]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,4]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 2<x1>[I] +2<x2>[I] >= 3") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x3>[I] +2<x4>[I] >= 3") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x1>[I] <= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: 1<x1>[I] +1<x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons[0] = SCIPfindCons(scip, "c3");
   mastercons[1] = SCIPfindCons(scip, "c4");
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, mastercons, 2) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(2, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
}

TEST_F(GcgAggregationTest, NonSetppcMasterTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons[2];
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,4]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,4]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 2<x1>[I] +2<x2>[I] >= 3") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x3>[I] +2<x4>[I] >= 3") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 2<x2>[I] +2<x4>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: 1<x1>[I] +1<x3>[I] <= 8") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons[0] = SCIPfindCons(scip, "c3");
   mastercons[1] = SCIPfindCons(scip, "c4");
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, mastercons, 2) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(2, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(2, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(0, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
   ASSERT_EQ(FALSE, GCGisPricingprobRelevant(scip, 1));
}

TEST_F(GcgAggregationTest, NonSetppcMasterWrongCoeffTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons[2];
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,4]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,4]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 2<x1>[I] +2<x2>[I] >= 3") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x3>[I] +2<x4>[I] >= 3") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 2<x2>[I] +3<x4>[I] <= 10") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: 1<x1>[I] +1<x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons[0] = SCIPfindCons(scip, "c3");
   mastercons[1] = SCIPfindCons(scip, "c4");
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, mastercons, 2) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(2, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
}

TEST_F(GcgAggregationTest, PresolvedMasterTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons[3];
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,2]") );

   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,2]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 1<x1>[I] +1<x3>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x2>[I] +1<x3>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: 1<x1>[I] +1<x2>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c5>: 1<x2>[I] +1<x3>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c6>: 1<x3>[I] +1<x1>[I] <= 2") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons[0] = SCIPfindCons(scip, "c4");
   mastercons[1] = SCIPfindCons(scip, "c5");
   mastercons[2] = SCIPfindCons(scip, "c6");
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, mastercons, 3) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(1, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
}

TEST_F(GcgAggregationTest, NonTriangleTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons[3];
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x5>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x6>: obj=2.0, original bounds=[0,2]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 1<x3>[I] +1<x4>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x5>[I] +1<x6>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: 1<x1>[I] +1<x3>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c5>: 1<x5>[I] +1<x3>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c6>: 1<x5>[I] +1<x1>[I] <= 2") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons[0] = SCIPfindCons(scip, "c4");
   mastercons[1] = SCIPfindCons(scip, "c5");
   mastercons[2] = SCIPfindCons(scip, "c6");
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, mastercons, 3) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(3, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 2));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 2));
}

TEST_F(GcgAggregationTest, NonExtendedTriangleTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons[3];
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x5>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x6>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x7>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x8>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x9>: obj=2.0, original bounds=[0,2]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 1<x3>[I] +1<x4>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x5>[I] +1<x6>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: 1<x1>[I] +1<x7>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c5>: 1<x3>[I] +1<x8>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c6>: 1<x5>[I] +1<x9>[I] <= 2") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons[0] = SCIPfindCons(scip, "c4");
   mastercons[1] = SCIPfindCons(scip, "c5");
   mastercons[2] = SCIPfindCons(scip, "c6");
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, mastercons, 3) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(3, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 2));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 2));
}

TEST_F(GcgAggregationTest, ExtendedMasterTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons[3];
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x5>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x6>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x7>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x8>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x9>: obj=2.0, original bounds=[0,2]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 1<x3>[I] +1<x4>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x5>[I] +1<x6>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: 1<x1>[I] +1<x3>[I] +1<x5>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c5>: 2<x1>[I] +2<x3>[I] +1<x7>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c6>: 1<x1>[I] +1<x3>[I] <= 2") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons[0] = SCIPfindCons(scip, "c4");
   mastercons[1] = SCIPfindCons(scip, "c5");
   mastercons[2] = SCIPfindCons(scip, "c6");
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, mastercons, 3) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(3, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(2, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(0, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 2));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
   ASSERT_EQ(FALSE, GCGisPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 2));
}

TEST_F(GcgAggregationTest, NonExtendedMasterTest) {
   GCG_DECOMP* decomp;
   SCIP_CONS* mastercons[3];
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x5>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x6>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x7>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x8>: obj=2.0, original bounds=[0,2]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x9>: obj=2.0, original bounds=[0,2]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 1<x1>[I] +1<x2>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 1<x3>[I] +1<x4>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 1<x5>[I] +1<x6>[I] >= 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: 1<x1>[I] +1<x3>[I] +1<x5>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c5>: 1<x1>[I] +1<x3>[I] <= 2") );
   SCIP_CALL_EXPECT( createCons("[linear] <c6>: 1<x1>[I] +1<x5>[I] <= 2") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons[0] = SCIPfindCons(scip, "c4");
   mastercons[1] = SCIPfindCons(scip, "c5");
   mastercons[2] = SCIPfindCons(scip, "c6");
   SCIP_CALL_EXPECT( GCGcreateDecompFromMasterconss(scip, &decomp, mastercons, 3) );
   SCIP_CALL_EXPECT( GCGconshdlrDecompAddDecomp(scip, decomp, false) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(3, GCGgetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(1, GCGgetNIdenticalBlocks(scip, 2));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 0));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGisPricingprobRelevant(scip, 2));
}
