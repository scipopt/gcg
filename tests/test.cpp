/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
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

/**@file   tests.cpp
 * @brief  generic unit tests
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "cons_decomp.h"
#include "reader_blk.h"
#include "reader_dec.h"
#include "pub_decomp.h"
#include "scip/cons_linear.h"
#include "gcg.h"

#include "test.h"

class GcgTest : public ::testing::Test {
 protected:
  static SCIP *scip;

  static void SetUpTestCase() {

       SCIP_CALL_ABORT( SCIPcreate(&scip) );
       SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
       SCIP_CALL_ABORT( SCIPcreateProb(scip, "test", NULL, NULL, NULL, NULL,NULL, NULL, NULL) );
       SCIP_CALL_ABORT( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
  }

  static void TearDownTestCase() {
     SCIP_CALL_ABORT( SCIPfree(&scip) );
  }

   virtual void SetUp() {

   }

   virtual void TearDown() {
      SCIPfreeTransform(scip);
   }
};

SCIP* GcgTest::scip = NULL;

class GcgResultTest : public ::testing::Test {
 protected:
  static SCIP *scip;

  static void SetUpTestCase() {
     SCIP_RESULT result;
     SCIP_CALL_ABORT( SCIPcreate(&scip) );
     SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
     SCIP_CALL_ABORT( SCIPcreateProb(scip, "test", NULL, NULL, NULL, NULL,NULL, NULL, NULL) );
     SCIP_CALL_ABORT( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrgpartition/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrcgpartition/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hcgpartition/enabled", FALSE) );

     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/staircase/enabled", FALSE) );


     SCIP_CALL_ABORT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
     SCIP_CALL_ABORT( SCIPpresolve(scip) );
     SCIP_CALL_ABORT( GCGdetectStructure(scip, &result) );
     SCIP_CALL_ABORT( SCIPsolve(scip) );

  }

  static void TearDownTestCase() {
     SCIP_CALL_ABORT( SCIPfree(&scip) );
  }

};

SCIP* GcgResultTest::scip = NULL;

TEST_F(GcgTest, StatusTest) {
   ASSERT_EQ(SCIP_STATUS_UNKNOWN, SCIPgetStatus(scip));
}

TEST_F(GcgTest, CreateTest) {
   ASSERT_FALSE(NULL == scip);
}

TEST_F(GcgTest, NameTest) {
   ASSERT_EQ(0, strcmp("test", SCIPgetProbName(scip)));
}

TEST_F(GcgTest, isGcgTest) {
   EXPECT_TRUE(GCGisOriginal(scip));
   EXPECT_TRUE(GCGisMaster(GCGgetMasterprob(scip)));
}


TEST_F(GcgTest, emptyProblem) {
   SCIP_SOL* bestsol;
   ASSERT_EQ(SCIP_OKAY, SCIPsolve(scip));
   ASSERT_EQ(0, SCIPgetNVars(scip));
   ASSERT_EQ(0, SCIPgetNConss(scip));
   ASSERT_EQ(SCIP_STATUS_OPTIMAL, SCIPgetStatus(scip));
   bestsol = SCIPgetBestSol(scip);
   ASSERT_TRUE(bestsol != NULL);
   ASSERT_FLOAT_EQ(0.0, SCIPgetSolTransObj(scip, bestsol));
}
TEST_F(GcgTest, detectEmptyProblem) {
   SCIP_RESULT result;
   ASSERT_EQ(SCIP_OKAY, GCGdetectStructure(scip, &result));
   ASSERT_EQ(SCIP_DIDNOTRUN, result);
}

TEST_F(GcgResultTest, numberOfBlocks) {
   ASSERT_EQ(50, GCGgetNPricingprobs(scip) );
}

TEST_F(GcgResultTest, optimalSolutionValue) {
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));
}

class GcgLibTest : public ::testing::Test {
 protected:
  static SCIP *scip;

  static void SetUpTestCase() {
  }

  static void TearDownTestCase() {
  }

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrgpartition/enabled", FALSE) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrcgpartition/enabled", FALSE) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hcgpartition/enabled", FALSE) );

      SCIP_CALL_ABORT( SCIPcreateProb(scip, "test", NULL, NULL, NULL, NULL,NULL, NULL, NULL) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
   }

   virtual void TearDown() {
     SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
};

SCIP* GcgLibTest::scip = NULL;

TEST_F(GcgLibTest, FreeTransformTest) {
   SCIP_RESULT result;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
   SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrgpartition/enabled", FALSE) );
    SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrcgpartition/enabled", FALSE) );
    SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hcgpartition/enabled", FALSE) );
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( GCGdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   SCIP_CALL_EXPECT( SCIPfreeTransform(scip) );

   ASSERT_EQ(SCIP_STAGE_PROBLEM, SCIPgetStage(scip));
   ASSERT_EQ(0, GCGconshdlrDecompGetNDecomps(scip));
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( GCGdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   ASSERT_EQ(SCIP_STATUS_OPTIMAL, SCIPgetStatus(scip));
}

TEST_F(GcgLibTest, FreeProbTest) {
   SCIP_RESULT result;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
   SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrgpartition/enabled", FALSE) );
    SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrcgpartition/enabled", FALSE) );
    SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hcgpartition/enabled", FALSE) );
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( GCGdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   SCIP_CALL_EXPECT( SCIPfreeProb(scip) );

   ASSERT_EQ(0, GCGconshdlrDecompGetNDecomps(scip));
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( GCGdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   ASSERT_EQ(SCIP_STATUS_OPTIMAL, SCIPgetStatus(scip));
}

TEST_F(GcgLibTest, FreeSolveTest) {
   SCIP_RESULT result;
   int nconss;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
   nconss = SCIPgetNConss(scip);
   SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrgpartition/enabled", FALSE) );
    SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrcgpartition/enabled", FALSE) );
    SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hcgpartition/enabled", FALSE) );
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( GCGdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   SCIP_CALL_EXPECT( SCIPfreeSolve(scip, FALSE) );

   ASSERT_EQ(nconss+1, SCIPgetNConss(scip));
   ASSERT_EQ(SCIP_STAGE_TRANSFORMED, SCIPgetStage(scip));
   ASSERT_LE(1, GCGconshdlrDecompGetNDecomps(scip));
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   ASSERT_EQ(SCIP_STATUS_OPTIMAL, SCIPgetStatus(scip));
   ASSERT_EQ(nconss+1, SCIPgetNConss(scip));
}
class GcgDecTest : public ::testing::Test {
 protected:
  static SCIP *scip;

   virtual void SetUp() {
     SCIP_CALL_ABORT( SCIPcreate(&scip) );
     SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
     SCIP_CALL_ABORT( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrgpartition/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrcgpartition/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hcgpartition/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/staircase/enabled", FALSE) );
   }

   virtual void TearDown() {
     SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
};

SCIP* GcgDecTest::scip = NULL;


TEST_F(GcgDecTest, ReadDecTest) {
   SCIP_RESULT result;
   GCG_DECOMP* decomp;
   int i;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/miplib/noswot.mps", "mps") );
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/miplib/noswot.dec", "dec") );
   /* ASSERT_EQ(SCIP_SUCCESS, result); */
   ASSERT_EQ(1, GCGconshdlrDecompGetNDecomps(scip));

   decomp = GCGconshdlrDecompGetDecomps(scip)[0];
   ASSERT_TRUE(decomp != NULL);
   EXPECT_EQ(5, GCGdecompGetNBlocks(decomp));
   EXPECT_EQ(17, GCGdecompGetNLinkingconss(decomp));
   EXPECT_EQ(3, GCGdecompGetNLinkingvars(decomp));
   ASSERT_TRUE(GCGdecompGetNSubscipconss(decomp) != NULL);

   for( i = 0; i < 5; ++i )
   {
      EXPECT_EQ(33, GCGdecompGetNSubscipconss(decomp)[i]);
      EXPECT_EQ(25, GCGdecompGetNSubscipvars(decomp)[i]);
   }
}

TEST_F(GcgDecTest, ReadBlkTest) {
   SCIP_RESULT result;
   GCG_DECOMP* decomp;
   int i;

   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   SCIP_CALL_EXPECT( SCIPreadBlk(scip, "check/instances/bpp/N1C3W1_A.blk", &result) );
   ASSERT_EQ(SCIP_SUCCESS, result);
   ASSERT_EQ(1, GCGconshdlrDecompGetNDecomps(scip));
   SCIP_CALL_EXPECT( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   decomp = GCGconshdlrDecompGetDecomps(scip)[0];
   ASSERT_TRUE(decomp != NULL);
   ASSERT_EQ(24, GCGdecompGetNBlocks(decomp));
   ASSERT_EQ(50, GCGdecompGetNLinkingconss(decomp));
   ASSERT_EQ(0, GCGdecompGetNLinkingvars(decomp));
   ASSERT_TRUE(GCGdecompGetNSubscipconss(decomp) != NULL);

   for( i = 0; i < 24; ++i )
   {
      ASSERT_EQ(1, GCGdecompGetNSubscipconss(decomp)[i]);
      ASSERT_EQ(51, GCGdecompGetNSubscipvars(decomp)[i]);
   }
}

TEST_F(GcgDecTest, NoDecTest) {
   GCG_DECOMP* decomp;

   SCIP_CALL_EXPECT( SCIPreadProb(scip, "lib/scip/check/instances/MIP/bell5.mps", "mps") );
   ASSERT_EQ(0, GCGconshdlrDecompGetNDecomps(scip));
   SCIP_CALL_EXPECT( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL_EXPECT( SCIPsetLongintParam(scip, "limits/nodes", 1L) );

   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_EQ(1, GCGconshdlrDecompGetNDecomps(scip));
   ASSERT_NEAR(+8.96640649152000e+06, SCIPgetLowerbound(scip), SCIPfeastol(scip));
   ASSERT_EQ(1, GCGconshdlrDecompGetNDecomps(scip));
   SCIP_CALL_EXPECT( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   decomp = GCGconshdlrDecompGetDecomps(scip)[0];
   ASSERT_TRUE(decomp != NULL);
   ASSERT_EQ(1, GCGdecompGetNBlocks(decomp));
   ASSERT_EQ(GCGdecompGetNLinkingconss(decomp), 0 );
}

TEST_F(GcgDecTest, WrongDecompTestBlk) {
   SCIP_RESULT result;
   SCIP_RETCODE retcode;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   retcode = SCIPreadBlk(scip, "check/instances/miplib/noswot.dec", &result);
   ASSERT_EQ(SCIP_READERROR, retcode);
}

TEST_F(GcgDecTest, WrongDecompTestDec) {
   SCIP_RETCODE retcode;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   retcode = SCIPreadProb(scip, "check/instances/cpmp/p2050-1.txt.dec", "dec");
   ASSERT_EQ(SCIP_READERROR, retcode);
}

TEST_F(GcgDecTest, MasterSpecificationTest) {
   SCIP_CONS** conss = NULL;
   GCG_DECOMP* decomp = NULL;
   int i = 0;
   char name[SCIP_MAXSTRLEN];

   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );

   SCIP_CALL_EXPECT( SCIPallocMemoryArray(scip, &conss, 50) );

   for( i = 0; i < 50; ++i )
   {
      SCIP_CONS* cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "Allocate_%d", i+1);
      cons = SCIPfindCons(scip, name);
      ASSERT_TRUE(cons != NULL);
      conss[i] = cons;
   }

   SCIP_CALL_EXPECT(GCGcreateDecompFromMasterconss(scip, &decomp, conss, 50) );

   ASSERT_TRUE(decomp != NULL);
   ASSERT_EQ(24, GCGdecompGetNBlocks(decomp));
   ASSERT_EQ(50, GCGdecompGetNLinkingconss(decomp));
   ASSERT_EQ(0, GCGdecompGetNLinkingvars(decomp));
   ASSERT_TRUE(GCGdecompGetNSubscipconss(decomp) != NULL);

   for( i = 0; i < 24; ++i )
   {
      ASSERT_EQ(1, GCGdecompGetNSubscipconss(decomp)[i]);
      ASSERT_EQ(51, GCGdecompGetNSubscipvars(decomp)[i]);
   }

   SCIP_CALL_EXPECT( GCGdecompFree(scip, &decomp) );
   SCIPfreeMemoryArray(scip, &conss);
}

TEST_F(GcgDecTest, EqualDecTest) {
   SCIP_CONS** conss = NULL;
   GCG_DECOMP* decomp1 = NULL;
   GCG_DECOMP* decomp2 = NULL;
   GCG_DECOMP* decomp3 = NULL;
   GCG_DECOMP* decomp4 = NULL;

   int i = 0;
   char name[SCIP_MAXSTRLEN];

   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );

   SCIP_CALL_EXPECT( SCIPallocMemoryArray(scip, &conss, 50) );

   for( i = 0; i < 50; ++i )
   {
      SCIP_CONS* cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "Allocate_%d", i+1);
      cons = SCIPfindCons(scip, name);
      ASSERT_TRUE(cons != NULL);
      conss[i] = cons;
   }

   SCIP_CALL_EXPECT(GCGcreateDecompFromMasterconss(scip, &decomp1, conss, 50) );

   for( i = 0; i < 24; ++i )
   {
      SCIP_CONS* cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "Capacity_%d", i+1);
      cons = SCIPfindCons(scip, name);
      ASSERT_TRUE(cons != NULL);
      conss[i] = cons;
   }

   SCIP_CALL_EXPECT(GCGcreateDecompFromMasterconss(scip, &decomp2, conss, 24) );
   SCIP_CALL_EXPECT(GCGcreateDecompFromMasterconss(scip, &decomp3, &conss[1], 1) );
   SCIP_CALL_EXPECT(GCGcreateDecompFromMasterconss(scip, &decomp4, &conss[0], 1) );

   ASSERT_TRUE(GCGdecompositionsAreEqual(scip, decomp1, decomp1));
   ASSERT_TRUE(GCGdecompositionsAreEqual(scip, decomp2, decomp2));
   ASSERT_FALSE(GCGdecompositionsAreEqual(scip, decomp2, decomp1));
   ASSERT_FALSE(GCGdecompositionsAreEqual(scip, decomp1, decomp2));
   ASSERT_FALSE(GCGdecompositionsAreEqual(scip, decomp3, decomp4));
   ASSERT_FALSE(GCGdecompositionsAreEqual(scip, decomp1, decomp3));

   SCIPfreeMemoryArray(scip, &conss);

   SCIP_CALL_EXPECT( GCGdecompFree(scip, &decomp1) );
   SCIP_CALL_EXPECT( GCGdecompFree(scip, &decomp2) );
   SCIP_CALL_EXPECT( GCGdecompFree(scip, &decomp3) );
   SCIP_CALL_EXPECT( GCGdecompFree(scip, &decomp4) );
}

TEST_F(GcgDecTest, FilterDecTest) {
   SCIP_CONS** conss = NULL;
   GCG_DECOMP* decomp1 = NULL;
   GCG_DECOMP* decomp2 = NULL;
   GCG_DECOMP* decomp3 = NULL;
   GCG_DECOMP* decomp4 = NULL;

   int i = 0;
   char name[SCIP_MAXSTRLEN];

   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );

   SCIP_CALL_EXPECT( SCIPallocMemoryArray(scip, &conss, 50) );

   for( i = 0; i < 50; ++i )
   {
      SCIP_CONS* cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "Allocate_%d", i+1);
      cons = SCIPfindCons(scip, name);
      ASSERT_TRUE(cons != NULL);
      conss[i] = cons;
   }

   SCIP_CALL_EXPECT(GCGcreateDecompFromMasterconss(scip, &decomp1, conss, 50) );

   for( i = 0; i < 24; ++i )
   {
      SCIP_CONS* cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "Capacity_%d", i+1);
      cons = SCIPfindCons(scip, name);
      ASSERT_TRUE(cons != NULL);
      conss[i] = cons;
   }

   SCIP_CALL_EXPECT(GCGcreateDecompFromMasterconss(scip, &decomp2, conss, 24) );
   SCIP_CALL_EXPECT(GCGcreateDecompFromMasterconss(scip, &decomp3, &conss[1], 1) );
   SCIP_CALL_EXPECT(GCGcreateDecompFromMasterconss(scip, &decomp4, &conss[0], 1) );

   GCG_DECOMP* decomps[5] = {decomp1, decomp2, decomp1, decomp3, decomp4};

   ASSERT_EQ(1, GCGfilterSimilarDecompositions(scip, decomps, 1) );
   ASSERT_EQ(decomp1, decomps[0]);
   ASSERT_EQ(2, GCGfilterSimilarDecompositions(scip, decomps, 2) );
   ASSERT_EQ(decomp2, decomps[1]);
   ASSERT_EQ(2, GCGfilterSimilarDecompositions(scip, decomps, 3) );
   ASSERT_EQ(decomp1, decomps[2]);
   ASSERT_EQ(3, GCGfilterSimilarDecompositions(scip, decomps, 4) );
   ASSERT_EQ(decomp3, decomps[2]);
   ASSERT_EQ(decomp1, decomps[3]);
   ASSERT_EQ(4, GCGfilterSimilarDecompositions(scip, decomps, 5) );
   ASSERT_EQ(decomp3, decomps[2]);
   ASSERT_EQ(decomp4, decomps[3]);
   ASSERT_EQ(decomp1, decomps[4]);

   SCIPfreeMemoryArray(scip, &conss);

   SCIP_CALL_EXPECT( GCGdecompFree(scip, &decomp1) );
   SCIP_CALL_EXPECT( GCGdecompFree(scip, &decomp2) );
   SCIP_CALL_EXPECT( GCGdecompFree(scip, &decomp3) );
   SCIP_CALL_EXPECT( GCGdecompFree(scip, &decomp4) );
}

class GcgMultProbsTest : public ::testing::Test {
 protected:
  static SCIP *scipone;
  static SCIP *sciptwo;

   virtual void SetUp() {
      SCIP_CALL_ABORT( SCIPcreate(&scipone) );
      SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scipone) );
      SCIP_CALL_ABORT( SCIPcreateProb(scipone, "testone", NULL, NULL, NULL, NULL,NULL, NULL, NULL) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipone, "display/verblevel", SCIP_VERBLEVEL_NONE) );

      SCIP_CALL_ABORT( SCIPcreate(&sciptwo) );
      SCIP_CALL_ABORT( SCIPincludeGcgPlugins(sciptwo) );
      SCIP_CALL_ABORT( SCIPfree(&sciptwo) );

   }

   virtual void TearDown() {
     SCIP_CALL_ABORT( SCIPfree(&scipone) );
   }
};

SCIP* GcgMultProbsTest::scipone = NULL;
SCIP* GcgMultProbsTest::sciptwo = NULL;

TEST_F(GcgMultProbsTest, FreeTransformTest) {
   SCIP_RESULT result;

   SCIP_CALL_EXPECT( SCIPreadProb(scipone, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
   SCIP_CALL_EXPECT( SCIPpresolve(scipone) );
   SCIP_CALL_EXPECT( GCGdetectStructure(scipone, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scipone) );
}



int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
