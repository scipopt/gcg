#include "gtest/gtest.h"

#include "gcgplugins.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "cons_decomp.h"
#include "reader_blk.h"
#include "reader_dec.h"

#define SCIP_CALL_EXPECT(x) do { EXPECT_EQ(SCIP_OKAY, (x)); } while(FALSE)

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
     SCIP_CALL_ABORT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
     SCIP_CALL_ABORT( SCIPpresolve(scip) );
     SCIP_CALL_ABORT( DECdetectStructure(scip, &result) );
     SCIP_CALL_ABORT( SCIPsolve(scip) );
         
  }

  static void TearDownTestCase() {
     SCIP_CALL_ABORT( SCIPfree(&scip) );
  }

   virtual void SetUp() {

   }

   virtual void TearDown() {
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
   EXPECT_TRUE(GCGisMaster(GCGrelaxGetMasterprob(scip)));
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
   ASSERT_EQ(SCIP_OKAY, DECdetectStructure(scip, &result));
   ASSERT_EQ(SCIP_DIDNOTRUN, result);
}

TEST_F(GcgResultTest, numberOfBlocks) {
   ASSERT_EQ(50, GCGrelaxGetNPricingprobs(scip) );
   ASSERT_EQ(50, GCGrelaxGetNIdenticalBlocks(scip, 0));
} 

TEST_F(GcgResultTest, optimalSolutionValue) {
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));
}

TEST_F(GcgResultTest, performanceTest) {
   // expect solving time less than 5 seconds in debug mode
   EXPECT_GT(5, SCIPgetSolvingTime(scip)); 
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
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( DECdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   SCIP_CALL_EXPECT( SCIPfreeTransform(scip) );

   ASSERT_EQ(SCIP_STAGE_PROBLEM, SCIPgetStage(scip));
   ASSERT_EQ(0, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( DECdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   ASSERT_EQ(SCIP_STATUS_OPTIMAL, SCIPgetStatus(scip));
}

TEST_F(GcgLibTest, FreeProbTest) {
   SCIP_RESULT result;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( DECdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   SCIP_CALL_EXPECT( SCIPfreeProb(scip) );

   ASSERT_EQ(0, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( DECdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   ASSERT_EQ(SCIP_STATUS_OPTIMAL, SCIPgetStatus(scip));
}

TEST_F(GcgLibTest, FreeSolveTest) {
   SCIP_RESULT result;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( DECdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   SCIP_CALL_EXPECT( SCIPfreeSolve(scip, FALSE) );

   ASSERT_EQ(SCIP_STAGE_TRANSFORMED, SCIPgetStage(scip));
   ASSERT_LE(1, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
  // SCIP_CALL_EXPECT( DECdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   ASSERT_EQ(SCIP_STATUS_OPTIMAL, SCIPgetStatus(scip));
}

class GcgDecTest : public ::testing::Test {
 protected:
  static SCIP *scip;

  static void SetUpTestCase() {
  }

  static void TearDownTestCase() {

  }

   virtual void SetUp() {
     SCIP_CALL_ABORT( SCIPcreate(&scip) );
     SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
     SCIP_CALL_ABORT( SCIPcreateProb(scip, "test", NULL, NULL, NULL, NULL,NULL, NULL, NULL) );
     SCIP_CALL_ABORT( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
   }

   virtual void TearDown() {
     SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
};

SCIP* GcgDecTest::scip = NULL;

TEST_F(GcgDecTest, ReadDecTest) {
   SCIP_RESULT result;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/miplib/noswot.mps", "mps") );
   SCIP_CALL_EXPECT( SCIPreadDec(scip, "check/instances/miplib/noswot.dec", &result) );
   ASSERT_EQ(SCIP_SUCCESS, result);
   ASSERT_EQ(1, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_EXPECT( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_EQ(1, SCIPconshdlrDecompGetNDecdecomps(scip));
   ASSERT_NEAR(-41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));
}

TEST_F(GcgDecTest, ReadBlkTest) {
   SCIP_RESULT result;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   SCIP_CALL_EXPECT( SCIPreadBlk(scip, "check/instances/bpp/N1C3W1_A.blk", &result) );
   ASSERT_EQ(SCIP_SUCCESS, result);
   ASSERT_EQ(1, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_EXPECT( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_EQ(1, SCIPconshdlrDecompGetNDecdecomps(scip));
   ASSERT_NEAR(16.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));
}

TEST_F(GcgDecTest, NoDecTest) {
   SCIP_RESULT result;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   ASSERT_EQ(0, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_EXPECT( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL_EXPECT( SCIPsetLongintParam(scip, "limits/nodes", 1L) );

   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_EQ(1, SCIPconshdlrDecompGetNDecdecomps(scip));
   ASSERT_NEAR(15.873333333333, SCIPgetLowerbound(scip), SCIPfeastol(scip));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
