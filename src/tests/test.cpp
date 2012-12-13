#include "gtest/gtest.h"

#include "gcgplugins.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "cons_decomp.h"

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
     SCIP_RESULT result;
     SCIP_CALL_ABORT( SCIPcreate(&scip) );
     SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
     SCIP_CALL_ABORT( SCIPcreateProb(scip, "test", NULL, NULL, NULL, NULL,NULL, NULL, NULL) );
     SCIP_CALL_ABORT( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/borderheur/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/arrowheur/enabled", FALSE) );
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

SCIP* GcgLibTest::scip = NULL;


TEST_F(GcgLibTest, FreeTransformTest) {
   SCIP_RESULT result;
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));
   SCIP_CALL_ABORT( SCIPfreeTransform(scip) );
   ASSERT_EQ(SCIP_STAGE_PROBLEM, SCIPgetStage(scip));
   ASSERT_EQ(0, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_ABORT( SCIPpresolve(scip) );
   SCIP_CALL_ABORT( DECdetectStructure(scip, &result) );
   SCIP_CALL_ABORT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   ASSERT_EQ(SCIP_STATUS_OPTIMAL, SCIPgetStatus(scip));
}

TEST_F(GcgLibTest, FreeProbTest) {
   SCIP_RESULT result;
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));
   SCIP_CALL_ABORT( SCIPfreeProb(scip) );
   ASSERT_EQ(0, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_ABORT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
   SCIP_CALL_ABORT( SCIPpresolve(scip) );
   SCIP_CALL_ABORT( DECdetectStructure(scip, &result) );
   SCIP_CALL_ABORT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   ASSERT_EQ(SCIP_STATUS_OPTIMAL, SCIPgetStatus(scip));
}

TEST_F(GcgLibTest, FreeSolveTest) {
   SCIP_RESULT result;
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));
   SCIP_CALL_ABORT( SCIPfreeSolve(scip, FALSE) );
   ASSERT_EQ(SCIP_STAGE_TRANSFORMED, SCIPgetStage(scip));
   ASSERT_LE(1, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_ABORT( SCIPpresolve(scip) );
  // SCIP_CALL_ABORT( DECdetectStructure(scip, &result) );
   SCIP_CALL_ABORT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   ASSERT_EQ(SCIP_STATUS_OPTIMAL, SCIPgetStatus(scip));
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
