
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
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/arrowheur/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/borderheur/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/random/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/staircase/enabled", FALSE) );

     SCIP_CALL_ABORT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
     SCIP_CALL_ABORT( SCIPpresolve(scip) );
     SCIP_CALL_ABORT( DECdetectStructure(scip, &result) );
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
      SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/arrowheur/enabled", FALSE) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/borderheur/enabled", FALSE) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/random/enabled", FALSE) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/staircase/enabled", FALSE) );
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
   SCIP_CALL_EXPECT( SCIPsetBoolParam(scip, "detectors/borderheur/enabled", FALSE) );
   SCIP_CALL_EXPECT( SCIPsetBoolParam(scip, "detectors/arrowheur/enabled", FALSE) );
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
   SCIP_CALL_EXPECT( SCIPsetBoolParam(scip, "detectors/borderheur/enabled", FALSE) );
   SCIP_CALL_EXPECT( SCIPsetBoolParam(scip, "detectors/arrowheur/enabled", FALSE) );
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( DECdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
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
   int nconss;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C1W4_M.BPP.lp", "lp") );
   nconss = SCIPgetNConss(scip);
   SCIP_CALL_EXPECT( SCIPsetBoolParam(scip, "detectors/borderheur/enabled", FALSE) );
   SCIP_CALL_EXPECT( SCIPsetBoolParam(scip, "detectors/arrowheur/enabled", FALSE) );
   SCIP_CALL_EXPECT( SCIPpresolve(scip) );
   SCIP_CALL_EXPECT( DECdetectStructure(scip, &result) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_NEAR(41.0, SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)), SCIPfeastol(scip));

   SCIP_CALL_EXPECT( SCIPfreeSolve(scip, FALSE) );

   ASSERT_EQ(nconss+1, SCIPgetNConss(scip));
   ASSERT_EQ(SCIP_STAGE_TRANSFORMED, SCIPgetStage(scip));
   ASSERT_LE(1, SCIPconshdlrDecompGetNDecdecomps(scip));
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
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/arrowheur/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/borderheur/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/random/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/staircase/enabled", FALSE) );
   }

   virtual void TearDown() {
     SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
};

SCIP* GcgDecTest::scip = NULL;

TEST_F(GcgDecTest, ReadDecTest) {
   SCIP_RESULT result;
   DEC_DECOMP* decomp;
   int i;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/miplib/noswot.mps", "mps") );
   SCIP_CALL_EXPECT( SCIPreadDec(scip, "check/instances/miplib/noswot.dec", &result) );
   ASSERT_EQ(SCIP_SUCCESS, result);
   ASSERT_EQ(1, SCIPconshdlrDecompGetNDecdecomps(scip));

   decomp = SCIPconshdlrDecompGetDecdecomps(scip)[0];
   ASSERT_TRUE(decomp != NULL);
   EXPECT_EQ(5, DECdecompGetNBlocks(decomp));
   EXPECT_EQ(17, DECdecompGetNLinkingconss(decomp));
   EXPECT_EQ(3, DECdecompGetNLinkingvars(decomp));
   ASSERT_TRUE(DECdecompGetNSubscipconss(decomp) != NULL);

   for( i = 0; i < 5; ++i )
   {
      EXPECT_EQ(33, DECdecompGetNSubscipconss(decomp)[i]);
      EXPECT_EQ(25, DECdecompGetNSubscipvars(decomp)[i]);
   }
}

TEST_F(GcgDecTest, ReadBlkTest) {
   SCIP_RESULT result;
   DEC_DECOMP* decomp;
   int i;

   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   SCIP_CALL_EXPECT( SCIPreadBlk(scip, "check/instances/bpp/N1C3W1_A.blk", &result) );
   ASSERT_EQ(SCIP_SUCCESS, result);
   ASSERT_EQ(1, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_EXPECT( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   decomp = SCIPconshdlrDecompGetDecdecomps(scip)[0];
   ASSERT_TRUE(decomp != NULL);
   ASSERT_EQ(24, DECdecompGetNBlocks(decomp));
   ASSERT_EQ(50, DECdecompGetNLinkingconss(decomp));
   ASSERT_EQ(0, DECdecompGetNLinkingvars(decomp));
   ASSERT_TRUE(DECdecompGetNSubscipconss(decomp) != NULL);

   for( i = 0; i < 24; ++i )
   {
      ASSERT_EQ(1, DECdecompGetNSubscipconss(decomp)[i]);
      ASSERT_EQ(51, DECdecompGetNSubscipvars(decomp)[i]);
   }
}

TEST_F(GcgDecTest, NoDecTest) {
   DEC_DECOMP* decomp;

   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   ASSERT_EQ(0, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_EXPECT( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL_EXPECT( SCIPsetBoolParam(scip, "constraints/decomp/createbasicdecomp", 1) );
   SCIP_CALL_EXPECT( SCIPsetLongintParam(scip, "limits/nodes", 1L) );

   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_EQ(1, SCIPconshdlrDecompGetNDecdecomps(scip));
   ASSERT_NEAR(15.873333333333, SCIPgetLowerbound(scip), SCIPfeastol(scip));
   SCIP_CALL_EXPECT( SCIPsetBoolParam(scip, "constraints/decomp/createbasicdecomp", 0) );
   ASSERT_EQ(1, SCIPconshdlrDecompGetNDecdecomps(scip));
   SCIP_CALL_EXPECT( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   decomp = SCIPconshdlrDecompGetDecdecomps(scip)[0];
   ASSERT_TRUE(decomp != NULL);
   ASSERT_EQ(0, DECdecompGetNBlocks(decomp));
   ASSERT_EQ(SCIPgetNOrigConss(scip), DECdecompGetNLinkingconss(decomp));
   ASSERT_EQ(SCIPgetNOrigVars(scip), DECdecompGetNLinkingvars(decomp));
   ASSERT_TRUE(DECdecompGetNSubscipconss(decomp) == NULL);
}

TEST_F(GcgDecTest, WrongDecompTestBlk) {
   SCIP_RESULT result;
   SCIP_RETCODE retcode;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   retcode = SCIPreadBlk(scip, "check/instances/miplib/noswot.dec", &result);
   ASSERT_EQ(SCIP_READERROR, retcode);
}

TEST_F(GcgDecTest, WrongDecompTestDec) {
   SCIP_RESULT result;
   SCIP_RETCODE retcode;
   SCIP_CALL_EXPECT( SCIPreadProb(scip, "check/instances/bpp/N1C3W1_A.lp", "lp") );
   retcode = SCIPreadDec(scip, "check/instances/cpmp/p2050-1.txt.dec", &result);
   ASSERT_EQ(SCIP_READERROR, retcode);
}

TEST_F(GcgDecTest, MasterSpecificationTest) {
   SCIP_CONS** conss = NULL;
   DEC_DECOMP* decomp = NULL;
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

   SCIP_CALL_EXPECT(DECcreateDecompFromMasterconss(scip, &decomp, conss, 50) );

   ASSERT_TRUE(decomp != NULL);
   ASSERT_EQ(24, DECdecompGetNBlocks(decomp));
   ASSERT_EQ(50, DECdecompGetNLinkingconss(decomp));
   ASSERT_EQ(0, DECdecompGetNLinkingvars(decomp));
   ASSERT_TRUE(DECdecompGetNSubscipconss(decomp) != NULL);

   for( i = 0; i < 24; ++i )
   {
      ASSERT_EQ(1, DECdecompGetNSubscipconss(decomp)[i]);
      ASSERT_EQ(51, DECdecompGetNSubscipvars(decomp)[i]);
   }

   SCIP_CALL_EXPECT( DECdecompFree(scip, &decomp) );
   SCIPfreeMemoryArray(scip, &conss);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
