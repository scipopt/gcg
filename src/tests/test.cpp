#include "gtest/gtest.h"

#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "cons_decomp.h"
#include "reader_blk.h"
#include "reader_dec.h"
#include "pub_decomp.h"
#include "scip/cons_linear.h"

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

  static void SetUpTestCase() {
  }

  static void TearDownTestCase() {

  }

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

class GcgDecStatisticTest : public ::testing::Test {
 protected:
  static SCIP *scip;
  SCIP_VAR* vars[5];
  SCIP_VAR* transvars[5];
  SCIP_CONS* conss[3];
  SCIP_CONS* transconss[3];
  SCIP_Real curvals[3];
  SCIP_VAR* curvars[3];
  DEC_DECOMP* decomp;

  static void SetUpTestCase() {
  }

  static void TearDownTestCase() {

  }

   virtual void SetUp() {
      int i;
      SCIP_HASHMAP* constoblock;
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );

      SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );
     /*
      * create problem
      *
      * min 3*x1 + x2 + 3*x3 + x4 + 3*x5
      * s.t.
      *
      * x1 -x2       +x5  = 1
      *       -x3+x4 -x5 <= 1
      * x1       -x4 +x5 >= 2
      */

      SCIP_CALL_ABORT( SCIPcreateVarBasic(scip, &(vars[0]), "x1", 0.0, 3.0, 1.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL_ABORT( SCIPcreateVarBasic(scip, &(vars[1]), "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL_ABORT( SCIPcreateVarBasic(scip, &(vars[2]), "x3", 0.0, 3.0, 1.0, SCIP_VARTYPE_IMPLINT) );
      SCIP_CALL_ABORT( SCIPcreateVarBasic(scip, &(vars[3]), "x4", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL_ABORT( SCIPcreateVarBasic(scip, &(vars[4]), "x5", 0.0, 3.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );

      curvars[0] = vars[0];
      curvars[1] = vars[1];
      curvars[2] = vars[4];
      curvals[0] = 1.0;
      curvals[1] = -1.0;
      curvals[2] = 1.0;
      SCIP_CALL_ABORT( SCIPcreateConsBasicLinear(scip, &(conss[0]), "c1", 3, curvars, curvals, 1.0, 1.0) );

      curvars[0] = vars[2];
      curvars[1] = vars[3];
      curvars[2] = vars[4];
      curvals[0] = -1.0;
      curvals[1] = 1.0;
      curvals[2] = -1.0;
      SCIP_CALL_ABORT( SCIPcreateConsBasicLinear(scip, &(conss[1]), "c2", 3, curvars, curvals, -SCIPinfinity(scip), 1.0) );

      curvars[0] = vars[0];
      curvars[1] = vars[3];
      curvars[2] = vars[4];
      curvals[0] = 1.0;
      curvals[1] = -1.0;
      curvals[2] = 1.0;
      SCIP_CALL_ABORT( SCIPcreateConsBasicLinear(scip, &(conss[2]), "c3", 3, curvars, curvals, 1.0, SCIPinfinity(scip)) );

      for( i = 0; i < 5; ++i )
      {
         SCIP_CALL_ABORT( SCIPaddVar(scip, vars[i]));
      }

      for( i = 0; i < 3; ++i )
      {
         SCIP_CALL_ABORT( SCIPaddCons(scip, conss[i]));
      }

      SCIP_CALL_ABORT( SCIPtransformProb(scip) );

      SCIP_CALL_ABORT( SCIPgetTransformedConss(scip, 3, conss, transconss) );
      SCIP_CALL_ABORT( SCIPgetTransformedVars(scip, 5, vars, transvars) );

      SCIP_CALL_ABORT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 3) );
      SCIP_CALL_ABORT( SCIPhashmapInsert(constoblock, transconss[0], (void*) (size_t)1) );
      SCIP_CALL_ABORT( SCIPhashmapInsert(constoblock, transconss[1], (void*) (size_t)2) );
      SCIP_CALL_ABORT( SCIPhashmapInsert(constoblock, transconss[2], (void*) (size_t)3) );
      SCIP_CALL_ABORT( DECdecompCreate(scip, &decomp) );

      SCIP_CALL_ABORT( DECfilloutDecdecompFromConstoblock(scip, decomp, constoblock, 2, transvars, 5, transconss, 3, FALSE) );
   }

   virtual void TearDown() {
      int i;
      for( i = 0; i < 5; ++i )
      {
         SCIP_CALL_ABORT( SCIPreleaseVar(scip, &(vars[i])));
      }

      for( i = 0; i < 3; ++i )
      {
         SCIP_CALL_ABORT( SCIPreleaseCons(scip, &(conss[i])));
      }
      SCIP_CALL_ABORT( DECdecompFree(scip, &decomp) );


      SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
};

SCIP* GcgDecStatisticTest::scip = NULL;

TEST_F(GcgDecStatisticTest, BlockTest) {
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
}

TEST_F(GcgDecStatisticTest, SubscipSizeTest) {
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
   ASSERT_EQ(2, DECdecompGetNSubscipvars(decomp)[0]);
   ASSERT_EQ(2, DECdecompGetNSubscipvars(decomp)[1]);
   ASSERT_EQ(1, DECdecompGetNLinkingvars(decomp));
   ASSERT_EQ(1, DECdecompGetNSubscipconss(decomp)[0]);
   ASSERT_EQ(1, DECdecompGetNSubscipconss(decomp)[0]);
   ASSERT_EQ(1, DECdecompGetNLinkingconss(decomp));
}

TEST_F(GcgDecStatisticTest, DensityTest) {
   int i;
   SCIP_VAR* densvars[5];
   SCIP_CONS* densconss[5];
   int varsubprobdens[5];
   int varmasterdens[5];
   int conssubprobdens[3];
   int consmasterdens[3];

   SCIP_CALL_EXPECT(DECgetDensityData(scip, decomp, densvars, 5, densconss, 3, varsubprobdens, varmasterdens, conssubprobdens, consmasterdens) );
   for( i = 0; i < 5; i++)
   {
      if( strcmp(SCIPvarGetName(densvars[i]), "t_x1") == 0)
      {
         ASSERT_EQ(1, varsubprobdens[i]);
         ASSERT_EQ(1, varmasterdens[i]);
      }
      else if( strcmp(SCIPvarGetName(densvars[i]), "t_x2") == 0)
      {
         ASSERT_EQ(1, varsubprobdens[i]);
         ASSERT_EQ(0, varmasterdens[i]);
      }
      else if( strcmp(SCIPvarGetName(densvars[i]), "t_x3") == 0)
      {
         ASSERT_EQ(1, varsubprobdens[i]);
         ASSERT_EQ(0, varmasterdens[i]);
      }
      else if( strcmp(SCIPvarGetName(densvars[i]), "t_x4") == 0)
      {
         ASSERT_EQ(1, varsubprobdens[i]);
         ASSERT_EQ(1, varmasterdens[i]);
      }
      else if( strcmp(SCIPvarGetName(densvars[i]), "t_x5") == 0)
      {
         ASSERT_EQ(2, varsubprobdens[i]);
         ASSERT_EQ(1, varmasterdens[i]);
      }
      else
      {
         ASSERT_FALSE(true);
      }
   }
   for( i = 0; i < 3; i++)
   {
      if( strcmp(SCIPconsGetName(densconss[i]), "c1") == 0)
      {
         ASSERT_EQ(2, conssubprobdens[i]);
         ASSERT_EQ(1, consmasterdens[i]);
      }
      else if( strcmp(SCIPconsGetName(densconss[i]), "c2") == 0)
      {
         ASSERT_EQ(2, conssubprobdens[i]);
         ASSERT_EQ(1, consmasterdens[i]);
      }
      else if( strcmp(SCIPconsGetName(densconss[i]), "c3") == 0)
      {
         ASSERT_EQ(0, conssubprobdens[i]);
         ASSERT_EQ(3, consmasterdens[i]);
      }
      else
      {
         ASSERT_FALSE(true);
      }
   }
}

TEST_F(GcgDecStatisticTest, VarsDataTest) {
   int nvars[2];
   int nbinvars[2];
   int nintvars[2];
   int nimplvars[2];
   int ncontvars[2];

   DECgetSubproblemVarsData(scip, decomp, nvars, nbinvars, nintvars, nimplvars, ncontvars, 2);
   ASSERT_EQ(2, nvars[0]);
   ASSERT_EQ(2, nvars[1]);
   ASSERT_EQ(1, nintvars[0]);
   ASSERT_EQ(0, nintvars[1]);
   ASSERT_EQ(1, nbinvars[0]);
   ASSERT_EQ(1, nbinvars[1]);
   ASSERT_EQ(0, nimplvars[0]);
   ASSERT_EQ(1, nimplvars[1]);
   ASSERT_EQ(0, ncontvars[0]);
   ASSERT_EQ(0, ncontvars[1]);

   DECgetLinkingVarsData(scip, decomp, nvars, nbinvars, nintvars, nimplvars, ncontvars);
   ASSERT_EQ(1, nvars[0]);
   ASSERT_EQ(0, nintvars[0]);
   ASSERT_EQ(0, nbinvars[0]);
   ASSERT_EQ(0, nimplvars[0]);
   ASSERT_EQ(1, ncontvars[0]);
}

TEST_F(GcgDecStatisticTest, VarlockTest) {
   SCIP_VAR* lockvars[5];
   int* sublockdown[2];
   int* sublockup[2];
   int masterlockdown[5];
   int masterlockup[5];
   int i;

   for( i = 0; i < 2; ++i)
   {
      sublockdown[i] = new int[5];
      sublockup[i] = new int[5];
   }

   SCIP_CALL_EXPECT( DECgetVarLockData(scip, decomp, lockvars, 5, 2, sublockdown, sublockup, masterlockdown, masterlockup) );
   for( i = 0; i < 5; i++)
   {
      if( strcmp(SCIPvarGetName(lockvars[i]), "t_x1") == 0)
      {
         ASSERT_EQ(1, sublockdown[0][i]);
         ASSERT_EQ(0, sublockdown[1][i]);
         ASSERT_EQ(1, sublockup[0][i]);
         ASSERT_EQ(0, sublockup[1][i]);
         ASSERT_EQ(1, masterlockdown[i]);
         ASSERT_EQ(0, masterlockup[i]);
      }
      else if( strcmp(SCIPvarGetName(lockvars[i]), "t_x2") == 0)
      {
         ASSERT_EQ(1, sublockdown[0][i]);
         ASSERT_EQ(0, sublockdown[1][i]);
         ASSERT_EQ(1, sublockup[0][i]);
         ASSERT_EQ(0, sublockup[1][i]);
         ASSERT_EQ(0, masterlockdown[i]);
         ASSERT_EQ(0, masterlockup[i]);
      }
      else if( strcmp(SCIPvarGetName(lockvars[i]), "t_x3") == 0)
      {
         ASSERT_EQ(0, sublockdown[0][i]);
         ASSERT_EQ(1, sublockdown[1][i]);
         ASSERT_EQ(0, sublockup[0][i]);
         ASSERT_EQ(0, sublockup[1][i]);
         ASSERT_EQ(0, masterlockdown[i]);
         ASSERT_EQ(0, masterlockup[i]);
      }
      else if( strcmp(SCIPvarGetName(lockvars[i]), "t_x4") == 0)
      {
         ASSERT_EQ(0, sublockdown[0][i]);
         ASSERT_EQ(0, sublockdown[1][i]);
         ASSERT_EQ(0, sublockup[0][i]);
         ASSERT_EQ(1, sublockup[1][i]);
         ASSERT_EQ(0, masterlockdown[i]);
         ASSERT_EQ(1, masterlockup[i]);
      }
      else if( strcmp(SCIPvarGetName(lockvars[i]), "t_x5") == 0)
      {
         ASSERT_EQ(1, sublockdown[0][i]);
         ASSERT_EQ(1, sublockdown[1][i]);
         ASSERT_EQ(1, sublockup[0][i]);
         ASSERT_EQ(0, sublockup[1][i]);
         ASSERT_EQ(1, masterlockdown[i]);
         ASSERT_EQ(0, masterlockup[i]);
      }
      else
      {
         ASSERT_FALSE(true);
      }
   }
   for( i = 0; i < 2; ++i)
   {
      delete[] sublockdown[i];
      delete[] sublockup[i];
   }
}

class GcgAggregationTest : public ::testing::Test {
 protected:
  static SCIP *scip;

  static void SetUpTestCase() {
  }

  static void TearDownTestCase() {

  }

   virtual void SetUp() {
     SCIP_CALL_ABORT( SCIPcreate(&scip) );
     SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
     SCIP_CALL_ABORT( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/arrowheur/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/borderheur/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/random/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detectors/staircase/enabled", FALSE) );
     SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
     SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );
   }

   virtual void TearDown() {
     SCIP_CALL_ABORT( SCIPfree(&scip) );
   }

   SCIP_RETCODE createVar(const char * str) {
      SCIP_VAR* var;
      SCIP_Bool success;
      SCIP_CALL( SCIPparseVar(scip, &var, str, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL, &success) );
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
   DEC_DECOMP* decomp;
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
   SCIP_CALL_EXPECT( DECcreateDecompFromMasterconss(scip, &decomp, &(mastercons), 1) );
   SCIP_CALL_EXPECT( SCIPconshdlrDecompAddDecdecomp(scip, decomp) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(2, GCGrelaxGetNPricingprobs(scip) );
   ASSERT_EQ(2, GCGrelaxGetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(0, GCGrelaxGetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(FALSE, GCGrelaxIsPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGrelaxIsPricingprobRelevant(scip, 0));
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
   SCIP_CALL_EXPECT( DECdetectStructure(scip, &result) );
   ASSERT_EQ(SCIP_SUCCESS, result);

   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_EQ(2, GCGrelaxGetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGrelaxGetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGrelaxGetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGrelaxIsPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGrelaxIsPricingprobRelevant(scip, 0));
}

TEST_F(GcgAggregationTest, WrongTypeTest) {
   DEC_DECOMP* decomp;
   SCIP_CONS* mastercons;
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[continuous] <x4>: obj=1.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 2<x1>[I] +2<x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x3>[I] +2<x4>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] +<x2>[C] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons = SCIPfindCons(scip, "c3");
   SCIP_CALL_EXPECT( DECcreateDecompFromMasterconss(scip, &decomp, &(mastercons), 1) );
   SCIP_CALL_EXPECT( SCIPconshdlrDecompAddDecdecomp(scip, decomp) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );
   ASSERT_EQ(2, GCGrelaxGetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGrelaxGetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGrelaxGetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGrelaxIsPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGrelaxIsPricingprobRelevant(scip, 0));
}

TEST_F(GcgAggregationTest, WrongBoundTest) {
   DEC_DECOMP* decomp;
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
   SCIP_CALL_EXPECT( DECcreateDecompFromMasterconss(scip, &decomp, &(mastercons), 1) );
   SCIP_CALL_EXPECT( SCIPconshdlrDecompAddDecdecomp(scip, decomp) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(2, GCGrelaxGetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGrelaxGetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGrelaxGetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGrelaxIsPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGrelaxIsPricingprobRelevant(scip, 0));
}

TEST_F(GcgAggregationTest, WrongCoeffSubproblemTest) {
   DEC_DECOMP* decomp;
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
   SCIP_CALL_EXPECT( DECcreateDecompFromMasterconss(scip, &decomp, &(mastercons), 1) );
   SCIP_CALL_EXPECT( SCIPconshdlrDecompAddDecdecomp(scip, decomp) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(2, GCGrelaxGetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGrelaxGetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGrelaxGetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGrelaxIsPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGrelaxIsPricingprobRelevant(scip, 0));
}

TEST_F(GcgAggregationTest, WrongCoeffMasterTest) {
   DEC_DECOMP* decomp;
   SCIP_CONS* mastercons;
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,4]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,4]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: 2<x1>[I] +2<x2>[I] >= 3") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: 2<x3>[I] +2<x4>[I] >= 3") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: 3<x1>[I] +2<x3>[I] <= 4") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   mastercons = SCIPfindCons(scip, "c3");
   SCIP_CALL_EXPECT( DECcreateDecompFromMasterconss(scip, &decomp, &(mastercons), 1) );
   SCIP_CALL_EXPECT( SCIPconshdlrDecompAddDecdecomp(scip, decomp) );
   SCIP_CALL_EXPECT( SCIPsolve(scip) );

   ASSERT_EQ(2, GCGrelaxGetNPricingprobs(scip) );
   ASSERT_EQ(1, GCGrelaxGetNIdenticalBlocks(scip, 0));
   ASSERT_EQ(1, GCGrelaxGetNIdenticalBlocks(scip, 1));
   ASSERT_EQ(TRUE, GCGrelaxIsPricingprobRelevant(scip, 1));
   ASSERT_EQ(TRUE, GCGrelaxIsPricingprobRelevant(scip, 0));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
