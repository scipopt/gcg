/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
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

/**@file   decomp_test.cpp
 * @brief  Description
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "test.h"
#include "pub_decomp.h"
#include "struct_decomp.h"

class GcgDecompTest : public ::testing::Test {

   virtual void SetUp() {
      decomp = NULL;
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
   }

   virtual void TearDown() {
      if( decomp != NULL)
      {
         SCIP_CALL_ABORT( DECdecompFree(scip, &decomp) );
      }

      SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
protected:
   DEC_DECOMP* decomp;

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

public:
   static SCIP* scip;
};

SCIP* GcgDecompTest::scip = NULL;

TEST_F(GcgDecompTest, CreateAndFreeTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );

   ASSERT_EQ(FALSE, decomp->presolved);
   ASSERT_EQ(0 , decomp->nblocks);
   ASSERT_EQ(NULL , decomp->subscipvars);
   ASSERT_EQ(NULL , decomp->nsubscipvars);
   ASSERT_EQ(NULL , decomp->subscipconss);
   ASSERT_EQ(NULL , decomp->nsubscipconss);
   ASSERT_EQ(NULL , decomp->linkingconss);
   ASSERT_EQ(0 , decomp->nlinkingconss);
   ASSERT_EQ(NULL , decomp->linkingvars);
   ASSERT_EQ(0 , decomp->nlinkingvars);
   ASSERT_EQ(NULL , decomp->stairlinkingvars);
   ASSERT_EQ(NULL , decomp->nstairlinkingvars);
   ASSERT_EQ(NULL , decomp->vartoblock);
   ASSERT_EQ(NULL , decomp->constoblock);
   ASSERT_EQ(NULL , decomp->varindex);
   ASSERT_EQ(NULL, decomp->consindex);
   ASSERT_EQ(DEC_DECTYPE_UNKNOWN , decomp->type);
   ASSERT_EQ(NULL, decomp->detector);

   SCIP_CALL_EXPECT( DECdecompFree(scip, &decomp) );
   ASSERT_EQ(NULL, decomp);
}

TEST_F(GcgDecompTest, GetDetectorTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, DECdecompGetDetector(decomp));
   decomp->detector = (DEC_DETECTOR*) 0xDEADBEEF;
   ASSERT_EQ((DEC_DETECTOR*) 0xDEADBEEF, DECdecompGetDetector(decomp));
}

TEST_F(GcgDecompTest, SetDetectorTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, decomp->detector);
   DECdecompSetDetector(decomp, (DEC_DETECTOR*) 0xDEADBEEF);
   ASSERT_EQ(decomp->detector, (DEC_DETECTOR*) 0xDEADBEEF);
}

TEST_F(GcgDecompTest, GetConsindexTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, DECdecompGetConsindex(decomp));
   decomp->consindex = (SCIP_HASHMAP*) 0xDEADBEEF;
   ASSERT_EQ((SCIP_HASHMAP*) 0xDEADBEEF, DECdecompGetConsindex(decomp));
   decomp->consindex = NULL;
}

TEST_F(GcgDecompTest, SetConsindexTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, decomp->consindex);
   DECdecompSetConsindex(decomp, (SCIP_HASHMAP*) 0xDEADBEEF);
   ASSERT_EQ(decomp->consindex, (SCIP_HASHMAP*) 0xDEADBEEF);
   decomp->consindex = NULL;
}

TEST_F(GcgDecompTest, GetVarindexTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, DECdecompGetVarindex(decomp));
   decomp->varindex = (SCIP_HASHMAP*) 0xDEADBEEF;
   ASSERT_EQ((SCIP_HASHMAP*) 0xDEADBEEF, DECdecompGetVarindex(decomp));
   decomp->varindex = NULL;
}

TEST_F(GcgDecompTest, SetVarindexTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, decomp->varindex);
   DECdecompSetVarindex(decomp, (SCIP_HASHMAP*) 0xDEADBEEF);
   ASSERT_EQ(decomp->varindex, (SCIP_HASHMAP*) 0xDEADBEEF);
   decomp->varindex = NULL;
   SCIP_CALL_EXPECT(DECdecompFree(scip, &decomp));
}

TEST_F(GcgDecompTest, SetTypeDiagonalTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   ASSERT_EQ(DEC_DECTYPE_UNKNOWN, decomp->type);
   ASSERT_EQ(SCIP_OKAY, DECdecompSetType(decomp, DEC_DECTYPE_DIAGONAL));
   decomp->nlinkingconss = 1;
   ASSERT_EQ(SCIP_INVALIDDATA, DECdecompSetType(decomp, DEC_DECTYPE_DIAGONAL));
   decomp->nlinkingconss = 0;
   decomp->linkingconss = (SCIP_CONS**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_INVALIDDATA, DECdecompSetType(decomp, DEC_DECTYPE_DIAGONAL));
   decomp->linkingconss = NULL;
   decomp->nlinkingvars = 1;
   ASSERT_EQ(SCIP_INVALIDDATA, DECdecompSetType(decomp, DEC_DECTYPE_DIAGONAL));
   decomp->nlinkingvars = 0;
   decomp->linkingvars = (SCIP_VAR**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_INVALIDDATA, DECdecompSetType(decomp, DEC_DECTYPE_DIAGONAL));
   decomp->linkingvars = NULL;
}

TEST_F(GcgDecompTest, SetTypeUnknownTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   ASSERT_EQ(DEC_DECTYPE_UNKNOWN, decomp->type);
   ASSERT_EQ(SCIP_INVALIDDATA, DECdecompSetType(decomp, DEC_DECTYPE_UNKNOWN));
}

TEST_F(GcgDecompTest, SetTypeArrowheadTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   ASSERT_EQ(DEC_DECTYPE_UNKNOWN, decomp->type);
   ASSERT_EQ(SCIP_OKAY, DECdecompSetType(decomp, DEC_DECTYPE_ARROWHEAD));
   decomp->nlinkingconss = 1;
   ASSERT_EQ(SCIP_OKAY, DECdecompSetType(decomp, DEC_DECTYPE_ARROWHEAD));
   decomp->nlinkingconss = 0;
   decomp->linkingconss = (SCIP_CONS**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_OKAY, DECdecompSetType(decomp, DEC_DECTYPE_ARROWHEAD));
   decomp->linkingconss = NULL;
   decomp->nlinkingvars = 1;
   ASSERT_EQ(SCIP_OKAY, DECdecompSetType(decomp, DEC_DECTYPE_ARROWHEAD));
   decomp->nlinkingvars = 0;
   decomp->linkingvars = (SCIP_VAR**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_OKAY, DECdecompSetType(decomp, DEC_DECTYPE_ARROWHEAD));
   decomp->linkingvars = NULL;
}

TEST_F(GcgDecompTest, SetTypeBorderedTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   ASSERT_EQ(DEC_DECTYPE_UNKNOWN, decomp->type);
   ASSERT_EQ(SCIP_OKAY, DECdecompSetType(decomp, DEC_DECTYPE_BORDERED));
   decomp->nlinkingconss = 1;
   ASSERT_EQ(SCIP_OKAY, DECdecompSetType(decomp, DEC_DECTYPE_BORDERED));
   decomp->nlinkingconss = 0;
   decomp->linkingconss = (SCIP_CONS**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_OKAY, DECdecompSetType(decomp, DEC_DECTYPE_BORDERED));
   decomp->linkingconss = NULL;
   decomp->nlinkingvars = 1;
   ASSERT_EQ(SCIP_INVALIDDATA, DECdecompSetType(decomp, DEC_DECTYPE_BORDERED));
   decomp->nlinkingvars = 0;
   decomp->linkingvars = (SCIP_VAR**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_INVALIDDATA, DECdecompSetType(decomp, DEC_DECTYPE_BORDERED));
   decomp->linkingvars = NULL;
}

TEST_F(GcgDecompTest, GetPresolvedTest) {
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   ASSERT_EQ(FALSE, decomp->presolved);
   ASSERT_EQ(FALSE, DECdecompGetPresolved(decomp));
   decomp->presolved = TRUE;
   ASSERT_EQ(TRUE, DECdecompGetPresolved(decomp));
}

TEST_F(GcgDecompTest, determineConsInPricing) {
   int block;
   SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
   SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] == 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: <x2>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( DECcreateDecompFromMasterconss(scip, &decomp, NULL, 0) );

   SCIP_CALL_EXPECT( DECdetermineConsBlock(scip, decomp, SCIPfindCons(scip, "c3"), &block));
   ASSERT_EQ(0, block);
   SCIP_CALL_EXPECT( DECdetermineConsBlock(scip, decomp, SCIPfindCons(scip, "c4"), &block));
   ASSERT_EQ(1, block);
}

TEST_F(GcgDecompTest, determineConsInMaster) {
   int block;
   SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
   SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );
   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] + <x2>[I] == 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: <x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CONS* masterconss[2] = {SCIPfindCons(scip, "c3"), SCIPfindCons(scip, "c4")};
   SCIP_CALL_EXPECT( DECcreateDecompFromMasterconss(scip, &decomp, masterconss, 2) );

   SCIP_CALL_EXPECT( DECdetermineConsBlock(scip, decomp, SCIPfindCons(scip, "c3"), &block));
   ASSERT_EQ(2, block);
   SCIP_CALL_EXPECT( DECdetermineConsBlock(scip, decomp, SCIPfindCons(scip, "c4"), &block));
   ASSERT_EQ(2, block);
}


TEST_F(GcgDecompTest, determineConsLinkingvarOnly) {
   SCIP_HASHMAP* constoblock;
   int block;
   SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
   SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 3) );
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecdecompFromConstoblock(scip, decomp, constoblock, 2, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetConss(scip), SCIPgetNConss(scip), FALSE) );

   SCIP_CALL_EXPECT( DECdetermineConsBlock(scip, decomp, SCIPfindCons(scip, "c3"), &block));
   ASSERT_EQ(2, block);
}

TEST_F(GcgDecompTest, determineConsNewPricingProblem) {
   SCIP_HASHMAP* constoblock;
   int block;
   SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
   SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,1]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x4>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 3) );
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecdecompFromConstoblock(scip, decomp, constoblock, 2, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetConss(scip), SCIPgetNConss(scip), FALSE) );

   SCIP_CALL_EXPECT( DECdetermineConsBlock(scip, decomp, SCIPfindCons(scip, "c3"), &block));
   ASSERT_EQ(-1, block);
}

TEST_F(GcgDecompTest, polishMasterconssToPricing) {
   SCIP_HASHMAP* constoblock;
   int transferred;

   SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
   SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x2>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 3) );

   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecdecompFromConstoblock(scip, decomp, constoblock, 2, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetConss(scip), SCIPgetNConss(scip), FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToPricing(scip, decomp, &transferred) );
   ASSERT_EQ(1, transferred);
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
}

TEST_F(GcgDecompTest, polishNoMasterconssToPricing) {
   SCIP_HASHMAP* constoblock;
   int transferred;

   SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
   SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I]  <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] + <x2>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 3) );
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecdecompFromConstoblock(scip, decomp, constoblock, 2, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetConss(scip), SCIPgetNConss(scip), FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToPricing(scip, decomp, &transferred) );
   ASSERT_EQ(0, transferred);
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
}

TEST_F(GcgDecompTest, polishLinkingVarsNoTransfer) {
   SCIP_HASHMAP* constoblock;
   int transferred;

   SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
   SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 3) );
   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecdecompFromConstoblock(scip, decomp, constoblock, 2, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetConss(scip), SCIPgetNConss(scip), FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToPricing(scip, decomp, &transferred) );
   ASSERT_EQ(0, transferred);
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
}

TEST_F(GcgDecompTest, polishLinkingVarsTransfer) {
   SCIP_HASHMAP* constoblock;
   int transferred;

   SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
   SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] + <x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 3) );

   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecdecompFromConstoblock(scip, decomp, constoblock, 2, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetConss(scip), SCIPgetNConss(scip), FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToPricing(scip, decomp, &transferred) );
   ASSERT_EQ(1, transferred);
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
}

TEST_F(GcgDecompTest, polishNewPricingProblem) {
   SCIP_HASHMAP* constoblock;
   int transferred;

   SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
   SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,8]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x4>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 3) );

   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecdecompFromConstoblock(scip, decomp, constoblock, 2, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetConss(scip), SCIPgetNConss(scip), FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToPricing(scip, decomp, &transferred) );
   ASSERT_EQ(1, transferred);
   ASSERT_EQ(3, DECdecompGetNBlocks(decomp));

}
