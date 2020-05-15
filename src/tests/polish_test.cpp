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

/**@file   decomp_test.cpp
 * @brief  Description
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "test.h"
#include "pub_decomp.h"
#include "struct_decomp.h"

class GcgPolishDecompTest : public ::testing::Test {

   virtual void SetUp() {
      decomp = NULL;
      newdecomp = NULL;

      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
      SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );

   }

   virtual void TearDown() {
      if( decomp != NULL)
      {
         SCIP_CALL_ABORT( DECdecompFree(scip, &decomp) );
      }
      if( newdecomp != NULL)
      {
         SCIP_CALL_ABORT( DECdecompFree(scip, &newdecomp) );
      }
      SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
protected:
   DEC_DECOMP* decomp;
   DEC_DECOMP* newdecomp;

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

SCIP* GcgPolishDecompTest::scip = NULL;

TEST_F(GcgPolishDecompTest, DetermineInPricing) {
   int block;

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

TEST_F(GcgPolishDecompTest, DetermineInMaster) {
   int block;

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] + <x2>[I] == 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: <x1>[I] + <x2>[I] + <x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CONS* masterconss[2] = {SCIPfindCons(scip, "c3"), SCIPfindCons(scip, "c4")};
   SCIP_CALL_EXPECT( DECcreateDecompFromMasterconss(scip, &decomp, masterconss, 2) );

   SCIP_CALL_EXPECT( DECdetermineConsBlock(scip, decomp, SCIPfindCons(scip, "c3"), &block));
   ASSERT_EQ(2, block);
   SCIP_CALL_EXPECT( DECdetermineConsBlock(scip, decomp, SCIPfindCons(scip, "c4"), &block));
   ASSERT_EQ(2, block);
}


TEST_F(GcgPolishDecompTest, DetermineLinkingvarOnly) {
   SCIP_HASHMAP* constoblock;
   int block;

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
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECdetermineConsBlock(scip, decomp, SCIPfindCons(scip, "c3"), &block));
   ASSERT_EQ(2, block);
}

TEST_F(GcgPolishDecompTest, DetermineNewPricingProblem) {
   SCIP_HASHMAP* constoblock;
   int block;

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
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECdetermineConsBlock(scip, decomp, SCIPfindCons(scip, "c3"), &block));
   ASSERT_EQ(-1, block);
}

TEST_F(GcgPolishDecompTest, TransferMasterconssToPricing) {
   SCIP_HASHMAP* constoblock;
   int transferred;
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
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToExistingPricing(scip, decomp, &transferred) );
   ASSERT_EQ(1, transferred);
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
}

TEST_F(GcgPolishDecompTest, DontTransferMasterconssToPricing) {
   SCIP_HASHMAP* constoblock;
   int transferred;

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
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToExistingPricing(scip, decomp, &transferred) );
   ASSERT_EQ(0, transferred);
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
}

TEST_F(GcgPolishDecompTest, DontTransferLinkingVarsOnly) {
   SCIP_HASHMAP* constoblock;
   int transferred;

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
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToExistingPricing(scip, decomp, &transferred) );
   ASSERT_EQ(0, transferred);
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
}

TEST_F(GcgPolishDecompTest, TransferLinkingVarToPricing) {
   SCIP_HASHMAP* constoblock;
   int transferred;

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
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToExistingPricing(scip, decomp, &transferred) );
   ASSERT_EQ(1, transferred);
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
}

TEST_F(GcgPolishDecompTest, TransferNewVarToPricing) {
   SCIP_HASHMAP* constoblock;
   int transferred;

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,1]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] + <x4>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 3) );

   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToExistingPricing(scip, decomp, &transferred) );
   ASSERT_EQ(1, transferred);
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
}

TEST_F(GcgPolishDecompTest, TransferNewVarToPricingWithLinking) {
   SCIP_HASHMAP* constoblock;
   int transferred;

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,1]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I] + <x3>[I] + <x4>[I] == 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: <x1>[I] + <x2>[I]  == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 4) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c4"), (void*) 3) );

   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToExistingPricing(scip, decomp, &transferred) );
   ASSERT_EQ(1, transferred);
   ASSERT_EQ(2, DECdecompGetNBlocks(decomp));
   ASSERT_EQ(1, DECdecompGetNLinkingvars(decomp));
   ASSERT_EQ(SCIPfindVar(scip, "t_x3"), DECdecompGetLinkingvars(decomp)[0]);
}

TEST_F(GcgPolishDecompTest, CreateNewPricingProblem) {
   SCIP_HASHMAP* constoblock;
   int transferred;

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
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToNewPricing(scip, decomp, &newdecomp, &transferred) );
   ASSERT_EQ(1, transferred);
   ASSERT_EQ(3, DECdecompGetNBlocks(newdecomp));
}

TEST_F(GcgPolishDecompTest, CreateNewPricingProblemWithLinking) {
   SCIP_HASHMAP* constoblock;
   int transferred;

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,8]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x4>[I] + <x3>[I]== 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 3) );

   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECtryAssignMasterconssToNewPricing(scip, decomp, &newdecomp, &transferred) );
   ASSERT_EQ(1, transferred);
   ASSERT_EQ(3, DECdecompGetNBlocks(newdecomp));
}

TEST_F(GcgPolishDecompTest, PolishDecompTransferAll) {
   SCIP_HASHMAP* constoblock;

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,8]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x5>: obj=2.0, original bounds=[0,8]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c1a>: <x1>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] + <x3>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2a>: <x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x4>[I] + <x3>[I]== 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: <x5>[I] + <x3>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 6) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1a"), (void*) 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2a"), (void*) 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c4"), (void*) 3) );

   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECcreatePolishedDecomp(scip, decomp, &newdecomp) );
   ASSERT_EQ(4, DECdecompGetNBlocks(newdecomp));
   ASSERT_EQ(0, DECdecompGetNLinkingconss(newdecomp));
}

TEST_F(GcgPolishDecompTest, PolishDecompOnlyNew) {
   SCIP_HASHMAP* constoblock;

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x4>: obj=2.0, original bounds=[0,8]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x5>: obj=2.0, original bounds=[0,8]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x3>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: <x4>[I]== 1") );
   SCIP_CALL_EXPECT( createCons("[linear] <c5>: <x5>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 5) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c4"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c5"), (void*) 1) );

   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 0, FALSE) );

   SCIP_CALL_EXPECT( DECcreatePolishedDecomp(scip, decomp, &newdecomp) );
   ASSERT_EQ(5, DECdecompGetNBlocks(newdecomp));
   ASSERT_EQ(0, DECdecompGetNLinkingconss(newdecomp));
}

TEST_F(GcgPolishDecompTest, PolishDecompNothingNew) {
   SCIP_HASHMAP* constoblock;

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x3>: obj=2.0, original bounds=[0,1]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x2>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x2>[I] <= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c3>: <x1>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c4>: <x1>[I]+ <x3>[I] + <x2>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 5) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 3) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c3"), (void*) 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c4"), (void*) 3) );

   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );

   SCIP_CALL_EXPECT( DECcreatePolishedDecomp(scip, decomp, &newdecomp) );
   ASSERT_EQ((DEC_DECOMP*) NULL, newdecomp);
}

TEST_F(GcgPolishDecompTest, DontPolishOneBlock) {
   SCIP_HASHMAP* constoblock;

   SCIP_CALL_EXPECT( createVar("[integer] <x1>: obj=2.0, original bounds=[0,1]") );
   SCIP_CALL_EXPECT( createVar("[integer] <x2>: obj=2.0, original bounds=[0,3]") );

   SCIP_CALL_EXPECT( createCons("[linear] <c1>: <x1>[I] + <x2>[I]<= 5") );
   SCIP_CALL_EXPECT( createCons("[linear] <c2>: <x1>[I] == 1") );

   SCIP_CALL_EXPECT( SCIPtransformProb(scip) );
   SCIP_CALL_EXPECT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), 2) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c1"), (void*) 1) );
   SCIP_CALL_EXPECT( SCIPhashmapInsert(constoblock, SCIPfindCons(scip, "c2"), (void*) 2) );

   SCIP_CALL_EXPECT( DECdecompCreate(scip, &decomp) );
   SCIP_CALL_EXPECT( DECfilloutDecompFromConstoblock(scip, decomp, constoblock, 1, FALSE) );

   SCIP_CALL_EXPECT( DECcreatePolishedDecomp(scip, decomp, &newdecomp) );
   ASSERT_EQ((DEC_DECOMP*) NULL, newdecomp);
}
