/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
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
      SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
   }

   virtual void TearDown() {
      if( decomp != NULL)
      {
         SCIP_CALL_ABORT( GCGdecompFree(scip, &decomp) );
      }

      SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
protected:
   GCG_DECOMP* decomp;

public:
   static SCIP* scip;
};

SCIP* GcgDecompTest::scip = NULL;

TEST_F(GcgDecompTest, CreateAndFreeTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );

   ASSERT_EQ((uint) FALSE, decomp->presolved);
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
   ASSERT_EQ(GCG_DECTYPE_UNKNOWN , decomp->type);
   ASSERT_EQ(NULL, decomp->detector);

   SCIP_CALL_EXPECT( GCGdecompFree(scip, &decomp) );
   ASSERT_EQ(NULL, decomp);
}

TEST_F(GcgDecompTest, GetDetectorTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, GCGdecompGetDetector(decomp));
   decomp->detector = (DEC_DETECTOR*) 0xDEADBEEF;
   ASSERT_EQ((DEC_DETECTOR*) 0xDEADBEEF, GCGdecompGetDetector(decomp));
}

TEST_F(GcgDecompTest, SetDetectorTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, decomp->detector);
   GCGdecompSetDetector(decomp, (DEC_DETECTOR*) 0xDEADBEEF);
   ASSERT_EQ(decomp->detector, (DEC_DETECTOR*) 0xDEADBEEF);
}

TEST_F(GcgDecompTest, GetConsindexTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, GCGdecompGetConsindex(decomp));
   decomp->consindex = (SCIP_HASHMAP*) 0xDEADBEEF;
   ASSERT_EQ((SCIP_HASHMAP*) 0xDEADBEEF, GCGdecompGetConsindex(decomp));
   decomp->consindex = NULL;
}

TEST_F(GcgDecompTest, SetConsindexTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, decomp->consindex);
   GCGdecompSetConsindex(decomp, (SCIP_HASHMAP*) 0xDEADBEEF);
   ASSERT_EQ(decomp->consindex, (SCIP_HASHMAP*) 0xDEADBEEF);
   decomp->consindex = NULL;
}

TEST_F(GcgDecompTest, GetVarindexTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, GCGdecompGetVarindex(decomp));
   decomp->varindex = (SCIP_HASHMAP*) 0xDEADBEEF;
   ASSERT_EQ((SCIP_HASHMAP*) 0xDEADBEEF, GCGdecompGetVarindex(decomp));
   decomp->varindex = NULL;
}

TEST_F(GcgDecompTest, SetVarindexTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );
   ASSERT_EQ(NULL, decomp->varindex);
   GCGdecompSetVarindex(decomp, (SCIP_HASHMAP*) 0xDEADBEEF);
   ASSERT_EQ(decomp->varindex, (SCIP_HASHMAP*) 0xDEADBEEF);
   decomp->varindex = NULL;
   SCIP_CALL_EXPECT(GCGdecompFree(scip, &decomp));
}

TEST_F(GcgDecompTest, SetTypeDiagonalTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );
   ASSERT_EQ(GCG_DECTYPE_UNKNOWN, decomp->type);
   ASSERT_EQ(SCIP_OKAY, GCGdecompSetType(decomp, GCG_DECTYPE_DIAGONAL));
   decomp->nlinkingconss = 1;
   ASSERT_EQ(SCIP_INVALIDDATA, GCGdecompSetType(decomp, GCG_DECTYPE_DIAGONAL));
   decomp->nlinkingconss = 0;
   decomp->linkingconss = (SCIP_CONS**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_INVALIDDATA, GCGdecompSetType(decomp, GCG_DECTYPE_DIAGONAL));
   decomp->linkingconss = NULL;
   decomp->nlinkingvars = 1;
   ASSERT_EQ(SCIP_INVALIDDATA, GCGdecompSetType(decomp, GCG_DECTYPE_DIAGONAL));
   decomp->nlinkingvars = 0;
   decomp->linkingvars = (SCIP_VAR**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_INVALIDDATA, GCGdecompSetType(decomp, GCG_DECTYPE_DIAGONAL));
   decomp->linkingvars = NULL;
}

TEST_F(GcgDecompTest, SetTypeUnknownTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );
   ASSERT_EQ(GCG_DECTYPE_UNKNOWN, decomp->type);
   ASSERT_EQ(SCIP_INVALIDDATA, GCGdecompSetType(decomp, GCG_DECTYPE_UNKNOWN));
}

TEST_F(GcgDecompTest, SetTypeArrowheadTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );
   ASSERT_EQ(GCG_DECTYPE_UNKNOWN, decomp->type);
   ASSERT_EQ(SCIP_OKAY, GCGdecompSetType(decomp, GCG_DECTYPE_ARROWHEAD));
   decomp->nlinkingconss = 1;
   ASSERT_EQ(SCIP_OKAY, GCGdecompSetType(decomp, GCG_DECTYPE_ARROWHEAD));
   decomp->nlinkingconss = 0;
   decomp->linkingconss = (SCIP_CONS**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_OKAY, GCGdecompSetType(decomp, GCG_DECTYPE_ARROWHEAD));
   decomp->linkingconss = NULL;
   decomp->nlinkingvars = 1;
   ASSERT_EQ(SCIP_OKAY, GCGdecompSetType(decomp, GCG_DECTYPE_ARROWHEAD));
   decomp->nlinkingvars = 0;
   decomp->linkingvars = (SCIP_VAR**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_OKAY, GCGdecompSetType(decomp, GCG_DECTYPE_ARROWHEAD));
   decomp->linkingvars = NULL;
}

TEST_F(GcgDecompTest, SetTypeBorderedTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );
   ASSERT_EQ(GCG_DECTYPE_UNKNOWN, decomp->type);
   ASSERT_EQ(SCIP_OKAY, GCGdecompSetType(decomp, GCG_DECTYPE_BORDERED));
   decomp->nlinkingconss = 1;
   ASSERT_EQ(SCIP_OKAY, GCGdecompSetType(decomp, GCG_DECTYPE_BORDERED));
   decomp->nlinkingconss = 0;
   decomp->linkingconss = (SCIP_CONS**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_OKAY, GCGdecompSetType(decomp, GCG_DECTYPE_BORDERED));
   decomp->linkingconss = NULL;
   decomp->nlinkingvars = 1;
   ASSERT_EQ(SCIP_INVALIDDATA, GCGdecompSetType(decomp, GCG_DECTYPE_BORDERED));
   decomp->nlinkingvars = 0;
   decomp->linkingvars = (SCIP_VAR**) 0xDEADBEEF;
   ASSERT_EQ(SCIP_INVALIDDATA, GCGdecompSetType(decomp, GCG_DECTYPE_BORDERED));
   decomp->linkingvars = NULL;
}

TEST_F(GcgDecompTest, GetPresolvedTest) {
   SCIP_CALL_EXPECT( GCGdecompCreate(scip, &decomp) );
   ASSERT_EQ((uint)FALSE, decomp->presolved);
   ASSERT_EQ((uint)FALSE, GCGdecompGetPresolved(decomp));
   decomp->presolved = TRUE;
   ASSERT_EQ((uint)TRUE, GCGdecompGetPresolved(decomp));
}
