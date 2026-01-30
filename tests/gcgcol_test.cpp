/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   gcgcol_test.cpp
 * @brief  Description
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "test.h"
#include "pub_gcgcol.h"
#include "scip/struct_var.h"
#include "gcg.h"
#include "type_gcgcol.h"
#include "struct_gcgcol.h"


class GcgColTest : public ::testing::Test {

   virtual void SetUp() {
      decomp = NULL;
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
   }

   virtual void TearDown() {
      SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
protected:
   GCG_DECOMP* decomp;

public:
   static SCIP* scip;
};

SCIP* GcgColTest::scip = NULL;

TEST_F(GcgColTest, CreateEmptyColTest) {

   GCG_COL* gcgcol;

   SCIP_CALL_EXPECT( GCGcreateGcgCol(scip, &gcgcol, 0, NULL, NULL, 0, FALSE, SCIPinfinity(scip)) );

   GCGfreeGcgCol(&gcgcol);
}

TEST_F(GcgColTest, CreateColTest) {

   GCG_COL* gcgcol;
   SCIP_VAR** vars;
   SCIP_Real vals[4] = {1.0, 2.0, 0.0, -1.0};

   int i;

   SCIP_CALL_EXPECT( SCIPallocMemoryArray(scip, &vars, 4) );

   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL_EXPECT( SCIPallocMemory(scip, &(vars[i])) );

      vars[i]->index = 4 - i;

      vars[i]->varstatus = SCIP_VARSTATUS_ORIGINAL;
   }

   SCIP_CALL_EXPECT( GCGcreateGcgCol(scip, &gcgcol, 0, vars, vals, 4, FALSE, SCIPinfinity(scip)) );

   ASSERT_EQ(gcgcol->nvars, 3);
   ASSERT_EQ(gcgcol->vars[0], vars[3]);
   ASSERT_EQ(gcgcol->vars[1], vars[1]);
   ASSERT_EQ(gcgcol->vars[2], vars[0]);

   ASSERT_EQ(gcgcol->vals[0], vals[3]);
   ASSERT_EQ(gcgcol->vals[1], vals[1]);
   ASSERT_EQ(gcgcol->vals[2], vals[0]);

   ASSERT_EQ(gcgcol->probnr, 0);
   ASSERT_EQ(gcgcol->isray, FALSE);
   ASSERT_EQ(gcgcol->redcost, SCIPinfinity(scip));

   GCGfreeGcgCol(&gcgcol);

   for( i = 0; i < 4; ++i )
   {
      SCIPfreeMemory(scip, &(vars[i]));
   }

   SCIPfreeMemoryArray(scip, &vars);
}

TEST_F(GcgColTest, CreateColFromSolTest) {

   SCIPinfoMessage(scip, NULL, "Cannot test GCGcreateGcgColFromSol(), because it uses GCG methods\n");
}

TEST_F(GcgColTest, EqColsColIsEqTest) {

   GCG_COL* gcgcol1;
   GCG_COL* gcgcol2;
   SCIP_VAR** vars;
   SCIP_Real vals[4] = {1.0, 2.0, 0.0, -1.0};

   int i;

   SCIP_CALL_EXPECT( SCIPallocMemoryArray(scip, &vars, 4) );

   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL_EXPECT( SCIPallocMemory(scip, &(vars[i])) );

      vars[i]->index = 4 - i;

      vars[i]->varstatus = SCIP_VARSTATUS_ORIGINAL;
   }

   SCIP_CALL_EXPECT( GCGcreateGcgCol(scip, &gcgcol1, 0, vars, vals, 4, FALSE, SCIPinfinity(scip)) );
   SCIP_CALL_EXPECT( GCGcreateGcgCol(scip, &gcgcol2, 0, vars, vals, 4, FALSE, 1.0) );

   ASSERT_EQ(GCGcolIsEq(gcgcol1, gcgcol2), TRUE);

   GCGfreeGcgCol(&gcgcol2);
   GCGfreeGcgCol(&gcgcol1);

   for( i = 0; i < 4; ++i )
   {
      SCIPfreeMemory(scip, &(vars[i]));
   }

   SCIPfreeMemoryArray(scip, &vars);
}

TEST_F(GcgColTest, NeqColsColIsEqTest) {

   GCG_COL* gcgcol1;
   GCG_COL* gcgcol2;
   SCIP_VAR** vars;
   SCIP_Real vals[4] = {1.0, 2.0, 0.0, -1.0};

   int i;

   SCIP_CALL_EXPECT( SCIPallocMemoryArray(scip, &vars, 4) );

   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL_EXPECT( SCIPallocMemory(scip, &(vars[i])) );

      vars[i]->index = 4 - i;

      vars[i]->varstatus = SCIP_VARSTATUS_ORIGINAL;
   }

   SCIP_CALL_EXPECT( GCGcreateGcgCol(scip, &gcgcol1, 0, vars, vals, 4, FALSE, SCIPinfinity(scip)) );

   vals[2] = 3.0;

   SCIP_CALL_EXPECT( GCGcreateGcgCol(scip, &gcgcol2, 0, vars, vals, 4, FALSE, 1.0) );

   ASSERT_EQ(GCGcolIsEq(gcgcol1, gcgcol2), FALSE);

   GCGfreeGcgCol(&gcgcol2);
   GCGfreeGcgCol(&gcgcol1);

   for( i = 0; i < 4; ++i )
   {
      SCIPfreeMemory(scip, &(vars[i]));
   }

   SCIPfreeMemoryArray(scip, &vars);
}

TEST_F(GcgColTest, GetSolValTest) {

   GCG_COL* gcgcol;
   SCIP_VAR** vars;
   SCIP_Real vals[4] = {1.0, 2.0, 0.0, -1.0};

   int i;

   SCIP_CALL_EXPECT( SCIPallocMemoryArray(scip, &vars, 4) );

   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL_EXPECT( SCIPallocMemory(scip, &(vars[i])) );

      vars[i]->index = 4 - i;

      vars[i]->varstatus = SCIP_VARSTATUS_ORIGINAL;
   }

   SCIP_CALL_EXPECT( GCGcreateGcgCol(scip, &gcgcol, 0, vars, vals, 4, FALSE, SCIPinfinity(scip)) );

   ASSERT_EQ(GCGcolGetSolVal(scip, gcgcol, vars[0]), 1.0);
   ASSERT_EQ(GCGcolGetSolVal(scip, gcgcol, vars[1]), 2.0);
   ASSERT_EQ(GCGcolGetSolVal(scip, gcgcol, vars[2]), 0.0);
   ASSERT_EQ(GCGcolGetSolVal(scip, gcgcol, vars[3]), -1.0);

   GCGfreeGcgCol(&gcgcol);

   for( i = 0; i < 4; ++i )
   {
      SCIPfreeMemory(scip, &(vars[i]));
   }

   SCIPfreeMemoryArray(scip, &vars);
}
