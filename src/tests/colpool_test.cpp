/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2014 Operations Research, RWTH Aachen University       */
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

/**@file   gcgcol_test.cpp
 * @brief  Description
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "type_gcgcol.h"
#include "pub_gcgcol.h"
#include "test.h"
#include "scip/struct_var.h"
#include "gcg.h"
#include "gcgpqueue.h"
#include "pub_gcgpqueue.h"
#include "class_colpool.h"


using gcg::Colpool;

class ColpoolTest : public ::testing::Test {

   virtual void SetUp() {
      decomp = NULL;
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
   }

   virtual void TearDown() {
      SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
protected:
   DEC_DECOMP* decomp;

public:
   static SCIP* scip;
};

SCIP* ColpoolTest::scip = NULL;

TEST_F(ColpoolTest, CreateEmptyColpoolTest) {

   Colpool* colpool;

   colpool = new Colpool(scip, 5, 10, 10);

   delete colpool;
}

TEST_F(ColpoolTest, CreateColpoolTest) {

   Colpool* colpool;

   GCG_COL** gcgcols;
   GCG_COL* gcgcol;
   SCIP_Real redcosts[4] = {1.0, 2.0, 0.0, -1.0};
   int probs[4] = {0, 1, 2, 3};

   SCIP_Bool success;

   int i;

   colpool = new Colpool(scip, 5, 10, 10);

   SCIP_CALL_EXPECT( SCIPallocMemoryArray(scip, &gcgcols, 4) );

   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL_EXPECT( SCIPallocMemory(scip, &(gcgcols[i])) );

      gcgcols[i]->redcost = redcosts[i];

      gcgcols[i]->probnr = probs[i];

      colpool->addCol(gcgcols[i], &success);

      ASSERT_EQ(success, TRUE);
   }

   colpool->getBestCol(&gcgcol);
   ASSERT_EQ(gcgcol, gcgcols[3]);
   SCIPfreeMemory(scip, &gcgcol);

   colpool->getBestCol(&gcgcol);
   ASSERT_EQ(gcgcol, gcgcols[2]);
   SCIPfreeMemory(scip, &gcgcol);

   colpool->getBestCol(&gcgcol);
   ASSERT_EQ(gcgcol, gcgcols[0]);
   SCIPfreeMemory(scip, &gcgcol);

   colpool->getBestCol(&gcgcol);
   ASSERT_EQ(gcgcol, gcgcols[1]);
   SCIPfreeMemory(scip, &gcgcol);

   SCIPfreeMemoryArray(scip, &gcgcols);

   delete colpool;
}

TEST_F(ColpoolTest, DeleteOldTest) {

   Colpool* colpool;

   GCG_COL** gcgcols;
   GCG_COL* gcgcol;
   SCIP_Real redcosts[4] = {1.0, 2.0, 0.0, -1.0};
   int probs[4] = {0, 1, 2, 3};
   int ages[4] = {4, 9, 2, 7};

   SCIP_Bool success;

   int i;

   colpool = new Colpool(scip, 5, 10, 10);

   SCIP_CALL_EXPECT( SCIPallocMemoryArray(scip, &gcgcols, 4) );

   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL_EXPECT( GCGcreateGcgCol(scip, &(gcgcols[i]), probs[i], NULL, NULL, 0, FALSE, redcosts[i] ) );

      gcgcols[i]->redcost = redcosts[i];

      gcgcols[i]->probnr = probs[i];

      gcgcols[i]->age = ages[i];

      colpool->addCol(gcgcols[i], &success);

      ASSERT_EQ(success, TRUE);
   }

   colpool->deleteOldColumns();

   ASSERT_EQ(colpool->getNCols(), 2);

   colpool->getBestCol(&gcgcol);
   ASSERT_EQ(gcgcol, gcgcols[2]);
   GCGfreeGcgCol(&gcgcol);

   colpool->getBestCol(&gcgcol);
   ASSERT_EQ(gcgcol, gcgcols[0]);
   GCGfreeGcgCol(&gcgcol);

   SCIPfreeMemoryArray(scip, &gcgcols);

   delete colpool;
}
