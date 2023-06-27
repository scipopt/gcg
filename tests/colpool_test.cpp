/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       */
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
#include "pub_colpool.h"

class ColpoolTest : public ::testing::Test {

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

SCIP* ColpoolTest::scip = NULL;

TEST_F(ColpoolTest, CreateEmptyColpoolTest) {

   GCG_COLPOOL* colpool;

   SCIP_CALL_EXPECT( GCGcolpoolCreate(scip, &colpool, 5) );

   SCIP_CALL_ABORT( GCGcolpoolFree(scip, &colpool) );
}
