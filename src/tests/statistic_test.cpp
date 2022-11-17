/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2022 Operations Research, RWTH Aachen University       */
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

/**@file   statistic_test.cpp
 * @brief  Description
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "test.h"
#include "pub_decomp.h"

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
      * min x1+x2+x3+x4+x5
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
      SCIP_CALL_ABORT( GCGdecompFreeCreate(scip, &decomp) );

      SCIP_CALL_ABORT( GCGfilloutDecompFromConstoblock(scip, decomp, constoblock, 2, FALSE) );
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
      SCIP_CALL_ABORT( GCGdecompFreeFree(scip, &decomp) );


      SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
};

SCIP* GcgDecStatisticTest::scip = NULL;

TEST_F(GcgDecStatisticTest, BlockTest) {
   ASSERT_EQ(2, GCGdecompFreeGetNBlocks(decomp));
}

TEST_F(GcgDecStatisticTest, SubscipSizeTest) {
   ASSERT_EQ(2, GCGdecompFreeGetNBlocks(decomp));
   ASSERT_EQ(2, GCGdecompGetNSubscipvars(decomp)[0]);
   ASSERT_EQ(2, GCGdecompGetNSubscipvars(decomp)[1]);
   ASSERT_EQ(1, GCGdecompFreeGetNLinkingvars(decomp));
   ASSERT_EQ(1, GCGdecompFreeGetNSubscipconss(decomp)[0]);
   ASSERT_EQ(1, GCGdecompFreeGetNSubscipconss(decomp)[0]);
   ASSERT_EQ(1, GCGdecompGetNLinkingconss(decomp));
}

TEST_F(GcgDecStatisticTest, DensityTest) {
   int i;
   SCIP_VAR* densvars[5];
   SCIP_CONS* densconss[5];
   int varsubprobdens[5];
   int varmasterdens[5];
   int conssubprobdens[3];
   int consmasterdens[3];

   SCIP_CALL_EXPECT(GCGgetDensityData(scip, decomp, densvars, 5, densconss, 3, varsubprobdens, varmasterdens, conssubprobdens, consmasterdens) );
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

   GCGgetSubproblemVarsData(scip, decomp, nvars, nbinvars, nintvars, nimplvars, ncontvars, 2);
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

   GCGgetLinkingVarsData(scip, decomp, nvars, nbinvars, nintvars, nimplvars, ncontvars);
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

   SCIP_CALL_EXPECT( GCGgetVarLockData(scip, decomp, lockvars, 5, 2, sublockdown, sublockup, masterlockdown, masterlockup) );
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
