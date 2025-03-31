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

/**@file   graphtest.cpp
 * @brief  Description
 * @author bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "graphtest.h"
#include <fstream>
#include <cerrno>
#include <cstdio>

void GraphTest::SetUp() {
  SCIP_CALL_ABORT( SCIPcreate(&scip) );
  SCIP_CALL_ABORT( SCIPincludeGcgPlugins(scip) );
  SCIP_CALL_ABORT( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
  SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrgpartition/enabled", FALSE) );
  SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hrcgpartition/enabled", FALSE) );
  SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/hcgpartition/enabled", FALSE) );
  SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/random/enabled", FALSE) );
  SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "detection/detectors/staircase/enabled", FALSE) );
  SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
  SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "prob") );
}

void GraphTest::TearDown() {
  SCIP_CALL_ABORT( SCIPfree(&scip) );
}

SCIP_RETCODE GraphTest::createVar(const char * str) {
   SCIP_VAR* var;
   char* endptr;
   SCIP_Bool success;
   SCIP_CALL( SCIPparseVar(scip, &var, str, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL, &endptr, &success) );
   assert(success);
   SCIP_CALL( SCIPaddVar(scip, var) );
   SCIP_CALL( SCIPreleaseVar(scip, &var) );
   return SCIP_OKAY;
}

SCIP_RETCODE GraphTest::createCons(const char * str) {
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_CALL( SCIPparseCons(scip, &cons, str, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   assert(success);
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   return SCIP_OKAY;
}

void GraphTest::parseFile(
   const char *str,
   std::vector<int> &array)
{
   std::ifstream stream(str);
   stream.exceptions( std::ios::failbit );

   ASSERT_TRUE(stream.good());

   int input;

   for( size_t i = 0; i < array.size(); ++i )
   {
      ASSERT_FALSE(stream.eof());
      ASSERT_TRUE(static_cast<bool>(stream >> input));
      ASSERT_EQ(array[i], input);
   }
   stream.close();
}
