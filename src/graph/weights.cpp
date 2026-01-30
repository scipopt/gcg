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

/**@file   weights.cpp
 * @brief  weight class for graphs
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "graph/weights.h"

namespace gcg {

Weights::Weights(
      int varweight_,
      int vbinary_,
      int vcontinous_,
      int vinteger_,
      int vimplint_,
      int consweight_
   ): varweight(varweight_),
      vbinary(vbinary_),
      vcontinous(vcontinous_),
      vinteger(vinteger_),
      vimplint(vimplint_),
      consweight(consweight_)
{
   // TODO Auto-generated constructor stub

}

Weights::Weights()
: varweight(1),
  vbinary(1),
  vcontinous(1),
  vinteger(1),
  vimplint(1),
  consweight(1)
{

}

Weights::~Weights()
{
   // TODO Auto-generated destructor stub
}

int Weights::calculate(SCIP_CONS* cons) const
{ /*lint -e715*/
   return consweight;
}
int Weights::calculate(SCIP_VAR* var) const

{
   int weight;

   assert(var != NULL);

   switch ( SCIPvarGetType(var) ) {
   case SCIP_VARTYPE_CONTINUOUS:
      if( SCIPvarGetImplType(var) != SCIP_IMPLINTTYPE_NONE )
         weight = vimplint;
      else
         weight = vcontinous;
      break;
   case SCIP_VARTYPE_INTEGER:
      weight = vinteger;
      break;
   case SCIP_VARTYPE_BINARY:
      weight = vbinary;
      break;
   default:
      weight = varweight;
      break;
   }

   return weight;
}
} /* namespace gcg */
