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

/**@file   weights.h
 * @brief  weight class for graphs
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef WEIGHTS_H_
#define WEIGHTS_H_
#include "objscip/objscip.h"

namespace gcg {

class Weights
{
protected:
   int varweight;                      /**< weight of a variable vertex */
   int vbinary;                        /**< weight of a binary variable vertex */
   int vcontinous;                     /**< weight of a continuous variable vertex */
   int vinteger;                       /**< weight of an integer variable vertex */
   int vimplint;                       /**< weight of an implicit integer variable vertex */
   int consweight;                     /**< weight of a constraint vertex */

public:
   Weights(
      int varweight_,
      int vbinary_,
      int vcontinous_,
      int vinteger_,
      int vimplint_,
      int consweight_
   );
   Weights();

   virtual ~Weights();
   int calculate(SCIP_CONS* cons) const;
   int calculate(SCIP_VAR* var) const;
};

} /* namespace gcg */
#endif /* WEIGHTS_H_ */
