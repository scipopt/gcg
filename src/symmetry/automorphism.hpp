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

/**@file    automorphism.hpp
 * @brief   automorphism recognition (C++ interface)
 *
 * @author  Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_AUTOMORPHISM_HPP__
#define GCG_AUTOMORPHISM_HPP__

#include "gcg/class_partialdecomp.h"


/** compare two graphs w.r.t. automorphism */
SCIP_RETCODE cmpGraphPair(
   SCIP*                   scip,               /**< SCIP data structure */
   gcg::PARTIALDECOMP*     partialdec,         /**< partialdec the graphs should be compared for */
   int                     block1,             /**< index of first pricing prob */
   int                     block2,             /**< index of second pricing prob */
   SCIP_RESULT*            result,             /**< result pointer to indicate success or failure */
   SCIP_HASHMAP*           varmap,             /**< hashmap to save permutation of variables */
   SCIP_HASHMAP*           consmap,            /**< hashmap to save permutation of constraints */
   unsigned int            searchnodelimit,    /**< search node limit */
   unsigned int            generatorlimit      /**< generator limit */
);

#endif //GCG_AUTOMORPHISM_HPP__
