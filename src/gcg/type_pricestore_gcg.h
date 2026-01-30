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

/**@file   type_pricestore_gcg.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for storing priced cols
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_TYPE_PRICESTORE_H__
#define __GCG_TYPE_PRICESTORE_H__

#ifdef __cplusplus
extern "C" {
#endif

/** possible settings for specifying the solution for which cuts are selected */
enum GCG_Efficiacychoice
{
   GCG_EFFICIACYCHOICE_DANTZIG = 0,          /**< use Dantzig's rule (reduced cost) to base efficacy on */
   GCG_EFFICIACYCHOICE_STEEPESTEDGE = 1,     /**< use steepest edge rule s( to base efficacy on */
   GCG_EFFICIACYCHOICE_LAMBDA = 2            /**< use lambda pricing to base efficacy on */
};
typedef enum GCG_Efficiacychoice GCG_EFFICIACYCHOICE;

typedef struct GCG_PriceStore GCG_PRICESTORE;     /**< storage for priced variables */

#ifdef __cplusplus
}
#endif

#endif
