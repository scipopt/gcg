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

/**@file   struct_pricingjob.h
 * @ingroup DATASTRUCTURES
 * @brief  data structure for pricing jobs
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_PRICINGJOB_H_
#define GCG_STRUCT_PRICINGJOB_H_

#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/scip.h"

#include "gcg/type_pricingjob.h"
#include "gcg/type_pricingprob.h"
#include "gcg/type_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

/** pricing job data structure */
struct GCG_PricingJob
{
   /* problem data */
   GCG_PRICINGPROB*     pricingprob;        /**< data structure of the corresponding pricing problem */
   GCG_SOLVER*          solver;             /**< solver with which to solve the pricing problem */

   /* strategic parameters */
   int                  chunk;              /**< chunk the pricing job belongs to */
   SCIP_Real            score;              /**< current score of the pricing job */
   SCIP_Bool            heuristic;          /**< shall the pricing problem be solved heuristically? */
   int                  nheuriters;         /**< number of times the pricing job was performed heuristically */
   SCIP_Bool            solverchanged;      /**< was the solver changed after the last solver call? */
};

#ifdef __cplusplus
}
#endif

#endif /* STRUCT_PRICINGJOB_H_ */
