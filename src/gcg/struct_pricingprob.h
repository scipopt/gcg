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

/**@file   struct_pricingprob.h
 * @ingroup DATASTRUCTURES
 * @brief  data structure to store pricing problem information
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_PRICINGPROB_H_
#define GCG_STRUCT_PRICINGPROB_H_

#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/scip.h"

#include "gcg/type_pricingprob.h"
#include "gcg/type_gcgcol.h"
#include "gcg/type_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

/** pricing problem data structure */
struct GCG_PricingProb
{
   /* problem data */
   SCIP*                pricingscip;        /**< SCIP data structure */
   int                  probnr;             /**< (block) index of the corresponding pricing problem */

   /* generic branching information */
   SCIP_CONS**          branchconss;        /**< stack of generic branching constraints */
   SCIP_Real*           branchduals;        /**< corresponding dual solution values */
   int                  nbranchconss;       /**< number of generic branching constraints */
   int                  branchconsssize;    /**< size of generic branching constraints array */
   int                  branchconsidx;      /**< lowest index generic branching constraint that is considered */
   SCIP_Bool            consisadded;        /**< flag to indicate whether this constraint has already been added */

   /* result values */
   GCG_PRICINGSTATUS    status;             /**< current status of the pricing problem */
   SCIP_Real            lowerbound;         /**< lower bound obtained by solving the pricing problem */
   int                  nimpcols;           /**< number of improving columns found in the current pricing round */

   /* statistics */
   int                  nsolves;            /**< number of times the pricing problem was solved during the loop */
   int*                 ncolsround;         /**< number of improving columns found in the last rounds */
   int                  maxcolsround;       /**< capacity of ncolsround */
};

#ifdef __cplusplus
}
#endif

#endif /* STRUCT_PRICINGPROB_H_ */
