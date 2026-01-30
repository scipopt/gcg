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

/**@file   type_pricingcb.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for the pricing callback plugin
 * @author Stephen J. Maher
 *
 *  This file defines the interface for pricing callback plugins implemented in C.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_PRICINGCB_H__
#define __SCIP_TYPE_PRICINGCB_H__

#include "scip/def.h"
#include "scip/type_pricer.h"
#include "scip/type_result.h"
#include "gcg/type_pricetype.h"
#include "gcg/type_gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GCG_Pricingcb GCG_PRICINGCB;              /**< the pricing callback structure */
typedef struct GCG_PricingcbData GCG_PRICINGCBDATA;      /**< locally defined pricing callback data */


/** copy method for the pricing callback plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - pricingcb       : the pricing callback structure
 */
#define GCG_DECL_PRICINGCBCOPY(x) SCIP_RETCODE x (GCG* gcg, GCG_PRICINGCB* pricingcb)

/** destructor of the pricing callback to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - pricingcb       : the pricing callback structure
 */
#define GCG_DECL_PRICINGCBFREE(x) SCIP_RETCODE x (GCG* gcg, GCG_PRICINGCB* pricingcb)

/** initialization method of the pricing callback (called after problem was transformed)
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - pricingcb       : the pricing callback structure
 */
#define GCG_DECL_PRICINGCBINIT(x) SCIP_RETCODE x (GCG* gcg, GCG_PRICINGCB* pricingcb)

/** deinitialization method of the pricing callback (called before transformed problem is freed)
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - pricingcb       : the pricing callback structure
 */
#define GCG_DECL_PRICINGCBEXIT(x) SCIP_RETCODE x (GCG* gcg, GCG_PRICINGCB* pricingcb)

/** solving process initialization method of the pricing callback (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - pricingcb       : the pricing callback structure
 */
#define GCG_DECL_PRICINGCBINITSOL(x) SCIP_RETCODE x (GCG* gcg, GCG_PRICINGCB* pricingcb)

/** solving process deinitialization method of the pricing callback (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The pricing callback should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - pricingcb       : the pricing callback structure
 */
#define GCG_DECL_PRICINGCBEXITSOL(x) SCIP_RETCODE x (GCG* gcg, GCG_PRICINGCB* pricingcb)

/** the pre-pricing method of the pricing callback technique
 *
 *  This method is called immediately before the pricing is performed in the GCG pricer. At this point, it is possible
 *  to modify the solving data is used within the pricing for new variables. Any data that is modified should be
 *  reverted in the post-pricing method.
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - pricingcb       : the pricing callback structure
 *  - pricer          : the pricer structure
 *  - type            : the pricing type, either farkas or reduced cost pricing
 *  - abort           : should the pricing be aborted. Care must be taken when setting this flag, since aborting the
 *                      pricing may lead to early branching and sub-optimal solution.
 *  - result          : pointer to store the result of the pre-pricing method
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_DIDNOTRUN  : if the pre-pricing method was not run.
 *  - SCIP_SUCCESS    : the pre-pricing method was executed successfully.
 */
#define GCG_DECL_PRICINGCBPREPRICING(x) SCIP_RETCODE x (GCG* gcg, GCG_PRICINGCB* pricingcb,\
   SCIP_PRICER* pricer, GCG_PRICETYPE type, SCIP_Bool* abort, SCIP_RESULT* result)

/** the post-pricing method of the pricing callback technique
 *
 *  This method is called immediately after the pricing is performed in the GCG pricer. This method should be used to
 *  revert any changes made to in the pre-pricing method.
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - pricingcb       : the pricing callback structure
 *  - pricer          : the pricer structure
 *  - type            : the pricing type, either farkas or reduced cost pricing
 *  - result          : pointer to store the result of the pre-pricing method
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_DIDNOTRUN  : if the post-pricing method was not run.
 *  - SCIP_SUCCESS    : the post-pricing method was executed successfully.
 */
#define GCG_DECL_PRICINGCBPOSTPRICING(x) SCIP_RETCODE x (GCG* gcg, GCG_PRICINGCB* pricingcb,\
   SCIP_PRICER* pricer, GCG_PRICETYPE type, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
