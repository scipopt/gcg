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

/**@file   pricingprob.h
 * @ingroup PRICING_PRIV
 * @brief  private methods for working with pricing problems, to be used by the pricing controller only
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PRICINGPROB_H__
#define GCG_PRICINGPROB_H__


#include "gcg/struct_pricingprob.h"
#include "gcg/type_pricingprob.h"

#include "gcg/pricer_gcg.h"
#include "gcg/type_colpool.h"
#include "gcg/type_pricestore_gcg.h"
#include "gcg/type_pricingjob.h"
#include "gcg/type_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create a pricing problem */
SCIP_RETCODE GCGpricingprobCreate(
   GCG*                  gcg,                /**< GCG data structure (master problem) */
   GCG_PRICINGPROB**     pricingprob,        /**< pricing problem to be created */
   SCIP*                 pricingscip,        /**< SCIP data structure of the corresponding pricing problem */
   int                   probnr,             /**< index of the corresponding pricing problem */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
);

/** free a pricing problem */
void GCGpricingprobFree(
   GCG*                  gcg,                /**< GCG data structure (master problem) */
   GCG_PRICINGPROB**     pricingprob         /**< pricing problem to be freed */
);

/** initialize pricing problem at the beginning of the pricing round */
void GCGpricingprobInitPricing(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** uninitialize pricing problem at the beginning of the pricing round */
void GCGpricingprobExitPricing(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
   );

/** add generic branching data (constraint and dual value) to the current pricing problem */
SCIP_RETCODE GCGpricingprobAddGenericBranchData(
   GCG*                  gcg,                /**< GCG data structure (master problem) */
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   SCIP_CONS*            branchcons,         /**< generic branching constraint */
   SCIP_Real             branchdual          /**< corresponding dual solution value */
   );

/** reset the pricing problem statistics for the current pricing round */
void GCGpricingprobReset(
   GCG*                  gcg,                /**< GCG data structure (master problem) */
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** update solution information of a pricing problem */
void GCGpricingprobUpdate(
   GCG*                  gcg,                /**< GCG data structure (master problem) */
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   GCG_PRICINGSTATUS     status,             /**< status of last pricing job */
   SCIP_Real             lowerbound,         /**< new lower bound */
   int                   nimpcols            /**< number of new improving columns */
   );

/** add the information that the next branching constraint must be added */
void GCGpricingprobNextBranchcons(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** set the lower bound of a pricing job */
void GCGpricingjobSetLowerbound(
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   SCIP_Real             lowerbound          /**< new lower bound */
   );

#ifdef __cplusplus
}
#endif

#endif
