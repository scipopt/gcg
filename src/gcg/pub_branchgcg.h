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

/**@file    pub_branchgcg.h
 * @ingroup PUBLICCOREAPI
 * @brief   public methods for branching rules in GCG projects
 * @author  Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_PUB_BRANCHGCG_H_
#define GCG_PUB_BRANCHGCG_H_

#include "gcg/type_branchgcg.h"
#include "scip/scip.h"
#include "gcg/def.h"
#include "gcg/type_gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates an extended master cons for a cons created by a branch rule */
GCG_EXPORT
SCIP_RETCODE GCGbranchCreateExtendedmastercons(
    GCG*                             gcg,                             /**< GCG data structure */
    GCG_BRANCHRULE*                  branchrule,                      /**< GCG branch rule */
    GCG_EXTENDEDMASTERCONSDATA**     extendedmasterconsdata,          /**< pointer to store the extended master cons data */
    SCIP_CONS*                       cons,                            /**< constraint in the master problem that represents the extended master cons */
    GCG_PRICINGMODIFICATION**        pricingmodifications,            /**< pricing modifications for the extended master cons */
    int                              npricingmodifications,           /**< number of pricing modifications for the extended master cons */
    GCG_BRANCHDATA*                  branchdata                       /**< any data that might be required to calculate the coefficient of a column solution */
    );

/** calculate the coefficient of a column solution in the extended master cons */
GCG_EXPORT
SCIP_RETCODE GCGbranchGetExtendedmasterconsCoeff(
    GCG*                             gcg,                          /**< GCG data structure */
    GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata,       /**< extended master cons data */
    SCIP_VAR**                       solvars,                      /**< array of column solution variables */
    SCIP_Real*                       solvals,                      /**< array of column solution values */
    int                              nsolvars,                     /**< number of column solution variables and values */
    int                              probnr,                       /**< the pricing problem that the column belongs to */
    SCIP_Real*                       coeff                         /**< pointer to store the coefficient */
    );

/** frees the branch data stored in the extened master cons data */
GCG_EXPORT
SCIP_RETCODE GCGbranchFreeExtendedmasterconsBranchData(
    GCG*                             gcg,                          /**< GCG data structure */
    GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata        /**< extended master cons data */
    );

#ifdef __cplusplus
}
#endif

#endif
