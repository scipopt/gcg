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

/**@file    branchgcg.c
 * @brief   methods for branching rules in GCG projects
 * @author  Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/branchgcg.h"
#include "gcg/struct_branchgcg.h"

struct GCG_ExtendedMasterConsDataData
{
    GCG_BRANCHRULE* branchrule;
    GCG_BRANCHDATA* branchdata;
};

/** creates an extended master cons for a cons created by a branch rule */
SCIP_RETCODE GCGbranchCreateExtendedmastercons(
    GCG*                             gcg,                             /**< GCG data structure */
    GCG_BRANCHRULE*                  branchrule,                      /**< GCG branch rule */
    GCG_EXTENDEDMASTERCONSDATA**     extendedmasterconsdata,          /**< pointer to store the extended master cons data */
    SCIP_CONS*                       cons,                            /**< constraint in the master problem that represents the extended master cons */
    GCG_PRICINGMODIFICATION**        pricingmodifications,            /**< pricing modifications for the extended master cons */
    int                              npricingmodifications,           /**< number of pricing modifications for the extended master cons */
    GCG_BRANCHDATA*                  branchdata                       /**< any data that might be required to calculate the coefficient of a column solution */
    )
{
    GCG_EXTENDEDMASTERCONSDATADATA* data = NULL;

    SCIP_CALL( SCIPallocBlockMemory(GCGgetMasterprob(gcg), &data) );
    data->branchrule = branchrule;
    data->branchdata = branchdata;
    SCIP_CALL( GCGextendedmasterconsCreateFromCons(gcg, extendedmasterconsdata, GCG_EXTENDEDMASTERCONSTYPE_BRANCH_CONS, cons, pricingmodifications, npricingmodifications, data) );
    return SCIP_OKAY;
}

/** calculate the coefficient of a column solution in the extended master cons */
SCIP_RETCODE GCGbranchGetExtendedmasterconsCoeff(
    GCG*                             gcg,                          /**< GCG data structure */
    GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata,       /**< extended master cons data */
    SCIP_VAR**                       solvars,                      /**< array of column solution variables */
    SCIP_Real*                       solvals,                      /**< array of column solution values */
    int                              nsolvars,                     /**< number of column solution variables and values */
    int                              probnr,                       /**< the pricing problem that the column belongs to */
    SCIP_Real*                       coeff                         /**< pointer to store the coefficient */
    )
{
    GCG_EXTENDEDMASTERCONSDATADATA* data;
    assert(extendedmasterconsdata != NULL);
    assert(GCGextendedmasterconsGetType(extendedmasterconsdata) == GCG_EXTENDEDMASTERCONSTYPE_BRANCH_CONS);

    data = GCGextendedmasterconsGetData(extendedmasterconsdata);
    assert(data != NULL);
    SCIP_CALL( data->branchrule->branchgetextendedmasterconscoeff(gcg, data->branchdata, extendedmasterconsdata, solvars, solvals, nsolvars, probnr, coeff) );
    return SCIP_OKAY;
}

/** frees the branch data stored in the extened master cons data */
SCIP_RETCODE GCGbranchFreeExtendedmasterconsBranchData(
    GCG*                             gcg,                          /**< GCG data structure */
    GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata        /**< extended master cons data */
    )
{
    GCG_EXTENDEDMASTERCONSDATADATA* data;
    assert(extendedmasterconsdata != NULL);
    assert(GCGextendedmasterconsGetType(extendedmasterconsdata) == GCG_EXTENDEDMASTERCONSTYPE_BRANCH_CONS);

    data = GCGextendedmasterconsGetData(extendedmasterconsdata);
    assert(data != NULL);
    SCIPfreeBlockMemory(GCGgetMasterprob(gcg), &data);
    return SCIP_OKAY;
}
