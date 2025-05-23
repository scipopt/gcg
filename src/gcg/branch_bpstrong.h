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

/**@file    branch_bpstrong.h
 * @brief   generic branch and price strong branching as described in
 *          Pecin, D., Pessoa, A., Poggi, M., Uchoa, E. Improved branch-cut-and-price for capacitated vehicle routing.
 *          In: Math. Prog. Comp. 9:61-100. Springer (2017).
 * @ingroup BRANCHINGRULES-GCG
 * @author  Oliver Gaul
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_BRANCH_BPSTRONG_H__
#define GCG_BRANCH_BPSTRONG_H__

#include "scip/scip.h"
#include "gcg/gcg.h"
#include "gcg/type_branchgcg.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** creates the xyz branching rule and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeBranchruleBPStrong(
    GCG*                 gcg                 /**< SCIP data structure */
);

SCIP_RETCODE GCGbranchSelectCandidateStrongBranchingOrig(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE       *origbranchrule,    /**< pointer storing original branching rule */
   SCIP_VAR              **branchvar,        /**< pointer to store output var pointer */
   SCIP_Bool             *upinf,             /**< pointer to store whether strong branching detected infeasibility in
                                                * the upbranch */
   SCIP_Bool             *downinf,           /**< pointer to store whether strong branching detected infeasibility in
                                                * the downbranch */
   SCIP_RESULT           *result,            /**< pointer to store result */
   SCIP_Bool             *stillusestrong     /**< pointer to store whether strong branching has reached a permanent
                                                * stopping condition for orig */
);

SCIP_RETCODE GCGbranchSelectCandidateStrongBranchingRyanfoster(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE*      rfbranchrule,       /**< Ryan-Foster branchrule */
   SCIP_VAR              **ovar1s,           /**< first elements of candidate pairs */
   SCIP_VAR              **ovar2s,           /**< second elements of candidate pairs */
   int                   *nspricingblock,    /**< pricing block numbers corresponding to input pairs */
   int                   npairs,             /**< number of input pairs */
   SCIP_VAR              **ovar1,            /**< pointer to store output var 1 pointer */
   SCIP_VAR              **ovar2,            /**< pointer to store output var 2 pointer */
   int                   *pricingblock,      /**< pointer to store output pricing block number */
   SCIP_Bool             *sameinf,           /**< pointer to store whether strong branching detected infeasibility in
                                                * the same branch */
   SCIP_Bool             *differinf,         /**< pointer to store whether strong branching detected infeasibility in
                                                * the differ branch */
   SCIP_RESULT           *result,            /**< pointer to store result */
   SCIP_Bool             *stillusestrong     /**< pointer to store whether strong branching has reached a permanent
                                                * stopping condition for Ryan-Foster */
);

#ifdef __cplusplus
}
#endif

#endif
