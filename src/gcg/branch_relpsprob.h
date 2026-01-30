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

/**@file   branch_relpsprob.h
 * @ingroup BRANCHINGRULES-GCG
 * @brief  reliable pseudo costs branching rule
 * @author Tobias Achterberg
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_BRANCH_RELPSPROB_H__
#define GCG_BRANCH_RELPSPROB_H__


#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the reliable pseudo cost braching rule and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeBranchruleRelpsprob(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** execution reliability pseudo cost probing branching with the given branching candidates */
SCIP_RETCODE GCGgetRelpsprobBranchVar(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_VAR**            branchcands,        /**< brancing candidates */
   SCIP_Real*            branchcandssol,     /**< solution value for the branching candidates */
   int                   nbranchcands,       /**< number of branching candidates */
   int                   nvars,              /**< number of variables to be watched by bdchgdata */
   SCIP_RESULT*          result,             /**< pointer to the result of the execution */
   SCIP_VAR**            branchvar           /**< pointer to the variable to branch on */
   );

#ifdef __cplusplus
}
#endif

#endif
