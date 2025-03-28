/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
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

/**@file   branch_orig.h
 * @brief  branching rule for original problem in GCG
 * @author Gerald Gamrath
 * @ingroup BRANCHINGRULES-GCG
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_BRANCH_ORIG_H__
#define GCG_BRANCH_ORIG_H__

#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the branching on original variable branching rule and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeBranchruleOrig(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** get the original variable on which the branching was performed */
SCIP_VAR* GCGbranchOrigGetOrigvar(
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   );

/** get the type of the new bound which resulted of the performed branching */
GCG_BOUNDTYPE GCGbranchOrigGetBoundtype(
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   );

/** get the new bound which resulted of the performed branching */
SCIP_Real GCGbranchOrigGetNewbound(
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   );

/** updates extern branching candidates before branching */
SCIP_RETCODE GCGbranchOrigUpdateExternBranchcands(
   GCG*                  gcg                /**< GCG data structure */
);

#ifdef __cplusplus
}
#endif

#endif
