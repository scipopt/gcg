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

/**@file   struct_branchgcg.h
 * @ingroup DATASTRUCTURES
 * @brief  data structures for branching rules
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_BRANCHGCG_H__
#define GCG_STRUCT_BRANCHGCG_H__

#include "gcg/type_branchgcg.h"

#include <scip/type_branch.h>

#ifdef __cplusplus
extern "C" {
#endif

struct GCG_BranchConsData
{
    GCG_BRANCHRULE* branchrule;
    GCG_BRANCHDATA* branchdata;
};

/** branching rule */
struct GCG_Branchrule
{
   SCIP_BRANCHRULE*      branchrule;                        /**< pointer to the SCIP branching rule */
   GCG_DECL_BRANCHACTIVEMASTER ((*branchactivemaster));     /**< node activation method of branching rule */
   GCG_DECL_BRANCHDEACTIVEMASTER ((*branchdeactivemaster)); /**< node deactivation method of branching rule */
   GCG_DECL_BRANCHPROPMASTER ((*branchpropmaster));         /**< propagation method of branching rule */
   GCG_DECL_BRANCHMASTERSOLVED((*branchmastersolved));      /**< lp solved method of branching rule */
   GCG_DECL_BRANCHDATADELETE ((*branchdatadelete));         /**< deinitialization method of branching rule */
   GCG_DECL_BRANCHNEWCOL ((*branchnewcol));                 /**< new column handler method of branching rule */
   GCG_DECL_BRANCHGETEXTENDEDMASTERCONS ((*branchgetextendedmastercons));     /**< extended master cons getter of branching rule */
   GCG_DECL_BRANCHGETEXTENDEDMASTERCONSCOEFF ((*branchgetextendedmasterconscoeff));             /**< column coefficient calculation method for extended master conss */
};

#ifdef __cplusplus
}
#endif

#endif
