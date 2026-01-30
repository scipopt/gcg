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

/**@file   type_branchgcg.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for branching rules in GCG projects
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_TYPE_BRANCHGCG_H__
#define GCG_TYPE_BRANCHGCG_H__

#include "scip/def.h"
#include "scip/type_result.h"
#include "scip/type_scip.h"
#include "gcg/type_extendedmasterconsdata.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GCG_BranchData GCG_BRANCHDATA;         /**< branching data */
typedef struct GCG_Branchrule GCG_BRANCHRULE;         /**< branching rule */
typedef struct GCG_BranchConsData GCG_BRANCHCONSDATA; /**< branching cons data */

/** type of variable bound: lower or upper bound */
enum GCG_BoundType
{
   GCG_BOUNDTYPE_LOWER = 0,            /**< lower bound */
   GCG_BOUNDTYPE_UPPER = 1,            /**< upper bound */
   GCG_BOUNDTYPE_FIXED = 2,            /**< variable fixed */
   GCG_BOUNDTYPE_NONE = 3              /**< no bound */
};
typedef enum GCG_BoundType GCG_BOUNDTYPE;

/** activation method for branchrule, called when a node in the master problem is activated,
 *  should perform changes to the current node's problem due to the branchdata
 *
 *  input:
 *  - gcg             : GCG main data structure
 *  - branchdata      : the branching data
 */
#define GCG_DECL_BRANCHACTIVEMASTER(x) SCIP_RETCODE x (GCG* gcg, GCG_BRANCHDATA* branchdata)

/** deactivation method for branchrule, called when a node in the master problem is deactivated,
 *  should undo changes to the current node's problem due to the branchdata
 *
 *  input:
 *  - gcg             : GCG main data structure
 *  - branchdata      : the branching data
 */
#define GCG_DECL_BRANCHDEACTIVEMASTER(x) SCIP_RETCODE x (GCG* gcg, GCG_BRANCHDATA* branchdata)

/** propagation method for branchrule, called when a node in the master problem is propagated,
 *  should perform propagation at the current node due to the branchdata
 *
 *  input:
 *  - gcg             : GCG main data structure
 *  - branchdata      : the branching data
 *  - node            : the activated node
 *  - result          : pointer to store the result of the propagation call
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_REDUCEDDOM : at least one domain reduction was found
 *  - SCIP_DIDNOTFIND : the propagator searched but did not find any domain reductions
 *  - SCIP_DIDNOTRUN  : the propagator was skipped
 *  - SCIP_DELAYED    : the propagator was skipped, but should be called again

 */
#define GCG_DECL_BRANCHPROPMASTER(x) SCIP_RETCODE x (GCG* gcg, GCG_BRANCHDATA* branchdata, SCIP_RESULT* result)

/** method for branchrule, called when the master LP is solved at one node,
 *  can store pseudocosts for the branching decisions
 *
 *  input:
 *  - gcg             : GCG main data structure
 *  - branchdata      : the branching data
 *  - newlowerbound   : the new local lower bound
 *
 */
#define GCG_DECL_BRANCHMASTERSOLVED(x) SCIP_RETCODE x (GCG* gcg, GCG_BRANCHDATA* branchdata, SCIP_Real newlowerbound)

/** frees branching data of an origbranch constraint (called when the origbranch constraint is deleted)
 *
 *  input:
 *    gcg             : GCG main data structure
 *    branchdata      : pointer to the branching data to free
 *    origbranch      : true iff an origbranch triggered this call
 *    force           : branch data must be deleted if true
 */
#define GCG_DECL_BRANCHDATADELETE(x) SCIP_RETCODE x (GCG* gcg, GCG_BRANCHDATA** branchdata, SCIP_Bool origbranch, SCIP_Bool force)

/** notify the branching rule that a new mastervariable was created while this node was active
 *
 *  input:
 *    gcg             : GCG main data structure
 *    branchdata      : the branching data
 *    mastervar       : pointer to the new master variable
 */
#define GCG_DECL_BRANCHNEWCOL(x) SCIP_RETCODE x (GCG* gcg, GCG_BRANCHDATA* branchdata, SCIP_VAR* mastervar)

/** get the extendedmasterconsdata created by this branching rule, if any
 *
 *  input:
 *    gcg             : GCG main data structure
 *    branchdata      : the branching data
 */
#define GCG_DECL_BRANCHGETEXTENDEDMASTERCONS(x) SCIP_RETCODE x (GCG* gcg, GCG_BRANCHDATA* branchdata, GCG_EXTENDEDMASTERCONSDATA** extendedmasterconsdata)

/** determine the coefficient of a column solution in the extended master cons
 *
 *  input:
 *    gcg             : GCG data structure
 *    branchdata      : branch data created by the branch rule
 *    extendedmasterconsdata   : the extended master cons data
 *    solvars         : array of column solution variables
 *    solvals         : array of column solution values
 *    nsolvars        : number of column solution variables and values
 *    probnr          : the pricing problem that the column belongs to
 *    gcgcol          : gcg column if available (or NULL)
 *    coef            : the calculated coefficient
 */
#define GCG_DECL_BRANCHGETEXTENDEDMASTERCONSCOEFF(x) SCIP_RETCODE x (GCG* gcg, GCG_BRANCHDATA* branchdata, GCG_EXTENDEDMASTERCONSDATA* extendedmasterconsdata, SCIP_VAR** solvars, SCIP_Real* solvals, int nsolvars, int probnr, GCG_COL* gcgcol, SCIP_Real* coef)

#ifdef __cplusplus
}
#endif

#endif
