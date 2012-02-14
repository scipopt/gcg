/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_branchgcg.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for branching rules in gcg projects
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_BRANCHGCG_H__
#define __SCIP_TYPE_BRANCHGCG_H__

#include "scip/def.h"
#include "scip/type_result.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GCG_BranchData GCG_BRANCHDATA;   /**< branching data */
typedef struct GCG_Branchrule GCG_BRANCHRULE;   /**< branching rule */


/** activation method for branchrule, called when a node in the master problem is activated,
 *  should perform changes to the current node's problem due to the branchdata
 *
 *  input:
 *  - scip            : SCIP main data structure of the master problem
 *  - branchdata      : the branching data
 */
#define GCG_DECL_BRANCHACTIVEMASTER(x) SCIP_RETCODE x (SCIP* scip, GCG_BRANCHDATA* branchdata)

/** deactivation method for branchrule, called when a node in the master problem is deactivated,
 *  should undo changes to the current node's problem due to the branchdata
 *
 *  input:
 *  - scip            : SCIP main data structure of the master problem
 *  - branchdata      : the branching data
 */
#define GCG_DECL_BRANCHDEACTIVEMASTER(x) SCIP_RETCODE x (SCIP* scip, GCG_BRANCHDATA* branchdata)

/** propagation method for branchrule, called when a node in the master problem is propagated,
 *  should perform propagation at the current node due to the branchdata
 *
 *  input:
 *  - scip            : SCIP main data structure of the master problem
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
#define GCG_DECL_BRANCHPROPMASTER(x) SCIP_RETCODE x (SCIP* scip, GCG_BRANCHDATA* branchdata, SCIP_RESULT* result)

/** method for branchrule, called when the master LP is solved at one node,
 *  can store pseudocosts for the branching decisions
 *
 *  input:
 *  - scip            : SCIP main data structure of the original problem
 *  - branchdata      : the branching data
 *  - newlowerbound   : the new local lower bound
 *
 */
#define GCG_DECL_BRANCHMASTERSOLVED(x) SCIP_RETCODE x (SCIP* scip, GCG_BRANCHDATA* branchdata, SCIP_Real newlowerbound)

/** frees branching data of an origbranch constraint (called when the origbranch constraint is deleted)
 *
 *  input:
 *    scip            : SCIP main data structure of the original problem
 *    branchdata      : pointer to the branching data to free
 */
#define GCG_DECL_BRANCHDATADELETE(x) SCIP_RETCODE x (SCIP* scip, GCG_BRANCHDATA** branchdata)

#ifdef __cplusplus
}
#endif

#endif
