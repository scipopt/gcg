/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_branchgcg.h
 * @brief  data structures for branching rules
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_BRANCHGCG_H__
#define __SCIP_STRUCT_BRANCHGCG_H__

#include "type_branchgcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** branching rule */
struct GCG_Branchrule
{
   SCIP_BRANCHRULE*      branchrule;         /**< pointer to the SCIP branching rule */
   GCG_DECL_BRANCHACTIVEMASTER ((*branchactivemaster));     /**< node activation method of branching rule */
   GCG_DECL_BRANCHDEACTIVEMASTER ((*branchdeactivemaster)); /**< node deactivation method of branching rule */
   GCG_DECL_BRANCHPROPMASTER ((*branchpropmaster));         /**< propagation method of branching rule */
   GCG_DECL_BRANCHMASTERSOLVED((*branchmastersolved));      /**< lp solved method of branching rule */
   GCG_DECL_BRANCHDATADELETE ((*branchdatadelete));         /**< deinitialization method of branching rule */
};


#ifdef __cplusplus
}
#endif

#endif
