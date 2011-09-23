/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   struct_branchgcg.h
 * @brief  datastructures for branching rules
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
   SCIP_BRANCHRULE*       branchrule;
   GCG_DECL_BRANCHACTIVEMASTER ((*branchactivemaster));
   GCG_DECL_BRANCHDEACTIVEMASTER ((*branchdeactivemaster));
   GCG_DECL_BRANCHPROPMASTER ((*branchpropmaster));
   GCG_DECL_BRANCHMASTERSOLVED((*branchmastersolved));
   GCG_DECL_BRANCHDATADELETE ((*branchdatadelete));
};


#ifdef __cplusplus
}
#endif

#endif
