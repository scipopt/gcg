/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
   GCG_DECL_BRANCHDATADELETE ((*branchdatadelete));
};


#ifdef __cplusplus
}
#endif

#endif
