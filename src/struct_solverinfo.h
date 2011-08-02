/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_solverinfo.h 198 2011-01-06 16:58:56Z ggamrath $"

/**@file   struct_solverinfo.h
 * @brief  datastructures for solverinfos
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SOLVERINFO_H__
#define __SCIP_STRUCT_SOLVERINFO_H__

#include "type_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

/** branching rule */
struct GCG_SolverInfo
{
   pthread_mutex_t access_masterscip;
   pthread_mutex_t update_count;
   pthread_cond_t update_cond;
   int* queue;
   int nqueueentries;

   int count;
};


#ifdef __cplusplus
}
#endif

#endif
