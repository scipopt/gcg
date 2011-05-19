/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_xxx.h,v 1.15 2010/01/04 20:35:41 bzfheinz Exp $"

/**@file   heur_restmaster.h
 * @brief  restricted master primal heuristic
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_RESTMASTER_H__
#define __SCIP_HEUR_RESTMASTER_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the restricted master primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurRestmaster(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
