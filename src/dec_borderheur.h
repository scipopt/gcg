/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_borderheur.h
 * @brief  borderheur detector
 * @author Martin Bergner
 * @ingroup DETECTORS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DEC_BORDERHEUR_H__
#define __SCIP_DEC_BORDERHEUR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

extern
/** creates the borderheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionBorderheur(
      SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
