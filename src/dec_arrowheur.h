/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_arrowheur.h
 * @brief  arrowheur presolver
 * @author Martin Bergner
 * @ingroup DETECTORS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DEC_ARROWHEUR_H__
#define __SCIP_DEC_ARROWHEUR_H__

#include "scip/scip.h"
#include "type_decomp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the arrowheur presolver and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeDetectionArrowheur(
   SCIP* scip                 /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
