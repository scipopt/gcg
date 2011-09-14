/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: presol_arrowheur.h,v 1.16 2010/01/04 20:35:45 bzfheinz Exp $"

/**@file   dec_arrowheur.h
 * @brief  arrowheur presolver
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DEC_ARROWHEUR_H__
#define __SCIP_DEC_ARROWHEUR_H__


#include "scip/scip.h"
#include "type_decomp.h"
#ifdef __cplusplus
extern "C" {
#endif

extern
/** creates the arrowheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionArrowheur(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
