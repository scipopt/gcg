/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_stairheur.h
 * @ingroup DETECTORS
 * @brief  stairheur presolver
 * @author Martin Bergner
 * @author Mathias Luers
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DEC_STAIRHEUR_H__
#define __SCIP_DEC_STAIRHEUR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

enum Blockingtype
{
   DYNAMIC           = 1,     /**< Tries to minimize the number of linking variables */
   STATIC            = 2,     /**< Creates blocks with the same number of rows */
   ASSOONASPOSSIBLE  = 3      /**< Blocking is done in a way such that three adjacent blocks just do not overlap. This results in (almost) exclusively linking variables. */
};
typedef enum Blockingtype BLOCKINGTYPE;

extern
/** creates the stairheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionStairheur(
      SCIP*              scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif