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
#pragma ident "@(#) $Id: presol_arrowheur.h,v 1.16 2010/01/04 20:35:45 bzfheinz Exp $"

/**@file   dec_arrowheur.h
 * @brief  arrowheur presolver
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DEC_ARROWHEUR_H__
#define __SCIP_DEC_ARROWHEUR_H__


#include "scip/scip.h"
#include "struct_decomp.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct SCIP_ArrowheurData SCIP_ARROWHEURDATA;
/** creates the arrowheur presolver and includes it in SCIP */
extern
SCIP_RETCODE SCIPArrowHeurSetDecomp(
   SCIP* scip,
   SCIP_ARROWHEURDATA* arrowheurdata,
   DECDECOMP* decdecomp
   );

extern
SCIP_RETCODE createArrowheurData(
      SCIP*         scip,
      SCIP_ARROWHEURDATA** arrowheurdata
      );

extern
void freeArrowheurData(
      SCIP* scip,
      SCIP_ARROWHEURDATA** arrowheurdata
   );

extern
/** creates the arrowheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionArrowheur(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_ARROWHEURDATA*   arrowheurdata
   );

extern
SCIP_RETCODE detectAndBuildArrowHead(
      SCIP*                scip,          /**< SCIP data structure */
      SCIP_ARROWHEURDATA*  arrowheurdata, /**< presolver data data structure */
      SCIP_RESULT*         result
      );
#ifdef __cplusplus
}
#endif

#endif
