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
#pragma ident "@(#) $Id: presol_borderheur.h,v 1.16 2010/01/04 20:35:45 bzfheinz Exp $"

/**@file   dec_borderheur.h
 * @brief  borderheur presolver
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DEC_BORDERHEUR_H__
#define __SCIP_DEC_BORDERHEUR_H__


#include "scip/scip.h"
#include "struct_decomp.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct SCIP_BorderheurData SCIP_BORDERHEURDATA;
/** creates the borderheur presolver and includes it in SCIP */
extern
SCIP_RETCODE SCIPBorderheurSetDecomp(
   SCIP* scip,
   SCIP_BORDERHEURDATA* borderheurdata,
   DECDECOMP* decdecomp
   );

extern
SCIP_RETCODE createBorderheurData(
      SCIP*         scip,
      SCIP_BORDERHEURDATA** borderheurdata
      );

extern
void freeBorderheurData(
      SCIP* scip,
      SCIP_BORDERHEURDATA** borderheurdata
   );

extern
/** creates the borderheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionBorderheur(
      SCIP*                 scip,                /**< SCIP data structure */
      SCIP_BORDERHEURDATA*  borderheurdata
   );

extern
SCIP_RETCODE detectAndBuildBordered(
      SCIP*                 scip,          /**< SCIP data structure */
      SCIP_BORDERHEURDATA*  borderheurdata, /**< presolver data data structure */
      SCIP_RESULT*          result
      );
#ifdef __cplusplus
}
#endif

#endif
