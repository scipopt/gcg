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
#pragma ident "@(#) $Id: presol_stairheur.h,v 1.16 2010/01/04 20:35:45 bzfheinz Exp $"

/**@file   presol_stairheur.h
 * @brief  stairheur presolver
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_STAIRHEUR_H__
#define __SCIP_PRESOL_STAIRHEUR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_StairheurData SCIP_STAIRHEURDATA;

/** execution method of presolver */
extern
SCIP_RETCODE detectStructureStairheur(
      SCIP*                 scip,           /**< SCIP data structure       */
      SCIP_STAIRHEURDATA*   stairheurdata,  /**< stairheur data structure  */
      SCIP_RESULT*          result          /**< pointer to hold result    */
      );

/** creates the stairheur detection and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeDetectionStairheur(
      SCIP*                 scip, /**< SCIP data structure */
      SCIP_STAIRHEURDATA*   stairheurdata
   );
/** allocates and initializes the stairheurdata */
extern
SCIP_RETCODE createStairheurData(
      SCIP*                scip,           /**< SCIP data structure */
      SCIP_STAIRHEURDATA** stairheurdata   /**< stairheur data structure */
   );

/** frees the stairheurdata */
extern
void freeStairheurData(
      SCIP*                scip,           /**< SCIP data structure */
      SCIP_STAIRHEURDATA** stairheurdata   /**< stairheur data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
