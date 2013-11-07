/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_base.h
 * @ingroup SEPARATORS
 * @brief  base separator
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_BASE_H__
#define __SCIP_SEPA_BASE_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the base separator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeSepaBase(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array of original cuts saved in the separator data */
extern
SCIP_ROW** GCGsepaBaseGetOrigcuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of original cuts saved in the separator data */
extern
int GCGsepaBaseGetNOrigcuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array of master cuts saved in the separator data */
extern
SCIP_ROW** GCGsepaBaseGetMastercuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of master cuts saved in the separator data */
extern
int GCGsepaBaseGetNMastercuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transforms cut in pricing variables to cut in original variables and adds it to newcuts array */
extern
SCIP_RETCODE GCGsepaBaseAddPricingCut(
   SCIP*                scip,
   int                  ppnumber,
   SCIP_ROW*            cut
   );

#ifdef __cplusplus
}
#endif

#endif
