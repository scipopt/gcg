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

/**@file   sepa_basis.h
 * @ingroup SEPARATORS
 * @brief  basis separator
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_BASIS_H__
#define __SCIP_SEPA_BASIS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the basis separator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeSepaBasis(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array of original cuts saved in the separator data */
extern
SCIP_ROW** GCGsepaBasisGetOrigcuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of original cuts saved in the separator data */
extern
int GCGsepaBasisGetNOrigcuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array of master cuts saved in the separator data */
extern
SCIP_ROW** GCGsepaBasisGetMastercuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of master cuts saved in the separator data */
extern
int GCGsepaBasisGetNMastercuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transforms cut in pricing variables to cut in original variables and adds it to newcuts array */
extern
SCIP_RETCODE GCGsepaBasisAddPricingCut(
   SCIP*                scip,               /**< SCIP data structure */
   int                  ppnumber,           /**< number of pricing problem */
   SCIP_ROW*            cut                 /**< cut to be added */
   );

/** Add cuts which are due to the latest objective function of the pricing problems
 *  (reduced cost non-negative) */
extern
SCIP_RETCODE SCIPsepaBasisAddPPObjConss(
   SCIP*                scip,               /**< SCIP data structure */
   int                  ppnumber,           /**< number of pricing problem */
   SCIP_Real            dualsolconv,        /**< dual solution corresponding to convexity constraint */
   SCIP_Bool            newcuts             /**< add cut to newcuts in sepadata? (otherwise add it just to the cutpool) */
   );

#ifdef __cplusplus
}
#endif

#endif
