/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
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
#include "gcg/gcg.h"


#ifdef __cplusplus
extern "C" {
#endif

/** creates the basis separator and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeSepaBasis(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the array of original cuts saved in the separator data */
GCG_EXPORT
SCIP_ROW** GCGsepaBasisGetOrigcuts(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the number of original cuts saved in the separator data */
GCG_EXPORT
int GCGsepaBasisGetNOrigcuts(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the array of master cuts saved in the separator data */
GCG_EXPORT
SCIP_ROW** GCGsepaBasisGetMastercuts(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the number of master cuts saved in the separator data */
GCG_EXPORT
int GCGsepaBasisGetNMastercuts(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** transforms cut in pricing variables to cut in original variables and adds it to newcuts array */
GCG_EXPORT
SCIP_RETCODE GCGsepaBasisAddPricingCut(
   GCG*                 gcg,                /**< GCG data structure */
   int                  ppnumber,           /**< number of pricing problem */
   SCIP_ROW*            cut                 /**< cut to be added */
   );

/** add cuts which are due to the latest objective function of the pricing problems
 *  (reduced cost non-negative) */
GCG_EXPORT
SCIP_RETCODE SCIPsepaBasisAddPPObjConss(
   GCG*                 gcg,                /**< GCG data structure */
   int                  ppnumber,           /**< number of pricing problem */
   SCIP_Real            dualsolconv,        /**< dual solution corresponding to convexity constraint */
   SCIP_Bool            newcuts             /**< add cut to newcuts in sepadata? (otherwise add it just to the cutpool) */
   );

#ifdef __cplusplus
}
#endif

#endif
