/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   sepa_original.h
 * @ingroup SEPARATORS
 * @brief  original separator
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_SEPA_ORIGINAL_H__
#define GCG_SEPA_ORIGINAL_H__


#include "scip/scip.h"
#include "gcg/gcg.h"


#ifdef __cplusplus
extern "C" {
#endif

/** creates the original separator and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeSepaOriginal(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the array of original cuts in the original problem saved in the separator data */
GCG_EXPORT
SCIP_ROW** GCGsepaGetOriginalSepaOrigcuts(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the number of cuts saved in the separator data */
GCG_EXPORT
int GCGsepaGetNOriginalSepaCuts(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the array of original cuts in the master problem saved in the separator data */
GCG_EXPORT
SCIP_ROW** GCGsepaGetOriginalSepaMastercuts(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** adds given original cut in both the original and master problem to master separator data */
GCG_EXPORT
SCIP_RETCODE GCGsepaAddOriginalSepaCuts(
   GCG*                 gcg,                /**< GCG data structure */
   SCIP_ROW*            origcut,            /**< pointer to orginal cut in the original problem */
   SCIP_ROW*            mastercut           /**< pointer to original cut in the master problem */
   );

/** checks whether a given original cut in the original problem is already known */
SCIP_Bool GCGsepaOriginalSepaOrigcutExists(
   GCG*                 gcg,             /**< GCG data structure */
   SCIP_ROW*            origcut          /**< pointer to orginal cut in the original problem */
   );

#ifdef __cplusplus
}
#endif

#endif
