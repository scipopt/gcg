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

/**@file   decomp.h
 * @ingroup DECOMP
 * @brief  private methods for working with decomp structures
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_DECOMP_H__
#define GCG_DECOMP_H__

#include "gcg/type_decomp.h"
#include "gcg/type_detector.h"
#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets the detectors for the given decomposition 
 * 
 * @note make sure you know what you are doing, only use at initialization
*/
GCG_EXPORT
SCIP_RETCODE GCGdecompSetDetectorChain(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   GCG_DETECTOR**        detectors,          /**< new detector chain */
   int                   ndetectors          /**< number of new detectors (i.e. length of the detector array) */
   );

/** sets the id of the original partialdec */
GCG_EXPORT
void GCGdecompSetPartialdecID(
   GCG_DECOMP*           decomp,              /**< decomposition data structure */
   int                   id                   /**< ID of partialdec */
   );

#ifdef __cplusplus
}
#endif

#endif
