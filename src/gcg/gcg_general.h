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

/**@file   gcg_general.h
 * @brief  gcg general public methods
 * @author Steffan Schlein
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_GENERAL_H__
#define GCG_GENERAL_H__



#include "scip/def.h"
#include "scip/type_scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns complete GCG version number in the format "major . minor tech"
 *
 *  @return complete GCG version
 */
GCG_EXPORT
SCIP_Real GCGversion(
   void
   );

/** returns GCG major version
 *
 *  @return major GCG version
 */
GCG_EXPORT
int GCGmajorVersion(
   void
   );

/** returns GCG minor version
 *
 *  @return minor GCG version
 */
GCG_EXPORT
int GCGminorVersion(
   void
   );

/** returns GCG technical version
 *
 *  @return technical GCG version
 */
GCG_EXPORT
int GCGtechVersion(
   void
   );

/** returns GCG sub version number
 *
 *  @return subversion GCG version
 */
GCG_EXPORT
int GCGsubversion(
   void
   );

/** prints out GCG version
 */
GCG_EXPORT
void GCGprintVersion(
   GCG*                  gcg,                /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

#ifdef __cplusplus
}
#endif

#endif
