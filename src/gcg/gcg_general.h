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
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   gcg_general.h
 * @brief  gcg general public methods
 * @author Steffan Schlein
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_GENERAL_H__
#define GCG_GENERAL_H__


#include "def.h"
#include "scip/def.h"
#include "scip/type_scip.h"

#include "def.h"

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
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

#ifdef __cplusplus
}
#endif

#endif
