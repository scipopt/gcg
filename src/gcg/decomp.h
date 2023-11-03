/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       */
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

/**@file   decomp.h
 * @ingroup DECOMP
 * @brief  private methods for working with decomp structures
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_DECOMP_H__
#define GCG_DECOMP_H__

#include "type_decomp.h"
#include "type_detector.h"
#include "scip/scip.h"
#include "gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets the detectors for the given decomposition 
 * 
 * @note make sure you know what you are doing, only use at initialization
*/
GCG_EXPORT
SCIP_RETCODE GCGdecompSetDetectorChain(
   SCIP*                 scip,               /**< SCIP data structure */
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
