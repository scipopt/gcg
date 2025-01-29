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

/**@file   benders_gcg.h
 * @ingroup BENDERS-GCG
 * @brief  GCG Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERS_GCG_H__
#define __SCIP_BENDERS_GCG_H__


#include "def.h"
#include "scip/scip.h"
#include "scip/bendersdefcuts.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the GCG Benders' decomposition and includes it in SCIP
 *
 *  @ingroup BendersIncludes
 */
GCG_EXPORT
SCIP_RETCODE SCIPincludeBendersGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origprob            /**< the SCIP instance of the original problem */
   );

/**@addtogroup BENDERS-GCG
 *
 * @{
 */

/** returns the last relaxation solution */
GCG_EXPORT
SCIP_SOL* SCIPbendersGetRelaxSol(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   );

/** returns the original problem for the given master problem */
GCG_EXPORT
SCIP* GCGbendersGetOrigprob(
   SCIP*                 masterprob          /**< the master problem SCIP instance */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
