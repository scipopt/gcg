/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file   sepa_original.h
 * @ingroup SEPARATORS
 * @brief  original separator
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_SEPA_ORIGINAL_H__
#define GCG_SEPA_ORIGINAL_H__


#include "scip/scip.h"
#include "def.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the original separator and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE SCIPincludeSepaOriginal(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array of original cuts in the original problem saved in the separator data */
GCG_EXPORT
SCIP_ROW** GCGsepaGetOriginalSepaOrigcuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of cuts saved in the separator data */
GCG_EXPORT
int GCGsepaGetNOriginalSepaCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array of original cuts in the master problem saved in the separator data */
GCG_EXPORT
SCIP_ROW** GCGsepaGetOriginalSepaMastercuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds given original cut in both the original and master problem to master separator data */
GCG_EXPORT
SCIP_RETCODE GCGsepaAddOriginalSepaCuts(
   SCIP*                scip,               /**< SCIP data structure */
   SCIP_ROW*            origcut,            /**< pointer to orginal cut in the original problem */
   SCIP_ROW*            mastercut           /**< pointer to original cut in the master problem */
   );

/** checks whether a given original cut in the original problem is already known */
SCIP_Bool GCGsepaOriginalSepaOrigcutExists(
   SCIP*                scip,            /**< SCIP data structure */
   SCIP_ROW*            origcut          /**< pointer to orginal cut in the original problem */
   );

#ifdef __cplusplus
}
#endif

#endif
