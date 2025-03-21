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

/**@file   dec_stairheur.h
 * @ingroup DETECTORS
 * @brief  detector for staircase structures via ROC algorithms
 * @author Martin Bergner
 * @author Mathias Luers
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_DEC_STAIRHEUR_H__
#define GCG_DEC_STAIRHEUR_H__


#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

enum Blockingtype
{
   DYNAMIC           = 1,     /**< Tries to minimize the number of linking variables */
   STATIC            = 2,     /**< Creates blocks with the same number of rows */
   ASSOONASPOSSIBLE  = 3      /**< Blocking is done in a way such that three adjacent blocks just do not overlap. This results in (almost) exclusively linking variables. */
};
typedef enum Blockingtype BLOCKINGTYPE;

/** creates the stairheur detector and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeDetectorStairheur(
      GCG*               gcg                 /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
