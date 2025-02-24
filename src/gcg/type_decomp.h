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

/**@file   type_decomp.h
 * @ingroup DECOMP
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for decomposition information in GCG projects
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_TYPE_DECOMP_H__
#define GCG_TYPE_DECOMP_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GCG_Decomp GCG_DECOMP; /**< decomposition structure */

/** type of the decomposition */
enum Dectype
{
   GCG_DECTYPE_UNKNOWN   = 0,                /**< unknown structure (used for initialization) */
   GCG_DECTYPE_ARROWHEAD = 1,                /**< arrowhead structure (linking variables and constraints) */
   GCG_DECTYPE_STAIRCASE = 2,                /**< staircase structure (linking variables between consecutive blocks) */
   GCG_DECTYPE_DIAGONAL  = 3,                /**< block diagonal structure (no linking variables and constraints) */
   GCG_DECTYPE_BORDERED  = 4                 /**< bordered block diagonal structure (linking constraints only) */
};

typedef enum Dectype GCG_DECTYPE; /**< decomposition type */

/** the decomposition mode */
enum Decmode
{
   GCG_DECMODE_DANTZIGWOLFE = 0,             /**< Datizig-Wolfe reformulation */
   GCG_DECMODE_BENDERS      = 1,             /**< Benders' decomposition */
   GCG_DECMODE_ORIGINAL     = 2,             /**< the original problem will be solved without decomposition */
   GCG_DECMODE_AUTO         = 3,             /**< the best of either Dantzig-Wolfe or Benders' will be applied */
   GCG_DECMODE_UNKNOWN      = 4              /**< the mode can not be determined from the given information */
};

typedef enum Decmode GCG_DECMODE; /**< decomposition mode */

#ifdef __cplusplus
}
#endif

#endif
