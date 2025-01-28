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

/**@file   struct_score.h
 * @ingroup INTERNALAPI-GCG
 * @brief  data structures for score
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_STRUCT_SCORE_H__
#define __GCG_STRUCT_SCORE_H__

#include "scip/def.h"
#include "type_score.h"

#ifdef __cplusplus
extern "C" {
#endif

/** score */
struct GCG_Score
{
   const char*           name;               /**< name of the score */
   const char*           shortname;          /**< shortname of the score */
   const char*           description;        /**< description of the score */
   GCG_DECL_SCOREFREE    ((*scorefree));     /**< destructor of score */
   GCG_DECL_SCORECALC    ((*scorecalc));     /**< calculate method of score */
   GCG_SCOREDATA*        scoredata;          /**< custom data structure of the score */
};

#ifdef __cplusplus
}
#endif

#endif
