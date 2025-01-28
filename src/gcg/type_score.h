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

/**@file   type_score.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for score
 * @author Jurgen Lentz
 */

#ifndef __GCG_TYPE_SCORE_H__
#define __GCG_TYPE_SCORE_H__

#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GCG_Score GCG_SCORE;               /**< score data structure */
typedef struct GCG_ScoreData GCG_SCOREDATA;       /**< score specific data */


/** destructor of score to free score data (called when GCG is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - score           : score data structure
 */
#define GCG_DECL_SCOREFREE(x) SCIP_RETCODE x (SCIP* scip, GCG_SCORE* score)

/**
 * calculates the score value of partialdecid and stores it in scorevalue
 *
 * input:
 *  - scip                 : SCIP data structure
 *  - score                : score data structure
 *  - partialdecid         : id of partialdec to calculate its score
 *  - scorevalue           : value of the score
 */
#define GCG_DECL_SCORECALC(x) SCIP_RETCODE x (SCIP* scip, GCG_SCORE* score, int partialdecid, SCIP_Real* scorevalue)

#ifdef __cplusplus
}
#endif

#endif
