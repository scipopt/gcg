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

/**@file   type_score.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for score
 * @author Jurgen Lentz
 */

#ifndef __GCG_TYPE_SCORE_H__
#define __GCG_TYPE_SCORE_H__

#include "scip/def.h"
#include "type_gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GCG_Score GCG_SCORE;               /**< score data structure */
typedef struct GCG_ScoreData GCG_SCOREDATA;       /**< score specific data */


/** destructor of score to free score data (called when GCG is exiting)
 *
 *  input:
 *  - gcg             : GCG main data structure
 *  - score           : score data structure
 */
#define GCG_DECL_SCOREFREE(x) SCIP_RETCODE x (GCG* gcg, GCG_SCORE* score)

/**
 * calculates the score value of partialdecid and stores it in scorevalue
 *
 * input:
 *  - gcg                  : GCG data structure
 *  - score                : score data structure
 *  - partialdecid         : id of partialdec to calculate its score
 *  - scorevalue           : value of the score
 */
#define GCG_DECL_SCORECALC(x) SCIP_RETCODE x (GCG* gcg, GCG_SCORE* score, int partialdecid, SCIP_Real* scorevalue)

#ifdef __cplusplus
}
#endif

#endif
