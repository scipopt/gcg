/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2022 Operations Research, RWTH Aachen University       */
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

#ifndef GCG_TYPE_SCORE_H__
#define GCG_TYPE_SCORE_H__


typedef struct DEC_Score DEC_SCORE;
typedef struct DEC_ScoreData DEC_SCOREDATA;       /**< score specific data */


/** destructor of classifier to free score data (called when GCG is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - score           : score data structure
 */
#define DEC_DECL_SCOREFREE(x) SCIP_RETCODE x (SCIP* scip, DEC_SCORE* score)

/**
 * Tries to sco with data of the according detprobdata and store the classification in the detprobdata
 *
 * input:
 *  - scip                 : SCIP data structure
 *  - score                : score data structure
 *  - partialdecid         : id of partialdec to calculate its score
 *  - scorevalue           : value of the score
 */
#define DEC_DECL_SCORECALC(x) SCIP_RETCODE x (SCIP* scip, DEC_SCORE* score, int partialdecid, SCIP_Real* scorevalue)


#endif //GCG_TYPE_CONSCLASSIFIER_H__
