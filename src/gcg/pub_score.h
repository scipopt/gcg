/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   pub_score.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for score
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_PUB_SCORE_H__
#define __GCG_PUB_SCORE_H__


#include "scip/scip.h"
#include "gcg/def.h"
#include "gcg/type_score.h"
#include "gcg/type_gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicScoreMethods
 *
 * @{
 */

/**
 * @brief creates a score and includes it in GCG
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGincludeScore(
   GCG*                  gcg,                /**< GCG data structure */
   const char*           name,               /**< name of score */
   const char*           shortname,          /**< shortname of score */
   const char*           description,        /**< description of score */
   GCG_SCOREDATA*        scoredata,          /**< score data */
   GCG_DECL_SCOREFREE    ((*scorefree)),     /**< destructor of score */
   GCG_DECL_SCORECALC    ((*scorecalc))      /**< score calculation method of score */
   );

/**
 * @brief searches for the score with the given name and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given name is not found
 */
GCG_EXPORT
GCG_SCORE* GCGfindScore(
   GCG*                  gcg,                /**< GCG data structure */
   const char*           name                /**< name of score */
   );

/**
 * @brief searches for the score with the given shortname and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given shortname is not found
 */
GCG_EXPORT
GCG_SCORE* GCGfindScoreByShortname(
   GCG*                  gcg,                /**< GCG data structure */
   const char*           shortname           /**< shortname of score */
   );

/** returns the array of currently available scores */
GCG_EXPORT
GCG_SCORE** GCGgetScores(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the number of currently available scores */
GCG_EXPORT
int GCGgetNScores(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** gets user data of score */
GCG_EXPORT
GCG_SCOREDATA* GCGscoreGetData(
   GCG_SCORE*            score               /**< score */
   );

/** sets user data of score; user has to free old data in advance! */
GCG_EXPORT
void GCGscoreSetData(
   GCG_SCORE*            score,              /**< score */
   GCG_SCOREDATA*        scoredata           /**< new score user data */
   );

/** gets name of score */
GCG_EXPORT
const char* GCGscoreGetName(
   GCG_SCORE*            score               /**< score */
   );

/** gets shortname of score */
GCG_EXPORT
const char* GCGscoreGetShortname(
   GCG_SCORE*            score               /**< score */
   );

/** gets description of score */
GCG_EXPORT
const char* GCGscoreGetDesc(
   GCG_SCORE*            score               /**< score */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
