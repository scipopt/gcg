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

/**@file   score.c
 * @ingroup DECOMP
 * @brief  interface for score
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/pub_score.h"
#include "gcg/gcg.h"
#include "scip/scip.h"
#include "gcg/struct_score.h"
#include "gcg/cons_decomp.h"

#include <assert.h>


/**
 * @brief creates a score and includes it in GCG
 * @returns scip return code
 */
SCIP_RETCODE GCGincludeScore(
   GCG*                  gcg,                /**< GCG data structure */
   const char*           name,               /**< name of score */
   const char*           shortname,          /**< shortname of score */
   const char*           description,        /**< description of score */
   GCG_SCOREDATA*        scoredata,          /**< score data */
   GCG_DECL_SCOREFREE    ((*scorefree)),     /**< destructor of score */
   GCG_DECL_SCORECALC    ((*scorecalc))      /**< score calculation method of score */
   )
{
   /* check whether score is already present */
   if( GCGfindScore(gcg, name) != NULL || GCGfindScoreByShortname(gcg, shortname) != NULL )
   {
      SCIPerrorMessage("Score <%s> is already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( GCGconshdlrDecompIncludeScore(gcg, name, shortname, description, scoredata, scorefree, scorecalc) );

   return SCIP_OKAY;
}

/**
 * @brief searches for the score with the given name and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given name is not found
 */
GCG_SCORE* GCGfindScore(
   GCG*                  gcg,                /**< GCG data structure */
   const char*           name                /**< name of score */
   )
{
   assert(gcg != NULL);
   assert(name != NULL);

   return GCGconshdlrDecompFindScore(gcg, name);
}

/**
 * @brief searches for the score with the given shortname and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given shortname is not found
 */
GCG_SCORE* GCGfindScoreByShortname(
   GCG*                  gcg,                /**< GCG data structure */
   const char*           shortname           /**< shortname of score */
   )
{
   assert(gcg != NULL);
   assert(shortname != NULL);

   return GCGconshdlrDecompFindScoreByShortname(gcg, shortname);
}

/** returns the array of currently available scores */
GCG_SCORE** GCGgetScores(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   assert(gcg != NULL);

   return GCGconshdlrDecompGetScores(gcg);
}

/** returns the number of currently available scores */
int GCGgetNScores(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   assert(gcg != NULL);

   return GCGconshdlrDecompGetNScores(gcg);
}

/** gets user data of score */
GCG_SCOREDATA* GCGscoreGetData(
   GCG_SCORE*            score               /**< score */
   )
{
   assert(score != NULL);

   return score->scoredata;
}

/** sets user data of score; user has to free old data in advance! */
void GCGscoreSetData(
   GCG_SCORE*            score,              /**< score */
   GCG_SCOREDATA*        scoredata           /**< new score user data */
   )
{
   assert(score != NULL);

   score->scoredata = scoredata;
}


/** gets name of score */
const char* GCGscoreGetName(
   GCG_SCORE*            score               /**< score */
   )
{
   assert(score != NULL);

   return score->name;
}

/** gets shortname of score */
const char* GCGscoreGetShortname(
   GCG_SCORE*            score               /**< score */
   )
{
   assert(score != NULL);

   return score->shortname;
}

/** gets description of score */
const char* GCGscoreGetDesc(
   GCG_SCORE*            score               /**< score */
   )
{
   assert(score != NULL);

   return score->description;
}
