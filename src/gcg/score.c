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

/**@file   score.c
 * @ingroup DECOMP
 * @brief  interface for score
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "pub_score.h"
#include "gcg.h"
#include "scip/scip.h"
#include "struct_score.h"
#include "cons_decomp.h"

#include <assert.h>


/**
 * @brief creates a score and includes it in GCG
 * @returns scip return code
 */
SCIP_RETCODE GCGincludeScore(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of score */
   const char*           shortname,          /**< shortname of score */
   const char*           description,        /**< description of score */
   GCG_SCOREDATA*        scoredata,          /**< score data */
   GCG_DECL_SCOREFREE    ((*scorefree)),     /**< destructor of score */
   GCG_DECL_SCORECALC    ((*scorecalc))      /**< score calculation method of score */
   )
{
   /* check whether score is already present */
   if( GCGfindScore(scip, name) != NULL || GCGfindScoreByShortname(scip, shortname) != NULL )
   {
      SCIPerrorMessage("Score <%s> is already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( GCGconshdlrDecompIncludeScore(scip, name, shortname, description, scoredata, scorefree, scorecalc) );

   return SCIP_OKAY;
}

/**
 * @brief searches for the score with the given name and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given name is not found
 */
GCG_SCORE* GCGfindScore(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of score */
   )
{
   assert(scip != NULL);
   assert(name != NULL);

   return GCGconshdlrDecompFindScore(scip, name);
}

/**
 * @brief searches for the score with the given shortname and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given shortname is not found
 */
GCG_SCORE* GCGfindScoreByShortname(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           shortname           /**< shortname of score */
   )
{
   assert(scip != NULL);
   assert(shortname != NULL);

   return GCGconshdlrDecompFindScoreByShortname(scip, shortname);
}

/** returns the array of currently available scores */
GCG_SCORE** GCGgetScores(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return GCGconshdlrDecompGetScores(scip);
}

/** returns the number of currently available scores */
int GCGgetNScores(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return GCGconshdlrDecompGetNScores(scip);
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
