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

/**@file   gcg_score.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for score plugins
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg_score.h"
#include "cons_decomp.h"

/**
 * @brief includes one score
 * @returns scip return code
 */
SCIP_RETCODE GCGincludeScore(
   SCIP*                 scip,
   const char*           name,               /**< name of the score*/
   const char*           shortname,          /**< shortname of the score */
   const char*           description,        /**< description of the score */
   DEC_SCOREDATA*        scoredata,
   DEC_DECL_SCOREFREE    ((*scorefree)),
   DEC_DECL_SCORECALC    ((*scorecalc))
   )
{
    /* check whether score is already present */
    if( GCGfindScore(scip, name) != NULL )
    {
       SCIPerrorMessage("score <%s> already included.\n", name);
       return SCIP_INVALIDDATA;
    }

    SCIP_CALL( DECincludeScore(scip, name, shortname, description, scoredata, scorefree, scorecalc) );

    return SCIP_OKAY;
}

/**
 * @brief searches for the score with the given name and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given name is not found
 */
DEC_SCORE* GCGfindScore(
   SCIP*                 scip,
   const char*           name
   )
{
    assert(scip != NULL);
    assert(name != NULL);

    return DECfindScore(scip, name);
}

/**
 * @brief searches for the score with the given shortname and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given shortname is not found
 */
DEC_SCORE* GCGfindScoreByShortname(
   SCIP*                 scip,
   const char*           shortname
   )
{
    assert(scip != NULL);
    assert(shortname != NULL);

    return DECfindScoreByShortname(scip, shortname);
}

/** returns the array of currently available scores */
DEC_SCORE** GCGgetScores(
   SCIP*                 scip
   )
{
    assert(scip != NULL);

    return GCGconshdlrDecompGetScores(scip);
}

/** returns the number of currently available scores */
int GCGgetNScores(
   SCIP*                 scip
   )
{
    assert(scip != NULL);

    return GCGconshdlrDecompGetNScores(scip);
}
