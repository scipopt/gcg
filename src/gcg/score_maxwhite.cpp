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

/**@file    score_maxwhite.cpp
 * @brief   max white score
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_partialdecomp.h"
#include "class_detprobdata.h"
#include "cons_decomp.h"
#include "cons_decomp.hpp"
#include "score_maxwhite.h"

/* score properties */
#define SCORE_NAME                "max white"
#define SCORE_SHORTNAME           "maxwhi"
#define SCORE_DESC                "maximum white area score (white area is nonblock and nonborder area)"


/*
 * Data structures
 */
struct GCG_ScoreData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * score callback methods
 */

/** destructor of score to free user data (called when GCG is exiting) */
#define scoreFreeMaxwhite NULL

static
GCG_DECL_SCORECALC(scoreCalcMaxwhite)
{
   SCIP_Real borderareascore;

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);

   SCIP_Real blockareascore = partialdec->calcBlockAreaScore(scip);
   borderareascore = partialdec->getScore(GCGconshdlrDecompFindScore(scip, "border area"));

   SCIP_Real maxwhitescore = blockareascore + borderareascore - 1.;

   if( maxwhitescore < 0. )
     maxwhitescore = 0.;

   *scorevalue = maxwhitescore;

   return SCIP_OKAY;
}




/*
 * score specific interface methods
 */

/** creates the max white score and includes it in SCIP */
SCIP_RETCODE GCGincludeScoreMaxwhite(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   GCG_SCOREDATA* scoredata = NULL;

   SCIP_CALL( 
      GCGincludeScore(scip, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata, 
         scoreFreeMaxwhite, scoreCalcMaxwhite) );

   return SCIP_OKAY;
}
