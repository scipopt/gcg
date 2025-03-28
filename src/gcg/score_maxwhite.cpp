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

/**@file    score_maxwhite.cpp
 * @brief   max white score
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include "gcg/score_maxwhite.h"

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

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(gcg, partialdecid);

   SCIP_Real blockareascore = partialdec->calcBlockAreaScore();
   borderareascore = partialdec->getScore(GCGconshdlrDecompFindScore(gcg, "border area"));

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
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_SCOREDATA* scoredata = NULL;

   SCIP_CALL( 
      GCGincludeScore(gcg, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata, 
         scoreFreeMaxwhite, scoreCalcMaxwhite) );

   return SCIP_OKAY;
}
