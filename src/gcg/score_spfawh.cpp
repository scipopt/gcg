/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
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

/**@file    score_spfawh.cpp
 * @brief   setpartitioning maximum foreseeing white area score with aggregation information
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include "gcg/score_spfawh.h"


/* score properties */
#define SCORE_NAME                "ppc-max-white with aggregation info"
#define SCORE_SHORTNAME           "spfawh"
#define SCORE_DESC                "setpartitioning maximum foreseeing white area score with aggregation information (convex combination of maximum foreseeing white area score and rewarding if a master contains only setppc and cardinality constraints)"


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
#define scoreFreeSpfawh NULL

static
GCG_DECL_SCORECALC(scoreCalcSpfawh)
{
   SCIP_Real maxforeseeingwhitescoreagg = 0;

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(gcg, partialdecid);

   maxforeseeingwhitescoreagg = partialdec->getScore(GCGconshdlrDecompFindScore(gcg, "max foreseeing white with aggregation info"));

   if( partialdec->hasSetppccardMaster() && !partialdec->isTrivial() && partialdec->getNBlocks() > 1 )
   {
      *scorevalue = 0.5 * maxforeseeingwhitescoreagg + 0.5;
   }
   else
   {
      *scorevalue = 0.5 * maxforeseeingwhitescoreagg;
   }

   return SCIP_OKAY;
}




/*
 * score specific interface methods
 */

/** creates the setpartitioning maximum foreseeing white area score with aggregation information score and includes it in SCIP */
SCIP_RETCODE GCGincludeScoreSpfawh(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_SCOREDATA* scoredata = NULL;

   SCIP_CALL(
      GCGincludeScore(gcg, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata,
         scoreFreeSpfawh, scoreCalcSpfawh) );

   return SCIP_OKAY;
}
