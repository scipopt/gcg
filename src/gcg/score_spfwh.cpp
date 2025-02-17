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

/**@file    score_spfwh.cpp
 * @brief   setpartitioning maximum foreseeing white area score
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_partialdecomp.h"
#include "class_detprobdata.h"
#include "cons_decomp.h"
#include "cons_decomp.hpp"
#include "score_spfwh.h"
#include "struct_score.h"


/* score properties */
#define SCORE_NAME                "ppc-max-white"
#define SCORE_SHORTNAME           "spfwh"
#define SCORE_DESC                "setpartitioning maximum foreseeing white area score (convex combination of maximum foreseeing white area score and rewarding if master contains only setppc and cardinality constraints)"


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
#define scoreFreeSpfwh NULL

static
GCG_DECL_SCORECALC(scoreCalcSpfwh)
{
   SCIP_Real maxforeseeingwhitescore = 0;

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);

   maxforeseeingwhitescore = partialdec->getScore(GCGconshdlrDecompFindScore(scip, "max foreseeing white"));

   if( partialdec->hasSetppccardMaster() && !partialdec->isTrivial() && partialdec->getNBlocks() > 1 )
   {
      *scorevalue = 0.5 * maxforeseeingwhitescore + 0.5;
   }
   else
   {
      *scorevalue = 0.5 * maxforeseeingwhitescore;
   }

   return SCIP_OKAY;
}




/*
 * score specific interface methods
 */

/** creates the setpartitioning maximum foreseeing white area score and includes it in SCIP */
SCIP_RETCODE GCGincludeScoreSpfwh(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   GCG_SCOREDATA* scoredata = NULL;

   SCIP_CALL(
      GCGincludeScore(scip, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata,
         scoreFreeSpfwh, scoreCalcSpfwh) );

   return SCIP_OKAY;
}
