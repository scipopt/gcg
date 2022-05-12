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

/**@file    score_maxwhite.cpp
 * @ingroup  SCORES
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
struct DEC_ScoreData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/**
 * @brief gets an intermediate score value for the blocks of a partialdec
 *
 * Used by several score calculations,
 * computed as (1 - fraction of block area to complete area)
 * 
 * @returns intermediate score value
 */
static
SCIP_Real calcBlockAreaScore(
   SCIP* scip,                /**< SCIP data structure */
   gcg::PARTIALDECOMP* partialdec  /**< compute for this partialdec */
   )
{
   unsigned long matrixarea;
   unsigned long blockarea;

   SCIP_CLOCK* clock;
   SCIP_CALL_ABORT( SCIPcreateClock( scip, &clock) );
   SCIP_CALL_ABORT( SCIPstartClock( scip, clock) );

   matrixarea = (unsigned long) partialdec->getNVars() * (unsigned long) partialdec->getNConss() ;
   blockarea = 0;

   for( int i = 0; i < partialdec->getNBlocks(); ++ i )
   {
      blockarea += (unsigned long) partialdec->getNConssForBlock(i) * ( (unsigned long) partialdec->getNVarsForBlock(i) );
   }

   SCIP_Real blockareascore = 1. - (matrixarea == 0 ? 0 : ( (SCIP_Real) blockarea / (SCIP_Real) matrixarea ));

   SCIP_CALL_ABORT( SCIPstopClock(scip, clock) );
   GCGconshdlrDecompAddScoreTime(scip, SCIPgetClockTime(scip, clock));
   SCIP_CALL_ABORT( SCIPfreeClock(scip, &clock) );

   return blockareascore;
}

/**
 * @brief calculates the border area score of a partialdec
 *
 * 1 - fraction of border area to complete area
 * @return scip return code
 */
static
SCIP_RETCODE GCGconshdlrDecompCalcBorderAreaScore(
   SCIP* scip,
   int partialdecid,
   SCIP_Real* score
   )
{
   unsigned long matrixarea;
   unsigned long borderarea;

   SCIP_CLOCK* clock;
   SCIP_CALL_ABORT( SCIPcreateClock( scip, &clock) );
   SCIP_CALL_ABORT( SCIPstartClock( scip, clock) );

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);

   matrixarea = (unsigned long) partialdec->getNVars() * (unsigned long)partialdec->getNConss();
   borderarea = 0;

   borderarea += (unsigned long) ( partialdec->getNLinkingvars() + partialdec->getNTotalStairlinkingvars() ) * (unsigned long) partialdec->getNConss();
   borderarea += (unsigned long) partialdec->getNMasterconss() * ( (unsigned long) partialdec->getNVars() - ( partialdec->getNLinkingvars() + partialdec->getNTotalStairlinkingvars() ) ) ;

   *score = 1. - (matrixarea == 0 ? 0 : ( (SCIP_Real) borderarea / (SCIP_Real) matrixarea ));

   partialdec->setBorderAreaScore(*score);

   SCIP_CALL_ABORT(SCIPstopClock( scip, clock) );
   GCGconshdlrDecompAddScoreTime(scip, SCIPgetClockTime( scip, clock));
   SCIP_CALL_ABORT(SCIPfreeClock( scip, &clock) );

   return SCIP_OKAY;
}




/*
 * score callback methods
 */

/** destructor of score to free user data (called when GCG is exiting) */
#define scoreFreeMaxwhite NULL

static
DEC_DECL_SCORECALC(scoreCalcMaxwhite)
{
   SCIP_Real borderareascore = 0;
   SCIP_CLOCK* clock;
   SCIP_CALL_ABORT( SCIPcreateClock(scip, &clock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, clock) );

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);

   SCIP_CALL_ABORT(SCIPstopClock(scip, clock) );
   GCGconshdlrDecompAddScoreTime(scip, SCIPgetClockTime(scip, clock));

   SCIP_Real blockareascore = calcBlockAreaScore(scip, partialdec);
   borderareascore = partialdec->getBorderAreaScore();
   if(borderareascore == -1)
      GCGconshdlrDecompCalcBorderAreaScore(scip, partialdecid, &borderareascore);

   SCIP_CALL_ABORT(SCIPresetClock( scip, clock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, clock) );

   SCIP_Real maxwhitescore = blockareascore + borderareascore - 1.;

   if( maxwhitescore < 0. )
     maxwhitescore = 0.;

   *scorevalue = maxwhitescore;
   partialdec->setMaxWhiteScore(*scorevalue);

   SCIP_CALL_ABORT(SCIPstopClock(scip, clock) );
   GCGconshdlrDecompAddScoreTime(scip, SCIPgetClockTime(scip, clock));
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &clock) );

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
   DEC_SCOREDATA* scoredata = NULL;

   SCIP_CALL( 
      GCGincludeScore(scip, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata, 
         scoreFreeMaxwhite, scoreCalcMaxwhite) );

   return SCIP_OKAY;
}
