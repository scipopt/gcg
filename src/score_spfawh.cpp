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

/**@file    score_spfawh.cpp
 * @ingroup SCORES
 * @brief   setpartitioning maximum foreseeing white area score with aggregation information
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_partialdecomp.h"
#include "class_detprobdata.h"
#include "cons_decomp.h"
#include "cons_decomp.hpp"
#include "score_spfawh.h"


/* score properties */
#define SCORE_NAME                "ppc-max-white with aggregation info"
#define SCORE_SHORTNAME           "spfawh"
#define SCORE_DESC                "setpartitioning maximum foreseeing white area score with aggregation information (convex combination of maximum foreseeing white area score and rewarding if a master contains only setppc and cardinality constraints)"


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
 * @brief calculates the maxforeseeingwhiteagg score of a partialdec
 *
 * maximum foreseeing white area score with respect to aggregatable blocks
 * (i.e. maximize fraction of white area score considering problem with copied linking variables
 * and corresponding master constraints;
 * white area is nonblock and nonborder area, stairlinking variables count as linking)
 * @return scip return code
 */
SCIP_RETCODE GCGconshdlrDecompCalcMaxForeseeingWhiteAggScore(
   SCIP* scip,
   int partialdecid,
   SCIP_Real* score
   )
{
   unsigned long sumblockshittinglinkingvar;
   unsigned long sumlinkingvarshittingblock;
   unsigned long newheight;
   unsigned long newwidth;
   unsigned long newmasterarea;
   unsigned long newblockareaagg;

   SCIP_CLOCK* clock;
   SCIP_CALL_ABORT( SCIPcreateClock(scip, &clock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, clock) );

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);

   std::vector<int> nlinkingvarsforblock(partialdec->getNBlocks(), 0);
   std::vector<int> nblocksforlinkingvar(partialdec->getNLinkingvars() + partialdec->getNTotalStairlinkingvars(), 0);

   gcg::DETPROBDATA* detprobdata = partialdec->getDetprobdata();

   partialdec->calcAggregationInformation(false);

   for( int lv = 0; lv < partialdec->getNLinkingvars(); ++lv )
   {
      int linkingvarid = partialdec->getLinkingvars()[lv];

      for( int b = 0; b < partialdec->getNBlocks(); ++b )
      {
         for ( int blc = 0; blc < partialdec->getNConssForBlock(b); ++blc )
         {
            int blockcons = partialdec->getConssForBlock(b)[blc];
            if( !SCIPisZero( scip, detprobdata->getVal(blockcons, linkingvarid) ) )
            {
               /* linking var hits block */
               ++nlinkingvarsforblock[b];
               ++nblocksforlinkingvar[lv];
               break;
            }
         }
      }
   }

   for( int b = 0; b < partialdec->getNBlocks(); ++b)
   {
      for( int slv = 0; slv < partialdec->getNStairlinkingvars(b); ++slv )
      {
         ++nlinkingvarsforblock[b];
         ++nlinkingvarsforblock[b+1];
         ++nblocksforlinkingvar[partialdec->getNLinkingvars() + slv];
         ++nblocksforlinkingvar[partialdec->getNLinkingvars() + slv];
      }
   }

   sumblockshittinglinkingvar = 0;
   sumlinkingvarshittingblock = 0;
   for( int b = 0; b < partialdec->getNBlocks(); ++b )
   {
      sumlinkingvarshittingblock += nlinkingvarsforblock[b];
   }
   for( int lv = 0; lv < partialdec->getNLinkingvars(); ++lv )
   {
      sumblockshittinglinkingvar += nblocksforlinkingvar[lv];
   }

   for( int slv = 0; slv < partialdec->getNTotalStairlinkingvars(); ++slv )
   {
      sumblockshittinglinkingvar += nblocksforlinkingvar[partialdec->getNLinkingvars() + slv];
   }

   newheight = partialdec->getNConss() + sumblockshittinglinkingvar;
   newwidth = partialdec->getNVars() + sumlinkingvarshittingblock;

   newmasterarea = ( partialdec->getNMasterconss() + sumblockshittinglinkingvar) * ( partialdec->getNVars() + sumlinkingvarshittingblock );
   newblockareaagg = 0;

   for( int br = 0; br < partialdec->getNReps(); ++br )
   {
      newblockareaagg += partialdec->getNConssForBlock( partialdec->getBlocksForRep(br)[0] ) * ( partialdec->getNVarsForBlock( partialdec->getBlocksForRep(br)[0] ) + nlinkingvarsforblock[partialdec->getBlocksForRep(br)[0]] );
   }

   SCIP_Real maxforeseeingwhitescoreagg = ((SCIP_Real ) newblockareaagg + (SCIP_Real) newmasterarea) / (SCIP_Real) newwidth;
   maxforeseeingwhitescoreagg =  maxforeseeingwhitescoreagg / (SCIP_Real) newheight ;

   maxforeseeingwhitescoreagg = 1. - maxforeseeingwhitescoreagg;
   *score = maxforeseeingwhitescoreagg;

   partialdec->setMaxForWhiteAggScore(*score);

   SCIP_CALL_ABORT(SCIPstopClock(scip, clock) );
   GCGconshdlrDecompAddScoreTime(scip, SCIPgetClockTime(scip, clock));
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &clock) );

   return SCIP_OKAY;
}




/*
 * score callback methods
 */

/** destructor of score to free user data (called when GCG is exiting) */
#define scoreFreeSpfawh NULL

static
DEC_DECL_SCORECALC(scoreCalcSpfawh)
{
   SCIP_Real maxforeseeingwhitescoreagg = 0;

   SCIP_CLOCK* clock;
   SCIP_CALL_ABORT( SCIPcreateClock(scip, &clock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, clock) );

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);

   SCIP_CALL_ABORT(SCIPstopClock(scip, clock) );
   GCGconshdlrDecompAddScoreTime(scip, SCIPgetClockTime(scip, clock));

   maxforeseeingwhitescoreagg = partialdec->getMaxForWhiteAggScore();
   if(maxforeseeingwhitescoreagg == -1)
      GCGconshdlrDecompCalcMaxForeseeingWhiteAggScore(scip, partialdecid, &maxforeseeingwhitescoreagg);

   SCIP_CALL_ABORT(SCIPresetClock( scip, clock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, clock) );

   if( partialdec->hasSetppccardMaster() && !partialdec->isTrivial() && partialdec->getNBlocks() > 1 )
   {
      *scorevalue = 0.5 * maxforeseeingwhitescoreagg + 0.5;
   }
   else
   {
      *scorevalue = 0.5 * maxforeseeingwhitescoreagg;
   }

   partialdec->setSetPartForWhiteAggScore(*scorevalue);

   SCIP_CALL_ABORT(SCIPstopClock(scip, clock) );
   GCGconshdlrDecompAddScoreTime(scip, SCIPgetClockTime(scip, clock));
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &clock) );

   return SCIP_OKAY;
}




/*
 * score specific interface methods
 */

/** creates the setpartitioning maximum foreseeing white area score with aggregation information score and includes it in SCIP */
SCIP_RETCODE GCGincludeScoreSpfawh(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   DEC_SCOREDATA* scoredata = NULL;

   SCIP_CALL(
      GCGincludeScore(scip, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata,
         scoreFreeSpfawh, scoreCalcSpfawh) );

   return SCIP_OKAY;
}
