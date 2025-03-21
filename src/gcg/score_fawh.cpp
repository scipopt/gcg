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

/**@file    score_fawh.cpp
 * @brief   maximum foreseeing white area score with aggregation info score
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include "gcg/score_fawh.h"


/* score properties */
#define SCORE_NAME                "max foreseeing white with aggregation info"
#define SCORE_SHORTNAME           "fawh"
#define SCORE_DESC                "maximum foreseeing white area score with aggregation info (considering copied linking vars and their master conss; white area is nonblock and nonborder area)"


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
#define scoreFreeFawh NULL

static
GCG_DECL_SCORECALC(scoreCalcFawh)
{
   unsigned long sumblockshittinglinkingvar;
   unsigned long sumlinkingvarshittingblock;
   unsigned long newheight;
   unsigned long newwidth;
   unsigned long newmasterarea;
   unsigned long newblockareaagg;

   SCIP* scip = GCGgetOrigprob(gcg);
   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(gcg, partialdecid);

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

   for( int br = 0; br < partialdec->getNEquivalenceClasses(); ++br )
   {
      newblockareaagg += partialdec->getNConssForBlock(partialdec->getBlocksForEqClass(br)[0])
            * ( partialdec->getNVarsForBlock(partialdec->getBlocksForEqClass(br)[0]) + nlinkingvarsforblock[partialdec->getBlocksForEqClass(br)[0]] );
   }

   SCIP_Real maxforeseeingwhitescoreagg = newwidth == 0 ? 1. : ((SCIP_Real ) newblockareaagg + (SCIP_Real) newmasterarea) / (SCIP_Real) newwidth;
   maxforeseeingwhitescoreagg =  newheight == 0 ? 1. : maxforeseeingwhitescoreagg / (SCIP_Real) newheight;

   maxforeseeingwhitescoreagg = 1. - maxforeseeingwhitescoreagg;
   assert(maxforeseeingwhitescoreagg == SCIP_INVALID || (maxforeseeingwhitescoreagg >= 0.0 && maxforeseeingwhitescoreagg <= 1.0));
   *scorevalue = maxforeseeingwhitescoreagg;

   return SCIP_OKAY;
}




/*
 * score specific interface methods
 */

/** creates the maximum foreseeing white area score with aggregation info score and includes it in SCIP */
SCIP_RETCODE GCGincludeScoreFawh(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_SCOREDATA* scoredata = NULL;

   SCIP_CALL(
      GCGincludeScore(gcg, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata,
         scoreFreeFawh, scoreCalcFawh) );

   return SCIP_OKAY;
}
