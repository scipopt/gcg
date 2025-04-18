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

/**@file    score_forswh.cpp
 * @brief   max foreseeing white score
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include "gcg/score_forswh.h"


/* score properties */
#define SCORE_NAME                "max foreseeing white"
#define SCORE_SHORTNAME           "forswh"
#define SCORE_DESC                "maximum foreseeing white area score (considering copied linking vars and their master conss; white area is nonblock and nonborder area)"


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
#define scoreFreeForswh NULL

static
GCG_DECL_SCORECALC(scoreCalcForswh)
{
   unsigned long sumblockshittinglinkingvar;
   unsigned long sumlinkingvarshittingblock;
   unsigned long newheight;
   unsigned long newwidth;
   unsigned long newmasterarea;
   unsigned long newblockarea;

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(gcg, partialdecid);
   gcg::DETPROBDATA* detprobdata = partialdec->getDetprobdata();

   std::vector<int> blockforconss(partialdec->getNConss(), -1);
   std::vector<int> nlinkingvarsforblock(partialdec->getNBlocks(), 0);
   std::vector<int> nblocksforlinkingvar(partialdec->getNLinkingvars() + partialdec->getNTotalStairlinkingvars(), 0);

   /* store for each cons to which block it is assigned (or -1 if border or unassigned) */
   for( int b = 0; b < partialdec->getNBlocks(); ++b )
   {
      std::vector<int>& blockconss = partialdec->getConssForBlock(b);
      for ( int blc = 0; blc < partialdec->getNConssForBlock(b); ++blc )
      {
         int blockcons = blockconss[blc];
         blockforconss[blockcons] = b;
      }
   }

   /* iterate linking vars and corresponding conss to recognize hit blocks */
   for( int lv = 0; lv < partialdec->getNLinkingvars(); ++lv )
   {
      int linkingvarid = partialdec->getLinkingvars()[lv];
      int nhittingconss = detprobdata->getNConssForVar(linkingvarid);
      std::vector<int>& hittingconss = detprobdata->getConssForVar(linkingvarid);

      std::vector<bool> hitblock(partialdec->getNBlocks(), false);

      /* find out which blocks the linking var is hitting */
      for ( int hittingcons = 0; hittingcons < nhittingconss; ++hittingcons )
      {
         int blockforcons = blockforconss[hittingconss[hittingcons]];
         if( blockforcons != -1 )
            hitblock[blockforcons] = true;
      }

      /* count for each block and each linking var how many linking vars or blocks, respectively, they hit */
      for( int b = 0; b < partialdec->getNBlocks(); ++b )
      {
         if ( hitblock[b] )
         {
            /* linking var hits block, so count it */
            ++nlinkingvarsforblock[b];
            ++nblocksforlinkingvar[lv];
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

   newmasterarea = ( (SCIP_Real) partialdec->getNMasterconss() + sumblockshittinglinkingvar) * ( (SCIP_Real) partialdec->getNVars() + sumlinkingvarshittingblock );
   newblockarea = 0;

   for( int b = 0; b < partialdec->getNBlocks(); ++b )
   {
      newblockarea += ((SCIP_Real) partialdec->getNConssForBlock(b) ) * ( (SCIP_Real) partialdec->getNVarsForBlock(b) + nlinkingvarsforblock[b] );
   }

   SCIP_Real maxforeseeingwhitescore = newwidth == 0 ? 1. : ((SCIP_Real ) newblockarea + (SCIP_Real) newmasterarea) / (SCIP_Real) newwidth;
   maxforeseeingwhitescore =  newheight == 0 ? 1. : maxforeseeingwhitescore / (SCIP_Real) newheight;

   maxforeseeingwhitescore = 1. - maxforeseeingwhitescore;
   assert(maxforeseeingwhitescore == SCIP_INVALID || (maxforeseeingwhitescore >= 0. && maxforeseeingwhitescore <= 1.));
   *scorevalue = maxforeseeingwhitescore;

   return SCIP_OKAY;
}




/*
 * score specific interface methods
 */

/** creates the max foreseeing white score and includes it in SCIP */
SCIP_RETCODE GCGincludeScoreForswh(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_SCOREDATA* scoredata = NULL;

   SCIP_CALL(
      GCGincludeScore(gcg, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata,
         scoreFreeForswh, scoreCalcForswh) );

   return SCIP_OKAY;
}
