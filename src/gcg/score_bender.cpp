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

/**@file    score_bender.cpp
 * @brief   benders score
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include "gcg/score_bender.h"

/* score properties */
#define SCORE_NAME                "experimental benders score"
#define SCORE_SHORTNAME           "bender"
#define SCORE_DESC                "experimental score to evaluate benders decompositions"


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
 * classifier callback methods
 */

/** destructor of classifier to free user data (called when GCG is exiting) */
#define scoreFreeBender NULL

static
GCG_DECL_SCORECALC(scoreCalcBender)
{
   SCIP_Real benderareascore = 0;

   unsigned long nrelevantconss;
   unsigned long nrelevantvars;
   unsigned long nrelevantconss2;
   unsigned long nrelevantvars2;
   unsigned long badblockvararea;
   long benderborderarea;
   unsigned long totalarea;

   nrelevantconss = 0;
   nrelevantvars = 0;
   nrelevantconss2 = 0;
   nrelevantvars2 = 0;
   badblockvararea = 0;

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(gcg, partialdecid);

   gcg::DETPROBDATA* detprobdata = partialdec->getDetprobdata();

   /* calc bender area score  (1 - fraction of white area in master constraints to complete area) */
   for( int  c = 0; c < partialdec->getNMasterconss(); ++c )
   {
      bool relevant = true;
      int cons = partialdec->getMasterconss()[c];
      for( int v = 0; v < detprobdata->getNVarsForCons(cons); ++v )
      {
         int var = detprobdata->getVarsForCons(cons)[v];
         if ( partialdec->isVarOpenvar(var) || partialdec->isVarMastervar(var) || partialdec->isVarLinkingvar(var) )
         {
            relevant = false;
            break;
         }

      }
      if( relevant )
         ++nrelevantconss;
   }

   for( int b = 0; b < partialdec->getNBlocks(); ++b )
   {
      for(int v = 0; v < partialdec->getNVarsForBlock(b); ++v )
      {
         bool relevant = true;
         int var = partialdec->getVarsForBlock(b)[v];
         for( int c = 0; c < detprobdata->getNConssForVar(var); ++c )
         {
            int cons = detprobdata->getConssForVar(var)[c];
            if( partialdec->isConsMastercons(cons) || partialdec->isConsOpencons(cons)  )
            {
               relevant  = false;
               for( int b2 = 0; b2 < partialdec->getNBlocks(); ++b2 )
               {
                  if( b2 != b )
                     badblockvararea += partialdec->getNConssForBlock(b2);
               }
               break;
            }
         }
         if( relevant )
            ++nrelevantvars;
      }
   }

   for( int  v = 0; v < partialdec->getNLinkingvars(); ++v )
   {
      bool relevant = true;
      int var = partialdec->getLinkingvars()[v];
      for( int c = 0; c < detprobdata->getNConssForVar(var); ++c )
      {
         int cons = detprobdata->getConssForVar(var)[c];
         if ( partialdec->isConsOpencons(cons) || partialdec->isConsMastercons(cons) )
         {
            relevant = false;
            break;
         }

      }
      if( relevant )
         ++nrelevantvars2;
   }

   for( int b = 0; b < partialdec->getNBlocks(); ++b )
   {
      for(int c = 0; c < partialdec->getNConssForBlock(b); ++c )
      {
         bool relevant = true;
         int cons = partialdec->getConssForBlock(b)[c];
         for( int v = 0; v < detprobdata->getNVarsForCons(cons); ++v )
         {
            int var = detprobdata->getVarsForCons(cons)[v];
            if( partialdec->isVarLinkingvar(var) || partialdec->isVarOpenvar(var)  )
            {
               relevant  = false;
               break;
            }
         }
         if( relevant )
            ++nrelevantconss2;
      }
   }

   benderborderarea = ( nrelevantconss * nrelevantvars  ) + ( nrelevantconss2 * nrelevantvars2  ) - badblockvararea;
   totalarea = ( (unsigned long) partialdec->getNConss() * (unsigned long) partialdec->getNVars() );
   benderareascore =  ( SCIP_Real) benderborderarea / totalarea;

   /* get block & border area score (note: these calculations have their own clock) */
   SCIP_Real borderareascore;

   SCIP_Real blockareascore = partialdec->calcBlockAreaScore();
   borderareascore = partialdec->getScore(GCGconshdlrDecompFindScore(gcg, "border area"));

   *scorevalue = blockareascore + benderareascore + borderareascore - 1.;

   if( *scorevalue < 0. )
      *scorevalue = 0.;

   return SCIP_OKAY;
}

/*
 * score specific interface methods
 */

/** creates the bender score and includes it in SCIP */
SCIP_RETCODE GCGincludeScoreBender(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_SCOREDATA* scoredata = NULL;

   SCIP_CALL(
      GCGincludeScore(gcg, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata,
         scoreFreeBender, scoreCalcBender) );

   return SCIP_OKAY;
}
