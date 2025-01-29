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

/**@file    score_bender.cpp
 * @brief   benders score
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_partialdecomp.h"
#include "class_detprobdata.h"
#include "cons_decomp.h"
#include "cons_decomp.hpp"
#include "score_bender.h"

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

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);

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

   SCIP_Real blockareascore = partialdec->calcBlockAreaScore(scip);
   borderareascore = partialdec->getScore(GCGconshdlrDecompFindScore(scip, "border area"));

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
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   GCG_SCOREDATA* scoredata = NULL;

   SCIP_CALL(
      GCGincludeScore(scip, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata,
         scoreFreeBender, scoreCalcBender) );

   return SCIP_OKAY;
}
