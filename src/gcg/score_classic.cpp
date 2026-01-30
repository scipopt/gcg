/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file    score_classic.cpp
 * @brief   classic score
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include "gcg/score_classic.h"

/* score properties */
#define SCORE_NAME                "classic"
#define SCORE_SHORTNAME           "classi"
#define SCORE_DESC                "classical score"


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
#define scoreFreeClassic NULL

static
GCG_DECL_SCORECALC(scoreCalcClassic)
{
   int i;
   int j;
   int k;

   unsigned long matrixarea;
   unsigned long borderarea;
   SCIP_Real borderscore; /* score of the border */
   SCIP_Real densityscore; /* score of block densities */
   SCIP_Real linkingscore; /* score related to interlinking blocks */
   SCIP_Real totalscore; /* accumulated score */

   SCIP_Real varratio;
   int* nzblocks = NULL;
   int* nlinkvarsblocks = NULL;
   int* nvarsblocks = NULL;
   SCIP_Real* blockdensities = NULL;
   int* blocksizes = NULL;
   SCIP_Real density;

   SCIP_Real alphaborderarea;
   SCIP_Real alphalinking;
   SCIP_Real alphadensity;

   alphaborderarea = 0.6;
   alphalinking = 0.2;
   alphadensity = 0.2;

   SCIP* scip = GCGgetOrigprob(gcg);
   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(gcg, partialdecid);
   gcg::DETPROBDATA* detprobdata = partialdec->getDetprobdata();

   SCIP_CALL( SCIPallocBufferArray(scip, &nzblocks, partialdec->getNBlocks()) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlinkvarsblocks, partialdec->getNBlocks()) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blockdensities, partialdec->getNBlocks()) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocksizes, partialdec->getNBlocks()) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nvarsblocks, partialdec->getNBlocks()) );

   /*
    * 3 Scores
    *
    * - Area percentage (min)
    * - block density (max)
    * - \pi_b {v_b|v_b is linking}/#vb (min)
    */

   /* calculate slave sizes, nonzeros and linkingvars */
   for( i = 0; i < partialdec->getNBlocks(); ++ i )
   {
      int ncurconss;
      int nvarsblock;
      SCIP_Bool *ishandled = NULL;

      SCIP_CALL( SCIPallocBufferArray(scip, &ishandled, partialdec->getNVars()) );
      nvarsblock = 0;
      nzblocks[i] = 0;
      nlinkvarsblocks[i] = 0;

      for( j = 0; j < partialdec->getNVars(); ++ j )
      {
         ishandled[j] = FALSE;
      }
      ncurconss = partialdec->getNConssForBlock(i);

      for( j = 0; j < ncurconss; ++ j )
      {
         int cons = partialdec->getConssForBlock(i)[j];
         int ncurvars;
         ncurvars = detprobdata->getNVarsForCons(cons);
         for( k = 0; k < ncurvars; ++ k )
         {
            int var = detprobdata->getVarsForCons(cons)[k];
            int block = -3;
            if( partialdec->isVarBlockvarOfBlock(var, i) )
               block = i + 1;
            else if( partialdec->isVarLinkingvar(var) || partialdec->isVarStairlinkingvar(var) )
               block = partialdec->getNBlocks() + 2;
            else if( partialdec->isVarMastervar(var) )
               block = partialdec->getNBlocks() + 1;

            ++(nzblocks[i]);

            if( block == partialdec->getNBlocks() + 1 && ishandled[var] == FALSE )
            {
               ++(nlinkvarsblocks[i]);
            }
            ishandled[var] = TRUE;
         }
      }

      for( j = 0; j < partialdec->getNVars(); ++ j )
      {
         if( ishandled[j] )
         {
            ++nvarsblock;
         }
      }

      blocksizes[i] = nvarsblock * ncurconss;
      nvarsblocks[i] = nvarsblock;
      if( blocksizes[i] > 0 )
      {
         blockdensities[i] = 1.0 * nzblocks[i] / blocksizes[i];
      }
      else
      {
         blockdensities[i] = 0.0;
      }

      assert( blockdensities[i] >= 0 && blockdensities[i] <= 1.0 );
      SCIPfreeBufferArray( scip, &ishandled );
   }

   borderarea = ((unsigned long) partialdec->getNMasterconss() * partialdec->getNVars() ) + (((unsigned long) partialdec->getNLinkingvars() + partialdec->getNMastervars() + partialdec->getNTotalStairlinkingvars())) * (partialdec->getNConss() - partialdec->getNMasterconss());

   matrixarea = ((unsigned long) partialdec->getNVars()) * ((unsigned long) partialdec->getNConss());

   density = 1E20;
   varratio = 1.0;
   linkingscore = 1.;
   borderscore =  1.;
   densityscore = 1.;

   for( i = 0; i < partialdec->getNBlocks(); ++ i )
   {
      density = MIN(density, blockdensities[i]);

      if( (partialdec->getNLinkingvars() + partialdec->getNMastervars() + partialdec->getNTotalStairlinkingvars()) > 0 )
      {
         varratio *= 1.0 * nlinkvarsblocks[i] / (partialdec->getNLinkingvars() + partialdec->getNMastervars() + partialdec->getNTotalStairlinkingvars());
      }
      else
      {
         varratio = 0.;
      }
   }
   linkingscore = 0.5 + 0.5 * varratio;

   densityscore =  1. - density;

   borderscore = 1.0 * ( borderarea ) / matrixarea;

   totalscore = 1. - (alphaborderarea * borderscore + alphalinking * linkingscore + alphadensity * densityscore);

   if(totalscore > 1)
      totalscore = 1;
   if(totalscore < 0)
      totalscore = 0;

   *scorevalue = totalscore;

   SCIPfreeBufferArray( scip, & nzblocks );
   SCIPfreeBufferArray(  scip, & nlinkvarsblocks) ;
   SCIPfreeBufferArray(  scip, & blockdensities);
   SCIPfreeBufferArray(  scip, & blocksizes);
   SCIPfreeBufferArray(  scip, & nvarsblocks);

   return SCIP_OKAY;
}




/*
 * score specific interface methods
 */

/** creates the classic score and includes it in SCIP */
SCIP_RETCODE GCGincludeScoreClassic(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_SCOREDATA* scoredata = NULL;

   SCIP_CALL(
      GCGincludeScore(gcg, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata,
         scoreFreeClassic, scoreCalcClassic) );

   return SCIP_OKAY;
}
