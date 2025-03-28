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

/**@file    score_border.cpp
 * @brief   border area score
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include "gcg/score_border.h"

/* score properties */
#define SCORE_NAME                "border area"
#define SCORE_SHORTNAME           "border"
#define SCORE_DESC                "minimum border score (i.e. minimizes fraction of border area score)"


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
#define scoreFreeBorder NULL

static
GCG_DECL_SCORECALC(scoreCalcBorder)
{
   unsigned long matrixarea;
   unsigned long borderarea;

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(gcg, partialdecid);

   matrixarea = (unsigned long) partialdec->getNVars() * (unsigned long)partialdec->getNConss();
   borderarea = 0;

   borderarea += (unsigned long) ( partialdec->getNLinkingvars() + partialdec->getNTotalStairlinkingvars() ) * (unsigned long) partialdec->getNConss();
   borderarea += (unsigned long) partialdec->getNMasterconss() * ( (unsigned long) partialdec->getNVars() - ( partialdec->getNLinkingvars() + partialdec->getNTotalStairlinkingvars() ) ) ;

   *scorevalue = 1. - (matrixarea == 0 ? 0 : ( (SCIP_Real) borderarea / (SCIP_Real) matrixarea ));

   return SCIP_OKAY;
}




/*
 * score specific interface methods
 */

/** creates the border score and includes it in SCIP */
SCIP_RETCODE GCGincludeScoreBorder(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_SCOREDATA* scoredata = NULL;

   SCIP_CALL( 
      GCGincludeScore(gcg, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata, 
         scoreFreeBorder, scoreCalcBorder) );

   return SCIP_OKAY;
}
