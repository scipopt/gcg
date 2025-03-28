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
