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

/**@file   score_xyz.cpp
 * @ingroup DEFPLUGINS_SCORE
 * @brief  xyz score (put your description here)
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_partialdecomp.h"
#include "class_detprobdata.h"
#include "cons_decomp.h"
#include "cons_decomp.hpp"
#include "score_xyz.h"


/* score properties */
#define SCORE_NAME             "xyz score"                           /**< name of score */
#define SCORE_SHORTNAME        "xyz"                                 /**< shortname of score*/
#define SCORE_DESC             "score template"                      /**< short description of score */


/*
 * Data structures
 */

/* TODO: (optional) fill in the necessary score data */

/** data for xyz score */
struct DEC_ScoreData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of score
 */

/* TODO: Implement all necessary score methods. The methods with an #if SCIP_DISABLED_CODE ... #else 
 * #define ... are optional 
 */

/** destructor of classifier to free user data (called when GCG is exiting) */
#ifdef SCIP_DISABLED_CODE
static
DEC_DECL_SCOREFREE(scoreFreeXyz)
{  /*lint --e{715}*/

   SCIPerrorMessage("Free function of score <%s> not implemented!\n", DEC_SCORENAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define scoreFreeXyz NULL
#endif

/** calculate score method of score */
static
GCG_DECL_SCORECALC(scoreCalcXyz)
{  /*lint --e{715}*/

   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);

   SCIPerrorMessage("method of xyz score not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}




/*
 * score specific interface methods
 */

/** creates the handler for XYZ score and includes it in SCIP */
SCIP_RETCODE GCGincludeScoreXyz(
   SCIP*                 scip                /**< SCIP data structure */
   ) 
{
   GCG_SCOREDATA* scoredata;

   /**@todo create xyz score data here*/
   scoredata = NULL;

   SCIP_CALL( GCGincludeScore(scip, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata, scoreFreeXyz, scoreCalcXyz) );

   return SCIP_OKAY;
}
