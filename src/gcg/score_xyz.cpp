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

/**@file   score_xyz.cpp
 * @ingroup DEFPLUGINS_SCORE
 * @brief  xyz score (put your description here)
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include "gcg/score_xyz.h"


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

   SCIP_CALL( GCGincludeScore(gcg, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata, scoreFreeXyz, scoreCalcXyz) );

   return SCIP_OKAY;
}
