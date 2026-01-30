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

/**@file   dec_generalmastersetcover.cpp
 * 
 * @brief  detector for set covering constraints
 * @author Michael Bastubbe
 *
 * This detector sets the following constraints to master:
 * - set covering constraints
 * - logical OR constraints
 * - constraints with infinity rhs and nonnegative lhs
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/dec_generalmastersetcover.h"
#include "gcg/cons_decomp.h"
#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/gcg.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "gcg/scip_misc.h"
#include "scip/clock.h"

#include <iostream>

/* constraint handler properties */
#define DEC_NAME                  "generalmastersetcover"       /**< name of detector */
#define DEC_DESC                  "detector generalmastersetcover" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          0     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  0     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               '?'         /**< display character of detector */
#define DEC_ENABLED               TRUE        /**< should the detection be enabled */
#define DEC_ENABLEDFINISHING      FALSE       /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated partialdec */

/*
 * Data structures
 */

/** @todo fill in the necessary detector data */

/** detector handler data */
struct GCG_DetectorData
{
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
#define freeGeneralmastersetcover NULL

/** destructor of detector to free detector data (called before the solving process begins) */
#define exitGeneralmastersetcover NULL

/** detection initialization function of detector (called before solving is about to begin) */
#define initGeneralmastersetcover NULL

static GCG_DECL_PROPAGATEPARTIALDEC(propagatePartialdecGeneralmastersetcover)
{
   *result = SCIP_DIDNOTFIND;
   SCIP* origprob = GCGgetOrigprob(gcg);
   char decinfo[SCIP_MAXSTRLEN];
   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(origprob, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(origprob, temporaryClock) );

   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nvars;
   bool relevant = true;


   gcg::PARTIALDECOMP* partialdec = partialdecdetectiondata->workonpartialdec;
   auto& openconss = partialdec->getOpenconss();
   for( auto itr = openconss.cbegin(); itr != openconss.cend(); )
   {
      bool found = false;
      cons = partialdecdetectiondata->detprobdata->getCons(*itr);

      /* set open setcovering and logicor constraints to master */
      if( GCGconsGetType(origprob, cons) == setcovering || GCGconsGetType(origprob, cons) == logicor )
      {
         itr = partialdec->fixConsToMaster(itr);
         found = true;
      }
      /* set constraints with infinity rhs and nonnegative lhs to master */
      else if(GCGconsGetType(origprob, cons) != logicor && GCGconsGetType(origprob, cons) != setpacking && GCGconsGetType(origprob, cons) != setpartitioning )
      {
         nvars = GCGconsGetNVars(origprob, cons);
         vars = NULL;
         vals = NULL;
         if( !SCIPisInfinity(origprob, GCGconsGetRhs(origprob, cons)) )
            relevant = false;
         if( SCIPisNegative(origprob, GCGconsGetLhs(origprob, cons)) )
            relevant = false;
         if( nvars > 0 )
         {
            SCIP_CALL( SCIPallocMemoryArray(origprob, &vars, nvars) );
            SCIP_CALL( SCIPallocMemoryArray(origprob, &vals, nvars) );
            SCIP_CALL( GCGconsGetVars(origprob, cons, vars, nvars) );
            SCIP_CALL( GCGconsGetVals(origprob, cons, vals, nvars) );
         }
         for( int j = 0; j < nvars && relevant; ++j )
         {
            assert(vars != NULL);
            assert(vals != NULL);

            if( !SCIPvarIsIntegral(vars[j]) && !SCIPvarIsBinary(vars[j]) )
            {
               SCIPdebugPrintf("(%s is not integral) ", SCIPvarGetName(vars[j]) );
               relevant = FALSE;
            }
            if( !SCIPisEQ(origprob, vals[j], 1.0) )
            {
               SCIPdebugPrintf("(coeff for var %s is %.2f != 1.0) ", SCIPvarGetName(vars[j]), vals[j] );
               relevant = FALSE;
            }
         }
         SCIPfreeMemoryArrayNull(origprob, &vals);
         SCIPfreeMemoryArrayNull(origprob, &vars);

         if(relevant)
         {
            itr = partialdec->fixConsToMaster(itr);
            found = true;
         }
      }
      if( !found )
      {
         ++itr;
      }
   }

   partialdec->sort();
   SCIP_CALL_ABORT( SCIPstopClock(origprob, temporaryClock) );

   partialdecdetectiondata->detectiontime = SCIPgetClockTime(origprob, temporaryClock);
   SCIP_CALL( SCIPallocMemoryArray(origprob, &(partialdecdetectiondata->newpartialdecs), 1) );
   partialdecdetectiondata->newpartialdecs[0] = partialdec;
   partialdecdetectiondata->nnewpartialdecs = 1;
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "genmastersetcover");
   partialdecdetectiondata->newpartialdecs[0]->addDetectorChainInfo(decinfo);
   partialdecdetectiondata->newpartialdecs[0]->addClockTime(SCIPgetClockTime(origprob, temporaryClock));
   // we used the provided partialdec -> prevent deletion
   partialdecdetectiondata->workonpartialdec = NULL;
   SCIP_CALL_ABORT(SCIPfreeClock(origprob, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

#define finishPartialdecGeneralmastersetcover NULL
#define detectorPostprocessPartialdecGeneralmastersetcover NULL

static
GCG_DECL_SETPARAMAGGRESSIVE(setParamAggressiveGeneralmastersetcover)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   int newval;
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE ) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxcallround", name);
   SCIP_CALL( SCIPgetIntParam(origprob, setstr, &newval) );
   ++newval;
   SCIP_CALL( SCIPsetIntParam(origprob, setstr, newval ) );
   SCIPinfoMessage(origprob, NULL, "After Setting %s = %d\n", setstr, newval);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origmaxcallround", name);
   SCIP_CALL( SCIPgetIntParam(origprob, setstr, &newval) );
   ++newval;
   SCIP_CALL( SCIPsetIntParam(origprob, setstr, newval ) );
   SCIPinfoMessage(origprob, NULL, "%s = %d\n", setstr, newval);

   return SCIP_OKAY;
}


static
GCG_DECL_SETPARAMDEFAULT(setParamDefaultGeneralmastersetcover)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, DEC_ENABLED) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, DEC_ENABLEDFINISHING ) );

   return SCIP_OKAY;
}

static
GCG_DECL_SETPARAMFAST(setParamFastGeneralmastersetcover)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE ) );

   return SCIP_OKAY;
}



/*
 * detector specific interface methods
 */

/** creates the handler for generalmastersetcover detector and includes it in SCIP */
SCIP_RETCODE GCGincludeDetectorGeneralmastersetcover(
   GCG*                 gcg                  /**< GCG data structure */
   )
{
   GCG_DETECTORDATA* detectordata;

   /**@todo create generalmastersetcover detector data here*/
   detectordata = NULL;

   SCIP_CALL(
      GCGincludeDetector(gcg, DEC_NAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND,
                         DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY,
                         DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL,
                         detectordata, freeGeneralmastersetcover, initGeneralmastersetcover, exitGeneralmastersetcover, propagatePartialdecGeneralmastersetcover, finishPartialdecGeneralmastersetcover, detectorPostprocessPartialdecGeneralmastersetcover, setParamAggressiveGeneralmastersetcover, setParamDefaultGeneralmastersetcover, setParamFastGeneralmastersetcover));

   /**@todo add generalmastersetcover detector parameters */

   return SCIP_OKAY;
}









