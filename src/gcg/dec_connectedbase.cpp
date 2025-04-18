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

/**@file   dec_connectedbase.cpp
 * 
 * @brief  detector connectedbase (completes the partialdec by bfs)
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/dec_connectedbase.h"
#include "gcg/cons_decomp.h"
#include "gcg/gcg.h"
#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "scip/scip.h"
#include "gcg/scip_misc.h"
#include "scip/clock.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>

/* constraint handler properties */
#define DEC_NAME                  "connectedbase"       /**< name of detector */
#define DEC_DESC                  "detector connectedbase" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_DECCHAR               'C'         /**< display character of detector */
#define DEC_ENABLED               FALSE       /**< should the detection be enabled */
#define DEC_ENABLEDFINISHING      TRUE        /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE       /**< should the finishing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated partialdec */
#define DEFAULT_USECONSSADJ       TRUE
/*
 * Data structures
 */

/** @todo fill in the necessary detector data */

/** detector handler data */
struct GCG_DetectorData
{
   SCIP_Bool useconssadj;
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
/** destructor of detector to free detector data (called when SCIP is exiting) */
static
GCG_DECL_FREEDETECTOR(freeConnectedbase)
{  /*lint --e{715}*/
   GCG_DETECTORDATA *detectordata;
   assert(detector != NULL);

   assert(strcmp(GCGdetectorGetName(detector), DEC_NAME) == 0);

   detectordata = GCGdetectorGetData(detector);
   assert(detectordata != NULL);

   SCIPfreeMemory(GCGgetOrigprob(gcg), &detectordata);

   return SCIP_OKAY;
}


/** destructor of detector to free detector data (called before the solving process begins) */
#define exitConnectedbase NULL

/** detection initialization function of detector (called before solving is about to begin) */
#define initConnectedbase NULL

#define propagatePartialdecConnectedbase NULL


static
GCG_DECL_FINISHPARTIALDEC(finishPartialdecConnectedbase)
{
   *result = SCIP_DIDNOTFIND;

   SCIP_CLOCK* temporaryClock;
   char decinfo[SCIP_MAXSTRLEN];

   SCIP_Bool byconssadj;

   gcg::PARTIALDECOMP* partialdec = partialdecdetectiondata->workonpartialdec;
   gcg::DETPROBDATA* detprobdata = partialdecdetectiondata->detprobdata;
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);
   assert(partialdecdetectiondata->workonpartialdec->getDetprobdata() == detprobdata);

   SCIPgetBoolParam(origprob, "detection/detectors/connectedbase/useconssadj", &byconssadj);

   if ( byconssadj && !detprobdata->isConssAdjInitialized() )
      detprobdata->createConssAdjacency();

   SCIP_CALL_ABORT(SCIPcreateClock(origprob, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(origprob, temporaryClock) );

   //complete the partialdec by bfs
   if( byconssadj && partialdec->getNLinkingvars() == 0 )
      partialdec->completeByConnectedConssAdjacency();
   else
      partialdec->completeByConnected();

   SCIP_CALL_ABORT( SCIPstopClock(origprob, temporaryClock) );

   partialdecdetectiondata->detectiontime = SCIPgetClockTime(origprob, temporaryClock);
   SCIP_CALL( SCIPallocMemoryArray(origprob, &(partialdecdetectiondata->newpartialdecs), 1) );
   partialdecdetectiondata->newpartialdecs[0] = partialdec;
   partialdecdetectiondata->nnewpartialdecs = 1;
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "connected");
   partialdecdetectiondata->newpartialdecs[0]->addDetectorChainInfo(decinfo);

   partialdecdetectiondata->newpartialdecs[0]->addClockTime( SCIPgetClockTime(origprob, temporaryClock)  );
   SCIP_CALL_ABORT(SCIPfreeClock(origprob, &temporaryClock) );
   // we used the provided partialdec -> prevent deletion
   partialdecdetectiondata->workonpartialdec = NULL;

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


#define detectorPostprocessPartialdecConnectedbase NULL


static
GCG_DECL_SETPARAMAGGRESSIVE(setParamAggressiveConnectedbase)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, TRUE ) );


   return SCIP_OKAY;

}


static
GCG_DECL_SETPARAMDEFAULT(setParamDefaultConnectedbase)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, DEC_ENABLED) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, DEC_ENABLEDFINISHING ) );

   return SCIP_OKAY;

}


static
GCG_DECL_SETPARAMFAST(setParamFastConnectedbase)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, TRUE ) );

   return SCIP_OKAY;

}


/*
 * detector specific interface methods
 */

/** creates the handler for connectedbase detector and includes it in SCIP */
SCIP_RETCODE GCGincludeDetectorConnectedbase(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_DETECTORDATA* detectordata;
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   /**@todo create connectedbase detector data here*/
   detectordata = NULL;
   SCIP_CALL( SCIPallocMemory(origprob, &detectordata) );
   assert(detectordata != NULL);

   detectordata->useconssadj = TRUE;

   SCIP_CALL( GCGincludeDetector(gcg, DEC_NAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata, freeConnectedbase, initConnectedbase, exitConnectedbase, propagatePartialdecConnectedbase, finishPartialdecConnectedbase, detectorPostprocessPartialdecConnectedbase, setParamAggressiveConnectedbase, setParamDefaultConnectedbase, setParamFastConnectedbase) );

   /* add consname detector parameters */
      /**@todo add connectedbase detector parameters */
   SCIP_CALL( SCIPaddBoolParam(origprob, "detection/detectors/connectedbase/useconssadj", "should the constraint adjacency be used", &detectordata->useconssadj, FALSE, DEFAULT_USECONSSADJ, NULL, NULL) );

   return SCIP_OKAY;
}
