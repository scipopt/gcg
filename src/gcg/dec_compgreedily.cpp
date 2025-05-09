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

/**@file   dec_compgreedily.cpp
 * 
 * @brief  detector compgreedily (assigns the open cons and open vars of the partialdec greedily)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/dec_compgreedily.h"
#include "gcg/cons_decomp.h"
#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "scip/clock.h"
#include <iostream>

/* constraint handler properties */
#define DEC_NAME                  "compgreedily"       /**< name of detector */
#define DEC_DESC                  "detector compgreedily" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               'g'         /**< display character of detector */
#define DEC_ENABLED               FALSE       /**< should the detection be enabled */
#define DEC_ENABLEDFINISHING      FALSE       /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE       /**< should the finishing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated partialdec */

/** parameter limits for emphasis default */

#define DEFAULT_LIMITHALFPERIMETERENABLEDFINISHING    20000   /** limit in terms of nrows + ncols for enabling finishing */
#define DEFAULT_LIMITHALFPERIMETERENABLEDORIGINAL     10000   /** limit in terms of nrows + ncols for enabling in detecting for unpresolved problem */


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

#define freeCompgreedily NULL

/** destructor of detector to free detector data (called before the solving process begins) */

#define exitCompgreedily NULL

#define initCompgreedily NULL

static
GCG_DECL_PROPAGATEPARTIALDEC(propagatePartialdecCompgreedily)
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   *result = SCIP_DIDNOTFIND;

   char decinfo[SCIP_MAXSTRLEN];
   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT( SCIPcreateClock(origprob, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(origprob, temporaryClock) );

   gcg::PARTIALDECOMP* partialdec = partialdecdetectiondata->workonpartialdec;

   //assign open conss and vars greedily
   partialdec->completeGreedily();

   SCIP_CALL_ABORT( SCIPstopClock(origprob, temporaryClock) );

   partialdecdetectiondata->detectiontime =  SCIPgetClockTime(origprob, temporaryClock);
   SCIP_CALL( SCIPallocMemoryArray(origprob, &(partialdecdetectiondata->newpartialdecs), 1) );
   partialdecdetectiondata->newpartialdecs[0] = partialdec;
   partialdecdetectiondata->nnewpartialdecs = 1;

   partialdecdetectiondata->newpartialdecs[0]->addClockTime(SCIPgetClockTime(origprob, temporaryClock));
   SCIP_CALL_ABORT(SCIPfreeClock(origprob, &temporaryClock) );
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "compgreed");
   partialdecdetectiondata->newpartialdecs[0]->addDetectorChainInfo(decinfo);
   // we used the provided partialdec -> prevent deletion
   partialdecdetectiondata->workonpartialdec = NULL;

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

static
GCG_DECL_FINISHPARTIALDEC(finishPartialdecCompgreedily)
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   *result = SCIP_DIDNOTFIND;
   char decinfo[SCIP_MAXSTRLEN];

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT (SCIPcreateClock(origprob, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(origprob, temporaryClock) );

   gcg::PARTIALDECOMP* partialdec = partialdecdetectiondata->workonpartialdec;

   //assign open conss and vars greedily
   partialdec->completeGreedily();

   SCIP_CALL_ABORT( SCIPstopClock(origprob, temporaryClock) );

   partialdecdetectiondata->detectiontime =  SCIPgetClockTime(origprob, temporaryClock);
   SCIP_CALL( SCIPallocMemoryArray(origprob, &(partialdecdetectiondata->newpartialdecs), 1) );
   partialdecdetectiondata->newpartialdecs[0] = partialdec;
   partialdecdetectiondata->nnewpartialdecs = 1;
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "compgreed");
   partialdecdetectiondata->newpartialdecs[0]->addDetectorChainInfo(decinfo);
   partialdecdetectiondata->newpartialdecs[0]->addClockTime(SCIPgetClockTime(origprob, temporaryClock));
   // we used the provided partialdec -> prevent deletion
   partialdecdetectiondata->workonpartialdec = NULL;

   SCIP_CALL_ABORT(SCIPfreeClock(origprob, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

#define detectorPostprocessPartialdecCompgreedily NULL


static
GCG_DECL_SETPARAMAGGRESSIVE(setParamAggressiveCompgreedily)
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


static
GCG_DECL_SETPARAMDEFAULT(setParamDefaultCompgreedily)
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
GCG_DECL_SETPARAMFAST(setParamFastCompgreedily)
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

/** creates the handler for compgreedily detector and includes it in SCIP */
SCIP_RETCODE GCGincludeDetectorCompgreedily(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_DETECTORDATA* detectordata;

   /**@todo create compgreedily detector data here*/
   detectordata = NULL;

   SCIP_CALL( GCGincludeDetector(gcg, DEC_NAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata, freeCompgreedily, initCompgreedily, exitCompgreedily, propagatePartialdecCompgreedily, finishPartialdecCompgreedily, detectorPostprocessPartialdecCompgreedily, setParamAggressiveCompgreedily, setParamDefaultCompgreedily, setParamFastCompgreedily) );

   /**@todo add compgreedily detector parameters */

   return SCIP_OKAY;
}
