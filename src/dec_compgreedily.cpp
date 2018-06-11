
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file   dec_compgreedily.cpp
 * @ingroup DETECTORS
 * @brief  detector compgreedily (assigns the open cons and open vars of the seeed greedily)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_compgreedily.h"
#include "cons_decomp.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "scip/clock.h"
#include <iostream>

/* constraint handler properties */
#define DEC_DETECTORNAME          "compgreedily"       /**< name of detector */
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
#define DEC_ENABLEDORIGINAL       FALSE        /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING      FALSE        /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the finishing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */
#define DEC_LEGACYMODE            FALSE       /**< should (old) DETECTSTRUCTURE method also be used for detection */

/** parameter limits for emphasis default */

#define DEFAULT_LIMITHALFPERIMETERENABLEDFINISHING    20000   /** limit in terms of nrows + ncols for enabling finishing */
#define DEFAULT_LIMITHALFPERIMETERENABLEDORIGINAL     10000   /** limit in terms of nrows + ncols for enabling in detecting for unpresolved problem */


/*
 * Data structures
 */

/** @todo fill in the necessary detector data */

/** detector handler data */
struct DEC_DetectorData
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
#if 0
static
DEC_DECL_EXITDETECTOR(exitCompgreedily)
{  /*lint --e{715}*/

   SCIPerrorMessage("Exit function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define exitCompgreedily NULL
#endif

/** detection initialization function of detector (called before solving is about to begin) */
#if 0
static
DEC_DECL_INITDETECTOR(initCompgreedily)
{  /*lint --e{715}*/

   SCIPerrorMessage("Init function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define initCompgreedily NULL
#endif

/** detection function of detector */
//static
//DEC_DECL_DETECTSTRUCTURE(detectCompgreedily)
//{ /*lint --e{715}*/
//   *result = SCIP_DIDNOTFIND;
//
//   SCIPerrorMessage("Detection function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
//   SCIPABORT();  /*lint --e{527}*/
//
//   return SCIP_OKAY;
//}

#define detectCompgreedily NULL

static
DEC_DECL_PROPAGATESEEED(propagateSeeedCompgreedily)
{
   *result = SCIP_DIDNOTFIND;

   char decinfo[SCIP_MAXSTRLEN];
   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );

   gcg::Seeed* seeed = seeedPropagationData->seeedToPropagate;
   //seeed = new gcg::Seeed(seeedPropagationData->seeedToPropagate);

   //assign open conss and vars greedily
   seeed->completeGreedily();


   seeed->setDetectorPropagated(detector);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), 1) );
   seeedPropagationData->newSeeeds[0] = seeed;
   seeedPropagationData->nNewSeeeds = 1;

   SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
   seeedPropagationData->newSeeeds[0]->addClockTime( SCIPclockGetTime(temporaryClock )  );
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "compgreed");
   seeedPropagationData->newSeeeds[0]->addDetectorChainInfo(decinfo);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

static
DEC_DECL_FINISHSEEED(finishSeeedCompgreedily)
{
   *result = SCIP_DIDNOTFIND;
   char decinfo[SCIP_MAXSTRLEN];

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );

   gcg::Seeed* seeed = NULL;
   seeed = new gcg::Seeed(seeedPropagationData->seeedToPropagate);

   //assign open conss and vars greedily
   seeed->completeGreedily();

   SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), 1) );
   seeedPropagationData->newSeeeds[0] = seeed;
   seeedPropagationData->nNewSeeeds = 1;
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "compgreed");
   seeedPropagationData->newSeeeds[0]->addDetectorChainInfo(decinfo);

   SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
   seeedPropagationData->newSeeeds[0]->addClockTime( SCIPclockGetTime(temporaryClock )  );
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

#define detectorPostprocessSeeedCompgreedily NULL


static
DEC_DECL_SETPARAMAGGRESSIVE(setParamAggressiveCompgreedily)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );


   return SCIP_OKAY;

}


static
DEC_DECL_SETPARAMDEFAULT(setParamDefaultCompgreedily)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM && SCIPgetNOrigConss(scip) + SCIPgetNOrigVars(scip) < DEFAULT_LIMITHALFPERIMETERENABLEDORIGINAL )
      SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );
   else
      SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );


   return SCIP_OKAY;

}

static
DEC_DECL_SETPARAMFAST(setParamFastCompgreedily)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );

   return SCIP_OKAY;

}





/*
 * detector specific interface methods
 */

/** creates the handler for compgreedily detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorCompgreedily(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /**@todo create compgreedily detector data here*/
   detectordata = NULL;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDORIGINAL, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, DEC_LEGACYMODE, detectordata, detectCompgreedily, freeCompgreedily,initCompgreedily, exitCompgreedily, propagateSeeedCompgreedily, NULL, NULL, finishSeeedCompgreedily, detectorPostprocessSeeedCompgreedily, setParamAggressiveCompgreedily, setParamDefaultCompgreedily, setParamFastCompgreedily) );

   /**@todo add compgreedily detector parameters */

   return SCIP_OKAY;
}
