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

/**@file   dec_postprocess.c
 * @ingroup DETECTORS
 * @brief  checks if there are master constraints that can be assigned to one block (without any other changes)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_postprocess.h"
#include "cons_decomp.h"
#include "gcg.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "scip/scip.h"
#include "scip_misc.h"
#include "scip/clock.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>

/* constraint handler properties */
#define DEC_DETECTORNAME          "postprocess"       /**< name of detector */
#define DEC_DESC                  "detector postprocess" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_PRIORITY              1000000     /**< priority of the constraint handler for separation */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_DECCHAR               'p'         /**< display character of detector */
#define DEC_ENABLED               FALSE        /**< should the detection be enabled */
#define DEC_ENABLEDORIGINAL       FALSE  /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING      FALSE        /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING TRUE          /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */
#define DEC_LEGACYMODE            FALSE       /**< should (old) DETECTSTRUCTURE method also be used for detection */
#define DEFAULT_USECONSSADJ       TRUE
/*
 * Data structures
 */

/** @todo fill in the necessary detector data */

/** detector handler data */
struct DEC_DetectorData
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
DEC_DECL_FREEDETECTOR(freePostprocess)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}




/** destructor of detector to free detector data (called before the solving process begins) */
#if 0
static
DEC_DECL_EXITDETECTOR(exitPostprocess)
{  /*lint --e{715}*/

   SCIPerrorMessage("Exit function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
#else
#define exitPostprocess NULL
#endif

/** detection initialization function of detector (called before solving is about to begin) */
#if 0
static
DEC_DECL_INITDETECTOR(initPostprocess)
{  /*lint --e{715}*/

   SCIPerrorMessage("Init function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define initPostprocess NULL
#endif

/** detection function of detector */
//static
//DEC_DECL_DETECTSTRUCTURE(detectPostprocess)
//{ /*lint --e{715}*/
//   *result = SCIP_DIDNOTFIND;
//
//   SCIPerrorMessage("Detection function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
//   SCIPABORT();  /*lint --e{527}*/
//
//   return SCIP_OKAY;
//}

#define detectPostprocess NULL
#define propagateSeeedPostprocess NULL
#define finishSeeedPostprocess NULL

static
DEC_DECL_POSTPROCESSSEEED(postprocessSeeedPostprocess)
{
   *result = SCIP_DIDNOTFIND;

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
   char decinfo[SCIP_MAXSTRLEN];
   SCIP_Bool success;
   SCIP_Bool byconssadj;
   SCIP_Bool conssadjcalculated;

   gcg::Seeed* seeed;

   assert(seeedPropagationData->seeedToPropagate->getSeeedpool() == seeedPropagationData->seeedpool);
   seeed  = new gcg::Seeed(seeedPropagationData->seeedToPropagate);
   assert(scip == seeedPropagationData->seeedpool->getScip() );

   SCIPgetBoolParam(scip, "detection/detectors/postprocess/useconssadj", &byconssadj);
   SCIPgetBoolParam(scip, "detection/conssadjcalculated", &conssadjcalculated);


   //complete the seeed by bfs
   if ( byconssadj && conssadjcalculated)
      seeed->postprocessMasterToBlocksConssAdjacency( &success );
   else
      seeed->postprocessMasterToBlocks( &success );
  

   if ( !success )
   {
     seeedPropagationData->nNewSeeeds = 0;
     delete seeed;
     *result = SCIP_DIDNOTFIND;
     SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
     SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );
     return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), 1) );
   seeedPropagationData->newSeeeds[0] = seeed;
   seeedPropagationData->nNewSeeeds = 1;
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "postprocess");
   seeedPropagationData->newSeeeds[0]->addDetectorChainInfo(decinfo);

   seeedPropagationData->newSeeeds[0]->buildDecChainString();

   SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
   seeedPropagationData->newSeeeds[0]->addClockTime( SCIPclockGetTime(temporaryClock )  );
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


static
DEC_DECL_SETPARAMAGGRESSIVE(setParamAggressivePostprocess)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE ) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/postprocessingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );


   return SCIP_OKAY;

}


static
DEC_DECL_SETPARAMDEFAULT(setParamDefaultPostprocess)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLED) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLEDORIGINAL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLEDFINISHING) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/postprocessingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLEDPOSTPROCESSING ) );


   return SCIP_OKAY;

}

static
DEC_DECL_SETPARAMFAST(setParamFastPostprocess)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/postprocessingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );


   return SCIP_OKAY;

}



/*
 * detector specific interface methods
 */

/** creates the handler for postprocess detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorPostprocess(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /**@todo create postprocess detector data here*/
   detectordata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   detectordata->useconssadj = TRUE;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDORIGINAL, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, DEC_LEGACYMODE, detectordata, detectPostprocess, freePostprocess,
      initPostprocess, exitPostprocess, propagateSeeedPostprocess, NULL, NULL, finishSeeedPostprocess,
      postprocessSeeedPostprocess, setParamAggressivePostprocess, setParamDefaultPostprocess, setParamFastPostprocess) );

   /* add consname detector parameters */
      /**@todo add postprocess detector parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/postprocess/useconssadj", "should the constraint adjacency be used", &detectordata->useconssadj, FALSE, DEFAULT_USECONSSADJ, NULL, NULL) );


   return SCIP_OKAY;
}
