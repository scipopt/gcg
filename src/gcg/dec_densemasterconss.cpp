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

/**@file   dec_densemasterconss.cpp
 * 
 * @brief  detector densemasterconss (put your description here)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/dec_densemasterconss.h"
#include "gcg/cons_decomp.h"
#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/gcg.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "gcg/scip_misc.h"
#include "scip/clock.h"

#include <sstream>

#include <iostream>
#include <algorithm>


/* constraint handler properties */
#define DEC_NAME                  "densemasterconss"       /**< name of detector */
#define DEC_DESC                  "detector densemasterconss" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          0           /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               'd'         /**< display character of detector */
#define DEC_ENABLED               TRUE        /**< should the detection be enabled */
#define DEC_ENABLEDFINISHING      FALSE        /**< should the detection be enabled */
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

struct sort_pred {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.first > right.first;
    }
};

/* put your local methods here, and declare them static */

/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
#define freeDensemasterconss NULL

/** destructor of detector to free detector data (called before the solving process begins) */
#define exitDensemasterconss NULL

/** detection initialization function of detector (called before solving is about to begin) */
#define initDensemasterconss NULL

#define detectDensemasterconss NULL

#define finishPartialdecDensemasterconss NULL

static GCG_DECL_PROPAGATEPARTIALDEC(propagatePartialdecDensemasterconss)
{
   *result = SCIP_DIDNOTFIND;
   char decinfo[SCIP_MAXSTRLEN];

   SCIP_CLOCK* temporaryClock;
   SCIP* origprob = GCGgetOrigprob(gcg);
   gcg::DETPROBDATA* detprobdata;
   gcg::PARTIALDECOMP* partialdec = partialdecdetectiondata->workonpartialdec;
   std::stringstream decdesc;

   SCIP_Real maxratio = 0.2;
   int maxdiff = -1;
   int maxdiffindex = -1;
   int lastindex = -1;

   detprobdata = partialdecdetectiondata->detprobdata;

   SCIP_CALL_ABORT( SCIPcreateClock(origprob, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(origprob, temporaryClock) );

   lastindex =  maxratio * detprobdata->getNConss();
   /* fix open conss that have a) type of the current subset or b) decomp info ONLY_MASTER as master conss */
   std::vector<std::pair<int,int>> nnonzeros(detprobdata->getNConss(), std::pair<int, int>(0, -1));
   for( int cons: partialdec->getOpenconss() )
   {
      nnonzeros[cons] = std::pair<int,int>(detprobdata->getNVarsForCons(cons), cons);
   }

   std::sort(nnonzeros.begin(), nnonzeros.end(), sort_pred() );

   for( int i = 0; i < lastindex && i < (int) nnonzeros.size() - 1; ++i )
   {
      if( maxdiff < nnonzeros[i].first - nnonzeros[i+1].first )
      {
         maxdiff = nnonzeros[i].first - nnonzeros[i+1].first;
         maxdiffindex = i;
      }
   }

   for( int i = 0; i < maxdiffindex; ++i )
   {
      partialdec->fixConsToMaster(nnonzeros[i].second);
   }

   decdesc << "densemasterconss" << "\\_" << maxdiffindex ;

   partialdec->sort();
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, decdesc.str().c_str());
   partialdec->addDetectorChainInfo(decinfo);

   SCIP_CALL_ABORT( SCIPstopClock(origprob, temporaryClock ) );

   partialdecdetectiondata->detectiontime = SCIPgetClockTime(origprob, temporaryClock);
   SCIP_CALL( SCIPallocMemoryArray(origprob, &(partialdecdetectiondata->newpartialdecs), 1) );
   partialdecdetectiondata->nnewpartialdecs  = 1;

   SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "dec_densemasterconss found %d new partialdec \n", partialdecdetectiondata->nnewpartialdecs  );

   partialdecdetectiondata->newpartialdecs[0] = partialdec;
   partialdecdetectiondata->newpartialdecs[0]->addClockTime(SCIPgetClockTime(origprob, temporaryClock));
   // we used the provided partialdec -> prevent deletion
   partialdecdetectiondata->workonpartialdec = NULL;

   SCIP_CALL_ABORT(SCIPfreeClock(origprob, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

#define detectorPostprocessPartialdecDensemasterconss NULL

static
GCG_DECL_SETPARAMAGGRESSIVE(setParamAggressiveDensemasterconss)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE ) );

   return SCIP_OKAY;
}


static
GCG_DECL_SETPARAMDEFAULT(setParamDefaultDensemasterconss)
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
GCG_DECL_SETPARAMFAST(setParamFastDensemasterconss)
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

/** creates the handler for densemasterconss detector and includes it in SCIP */
SCIP_RETCODE GCGincludeDetectorDensemasterconss(
   GCG*                 gcg                  /**< GCG data structure */
   )
{
   GCG_DETECTORDATA* detectordata;

   /**@todo create densemasterconss detector data here*/
   detectordata = NULL;

   SCIP_CALL(
      GCGincludeDetector(gcg, DEC_NAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND,
                         DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata,
                         freeDensemasterconss, initDensemasterconss, exitDensemasterconss, propagatePartialdecDensemasterconss, finishPartialdecDensemasterconss, detectorPostprocessPartialdecDensemasterconss, setParamAggressiveDensemasterconss, setParamDefaultDensemasterconss, setParamFastDensemasterconss));

   /**@todo add densemasterconss detector parameters */

   return SCIP_OKAY;
}
