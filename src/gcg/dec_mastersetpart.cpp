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

/**@file   dec_mastersetpart.cpp
 * 
 * @brief  detector mastersetpart (set setpartitioning constraints to master)
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/dec_mastersetpart.h"
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
#define DEC_NAME                  "mastersetpart"       /**< name of detector */
#define DEC_DESC                  "detector mastersetpart" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /**< frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /**< last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /**< first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /**< frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /**< last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /**< first round the detector gets called while detecting the original problem    */
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
#define freeMastersetpart NULL

/** destructor of detector to free detector data (called before the solving process begins) */
#define exitMastersetpart NULL

/** detection initialization function of detector (called before solving is about to begin) */
#define initMastersetpart NULL


static GCG_DECL_PROPAGATEPARTIALDEC(propagatePartialdecMastersetpart)
{
   *result = SCIP_DIDNOTFIND;
   SCIP* origprob = GCGgetOrigprob(gcg);
   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT( SCIPcreateClock(origprob, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(origprob, temporaryClock) );

   SCIP_CONS* cons;

   gcg::PARTIALDECOMP* partialdec = partialdecdetectiondata->workonpartialdec;

   /* set open setpartitioning constraints to Master */
   auto& openconss = partialdec->getOpenconss();
   for( auto itr = openconss.cbegin(); itr != openconss.cend(); )
   {
      cons = partialdecdetectiondata->detprobdata->getCons(*itr);
      if( GCGconsGetType(origprob, cons) == setpartitioning )
      {
          itr = partialdec->fixConsToMaster(itr);
      }
      else
      {
         ++itr;
      }
   }

   partialdec->sort();
   SCIP_CALL_ABORT( SCIPstopClock(origprob, temporaryClock ) );

   partialdecdetectiondata->detectiontime = SCIPgetClockTime(origprob, temporaryClock);
   SCIP_CALL( SCIPallocMemoryArray(origprob, &(partialdecdetectiondata->newpartialdecs), 1) );
   partialdecdetectiondata->newpartialdecs[0] = partialdec;
   partialdecdetectiondata->nnewpartialdecs = 1;
   partialdecdetectiondata->newpartialdecs[0]->addClockTime(SCIPgetClockTime(origprob, temporaryClock));
   partialdecdetectiondata->newpartialdecs[0]->addDetectorChainInfo(DEC_NAME);
   // we used the provided partialdec -> prevent deletion
   partialdecdetectiondata->workonpartialdec = NULL;
   SCIP_CALL_ABORT(SCIPfreeClock(origprob, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

#define finishPartialdecMastersetpart NULL

#define detectorPostprocessPartialdecMastersetpart NULL

#define setParamAggressiveMastersetpart NULL
#define setParamDefaultMastersetpart NULL
#define setParamFastMastersetpart NULL


/*
 * detector specific interface methods
 */

/** creates the handler for mastersetpart detector and includes it in SCIP */
SCIP_RETCODE GCGincludeDetectorMastersetpart(
   GCG*                    gcg                     /**< GCG data structure */
   )
{
   GCG_DETECTORDATA* detectordata;

   /**@todo create mastersetpart detector data here*/
   detectordata = NULL;

   SCIP_CALL(
      GCGincludeDetector(gcg, DEC_NAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata, freeMastersetpart, initMastersetpart, exitMastersetpart, propagatePartialdecMastersetpart, finishPartialdecMastersetpart, detectorPostprocessPartialdecMastersetpart, setParamAggressiveMastersetpart, setParamDefaultMastersetpart, setParamFastMastersetpart));

   /**@todo add mastersetpart detector parameters */

   return SCIP_OKAY;
}




