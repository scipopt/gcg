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

/**@file   dec_postprocess.cpp
 * 
 * @brief  checks if there are master constraints that can be assigned to one block (without any other changes)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/dec_postprocess.h"
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
#define DEC_NAME                  "postprocess"       /**< name of detector */
#define DEC_DESC                  "detector postprocess" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /**< frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /**< last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /**< first round the detector gets called                              */
#define DEC_PRIORITY              1000000     /**< priority of the constraint handler for separation */
#define DEC_FREQCALLROUNDORIGINAL 1           /**< frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /**< last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /**< first round the detector gets called while detecting the original problem    */
#define DEC_DECCHAR               'p'         /**< display character of detector */
#define DEC_ENABLED               FALSE        /**< should the detection be enabled */
#define DEC_ENABLEDFINISHING      FALSE        /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING TRUE          /**< should the postprocessing be enabled */
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
GCG_DECL_FREEDETECTOR(freePostprocess)
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
#define exitPostprocess NULL

/** detection initialization function of detector (called before solving is about to begin) */
#define initPostprocess NULL

#define propagatePartialdecPostprocess NULL
#define finishPartialdecPostprocess NULL

static
GCG_DECL_POSTPROCESSPARTIALDEC(postprocessPartialdecPostprocess)
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   *result = SCIP_DIDNOTFIND;

   assert(partialdecdetectiondata->workonpartialdec->isComplete());

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT( SCIPcreateClock(origprob, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(origprob, temporaryClock) );
   char decinfo[SCIP_MAXSTRLEN];
   SCIP_Bool success;
   SCIP_Bool byconssadj;

   gcg::PARTIALDECOMP* partialdec = partialdecdetectiondata->workonpartialdec;
   gcg::DETPROBDATA* detprobdata = partialdecdetectiondata->detprobdata;
   assert(partialdecdetectiondata->workonpartialdec->getDetprobdata() == detprobdata);

   SCIPgetBoolParam(origprob, "detection/detectors/postprocess/useconssadj", &byconssadj);

   if ( byconssadj && !detprobdata->isConssAdjInitialized() )
      detprobdata->createConssAdjacency();

   //complete the partialdec by bfs
   if ( byconssadj )
   {
      success = FALSE;
      std::vector<int> constoreassign;
      std::vector<int> blockforconstoreassign;

      partialdec->sort();

      std::vector<int> blockforvar(partialdec->getNVars(), -1 );

      for( int b = 0; b < partialdec->getNBlocks(); ++b )
      {
         for( size_t j  = 0; j < (size_t) partialdec->getNVarsForBlock(b); ++j )
         {
            blockforvar[partialdec->getVarsForBlock(b)[j]] = b;
         }
      }

      for( int mc = 0; mc < partialdec->getNMasterconss(); ++mc )
      {
         int masterconsid = partialdec->getMasterconss()[mc];
         int hittenblock  = -1;

         SCIP_Bool lockedcons = FALSE;

         for( int var = 0; var < detprobdata->getNVarsForCons(masterconsid); ++var )
         {
            int varid = detprobdata->getVarsForCons(masterconsid)[var];
            /* do not reassign cons that contain static master vars or stairlinking vars */
            /* @todo: we could check if cons can be moved without destroying the staircase structure*/
            if( partialdec->isVarMastervar(varid) || partialdec->isVarStairlinkingvar(varid) )
            {
               lockedcons = TRUE;
               break;
            }

            if ( blockforvar[varid] != -1 )
            {
               if( hittenblock == -1 )
                  hittenblock = blockforvar[varid];
               else if( hittenblock != blockforvar[varid] )
               {
                  lockedcons = TRUE;
                  break;
               }
            }
         }

         if( lockedcons )
            continue;

         if ( hittenblock != -1 )
         {
            constoreassign.push_back(masterconsid);
            blockforconstoreassign.push_back(hittenblock);
         }
      }

      for( size_t i = 0; i < constoreassign.size() ; ++i )
      {
         assert(partialdec->isConsMastercons(constoreassign[i]));
         partialdec->removeMastercons(constoreassign[i]);
         partialdec->setConsToBlock(constoreassign[i], blockforconstoreassign[i]);
      }

      if( !constoreassign.empty() )
         success = TRUE;

      partialdec->prepare();
   }
   else
      success = FALSE;

   SCIP_CALL_ABORT( SCIPstopClock(origprob, temporaryClock ) );

   assert(partialdec->checkConsistency());

   if ( !success )
   {
      partialdecdetectiondata->nnewpartialdecs = 0;
      *result = SCIP_DIDNOTFIND;
   }
   else
   {
      partialdecdetectiondata->detectiontime = SCIPgetClockTime(origprob, temporaryClock);
      SCIP_CALL( SCIPallocMemoryArray(origprob, &(partialdecdetectiondata->newpartialdecs), 1) );
      partialdecdetectiondata->newpartialdecs[0] = partialdec;
      partialdecdetectiondata->nnewpartialdecs = 1;
      (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "postprocess");
      partialdecdetectiondata->newpartialdecs[0]->addDetectorChainInfo(decinfo);
      partialdecdetectiondata->newpartialdecs[0]->addClockTime(SCIPgetClockTime(origprob, temporaryClock));
      // we used the provided partialdec -> prevent deletion
      partialdecdetectiondata->workonpartialdec = NULL;
      *result = SCIP_SUCCESS;
   }

   SCIP_CALL_ABORT( SCIPfreeClock(origprob, &temporaryClock) );

   return SCIP_OKAY;
}


static
GCG_DECL_SETPARAMAGGRESSIVE(setParamAggressivePostprocess)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, TRUE ) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/postprocessingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, TRUE) );


   return SCIP_OKAY;

}


static
GCG_DECL_SETPARAMDEFAULT(setParamDefaultPostprocess)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, DEC_ENABLED) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, DEC_ENABLEDFINISHING) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/postprocessingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, DEC_ENABLEDPOSTPROCESSING ) );


   return SCIP_OKAY;

}

static
GCG_DECL_SETPARAMFAST(setParamFastPostprocess)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE ) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/postprocessingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE ) );


   return SCIP_OKAY;

}



/*
 * detector specific interface methods
 */

/** creates the handler for postprocess detector and includes it in SCIP */
SCIP_RETCODE GCGincludeDetectorPostprocess(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_DETECTORDATA* detectordata;
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   /**@todo create postprocess detector data here*/
   detectordata = NULL;
   SCIP_CALL( SCIPallocMemory(origprob, &detectordata) );
   assert(detectordata != NULL);

   detectordata->useconssadj = TRUE;

   SCIP_CALL( GCGincludeDetector(gcg, DEC_NAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata, freePostprocess,
                                 initPostprocess, exitPostprocess, propagatePartialdecPostprocess, finishPartialdecPostprocess,
                                 postprocessPartialdecPostprocess, setParamAggressivePostprocess, setParamDefaultPostprocess, setParamFastPostprocess) );

   /* add consname detector parameters */
      /**@todo add postprocess detector parameters */
   SCIP_CALL( SCIPaddBoolParam(origprob, "detection/detectors/postprocess/useconssadj", "should the constraint adjacency be used", &detectordata->useconssadj, FALSE, DEFAULT_USECONSSADJ, NULL, NULL) );


   return SCIP_OKAY;
}
