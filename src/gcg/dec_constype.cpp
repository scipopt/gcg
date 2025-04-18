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

/**@file   dec_constype.cpp
 * 
 * @brief  detector constype (put your description here)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/dec_constype.h"
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
#define DEC_NAME                  "constype"       /**< name of detector */
#define DEC_DESC                  "detector constype" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          0           /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  0     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               't'         /**< display character of detector */
#define DEC_ENABLED               FALSE        /**< should the detection be enabled */
#define DEC_ENABLEDFINISHING      FALSE       /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the finishing be enabled */
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

/** method to enumerate all subsets */
std::vector< std::vector<int> > getSubsets(std::vector<int> set)
{
    std::vector< std::vector<int> > subset;
    std::vector<int> empty;
    subset.push_back( empty );

    for ( size_t i = 0; i < set.size(); ++i )
    {
        std::vector< std::vector<int> > subsetTemp = subset;

        for (size_t j = 0; j < subsetTemp.size(); ++j)
            subsetTemp[j].push_back( set[i] );

        for (size_t j = 0; j < subsetTemp.size(); ++j)
            subset.push_back( subsetTemp[j] );
    }
    return subset;
}

/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
#define freeConstype NULL

/** destructor of detector to free detector data (called before the solving process begins) */
#define exitConstype NULL

/** detection initialization function of detector (called before solving is about to begin) */
#define initConstype NULL

static GCG_DECL_PROPAGATEPARTIALDEC(propagatePartialdecConstype)
{
   *result = SCIP_DIDNOTFIND;
   char decinfo[SCIP_MAXSTRLEN];
   SCIP* origprob = GCGgetOrigprob(gcg);

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT( SCIPcreateClock(origprob, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(origprob, temporaryClock) );

   SCIP_CONS* cons;

   int partialdecCounter = 0;
   gcg::PARTIALDECOMP* partialdecOrig;
   gcg::PARTIALDECOMP* partialdec;

   std::vector<consType> foundConstypes(0);
   std::vector<int> constypesIndices(0);

   partialdecOrig = partialdecdetectiondata->workonpartialdec;

   for( int i = 0; i < partialdecOrig->getNOpenconss(); ++i)
   {
      cons = partialdecdetectiondata->detprobdata->getCons(partialdecOrig->getOpenconss()[i]);
      consType cT = GCGconsGetType(origprob, cons);

      /* find constype or not */
      std::vector<consType>::const_iterator constypeIter = foundConstypes.begin();
      for(; constypeIter != foundConstypes.end(); ++constypeIter)
      {
    	  if(*constypeIter == cT)
    		  break;
      }

      if( constypeIter  == foundConstypes.end()  )
      {
         foundConstypes.push_back(GCGconsGetType(origprob, cons) );
      }
   }

   for(size_t i = 0; i < foundConstypes.size(); ++i)
   {
      constypesIndices.push_back((int) i);
   }

   std::vector< std::vector<int> > subsetsOfConstypes = getSubsets(constypesIndices);

   SCIP_CALL( SCIPallocMemoryArray(origprob, &(partialdecdetectiondata->newpartialdecs), subsetsOfConstypes.size() - 1) );
   partialdecdetectiondata->nnewpartialdecs = (int) subsetsOfConstypes.size() - 1;

   for(size_t subset = 0; subset < subsetsOfConstypes.size(); ++subset)
   {
      if(subsetsOfConstypes[subset].size() == 0)
          continue;

      partialdec = new gcg::PARTIALDECOMP(partialdecOrig);
      /* set open cons that have type of the current subset to Master */
      auto& openconss = partialdec->getOpenconss();
      for( auto itr = openconss.cbegin(); itr != openconss.cend(); )
      {
         cons = partialdecdetectiondata->detprobdata->getCons(*itr);
         bool found = false;
         for(size_t constypeId = 0; constypeId < subsetsOfConstypes[subset].size(); ++constypeId )
         {
            if( GCGconsGetType(origprob, cons) == foundConstypes[subsetsOfConstypes[subset][constypeId]] )
            {
               itr = partialdec->fixConsToMaster(itr);
               found = true;
               break;
            }
         }
         if( !found )
         {
            ++itr;
         }
      }
      partialdec->sort();
      (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "constype-%llu", subset);
      partialdec->addDetectorChainInfo(decinfo);
      partialdecdetectiondata->newpartialdecs[partialdecCounter] = partialdec;
      partialdecCounter++;
   }

   SCIP_CALL_ABORT( SCIPstopClock(origprob, temporaryClock) );
   partialdecdetectiondata->detectiontime = SCIPgetClockTime(origprob, temporaryClock);
   for( int s = 0; s < partialdecdetectiondata->nnewpartialdecs; ++s )
   {
      partialdecdetectiondata->newpartialdecs[s]->addClockTime(partialdecdetectiondata->detectiontime / partialdecdetectiondata->nnewpartialdecs);
   }
   SCIP_CALL_ABORT(SCIPfreeClock(origprob, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

#define finishPartialdecConstype NULL
#define detectorPostprocessPartialdecConstype NULL
static
GCG_DECL_SETPARAMAGGRESSIVE(setParamAggressiveConstype)
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
GCG_DECL_SETPARAMDEFAULT(setParamDefaultConstype)
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
GCG_DECL_SETPARAMFAST(setParamFastConstype)
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

/** creates the handler for constype detector and includes it in SCIP */
SCIP_RETCODE GCGincludeDetectorConstype(
   GCG*                 gcg                  /**< SCIP data structure */
)
{
   GCG_DETECTORDATA* detectordata;

   /**@todo create constype detector data here*/
   detectordata = NULL;

   SCIP_CALL(
      GCGincludeDetector(gcg, DEC_NAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND,
                         DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata,
                         freeConstype, initConstype, exitConstype, propagatePartialdecConstype, finishPartialdecConstype, detectorPostprocessPartialdecConstype, setParamAggressiveConstype, setParamDefaultConstype, setParamFastConstype));

   /**@todo add constype detector parameters */

   return SCIP_OKAY;
}
