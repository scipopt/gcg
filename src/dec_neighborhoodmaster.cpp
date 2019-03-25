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

/**@file   dec_neighborhoodmaster.cpp
 * @ingroup DETECTORS
 * @brief  detector neighborhoodmaster (This detector calculates cons-cons adjacency (if not already done), sorts constraints according size of neighborhood. Searching two consecutive constraints with largest size difference (according neighborhood size) in sorted constraints. All constraints having a larger neighborhood than the second one are assigned to the master)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_neighborhoodmaster.h"
#include "cons_decomp.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "gcg.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "scip_misc.h"
#include "scip/clock.h"

#include <sstream>

#include <iostream>
#include <algorithm>

/**
This detector calculates cons-cons adjacency (if not already done), and sorts constraints according size of neighborhood. Searching two consecutive constraints with largest size difference (according neighborhood size) in sorted constraints. All constraints having a larger neighborhood than the second one are assigned to the master
*/

/* constraint handler properties */
#define DEC_DETECTORNAME          "neighborhoodmaster"       /**< name of detector */
#define DEC_DESC                  "detector neighborhoodmaster" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          0           /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               'n'         /**< display character of detector */
#define DEC_ENABLED               TRUE        /**< should the detection be enabled */
#define DEC_ENABLEDORIGINAL       FALSE       /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING      FALSE       /**< should the detection be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE       /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */
#define DEC_LEGACYMODE            FALSE       /**< should (old) DETECTSTRUCTURE method also be used for detection */

#define DEFAULT_MAXRATIO          0.2

/*
 * Data structures
 */

/** detector handler data */
struct DEC_DetectorData
{
   SCIP_Real maxratio;
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
DEC_DECL_FREEDETECTOR(freeNeighborhoodmaster)
{ /*lint --e{715}*/

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
#define exitNeighborhoodmaster NULL

/** detection initialization function of detector (called before solving is about to begin) */
#define initNeighborhoodmaster NULL

#define detectNeighborhoodmaster NULL

#define finishSeeedNeighborhoodmaster NULL

static DEC_DECL_PROPAGATESEEED(propagateSeeedNeighborhoodmaster)
{
  *result = SCIP_DIDNOTFIND;
  char decinfo[SCIP_MAXSTRLEN];

  SCIP_CLOCK* temporaryClock;

  gcg::Seeedpool* seeedpool;
  std::vector<gcg::Seeed*> foundseeeds(0);

  gcg::Seeed* seeedOrig;
  gcg::Seeed* seeed;

  DEC_DetectorData* detectorData;

  detectorData = DECdetectorGetData(detector);

  seeedOrig = seeedPropagationData->seeedToPropagate;
  std::stringstream decdesc;

  int maxdiff = -1;
  int maxdiffindex = -1;
  int lastindex = -1;

  seeedpool = seeedPropagationData->seeedpool;

  if ( !seeedpool->isConssAdjInitilized() )
     seeedpool->createConssAdjacency();

  seeed = new gcg::Seeed(seeedOrig);
  std::vector<std::pair<int,int>> neighborhoodsize = std::vector<std::pair<int,int>>(seeed->getNOpenconss(), std::pair<int, int>(0,-1)  );

  SCIP_CALL_ABORT( SCIPcreateClock(scip, &temporaryClock) );
  SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );

  lastindex = (int) (detectorData->maxratio * seeed->getNOpenconss());

  /** book open conss that have a) type of the current subset or b) decomp info ONLY_MASTER as master conss */
  for( int i = 0; i < seeed->getNOpenconss(); ++i )
  {
     int cons = seeed->getOpenconss()[i];
     int neighborhoodsizecons = seeedpool->getNConssForCons(cons);

     neighborhoodsize[i] = std::pair<int,int>(neighborhoodsizecons, i);
  }

  std::sort(neighborhoodsize.begin(), neighborhoodsize.end(), sort_pred() );

  for( int i = 0; i < lastindex && i < (int) neighborhoodsize.size() - 1; ++i )
  {
     if( maxdiff < neighborhoodsize[i].first - neighborhoodsize[i+1].first )
     {
        maxdiff = neighborhoodsize[i].first - neighborhoodsize[i+1].first;
        maxdiffindex = i;
     }
  }


  for( int i = 0; i <= maxdiffindex; ++i )
  {
     seeed->bookAsMasterCons(seeed->getOpenconss()[neighborhoodsize[i].second]);
  }

  decdesc << "neighborhoodmaster" << "\\_" << maxdiffindex ;

  seeed->flushBooked();
  (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, decdesc.str().c_str());
  seeed->addDetectorChainInfo(decinfo);

  foundseeeds.push_back(seeed);


  SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );

  SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), foundseeeds.size() ) );
  seeedPropagationData->nNewSeeeds  = foundseeeds.size();

  SCIPinfoMessage(scip, NULL, "dec_neighborhoodmaster found %d new seeed \n", seeedPropagationData->nNewSeeeds  );

  for( int s = 0; s < seeedPropagationData->nNewSeeeds; ++s )
  {
     seeedPropagationData->newSeeeds[s] = foundseeeds[s];
     seeedPropagationData->newSeeeds[s]->addClockTime(SCIPgetClockTime(scip, temporaryClock ) );
  }

  SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );

  *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

#define detectorPostprocessSeeedNeighborhoodmaster NULL

static
DEC_DECL_SETPARAMAGGRESSIVE(setParamAggressiveNeighborhoodmaster)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );


   return SCIP_OKAY;

}


static
DEC_DECL_SETPARAMDEFAULT(setParamDefaultNeighborhoodmaster)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLED) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLEDORIGINAL ) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLEDFINISHING ) );

   return SCIP_OKAY;

}

static
DEC_DECL_SETPARAMFAST(setParamFastNeighborhoodmaster)
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

/** creates the handler for neighborhoodmaster detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorNeighborhoodmaster(SCIP* scip /**< SCIP data structure */
)
{
   DEC_DETECTORDATA* detectordata;

   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   SCIP_CALL(
      DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND,
         DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY,
         DEC_ENABLED, DEC_ENABLEDORIGINAL, DEC_ENABLEDFINISHING,DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL,
         DEC_LEGACYMODE, detectordata, detectNeighborhoodmaster, freeNeighborhoodmaster, initNeighborhoodmaster,
         exitNeighborhoodmaster, propagateSeeedNeighborhoodmaster, NULL, NULL, finishSeeedNeighborhoodmaster,
         detectorPostprocessSeeedNeighborhoodmaster, setParamAggressiveNeighborhoodmaster,
         setParamDefaultNeighborhoodmaster, setParamFastNeighborhoodmaster));

   SCIP_CALL( SCIPaddRealParam(scip, "detection/detectors/neighborhoodmaster/maxratio",
         "the maximal ratio of open constraints that are assigned to the master problem",
         &detectordata->maxratio, FALSE, DEFAULT_MAXRATIO, 0.0, 1.0, NULL, NULL ) );

   return SCIP_OKAY;
}
