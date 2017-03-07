/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2015 Operations Research, RWTH Aachen University       */
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

/**@file   dec_consclass.cpp
 * @ingroup DETECTORS
 * @brief  detector consclass (put your description here)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_consclass.h"
#include "cons_decomp.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "gcg.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "scip_misc.h"
#include "scip/clock.h"

#include <iostream>

/* constraint handler properties */
#define DEC_DETECTORNAME          "consclass"       /**< name of detector */
#define DEC_DESC                  "detector consclass" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          0           /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               'c'         /**< display character of detector */
#define DEC_ENABLED               TRUE        /**< should the detection be enabled */
#define DEC_ENABLEDORIGINAL       TRUE        /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING      FALSE        /**< should the detection be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */

#define DEFAULT_MAXIMUMNCLASSES     8
#define AGGRESSIVE_MAXIMUMNCLASSES  10
#define FAST_MAXIMUMNCLASSES        6

#define SET_MULTIPLEFORSIZETRANSF   12500

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

/** method to enumerate all subsets */
std::vector< std::vector<int> > getAllSubsets(std::vector<int> set)
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
#define freeConsclass NULL

/** destructor of detector to free detector data (called before the solving process begins) */
#if 0
static
DEC_DECL_EXITDETECTOR(exitConsclass)
{ /*lint --e{715}*/

   SCIPerrorMessage("Exit function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define exitConsclass NULL
#endif

/** detection initialization function of detector (called before solving is about to begin) */
#if 0
static
DEC_DECL_INITDETECTOR(initConsclass)
{ /*lint --e{715}*/

   SCIPerrorMessage("Init function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define initConsclass NULL
#endif

/** detection function of detector */
static DEC_DECL_DETECTSTRUCTURE(detectConsclass)
{ /*lint --e{715}*/
   *result = SCIP_DIDNOTFIND;

   SCIPerrorMessage("Detection function of detector <%s> not implemented!\n", DEC_DETECTORNAME)
;   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


#define finishSeeedConsclass NULL

static DEC_DECL_PROPAGATESEEED(propagateSeeedConsclass)
{
  *result = SCIP_DIDNOTFIND;

  SCIP_CLOCK* temporaryClock;

  if (seeedPropagationData->seeedToPropagate->getNOpenconss() != seeedPropagationData->seeedpool->getNConss() ||  seeedPropagationData->seeedToPropagate->getNOpenvars() != seeedPropagationData->seeedpool->getNVars() )
  {
    *result = SCIP_SUCCESS;
     return SCIP_OKAY;
  }

  SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );
  SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );

  std::vector<gcg::Seeed*> foundseeeds(0);

  gcg::Seeed* seeedOrig;
  gcg::Seeed* seeed;

  int maximumnclasses;

  SCIPgetIntParam(scip, "detectors/consclass/maxnclasses", &maximumnclasses); /* if  distribution of classes exceed this number its skipped */

  for( int conssclass = 0; conssclass < seeedPropagationData->seeedpool->getNConssClassDistributions(); ++conssclass )
  {
    int nclasses = seeedPropagationData->seeedpool->getNClassesOfDistribution(conssclass);
    std::vector<int> classforcons = seeedPropagationData->seeedpool->getConssClassDistributionVector(conssclass);
    std::vector<int> consclassindices = std::vector<int>(0);

    /** check if there are to  many classes in this distribution and skip it if so */

    if ( nclasses > maximumnclasses )
    {
       std::cout << " the current consclass distribution includes " <<  nclasses << " classes but only " << maximumnclasses << " are allowed for propagateSeeed() of cons class detector" << std::endl;
       continue;
    }

  seeedOrig = new gcg::Seeed(seeedPropagationData->seeedToPropagate, seeedPropagationData->seeedpool);
  seeedOrig->setDetectorPropagated(detector);

  if(!seeedOrig->areOpenVarsAndConssCalculated())
  {
      seeedOrig->calcOpenconss();
      seeedOrig->calcOpenvars();
      seeedOrig->setOpenVarsAndConssCalculated(true);
  }

  for( int i = 0; i < nclasses; ++ i)
      consclassindices.push_back(i);

  std::vector< std::vector<int> > subsetsOfConsclasses = getAllSubsets(consclassindices);




  for(size_t subset = 0; subset < subsetsOfConsclasses.size(); ++subset)
  {
      if(subsetsOfConsclasses[subset].size() == 0)
          continue;

      seeed = new gcg::Seeed(seeedOrig, seeedPropagationData->seeedpool);
         /** set open cons that have type of the current subset to Master */
      for( int i = 0; i < seeed->getNOpenconss(); ++i)
      {
          for(size_t consclassId = 0; consclassId < subsetsOfConsclasses[subset].size(); ++consclassId )
          {
              if( classforcons[seeed->getOpenconss()[i]] == subsetsOfConsclasses[subset][consclassId] )
              {
                  seeed->bookAsMasterCons(seeed->getOpenconss()[i]);
                  break;
              }
          }
      }
      seeed->flushBooked();

      foundseeeds.push_back(seeed);
  }
  delete seeedOrig;
 }

  SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );

  SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), foundseeeds.size() ) );
  seeedPropagationData->nNewSeeeds = foundseeeds.size();

  for( int s = 0; s < seeedPropagationData->nNewSeeeds; ++s )
  {
     seeedPropagationData->newSeeeds[s] = foundseeeds[s];
     seeedPropagationData->newSeeeds[s]->addClockTime(SCIPclockGetTime(temporaryClock )  );
  }

  SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );

  *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

static
DEC_DECL_SETPARAMAGGRESSIVE(setParamAggressiveConsclass)
{
   char setstr[SCIP_MAXSTRLEN];
   SCIP_Real modifier;

   int newval;
   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );


   modifier = ((SCIP_Real)SCIPgetNConss(scip) + (SCIP_Real)SCIPgetNVars(scip) ) / SET_MULTIPLEFORSIZETRANSF;
   modifier = log(modifier) / log(2.);

   if (!SCIPisFeasPositive(scip, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(scip, modifier);

   newval = MAX( 2, AGGRESSIVE_MAXIMUMNCLASSES - modifier );
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/maxnclasses", name);

   SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );
   SCIPinfoMessage(scip, NULL, "\n%s = %d\n", setstr, newval);


   return SCIP_OKAY;

}


static
DEC_DECL_SETPARAMDEFAULT(setParamDefaultConsclass)
{
   char setstr[SCIP_MAXSTRLEN];
   SCIP_Real modifier;

   int newval;
   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE ) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );

   modifier = ( (SCIP_Real)SCIPgetNConss(scip) + (SCIP_Real)SCIPgetNVars(scip) ) / SET_MULTIPLEFORSIZETRANSF;
   modifier = log(modifier) / log(2);

   if (!SCIPisFeasPositive(scip, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(scip, modifier);

   newval = MAX( 2, DEFAULT_MAXIMUMNCLASSES - modifier );
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/maxnclasses", name);

   SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );
   SCIPinfoMessage(scip, NULL, "\n%s = %d\n", setstr, newval);

   return SCIP_OKAY;

}

static
DEC_DECL_SETPARAMFAST(setParamFastConsclass)
{
   char setstr[SCIP_MAXSTRLEN];
   SCIP_Real modifier;
   int newval;

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );

   modifier = ( (SCIP_Real)SCIPgetNConss(scip) + (SCIP_Real)SCIPgetNVars(scip) ) / SET_MULTIPLEFORSIZETRANSF;

   modifier = log(modifier) / log(2);

   if (!SCIPisFeasPositive(scip, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(scip, modifier);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/maxnclasses", name);

   newval = MAX( 2, FAST_MAXIMUMNCLASSES - modifier );

   SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );
   SCIPinfoMessage(scip, NULL, "\n%s = %d\n", setstr, newval);

   return SCIP_OKAY;

}



/*
 * detector specific interface methods
 */

/** creates the handler for consclass detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorConsclass(SCIP* scip /**< SCIP data structure */
)
{
   DEC_DETECTORDATA* detectordata;
   char setstr[SCIP_MAXSTRLEN];

   /**@todo create consclass detector data here*/
   detectordata = NULL;

   SCIP_CALL(
      DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND,
         DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDORIGINAL, DEC_ENABLEDFINISHING, DEC_SKIP, DEC_USEFULRECALL, detectordata, detectConsclass,
         freeConsclass, initConsclass, exitConsclass, propagateSeeedConsclass, finishSeeedConsclass, setParamAggressiveConsclass, setParamDefaultConsclass, setParamFastConsclass));

   /**@todo add consclass detector parameters */

   const char* name = DEC_DETECTORNAME;
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/maxnclasses", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, "maximum number of classes ",  NULL, FALSE, DEFAULT_MAXIMUMNCLASSES, 1, INT_MAX, NULL, NULL ) );

   return SCIP_OKAY;
}
