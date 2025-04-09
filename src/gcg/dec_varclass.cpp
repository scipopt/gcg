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

/**@file   dec_varclass.cpp
 * 
 * @brief  detector varclass
 * @author Julius Hense
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/dec_varclass.h"
#include "gcg/cons_decomp.h"
#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/class_varpartition.h"
#include "gcg/gcg.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "gcg/scip_misc.h"
#include "scip/clock.h"

#include <sstream>

#include <iostream>
#include <algorithm>

/* constraint handler properties */
#define DEC_NAME                  "varclass"       /**< name of detector */
#define DEC_DESC                  "detector varclass" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /**< frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          0           /**< last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /**< first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /**< frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /**< last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /**< first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               'v'         /**< display character of detector */
#define DEC_ENABLED               TRUE        /**< should the detection be enabled */
#define DEC_ENABLEDFINISHING      FALSE        /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated partialdec */

#define DEFAULT_MAXIMUMNCLASSES     8
#define AGGRESSIVE_MAXIMUMNCLASSES  10
#define FAST_MAXIMUMNCLASSES        6

#define SET_MULTIPLEFORSIZETRANSF   12500

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
#define freeVarclass NULL

/** destructor of detector to free detector data (called before the solving process begins) */
#define exitVarclass NULL

/** detection initialization function of detector (called before solving is about to begin) */
#define initVarclass NULL

#define finishPartialdecVarclass NULL

static GCG_DECL_PROPAGATEPARTIALDEC(propagatePartialdecVarclass)
{
   *result = SCIP_DIDNOTFIND;
   char decinfo[SCIP_MAXSTRLEN];
   SCIP* origprob = GCGgetOrigprob(gcg);

   SCIP_CLOCK* temporaryClock;

   if (partialdecdetectiondata->workonpartialdec->getNOpenconss() != partialdecdetectiondata->detprobdata->getNConss() ||  partialdecdetectiondata->workonpartialdec->getNOpenvars() != partialdecdetectiondata->detprobdata->getNVars() )
   {
    *result = SCIP_SUCCESS;
     return SCIP_OKAY;
   }

   SCIP_CALL_ABORT( SCIPcreateClock(origprob, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(origprob, temporaryClock) );

   std::vector<gcg::PARTIALDECOMP*> foundpartialdecs(0);

   gcg::PARTIALDECOMP* partialdecOrig;
   gcg::PARTIALDECOMP* partialdec;

   int maximumnclasses;

   if( partialdecdetectiondata->detprobdata->getNConss() + partialdecdetectiondata->detprobdata->getNVars() >= 50000 )
      SCIPgetIntParam(origprob, "detection/classification/maxnclassesperpartitionforlargeprobs", &maximumnclasses);
   else
      SCIPgetIntParam(origprob, "detection/classification/maxnclassesperpartition", &maximumnclasses);

   for( int classifierIndex = 0; classifierIndex < partialdecdetectiondata->detprobdata->getNVarPartitions(); ++classifierIndex )
   {
      gcg::VarPartition* classifier = partialdecdetectiondata->detprobdata->getVarPartition(classifierIndex);
      std::vector<int> varclassindices_master = std::vector<int>(0);
      std::vector<int> varclassindices_linking = std::vector<int>(0);

      if ( classifier->getNClasses() > maximumnclasses )
      {
         std::cout << " the current varclass distribution includes " <<  classifier->getNClasses() << " classes but only " << maximumnclasses << " are allowed for propagatePartialdec() of var class detector" << std::endl;
         continue;
      }

      partialdecOrig = partialdecdetectiondata->workonpartialdec;

      for( int i = 0; i < classifier->getNClasses(); ++ i )
      {
         switch( classifier->getClassDecompInfo(i) )
         {
            case gcg::ALL:
               break;
            case gcg::LINKING:
               varclassindices_linking.push_back(i);
               break;
            case gcg::MASTER:
               varclassindices_master.push_back(i);
               break;
            case gcg::BLOCK:
               break;
         }
      }

      std::vector< std::vector<int> > subsetsOfVarclasses = classifier->getAllSubsets( true, false, false, false );

      for( auto& subset : subsetsOfVarclasses )
      {
         if( subset.empty() && varclassindices_master.empty() && varclassindices_linking.empty() )
            continue;

         partialdec = new gcg::PARTIALDECOMP(partialdecOrig);

         /* fix open vars that have a) type of the current subset or b) decomp info LINKING as linking vars */
         auto& openvars = partialdec->getOpenvars();
         for( auto itr = openvars.cbegin(); itr != openvars.cend(); )
         {
            bool foundVar = false;
            for( int varclassId : subset )
            {
               if( classifier->getClassOfVar(*itr) == varclassId )
               {
                  itr = partialdec->fixVarToLinking(itr);
                  foundVar = true;
                  break;
               }
            }
            /* only check varclassindices_linking if current var has not already been found in a subset */
            if ( !foundVar )
            {
               for( int varclassId : varclassindices_linking )
               {
                  if( classifier->getClassOfVar(*itr) == varclassId )
                  {
                     itr = partialdec->fixVarToLinking(itr);
                     foundVar = true;
                     break;
                  }
               }
            }
            /* only check varclassindices_master if current var has not already been found in a subset */
            if ( !foundVar )
            {
               for( int varclassId : varclassindices_master )
               {
                  if( classifier->getClassOfVar(*itr) == varclassId )
                  {
                     itr = partialdec->fixVarToMaster(itr);
                     foundVar = true;
                     break;
                  }
               }
            }
            if( !foundVar )
            {
               ++itr;
            }
         }

         /* set decinfo to: varclass_<classfier_name>:<linking_class_name#1>-...-<linking_class_name#n> */
         std::stringstream decdesc;
         decdesc << "varclass" << "\\_" << classifier->getName() << ": \\\\ ";
         std::vector<int> curlinkingclasses( varclassindices_linking );
         for ( size_t varclassId = 0; varclassId < subset.size(); ++varclassId )
         {
            if ( varclassId > 0 )
            {
               decdesc << "-";
            }
            decdesc << classifier->getClassName(subset[varclassId]);

            if( std::find( varclassindices_linking.begin(), varclassindices_linking.end(),
               subset[varclassId] ) == varclassindices_linking.end() )
            {
               curlinkingclasses.push_back(subset[varclassId]);
            }
         }
         for( int varclassId : varclassindices_linking )
         {
            if( varclassId > 0 || !subset.empty() )
            {
               decdesc << "-";
            }
            decdesc << classifier->getClassName(varclassId);
         }

         partialdec->sort();
         (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, decdesc.str().c_str());
         partialdec->addDetectorChainInfo(decinfo);
         partialdec->setVarPartitionStatistics(partialdec->getNDetectors(), classifier, curlinkingclasses,
                                               varclassindices_master);

         foundpartialdecs.push_back(partialdec);
      }
   }

   SCIP_CALL_ABORT( SCIPstopClock(origprob, temporaryClock) );

   partialdecdetectiondata->detectiontime = SCIPgetClockTime(origprob, temporaryClock);
   SCIP_CALL( SCIPallocMemoryArray(origprob, &(partialdecdetectiondata->newpartialdecs), foundpartialdecs.size()) );
   partialdecdetectiondata->nnewpartialdecs = (int) foundpartialdecs.size();

   for( int s = 0; s < partialdecdetectiondata->nnewpartialdecs; ++s )
   {
      partialdecdetectiondata->newpartialdecs[s] = foundpartialdecs[s];
      partialdecdetectiondata->newpartialdecs[s]->addClockTime(partialdecdetectiondata->detectiontime / partialdecdetectiondata->nnewpartialdecs);
   }

   SCIP_CALL_ABORT(SCIPfreeClock(origprob, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


#define detectorPostprocessPartialdecVarclass NULL

static
GCG_DECL_SETPARAMAGGRESSIVE(setParamAggressiveVarclass)
{
   char setstr[SCIP_MAXSTRLEN];
   SCIP_Real modifier;

   int newval;
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE ) );

   if( SCIPgetStage(origprob) < SCIP_STAGE_PROBLEM )
      return SCIP_OKAY;

   modifier = ((SCIP_Real)SCIPgetNConss(origprob) + (SCIP_Real)SCIPgetNVars(origprob) ) / SET_MULTIPLEFORSIZETRANSF;
   modifier = log(modifier) / log(2.);

   if (!SCIPisFeasPositive(origprob, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(origprob, modifier);

   newval = MAX( 2, AGGRESSIVE_MAXIMUMNCLASSES - modifier );
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnclasses", name);

   SCIP_CALL( SCIPsetIntParam(origprob, setstr, newval ) );
   SCIPverbMessage(origprob, SCIP_VERBLEVEL_DIALOG, NULL, "\n%s = %d\n", setstr, newval);

   return SCIP_OKAY;
}


static
GCG_DECL_SETPARAMDEFAULT(setParamDefaultVarclass)
{
   char setstr[SCIP_MAXSTRLEN];
   SCIP_Real modifier;

   int newval;
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, DEC_ENABLED) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, DEC_ENABLEDFINISHING ) );

   if( SCIPgetStage(origprob) < SCIP_STAGE_PROBLEM )
         return SCIP_OKAY;

   modifier = ( (SCIP_Real)SCIPgetNConss(origprob) + (SCIP_Real)SCIPgetNVars(origprob) ) / SET_MULTIPLEFORSIZETRANSF;
   modifier = log(modifier) / log(2);

   if (!SCIPisFeasPositive(origprob, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(origprob, modifier);

   newval = MAX( 2, DEFAULT_MAXIMUMNCLASSES - modifier );
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnclasses", name);

   SCIP_CALL( SCIPsetIntParam(origprob, setstr, newval ) );
   SCIPverbMessage(origprob, SCIP_VERBLEVEL_DIALOG, NULL, "\n%s = %d\n", setstr, newval);

   return SCIP_OKAY;

}

static
GCG_DECL_SETPARAMFAST(setParamFastVarclass)
{
   char setstr[SCIP_MAXSTRLEN];
   SCIP_Real modifier;
   int newval;
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE ) );

   if( SCIPgetStage(origprob) < SCIP_STAGE_PROBLEM )
         return SCIP_OKAY;

   modifier = ( (SCIP_Real)SCIPgetNConss(origprob) + (SCIP_Real)SCIPgetNVars(origprob) ) / SET_MULTIPLEFORSIZETRANSF;

   modifier = log(modifier) / log(2);

   if (!SCIPisFeasPositive(origprob, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(origprob, modifier);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnclasses", name);

   newval = MAX( 2, FAST_MAXIMUMNCLASSES - modifier );

   SCIP_CALL( SCIPsetIntParam(origprob, setstr, newval ) );
   SCIPverbMessage(origprob, SCIP_VERBLEVEL_DIALOG, NULL, "\n%s = %d\n", setstr, newval);

   return SCIP_OKAY;
}



/*
 * detector specific interface methods
 */

/** creates the handler for varclass detector and includes it in SCIP */
SCIP_RETCODE GCGincludeDetectorVarclass(
   GCG*                 gcg                  /**< GCG data structure */
   )
{
   GCG_DETECTORDATA* detectordata;
   char setstr[SCIP_MAXSTRLEN];
   SCIP* origprob = GCGgetOrigprob(gcg);

   /**@todo create varclass detector data here*/
   detectordata = NULL;

   SCIP_CALL(
      GCGincludeDetector(gcg, DEC_NAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND,
                         DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata,
                         freeVarclass, initVarclass, exitVarclass, propagatePartialdecVarclass, finishPartialdecVarclass, detectorPostprocessPartialdecVarclass, setParamAggressiveVarclass, setParamDefaultVarclass, setParamFastVarclass));

   /**@todo add varclass detector parameters */

   const char* name = DEC_NAME;
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnclasses", name);
   SCIP_CALL( SCIPaddIntParam(origprob, setstr, "maximum number of classes ",  NULL, FALSE, DEFAULT_MAXIMUMNCLASSES, 1, INT_MAX, NULL, NULL ) );

   return SCIP_OKAY;
}
