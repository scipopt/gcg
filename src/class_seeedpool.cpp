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

/**@file   class_seeedpool.cpp
 * @brief  class with functions for seeedpool
 * @author Michael Bastubbe
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg.h"
#include "objscip/objscip.h"
#include "class_seeedpool.h"
#include "struct_detector.h"
#include "struct_decomp.h"
#include "cons_decomp.h"
#include "decomp.h"
#include "scip_misc.h"
#include "scip/clock.h"

#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <queue>

#include <exception>

#if defined(_WIN32) || defined(_WIN64)
#define LINEBREAK "\r\n"
#else
#define LINEBREAK "\n"
#endif



#define SCIP_CALL_EXC(x)   do                                                                                  \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) !=  SCIP_OKAY )                                                \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             throw std::exception();                                                          \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )



/** constraint handler data */
struct SCIP_ConshdlrData
{
   DEC_DECOMP**          decdecomps;         /**< array of decomposition structures */
   DEC_DETECTOR**        detectors;          /**< array of structure detectors */
   int*                  priorities;         /**< priorities of the detectors */
   int                   ndetectors;         /**< number of detectors */
   SCIP_CLOCK*           detectorclock;      /**< clock to measure detection time */
   SCIP_Bool             hasrun;             /**< flag to indicate whether we have already detected */
   int                   ndecomps;           /**< number of decomposition structures  */
   SCIP_Bool             createbasicdecomp;  /**< indicates whether to create a decomposition with all constraints in the master if no other specified */
};

namespace gcg {

/** local methods */

int calcLevenshteinDistance(std::string s, std::string t)
{
    // trivial cases
    if (s.compare(t) == 0) return 0;
    if (s.length() == 0) return t.length();
    if (t.length() == 0) return s.length();

    // create two work vectors of integer distances
    std::vector<int> v0 (t.length() + 1);
    std::vector<int> v1 (t.length() + 1);

    // initialize v0 (the previous row of distances)
    // this row is A[0][i]: edit distance for an empty s
    // the distance is just the number of characters to delete from t
    for ( size_t i = 0; i < v0.size(); i++ )
        v0[i] = i;

    for (size_t i = 0; i < s.length(); i++)
    {
        // calculate v1 (current row distances) from the previous row v0

        // first element of v1 is A[i+1][0]
        //   edit distance is delete (i+1) chars from s to match empty t
        v1[0] = i + 1;

        // use formula to fill in the rest of the row
        for (size_t j = 0; j < t.length(); j++)
        {
            int cost = (s[i] == t[j]) ? 0 : 1;
            v1[j + 1] = std::min(v1[j] + 1, std::min( v0[j + 1] + 1, v0[j] + cost ) );
        }

        // copy v1 (current row) to v0 (previous row) for next iteration
        for (size_t j = 0; j < v0.size(); j++)
            v0[j] = v1[j];
    }

    return v1[t.length()];
}


void removeDigits(char *str, int *nremoved) {

    char digits[11] = "0123456789";
    *nremoved = 0;

    for(int i = 0; i < 10; ++i )
    {
       char digit = digits[i];
       size_t j = 0;
       while ( j < strlen(str) )
       {
          if (str[j] == digit)
          {
             *nremoved = *nremoved + 1;
             for(size_t k = j; k < strlen(str); ++k)
             {
                str[k] = str[k+1];
             }
          }
          else ++j;
       }
    }
}


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

/** methode to calculate the greates common divisor */

int gcd(int a, int b) {
    return b == 0 ? a : gcd(b, a % b);
}


SCIP_CONS* consGetRelevantRepr(SCIP* scip, SCIP_CONS* cons){

   return cons;
}

SCIP_VAR* varGetRelevantRepr(SCIP* scip, SCIP_VAR* var){

        return SCIPvarGetProbvar(var);
}

SCIP_Bool seeedIsNoDuplicateOfSeeeds(SeeedPtr compseeed, std::vector<SeeedPtr> const & seeeds, bool sort){

   assert(compseeed != NULL);
   SCIP_Bool isDuplicate;

   for( size_t i = 0; i < seeeds.size(); ++i )
   {
      assert(seeeds[i] != NULL);

      compseeed->isEqual(seeeds[i], &isDuplicate, sort );
      if ( isDuplicate )
         return FALSE;
   }
   return TRUE;
}

SCIP_Bool seeedIsNoDuplicate(SeeedPtr seeed, std::vector<SeeedPtr> const & currSeeeds, std::vector<SeeedPtr> const & finishedSeeeds, bool sort){

   SCIP_Bool bool1 = seeedIsNoDuplicateOfSeeeds(seeed, currSeeeds, sort);
   SCIP_Bool bool2 = seeedIsNoDuplicateOfSeeeds(seeed, finishedSeeeds, sort);
   return ( bool1 && bool2 );
}


/** constructor */
 Seeedpool::Seeedpool(
    SCIP*               givenScip, /**< SCIP data structure */
        const char*             conshdlrName
    ):scip(givenScip), currSeeeds(0), nTotalSeeeds(0),nVars(SCIPgetNVars(givenScip) ), nConss(SCIPgetNConss(givenScip) ), nDetectors(0), ndecompositions(0), candidatesNBlocks(0)
 {
         SCIP_CONS** conss;
         SCIP_VAR** vars;

         SCIP_CONSHDLR* conshdlr;  /** cons_decomp to get detectors */
         SCIP_CONSHDLRDATA* conshdlrdata;


         int relevantVarCounter = 0;
         int relevantConsCounter = 0;

         /** store all enabled detectors */

         conshdlr = SCIPfindConshdlr(scip, conshdlrName);
         assert(conshdlr != NULL);
         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert(conshdlrdata != NULL);

         /** set detection data */
         SCIP_CALL_ABORT( SCIPgetIntParam(givenScip, "detection/maxrounds", &maxndetectionrounds) );

         /** store priorities of the detectors */
         for(int d = 0; d < conshdlrdata->ndetectors; ++d )
         {
                 DEC_DETECTOR *detector;
                 detector = conshdlrdata->detectors[d];
                 assert(detector != NULL);
                 conshdlrdata->priorities[d] = detector->priority;
         }

         SCIPdebugMessage("Sorting %i detectors\n", conshdlrdata->ndetectors);

         /** sort the detectors according their priorities */
         SCIPsortIntPtr(conshdlrdata->priorities, (void**)conshdlrdata->detectors, conshdlrdata->ndetectors);

         SCIPdebugMessage("Trying %d detectors.\n", conshdlrdata->ndetectors);

         for(int d = 0; d < conshdlrdata->ndetectors; ++d )
         {
                 DEC_DETECTOR* detector;

                 detector = conshdlrdata->detectors[d];
                 assert(detector != NULL);
                 if( !detector->enabled || detector->propagateSeeed == NULL)
                         continue;

                 scipDetectorToIndex[detector] = nDetectors;
                 detectorToScipDetector.push_back(detector);
                 ++nDetectors;

         }

         /** initilize matrix datastructures */
         conss = SCIPgetConss(scip);
         vars = SCIPgetVars(scip);

         /** assign an index to every cons and var
          * @TODO: are all constraints/variables relevant? (probvars etc)  */
         for(int i = 0; i < nConss; ++i)
         {
                 SCIP_CONS* relevantCons;

                 relevantCons = consGetRelevantRepr(scip, conss[i]);
                 if( relevantCons != NULL )
                 {
                         scipConsToIndex[relevantCons] = relevantConsCounter ;
                         consToScipCons.push_back(relevantCons);
                         ++relevantConsCounter;
                 }
         }

         for(int i = 0; i < nVars; ++i)
         {
                 SCIP_VAR* relevantVar;

                 relevantVar = varGetRelevantRepr(scip, vars[i]);

                 if( relevantVar != NULL )
                 {
                         scipVarToIndex[relevantVar] = relevantVarCounter ;
                         varToScipVar.push_back(relevantVar);
                         ++relevantVarCounter;
                 }
         }

         /** from here on nVars and nConss represents the relevant numbers */
         nVars = relevantVarCounter;
         nConss = relevantConsCounter;
         varsForConss = std::vector<std::vector<int>>(nConss);
         valsForConss = std::vector<std::vector<SCIP_Real>>(nConss);
         conssForVars = std::vector<std::vector<int>>(nVars);

         assert((int) varToScipVar.size() == nVars);
         assert((int) consToScipCons.size() == nConss);

         /** assumption: now every relevant constraint and variable has its index and is stored in the corresponding unordered_map */
         /** find constraint <-> variable relationships and store them in both directions */
         for( int i = 0; i < (int)consToScipCons.size() ; ++i )
         {
                 SCIP_CONS* cons;
                 SCIP_VAR** currVars;
                 SCIP_Real* currVals;
                 int            nCurrVars;

                 cons = consToScipCons[i];

                 nCurrVars = GCGconsGetNVars(scip, cons);

                 SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &currVars, nCurrVars) ); /** free in line 321 */
                 SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &currVals, nCurrVars) ); /** free in line 321 */
                 SCIP_CALL_ABORT(GCGconsGetVars(scip, cons, currVars, nCurrVars));
                 SCIP_CALL_ABORT(GCGconsGetVals(scip, cons, currVals, nCurrVars));

                 for(int currVar = 0; currVar < nCurrVars; ++currVar)
                 {
                     int varIndex;
                     std::tr1::unordered_map<SCIP_VAR*, int>::const_iterator iterVar;

                     /** because of the bug of GCGconsGet*()-methods some variables have to be negated */
                     if(!SCIPvarIsNegated(currVars[currVar]))
                        iterVar = scipVarToIndex.find(currVars[currVar]);
                     else
                        iterVar = scipVarToIndex.find(SCIPvarGetNegatedVar(currVars[currVar]));

                     if(iterVar == scipVarToIndex.end() )
                        continue;


                         varIndex = iterVar->second;

                         varsForConss[i].push_back(varIndex);
                         conssForVars[varIndex].push_back(i);
                         valsForConss[i].push_back(currVals[currVar]);
                         valsMap[std::pair<int,int>(i, varIndex)] =  currVals[currVar] ;
                 }
                 SCIPfreeBufferArrayNull(scip, &currVars);
                 SCIPfreeBufferArrayNull(scip, &currVals);
         }

         /* populate seeedpool with empty seeed */

         currSeeeds.push_back(new Seeed(scip, nTotalSeeeds,nDetectors,nConss,nVars) );

         nTotalSeeeds++;

         decompositions = NULL;

         addConssClassesForSCIPConstypes();
         addConssClassesForConsnamesDigitFreeIdentical();
         addConssClassesForConsnamesLevenshteinDistanceConnectivity(1);
         calcCandidatesNBlocks();


 }//end constructor

 Seeedpool::~Seeedpool(){

 }


 /** finds decompositions  */
  /** access coefficient matrlix constraint-wise */
 void    Seeedpool::findDecompositions(
 ){

         /** 1) read parameter, as there are: maxrounds
          *  2) loop rounds
          *  3) every seeed in seeeds
          *  4) every detector not registered yet propagates seeed
          *  5)  */


         SEEED_PROPAGATION_DATA* seeedPropData;
//         SCIP_VAR* probvar;
//         SCIP_CONS* cons;
//         int cindex = 0;
//         int vindex = 0;
//         int currblock;
         bool displaySeeeds = false;
         int verboseLevel;
         std::vector<int> successDetectors;
         std::vector<SeeedPtr> delSeeeds;
         bool duplicate;

//         SCIP_CLOCK* temporaryClock; /* @TODO replace with finer measurement in detectors */
//
//         SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );

         successDetectors = std::vector<int>(nDetectors, 0);
         ndecompositions = 0;
         seeedPropData = new SEEED_PROPAGATION_DATA();
         seeedPropData->seeedpool = this;
         seeedPropData->nNewSeeeds = 0;
         delSeeeds = std::vector<SeeedPtr>(0);

         verboseLevel = 0;

         for(size_t s = 0; s < currSeeeds.size(); ++s)
         {
            currSeeeds[s]->sort();
            currSeeeds[s]->considerImplicits(this);
            currSeeeds[s]->calcHashvalue();
         }

         for(int round = 0; round < maxndetectionrounds; ++round)
         {
                 std::cout << "currently in detection round " << round << std::endl;
                 std::vector<SeeedPtr> nextSeeeds = std::vector<SeeedPtr>(0);
                 std::vector<SeeedPtr> currSeeedsToDelete = std::vector<SeeedPtr>(0);

                 for(size_t s = 0; s < currSeeeds.size(); ++s )
                 {
                         SeeedPtr seeedPtr;
                         seeedPtr= currSeeeds[s];
                         if(displaySeeeds)
                         {
                            std::cout << "Start to propagate seeed " << seeedPtr->getID() << " in round " << round << ":" << std::endl;
                            seeedPtr->displaySeeed();
                         }

                         /** the current seeed is handled by all detectors */
                         for(int d = 0; d < nDetectors; ++d)
                         {


                                 DEC_DETECTOR* detector;
                                 std::vector<SeeedPtr>::const_iterator newSIter;
                                 std::vector<SeeedPtr>::const_iterator newSIterEnd;


                                 SCIP_RESULT result = SCIP_DIDNOTFIND;
                                 detector = detectorToScipDetector[d];

                                 /** if the seeed is also propagated by the detector go on with the next detector */
                                 if(seeedPtr->isPropagatedBy(d) && !detector->usefulRecall )
                                         continue;

                                 /** check if detector is callable in current detection round */
                                 if(detector->maxCallRound < round || detector->minCallRound > round)
                                     continue;

                                 if( (round - detector->minCallRound) % detector->freqCallRound != 0 )
                                     continue;

                                 seeedPropData->seeedToPropagate = seeedPtr;

                                 /** new seeeds are created by the current detector */
                                 SCIP_CALL_ABORT( SCIPstartClock(scip, detectorToScipDetector[d]->dectime) );
                                 //SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
                                 if(verboseLevel > 2)
                                     std::cout << "detector " << DECdetectorGetName(detectorToScipDetector[d]) << " started to propagate the " << s+1 << ". seeed (ID " << seeedPtr->getID() << ") in round " << round+1 << std::endl;

                                 SCIP_CALL_ABORT(detectorToScipDetector[d]->propagateSeeed(scip, detectorToScipDetector[d],seeedPropData, &result) );



                                 for( int j = 0; j < seeedPropData->nNewSeeeds; ++j )
                                 {
                                    seeedPropData->newSeeeds[j]->considerImplicits(this);
                                    seeedPropData->newSeeeds[j]->sort();
                                    if(!seeedPropData->newSeeeds[j]->checkConsistency())
                                    {
                                        seeedPropData->newSeeeds[j]->displaySeeed();
                                        assert(false);
                                    }
                                    seeedPropData->newSeeeds[j]->calcHashvalue();
                                    seeedPropData->newSeeeds[j]->addDecChangesFromAncestor(seeedPtr);
                                    //seeedPropData->newSeeeds[j]->addClockTime( SCIPclockGetTime(temporaryClock )  );
                                 }

                                 SCIP_CALL_ABORT( SCIPstopClock(scip, detectorToScipDetector[d]->dectime) );
                                 //SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
                                 //SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );


                                 if(seeedPropData->nNewSeeeds != 0 && (displaySeeeds ) )
                                 {
                                    std::cout << "detector " << DECdetectorGetName(detectorToScipDetector[d] ) << " found " << seeedPropData->nNewSeeeds << " new seeed(s): ";
                                    std::cout << seeedPropData->newSeeeds[0]->getID();
                                    for( int j = 1; j < seeedPropData->nNewSeeeds; ++j )
                                       std::cout << ", " << seeedPropData->newSeeeds[j]->getID();
                                    std::cout << "\n";

                                    if(displaySeeeds)
                                    {
                                       for( int j = 0; j < seeedPropData->nNewSeeeds; ++j )
                                          seeedPropData->newSeeeds[j]->displaySeeed();
                                    }
                                 }
                                 else
                                     if(displaySeeeds)
                                         std::cout << "detector " << DECdetectorGetName(detectorToScipDetector[d] ) << " found 0 new seeeds" << std::endl;





                                 /** if the new seeeds are no duplicate they're added to the currSeeeds */
                                 for( int seeed = 0; seeed < seeedPropData->nNewSeeeds; ++seeed )
                                 {
                                         if( !seeedPropData->newSeeeds[seeed]->isTrivial() && seeedIsNoDuplicate(seeedPropData->newSeeeds[seeed], nextSeeeds, finishedSeeeds, false) )
                                         {
                                            seeedPropData->newSeeeds[seeed]->calcOpenconss();
                                            seeedPropData->newSeeeds[seeed]->calcOpenvars();
                                            if(seeedPropData->newSeeeds[seeed]->getNOpenconss() == 0 && seeedPropData->newSeeeds[seeed]->getNOpenvars() == 0)
                                            {
                                               if(verboseLevel > 2)
                                               {
                                                   std::cout << "seeed " << seeedPropData->newSeeeds[seeed]->getID() << " is addded to finished seeeds!" << std::endl;
                                                   seeedPropData->newSeeeds[seeed]->showScatterPlot(this);
                                               }
                                                   finishedSeeeds.push_back(seeedPropData->newSeeeds[seeed]);
                                            }
                                            else
                                            {
                                               if(verboseLevel > 2)
                                               {
                                                   std::cout << "seeed " << seeedPropData->newSeeeds[seeed]->getID() << " is addded to next round seeeds!" << std::endl;
                                                   seeedPropData->newSeeeds[seeed]->showScatterPlot(this);
                                               }
                                               nextSeeeds.push_back(seeedPropData->newSeeeds[seeed]);
                                            }
                                         }
                                         else
                                         {
                                            delete seeedPropData->newSeeeds[seeed];
                                            seeedPropData->newSeeeds[seeed] = NULL;
                                         }
                                 }
                                 /** cleanup propagation data structure */
                                 SCIPfreeMemoryArrayNull(scip, &seeedPropData->newSeeeds);
                                 seeedPropData->newSeeeds = NULL;
                                 seeedPropData->nNewSeeeds = 0;
                         }

//                         SCIP_CALL_ABORT(seeedPtr->completeByConnected( seeedPropData->seeedpool ) );
                         for(int d = 0; d < nDetectors; ++d)
                         {
                            DEC_DETECTOR* detector = detectorToScipDetector[d];
                            SCIP_RESULT result = SCIP_DIDNOTFIND;
                            seeedPropData->seeedToPropagate = seeedPtr;

                            /** if the finishing of the detector is not enabled go on with the next detector */
                            if( !detector->enabledFinishing )
                                    continue;

                            SCIP_CALL_ABORT(detectorToScipDetector[d]->finishSeeed(scip, detectorToScipDetector[d],seeedPropData, &result) );

                            for(int finished = 0; finished < seeedPropData->nNewSeeeds; ++finished)
                            {
                               SeeedPtr seeed = seeedPropData->newSeeeds[finished];
                               seeed->calcHashvalue();
                               if( seeedIsNoDuplicateOfSeeeds(seeed, finishedSeeeds, false) )
                               {
                                  finishedSeeeds.push_back(seeed);
                               }
                               else
                               {
                                  bool isIdentical = false;
                                  for (size_t h = 0; h < finishedSeeeds.size(); ++h )
                                  {
                                     if( seeed == finishedSeeeds[h] )
                                     {
                                        isIdentical = true;
                                        break;
                                     }
                                  }

                                  if( !isIdentical )
                                  {
                                     currSeeedsToDelete.push_back(seeed);
                                  }
                               }
                            }
                            SCIPfreeMemoryArrayNull(scip, &seeedPropData->newSeeeds);
                            seeedPropData->newSeeeds = NULL;
                            seeedPropData->nNewSeeeds = 0;
                         }
//                         seeedPtr->calcHashvalue();


//                   if( seeedIsNoDuplicateOfSeeeds(seeedPtr, finishedSeeeds, false) )
//                   {
//                      finishedSeeeds.push_back(seeedPtr);
//                   }
//                   else
//                   {
//                      bool isIdentical = false;
//                      for (size_t h = 0; h < finishedSeeeds.size(); ++h )
//                      {
//                         if( seeedPtr == finishedSeeeds[h] )
//                         {
//                            isIdentical = true;
//                            break;
//                         }
//                      }
//
//                      if( !isIdentical )
//                      {
//                         currSeeedsToDelete.push_back(seeedPtr);
//                      }
//                   }
                 }

                 for(size_t s = 0; s < currSeeedsToDelete.size(); ++s )
                    delete currSeeedsToDelete[s];

                 currSeeeds = nextSeeeds;

         }

//         /* completeByconnected() on  currseeeds (from last round) and add them to finished seeeds */
//
//         for(size_t i = 0; i < currSeeeds.size(); ++i)
//         {
//             SeeedPtr seeedPtr = currSeeeds[i];
//
//             SCIP_CALL_ABORT(seeedPtr->completeByConnected( seeedPropData->seeedpool ) );
//             seeedPtr->calcHashvalue();
//             /* currseeeds are freed later */
//             if(seeedIsNoDuplicateOfSeeeds(seeedPtr, finishedSeeeds, false))
//             {
//                if(verboseLevel > 2)
//                {
//                   std::cout << "seeed " << seeedPtr->getID() << " is finished from next round seeeds!" << std::endl;
//                   seeedPtr->showScatterPlot(this);
//                }
//                finishedSeeeds.push_back(seeedPtr);
//             }
//
//         }

         /** complete the currseeeds with finishing detectors and add them to finished seeeds */
         for(size_t i = 0; i < currSeeeds.size(); ++i)
         {
            SeeedPtr seeedPtr = currSeeeds[i];
            for(int d = 0; d < nDetectors; ++d)
            {
               DEC_DETECTOR* detector = detectorToScipDetector[d];
               SCIP_RESULT result = SCIP_DIDNOTFIND;
               seeedPropData->seeedToPropagate = seeedPtr;

               /** if the finishing of the detector is not enabled go on with the next detector */
               if( !detector->enabledFinishing )
                  continue;

               SCIP_CALL_ABORT(detectorToScipDetector[d]->finishSeeed(scip, detectorToScipDetector[d],seeedPropData, &result) );

               for(int finished = 0; finished < seeedPropData->nNewSeeeds; ++finished)
               {
                  SeeedPtr seeed = seeedPropData->newSeeeds[finished];
                  seeed->calcHashvalue();
                  if( seeedIsNoDuplicateOfSeeeds(seeed, finishedSeeeds, false) )
                  {
                     if(verboseLevel > 2)
                     {
                        std::cout << "seeed " << seeed->getID() << " is finished from next round seeeds!" << std::endl;
                        seeed->showScatterPlot(this);
                     }
                     finishedSeeeds.push_back(seeed);
                  }

                  SCIPfreeMemoryArrayNull(scip, &seeedPropData->newSeeeds);
                  seeedPropData->newSeeeds = NULL;
                  seeedPropData->nNewSeeeds = 0;
               }
            }
         }

         std::cout << (int) finishedSeeeds.size() << " finished seeeds are found." << std::endl;



         if(displaySeeeds)
         {
            for(size_t i = 0; i < finishedSeeeds.size(); ++i)
            {
               std::cout << i+1 << "th finished seeed: " << std::endl;
               finishedSeeeds[i]->displaySeeed();
            }
         }

         /** count the successful refinement calls for each detector */

         for(size_t i = 0; i < finishedSeeeds.size(); ++i)
         {
            assert(finishedSeeeds[i]->checkConsistency() );
            assert(finishedSeeeds[i]->getNOpenconss() == 0);
            assert(finishedSeeeds[i]->getNOpenvars() == 0);


             for( int d = 0; d < nDetectors; ++d )
             {
                 if(finishedSeeeds[i]->isPropagatedBy(d))
                     successDetectors[d] += 1;
             }
         }

         /** preliminary output detector stats */

         std::cout << "Begin preliminary detector times: " << std::endl;

         for( int i = 0; i < nDetectors; ++i )
         {
             std::cout << "Detector " << std::setw(25) << std::setiosflags(std::ios::left) << DECdetectorGetName(detectorToScipDetector[i] ) << " \t worked on \t " << successDetectors[i] << " of " << finishedSeeeds.size() << "\t and took a total time of \t" << SCIPgetClockTime(scip, detectorToScipDetector[i]->dectime)  << std::endl;
         }

         if((int)finishedSeeeds.size() != 0)
         {
            SCIP_Real minscore = finishedSeeeds[0]->evaluate(this);
            SeeedPtr bestSeeed = finishedSeeeds[0];
            for( size_t i = 1; i < finishedSeeeds.size(); ++i )
            {
               SCIP_Real score = finishedSeeeds[i]->evaluate(this);
               if (score < minscore)
               {
                  minscore = score;
                  bestSeeed = finishedSeeeds[i];
               }
            }
            bestSeeed->showScatterPlot(this);
         }


         /** fill out the decompositions */

         SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &decompositions, (int) finishedSeeeds.size())); /** free in decomp.c:470 */
         for( size_t i = 0; i < finishedSeeeds.size(); ++i )
         {
            SeeedPtr seeed = finishedSeeeds[i];

            SCIP_HASHMAP* vartoblock;
            SCIP_HASHMAP* constoblock;
            SCIP_HASHMAP* varindex;
            SCIP_HASHMAP* consindex;

            SCIP_VAR*** stairlinkingvars;
            SCIP_VAR*** subscipvars;
            SCIP_VAR**  linkingvars;
            SCIP_CONS**  linkingconss;
            SCIP_CONS*** subscipconss;

            int* nsubscipconss;
            int* nsubscipvars;
            int* nstairlinkingvars;
            int  nlinkingvars;

            int varcounter = 1;  /* in varindex counting starts with 1 */
            int conscounter = 1; /* in consindex counting starts with 1 */
            int counterstairlinkingvars = 0;

            int size;

            assert(seeed->checkConsistency() );

            /* create decomp data structure */
            SCIP_CALL_ABORT( DECdecompCreate(scip, &(decompositions[i])) );

 //           seeed->displayConss();
 //           seeed->showScatterPlot(this);


            /** set nblocks */
            DECdecompSetNBlocks(decompositions[i], seeed->getNBlocks() );

            /** set constraints */
            if( seeed->getNMasterconss( )  != 0 )
               SCIP_CALL_ABORT (SCIPallocBufferArray(scip, &linkingconss, seeed->getNMasterconss() ) );
            else  linkingconss = NULL;

            SCIP_CALL_ABORT (SCIPallocBufferArray(scip, &nsubscipconss, seeed->getNBlocks() ) );
            SCIP_CALL_ABORT (SCIPallocBufferArray(scip, &subscipconss, seeed->getNBlocks() ) );

            SCIP_CALL_ABORT( SCIPhashmapCreate( &constoblock, SCIPblkmem(scip), seeed->getNConss() ) );
            SCIP_CALL_ABORT( SCIPhashmapCreate( &consindex, SCIPblkmem(scip), seeed->getNConss() ) );

            /* set linking constraints */
            for (int c = 0; c < seeed->getNMasterconss() ; ++c)
            {
               int consid = seeed->getMasterconss()[c];
               SCIP_CONS* scipcons = consToScipCons[consid];
               linkingconss[c] = scipcons;
               SCIP_CALL_ABORT( SCIPhashmapInsert(constoblock, scipcons, (void*) (size_t) (seeed->getNBlocks() + 1) ) );
               SCIP_CALL_ABORT( SCIPhashmapInsert(consindex, scipcons, (void*) (size_t) conscounter) );
               conscounter++;
            }

            if (seeed->getNMasterconss() != 0 )
               DECdecompSetLinkingconss(scip, decompositions[i], linkingconss, seeed->getNMasterconss());
            else
               linkingconss = NULL;
            /* set block constraints */
            for ( int b = 0; b < seeed->getNBlocks(); ++b )
            {
               SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &subscipconss[b], seeed->getNConssForBlock(b) ) );
               nsubscipconss[b] = seeed->getNConssForBlock(b);
               for ( int c = 0; c < seeed->getNConssForBlock(b); ++c )
               {
                  int consid  = seeed->getConssForBlock(b)[c];
                  SCIP_CONS* scipcons = consToScipCons[consid];

                  assert(scipcons != NULL);
                  subscipconss[b][c] = scipcons;
                  SCIP_CALL_ABORT( SCIPhashmapInsert(constoblock, scipcons, (void*) (size_t) (b + 1 ) ) ) ;
                  SCIP_CALL_ABORT( SCIPhashmapInsert(consindex, scipcons, (void*) (size_t) conscounter) );
                  conscounter++;
               }
            }


            DECdecompSetSubscipconss(scip, decompositions[i], subscipconss, nsubscipconss );

            DECdecompSetConstoblock(decompositions[i], constoblock);
            DECdecompSetConsindex(decompositions[i], consindex);

            /* finished setting constraint data structures */
            /** now: set variables */


            SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &nsubscipvars, seeed->getNBlocks() ) );
            SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &subscipvars, seeed->getNBlocks() ) );
            SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &stairlinkingvars, seeed->getNBlocks() ) );
            SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &nstairlinkingvars, seeed->getNBlocks() ) );

            SCIP_CALL_ABORT( SCIPhashmapCreate( &vartoblock, SCIPblkmem(scip), seeed->getNVars() ) );
            SCIP_CALL_ABORT( SCIPhashmapCreate( &varindex, SCIPblkmem(scip), seeed->getNVars() ) );

             /** set linkingvars */

            nlinkingvars = seeed->getNLinkingvars() + seeed->getNMastervars() + seeed->getNTotalStairlinkingvars();

            if( nlinkingvars != 0 )
               SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &linkingvars, nlinkingvars) );
            else
               linkingvars = NULL;

            for( int v = 0; v < seeed->getNLinkingvars(); ++v )
            {
               int var = seeed->getLinkingvars()[v];
               SCIP_VAR* scipvar = SCIPvarGetProbvar( varToScipVar[var] );
               assert(scipvar != NULL);

               linkingvars[v] = scipvar;
               SCIP_CALL_ABORT( SCIPhashmapInsert(vartoblock, scipvar, (void*) (size_t) (seeed->getNBlocks() + 2) ) );
               SCIP_CALL_ABORT( SCIPhashmapInsert(varindex, scipvar, (void*) (size_t) varcounter) );
               varcounter++;
            }

            for( int v = 0; v < seeed->getNMastervars(); ++v )
            {
               int var = seeed->getMastervars()[v];
               SCIP_VAR* scipvar = SCIPvarGetProbvar( varToScipVar[var] );
               linkingvars[v+seeed->getNLinkingvars()] = scipvar;
               SCIP_CALL_ABORT( SCIPhashmapInsert(vartoblock, scipvar, (void*) (size_t) (seeed->getNBlocks() + 1) ) );
               SCIP_CALL_ABORT( SCIPhashmapInsert(consindex, scipvar, (void*) (size_t) varcounter) );
               varcounter++;
            }


            /* set block variables */
            for ( int b = 0; b < seeed->getNBlocks(); ++b )
            {

               if(seeed->getNVarsForBlock(b) > 0)
                  SCIP_CALL_ABORT(SCIPallocBufferArray(scip, &subscipvars[b], seeed->getNVarsForBlock(b) ) );
               else subscipvars[b] = NULL;

               if(seeed->getNStairlinkingvars(b) > 0)
                  SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &stairlinkingvars[b], seeed->getNStairlinkingvars(b) ) );
               else stairlinkingvars[b] = NULL;

               nsubscipvars[b] = seeed->getNVarsForBlock(b);
               nstairlinkingvars[b] = seeed->getNStairlinkingvars(b);

               for ( int v = 0; v < seeed->getNVarsForBlock(b); ++v )
               {
                  int var = seeed->getVarsForBlock(b)[v];
                  SCIP_VAR* scipvar = SCIPvarGetProbvar( varToScipVar[var] );
                  assert(scipvar != NULL);

                  subscipvars[b][v] = scipvar;
                  SCIP_CALL_ABORT( SCIPhashmapInsert(vartoblock, scipvar, (void*) (size_t) (b + 1) ) );
                  SCIP_CALL_ABORT( SCIPhashmapInsert(varindex, scipvar, (void*) (size_t) varcounter) );
                  varcounter++;
               }

               for ( int v = 0; v < seeed->getNStairlinkingvars(b); ++v )
                {
                   int var = seeed->getStairlinkingvars(b)[v];
                   SCIP_VAR* scipvar = SCIPvarGetProbvar( varToScipVar[var] );
                   assert(scipvar != NULL);

                   stairlinkingvars[b][v] = scipvar;
                   linkingvars[seeed->getNLinkingvars() + seeed->getNMastervars() + counterstairlinkingvars ] = scipvar;
                   SCIP_CALL_ABORT( SCIPhashmapInsert(vartoblock, scipvar, (void*) (size_t) (seeed->getNBlocks() + 2) ) );
                   SCIP_CALL_ABORT( SCIPhashmapInsert(varindex, scipvar, (void*) (size_t) varcounter) );
                   varcounter++;
                   counterstairlinkingvars++;
                }
            }

            DECdecompSetSubscipvars(scip, decompositions[i], subscipvars, nsubscipvars);
            DECdecompSetStairlinkingvars(scip, decompositions[i], stairlinkingvars, nstairlinkingvars);
            DECdecompSetLinkingvars(scip, decompositions[i], linkingvars, nlinkingvars);
            DECdecompSetVarindex(decompositions[i], varindex);
            DECdecompSetVartoblock(decompositions[i], vartoblock) ;

            /** free stuff */

            /** free constraints */

            SCIPfreeBufferArrayNull(scip, &(linkingconss));
            SCIPfreeBufferArrayNull(scip, &(nsubscipconss));
            for( int b = seeed->getNBlocks()-1; b >= 0; --b )
            {
               SCIPfreeBufferArrayNull(scip, &(subscipconss[b]));
            }
            SCIPfreeBufferArrayNull(scip, &(subscipconss));

            /** free vars stuff */

            SCIPfreeBufferArrayNull(scip, &(linkingvars) );
            for( int b = seeed->getNBlocks()-1; b >= 0; --b )
            {
                  if( nsubscipvars[b] != 0 )
                  {
                     SCIPfreeBufferArrayNull(scip, &(subscipvars[b]));
                  }
            }

            SCIPfreeBufferArrayNull(scip, &(subscipvars) );
            SCIPfreeBufferArrayNull(scip, &(nsubscipvars));

            for( int b = seeed->getNBlocks()-1; b >= 0; --b )
             {
                if( nstairlinkingvars[b] != 0 )
                {
                   SCIPfreeBufferArrayNull(scip, &(stairlinkingvars[b]));
                }
             }
            SCIPfreeBufferArrayNull(scip, &(stairlinkingvars) );
            SCIPfreeBufferArrayNull(scip, &(nstairlinkingvars));





            /**** OLD STUFF below */
//
//            for( int j = 0; j < seeed->getNMastervars(); ++j )
//            {
//               decompositions[i]->linkingvars[j + seeed->getNLinkingvars()] = SCIPvarGetProbvar( getVarForIndex(seeed->getMastervars()[j]));
//               SCIPcaptureVar(scip, decompositions[i]->linkingvars[j + seeed->getNLinkingvars()]);
//            }
//
//            /** set stairlinkingvars */
//            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->stairlinkingvars, nblocks) ); /** free in decomp.c:466 */
//            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->nstairlinkingvars, nblocks) ); /** free in decomp.c:467 */
//
//            for ( int j = 0; j < nblocks; ++j )
//            {
//              decompositions[i]->nstairlinkingvars[j] = seeed->getNStairlinkingvars(j);
//              SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->stairlinkingvars[j], decompositions[i]->nstairlinkingvars[j]) ); /** free in decomp.c:444 */
//               for ( int k = 0; k < decompositions[i]->nstairlinkingvars[j]; ++k )
//               {
//                  decompositions[i]->stairlinkingvars[j][k] = SCIPvarGetProbvar( getVarForIndex(seeed->getStairlinkingvars(j)[k]) );
//                  SCIPcaptureVar(scip, decompositions[i]->stairlinkingvars[j][k]);
//               }
//            }
//
//            /** create hashmap constoblock and consindex */
//
//
//            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->subscipconss, nblocks) ); /** free in decomp.c:463 */
//            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->nsubscipconss, nblocks) ); /** free in decomp.c:464 */
//            /** add block conss */
//            for( int b = 0; b < seeed->getNBlocks(); ++b )
//            {
//               currblock = b+1;
//               decompositions[i]->nsubscipconss[b] = seeed->getNConssForBlock(b);
//               SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->subscipconss[b], decompositions[i]->nsubscipconss[b]) ); /** free in decomp.c:416 */
//               for ( int j = 0; j < seeed->getNConssForBlock(b); ++j)
//               {
//                  cindex++;
//                  cons = getConsForIndex((seeed->getConssForBlock(b))[j]);
//                  decompositions[i]->subscipconss[b][j] = cons;
//                  SCIPdebugMessage("cons %d is cons of block %d\n", cindex, currblock );
//                  SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->constoblock, cons, (void*) (size_t) currblock) );
//                  SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->consindex, (void*) (size_t) cindex, (void*) (size_t) currblock) );
//                  SCIP_CALL_ABORT( SCIPcaptureCons(scip, cons) );
//               }
//            }
//
//            /** add master conss */
//            for( int j = 0; j < seeed->getNMasterconss(); ++j )
//            {
//               cindex++;
//               cons = getConsForIndex(seeed->getMasterconss()[j]);
//               SCIPdebugMessage("cons %d is mastercons\n", cindex);
//               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->constoblock, cons, (void*) (size_t) (nblocks+1)) );
//               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->consindex, (void*) (size_t) cindex, (void*) (size_t) (nblocks+1)) );
//               SCIP_CALL_ABORT( SCIPcaptureCons(scip, cons) );
//            }
//
//            /** create hashmap vartoblock and varindex */
//            SCIP_CALL_ABORT( SCIPhashmapCreate( &decompositions[i]->vartoblock, SCIPblkmem(scip), nVars ) );
//            SCIP_CALL_ABORT( SCIPhashmapCreate( &decompositions[i]->varindex, SCIPblkmem(scip), nVars ) );
//
//            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->subscipvars, nblocks) ); /** free in decomp.c:461 */
//            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->nsubscipvars, nblocks) ); /** free in decomp.c:462 */
//
//            SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &stairlinkingvars, nblocks) );
//            SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &nstairlinkingvars, nblocks) );
//
//
//            /** add block vars and stairlinkingvars */
//            for( int b = 0; b < nblocks; ++b)
//            {
//               currblock = b+1;
//               decompositions[i]->nsubscipvars[b] = seeed->getNVarsForBlock(b);
//               if ( decompositions[i]->nsubscipvars[b] != 0)
//                  SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->subscipvars[b], decompositions[i]->nsubscipvars[b]) ); /** free in decomp.c:405 */
//               else decompositions[i]->subscipvars[b] = NULL;
//
//               SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &(stairlinkingvars[b]), seeed->getNStairlinkingvars(b)) ); /*lint !e866*/
//               nstairlinkingvars[b] = 0;
//
//
//               for ( int j = 0; j < seeed->getNVarsForBlock(b); ++j )
//               {
//                  vindex++;
//                  probvar = SCIPvarGetProbvar( getVarForIndex(seeed->getVarsForBlock(b)[j]) );
//                  decompositions[i]->subscipvars[b][j] = probvar;
//                  if( SCIPhashmapExists(decompositions[i]->vartoblock, probvar) )
//                  {
//                     SCIPdebugMessage("var <%s> has been handled before, it should not been add to block %d\n", SCIPvarGetName(probvar), currblock);
//                  }
//                  else
//                  {
//                     SCIPdebugMessage("var <%s> has not been handled before, adding to block %d\n", SCIPvarGetName(probvar), currblock );
//                     SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->vartoblock, probvar, (void*) (size_t) currblock) );
//                     SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->varindex, (void*) (size_t) vindex, (void*) (size_t) currblock) );
//                     SCIP_CALL_ABORT( SCIPcaptureVar(scip, probvar) );
//                  }
//               }
//
//               for( int k = 0; k < seeed->getNStairlinkingvars(b); ++k )
//               {
//                  probvar = SCIPvarGetProbvar( getVarForIndex((seeed->getStairlinkingvars(b))[k]));
//                  if( !SCIPhashmapExists(decompositions[i]->vartoblock, probvar) )
//                  {
//                     vindex++;
//                     SCIPdebugMessage("var <%s> is stairlinkingvar\n", SCIPvarGetName(probvar));
//                     SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->vartoblock, probvar, (void*) (size_t) (nblocks+2)) );
//                     SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->varindex, (void*) (size_t) vindex, (void*) (size_t) (nblocks+2)) );
//                     SCIP_CALL_ABORT( SCIPcaptureVar(scip, probvar) );
// //                    stairlinkingvars[b][k] = probvar;
//                  }
//                  assert( ((size_t) SCIPhashmapGetImage(decompositions[i]->vartoblock, probvar)) == (size_t)nblocks+2);
//               }
//            }
//
//
//            SCIPfreeBufferArrayNull(scip, &(stairlinkingvars));
//            SCIPfreeBufferArrayNull(scip, &(nstairlinkingvars));
//
//
//              /** add linking vars */
//            for( int j = 0; j < seeed->getNLinkingvars(); ++j )
//            {
//               vindex++;
//               probvar = SCIPvarGetProbvar( getVarForIndex(seeed->getLinkingvars()[j]) );
//               SCIPdebugMessage("var <%s> is linkingvar\n", SCIPvarGetName(probvar));
//               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->vartoblock, probvar, (void*) (size_t) (nblocks+2)) );
//               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->varindex, (void*) (size_t) vindex, (void*) (size_t) (nblocks+2)) );
//               SCIP_CALL_ABORT( SCIPcaptureVar(scip, probvar) );
//            }
//
//            /** add master vars */
//            for( int j = 0; j < seeed->getNMastervars(); ++j )
//            {
//               vindex++;
//               probvar = SCIPvarGetProbvar( getVarForIndex(seeed->getMastervars()[j]) );
//               SCIPdebugMessage("var <%s> is mastervar\n", SCIPvarGetName(probvar));
//               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->vartoblock, probvar, (void*) (size_t) (nblocks+1)) );
//               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->varindex, (void*) (size_t) vindex, (void*) (size_t) (nblocks+1)) );
//               SCIP_CALL_ABORT( SCIPcaptureVar(scip, probvar) );
//            }


//            /** test detector chain output */
//            char detectorchainstring[SCIP_MAXSTRLEN];
//
//            sprintf(detectorchainstring, "%s", DECdetectorGetName(decompositions[i]->detectorchain[0]));
//
//              for( i=1; i < ndetectors; ++i )
//              {
//                 sprintf(detectorchainstring, "%s-%s",detectorchainstring, DECdetectorGetName(decompositions[i]->detectorchain[i]) );
//              }
//
//              SCIPinfoMessage(scip, NULL, "%s %s", detectorchainstring, LINEBREAK);


            /*** OLD stuff above */


            /** set detectorchain */
            int ndetectors = seeed->getNDetectors();
            decompositions[i]->sizedetectorchain = ndetectors;
            size = SCIPcalcMemGrowSize(scip, decompositions[i]->sizedetectorchain);
            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->detectorchain, size) ); /** free in decomp.c:469 */
            for( int k = 0; k < ndetectors; ++k )
            {
               //          std::cout << " added detector of " << i << "-th seeed to its detetcor chain" << std::endl;
               decompositions[i]->detectorchain[k] = getDetectorForIndex(seeed->getDetectorchain()[k]);
            }


            /** set statistical detector chain data */

            DECdecompSetSeeedID(decompositions[i], seeed->getID() );
            if(seeed->getNDetectors() > 0 )
            {
               DECdecompSetDetectorClockTimes(scip, decompositions[i], &(seeed->detectorClockTimes[0]) );
               DECdecompSetDetectorPctVarsToBorder(scip, decompositions[i], &(seeed->pctVarsToBorder[0] ) );
               DECdecompSetDetectorPctVarsToBlock(scip, decompositions[i], &(seeed->pctVarsToBlock[0] ) );
               DECdecompSetDetectorPctVarsFromOpen(scip, decompositions[i], &(seeed->pctVarsFromFree[0] ) );
               DECdecompSetDetectorPctConssToBorder(scip, decompositions[i], &(seeed->pctConssToBorder[0] ) );
               DECdecompSetDetectorPctConssToBlock(scip, decompositions[i], &(seeed->pctConssToBlock[0] ) );
               DECdecompSetDetectorPctConssFromOpen(scip, decompositions[i], &(seeed->pctConssFromFree[0] ) );
               DECdecompSetNNewBlocks(scip, decompositions[i], &(seeed->nNewBlocks[0] ) );
            }
            /** set dectype */
            if(decompositions[i]->nlinkingvars == seeed->getNTotalStairlinkingvars() && decompositions[i]->nlinkingconss == 0 && DECdecompGetNLinkingvars(decompositions[i]) > 0)
            {
               decompositions[i]->type = DEC_DECTYPE_STAIRCASE;
            }
            else if(decompositions[i]->nlinkingvars > 0 || seeed->getNTotalStairlinkingvars() )
            {
               decompositions[i]->type = DEC_DECTYPE_ARROWHEAD;
            }
            else if(decompositions[i]->nlinkingconss > 0)
            {
               decompositions[i]->type = DEC_DECTYPE_BORDERED;
            }
            else if(decompositions[i]->nlinkingconss == 0 && seeed->getNTotalStairlinkingvars() == 0)
            {
               decompositions[i]->type = DEC_DECTYPE_DIAGONAL;
            }
            else
            {
               decompositions[i]->type = DEC_DECTYPE_UNKNOWN;
            }

            ndecompositions++;

            assert(DECdecompCheckConsistency(scip, decompositions[i] ) );

            assert(!SCIPhashmapIsEmpty(decompositions[i]->constoblock));
            assert(!SCIPhashmapIsEmpty(decompositions[i]->vartoblock));


         }

         //SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );

         /** delete the seeeds */
         for(size_t c = 0; c < currSeeeds.size(); ++c)
         {
            duplicate = false;
            for(size_t d = 0; d < delSeeeds.size(); ++d)
            {
               if(currSeeeds[c]==delSeeeds[d])
               {
                  duplicate=true;
                  break;
               }
            }
            if(!duplicate)
            {
               delSeeeds.push_back(currSeeeds[c]);
            }
         }

         for(size_t f = 0; f < finishedSeeeds.size(); ++f)
         {
            duplicate = false;
            for(size_t d = 0; d < delSeeeds.size(); ++d)
            {
               if(finishedSeeeds[f]==delSeeeds[d])
               {
                  duplicate=true;
                  break;
               }
            }
            if(!duplicate)
            {
               delSeeeds.push_back(finishedSeeeds[f]);
            }
         }


         for( size_t d =  delSeeeds.size(); d > 0; d--)
         {
            delete delSeeeds[d-1];
         }

         delSeeeds.clear();

         delete seeedPropData;

         return;

 }

/*SCIP_RETCODE DECdecompCheckConsistency(DEC_DECOMP* decomp)
{
   int c;
   int b;
   int v;

   for( v = 0; v < SCIPgetNVars(scip); ++v )
   {
      assert(SCIPhashmapExists(DECdecompGetVartoblock(decomp), SCIPgetVars(scip)[v]));

   }
}*/

const  int * Seeedpool::getVarsForCons(int cons){
         return &varsForConss[cons][0];
 }

const  SCIP_Real * Seeedpool::getValsForCons(int cons){
         return &valsForConss[cons][0];
 }


 /** access coefficient matrix variable-wise */
 const  int * Seeedpool::getConssForVar(int var){
         return &conssForVars[var][0];
 }

 int Seeedpool::getNVarsForCons(int cons){
    return varsForConss[cons].size();
 }

 int Seeedpool::getNConssForVar(int var){
    return conssForVars[var].size();
 }

 SCIP_VAR* Seeedpool::getVarForIndex(int varIndex){
         return varToScipVar[varIndex];
 }

 SCIP_CONS* Seeedpool::getConsForIndex(int consIndex){
         return consToScipCons[consIndex];
 }

 DEC_DETECTOR* Seeedpool::getDetectorForIndex(int detectorIndex){
    return detectorToScipDetector[detectorIndex];
 }

 SCIP_Real Seeedpool::getVal(int row, int col){

    std::tr1::unordered_map< std::pair<int, int>, SCIP_Real, pair_hash>::const_iterator iter =  valsMap.find(std::pair<int, int>(row, col) ) ;

    if ( iter == valsMap.end()  )
       return 0;

    return iter->second;
 }

 int Seeedpool::getIndexForVar(SCIP_VAR* var){
         return scipVarToIndex[var];
 }

 int Seeedpool::getIndexForCons(SCIP_CONS* cons){
         return scipConsToIndex[cons];
 }

 int Seeedpool::getIndexForDetector(DEC_DETECTOR* detector){
    return scipDetectorToIndex[detector];
 }

 int Seeedpool::getNewIdForSeeed(){
    nTotalSeeeds++;
    return (nTotalSeeeds-1);
 }

 void Seeedpool::decrementSeeedcount(){
     nTotalSeeeds--;
     return;
  }


 DEC_DECOMP** Seeedpool::getDecompositions(){
    return decompositions;
 }

 int Seeedpool::getNDecompositions(){
    return ndecompositions;
 }

 int Seeedpool::getNDetectors(){
    return nDetectors;
 }

 int Seeedpool::getNVars(){
    return nVars;
 }

 int Seeedpool::getNConss(){
    return nConss;
 }

 std::vector<int> Seeedpool::getCandidatesNBlocks() const
 {
    return candidatesNBlocks;
 }

 void Seeedpool::addCandidatesNBlocks(
    int                 candidate            /**< candidate for block size */
    )
 {

    if(candidate > 1)
    {
       bool alreadyIn = false;
       for(size_t i = 0; i < candidatesNBlocks.size(); ++i )
       {
          if(candidatesNBlocks[i] == candidate)
          {
             alreadyIn = true;
             break;
          }
       }
       if(!alreadyIn)
       {
          std::cout << "added block number candidate : " << candidate << std::endl;
          candidatesNBlocks.push_back(candidate);
       }
    }
 }

 void Seeedpool::calcCandidatesNBlocks()
  {
    /**
     * for every subset of constraint classes calculate gcd (greatest common divisors) of the corresponding number of occurrences
     */


    for( size_t conssclass = 0; conssclass < consclassescollection.size(); ++conssclass )
    {
       std::vector<int> nconssofclass = std::vector<int>(consclassesnclasses[conssclass], 0);
       std::vector<int> consclassindices = std::vector<int>(0);
       for( int i = 0; i < consclassesnclasses[conssclass]; ++ i)
          consclassindices.push_back(i);
       std::vector< std::vector<int> > subsetsOfConstypes = getAllSubsets(consclassindices);

       for ( int i =0; i < getNConss(); ++i)
          ++nconssofclass[ consclassescollection[conssclass][i] ];

       for(size_t subset = 0; subset < subsetsOfConstypes.size(); ++subset)
       {
          int greatestCD = 1;

          if(subsetsOfConstypes[subset].size() == 0 || subsetsOfConstypes[subset].size() == 1)
               continue;

          greatestCD = gcd(nconssofclass[subsetsOfConstypes[subset][0]], nconssofclass[subsetsOfConstypes[subset][1]]  );

          for( size_t i = 2; i < subsetsOfConstypes[subset].size() ; ++i)
          {
             greatestCD = gcd( greatestCD, nconssofclass[subsetsOfConstypes[subset][i]] );
          }

          if( greatestCD > 1 )
          {
             bool alreadyIn = false;
             for( size_t i = 0; i < candidatesNBlocks.size(); ++i )
             {
                if( candidatesNBlocks[i] == greatestCD )
                {
                   alreadyIn = true;
                   break;
                }
             }
             if( !alreadyIn )
             {
                std::cout << "added block number candidate : " << greatestCD << std::endl;
                candidatesNBlocks.push_back(greatestCD);
             }
          }
       }
    }

    return ;
  }

 int Seeedpool::getNConssClassDistributions(){
    return consclassescollection.size();
 }

 int* Seeedpool::getConssClassDistribution(
    int consclassdistr
    ){
    return &(consclassescollection[consclassdistr][0]);
 }

 int Seeedpool::getNClassesOfDistribution(
    int consclassdistr
    )
 {
    return consclassesnclasses[consclassdistr];
 }

 void Seeedpool::addConssClassesForSCIPConstypes()
 {
    /**
     * at first for every subset of constypes calculate gcd (greatest common divisors) of the corresponding number of occurrences
     */
    std::vector<consType> foundConstypes(0);
    std::vector<int> constypesIndices(0);
    std::vector<int> nConssConstype(0);
    std::vector<int> classForCons = std::vector<int>(getNConss(), -1);

    for( int i = 0; i < getNConss(); ++i)
    {
       SCIP_CONS* cons;
       bool found = false;
       cons = getConsForIndex(i);
       consType cT = GCGconsGetType(cons);
       size_t constype;

       /** find constype or not */

       for( constype = 0; constype < foundConstypes.size(); ++constype)
       {
          if( foundConstypes[constype] == cT )
          {
             found = true;
             break;
          }
       }
       if( !found )
       {
          foundConstypes.push_back(GCGconsGetType(cons) );
          classForCons[i] = foundConstypes.size() - 1;
       }
       else
          classForCons[i] = constype;
     }

    consclassescollection.push_back(classForCons);
    consclassesnclasses.push_back(foundConstypes.size() );

    return;
 }

 void Seeedpool::addConssClassesForConsnamesDigitFreeIdentical()
  {
     /**
      * at first remove all digits from the consnames
      */
     std::vector<std::string> consnamesToCompare(getNConss(), "");
     std::vector<int> nConssConstype(0);
     std::vector<int> classForCons = std::vector<int>(getNConss(), -1);
     std::vector<std::string> nameClasses(0);


     for( int i = 0; i < getNConss(); ++i )
     {
        int nremoved;
        char consname[SCIP_MAXSTRLEN];
        strcpy(consname, SCIPconsGetName(getConsForIndex(i) ) );

        removeDigits(consname, &nremoved);
        consnamesToCompare[i] = std::string(consname);
     }

     /** test of reduced consnames */

     if( false )
     {
        for( int i = 0; i < getNConss(); ++i )
        {
          std::cout << " old consname : " << SCIPconsGetName(getConsForIndex(i) ) << std::endl;
//           std::cout << " new consname : " << consnamesToCompare[i]  << std::endl;
//           std::cout << std::endl;
        }
     }

     for( int i = 0; i < getNConss(); ++i )
     {
        bool belongstoexistingclass = false;
        /** test if string belongs to an existing name class */
        for ( size_t j = 0; j < nameClasses.size(); ++j )
        {
           if ( nameClasses[j].compare(consnamesToCompare[i]) == 0 )
           {
              belongstoexistingclass = true;
              classForCons[i] = j;
              nConssConstype[j]++;
              break;
           }
        }
        if ( !belongstoexistingclass )
        {
           nameClasses.push_back(consnamesToCompare[i] );
           nConssConstype.push_back(1);
           classForCons[i] = nameClasses.size()-1;

        }
     }

     consclassescollection.push_back(classForCons);
     consclassesnclasses.push_back(nameClasses.size() );

     std::cout << " comparison of digit-reduced consnames for identity  " << " yields a distribution with " << nameClasses.size()  << " different constraint classes" << std::endl;

     return;


  }

 void Seeedpool::addConssClassesForConsnamesLevenshteinDistanceConnectivity(
    int connectivity
    )
  {
     /**
      * at first remove all digits from the consnames
      */
     std::vector<std::string> consnamesToCompare(getNConss(), "");
     std::vector<int> nConssConstype(0);
     std::vector<int> classForCons = std::vector<int>(getNConss(), -1);
     std::vector<bool> alreadyReached(getNConss(), false);
     std::queue<int> helpqueue = std::queue<int>();
     //std::vector<int> neighborConss(0);
     int nUnreachedConss = getNConss();
     int currentClass = -1;

     std::vector< std::vector<int> > levenshteindistances(getNConss(), std::vector<int>(getNConss(), -1) );

     for( int i = 0; i < getNConss(); ++i )
     {
        consnamesToCompare[i] = std::string(SCIPconsGetName(getConsForIndex(i) ));
     }

     for( int i = 0; i < getNConss(); ++i )
     {
        for ( int j = i+1; j < getNConss(); ++j)
        {
           levenshteindistances[i][j] = calcLevenshteinDistance(consnamesToCompare[i], consnamesToCompare[j]);
           levenshteindistances[j][i] = levenshteindistances[i][j];
//           if(levenshteindistances[i][j] == 1)
//              std::cout << " string1 : " << consnamesToCompare[i] << " string2 : " << consnamesToCompare[j] << " distance: " << levenshteindistances[i][j] << std::endl;
        }
     }

     /** do breadth first search */
        while( nUnreachedConss > 0 )
        {
           int firstUnreached = -1;
           currentClass++;
           assert(helpqueue.empty());
           for( int i = 0; i < getNConss(); ++i )
           {
              if( classForCons[i] == -1 )
              {
                 firstUnreached = i;
                 break;
              }
           }

           helpqueue.push(firstUnreached);
//           neighborConss.clear();
//           neighborConss.push_back(firstUnreached);
           alreadyReached[firstUnreached] = true;
           classForCons[firstUnreached] = currentClass;
           --nUnreachedConss;

           while( !helpqueue.empty() )
           {
              int nodecons = helpqueue.front();
              helpqueue.pop();
              for( int j = 0; j < getNConss() ; ++j )
              {

                 if( alreadyReached[j] )
                    continue;

                 if(j == nodecons)
                    continue;

                 if(levenshteindistances[j][nodecons] > connectivity)
                    continue;

                    alreadyReached[j] = true;
                    classForCons[j] = currentClass;
                    --nUnreachedConss;
                    helpqueue.push(j);
                 }
           } //endwhile(!queue.empty() )

        } // endwhile( !openConss.empty() )

        consclassescollection.push_back(classForCons);
        consclassesnclasses.push_back(currentClass+1 );

        std::cout << " levenshtein connectivity of " << connectivity << " yields a distribution with " << currentClass+1  << " different constraint classes" << std::endl;


     return;


  }



} /* namespace gcg */
