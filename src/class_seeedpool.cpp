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
#include "scip/cons.h"


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

/** method to calculate the greatest common divisor */

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
        const char*             conshdlrName,
        SCIP_Bool                transformed
    ):scip(givenScip), currSeeeds(0), nTotalSeeeds(0),nVars(SCIPgetNVars(givenScip) ), nConss(SCIPgetNConss(givenScip) ), nDetectors(0), nFinishingDetectors(0),ndecompositions(0), candidatesNBlocks(0), transformed(transformed)
 {
         SCIP_CONS** conss;
         SCIP_VAR** vars;

         SCIP_CONSHDLR* conshdlr;  /** cons_decomp to get detectors */
         SCIP_CONSHDLRDATA* conshdlrdata;


         if( !transformed )
         {
            nVars = SCIPgetNOrigVars(scip);
            nConss = SCIPgetNOrigConss(scip);
         }




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

         /** set up enabled detectors */
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

         /** set up enabled finishing detectors */
         for(int d = 0; d < conshdlrdata->ndetectors; ++d )
         {
                 DEC_DETECTOR* detector;

                 detector = conshdlrdata->detectors[d];
                 assert(detector != NULL);
                 if( !detector->enabledFinishing || detector->finishSeeed == NULL)
                         continue;

                 scipFinishingDetectorToIndex[detector] = nFinishingDetectors;
                 detectorToFinishingScipDetector.push_back(detector);
                 ++nFinishingDetectors;
         }



         /** initilize matrix datastructures */
         if( transformed )
         {
            conss = SCIPgetConss(scip);
            vars = SCIPgetVars(scip);
         }
         else
         {
            conss = SCIPgetOrigConss(scip);
            vars = SCIPgetOrigVars(scip);
         }

         /** assign an index to every cons and var
          * @TODO: are all constraints/variables relevant? (probvars etc)  */

         for(int i = 0; i < nConss; ++i)
         {
                 SCIP_CONS* relevantCons;

                 relevantCons = transformed ? consGetRelevantRepr(scip, conss[i]) : conss[i];

                 if( relevantCons != NULL )
                 {
                         scipConsToIndex[relevantCons] = relevantConsCounter ;
                         consToScipCons.push_back(relevantCons);
                         ++relevantConsCounter;
                 }
                 else
                 {
                    std::cout << "NULL" << std::endl;
                 }
         }

         for(int i = 0; i < nVars; ++i)
         {
                 SCIP_VAR* relevantVar;


                 if( transformed )
                    relevantVar = varGetRelevantRepr(scip, vars[i]);
                 else relevantVar = vars[i];

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
         std::cout << "nVars: " << nVars << " / nConss: " << nConss << std::endl;
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

         /*    seeedpool with empty seeed and translated original seeeds*/

         currSeeeds.push_back(new Seeed(scip, nTotalSeeeds,nDetectors,nConss,nVars) );

         for (size_t i = 0; i < translatedOrigSeeeds.size(); ++i)
         {
            translatedOrigSeeeds[i]->calcHashvalue();
            if( seeedIsNoDuplicateOfSeeeds(translatedOrigSeeeds[i], currSeeeds, true) )
                currSeeeds.push_back(translatedOrigSeeeds[i] );
         }

         nTotalSeeeds++;

         decompositions = NULL;

         addConssClassesForSCIPConstypes();
         addConssClassesForConsnamesDigitFreeIdentical();
     //    addConssClassesForConsnamesLevenshteinDistanceConnectivity(1);
         addConssClassesForNNonzeros();

         calcCandidatesNBlocks();


 }//end constructor

 Seeedpool::~Seeedpool(){

 }


 /** finds decompositions  */
  /** access coefficient matrlix constraint-wise */
 std::vector<SeeedPtr>    Seeedpool::findSeeeds(
 ){

         /** 1) read parameter, as there are: maxrounds
          *  2) loop rounds
          *  3) every seeed in seeeds
          *  4) every detector not registered yet propagates seeed
          *  5)  */


         SEEED_PROPAGATION_DATA* seeedPropData;
         bool displaySeeeds = false;
         int verboseLevel;
         std::vector<int> successDetectors;
         std::vector<SeeedPtr> delSeeeds;
         bool duplicate;

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

         for( int round = 0; round < maxndetectionrounds; ++round )
         {
                 std::cout << "currently in detection round " << round << std::endl;
                 std::vector<SeeedPtr> nextSeeeds = std::vector<SeeedPtr>(0);
                 std::vector<SeeedPtr> currSeeedsToDelete = std::vector<SeeedPtr>(0);

                 for( size_t s = 0; s < currSeeeds.size(); ++s )
                 {
                         SeeedPtr seeedPtr;
                         seeedPtr= currSeeeds[s];
                         if( displaySeeeds || verboseLevel >= 1 )
                         {
                            std::cout << "Start to propagate seeed " << seeedPtr->getID() << " in round " << round << ":" << std::endl;
                            if (displaySeeeds)
                               seeedPtr->displaySeeed();
                         }

                         /** the current seeed is handled by all detectors */
                         for( int d = 0; d < nDetectors; ++d )
                         {
                                 DEC_DETECTOR* detector;
                                 std::vector<SeeedPtr>::const_iterator newSIter;
                                 std::vector<SeeedPtr>::const_iterator newSIterEnd;

                                 SCIP_RESULT result = SCIP_DIDNOTFIND;
                                 detector = detectorToScipDetector[d];

                                 /** if the seeed is also propagated by the detector go on with the next detector */
                                 if(seeedPtr->isPropagatedBy(detector) && !detector->usefulRecall )
                                         continue;

                                 /** check if detector is callable in current detection round */
                                 if(detector->maxCallRound < round || detector->minCallRound > round)
                                     continue;

                                 if( (round - detector->minCallRound) % detector->freqCallRound != 0 )
                                     continue;

                                 seeedPropData->seeedToPropagate = seeedPtr;

                                 /** new seeeds are created by the current detector */
                                 SCIP_CALL_ABORT( SCIPstartClock(scip, detectorToScipDetector[d]->dectime) );
                                 if(verboseLevel >= 1)
                                     std::cout << "detector " << DECdetectorGetName(detectorToScipDetector[d]) << " started to propagate the " << s+1 << ". seeed (ID " << seeedPtr->getID() << ") in round " << round << std::endl;

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
                         } // end for detectors

                         for(int d = 0; d < nFinishingDetectors; ++d)
                         {
                            DEC_DETECTOR* detector = detectorToFinishingScipDetector[d];
                            SCIP_RESULT result = SCIP_DIDNOTFIND;
                            seeedPropData->seeedToPropagate = seeedPtr;

                            if(verboseLevel > 2 )
                               std::cout << "check if finisher of detector " << DECdetectorGetName(detectorToFinishingScipDetector[d] ) << " is enabled " << std::endl;

                            /** if the finishing of the detector is not enabled go on with the next detector */
                            if( !detector->enabledFinishing )
                                    continue;

                            if(verboseLevel > 2 )
                               std::cout << "call finisher for detector " << DECdetectorGetName(detectorToFinishingScipDetector[d] ) << std::endl;

                            SCIP_CALL_ABORT(detectorToFinishingScipDetector[d]->finishSeeed(scip, detectorToFinishingScipDetector[d], seeedPropData, &result) );

                            for(int finished = 0; finished < seeedPropData->nNewSeeeds; ++finished)
                            {
                               SeeedPtr seeed = seeedPropData->newSeeeds[finished];
                               seeedPropData->newSeeeds[finished]->sort();
                               seeed->calcHashvalue();
                               seeedPropData->newSeeeds[finished]->addDecChangesFromAncestor(seeedPtr);
                               seeed->setFinishedByFinisher(true);

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
                 }// end for currseeeds

                 for(size_t s = 0; s < currSeeedsToDelete.size(); ++s )
                    delete currSeeedsToDelete[s];

                 currSeeeds = nextSeeeds;
         } // end for rounds

         /** complete the currseeeds with finishing detectors and add them to finished seeeds */
         for(size_t i = 0; i < currSeeeds.size(); ++i)
         {
            SeeedPtr seeedPtr = currSeeeds[i];
            for(int d = 0; d < nFinishingDetectors; ++d)
            {
               DEC_DETECTOR* detector = detectorToFinishingScipDetector[d];
               SCIP_RESULT result = SCIP_DIDNOTFIND;
               seeedPropData->seeedToPropagate = seeedPtr;

               std::cout << "check if finisher of detector " << DECdetectorGetName(detectorToScipDetector[d] ) << " is enabled " << std::endl;

               /** if the finishing of the detector is not enabled go on with the next detector */
               if( !detector->enabledFinishing )
                  continue;

               std::cout << "call finisher for detector " << DECdetectorGetName(detectorToFinishingScipDetector[d] ) << std::endl;

               SCIP_CALL_ABORT(detectorToFinishingScipDetector[d]->finishSeeed(scip, detectorToFinishingScipDetector[d],seeedPropData, &result) );

               for(int finished = 0; finished < seeedPropData->nNewSeeeds; ++finished)
               {
                  SeeedPtr seeed = seeedPropData->newSeeeds[finished];
                  seeed->calcHashvalue();
                  seeed->addDecChangesFromAncestor(seeedPtr);
                  seeed->setFinishedByFinisher(true);

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
                 if(finishedSeeeds[i]->isPropagatedBy(detectorToScipDetector[d]))
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
//            bestSeeed->showScatterPlot(this);
         }


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


         for( size_t d =  delSeeeds.size(); d > 0; d--)
         {
            delete delSeeeds[d-1];
         }

         delSeeeds.clear();

         delete seeedPropData;

         return finishedSeeeds;

}

 void    Seeedpool::findDecompositions(
 ){

    /**
     * finds seeeds and translates them to decompositions
     *   */

    SEEED_PROPAGATION_DATA* seeedPropData;
    bool displaySeeeds = false;
    int verboseLevel;
    std::vector<int> successDetectors;
    std::vector<SeeedPtr> delSeeeds;
    std::vector<SeeedPtr> finishedSeeeds;
    bool duplicate;

    successDetectors = std::vector<int>(nDetectors, 0);
    ndecompositions = 0;
    seeedPropData = new SEEED_PROPAGATION_DATA();
    seeedPropData->seeedpool = this;
    seeedPropData->nNewSeeeds = 0;
    delSeeeds = std::vector<SeeedPtr>(0);

    verboseLevel = 0;

    finishedSeeeds = findSeeeds();

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
               if(k != ndetectors-1 || !seeed->getFinishedByFinisher() )
               {
               //          std::cout << " added detector of " << i << "-th seeed to its detetcor chain" << std::endl;
                  decompositions[i]->detectorchain[k] = seeed->getDetectorchain()[k];
               }else
                  decompositions[i]->detectorchain[k] = seeed->getDetectorchain()[k];
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


std::vector<Seeed*> Seeedpool::translateSeeeds( Seeedpool* origpool, std::vector<Seeed*> origseeeds )
{
   int nrowsother  = origpool->nConss;
   int nrowsthis  = nConss;
   int ncolsother  = origpool->nVars;
   int ncolsthis  = nVars;

   std::vector<Seeed*> newseeeds(0);
   std::vector<int> rowothertothis(nrowsother, -1);
   std::vector<int> rowthistoother(nrowsthis, -1);
   std::vector<int> colothertothis(ncolsother, -1);
   std::vector<int> colthistoother(ncolsthis, -1);

   std::vector<int> missingrowinthis(0);
   std::vector<int> newrowsthis(0);
   std::vector<int> missingcolsinthis(0);
   std::vector<int> newcolsthis(0);

   SCIPdebugMessagePrint(this->scip, " started translate seeed method \n" );


   /* identify new and deleted rows and cols; and identify bijection between maintained variables */

   for( int i = 0; i < nrowsother; ++i  )
   {
      SCIP_CONS* otherrow = origpool->getConsForIndex(i);
      assert(otherrow != NULL);
      if (SCIPconsGetTransformed(otherrow ) == NULL )
         continue;
      SCIP_Bool foundmaintained = FALSE;

//      SCIPdebugMessagePrint(this->scip, " otherrow: ptr %s ; name: %s  \n", SCIPconsGetTransformed(otherrow),  SCIPconsGetName(SCIPconsGetTransformed(otherrow))  );
      for( int j = 0; j < nrowsthis; ++j  )
      {
//         SCIPdebugMessagePrint(this->scip, " %s vs %s \n", SCIPconsGetName(otherrow), SCIPconsGetName(this->getConsForIndex(j))  );
         if( SCIPconsGetTransformed(otherrow) == this->getConsForIndex(j) )
         {
//            std::cout << " EQUAL! " << std::endl;
//            SCIPdebugMessagePrint(this->scip, " EQUAL \n" );
            rowothertothis[i] = j;
            rowthistoother[j] = i;
            foundmaintained = TRUE;
            break;
         }
      }
      if (!foundmaintained)
         missingrowinthis.push_back(i);
   }

   for( int i = 0; i < ncolsother; ++i  )
   {
      SCIP_VAR* othervar = origpool->getVarForIndex(i);
      SCIP_Bool foundmaintained = FALSE;
      for( int j = 0; j < ncolsthis; ++j  )
      {
         if( othervar == this->getVarForIndex(j) )
         {
            colothertothis[i] = j;
            colthistoother[j] = i;
            foundmaintained = TRUE;
            break;
         }
      }
   }

   SCIPdebugMessagePrint(this->scip, " calculated translation; number of missing constraints: %d; number of other seeeds: %d \n", missingrowinthis.size(), origseeeds.size() );

   /** constructing seeeds for this seeedpool */

   for( size_t s = 0; s < origseeeds.size(); ++s )
   {
      SeeedPtr otherseeed;
      SeeedPtr newseeed;

      otherseeed = origseeeds[s];

      SCIPdebugMessagePrint(this->scip, " otherseeed seeed %d has %d many blocks \n", otherseeed->getID(), otherseeed->getNBlocks() );

      /** ignore seeeds with one block or no block, they are supposed to be find anyway */
      if( otherseeed->getNBlocks() == 1 || otherseeed->getNBlocks() == 0  )
         continue;

      SCIPdebugMessagePrint(this->scip, " transform seeed %d \n", otherseeed->getID() );

      newseeed = new Seeed(scip, this->getNewIdForSeeed(), this->getNDetectors(), this->getNConss(), this->getNVars() );

      /** prepare new seeed */
      newseeed->calcOpenconss();
      newseeed->calcOpenvars();
      newseeed->setOpenVarsAndConssCalculated(true);

      newseeed->setNBlocks(otherseeed->getNBlocks() );


      /** set all (which have representative in the unpresolved seeed) constraints according to their representatives in the unpresolved seeed */
      for(int b = 0; b < otherseeed->getNBlocks() ; ++b )
      {
         for ( int i = 0; i < otherseeed->getNConssForBlock(b); i++ )
         {
            int thiscons = rowothertothis[otherseeed->getConssForBlock(b)[i] ];
            if( thiscons != -1 )
            {
               newseeed->setConsToBlock(thiscons, b);
               newseeed->deleteOpencons(thiscons);
            }
         }

/*         for ( int j = 0; j < otherseeed->getNVarsForBlock(b); j++ )
         {
            int thisvar = colothertothis[otherseeed->getVarsForBlock(b)[j] ];
            if( thisvar != -1 )
            {
               newseeed->setVarToBlock(thisvar, b);
               newseeed->deleteOpenvar(thisvar);
            }
         }*/
      }

      for ( int i = 0; i < otherseeed->getNMasterconss(); i++ )
      {
         int thiscons = rowothertothis[otherseeed->getMasterconss()[i] ];
         if( thiscons != -1 )
         {
            newseeed->setConsToMaster(thiscons);
            newseeed->deleteOpencons(thiscons);
         }
      }

      /** set linking and master vars according to their representatives in the unpresolved seeed */

      for ( int j = 0; j < otherseeed->getNLinkingvars(); j++ )
      {
         int thisvar = colothertothis[otherseeed->getLinkingvars()[j] ];
         if( thisvar != -1 )
         {
            newseeed->setVarToLinking(thisvar);
            newseeed->deleteOpenvar(thisvar);
         }
      }

      for ( int j = 0; j < otherseeed->getNMastervars(); j++ )
      {
         int thisvar = colothertothis[otherseeed->getMastervars()[j] ];
         if( thisvar != -1 )
         {
            newseeed->setVarToMaster(thisvar);
            newseeed->deleteOpenvar(thisvar);
         }
      }

      newseeed->detectorChain = otherseeed->detectorChain;

      newseeed->stemsFromUnpresolved = true;
      newseeed->isFinishedByFinisherUnpresolved  = otherseeed->isFinishedByFinisher;

      if ( otherseeed->isFinishedByFinisher )
         newseeed->finishedUnpresolvedBy = otherseeed->detectorChain[otherseeed->detectorChain.size()-1];


      newseeed->setFinishedByFinisher(otherseeed->isFinishedByFinisher);
      newseeed->sort();
      newseeed->considerImplicits(this);
      newseeed->deleteEmptyBlocks();
      newseeed->checkConsistency();

      if (newseeed->getNOpenconss() == getNConss() && newseeed->getNOpenvars() == newseeed->getNVars())
      {
         delete newseeed;
         continue;
      }

      std::cout << "unpresolved seeed " << std::endl;
      otherseeed->showScatterPlot(origpool);
      std::cout << "has become " << std::endl;
      newseeed->showScatterPlot(this);

      if(newseeed->checkConsistency() )
         newseeeds.push_back(newseeed);
      else delete newseeed;
   }

   return newseeeds;
}

void Seeedpool::populate(std::vector<SeeedPtr> seeeds){
   translatedOrigSeeeds = seeeds;
}



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

 DEC_DETECTOR* Seeedpool::getFinishingDetectorForIndex(int detectorIndex){
    return detectorToFinishingScipDetector[detectorIndex];
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

 int Seeedpool::getIndexForFinishingDetector(DEC_DETECTOR* detector){
     return scipFinishingDetectorToIndex[detector];
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

 int Seeedpool::getNFinishingDetectors(){
     return nFinishingDetectors;
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

    int maximumnclasses = 18; /* if  distribution of classes exceed this number its skipped */

    for( size_t conssclass = 0; conssclass < consclassescollection.size(); ++conssclass )
    {
       std::vector<int> nconssofclass = std::vector<int>(consclassesnclasses[conssclass], 0);
       std::vector<int> consclassindices = std::vector<int>(0);

       /** check if there are to  many classes in this distribution and skip it if so */

       if ( consclassesnclasses[conssclass] > maximumnclasses)
       {
          std::cout << " the current consclass distribution includes " <<  consclassesnclasses[conssclass] << " classes but only " << maximumnclasses << " are allowed for calcCandidatesNBlocks()" << std::endl;
          continue;
       }


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

 std::vector<int> Seeedpool::getConssClassDistributionVector(
     int consclassdistr
     ){
     return (consclassescollection[consclassdistr]);
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
     int nmaxconss = 5000;

     if (getNConss() > nmaxconss)
     {
        std::cout << " skipped levenshtein distance based constraint classes calculating since number of constraints " << getNConss() << " exceeds limit " << nmaxconss   << std::endl;
        return;
     }


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


 void Seeedpool::addConssClassesForNNonzeros()
  {
     /**
      * at first remove all digits from the consnames
      */
     std::vector<int> nconssforclass(0);
     std::vector<int> differentNNonzeros(0);
     std::vector<int> classForCons (getNConss(), -1);

     int counterClasses = 0;

     for( int i = 0; i < getNConss(); ++i )
     {
        int nnonzeros = getNVarsForCons(i);
        bool nzalreadyfound = false;

        for ( size_t nzid = 0; nzid < differentNNonzeros.size(); ++nzid )
        {
           if ( nnonzeros == differentNNonzeros[nzid] )
           {
              nzalreadyfound = true;
              classForCons[i] = nzid;
              ++nconssforclass[nzid];
              break;
           }
        }

        if(!nzalreadyfound)
        {
           classForCons[i] = counterClasses;
           ++counterClasses;
           differentNNonzeros.push_back(nnonzeros);
           nconssforclass.push_back(1);
        }
     }

     /** test of reduced consnames */

     if( true )
     {
        std::cout << " nNonzero : nConsWithNNonzero"  << std::endl;
        for( size_t i = 0; i < differentNNonzeros.size(); ++i )
        {
           std::cout << differentNNonzeros[i] << " : " << nconssforclass[i] << std::endl;
//           std::cout << " new consname : " << consnamesToCompare[i]  << std::endl;
//           std::cout << std::endl;
        }
     }

     consclassescollection.push_back(classForCons);
     consclassesnclasses.push_back(differentNNonzeros.size() );

     std::cout << " comparison of number of nonzeros  " << " yields a distribution with " << differentNNonzeros.size()  << " different constraint classes" << std::endl;

     return;


  }

 void Seeedpool::addConssClassDistribution(
     std::vector<int>                             conssClassDistribution,     /**< distribution to add */
     std::vector<SCIP_CONS*>                      indexToCons                 /**< stores the corresponding scip constraints pointer */
     )
  {
     int consindex;
     SCIP_CONS* cons;
     int classCounter = 0;
     std::vector<int> classes;
     int consCounter = 0;

     //fillout the distributionvector
     std::vector<int> distribution (nConss);

     for( size_t c = 0; c < indexToCons.size(); ++c )
     {
        cons = indexToCons[c];
        if( find(consToScipCons.begin(), consToScipCons.end(), cons) == consToScipCons.end() )
           continue;
        consindex = getIndexForCons(cons);
        distribution[consindex] = conssClassDistribution[c];
        consCounter ++;
        if( find( classes.begin(), classes.end(), conssClassDistribution[c]) == classes.end() )
        {
           classCounter++;
           classes.push_back(conssClassDistribution[c]);
        }

     }

     if( consCounter != nConss )
     {
        std::cout << "Distribution can not be added because of missing constraints!" << std::endl;
        return;
     }

     //the different classes should be consecutively numbered
     int maximum = *std::max_element(distribution.begin(), distribution.end());
     assert( maximum >= classCounter - 1);
     while( maximum != classCounter - 1 )
     {
        bool found;
        int number = -1;

        //search for a number between 0 and (classCounter - 1) which is not assigned to a class yet
        do
        {
           found = false;
           number++;
           for( size_t j = 0; j < distribution.size(); ++j )
           {
              if( number == distribution[j] )
              {
                 found = true;
                 break;
              }
           }
        }while( found && number < classCounter );

        assert( number != classCounter );

        for( size_t j = 0; j < distribution.size(); ++j )
        {
           if( distribution[j] == maximum )
              distribution[j] = number;
        }
     }

     if(distributionIsNoDuplicateOfDistributions(distribution, classCounter, consclassescollection))
     {
        consclassescollection.push_back( distribution );
        consclassesnclasses.push_back( classCounter );
     }



 }

 bool Seeedpool::distributionIsNoDuplicateOfDistributions(
    std::vector<int>              compDistribution,
    int                           nClasses,
    std::vector<std::vector<int>> distributions
    )
 {
    assert( (int) compDistribution.size() == nConss );
    bool equal = false;
    for( size_t j = 0; j < distributions.size(); ++j )
    {
       equal = true;
       for( int k = 0; k < nClasses && equal; ++k )
       {
          int classnumber = -1;
          for( size_t l = 0; l < compDistribution.size() && equal; ++l )
          {
             if( compDistribution[l] == k )
             {
                if( classnumber == -1 )
                   classnumber = consclassescollection[j][l];
                if( consclassescollection[j][l] != classnumber )
                   equal = false;
             }
          }
       }
       if(equal)
          return false;
    }
    return true;
 }




} /* namespace gcg */
