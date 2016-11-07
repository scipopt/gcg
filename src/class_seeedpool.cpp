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

#include <algorithm>
#include <iostream>


#include <exception>

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

SCIP_CONS* consGetRelevantRepr(SCIP* scip, SCIP_CONS* cons){

                return cons;
}

SCIP_VAR* varGetRelevantRepr(SCIP* scip, SCIP_VAR* var){

        return SCIPvarGetProbvar(var);
}

SCIP_Bool seeedIsNoDuplicateOfSeeeds(SeeedPtr compseeed, std::vector<SeeedPtr> const & seeeds, bool sort){

   for( size_t i = 0; i < seeeds.size(); ++i )
   {
      bool noDuplicate = false;

      /** compares the number of master conss, master vars, blocks, linking vars and stairlinking vars */
      if( compseeed->getNMasterconss() != seeeds[i]->getNMasterconss() || compseeed->getNMastervars() != seeeds[i]->getNMastervars() ||
         compseeed->getNBlocks() != seeeds[i]->getNBlocks() || compseeed->getNLinkingvars() != seeeds[i]->getNLinkingvars())
         continue;

      /** compares the number of stairlinking vars */
      for( int b = 0; b < compseeed->getNBlocks(); ++b)
      {
         if( compseeed->getNStairlinkingvars(b) != seeeds[i]->getNStairlinkingvars(b))
            noDuplicate = true;;
      }

      /** compares the number of constraints and variables in the blocks*/
      for( int j = 0; j < compseeed->getNBlocks() && !noDuplicate ; ++j )
      {
         if( (compseeed->getNVarsForBlock(j) != seeeds[i]->getNVarsForBlock(j)) || (compseeed->getNConssForBlock(j) != seeeds[i]->getNConssForBlock(j)) )
            noDuplicate = true;
      }

      /** sorts the the master conss, master vars, conss in blocks, vars in blocks, linking vars and stairlinking vars */
      if( sort && !noDuplicate)
      {
         compseeed->sort();
         seeeds[i]->sort();
      }

      /** compares the master cons */
      for( int j = 0; j < compseeed->getNMasterconss() && !noDuplicate; ++j)
      {
         if( compseeed->getMasterconss()[j] != seeeds[i]->getMasterconss()[j] )
            noDuplicate = true;
      }

      /** compares the master vars */
      for( int j = 0; j < compseeed->getNMastervars() && !noDuplicate; ++j)
      {
         if( compseeed->getMastervars()[j] != seeeds[i]->getMastervars()[j] )
            noDuplicate = true;
      }

      /** compares the constrains and variables in the blocks */
      for( int j = 0; j < compseeed->getNBlocks() && !noDuplicate; ++j )
      {
         for( int k = 0; k < compseeed->getNConssForBlock(j) && !noDuplicate; ++k)
         {
            if( compseeed->getConssForBlock(j)[k] != seeeds[i]->getConssForBlock(j)[k] )
               noDuplicate = true;
         }
         for( int k = 0; k < compseeed->getNVarsForBlock(j) && !noDuplicate; ++k)
         {
            if( compseeed->getVarsForBlock(j)[k] != seeeds[i]->getVarsForBlock(j)[k] )
               noDuplicate = true;
         }
      }

      /** compares the linking vars */
      for( int j = 0; j < compseeed->getNLinkingvars() && !noDuplicate; ++j)
      {
         if( compseeed->getLinkingvars()[j] != seeeds[i]->getLinkingvars()[j] )
            noDuplicate = true;
      }

      /** compares the stairlinking vars */
      for( int b = 0; b < compseeed->getNBlocks() && !noDuplicate; ++b)
      {
         for( int j = 0; j < compseeed->getNStairlinkingvars(b) && !noDuplicate; ++j)
         {
            if( compseeed->getStairlinkingvars(b)[j] != seeeds[i]->getStairlinkingvars(b)[j] )
               noDuplicate = true;
         }
      }

      if(!noDuplicate)
      {
         //std::cout << "seeed " << compseeed->getID() << " is a duplicate of seeed " << seeeds[i]->getID() << std::endl;
         if(compseeed->getHashValue() != seeeds[i]->getHashValue() )
         {
             compseeed->displaySeeed();
             seeeds[i]->displaySeeed();
         }
         assert(compseeed->getHashValue() == seeeds[i]->getHashValue() );
         return FALSE;
      }
      assert(compseeed->getHashValue() != seeeds[i]->getHashValue() );
   }
   return TRUE;
}

SCIP_Bool seeedIsNoDuplicate(SeeedPtr seeed, std::vector<SeeedPtr> const & currSeeeds, std::vector<SeeedPtr> const & finishedSeeeds, bool sort){
   return ( seeedIsNoDuplicateOfSeeeds(seeed, currSeeeds, sort) && seeedIsNoDuplicateOfSeeeds(seeed, finishedSeeeds, sort) );
}


/** constructor */
 Seeedpool::Seeedpool(
    SCIP*               givenScip, /**< SCIP data structure */
        const char*             conshdlrName
    ):scip(givenScip), currSeeeds(0), nTotalSeeeds(0),nVars(SCIPgetNVars(givenScip) ), nConss(SCIPgetNConss(givenScip) ), nDetectors(0), ndecompositions(0)
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
         conssForVars = std::vector<std::vector<int>>(nVars);

         assert(varToScipVar.size() == nVars);
         assert(consToScipCons.size() == nConss);

         /** assumption: now every relevant constraint and variable has its index and is stored in the corresponding unordered_map */
         /** find constraint <-> variable relationships and store them in both directions */
         for(int i = 0; i < (int)consToScipCons.size() ; ++i)
         {
                 SCIP_CONS* cons;
                 SCIP_VAR** currVars;
                 int            nCurrVars;
                 SCIP_Bool  success;

                 cons = consToScipCons[i];

//                 SCIP_CALL_ABORT( SCIPgetConsNVars(scip, cons, &nCurrVars, &success ) );
                 nCurrVars = GCGconsGetNVars(scip, cons);
                 std::cout << "\n\nConstraint: " << SCIPconsGetName(cons) << " with " << nCurrVars << " variables" << std::endl;
//                 assert(success);

                 SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &currVars, nCurrVars) ); /** free in line 321 */
//                 SCIPgetConsVars(scip, cons, currVars, nCurrVars, &success );
                 SCIP_CALL_ABORT(GCGconsGetVars(scip, cons, currVars, nCurrVars));

                 for(int currVar = 0; currVar < nCurrVars; ++currVar)
                 {

                     int varIndex;
                     std::tr1::unordered_map<SCIP_VAR*, int>::const_iterator iterVar;

                         std::cout << " try ("<< currVar << ")"<<varIndex << "/" << SCIPvarGetName(currVars[currVar]) << "\t";
                         if(!SCIPvarIsNegated(currVars[currVar]))
                         {
                             iterVar = scipVarToIndex.find(currVars[currVar]);
                         }
                         else
                         {
                            iterVar = scipVarToIndex.find(SCIPvarGetNegatedVar(currVars[currVar]));
                         }



                         if(iterVar == scipVarToIndex.end() )
                                 continue;


                         varIndex = iterVar->second;

                         varsForConss[i].push_back(varIndex);
                         conssForVars[varIndex].push_back(i);
                         std::cout << "("<< currVar << ")"<<varIndex << "/" << SCIPvarGetName(currVars[currVar]) << "\t";

                 }
                 std::cout << "\n" << std::endl;
                 SCIPfreeBufferArray(scip, &currVars) ;
         }

         /* populate seeed vector with empty seeed */

         currSeeeds.push_back(new Seeed(scip, nTotalSeeeds,nDetectors,nConss,nVars) );

         nTotalSeeeds++;

         decompositions = NULL;


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

         int maxRounds;
         SEEED_PROPAGATION_DATA* seeedPropData;
         SCIP_VAR* probvar;
         SCIP_CONS* cons;
         int cindex = 0;
         int vindex = 0;
         int currblock;
         bool displaySeeeds = false;
         int verboseLevel;
         std::vector<int> successDetectors;
         std::vector<SeeedPtr> delSeeeds;
         bool duplicate;

         successDetectors = std::vector<int>(nDetectors, 0);
         ndecompositions = 0;
         maxRounds = 2;
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

         for(int round = 0; round < maxRounds; ++round)
         {
                 std::cout << "currently in detection round " << round << std::endl;
                 std::vector<SeeedPtr> nextSeeeds = std::vector<SeeedPtr>(0);

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
                                 std::vector<SeeedPtr>::const_iterator newSIter;
                                 std::vector<SeeedPtr>::const_iterator newSIterEnd;

                                 SCIP_RESULT result = SCIP_DIDNOTFIND;

                                 /** if the seeed is also propageted by the detector go on with the next detector */
                                 if(seeedPtr->isPropagatedBy(d)  )
                                         continue;

                                 seeedPropData->seeedToPropagate = seeedPtr;

                                 /** new seeeds are created by the current detector */
                                 SCIP_CALL_ABORT( SCIPstartClock(scip, detectorToScipDetector[d]->dectime) );
                                 if(verboseLevel > 2)
                                     std::cout << "detector " << DECdetectorGetName(detectorToScipDetector[d]) << " started to propagate the " << s+1 << ". seeed (ID " << seeedPtr->getID() << ") in round " << round+1 << std::endl;
                                 SCIP_CALL_ABORT(detectorToScipDetector[d]->propagateSeeed(scip, detectorToScipDetector[d],seeedPropData, &result) );

                                 for( int j = 0; j < seeedPropData->nNewSeeeds; ++j )
                                 {
                                    seeedPropData->newSeeeds[j]->considerImplicits(this);
                                    seeedPropData->newSeeeds[j]->sort();
                                    seeedPropData->newSeeeds[j]->checkConsistency();
                                    seeedPropData->newSeeeds[j]->calcHashvalue();
                                 }

                                 if(seeedPropData->nNewSeeeds != 0 && displaySeeeds)
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

                                 SCIP_CALL_ABORT( SCIPstopClock(scip, detectorToScipDetector[d]->dectime) );


                                 assert(seeedPtr->isPropagatedBy(d));

                                 /** if the new seeeds are no duplicate they're added to the currSeeeds */
                                 for(int seeed = 0; seeed<seeedPropData->nNewSeeeds; ++seeed)
                                 {
                                         if(seeedIsNoDuplicate(seeedPropData->newSeeeds[seeed], currSeeeds, finishedSeeeds, false))
                                         {
                                            seeedPropData->newSeeeds[seeed]->calcOpenconss();
                                            seeedPropData->newSeeeds[seeed]->calcOpenvars();
                                            if(seeedPropData->newSeeeds[seeed]->getNOpenconss() == 0 && seeedPropData->newSeeeds[seeed]->getNOpenvars() == 0)
                                            {
                                               if(verboseLevel > 2)
                                                   std::cout << "seeed " << seeedPropData->newSeeeds[seeed]->getID() << " is addded to finished seeeds!" << std::endl;
                                               finishedSeeeds.push_back(seeedPropData->newSeeeds[seeed]);
                                            }
                                            else
                                            {
                                               if(verboseLevel > 2)
                                                   std::cout << "seeed " << seeedPropData->newSeeeds[seeed]->getID() << " is addded to next round seeeds!" << std::endl;
                                               nextSeeeds.push_back(seeedPropData->newSeeeds[seeed]);
                                            }
                                         }
                                         else
                                         {
                                            delete seeedPropData->newSeeeds[seeed];
                                         }
                                 }
                                 SCIPfreeMemoryArrayNull(scip, &seeedPropData->newSeeeds);

                         }

                         SCIP_CALL_ABORT(seeedPtr->completeGreedily( seeedPropData->seeedpool ) );
                         seeedPtr->calcHashvalue();


                   if(seeedIsNoDuplicateOfSeeeds(seeedPtr, finishedSeeeds, false))
                   {
                      finishedSeeeds.push_back(seeedPtr);
                   }
                 }

                 currSeeeds = nextSeeeds;

         }

         /* completeGreedily() on  currseeeds (from last round) and add them to finished seeeds */

         for(size_t i = 0; i < currSeeeds.size(); ++i)
         {
             SeeedPtr seeedPtr = currSeeeds[i];
             SCIP_CALL_ABORT(seeedPtr->completeGreedily( seeedPropData->seeedpool ) );
             seeedPtr->calcHashvalue();
             if(seeedIsNoDuplicateOfSeeeds(seeedPtr, finishedSeeeds, false))
             {
                finishedSeeeds.push_back(seeedPtr);
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
          //   std::vector<int>::const_iterator detectorIter =  finishedSeeeds.at(i)->detectorChain.begin();
          //   std::vector<int>::const_iterator detectorIterEnd =  finishedSeeeds.at(i)->detectorChain.end();

          //   for (; detectorIter != detectorIterEnd; ++detectorIter)
          //   {
          //       successDetectors[*detectorIter]++;
          //   }

             for(size_t d = 0; d < nDetectors; ++d)
             {
                 if(finishedSeeeds[i]->isPropagatedBy(d))
                     successDetectors[d] += 1;
             }
         }

         /** preliminary output detector stats */

         std::cout << "Begin preliminary detector times: " << std::endl;

         for(size_t i = 0; i < nDetectors; ++i)
         {
             std::cout << "Detector " << DECdetectorGetName(detectorToScipDetector[i] ) << " worked on " << successDetectors[i] << " of " << finishedSeeeds.size() << " and took a total time of " << SCIPgetClockTime(scip, detectorToScipDetector[i]->dectime)  << std::endl;
         }


         /** fill out the decompositions */

         SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &decompositions, (int) finishedSeeeds.size())); /** free in decomp.c:470 */
         for( size_t i = 0; i < finishedSeeeds.size(); ++i )
         {
            SeeedPtr seeed = finishedSeeeds[i];

            /** set nblocks */
            int nblocks = seeed->getNBlocks();
            SCIP_CALL_ABORT( DECdecompCreate(scip, &(decompositions[i])) );
            decompositions[i]->nblocks = nblocks;

            /** set linkingconss */
            int nlinkingconss = seeed->getNMasterconss();
            decompositions[i]->nlinkingconss = nlinkingconss;
            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->linkingconss, nlinkingconss) ); /** free in decomp.c:468 */
            for( int j = 0; j < nlinkingconss; ++j )
            {
               decompositions[i]->linkingconss[j] = getConsForIndex(seeed->getMasterconss()[j]);
            }

            /** set linkingvars */
            int nlinkingvars = seeed->getNLinkingvars() + seeed->getNMastervars();
            decompositions[i]->nlinkingvars = nlinkingvars;
            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->linkingvars, nlinkingvars) ); /** free in decomp.c:465 */
            for( int j = 0; j < seeed->getNLinkingvars(); ++j )
            {
               decompositions[i]->linkingvars[j] = SCIPvarGetProbvar( getVarForIndex(seeed->getLinkingvars()[j]) );
            }
            for( int j = 0; j < seeed->getNMastervars(); ++j )
            {
               decompositions[i]->linkingvars[j + seeed->getNLinkingvars()] = SCIPvarGetProbvar( getVarForIndex(seeed->getMastervars()[j]));
            }


            /** set stairlinkingvars */
            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->stairlinkingvars, nblocks) ); /** free in decomp.c:466 */
            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->nstairlinkingvars, nblocks) ); /** free in decomp.c:467 */
            for ( int j = 0; j < nblocks; ++j )
            {
               int nstairlinkingvars = seeed->getNStairlinkingvars(j);
               decompositions[i]->nstairlinkingvars[j]=nstairlinkingvars;
               SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->stairlinkingvars[j], nstairlinkingvars) ); /** free in decomp.c:444 */
               for ( int k = 0; k < nstairlinkingvars; ++k )
               {
                  decompositions[i]->stairlinkingvars[j][k] = SCIPvarGetProbvar( getVarForIndex(seeed->getStairlinkingvars(j)[k]) );
               }
            }

            /** create hashmap constoblock and consindex */
            SCIP_CALL_ABORT( SCIPhashmapCreate( &decompositions[i]->constoblock, SCIPblkmem(scip), nVars ) );
            SCIP_CALL_ABORT( SCIPhashmapCreate( &decompositions[i]->consindex, SCIPblkmem(scip), nVars ) );


            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->subscipconss, nblocks) ); /** free in decomp.c:463 */
            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->nsubscipconss, nblocks) ); /** free in decomp.c:464 */
            /** add block conss */
            for( int b = 0; b < nblocks; ++b )
            {
               currblock = b+1;
               decompositions[i]->nsubscipconss[b] = seeed->getNConssForBlock(b);
               SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->subscipconss[b], decompositions[i]->nsubscipconss[b]) ); /** free in decomp.c:416 */
               for ( int j = 0; j < seeed->getNConssForBlock(b); ++j)
               {
                  cindex++;
                  cons = getConsForIndex((seeed->getConssForBlock(b))[j]);
                  decompositions[i]->subscipconss[b][j] = cons;
                  SCIPdebugMessage("cons %d is cons of block %d\n", cindex, currblock );
                  SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->constoblock, cons, (void*) (size_t) currblock) );
                  SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->consindex, (void*) (size_t) cindex, (void*) (size_t) currblock) );
                  SCIP_CALL_ABORT( SCIPcaptureCons(scip, cons) );
               }
            }

            /** add master conss */
            for( int j = 0; j < seeed->getNMasterconss(); ++j )
            {
               cindex++;
               cons = getConsForIndex(seeed->getMasterconss()[j]);
               SCIPdebugMessage("cons %d is mastercons\n", cindex);
               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->constoblock, cons, (void*) (size_t) (nblocks+1)) );
               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->consindex, (void*) (size_t) cindex, (void*) (size_t) (nblocks+1)) );
               SCIP_CALL_ABORT( SCIPcaptureCons(scip, cons) );
            }

            /** create hashmap vartoblock and varindex */
            SCIP_CALL_ABORT( SCIPhashmapCreate( &decompositions[i]->vartoblock, SCIPblkmem(scip), nVars ) );
            SCIP_CALL_ABORT( SCIPhashmapCreate( &decompositions[i]->varindex, SCIPblkmem(scip), nVars ) );

            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->subscipvars, nblocks) ); /** free in decomp.c:461 */
            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->nsubscipvars, nblocks) ); /** free in decomp.c:462 */
            /** add block vars and stairlinkingvars */
            for( int b = 0; b < nblocks; ++b)
            {
               currblock = b+1;
               decompositions[i]->nsubscipvars[b] = seeed->getNVarsForBlock(b);
               SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->subscipvars[b], decompositions[i]->nsubscipvars[b]) ); /** free in decomp.c:405 */
               for ( int j = 0; j < seeed->getNVarsForBlock(b); ++j )
               {
                  vindex++;
                  probvar = SCIPvarGetProbvar( getVarForIndex(seeed->getVarsForBlock(b)[j]) );
                  decompositions[i]->subscipvars[b][j] = probvar;
                  if( SCIPhashmapExists(decompositions[i]->vartoblock, probvar) )
                  {
                     SCIPdebugMessage("var <%s> has been handled before, it should not been add to block %d\n", SCIPvarGetName(probvar), currblock);
                  }
                  else
                  {
                     SCIPdebugMessage("var <%s> has not been handled before, adding to block %d\n", SCIPvarGetName(probvar), currblock );
                     SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->vartoblock, probvar, (void*) (size_t) currblock) );
                     SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->varindex, (void*) (size_t) vindex, (void*) (size_t) currblock) );
                     SCIP_CALL_ABORT( SCIPcaptureVar(scip, probvar) );
                  }
               }

               for( int k = 0; k < seeed->getNStairlinkingvars(b); ++k )
               {
                  probvar = SCIPvarGetProbvar( getVarForIndex((seeed->getStairlinkingvars(b))[k]));
                  if( !SCIPhashmapExists(decompositions[i]->vartoblock, probvar) )
                  {
                     vindex++;
                     SCIPdebugMessage("var <%s> is stairlinkingvar\n", SCIPvarGetName(probvar));
                     SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->vartoblock, probvar, (void*) (size_t) (nblocks+2)) );
                     SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->varindex, (void*) (size_t) vindex, (void*) (size_t) (nblocks+2)) );
                     SCIP_CALL_ABORT( SCIPcaptureVar(scip, probvar) );
                  }
                  assert( ((size_t) SCIPhashmapGetImage(decompositions[i]->vartoblock, probvar)) == nblocks+2);
               }
            }

            /** add linking vars */
            for( int j = 0; j < seeed->getNLinkingvars(); ++j )
            {
               vindex++;
               SCIPdebugMessage("var <%s> is linkingvar\n", SCIPvarGetName(probvar));
               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->vartoblock, probvar, (void*) (size_t) (nblocks+2)) );
               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->varindex, (void*) (size_t) vindex, (void*) (size_t) (nblocks+2)) );
               SCIP_CALL_ABORT( SCIPcaptureVar(scip, probvar) );
            }

            /** add master vars */
            for( int j = 0; j < seeed->getNMastervars(); ++j )
            {
               vindex++;
               SCIPdebugMessage("var <%s> is mastervar\n", SCIPvarGetName(probvar));
               probvar = SCIPvarGetProbvar( getVarForIndex(seeed->getMastervars()[j]) );
               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->vartoblock, probvar, (void*) (size_t) (nblocks+1)) );
               SCIP_CALL_ABORT( SCIPhashmapSetImage(decompositions[i]->varindex, (void*) (size_t) vindex, (void*) (size_t) (nblocks+1)) );
               SCIP_CALL_ABORT( SCIPcaptureVar(scip, probvar) );
            }

            /** set detectorchain */
            int ndetectors = seeed->getNDetectors();
            decompositions[i]->sizeDetectorchain = ndetectors;
            SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &decompositions[i]->detectorchain, ndetectors) ); /** free in decomp.c:469 */
            for( int k = 0; k < ndetectors; ++k )
            {
               decompositions[i]->detectorchain[k] = getDetectorForIndex(seeed->getDetectorchain()[k]);
            }

            /** set dectype */
            if(decompositions[i]->nlinkingvars > 0 && decompositions[i]->nlinkingconss == 0)
            {
               decompositions[i]->type = DEC_DECTYPE_STAIRCASE;
            }
            else if(decompositions[i]->nlinkingvars > 0)
            {
               decompositions[i]->type = DEC_DECTYPE_ARROWHEAD;
            }
            else if(decompositions[i]->nlinkingconss > 0)
            {
               decompositions[i]->type = DEC_DECTYPE_BORDERED;
            }
            else if(decompositions[i]->nlinkingconss == 0)
            {
               decompositions[i]->type = DEC_DECTYPE_DIAGONAL;
            }
            else
            {
               decompositions[i]->type = DEC_DECTYPE_UNKNOWN;
            }

            ndecompositions++;
            assert(!SCIPhashmapIsEmpty(decompositions[i]->constoblock));
            assert(!SCIPhashmapIsEmpty(decompositions[i]->vartoblock));
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

         delete seeedPropData;
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




} /* namespace gcg */
