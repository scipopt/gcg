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

/**@file   dec_varclass.cpp
 * @ingroup DETECTORS
 * @brief  detector varclass
 * @author Julius Hense
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_varclass.h"
#include "cons_decomp.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "class_varclassifier.h"
#include "gcg.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "scip_misc.h"
#include "scip/clock.h"

#include <sstream>

#include <iostream>
#include <algorithm>

/* constraint handler properties */
#define DEC_DETECTORNAME          "varclass"       /**< name of detector */
#define DEC_DESC                  "detector varclass" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          0           /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               'v'         /**< display character of detector */
#define DEC_ENABLED               TRUE        /**< should the detection be enabled */
#define DEC_ENABLEDORIGINAL       TRUE        /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING      FALSE        /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */
#define DEC_LEGACYMODE            FALSE       /**< should (old) DETECTSTRUCTURE method also be used for detection */

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

/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
#define freeVarclass NULL

/** destructor of detector to free detector data (called before the solving process begins) */
#if 0
static
DEC_DECL_EXITDETECTOR(exitVarclass)
{ /*lint --e{715}*/

   SCIPerrorMessage("Exit function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define exitVarclass NULL
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
#define initVarclass NULL
#endif

/** detection function of detector */
//static DEC_DECL_DETECTSTRUCTURE(detectVarclass)
//{ /*lint --e{715}*/
//   *result = SCIP_DIDNOTFIND;
//
//   SCIPerrorMessage("Detection function of detector <%s> not implemented!\n", DEC_DETECTORNAME)
//;   SCIPABORT(); /*lint --e{527}*/
//
//   return SCIP_OKAY;
//}

#define detectVarclass NULL

#define finishSeeedVarclass NULL

static DEC_DECL_PROPAGATESEEED(propagateSeeedVarclass)
{
  *result = SCIP_DIDNOTFIND;
  char decinfo[SCIP_MAXSTRLEN];

  SCIP_CLOCK* temporaryClock;

  if (seeedPropagationData->seeedToPropagate->getNOpenconss() != seeedPropagationData->seeedpool->getNConss() ||  seeedPropagationData->seeedToPropagate->getNOpenvars() != seeedPropagationData->seeedpool->getNVars() )
  {
    *result = SCIP_SUCCESS;
     return SCIP_OKAY;
  }

  SCIP_CALL_ABORT( SCIPcreateClock(scip, &temporaryClock) );
  SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );

  std::vector<gcg::Seeed*> foundseeeds(0);

  gcg::Seeed* seeedOrig;
  gcg::Seeed* seeed;

  int maximumnclasses;

  if( seeedPropagationData->seeedpool->getNConss() + seeedPropagationData->seeedpool->getNVars() >= 50000 )
        SCIPgetIntParam(scip, "detection/maxnclassesperclassifierforlargeprobs", &maximumnclasses);
     else
        SCIPgetIntParam(scip, "detection/maxnclassesperclassifier", &maximumnclasses);

  for( int classifierIndex = 0; classifierIndex < seeedPropagationData->seeedpool->getNVarClassifiers(); ++classifierIndex )
  {
    gcg::VarClassifier* classifier = seeedPropagationData->seeedpool->getVarClassifier( classifierIndex );
    std::vector<int> varclassindices_master = std::vector<int>(0);
    std::vector<int> varclassindices_linking = std::vector<int>(0);

    if ( classifier->getNClasses() > maximumnclasses )
    {
       std::cout << " the current varclass distribution includes " <<  classifier->getNClasses() << " classes but only " << maximumnclasses << " are allowed for propagateSeeed() of var class detector" << std::endl;
       continue;
    }

    seeedOrig = seeedPropagationData->seeedToPropagate;

    for( int i = 0; i < classifier->getNClasses(); ++ i )
    {
       switch( classifier->getClassDecompInfo( i ) )
       {
          case gcg::ALL:
             break;
          case gcg::LINKING:
             varclassindices_linking.push_back( i );
             break;
          case gcg::MASTER:
             varclassindices_master.push_back( i );
             break;
          case gcg::BLOCK:
             break;
       }
    }

    std::vector< std::vector<int> > subsetsOfVarclasses = classifier->getAllSubsets( true, false, false, false );

    for( size_t subset = 0; subset < subsetsOfVarclasses.size(); ++subset )
    {
       if( subsetsOfVarclasses[subset].size() == 0 && varclassindices_master.size() == 0 && varclassindices_linking.size() == 0 )
          continue;

       seeed = new gcg::Seeed(seeedOrig);

       /** book open vars that have a) type of the current subset or b) decomp info LINKING as linking vars */
       for( int i = 0; i < seeed->getNOpenvars(); ++i )
       {
          bool foundVar = false;
          for( size_t varclassId = 0; varclassId < subsetsOfVarclasses[subset].size(); ++varclassId )
          {
              if( classifier->getClassOfVar( seeed->getOpenvars()[i] ) == subsetsOfVarclasses[subset][varclassId] )
              {
                  seeed->bookAsLinkingVar(seeed->getOpenvars()[i]);
                  foundVar = true;
                  break;
              }
          }
          /** only check varclassindices_linking if current var has not already been found in a subset */
          if ( !foundVar )
          {
             for( size_t varclassId = 0; varclassId < varclassindices_linking.size(); ++varclassId )
             {
                if( classifier->getClassOfVar( seeed->getOpenvars()[i] ) == varclassindices_linking[varclassId] )
                {
                   seeed->bookAsLinkingVar(seeed->getOpenvars()[i]);
                   foundVar = true;
                   break;
                }
             }
          }
          /** only check varclassindices_master if current var has not already been found in a subset */
          if ( !foundVar )
          {
             for( size_t varclassId = 0; varclassId < varclassindices_master.size(); ++varclassId )
             {
                if( classifier->getClassOfVar( seeed->getOpenvars()[i] ) == varclassindices_master[varclassId] )
                {
                   seeed->bookAsMasterVar(seeed->getOpenvars()[i]);
                   break;
                }
             }
          }
       }

       /** set decinfo to: varclass_<classfier_name>:<linking_class_name#1>-...-<linking_class_name#n> */
       std::stringstream decdesc;
       decdesc << "varclass" << "\\_" << classifier->getName() << ": \\\\ ";
       std::vector<int> curlinkingclasses( varclassindices_linking );
       for ( size_t varclassId = 0; varclassId < subsetsOfVarclasses[subset].size(); ++varclassId )
       {
          if ( varclassId > 0 )
          {
             decdesc << "-";
          }
          decdesc << classifier->getClassName( subsetsOfVarclasses[subset][varclassId] );

          if( std::find( varclassindices_linking.begin(), varclassindices_linking.end(),
             subsetsOfVarclasses[subset][varclassId] ) == varclassindices_linking.end() )
          {
             curlinkingclasses.push_back( subsetsOfVarclasses[subset][varclassId] );
          }
       }
       for ( size_t varclassId = 0; varclassId < varclassindices_linking.size(); ++varclassId )
       {
          if ( varclassId > 0 || subsetsOfVarclasses[subset].size() > 0)
          {
             decdesc << "-";
          }
          decdesc << classifier->getClassName( varclassindices_linking[varclassId] );
       }

       seeed->flushBooked();
       (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, decdesc.str().c_str());
       seeed->addDetectorChainInfo(decinfo);
       seeed->setVarClassifierStatistics( seeed->getNDetectors(), classifier, curlinkingclasses,
          varclassindices_master );

       foundseeeds.push_back(seeed);
    }
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
DEC_DECL_PROPAGATEFROMTOOLBOX(propagateFromToolboxVarclass)
{
   *result = SCIP_DIDNOTFIND;
   char decinfo[SCIP_MAXSTRLEN];
   gcg::VarClassifier** classifiers;
   int nclassifiers;
   SCIP_Bool newclass, newclassifier;
   gcg::VarClassifier* selectedclassifier;
   std::vector<int> selectedclasses;
   int i, j;
   char stri[SCIP_MAXSTRLEN];
   SCIP_Bool finished;

   char* command;
   int commandlen;
   SCIP_Bool endoffile;

   if (seeedPropagationData->seeedToPropagate->getNOpenconss() != seeedPropagationData->seeedpool->getNConss() ||  seeedPropagationData->seeedToPropagate->getNOpenvars() != seeedPropagationData->seeedpool->getNVars() )
   {
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   if( seeedPropagationData->seeedpool->getNVarClassifiers() == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No VarClassifiers listed for propagation, starting classification.\n");
      seeedPropagationData->seeedpool->calcClassifierAndNBlockCandidates(scip);
      if( seeedPropagationData->seeedpool->getNVarClassifiers() == 0 )
      {
         SCIPinfoMessage(scip, NULL, "No VarClassifiers found after calculation, aborting!.\n");
         return SCIP_ERROR;
      }
   }

   std::vector<gcg::Seeed*> foundseeeds(0);

   gcg::Seeed* seeedOrig;
   gcg::Seeed* seeed;

   int maximumnclasses;

   if( seeedPropagationData->seeedpool->getNConss() + seeedPropagationData->seeedpool->getNVars() >= 50000 )
         SCIPgetIntParam(scip, "detection/maxnclassesperclassifierforlargeprobs", &maximumnclasses);
      else
         SCIPgetIntParam(scip, "detection/maxnclassesperclassifier", &maximumnclasses);

   SCIP_CALL( SCIPallocMemoryArray(scip, &classifiers, seeedPropagationData->seeedpool->getNVarClassifiers()) );

   if( seeedPropagationData->seeedpool->getNVarClassifiers() == 0 )
   {
       SCIPinfoMessage(scip, NULL, "No varclassifiers available for propagation, aborting!\n");
       *result = SCIP_DIDNOTFIND;
       return SCIP_OKAY;
   }

   SCIPinfoMessage(scip, NULL, "\n%d varclassifiers available for propagation.\n", seeedPropagationData->seeedpool->getNVarClassifiers() );

   nclassifiers = 0;
   for( int classifierIndex = 0; classifierIndex < seeedPropagationData->seeedpool->getNVarClassifiers(); ++classifierIndex )
   {
      gcg::VarClassifier* classifier = seeedPropagationData->seeedpool->getVarClassifier( classifierIndex );
      if( classifier->getNClasses() > maximumnclasses )
      {
         std::cout << " the current varclass distribution includes " << classifier->getNClasses() << " classes but only " << maximumnclasses << " are allowed for propagateSeeed() of cons class detector" << std::endl;
         continue;
      }
      newclassifier = TRUE;
      for( i = 0; i < nclassifiers; ++i )
      {
         if( classifiers[i] == classifier )
         {
            newclassifier = FALSE;
         }
      }

      if( newclassifier )
      {
         classifiers[nclassifiers] = classifier;
         ++nclassifiers;
      }
   }

   selectedclassifier = classifiers[0]; //default case to omit warnings
   /* user selects a consclassifier */
   finished = FALSE;
   while( !finished )
   {
      SCIPinfoMessage(scip, NULL, "Available consclassifiers:\n");
      for( i = 0; i < nclassifiers; ++i )
      {
         SCIPinfoMessage(scip, NULL, "%d) ", i+1); //+1 as we want the list to start with 1)
         SCIPinfoMessage(scip, NULL, "%s\n", classifiers[i]->getName() );
      }

      commandlen = 0;
      while( commandlen == 0 )
      {
         SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Type in the name or number of the varclassifier that you want to use (seperated by spaces) or \"done\", (use \"quit\" to exit detector): \nGCG/toolbox> ", &command, &endoffile) );
         commandlen = strlen(command);
      }

      if( !strncmp( command, "done", commandlen) == 0 && !strncmp( command, "quit", commandlen) == 0 )
      {
         for( i = 0; i < nclassifiers; ++i )
         {
            sprintf(stri, "%d", i+1); //used for matching numberings in the list, off-by-one since classifiers array starts with index 0 and our list with 1)
            if( strncmp( command, classifiers[i]->getName(), commandlen) == 0 || strncmp( command, stri, commandlen ) == 0 )
            {
               selectedclassifier = classifiers[i];
               finished = TRUE;
               continue;
            }
         }
      }
      else if( strncmp( command, "done", commandlen) == 0 )
      {
         finished = TRUE;
         continue;
      }
      else if( strncmp( command, "quit", commandlen) == 0 )
      {
         SCIPfreeMemoryArray(scip, &classifiers);
         *result = SCIP_DIDNOTFIND;
         return SCIP_OKAY;
      }
   }

   std::vector<int> varclassindices = std::vector<int>(0);
   for( i = 0; i < selectedclassifier->getNClasses(); ++i )
   {
      varclassindices.push_back(i);
   }

   SCIPinfoMessage(scip, NULL, "You will now be asked to enter a selection of classes iteratively. If you have finished your selection, enter \"done\".\n");
   finished = FALSE;
   while( !finished )
   {
      std::vector<int> nVarsOfClasses = selectedclassifier->getNVarsOfClasses();
      SCIPinfoMessage(scip, NULL, "The following classes are available for the selected varclassifier \"%s\":\n",selectedclassifier->getName());
      for( i = 0; i < static_cast<int>(varclassindices.size()); ++i )
      {
         SCIPinfoMessage(scip, NULL, "%d) ", i+1); //+1 as we want the list to start with 1)
         SCIPinfoMessage(scip, NULL, "%s || NConss: %d || %s\n", selectedclassifier->getClassName(varclassindices[i]), nVarsOfClasses[i], selectedclassifier->getClassDescription(varclassindices[i]));
      }

      commandlen = 0;
      while( commandlen == 0 )
      {
         SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Type in the name(s) or number(s) of classes (seperated by spaces) or \"done\", (use \"quit\" to exit detector): \nGCG/toolbox> ", &command, &endoffile) );
         commandlen = strlen(command);
      }

      if( !strncmp( command, "done", commandlen) == 0 && !strncmp( command, "quit", commandlen) == 0 )
      {
         for( i = 0; i < static_cast<int>(varclassindices.size()); ++i )
         {
            newclass = TRUE;
            sprintf(stri, "%d", i+1); //used for matching numberings in the list, off-by-one since classifiers array starts with index 0 and our list with 1)
            if( strncmp( command, selectedclassifier->getClassNameOfVar(varclassindices[i]), commandlen) == 0 || strncmp( command, stri, commandlen ) == 0 )
            {
               /* check that we do not select the same classifier multiple times*/
               for( j = 0; j < static_cast<int>(selectedclasses.size()); ++j )
               {
                  if( selectedclasses[j] == varclassindices[i] )
                  {
                     newclass = FALSE;
                  }
               }
               if( newclass )
               {
                  selectedclasses.push_back(varclassindices[i]);

                  SCIPinfoMessage(scip, NULL, "\nCurrently selected classifiers: ");
                  for( j = 0; j < static_cast<int>(selectedclasses.size()); ++j )
                  {
                     SCIPinfoMessage(scip, NULL, "\"%s\" ", selectedclassifier->getClassNameOfVar(selectedclasses[j]));
                  }
                  SCIPinfoMessage(scip, NULL, "\n\n");

                  if( selectedclasses.size() >= varclassindices.size() )
                  {
                     finished = TRUE;
                     break;
                  }
               }
               else
               {
                  SCIPinfoMessage(scip, NULL, "\n+++Class \"%s\" is already selected!+++\n\n", selectedclassifier->getClassNameOfVar(varclassindices[i]));
               }
            }
         }
      }
      else if( strncmp( command, "done", commandlen) == 0 )
      {
         finished = TRUE;
         continue;
      }
      else if( strncmp( command, "quit", commandlen) == 0 )
      {
         SCIPfreeMemoryArray(scip, &classifiers);
         *result = SCIP_DIDNOTFIND;
         return SCIP_OKAY;
      }
   }

   std::vector<int> varclassindices_master = std::vector<int>(0);
   std::vector<int> varclassindices_linking = std::vector<int>(0);

   seeedOrig = seeedPropagationData->seeedToPropagate;

   for( i = 0; i < selectedclassifier->getNClasses(); ++ i )
   {
      switch( selectedclassifier->getClassDecompInfo( i ) )
      {
         case gcg::ALL:
            break;
         case gcg::LINKING:
            varclassindices_linking.push_back( i );
            break;
         case gcg::MASTER:
            varclassindices_master.push_back( i );
            break;
         case gcg::BLOCK:
            break;
      }
   }

   if( selectedclasses.size() == 0 && varclassindices_master.size() == 0 && varclassindices_linking.size() == 0 )
   {
      *result = SCIP_DIDNOTFIND;
      SCIPfreeMemoryArray(scip, &classifiers);
      return SCIP_OKAY;
   }

   seeed = new gcg::Seeed(seeedOrig);

   /** book open vars that have a) type of the current subset or b) decomp info LINKING as linking vars */
   for( i = 0; i < seeed->getNOpenvars(); ++i )
   {
      bool foundVar = false;
      for( size_t varclassId = 0; varclassId < selectedclasses.size(); ++varclassId )
      {
         if( selectedclassifier->getClassOfVar( seeed->getOpenvars()[i] ) == selectedclasses[varclassId] )
         {
            seeed->bookAsLinkingVar(seeed->getOpenvars()[i]);
            foundVar = true;
            break;
         }
      }
      /** only check varclassindices_linking if current var has not already been found in a subset */
      if ( !foundVar )
      {
         for( size_t varclassId = 0; varclassId < varclassindices_linking.size(); ++varclassId )
         {
            if( selectedclassifier->getClassOfVar( seeed->getOpenvars()[i] ) == varclassindices_linking[varclassId] )
            {
               seeed->bookAsLinkingVar(seeed->getOpenvars()[i]);
               foundVar = true;
               break;
            }
         }
      }
      /** only check varclassindices_master if current var has not already been found in a subset */
      if ( !foundVar )
      {
         for( size_t varclassId = 0; varclassId < varclassindices_master.size(); ++varclassId )
         {
            if( selectedclassifier->getClassOfVar( seeed->getOpenvars()[i] ) == varclassindices_master[varclassId] )
            {
               seeed->bookAsMasterVar(seeed->getOpenvars()[i]);
               break;
            }
         }
      }
   }

   /** set decinfo to: varclass_<classfier_name>:<linking_class_name#1>-...-<linking_class_name#n> */
   std::stringstream decdesc;
   decdesc << "varclass" << "\\_" << selectedclassifier->getName() << ": \\\\ ";
   std::vector<int> curlinkingclasses( varclassindices_linking );
   for ( size_t varclassId = 0; varclassId < selectedclasses.size(); ++varclassId )
   {
      if ( varclassId > 0 )
      {
         decdesc << "-";
      }
      decdesc << selectedclassifier->getClassName( selectedclasses[varclassId] );

      if( std::find( varclassindices_linking.begin(), varclassindices_linking.end(),
         selectedclasses[varclassId] ) == varclassindices_linking.end() )
      {
         curlinkingclasses.push_back( selectedclasses[varclassId] );
      }
   }
   for ( size_t varclassId = 0; varclassId < varclassindices_linking.size(); ++varclassId )
   {
      if ( varclassId > 0 || selectedclasses.size() > 0)
      {
         decdesc << "-";
      }
      decdesc << selectedclassifier->getClassName( varclassindices_linking[varclassId] );
   }

   seeed->flushBooked();
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, decdesc.str().c_str());
   seeed->addDetectorChainInfo(decinfo);
   seeed->setDetectorPropagated(detector);
   seeed->setVarClassifierStatistics( seeed->getNDetectors(), selectedclassifier, curlinkingclasses,
      varclassindices_master );

   foundseeeds.push_back(seeed);



   //@TODO: This alloc is already done in cons_decomp:SCIPconshdlrDecompToolboxActOnSeeed(..). 
   //This is contrary to the behaviour of other detectors such as hrcg, hrg and hr but not connectedbase
   SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), foundseeeds.size() ) );
   seeedPropagationData->nNewSeeeds = foundseeeds.size();

   for( int s = 0; s < seeedPropagationData->nNewSeeeds; ++s )
   {
      seeedPropagationData->newSeeeds[s] = foundseeeds[s];
   }

   *result = SCIP_SUCCESS;
   SCIPfreeMemoryArray(scip, &classifiers);
   return SCIP_OKAY;
}

#define finishFromToolboxVarclass NULL

#define detectorPostprocessSeeedVarclass NULL

static
DEC_DECL_SETPARAMAGGRESSIVE(setParamAggressiveVarclass)
{
   char setstr[SCIP_MAXSTRLEN];
   SCIP_Real modifier;

   int newval;
   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
      return SCIP_OKAY;

   modifier = ((SCIP_Real)SCIPgetNConss(scip) + (SCIP_Real)SCIPgetNVars(scip) ) / SET_MULTIPLEFORSIZETRANSF;
   modifier = log(modifier) / log(2.);

   if (!SCIPisFeasPositive(scip, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(scip, modifier);

   newval = MAX( 2, AGGRESSIVE_MAXIMUMNCLASSES - modifier );
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnclasses", name);

   SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );
   SCIPinfoMessage(scip, NULL, "\n%s = %d\n", setstr, newval);


   return SCIP_OKAY;

}


static
DEC_DECL_SETPARAMDEFAULT(setParamDefaultVarclass)
{
   char setstr[SCIP_MAXSTRLEN];
   SCIP_Real modifier;

   int newval;
   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLED) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLEDORIGINAL ) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLEDFINISHING ) );

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
         return SCIP_OKAY;

   modifier = ( (SCIP_Real)SCIPgetNConss(scip) + (SCIP_Real)SCIPgetNVars(scip) ) / SET_MULTIPLEFORSIZETRANSF;
   modifier = log(modifier) / log(2);

   if (!SCIPisFeasPositive(scip, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(scip, modifier);

   newval = MAX( 2, DEFAULT_MAXIMUMNCLASSES - modifier );
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnclasses", name);

   SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );
   SCIPinfoMessage(scip, NULL, "\n%s = %d\n", setstr, newval);

   return SCIP_OKAY;

}

static
DEC_DECL_SETPARAMFAST(setParamFastVarclass)
{
   char setstr[SCIP_MAXSTRLEN];
   SCIP_Real modifier;
   int newval;

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
         return SCIP_OKAY;

   modifier = ( (SCIP_Real)SCIPgetNConss(scip) + (SCIP_Real)SCIPgetNVars(scip) ) / SET_MULTIPLEFORSIZETRANSF;

   modifier = log(modifier) / log(2);

   if (!SCIPisFeasPositive(scip, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(scip, modifier);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnclasses", name);

   newval = MAX( 2, FAST_MAXIMUMNCLASSES - modifier );

   SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );
   SCIPinfoMessage(scip, NULL, "\n%s = %d\n", setstr, newval);

   return SCIP_OKAY;

}



/*
 * detector specific interface methods
 */

/** creates the handler for varclass detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorVarclass(SCIP* scip /**< SCIP data structure */
)
{
   DEC_DETECTORDATA* detectordata;
   char setstr[SCIP_MAXSTRLEN];

   /**@todo create varclass detector data here*/
   detectordata = NULL;

   SCIP_CALL(
      DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND,
         DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDORIGINAL, DEC_ENABLEDFINISHING,DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, DEC_LEGACYMODE, detectordata, detectVarclass,
         freeVarclass, initVarclass, exitVarclass, propagateSeeedVarclass, propagateFromToolboxVarclass, finishFromToolboxVarclass, finishSeeedVarclass, detectorPostprocessSeeedVarclass, setParamAggressiveVarclass, setParamDefaultVarclass, setParamFastVarclass));

   /**@todo add varclass detector parameters */

   const char* name = DEC_DETECTORNAME;
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnclasses", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, "maximum number of classes ",  NULL, FALSE, DEFAULT_MAXIMUMNCLASSES, 1, INT_MAX, NULL, NULL ) );

   return SCIP_OKAY;
}
