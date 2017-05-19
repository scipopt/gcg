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

/*lint -e64 disable useless and wrong lint warning */

#ifdef __INTEL_COMPILER
#ifndef _OPENMP
#pragma warning disable 3180  /* disable wrong and useless omp warnings */
#endif
#endif

//#define SCIP_DEBUG

#include "gcg.h"
#include "objscip/objscip.h"
#include "class_seeedpool.h"
#include "struct_detector.h"
#include "pub_decomp.h"
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
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif


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


#define ENUM_TO_STRING(x) # x
#define DEFAULT_THREADS                  0          /**< number of threads (0 is OpenMP default) */


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
   int                   nthreads;
};



namespace gcg {


/** local methods */

struct sort_decr
{
   bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right)
   {
      return left.second > right.second;
   }
};


struct sort_pred {
   bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
      return left.second < right.second;
   }
};

std::string getSeeedFolderLatex( SeeedPtr seeed )
{
   std::stringstream decompfilename;
   decompfilename << "dec" << seeed->getID() << ".pdf";

   return decompfilename.str();
}


SCIP_Bool unfinishedchildexists(std::vector<SCIP_Bool> const& childsfinished)
{
   for( size_t s = 0; s < childsfinished.size(); ++s )
   {
      if( !childsfinished[s] )
         return true;
   }
   return false;
}

int getfirstunfinishedchild(std::vector<SCIP_Bool> const& childsfinished, std::vector<int> const& childs)
{
   for( size_t s = 0; s < childsfinished.size(); ++s )
   {
      if( !childsfinished[s] )
         return childs[s];
   }
   return -1;
}

int getfirstunfinishedchildid(std::vector<SCIP_Bool> const& childsfinished, std::vector<int> const& childs)
{
   for( size_t s = 0; s < childsfinished.size(); ++s )
   {
      if( !childsfinished[s] )
         return (int)s;
   }
   return -1;
}


/**
 * @return is nextchild the last unfinished child
 */
SCIP_Bool finishnextchild( std::vector<int>& childs, std::vector<SCIP_Bool>& childsfinished, int child )
{
   for( size_t s = 0; s < childsfinished.size(); ++s )
   {
      if( !childsfinished[s] )
      {
         assert(childs[s] == child);
         childsfinished[s] = TRUE;
         return s == childsfinished.size() - 1;
      }
   }
   return FALSE;
}

std::string writeSeeedDetectorChainInfoLatex( SeeedPtr seeed, int currheight, int visucounter )
{
   std::stringstream line;
   std::string relposition;
   int position = visucounter % 3;
   if( position == 0 )
      relposition = "above";
   else if ( position == 1)
      relposition = "";
   else if ( position == 2)
      relposition = "below";
   else
      relposition = "below left";

   if ( currheight != 1)
      relposition = "";

   if ( (size_t) currheight >  seeed->detectorchaininfo.size() )
      line << "edge from parent node [" << relposition << "] {no info" << seeed->getID() << "-" << currheight -1 << " } " ;
   else
   {
      std::string oldinfo = seeed->detectorchaininfo[ currheight - 1];
      /** take latexified detctorchaininfo */
      size_t index = 0;
      while (true) {
         /* Locate the substring to replace. */
         index = oldinfo.find("_", index);
         if (index == std::string::npos)
            break;
         if ( index > 0 &&   oldinfo.at(index-1) == '\\' )
         {
            ++index;
            continue;
         }

         /* Make the replacement. */
         oldinfo.replace(index, 1, "\\_");

         /* Advance index forward so the next iteration doesn't pick it up as well. */
         index += 2;
      }
      std::cout << "oldinfo: " << oldinfo << std::endl;

      line << "edge from parent node [" << relposition << "] {" << oldinfo <<"} " ;
   }


   return line.str();
}


std::string writeSeeedInfoLatex( SeeedPtr seeed )
{
   std::stringstream line;
   line << "\\node[below = \\belowcaptionskip of s" << seeed->getID() << "] (caps" << seeed->getID() << ") {\\scriptsize " << seeed->getShortCaption() << "}; " << std::endl;

   return line.str();
}


std::string writeSeeedIncludeLatex( SeeedPtr seeed, std::string workfolder )
{
   std::stringstream line;
   line << " (s" << seeed->getID() << ") { \\includegraphics[width=0.15\\textwidth]{" << getSeeedFolderLatex(seeed) << "} }" << std::endl;

   return line.str();
}

SCIP_RETCODE getDetectorCallRoundInfo(SCIP* scip, const char* detectorname, SCIP_Bool transformed, int* maxcallround, int* mincallround, int* freqcallround)
{
   char  setstr[SCIP_MAXSTRLEN];
   if(transformed)
   {
      (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/maxcallround", detectorname);
      SCIP_CALL( SCIPgetIntParam(scip, setstr, maxcallround) );
      (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/mincallround", detectorname);
      SCIP_CALL( SCIPgetIntParam(scip, setstr, mincallround) );
      (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/freqcallround", detectorname);
      SCIP_CALL_ABORT( SCIPgetIntParam(scip, setstr, freqcallround) );
   }
   else
   {
      (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/origmaxcallround", detectorname);
      SCIP_CALL( SCIPgetIntParam(scip, setstr, maxcallround) );
      (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/origmincallround", detectorname);
      SCIP_CALL( SCIPgetIntParam(scip, setstr, mincallround) );
      (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/origfreqcallround", detectorname);
      SCIP_CALL( SCIPgetIntParam(scip, setstr, freqcallround) );
   }

   return SCIP_OKAY;
}

SCIP_Bool cmpSeeedsMaxWhite(SeeedPtr i, SeeedPtr j)
{
   return (i->getMaxWhiteScore() < j->getMaxWhiteScore() );
}


/* method to thin out the vector of given seeeds */
std::vector<SeeedPtr> thinout( std::vector<SeeedPtr> finishedSeeeds, size_t nDecomps, SCIP_Bool addTrivialDecomp )
         {
   std::vector<SeeedPtr> justBest(0);
   for( size_t dec = 0; dec < nDecomps && dec < finishedSeeeds.size(); ++dec)
   {
      justBest.push_back(finishedSeeeds[dec]);
   }

   if(addTrivialDecomp)
   {
      for(size_t dec = 0; dec < finishedSeeeds.size(); ++dec)
      {
         if(finishedSeeeds[dec]->getNMasterconss() == 0 && finishedSeeeds[dec]->getNLinkingvars() == 0 && finishedSeeeds[dec]->getNBlocks() == 1)
         {
            justBest.push_back(finishedSeeeds[dec]);
         }
      }
   }
   return justBest;
      }


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


void removeDigits(char *str, int *nremoved)
{
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

/** method to calculate the greatest common divisor */

int gcd(int a, int b)
{
   return b == 0 ? a : gcd(b, a % b);
}


SCIP_CONS* consGetRelevantRepr(SCIP* scip, SCIP_CONS* cons)
{
   return cons;
}

SCIP_VAR* varGetRelevantRepr(SCIP* scip, SCIP_VAR* var)
{
   return SCIPvarGetProbvar(var);
}

SCIP_Bool seeedIsNoDuplicateOfSeeeds(SeeedPtr compseeed, std::vector<SeeedPtr> const & seeeds, bool sort)
{
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


SCIP_Bool seeedIsNoDuplicate(SeeedPtr seeed, std::vector<SeeedPtr> const & currSeeeds, std::vector<SeeedPtr> const & finishedSeeeds, bool sort)
{
   SCIP_Bool bool1 = seeedIsNoDuplicateOfSeeeds(seeed, currSeeeds, sort);
   SCIP_Bool bool2 = seeedIsNoDuplicateOfSeeeds(seeed, finishedSeeeds, sort);
   return ( bool1 && bool2 );
}


/** constructor */
Seeedpool::Seeedpool(
   SCIP*               givenScip, /**< SCIP data structure */
   const char*       conshdlrName,
   SCIP_Bool         _transformed
):scip(givenScip), currSeeeds(0), allrelevantseeeds(0), nTotalSeeeds(0), nVars(SCIPgetNVars(givenScip) ), nConss(SCIPgetNConss(givenScip) ), nDetectors(0), nFinishingDetectors(0),ndecompositions(0), candidatesNBlocks(0), transformed(_transformed), helpvisucounter(0)
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

//   SCIPdebugMessagePrint(scip, " conshdlrname: %s \n ", conshdlrName);
 //  SCIPdebugMessagePrint(scip, " number of detectors: %d \n ", conshdlrdata->ndetectors);

   assert(conshdlrdata->ndetectors > 0);

   /** store priorities of the detectors */
   for( int d = 0; d < conshdlrdata->ndetectors; ++d )
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
   for( int d = 0; d < conshdlrdata->ndetectors; ++d )
   {
      DEC_DETECTOR* detector;

      detector = conshdlrdata->detectors[d];
      assert(detector != NULL);
      if( transformed )
      {
         if( !detector->enabled || detector->propagateSeeed == NULL)
            continue;
      }
      else
      {
         if( !detector->enabledOrig || detector->propagateSeeed == NULL)
            continue;
      }

      scipDetectorToIndex[detector] = nDetectors;
      detectorToScipDetector.push_back(detector);
      ++nDetectors;
   }

   /** set up enabled finishing detectors */
   for( int d = 0; d < conshdlrdata->ndetectors; ++d )
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

   for( int i = 0; i < nConss; ++i )
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

   for( int i = 0; i < nVars; ++i )
   {
      SCIP_VAR* relevantVar;

      if( transformed )
         relevantVar = varGetRelevantRepr(scip, vars[i]);
      else
         relevantVar = vars[i];

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
      int        nCurrVars;

      cons = consToScipCons[i];

      nCurrVars = GCGconsGetNVars(scip, cons);

      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &currVars, nCurrVars) ); /** free in line 321 */
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &currVals, nCurrVars) ); /** free in line 321 */
      SCIP_CALL_ABORT(GCGconsGetVars(scip, cons, currVars, nCurrVars));
      SCIP_CALL_ABORT(GCGconsGetVals(scip, cons, currVals, nCurrVars));

      for( int currVar = 0; currVar < nCurrVars; ++currVar )
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
         valsMap[std::pair<int,int>(i, varIndex)] =  currVals[currVar];
      }
      SCIPfreeBufferArrayNull(scip, &currVars);
      SCIPfreeBufferArrayNull(scip, &currVals);
   }

   /*  init  seeedpool with empty seeed and translated original seeeds*/
   addSeeedToCurr( new Seeed( scip, -1, nDetectors, nConss, nVars) );

   decompositions = NULL;

}//end constructor

Seeedpool::~Seeedpool()
{
   for( size_t i = 0; i < allrelevantseeeds.size(); ++i )
   {
      size_t help = allrelevantseeeds.size() - i - 1;
      if( allrelevantseeeds[help] != NULL && allrelevantseeeds[help]->getID() >= 0 )
         delete allrelevantseeeds[help];
   }

   for ( size_t i = 0; i < consclassescollection.size(); ++i )
   {
      size_t help = consclassescollection.size() - i - 1;
      if( consclassescollection[help] != NULL )
         delete consclassescollection[help];
   }

   for ( size_t i = 0; i < varclassescollection.size(); ++i )
   {
      size_t help = varclassescollection.size() - i - 1;
      if( varclassescollection[help] != NULL )
         delete varclassescollection[help];
   }
}


/** calculates constraint and variables classes, and find block number candidates */
SCIP_RETCODE Seeedpool::calcConsClassifierAndNBlockCandidates(
   SCIP*               givenScip /**< SCIP data structure */
   ){

   SCIP_Bool conssclassnnonzeros;
   SCIP_Bool conssclassscipconstypes;
   SCIP_Bool conssclassconsnamenonumbers;
   SCIP_Bool conssclassconsnamelevenshtein;
   SCIP_Bool varclassscipvartypes;

   if( transformed )
     {
        SCIPgetBoolParam(scip, "detection/conssclassifier/nnonzeros/enabled", &conssclassnnonzeros);
        SCIPgetBoolParam(scip, "detection/conssclassifier/scipconstype/enabled", &conssclassscipconstypes);
        SCIPgetBoolParam(scip, "detection/conssclassifier/consnamenonumbers/enabled", &conssclassconsnamenonumbers);
        SCIPgetBoolParam(scip, "detection/conssclassifier/consnamelevenshtein/enabled", &conssclassconsnamelevenshtein);
        SCIPgetBoolParam(scip, "detection/varclassifier/scipvartype/enabled", &varclassscipvartypes);
     }
     else
     {
        SCIPgetBoolParam(scip, "detection/conssclassifier/nnonzeros/origenabled", &conssclassnnonzeros);
        SCIPgetBoolParam(scip, "detection/conssclassifier/scipconstype/origenabled", &conssclassscipconstypes);
        SCIPgetBoolParam(scip, "detection/conssclassifier/consnamenonumbers/origenabled", &conssclassconsnamenonumbers);
        SCIPgetBoolParam(scip, "detection/conssclassifier/consnamelevenshtein/origenabled", &conssclassconsnamelevenshtein);
        SCIPgetBoolParam(scip, "detection/varclassifier/scipvartype/origenabled", &varclassscipvartypes);
     }

     std::cout << "consclass nonzeros enabled: " <<conssclassnnonzeros << std::endl;

     if( conssclassnnonzeros )
        addConsClassifier( createConsClassifierForNNonzeros() );
     if( conssclassscipconstypes )
        addConsClassifier( createConsClassifierForSCIPConstypes() );
     if( conssclassconsnamenonumbers )
        addConsClassifier( createConsClassifierForConsnamesDigitFreeIdentical() );
     if( conssclassconsnamelevenshtein )
        addConsClassifier( createConsClassifierForConsnamesLevenshteinDistanceConnectivity(1) );

     if ( varclassscipvartypes )
        addVarClassifier( createVarClassifierForSCIPVartypes() );

     reduceConsclasses();
     reduceVarclasses();

     calcCandidatesNBlocks();

     return SCIP_OKAY;
}


 /** finds decompositions  */
  /** access coefficient matrix constraint-wise */
 std::vector<SeeedPtr>    Seeedpool::findSeeeds(
 ){
    /** 1) read parameter, as there are: maxrounds
     *  2) loop rounds
     *  3) every seeed in seeeds
     *  4) every detector not registered yet propagates seeed
     *  5)  */

    bool displaySeeeds = false;
    int verboseLevel;
    std::vector<int> successDetectors;
    std::vector<SeeedPtr> delSeeeds;
    bool duplicate;

    successDetectors = std::vector<int>(nDetectors, 0);
    ndecompositions = 0;

    delSeeeds = std::vector<SeeedPtr>(0);

    verboseLevel = 1;

    for( size_t i = 0; i < currSeeeds.size(); ++i )
    {
       SCIP_CALL_ABORT( prepareSeeed( currSeeeds[i]) );
       currSeeeds[i]->setID(getNewIdForSeeed() );
    }

    /** add translated original seeeds (of unpresolved problem) */
    for( size_t i = 0; i < translatedOrigSeeeds.size(); ++i )
    {
       SCIP_CALL_ABORT( prepareSeeed( translatedOrigSeeeds[i]) );
       translatedOrigSeeeds[i]->setID(getNewIdForSeeed() );
       if( seeedIsNoDuplicateOfSeeeds(translatedOrigSeeeds[i], currSeeeds, true) )
          currSeeeds.push_back(translatedOrigSeeeds[i] );
    }

    for( int round = 0; round < maxndetectionrounds; ++round )
    {
       std::cout << "currently in detection round " << round << std::endl;
       std::vector<SeeedPtr> nextSeeeds = std::vector<SeeedPtr>(0);
       std::vector<SeeedPtr> currSeeedsToDelete = std::vector<SeeedPtr>(0);

       #pragma omp parallel for schedule(static,1)
       for( size_t s = 0; s < currSeeeds.size(); ++s )
       {
          SeeedPtr seeedPtr;
          seeedPtr = currSeeeds[s];

          #pragma omp critical (ostream)
          {
             if( displaySeeeds || verboseLevel >= 1 )
             {
                std::cout << "Start to propagate seeed " << seeedPtr->getID() << " in round " << round << ":" << std::endl;
                if( displaySeeeds )
                   seeedPtr->displaySeeed();
             }
          }

//          SCIPdebugMessagePrint(scip, "number of detectors: %d \n", nDetectors);

          /** the current seeed is handled by all detectors */
          for( int d = 0; d < nDetectors; ++d )
          {
             SEEED_PROPAGATION_DATA* seeedPropData;
             seeedPropData = new SEEED_PROPAGATION_DATA();
             seeedPropData->seeedpool = this;
             seeedPropData->nNewSeeeds = 0;
             DEC_DETECTOR* detector;
             std::vector<SeeedPtr>::const_iterator newSIter;
             std::vector<SeeedPtr>::const_iterator newSIterEnd;
             int maxcallround;
             int mincallround;
             int freqcallround;
             const char* detectorname;
             SCIP_CLOCK* detectorclock;

             detector = detectorToScipDetector[d];
             detectorname = DECdetectorGetName(detector);
             SCIP_RESULT result = SCIP_DIDNOTFIND;

             /** if the seeed is also propagated by the detector go on with the next detector */
             if( seeedPtr->isPropagatedBy(detector) && !detector->usefulRecall )
                continue;

             /** check if detector is callable in current detection round */
             SCIP_CALL_ABORT( getDetectorCallRoundInfo( scip, detectorname, transformed, &maxcallround, &mincallround, &freqcallround) );

             if( maxcallround < round || mincallround > round || ( (round - mincallround) % freqcallround != 0 ) )
                continue;

             #pragma omp critical (seeedcount)
             seeedPropData->seeedToPropagate = new gcg::Seeed(seeedPtr, this );


             /** new seeeds are created by the current detector */
             #pragma omp critical (clock)
             SCIPcreateClock(scip, &detectorclock);
             SCIP_CALL_ABORT( SCIPstartClock(scip, detectorclock) );

             if( verboseLevel >= 1 )
             {
                #pragma omp critical (scipinfo)
                SCIPinfoMessage(scip, NULL, "detector %s started to propagate the %d-th seeed (ID %d ) in round %d \n", DECdetectorGetName(detectorToScipDetector[d]), s+1, seeedPtr->getID(), round);
             }

             SCIP_CALL_ABORT(detectorToScipDetector[d]->propagateSeeed(scip, detectorToScipDetector[d], seeedPropData, &result) );

             for( int j = 0; j < seeedPropData->nNewSeeeds; ++j )
             {
               #pragma omp critical (seeedcount)
                seeedPropData->newSeeeds[j]->setID(getNewIdForSeeed());
                prepareSeeed( seeedPropData->newSeeeds[j] );
                assert(seeedPropData->newSeeeds[j]->checkConsistency(this) );
                seeedPropData->newSeeeds[j]->addDecChangesFromAncestor(seeedPtr);
             }

             SCIP_CALL_ABORT( SCIPstopClock(scip, detectorclock) );

             #pragma omp critical (clockcount)
             detectorToScipDetector[d]->dectime += SCIPgetClockTime(scip, detectorclock);

             #pragma omp critical (clock)
             SCIPfreeClock(scip, &detectorclock);

             if(seeedPropData->nNewSeeeds != 0 && ( displaySeeeds ) )
             {
                #pragma omp critical (ostream)
                std::cout << "detector " << DECdetectorGetName(detectorToScipDetector[d] ) << " found " << seeedPropData->nNewSeeeds << " new seeed(s): ";
                #pragma omp critical (ostream)
                std::cout << seeedPropData->newSeeeds[0]->getID();
                for( int j = 1; j < seeedPropData->nNewSeeeds; ++j )
                {
                   #pragma omp critical (ostream)
                   std::cout << ", " << seeedPropData->newSeeeds[j]->getID();
                }
                #pragma omp critical (ostream)
                std::cout << "\n";

                if( displaySeeeds )
                {
                   for( int j = 0; j < seeedPropData->nNewSeeeds; ++j )
                   {
                      #pragma omp critical (ostream)
                      seeedPropData->newSeeeds[j]->displaySeeed();
                   }
                }
             }
             else
                if( displaySeeeds )
                {
                   #pragma omp critical (ostream)
                   std::cout << "detector " << DECdetectorGetName(detectorToScipDetector[d] ) << " found 0 new seeeds" << std::endl;
                }

             /** if the new seeeds are no duplicate they're added to the currSeeeds */
             for( int seeed = 0; seeed < seeedPropData->nNewSeeeds; ++seeed )
             {
                SCIP_Bool noduplicate;
                #pragma omp critical (seeedptrstore)
                noduplicate = seeedIsNoDuplicate(seeedPropData->newSeeeds[seeed], nextSeeeds, finishedSeeeds, false);
                if( !seeedPropData->newSeeeds[seeed]->isTrivial() && noduplicate )
                {
                   seeedPropData->newSeeeds[seeed]->calcOpenconss();
                   seeedPropData->newSeeeds[seeed]->calcOpenvars();
                   if(seeedPropData->newSeeeds[seeed]->getNOpenconss() == 0 && seeedPropData->newSeeeds[seeed]->getNOpenvars() == 0)
                   {
                      if(verboseLevel > 2)
                      {
                         #pragma omp critical (ostream)
                         {
                            std::cout << "seeed " << seeedPropData->newSeeeds[seeed]->getID() << " is addded to finished seeeds!" << std::endl;
                            seeedPropData->newSeeeds[seeed]->showScatterPlot(this);
                         }
                      }
                      #pragma omp critical (seeedptrstore)
                      {
                         assert(seeedPropData->newSeeeds[seeed]->getID() >= 0);
                         finishedSeeeds.push_back(seeedPropData->newSeeeds[seeed]);
                         allrelevantseeeds.push_back(seeedPropData->newSeeeds[seeed]);
                      }
                   }
                   else
                   {
                      if(verboseLevel > 2)
                      {
                         #pragma omp critical (ostream)
                         {
                            std::cout << "seeed " << seeedPropData->newSeeeds[seeed]->getID() << " is addded to next round seeeds!" << std::endl;
                            seeedPropData->newSeeeds[seeed]->showScatterPlot(this);
                         }
                      }
                      #pragma omp critical (seeedptrstore)
                      {
                         nextSeeeds.push_back(seeedPropData->newSeeeds[seeed]);
                         allrelevantseeeds.push_back(seeedPropData->newSeeeds[seeed]);
                      }
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
             delete seeedPropData->seeedToPropagate;
             seeedPropData->newSeeeds = NULL;
             seeedPropData->nNewSeeeds = 0;
             delete seeedPropData;
          } // end for detectors

          for( int d = 0; d < nFinishingDetectors; ++d )
          {
             DEC_DETECTOR* detector = detectorToFinishingScipDetector[d];
             SCIP_RESULT result = SCIP_DIDNOTFIND;
             SEEED_PROPAGATION_DATA* seeedPropData;
             seeedPropData = new SEEED_PROPAGATION_DATA();
             seeedPropData->seeedpool = this;
             seeedPropData->nNewSeeeds = 0;
             #pragma omp critical (seeedcount)
             seeedPropData->seeedToPropagate = new gcg::Seeed(seeedPtr, this );

             if(verboseLevel > 2 )
             #pragma omp critical (ostream)
             {
                std::cout << "check if finisher of detector " << DECdetectorGetName(detectorToFinishingScipDetector[d] ) << " is enabled " << std::endl;
             }

             /** if the finishing of the detector is not enabled go on with the next detector */
             if( !detector->enabledFinishing )
                continue;

             if(verboseLevel > 2 )
             #pragma omp critical (ostream)
             {
                std::cout << "call finisher for detector " << DECdetectorGetName(detectorToFinishingScipDetector[d] ) << std::endl;
             }
             SCIP_CALL_ABORT(detectorToFinishingScipDetector[d]->finishSeeed(scip, detectorToFinishingScipDetector[d], seeedPropData, &result) );

             for( int finished = 0; finished < seeedPropData->nNewSeeeds; ++finished )
             {
                SeeedPtr seeed = seeedPropData->newSeeeds[finished];
                #pragma omp critical (seeedcount)
                seeed->setID(getNewIdForSeeed());
                seeed->sort();
                seeed->calcHashvalue();
                seeed->addDecChangesFromAncestor(seeedPtr);
                seeed->setFinishedByFinisher(true);
                #pragma omp critical (seeedptrstore)
                {
                   if( seeedIsNoDuplicateOfSeeeds(seeed, finishedSeeeds, false) )
                   {
                      assert(seeed->getID() >= 0);
                      finishedSeeeds.push_back(seeed);
                      allrelevantseeeds.push_back(seeed);
                   }
                   else
                   {
                      bool isIdentical = false;
                      for ( size_t h = 0; h < finishedSeeeds.size(); ++h )
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
             }
             SCIPfreeMemoryArrayNull(scip, &seeedPropData->newSeeeds);
             delete seeedPropData->seeedToPropagate;
             seeedPropData->newSeeeds = NULL;
             seeedPropData->nNewSeeeds = 0;
             delete seeedPropData;
          }
       }// end for currseeeds

       for(size_t s = 0; s < currSeeedsToDelete.size(); ++s )
       {
          delete currSeeedsToDelete[s];
          currSeeedsToDelete[s] = NULL;
       }

       currSeeeds = nextSeeeds;
    } // end for rounds

    /** complete the currseeeds with finishing detectors and add them to finished seeeds */
    #pragma omp parallel for schedule(static,1)
    for( size_t i = 0; i < currSeeeds.size(); ++i )
    {
       SeeedPtr seeedPtr = currSeeeds[i];
       for( int d = 0; d < nFinishingDetectors; ++d )
       {
          DEC_DETECTOR* detector = detectorToFinishingScipDetector[d];
          SCIP_RESULT result = SCIP_DIDNOTFIND;
          SEEED_PROPAGATION_DATA* seeedPropData;
          seeedPropData = new SEEED_PROPAGATION_DATA();
          seeedPropData->seeedpool = this;
          seeedPropData->nNewSeeeds = 0;

          #pragma omp critical (seeedcount)
          seeedPropData->seeedToPropagate = new gcg::Seeed(seeedPtr, this );

          if( verboseLevel > 2)
             std::cout << "check if finisher of detector " << DECdetectorGetName(detectorToScipDetector[d] ) << " is enabled " << std::endl;

          /** if the finishing of the detector is not enabled go on with the next detector */
          if( !detector->enabledFinishing )
             continue;

          std::cout << "call finisher for detector " << DECdetectorGetName(detectorToFinishingScipDetector[d] ) << std::endl;

          SCIP_CALL_ABORT(detectorToFinishingScipDetector[d]->finishSeeed(scip, detectorToFinishingScipDetector[d],seeedPropData, &result) );

          for( int finished = 0; finished < seeedPropData->nNewSeeeds; ++finished )
          {
             SeeedPtr seeed = seeedPropData->newSeeeds[finished];
             #pragma omp critical (seeedcount)
             seeed->setID( getNewIdForSeeed() );

             seeed->calcHashvalue();
             seeed->addDecChangesFromAncestor(seeedPtr);
             seeed->setFinishedByFinisher(true);

             if( seeedIsNoDuplicateOfSeeeds(seeed, finishedSeeeds, false) )
             {
                if( verboseLevel > 2 )
                {
                   std::cout << "seeed " << seeed->getID() << " is finished from next round seeeds!" << std::endl;
                   seeed->showScatterPlot(this);
                }
                #pragma omp critical (seeedptrstore)
                {
                   assert(seeed->getID() >= 0);
                   finishedSeeeds.push_back(seeed);
                   allrelevantseeeds.push_back(seeed);
                }
             }

             SCIPfreeMemoryArrayNull(scip, &seeedPropData->newSeeeds);
             seeedPropData->newSeeeds = NULL;
             seeedPropData->nNewSeeeds = 0;
          }

          delete seeedPropData->seeedToPropagate;
          delete seeedPropData;
       }
    }// end for finishing curr seeeds

    std::cout << (int) finishedSeeeds.size() << " finished seeeds are found." << std::endl;

    if( displaySeeeds )
    {
       for( size_t i = 0; i < finishedSeeeds.size(); ++i )
       {
          std::cout << i+1 << "th finished seeed: " << std::endl;
          finishedSeeeds[i]->displaySeeed();
       }
    }

    /** count the successful refinement calls for each detector */

    for(size_t i = 0; i < finishedSeeeds.size(); ++i)
    {
       assert(finishedSeeeds[i]->checkConsistency(this) );
       assert(finishedSeeeds[i]->getNOpenconss() == 0);
       assert(finishedSeeeds[i]->getNOpenvars() == 0);

       SCIP_CALL_ABORT(finishedSeeeds[i]->buildDecChainString() );

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
       std::cout << "Detector " << std::setw(25) << std::setiosflags(std::ios::left) << DECdetectorGetName(detectorToScipDetector[i] ) << " \t worked on \t " << successDetectors[i] << " of " << finishedSeeeds.size() << "\t and took a total time of \t" << detectorToScipDetector[i]->dectime  << std::endl;
    }

    if( (int) finishedSeeeds.size() != 0)
    {
       SCIP_Real minscore = finishedSeeeds[0]->evaluate(this);
       //            SeeedPtr bestSeeed = finishedSeeeds[0];
       for( size_t i = 1; i < finishedSeeeds.size(); ++i )
       {
          SCIP_Real score = finishedSeeeds[i]->evaluate(this);
          if (score < minscore)
          {
             minscore = score;
             //                  bestSeeed = finishedSeeeds[i];
          }
       }
       //            bestSeeed->showScatterPlot(this);
    }


    /** delete the seeeds */
    for( size_t c = 0; c < currSeeeds.size(); ++c )
    {
       duplicate = false;
       SCIP_CALL_ABORT(currSeeeds[c]->buildDecChainString() );
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

    /* postpone deleting to destructor */
    //         for( size_t d =  delSeeeds.size(); d > 0; d--)
    //         {
    //            delete delSeeeds[d-1];
    //            delSeeeds[d-1] = NULL;
    //         }
    //
    //         delSeeeds.clear();



    sortAllRelevantSeeeds();

    return finishedSeeeds;

 }


 /** method to complete a set of incomplete seeeds with the help of all included detectors that implement a finishing method */
/*
 * @return set of completed decomposition
 * */

 std::vector<SeeedPtr>    Seeedpool::finishIncompleteSeeeds(
    std::vector<SeeedPtr> incompleteseeeds
  )
{
    std::vector<SeeedPtr> finisheds(0, NULL);
    int verboseLevel = 1;

    #pragma omp parallel for schedule(static,1)
    for( size_t i = 0; i < incompleteseeeds.size(); ++i )
    {
       SeeedPtr seeedPtr = incompleteseeeds[i];

       for( int d = 0; d < nFinishingDetectors; ++d )
       {
          DEC_DETECTOR* detector = detectorToFinishingScipDetector[d];
          SCIP_RESULT result = SCIP_DIDNOTFIND;
          SEEED_PROPAGATION_DATA* seeedPropData;
          seeedPropData = new SEEED_PROPAGATION_DATA();
          seeedPropData->seeedpool = this;
          seeedPropData->nNewSeeeds = 0;

          #pragma omp critical (seeedcount)
          seeedPropData->seeedToPropagate = new gcg::Seeed(seeedPtr, this );

          if( verboseLevel > 2)
             std::cout << "check if finisher of detector " << DECdetectorGetName(detectorToScipDetector[d] ) << " is enabled " << std::endl;

          /** if the finishing of the detector is not enabled go on with the next detector */
          if( !detector->enabledFinishing )
             continue;

          std::cout << "call finisher for detector " << DECdetectorGetName(detectorToFinishingScipDetector[d] ) << std::endl;

          SCIP_CALL_ABORT(detectorToFinishingScipDetector[d]->finishSeeed(scip, detectorToFinishingScipDetector[d],seeedPropData, &result) );

          for( int finished = 0; finished < seeedPropData->nNewSeeeds; ++finished )
          {
             SeeedPtr seeed = seeedPropData->newSeeeds[finished];
#pragma omp critical (seeedcount)
             seeed->setID( getNewIdForSeeed() );

             seeed->calcHashvalue();
             seeed->addDecChangesFromAncestor(seeedPtr);
             seeed->setFinishedByFinisher(true);

             if( seeedIsNoDuplicateOfSeeeds(seeed, finishedSeeeds, false) )
             {
                if( verboseLevel > 2 )
                {
                   std::cout << "seeed " << seeed->getID() << " is finished from next round seeeds!" << std::endl;
                   seeed->showScatterPlot(this);
                }
#pragma omp critical (seeedptrstore)
                {
                   assert(seeed->getID() >= 0);
                   finisheds.push_back(seeed);
       //            allrelevantseeeds.push_back(seeed);
                }
             }

             SCIPfreeMemoryArrayNull(scip, &seeedPropData->newSeeeds);
             seeedPropData->newSeeeds = NULL;
             seeedPropData->nNewSeeeds = 0;
          }

          delete seeedPropData->seeedToPropagate;
          delete seeedPropData;
       }
    }// end for finishing curr seeeds
    return finisheds;
}

 void    Seeedpool::findDecompositions(
 ){

    /**
     * finds seeeds and translates them to decompositions
     */

    std::vector<int> successDetectors;
    std::vector<SeeedPtr> delSeeeds;
    SCIP_Bool usemaxwhitescore;
    SCIP_Bool dothinout;


    size_t nDecomps = 6;

//    SCIPdebugMessage("Starting find decomps of seeedpool\n");

    SCIP_Bool addTrivialDecomp = FALSE;

    successDetectors = std::vector<int>(nDetectors, 0);
    ndecompositions = 0;
    delSeeeds = std::vector<SeeedPtr>(0);
    usemaxwhitescore = TRUE;
    dothinout = FALSE;

    finishedSeeeds = findSeeeds();

    /* sort the seeeds according to maximum white measure */

    std::sort (finishedSeeeds.begin(), finishedSeeeds.end(), cmpSeeedsMaxWhite);

    /** hack to just use max white seeed */
    if( usemaxwhitescore && dothinout )
       finishedSeeeds = thinout( finishedSeeeds, nDecomps, addTrivialDecomp );

    /** fill out the decompositions */

    SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &decompositions, (int) finishedSeeeds.size())); /** free in decomp.c:470 */
    for( size_t i = 0; i < finishedSeeeds.size(); ++i )
    {
       SeeedPtr seeed = finishedSeeeds[i];

       SCIP_CALL_ABORT( createDecompFromSeeed(seeed, &decompositions[i]) );
    }

    //SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );

    /** delete the seeeds */

    //         for(size_t f = 0; f < finishedSeeeds.size(); ++f)
    //         {
    //            duplicate = false;
    //            for(size_t d = 0; d < delSeeeds.size(); ++d)
    //            {
    //               if(finishedSeeeds[f]==delSeeeds[d])
    //               {
    //                  duplicate=true;
    //                  break;
    //               }
    //            }
    //            if(!duplicate)
    //            {
    //               delSeeeds.push_back(finishedSeeeds[f]);
    //            }
    //         }
    //
    //
    //         for( size_t d =  delSeeeds.size(); d > 0; d--)
    //         {
    //            delete delSeeeds[d-1];
    //            delSeeeds[d-1] = NULL;
    //         }

    delSeeeds.clear();

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


 SCIP_RETCODE Seeedpool::prepareSeeed( SeeedPtr seeed)
 {
   seeed->considerImplicits(this);
   seeed->sort();
   seeed->calcHashvalue();
   seeed->evaluate(this);

   return SCIP_OKAY;
 }


void Seeedpool::freeCurrSeeeds()
{
   for( size_t i = 0; i < currSeeeds.size(); ++i )
   {
      if ( currSeeeds[i] != NULL )
      {
         currSeeeds[i]->checkConsistency(this);
         delete currSeeeds[i];
         currSeeeds[i] = NULL;
      }
   }
   return;
}


void Seeedpool::addSeeedToCurr(SeeedPtr seeed){

   currSeeeds.push_back(seeed);
   allrelevantseeeds.push_back(seeed);
   return;
}

void Seeedpool::addSeeedToFinished(SeeedPtr seeed){

   finishedSeeeds.push_back(seeed);
   allrelevantseeeds.push_back(seeed);
   return;
}

void Seeedpool::sortAllRelevantSeeeds(){

   int maxid  = 0;
   std::vector<SeeedPtr> tmpAllRelevantSeeeds(0);

   for ( size_t i = 0; i < allrelevantseeeds.size(); ++i )
   {
      if( allrelevantseeeds[i]->getID() > maxid )
         maxid = allrelevantseeeds[i]->getID();
   }

   tmpAllRelevantSeeeds = std::vector<SeeedPtr>(maxid+1, NULL );

   for ( size_t i = 0; i < allrelevantseeeds.size(); ++i )
   {
      if ( allrelevantseeeds[i]->getID() < 0  )
         continue;
      tmpAllRelevantSeeeds[allrelevantseeeds[i]->getID()] = allrelevantseeeds[i];
   }

   allrelevantseeeds = tmpAllRelevantSeeeds;

}


/** translates seeeds and classifiers if the index structure of the problem has changed, e.g. due to presolving */
void Seeedpool::translateSeeedData( Seeedpool* origpool, std::vector<Seeed*> origseeeds, std::vector<Seeed*>& newseeeds,
   std::vector<ConsClassifier*> otherconsclassifiers, std::vector<ConsClassifier*>& newconsclassifiers,
   std::vector<VarClassifier*> othervarclassifiers, std::vector<VarClassifier*>& newvarclassifiers )
{
   assert( newseeeds.empty() );
   assert( newconsclassifiers.empty() );
   assert( newvarclassifiers.empty() );

   std::vector<int> rowothertothis;
   std::vector<int> rowthistoother;
   std::vector<int> colothertothis;
   std::vector<int> colthistoother;
   std::vector<int> missingrowinthis;

   SCIPdebugMessagePrint(this->scip, " started translate seeed method \n" );

   calcTranslationMapping( origpool, rowothertothis, rowthistoother, colothertothis, colthistoother, missingrowinthis );

   SCIPdebugMessagePrint(this->scip, " calculated translation; number of missing constraints: %d; number of other seeeds: %d \n", missingrowinthis.size(), origseeeds.size() );

   newseeeds = getTranslatedSeeeds( origseeeds, rowothertothis, rowthistoother, colothertothis, colthistoother );
   newconsclassifiers = getTranslatedConsClassifiers( otherconsclassifiers, rowothertothis, rowthistoother );
   newvarclassifiers = getTranslatedVarClassifiers( othervarclassifiers, colothertothis, colthistoother );
}

/** translates seeeds if the index structure of the problem has changed, e.g. due to presolving */
void Seeedpool::translateSeeeds( Seeedpool* origpool, std::vector<Seeed*> origseeeds, std::vector<Seeed*>& newseeeds )
{
   assert( newseeeds.empty() );

   std::vector<int> rowothertothis(0);
   std::vector<int> rowthistoother(0);
   std::vector<int> colothertothis(0);
   std::vector<int> colthistoother(0);
   std::vector<int> missingrowinthis(0);

   SCIPdebugMessagePrint(this->scip, " started translate seeed method \n" );

   calcTranslationMapping( origpool, rowothertothis, rowthistoother, colothertothis, colthistoother, missingrowinthis );

   SCIPdebugMessagePrint(this->scip, " calculated translation; number of missing constraints: %d; number of other seeeds: %d \n", missingrowinthis.size(), origseeeds.size() );

   newseeeds = getTranslatedSeeeds( origseeeds, rowothertothis, rowthistoother, colothertothis, colthistoother );
}


/** calculates necessary data for translating seeeds and classifiers */
void Seeedpool::calcTranslationMapping( Seeedpool* origpool, std::vector<int>& rowothertothis, std::vector<int>& rowthistoother,
   std::vector<int>& colothertothis, std::vector<int>& colthistoother, std::vector<int>& missingrowinthis )
{
   int nrowsother  = origpool->nConss;
   int nrowsthis  = nConss;
   int ncolsother  = origpool->nVars;
   int ncolsthis  = nVars;

   rowothertothis.assign(nrowsother, -1);
   rowthistoother.assign(nrowsthis, -1);
   colothertothis.assign(ncolsother, -1);
   colthistoother.assign(ncolsthis, -1);

   missingrowinthis.clear();


   /* identify new and deleted rows and cols; and identify bijection between maintained variables */
   for( int i = 0; i < nrowsother; ++i  )
   {
      SCIP_CONS* otherrow = origpool->getConsForIndex(i);
      assert(otherrow != NULL);
      SCIP_Bool foundmaintained = FALSE;
      for( int j = 0; j < nrowsthis; ++j  )
      {
         SCIP_CONS* thisrow = this->getConsForIndex(j);
         assert(SCIPconsIsTransformed(thisrow) );
         char buffer[SCIP_MAXSTRLEN];
         assert(this->scip != NULL);
         strcpy(buffer, SCIPconsGetName(thisrow) + 2);
         assert(this->scip != NULL);
         if( strcmp(SCIPconsGetName(otherrow), SCIPconsGetName(thisrow) ) == 0 )
         {
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
      for( int j = 0; j < ncolsthis; ++j  )
      {
         if( othervar == this->getVarForIndex(j) )
         {
            colothertothis[i] = j;
            colthistoother[j] = i;
            break;
         }
      }
   }
}

/** returns translated seeeds derived from given mapping data */
std::vector<Seeed*> Seeedpool::getTranslatedSeeeds( std::vector<Seeed*>& origseeeds, std::vector<int>& rowothertothis, std::vector<int>& rowthistoother,
   std::vector<int>& colothertothis, std::vector<int>& colthistoother )
{
   std::vector<Seeed*> newseeeds(0);

   for( size_t s = 0; s < origseeeds.size(); ++s )
   {
      SeeedPtr otherseeed;
      SeeedPtr newseeed;

      otherseeed = origseeeds[s];

 //     SCIPdebugMessagePrint(this->scip, " otherseeed seeed %d has %d many blocks \n", otherseeed->getID(), otherseeed->getNBlocks() );

      /** ignore seeeds with one block or no block, they are supposed to be found anyway */
      if( otherseeed->getNBlocks() == 1 || otherseeed->getNBlocks() == 0  )
         continue;

      SCIPdebugMessagePrint(this->scip, " transform seeed %d \n", otherseeed->getID() );

      newseeed = new Seeed(scip, this->getNewIdForSeeed(), this->getNDetectors(), this->getNConss(), this->getNVars() );

      /** prepare new seeed */
      newseeed->calcOpenconss();
      newseeed->calcOpenvars();
      newseeed->setOpenVarsAndConssCalculated(true);

      newseeed->setNBlocks(otherseeed->getNBlocks() );

      newseeed->usergiven = otherseeed->usergiven;

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

//         for ( int j = 0; j < otherseeed->getNVarsForBlock(b); j++ )
//         {
//            int thisvar = colothertothis[otherseeed->getVarsForBlock(b)[j] ];
//            if( thisvar != -1 )
//            {
//               newseeed->setVarToBlock(thisvar, b);
//               newseeed->deleteOpenvar(thisvar);
//            }
//         }
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
      newseeed->detectorClockTimes = otherseeed->detectorClockTimes;
      newseeed->pctVarsFromFree = otherseeed->pctVarsFromFree;
      newseeed->pctVarsToBlock = otherseeed->pctVarsToBlock;
      newseeed->pctVarsToBorder = otherseeed->pctVarsToBorder;
      newseeed->pctConssToBorder = otherseeed->pctConssToBorder;
      newseeed->pctConssFromFree = otherseeed->pctConssFromFree;
      newseeed->pctConssToBlock = otherseeed->pctConssToBlock;
      newseeed->nNewBlocks = otherseeed->nNewBlocks;
      newseeed->detectorchainstring = otherseeed->detectorchainstring;





      newseeed->stemsFromUnpresolved = true;
      newseeed->isFinishedByFinisherUnpresolved  = otherseeed->isFinishedByFinisher;

      if ( otherseeed->isFinishedByFinisher )
         newseeed->finishedUnpresolvedBy = otherseeed->detectorChain[otherseeed->detectorChain.size()-1];


      newseeed->setFinishedByFinisher(otherseeed->isFinishedByFinisher);
      newseeed->sort();
      newseeed->considerImplicits(this);
      newseeed->deleteEmptyBlocks();
      newseeed->evaluate(this);

      //assert(newseeed->checkConsistency(this) );



      /*std::cout << "unpresolved seeed " << std::endl;
         otherseeed->showScatterPlot(origpool);
         std::cout << "has become " << std::endl;
         newseeed->showScatterPlot(this);
       */

      if( newseeed->checkConsistency(this) )
         newseeeds.push_back(newseeed);
      else {
         delete newseeed;
         newseeed = NULL;
      }
   }

   return newseeeds;
}


/** returns translated ConsClassifiers derived from given mapping data */
std::vector<ConsClassifier*> Seeedpool::getTranslatedConsClassifiers( std::vector<ConsClassifier*>& otherclassifiers,
   std::vector<int>& rowothertothis, std::vector<int>& rowthistoother )
{
   std::vector<ConsClassifier*> newclassifiers(0);

   for ( size_t i = 0; i < otherclassifiers.size(); ++i )
   {
      ConsClassifier* oldclassifier = otherclassifiers[i];
      ConsClassifier* newclassifier;
      std::stringstream newname;

      newname << oldclassifier->getName() << "-origp";
      newclassifier = new ConsClassifier(scip, newname.str().c_str(), oldclassifier->getNClasses(), (int) rowthistoother.size());
      int bufferclassindex = -1;

      /** copy class information */
      for ( int j = 0; j < oldclassifier->getNClasses(); ++j )
      {
         newclassifier->setClassName( j, oldclassifier->getClassName(j) );
         newclassifier->setClassDescription( j, oldclassifier->getClassDescription(j) );
         newclassifier->setClassDecompInfo( j, oldclassifier->getClassDecompInfo(j) );
      }

      /** assign new conss to classes */
      for ( int c = 0; c < (int) rowthistoother.size(); ++c )
      {
         if ( rowthistoother[c] != -1 )
         {
            newclassifier->assignConsToClass( c, oldclassifier->getClassOfCons( rowthistoother[c] ) );
         }
         else
         {
            if ( bufferclassindex == -1)
            {
               bufferclassindex = newclassifier->addClass( "buffer", "This class contains constraints which are new in the presolved problem.", BOTH );
            }
            newclassifier->assignConsToClass( c, bufferclassindex );
         }
      }

      /** remove empty classes */
      newclassifier->removeEmptyClasses();

      newclassifiers.push_back(newclassifier);
   }

   return newclassifiers;
}


/** returns translated VarClassifiers derived from given mapping data */
std::vector<VarClassifier*> Seeedpool::getTranslatedVarClassifiers( std::vector<VarClassifier*>& otherclassifiers,
   std::vector<int>& colothertothis, std::vector<int>& colthistoother )
{
   std::vector<VarClassifier*> newclassifiers(0);

   for ( size_t i = 0; i < otherclassifiers.size(); ++i )
   {
      VarClassifier* oldclassifier = otherclassifiers[i];
      VarClassifier* newclassifier;
      std::stringstream newname;

      newname << oldclassifier->getName() << "-origp";
      newclassifier = new VarClassifier( scip, newname.str().c_str(), oldclassifier->getNClasses(), (int) colthistoother.size());
      int bufferclassindex = -1;

      /** copy class information */
      for ( int j = 0; j < oldclassifier->getNClasses(); ++j )
      {
         newclassifier->setClassName( j, oldclassifier->getClassName(j) );
         newclassifier->setClassDescription( j, oldclassifier->getClassDescription(j) );
         newclassifier->setClassDecompInfo( j, oldclassifier->getClassDecompInfo(j) );
      }

      /** assign new vars to classes */
      for ( int c = 0; c < (int) colthistoother.size(); ++c )
      {
         if ( colthistoother[c] != -1 )
         {
            newclassifier->assignVarToClass( c, oldclassifier->getClassOfVar( colthistoother[c] ) );
         }
         else
         {
            if ( bufferclassindex == -1)
            {
               bufferclassindex = newclassifier->addClass( "buffer", "This class contains variables which are new in the presolved problem.", ALL );
            }
            newclassifier->assignVarToClass( c, bufferclassindex );
         }
      }

      /** remove empty classes */
      newclassifier->removeEmptyClasses();

      newclassifiers.push_back(newclassifier);
   }

   return newclassifiers;
}

void Seeedpool::populate(std::vector<SeeedPtr> seeeds)
{
   translatedOrigSeeeds = seeeds;
}

const  int * Seeedpool::getVarsForCons(int cons)
{
   return &varsForConss[cons][0];
}

const  SCIP_Real * Seeedpool::getValsForCons(int cons)
{
   return &valsForConss[cons][0];
}

/** access coefficient matrix variable-wise */
const  int * Seeedpool::getConssForVar(int var)
{
   return &conssForVars[var][0];
}

int Seeedpool::getNVarsForCons(int cons)
{
   return varsForConss[cons].size();
}

int Seeedpool::getNConssForVar(int var)
{
   return conssForVars[var].size();
}

SCIP_VAR* Seeedpool::getVarForIndex(int varIndex)
{
   return varToScipVar[varIndex];
}

SCIP_CONS* Seeedpool::getConsForIndex(int consIndex)
{
   return consToScipCons[consIndex];
}

DEC_DETECTOR* Seeedpool::getDetectorForIndex(int detectorIndex)
{
   return detectorToScipDetector[detectorIndex];
}

DEC_DETECTOR* Seeedpool::getFinishingDetectorForIndex(int detectorIndex)
{
   return detectorToFinishingScipDetector[detectorIndex];
}


SCIP_Real Seeedpool::getVal(int row, int col)
{
   std::tr1::unordered_map< std::pair<int, int>, SCIP_Real, pair_hash>::const_iterator iter =  valsMap.find(std::pair<int, int>(row, col) ) ;

   if ( iter == valsMap.end()  )
      return 0;

   return iter->second;
}

int Seeedpool::getIndexForVar(SCIP_VAR* var)
{
   return scipVarToIndex[var];
}

int Seeedpool::getIndexForCons(SCIP_CONS* cons)
{
   return scipConsToIndex[cons];
}

int Seeedpool::getIndexForDetector(DEC_DETECTOR* detector)
{
   return scipDetectorToIndex[detector];
}

int Seeedpool::getIndexForFinishingDetector(DEC_DETECTOR* detector)
{
   return scipFinishingDetectorToIndex[detector];
}


int Seeedpool::getNewIdForSeeed()
{
   return SCIPconshdlrDecompGetNextSeeedID(scip) ;
}

DEC_DECOMP** Seeedpool::getDecompositions()
{
   return decompositions;
}

int Seeedpool::getNDecompositions()
{
   return ndecompositions;
}

int Seeedpool::getNDetectors()
{
   return nDetectors;
}

int Seeedpool::getNFinishingDetectors()
{
   return nFinishingDetectors;
}


int Seeedpool::getNVars()
{
   return nVars;
}

int Seeedpool::getNConss()
{
   return nConss;
}

std::vector<int> Seeedpool::getSortedCandidatesNBlocks()
{
   std::vector<int> toreturn(0);
   SCIP_Bool output = FALSE;

   /** first: sort the current candidates */
   std::sort(candidatesNBlocks.begin(), candidatesNBlocks.end(), sort_decr() );

   if( output )
   {
      std::cout << "nCandidates: " << candidatesNBlocks.size() << std::endl;
      for( size_t i = 0; i < candidatesNBlocks.size(); ++i )
         std::cout << "nblockcandides: " << candidatesNBlocks[i].first << " ; " << candidatesNBlocks[i].second << " times prop " << std::endl;
   }

   for( size_t i = 0; i < candidatesNBlocks.size(); ++i )
      toreturn.push_back(candidatesNBlocks[i].first);

   return toreturn;
}

void Seeedpool::addCandidatesNBlocks(
   int                 candidate            /**< candidate for block size */
)
{
   if( candidate > 1 )
   {
      bool alreadyIn = false;
      for(size_t i = 0; i < candidatesNBlocks.size(); ++i )
      {
         if(candidatesNBlocks[i].first == candidate)
         {
            alreadyIn = true;
            ++candidatesNBlocks[i].second;
            break;
         }
      }
      if(!alreadyIn)
      {
         std::cout << "added block number candidate : " << candidate << std::endl;
         candidatesNBlocks.push_back(std::pair<int,int>(candidate, 1) );
      }
   }
   return;
}

void Seeedpool::calcCandidatesNBlocks()
 {
   /**
    * for every subset of constraint classes calculate gcd (greatest common divisors) of the corresponding number of occurrences
    */

   int maximumnclasses = 18; /* if  distribution of classes exceed this number its skipped */

   /** firstly, iterate over all consclassifiers */
   for( size_t classifier = 0; classifier < consclassescollection.size(); ++classifier )
   {
      /** check if there are to  many classes in this distribution and skip it if so */
      if ( consclassescollection[classifier]->getNClasses() > maximumnclasses)
      {
         std::cout << " the current consclass distribution includes " <<  consclassescollection[classifier]->getNClasses() << " classes but only " << maximumnclasses << " are allowed for calcCandidatesNBlocks()" << std::endl;
         continue;
      }

      /** get necessary data of current classifier */
      std::vector< std::vector<int> > subsetsOfConstypes = consclassescollection[classifier]->getAllSubsets( true, true, true );
      std::vector<int> nConssOfClasses = consclassescollection[classifier]->getNConssOfClasses();

      /** start with the cardinalities of the consclasses as candidates */
      for( size_t i = 0; i < nConssOfClasses.size(); ++i)
      {
        addCandidatesNBlocks(nConssOfClasses[i]);
      }

      /** continue with gcd of all cardinalities in this subset */
      for(size_t subset = 0; subset < subsetsOfConstypes.size(); ++subset)
      {
         int greatestCD = 1;

         if( subsetsOfConstypes[subset].size() == 0 || subsetsOfConstypes[subset].size() == 1 )
              continue;

         greatestCD = gcd(nConssOfClasses[subsetsOfConstypes[subset][0]], nConssOfClasses[subsetsOfConstypes[subset][1]]  );

         for( size_t i = 2; i < subsetsOfConstypes[subset].size() ; ++i)
         {
            greatestCD = gcd( greatestCD, nConssOfClasses[subsetsOfConstypes[subset][i]] );
         }

         addCandidatesNBlocks(greatestCD);
      }
   }

   /** secondly, iterate over all varclassifiers */
   for( size_t classifier = 0; classifier < varclassescollection.size(); ++classifier )
   {
      /** check if there are to  many classes in this distribution and skip it if so */
      if ( varclassescollection[classifier]->getNClasses() > maximumnclasses)
      {
         std::cout << " the current varclass distribution includes " <<  varclassescollection[classifier]->getNClasses() << " classes but only " << maximumnclasses << " are allowed for calcCandidatesNBlocks()" << std::endl;
         continue;
      }

      /** get necessary data of current classifier */
      std::vector< std::vector<int> > subsetsOfVartypes = varclassescollection[classifier]->getAllSubsets( true, true, true, true );
      std::vector<int> nVarsOfClasses = varclassescollection[classifier]->getNVarsOfClasses();

      /** start with the cardinalities of the varclasses as candidates */
      for( size_t i = 0; i < nVarsOfClasses.size(); ++i)
      {
         addCandidatesNBlocks(nVarsOfClasses[i]);
      }

      /** continue with gcd of all cardinalities in this subset */
      for(size_t subset = 0; subset < subsetsOfVartypes.size(); ++subset)
      {
         int greatestCD = 1;

         if( subsetsOfVartypes[subset].size() == 0 || subsetsOfVartypes[subset].size() == 1 )
              continue;

         greatestCD = gcd(nVarsOfClasses[subsetsOfVartypes[subset][0]], nVarsOfClasses[subsetsOfVartypes[subset][1]]  );

         for( size_t i = 2; i < subsetsOfVartypes[subset].size() ; ++i)
         {
            greatestCD = gcd( greatestCD, nVarsOfClasses[subsetsOfVartypes[subset][i]] );
         }

         addCandidatesNBlocks(greatestCD);
      }
   }
}

int Seeedpool::getNConssClassDistributions(){
   return (int) consclassescollection.size();
}

int* Seeedpool::getConssClassDistribution(
   int consclassdistr
)
{
   int nconss = consclassescollection[consclassdistr]->getNConss();
   int* output = new int[nconss];
   for ( int i = 0; i < nconss; ++i )
      output[i] = consclassescollection[consclassdistr]->getClassOfCons( i );
   return &output[0];
}

std::vector<int> Seeedpool::getConssClassDistributionVector(
   int consclassdistr
)
{
   int nconss = consclassescollection[consclassdistr]->getNConss();
   std::vector<int> output(nconss, 0);
   for ( int i = 0; i < nconss; ++i )
      output[i] = consclassescollection[consclassdistr]->getClassOfCons( i );
   return output;
}

int Seeedpool::getNClassesOfDistribution( int consclassdistr )
{
   return consclassescollection[consclassdistr]->getNClasses();
}

/** returns number of different constraint classifiers */
int Seeedpool::getNConsClassifiers()
{
   return (int) consclassescollection.size();
}

 /** returns pointer to a constraint classifier */
ConsClassifier* Seeedpool::getConsClassifier( int givenClassifierIndex )
{
   assert( 0 <= givenClassifierIndex && givenClassifierIndex < (int) consclassescollection.size() );

   return consclassescollection[givenClassifierIndex];
}

ConsClassifier* Seeedpool::createConsClassifierForSCIPConstypes()
 {
    /**
     * at first for every subset of constypes calculate gcd (greatest common divisors) of the corresponding number of occurrences
     */
    std::vector<consType> foundConstypes(0);
    std::vector<int> constypesIndices(0);
    std::vector<int> classForCons = std::vector<int>(getNConss(), -1);
    ConsClassifier* classifier;

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

    classifier = new ConsClassifier( scip, "constypes", (int) foundConstypes.size(), getNConss() );

    for( int c = 0; c < classifier->getNClasses(); ++c )
    {
       std::string name;
       std::stringstream text;
       switch (foundConstypes[c])
       {
       case linear:
          name = "linear";
          break;
       case knapsack:
          name = "knapsack";
          break;
       case varbound:
          name = "varbound";
          break;
       case setpacking:
          name = "setpacking";
          break;
       case setcovering:
          name = "setcovering";
          break;
       case setpartitioning:
          name = "setpartitioning";
          break;
       case logicor:
          name= "logicor";
          break;
       case sos1:
          name = "sos1";
          break;
       case sos2:
          name = "sos2";
          break;
       case unknown:
          name = "unknown";
          break;
       case nconsTypeItems:
          name = "nconsTypeItems";
          break;
       default:
          name = "newConstype";
          break;
       }
       classifier->setClassName( c, name.c_str() );
       text << "This class contains all constraints that are of (SCIP) constype \"" << name << "\".";
       classifier->setClassDescription( c, text.str().c_str() );
    }
    for( int i = 0; i < classifier->getNConss(); ++i )
    {
       classifier->assignConsToClass( i, classForCons[i] );
    }

    std::cout << " consclassifier scipconstypes:" << " yields a classification with " << foundConstypes.size()  << " different constraint classes" << std::endl;

    return classifier;
}

 ConsClassifier* Seeedpool::createConsClassifierForConsnamesDigitFreeIdentical()
  {
     /**
      * at first remove all digits from the consnames
      */
     std::vector<std::string> consnamesToCompare(getNConss(), "");
     std::vector<int> nConssConstype(0);
     std::vector<int> classForCons = std::vector<int>(getNConss(), -1);
     std::vector<std::string> nameClasses(0);
     ConsClassifier* classifier;


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

     classifier = new ConsClassifier( scip, "consnames", (int) nameClasses.size(), getNConss() );

     for( int c = 0; c < classifier->getNClasses(); ++c )
     {
        std::stringstream text;
        classifier->setClassName( c, nameClasses[c].c_str() );
        text << "This class contains all constraints with name \"" << nameClasses[c] << "\".";
        classifier->setClassDescription( c, text.str().c_str() );
     }
     for( int i = 0; i < classifier->getNConss(); ++i )
     {
        classifier->assignConsToClass( i, classForCons[i] );
     }

     std::cout << " consclassifier digit-reduced consnames (check for identity):" << " yields a classificiation with " << nameClasses.size()  << " different constraint classes" << std::endl;

     return classifier;
}

 ConsClassifier* Seeedpool::createConsClassifierForConsnamesLevenshteinDistanceConnectivity(
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

     std::stringstream classifierName;
     classifierName << "lev-dist-" << connectivity;
     ConsClassifier* classifier = new ConsClassifier( scip, classifierName.str().c_str(), 0, getNConss() );


     if (getNConss() > nmaxconss)
     {
        std::cout << " skipped levenshtein distance based constraint classes calculating since number of constraints " << getNConss() << " exceeds limit " << nmaxconss   << std::endl;
        return NULL;
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

           std::stringstream text;
           text << "This class contains all constraints with a name similar to \"" << consnamesToCompare[firstUnreached] << "\".";
           int newClass = classifier->addClass( consnamesToCompare[firstUnreached].c_str(), text.str().c_str(), BOTH );
           assert( newClass == currentClass );

        } // endwhile( !openConss.empty() )

        for( int i = 0; i < classifier->getNConss(); ++i )
        {
           classifier->assignConsToClass( i, classForCons[i] );
        }

        std::cout << " consclassifier levenshtein: connectivity of " << connectivity << " yields a classification with " << currentClass+1  << " different constraint classes" << std::endl;

        return classifier;
}


 ConsClassifier* Seeedpool::createConsClassifierForNNonzeros()
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

     if( false )
     {
        std::cout << " nNonzero : nConsWithNNonzero"  << std::endl;
        for( size_t i = 0; i < differentNNonzeros.size(); ++i )
        {
           std::cout << differentNNonzeros[i] << " : " << nconssforclass[i] << std::endl;
//           std::cout << " new consname : " << consnamesToCompare[i]  << std::endl;
//           std::cout << std::endl;
        }
     }

     ConsClassifier* classifier = new ConsClassifier( scip, "nonzeros", (int) differentNNonzeros.size(), getNConss() );

     for( int c = 0; c < classifier->getNClasses(); ++c )
     {
        std::stringstream text;
        text << differentNNonzeros[c];
        classifier->setClassName( c, text.str().c_str() );
        text.str("");
        text.clear();
        text << "This class contains all constraints with " << differentNNonzeros[c] << " nonzero coefficients.";
        classifier->setClassDescription( c, text.str().c_str() );
     }
     for( int i = 0; i < classifier->getNConss(); ++i )
     {
        classifier->assignConsToClass( i, classForCons[i] );
     }

     std::cout << " consclassifier nonzeros: comparison of number of nonzeros" << " yields a classification with " << differentNNonzeros.size()  << " different constraint classes" << std::endl;

     return classifier;
 }

/** adds a constraint classifier if it is no duplicate of an existing classifier */
void Seeedpool::addConsClassifier( ConsClassifier* givenClassifier )
{
   if ( givenClassifier != NULL )
   {
      /** check whether there already exists an equivalent consclassifier */
      ConsClassifier* equiv = NULL;

      for ( size_t i = 0; i < consclassescollection.size(); ++i )
      {
         if ( givenClassifier->classifierIsDuplicateOfClassifier( consclassescollection[i] ) )
         {
            equiv = consclassescollection[i];
            break;
         }
      }

      if ( equiv == NULL )
         consclassescollection.push_back( givenClassifier );
      else
         std::cout << "  consclassifier " << givenClassifier->getName() << " not considered since it offers the same structure as " << equiv->getName() << std::endl;
   }
}

/** adds constraint classifiers with a reduced number of classes */
void Seeedpool::reduceConsclasses()
{
   int maxnclasses = 9;

   if( getNConss() + getNVars() > 50000 )
      maxnclasses = 3;

   for( size_t classifierid = 0; classifierid < consclassescollection.size(); ++classifierid )
   {
      ConsClassifier* newclassifier = consclassescollection[classifierid]->reduceClasses( maxnclasses );

      if ( newclassifier != NULL )
      {
         std::cout <<  "add reduced version of consclassifier " << consclassescollection[classifierid]->getName() << " with " << maxnclasses << " classes" << std::endl;
         addConsClassifier( newclassifier );
      }
   }
}

/** returns number of different variable classifiers */
int Seeedpool::getNVarClassifiers()
{
   return (int) varclassescollection.size();
}

/** returns pointer to a variable classifier */
VarClassifier* Seeedpool::getVarClassifier( int givenClassifierIndex )
{
   assert( 0 <= givenClassifierIndex && givenClassifierIndex < (int) varclassescollection.size() );

   return varclassescollection[givenClassifierIndex];
}

VarClassifier* Seeedpool::createVarClassifierForSCIPVartypes()
{
   std::vector<SCIP_VARTYPE> foundVartypes(0);
   std::vector<int> vartypesIndices(0);
   std::vector<int> classForVars = std::vector<int>(getNVars(), -1);
   VarClassifier* classifier;

   for( int i = 0; i < getNVars(); ++i)
   {
      SCIP_VAR* var;
      bool found = false;
      var = getVarForIndex(i);
      SCIP_VARTYPE vT = SCIPvarGetType( var );
      size_t vartype;

      /** find vartype or not */
      for( vartype = 0; vartype < foundVartypes.size(); ++vartype )
      {
         if( foundVartypes[vartype] == vT )
         {
            found = true;
            break;
         }
      }
      if( !found )
      {
         foundVartypes.push_back( vT );
         classForVars[i] = foundVartypes.size() - 1;
      }
      else
         classForVars[i] = vartype;
   }

   classifier = new VarClassifier( scip, "vartypes", (int) foundVartypes.size(), getNVars() );

   for( int c = 0; c < classifier->getNClasses(); ++c )
   {
      std::string name;
      std::stringstream text;
      switch (foundVartypes[c])
      {
          case SCIP_VARTYPE_BINARY:
             name = "bin";
             break;
          case SCIP_VARTYPE_INTEGER:
             name = "int";
             break;
          case SCIP_VARTYPE_IMPLINT:
             name = "impl";
             break;
          case SCIP_VARTYPE_CONTINUOUS:
             name = "cont";
             break;
          default:
             name = "newVartype";
             break;
      }
      classifier->setClassName( c, name.c_str() );
      text << "This class contains all variables that are of (SCIP) vartype \"" << name << "\".";
      classifier->setClassDescription( c, text.str().c_str() );
   }
   for( int i = 0; i < classifier->getNVars(); ++i )
   {
      classifier->assignVarToClass( i, classForVars[i] );
   }

   std::cout << " varclassifier scipvartypes:" << " yields a classification with " << foundVartypes.size()  << " different variable classes" << std::endl;

   return classifier;
}

/** adds a variable classifier if it is no duplicate of an existing variable classifier */
void Seeedpool::addVarClassifier( VarClassifier* givenClassifier )
{
   if ( givenClassifier != NULL )
   {
      /** check whether there already exists an equivalent varclassifier */
      VarClassifier* equiv = NULL;

      for ( size_t i = 0; i < varclassescollection.size(); ++i )
      {
         if ( givenClassifier->classifierIsDuplicateOfClassifier( varclassescollection[i] ) )
         {
            equiv = varclassescollection[i];
            break;
         }
      }

      if ( equiv == NULL )
         varclassescollection.push_back( givenClassifier );
      else
         std::cout << "  varclassifier " << givenClassifier->getName() << " not considered since it offers the same structure as " << equiv->getName() << std::endl;
   }
}

/** adds variable classifiers with a reduced number of classes */
void Seeedpool::reduceVarclasses()
{
   int maxnclasses = 9;

   if( getNConss() + getNVars() > 50000 )
      maxnclasses = 3;

   for( size_t classifierid = 0; classifierid < varclassescollection.size(); ++classifierid )
   {
      VarClassifier* newclassifier = varclassescollection[classifierid]->reduceClasses( maxnclasses );

      if ( newclassifier != NULL )
      {
         std::cout <<  "add reduced version of consclassifier " << varclassescollection[classifierid]->getName() << " with " << maxnclasses << " classes" << std::endl;
         addVarClassifier( newclassifier );
      }
   }
}

std::vector<SeeedPtr> Seeedpool::removeSomeOneblockDecomps(
   std::vector<SeeedPtr> seeeds)
{
   std::vector<SeeedPtr> remainingSeeeds(0);
   std::vector<SeeedPtr> oneBlockSeeeds(0);

   int nmasterconssfirst = 1000;
   int nmasterconsssecond = 1001;

   for( size_t i = 0; i < seeeds.size(); ++i )
   {
      if( seeeds[i]->getNBlocks() == 1)
      {
         if(seeeds[i]->getNMasterconss() < nmasterconssfirst )
         {
            nmasterconsssecond = nmasterconssfirst;
            nmasterconssfirst = seeeds[i]->getNMasterconss();
         }else if(seeeds[i]->getNMasterconss() < nmasterconsssecond )
            nmasterconsssecond = seeeds[i]->getNMasterconss();

      }
      else
         remainingSeeeds.push_back(seeeds[i]);
   }

   for(int i = 0; i < (int) seeeds.size(); ++i)
   {
      if( seeeds[i]->getNBlocks() == 1 && ( seeeds[i]->getNMasterconss() == nmasterconssfirst || seeeds[i]->getNMasterconss() == nmasterconsssecond ) )
         remainingSeeeds.push_back(seeeds[i]);
      else if (seeeds[i]->getNBlocks() == 1)
         oneBlockSeeeds.push_back(seeeds[i]);
   }

   for(int i = 0; i < (int) oneBlockSeeeds.size(); ++i)
   {
      delete oneBlockSeeeds[i];
      oneBlockSeeeds[i] = NULL;
   }

   return remainingSeeeds;
}

SCIP_RETCODE Seeedpool::writeFamilyTreeLatexFile(
   const char* filename,                                 /* filename the output should be written to (including directory) */
   const char* workfolder,                               /* directory in which should be worked */
   std::vector<SeeedPtr> seeeds,                         /* vector of seeed pointers the  family tree should be constructed for */
   SCIP_Bool draft
){

   std::ofstream ofs;
   int curr = -1;
   int currheight = 0;
   SCIP_Real firstsibldist = -1.;

   std::stringstream preambel;
   std::string closing = "\\end{tikzpicture}\n\\end{document}";

   /* collectiotn of treeseeds */
   std::vector<SeeedPtr> treeseeeds(0);
   std::vector<int> treeseeedids(0);
   std::vector<SCIP_Bool> isseeedintree(allrelevantseeeds.size(), FALSE );

   int root = -1;
   std::vector<int> parents(allrelevantseeeds.size(), -1);
   std::vector<std::vector<int> > childs (allrelevantseeeds.size(), std::vector<int>(0));
   std::vector<std::vector<SCIP_Bool> > childsfinished(allrelevantseeeds.size(), std::vector<SCIP_Bool>(0));
   std::vector<SCIP_Bool> visited(allrelevantseeeds.size(), FALSE);

   /** check allrelevant seeeds **/
   for( size_t s = 0; s < allrelevantseeeds.size(); ++s )
   {
      assert(allrelevantseeeds[s] == NULL || (int) s == allrelevantseeeds[s]->getID() );
   }

   /** 1) find relevant seeeds in tree and build tree */
   for( size_t s = 0; s < seeeds.size(); ++s )
   {
      int currid;
      if ( seeeds[s] == NULL )
         continue;
      currid = seeeds[s]->getID();
      if( !isseeedintree[seeeds[s]->getID()] )
      {
         isseeedintree[seeeds[s]->getID()] = TRUE;
         treeseeeds.push_back( seeeds[s]);
         treeseeedids.push_back(seeeds[s]->getID());
      }
      else
         break;

      for( size_t i = 0; i < seeeds[s]->listofancestorids.size(); ++i )
      {
         int ancestorid;
         ancestorid = seeeds[s]->listofancestorids[seeeds[s]->listofancestorids.size() - i -1];
         parents[currid] = ancestorid;
         childs[ancestorid].push_back(currid);
         childsfinished[ancestorid].push_back(FALSE);

         if( !isseeedintree[ancestorid] )
         {
            isseeedintree[ancestorid] = TRUE;
            treeseeeds.push_back( allrelevantseeeds[ancestorid] );
            treeseeedids.push_back(ancestorid);
            if( i == seeeds[s]->listofancestorids.size() -1 )
            {
               root = ancestorid;
            }
            currid = ancestorid;
         }
         else
            break;
      }
   }

   for( size_t i = 0; i < treeseeeds.size(); ++i )
   {
      SeeedPtr seeed = treeseeeds[i];
      std::stringstream decompfilename;

      seeed = treeseeeds[i];
      decompfilename << workfolder << "/" << getSeeedFolderLatex(seeed);

      seeed->showScatterPlot(this, TRUE, decompfilename.str().c_str(), draft );
   }

   //  finishedSeeeds[0]->showScatterPlot(this, TRUE, "./testdecomp/001.pdf") ;

   firstsibldist = 1. / (childs[root].size() - 1 );
   preambel.precision(2);

   preambel << "\\documentclass[a4paper,landscape]{scrartcl}\n\\usepackage{fancybox}\n\\usepackage{tikz}";
   preambel << "\n\\usetikzlibrary{positioning}\n\\title{Detection Tree}\n\\date{}\n\\begin{document}\n\n";
   preambel << "\\begin{tikzpicture}[level/.style={sibling distance=" << firstsibldist << "\\textwidth/#1}, level distance=12em, ->, dashed]\n\\node";



   /** start writing file */
   ofs.open (filename, std::ofstream::out );
   ofs << preambel.str();

   /** iterate tree and write file */
   curr = root;
   while ( curr != -1 )
   {
      if( !visited[curr] )
      {
         /** write node */
         ofs << writeSeeedIncludeLatex( allrelevantseeeds[curr], workfolder );
         /* set node visited */
         visited[curr] = TRUE;
         if( parents[curr] != -1 )
            finishnextchild(childs[parents[curr]], childsfinished[parents[curr]], curr);

      }
      if ( unfinishedchildexists(childsfinished[curr] ) )
      {
         int unfinishedchild = getfirstunfinishedchild(childsfinished[curr], childs[curr] );
         /* is first child unfinihsed? */
         //         if( unfinishedchild == childs[curr][0] )
         ofs << " child { node " ;
         curr = unfinishedchild;
         ++currheight;
      }
      else
      {
         if ( parents[curr] != -1 ){
            ofs << writeSeeedDetectorChainInfoLatex( allrelevantseeeds[curr], currheight, helpvisucounter);
            ++helpvisucounter;
         }
         --currheight;
         curr = parents[curr];
         if( curr != -1)
            ofs << " } " ;
      }
   }

   ofs << ";" << std::endl;
   for( size_t i = 0; i < treeseeeds.size(); ++i)
   {
      ofs << writeSeeedInfoLatex( treeseeeds[i] );
   }

   ofs << closing << std::endl;

   ofs.close();

   return SCIP_OKAY;
}



/**
 * creates a decomposition for a given seeed
 */
SCIP_RETCODE Seeedpool::createDecompFromSeeed(
   SeeedPtr       seeed,                                 /** seeed the decomposition is created for */
   DEC_DECOMP**   newdecomp                              /** the new decomp created from the seeed */
  )
{

   char detectorchaininfo[SCIP_MAXSTRLEN];

   SCIP_HASHMAP* vartoblock;
   SCIP_HASHMAP* constoblock;
   SCIP_HASHMAP* varindex;
   SCIP_HASHMAP* consindex;

   SCIP_VAR*** stairlinkingvars;
   SCIP_VAR*** subscipvars;
   SCIP_VAR**  linkingvars;
   SCIP_CONS**  linkingconss;
   SCIP_CONS*** subscipconss;

   DEC_SCORES scores;



   int* nsubscipconss;
   int* nsubscipvars;
   int* nstairlinkingvars;
   int  nlinkingvars;

   int varcounter = 1;  /* in varindex counting starts with 1 */
   int conscounter = 1; /* in consindex counting starts with 1 */
   int counterstairlinkingvars = 0;

   int size;

   assert(seeed->checkConsistency(this) );

   /* create decomp data structure */
   SCIP_CALL_ABORT( DECdecompCreate( scip, newdecomp) );

   //           seeed->displayConss();
   //     if(seeed->detectorChain.size() > 2)
   //    seeed->showScatterPlot(this);


   /** set nblocks */
   DECdecompSetNBlocks( *newdecomp, seeed->getNBlocks() );

   //detectorchaininfo ;
   /** set constraints */
   if( seeed->getNMasterconss( )  != 0 )
      SCIP_CALL_ABORT ( SCIPallocBufferArray(scip, &linkingconss, seeed->getNMasterconss() ) );
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
      DECdecompSetLinkingconss(scip, *newdecomp, linkingconss, seeed->getNMasterconss());
   else
      linkingconss = NULL;
   /* set block constraints */
   for ( int b = 0; b < seeed->getNBlocks(); ++b )
   {
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &subscipconss[b], seeed->getNConssForBlock(b) ) );
      nsubscipconss[b] = seeed->getNConssForBlock(b);
      for( int c = 0; c < seeed->getNConssForBlock(b); ++c )
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


   DECdecompSetSubscipconss(scip, *newdecomp, subscipconss, nsubscipconss );

   DECdecompSetConstoblock(*newdecomp, constoblock);
   DECdecompSetConsindex(*newdecomp, consindex);

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
      SCIP_CALL_ABORT( SCIPhashmapInsert(varindex, scipvar, (void*) (size_t) varcounter) );
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

   DECdecompSetSubscipvars(scip, *newdecomp, subscipvars, nsubscipvars);
   DECdecompSetStairlinkingvars(scip, *newdecomp, stairlinkingvars, nstairlinkingvars);
   DECdecompSetLinkingvars(scip, *newdecomp, linkingvars, nlinkingvars);
   DECdecompSetVarindex(*newdecomp, varindex);
   DECdecompSetVartoblock(*newdecomp, vartoblock) ;



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
   (*newdecomp)->sizedetectorchain = ndetectors;
   size = SCIPcalcMemGrowSize(scip, (*newdecomp)->sizedetectorchain);
   SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &(*newdecomp)->detectorchain, size) ); /** free in decomp.c:469 */
   for( int k = 0; k < ndetectors; ++k )
   {
      if(k != ndetectors-1 || !seeed->getFinishedByFinisher() )
      {
         //          std::cout << " added detector of " << i << "-th seeed to its detetcor chain" << std::endl;
         (*newdecomp)->detectorchain[k] = seeed->getDetectorchain()[k];
      }else
         (*newdecomp)->detectorchain[k] = seeed->getDetectorchain()[k];
   }


   /** set statistical detector chain data */

   DECdecompSetSeeedID(*newdecomp, seeed->getID() );
   if( seeed->getNDetectors() > 0 )
   {
      DECdecompSetDetectorClockTimes(scip, *newdecomp, &(seeed->detectorClockTimes[0]) );
      DECdecompSetDetectorPctVarsToBorder(scip, *newdecomp, &(seeed->pctVarsToBorder[0] ) );
      DECdecompSetDetectorPctVarsToBlock(scip, *newdecomp, &(seeed->pctVarsToBlock[0] ) );
      DECdecompSetDetectorPctVarsFromOpen(scip, *newdecomp, &(seeed->pctVarsFromFree[0] ) );
      DECdecompSetDetectorPctConssToBorder(scip, *newdecomp, &(seeed->pctConssToBorder[0] ) );
      DECdecompSetDetectorPctConssToBlock(scip, *newdecomp, &(seeed->pctConssToBlock[0] ) );
      DECdecompSetDetectorPctConssFromOpen(scip, *newdecomp, &(seeed->pctConssFromFree[0] ) );
      DECdecompSetNNewBlocks(scip, *newdecomp, &(seeed->nNewBlocks[0] ) );
   }

   if( seeed->detectorchainstring == NULL)
   {
      seeed->buildDecChainString();
   }
   SCIP_CALL(DECdecompSetDetectorChainString(scip, *newdecomp, seeed->detectorchainstring ) );

   /** set dectype */
   if( (*newdecomp)->nlinkingvars == seeed->getNTotalStairlinkingvars() && (*newdecomp)->nlinkingconss == 0 && DECdecompGetNLinkingvars((*newdecomp)) > 0)
   {
      (*newdecomp)->type = DEC_DECTYPE_STAIRCASE;
   }
   else if( (*newdecomp)->nlinkingvars > 0 || seeed->getNTotalStairlinkingvars() > 0 )
   {
      (*newdecomp)->type = DEC_DECTYPE_ARROWHEAD;
   }
   else if( (*newdecomp)->nlinkingconss > 0)
   {
      (*newdecomp)->type = DEC_DECTYPE_BORDERED;
   }
   else if( (*newdecomp)->nlinkingconss == 0 && seeed->getNTotalStairlinkingvars() == 0)
   {
      (*newdecomp)->type = DEC_DECTYPE_DIAGONAL;
   }
   else
   {
      (*newdecomp)->type = DEC_DECTYPE_UNKNOWN;
   }

   ndecompositions++;

   SCIP_CALL(DECevaluateDecomposition(scip, *newdecomp, &scores) );

   assert(scores.maxwhitescore == seeed->getMaxWhiteScore() );

   assert(DECdecompCheckConsistency(scip, (*newdecomp) ) );

   assert(!SCIPhashmapIsEmpty((*newdecomp)->constoblock));
   assert(!SCIPhashmapIsEmpty((*newdecomp)->vartoblock));


   return SCIP_OKAY;
}

/**
 * creates a seeed for a given decomposition
 */
SCIP_RETCODE Seeedpool::createSeeedFromDecomp(
      DEC_DECOMP* decomp,                                    /** decomposition the seeed is created for */
      SeeedPtr*   newseeed                                   /** the new seeed created from the decomp */
  )
{
   SeeedPtr seeed;

   seeed = NULL;

   return SCIP_OKAY;

}

/**
 * returns transformation information
 */
SCIP_Bool Seeedpool::getTransformedInfo(
   )
{
   return transformed;
}



} /* namespace gcg */
