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

/**@file   clscons_consnamelevenshtein.cpp
 * 
 * @brief classifies constraints according to levenshtein distance graph of their names
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/clscons_consnamelevenshtein.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include <vector>
#include <stdio.h>
#include <sstream>
#include <queue>

#include "gcg/class_detprobdata.h"

#include "gcg/class_conspartition.h"
#include "gcg/scip_misc.h"

/* classifier properties */
#define CLSCONS_NAME              "consnamelevenshtein"       /**< name of classifier */
#define CLSCONS_DESC              "constraint names (according to levenshtein distance graph)"     /**< short description of classification*/
#define CLSCONS_PRIORITY          0

#define CLSCONS_ENABLED           FALSE


/*
 * Data structures
 */

/** classifier handler data */
struct GCG_ClassifierData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * classifier callback methods
 */

/** destructor of classifier to free user data (called when GCG is exiting) */
#define classifierFree NULL

/** returns levenshtein distance between two strings */
int calcLevenshteinDistance(
   std::string s,
   std::string t
   )
{
   /* easy cases */
   if( s.compare( t ) == 0 )
      return 0;
   if( s.length() == 0 )
      return (int) t.length();
   if( t.length() == 0 )
      return (int) s.length();

   /* vectors to store integer distances */
   std::vector<int> prev( t.length() + 1 );
   std::vector<int> curr( t.length() + 1 );

   /* initialize prev (previous row of distances) */
   for( size_t i = 0; i < prev.size(); ++i )
   {
      prev[i] = (int) i;
   }
   for( size_t i = 0; i < s.length(); ++i )
   {
      /* calculate curr (row distances) from the previous one */

      curr[0] = (int) i + 1;

      /* fill remaining of row using 'Bellman' equality */
      for( size_t j = 0; j < t.length(); ++j )
      {
         int cost = ( s[i] == t[j] ) ? 0 : 1;
         curr[j + 1] = std::min( curr[j] + 1, std::min( prev[j + 1] + 1, prev[j] + cost ) );
      }

      /* copy curr to prev for next iteration */
      for( size_t j = 0; j < prev.size(); ++j )
         prev[j] = curr[j];
   }

   return curr[t.length()];
}



static
GCG_DECL_CONSCLASSIFY(classifierClassify) {
   gcg::DETPROBDATA* detprobdata;
   SCIP* origprob = GCGgetOrigprob(gcg);
   if( transformed )
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataPresolved(gcg);
   }
   else
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataOrig(gcg);
   }

   std::vector < std::string > consnamesToCompare(detprobdata->getNConss(), "");
   std::vector<int> nConssConstype;
   std::vector<int> classForCons(detprobdata->getNConss(), - 1);
   std::vector<bool> alreadyReached(detprobdata->getNConss(), false);
   std::queue<int> helpqueue;
   int nUnreachedConss = detprobdata->getNConss();
   int currentClass = - 1;
   int nmaxconss = 5000;
   int connectivity = 1;

   std::stringstream classifierName;
   classifierName << "lev-dist-" << connectivity;
   gcg::ConsPartition* classifier = new gcg::ConsPartition(gcg, classifierName.str().c_str(), 0, detprobdata->getNConss());

   /* if number of conss exceeds this number, skip calculating such a classifier */
   if( detprobdata->getNConss() > nmaxconss )
   {

      SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, " skipped levenshtein distance based constraint classes calculating since number of constraints  %d  exceeds limit %d \n", detprobdata->getNConss(), nmaxconss );
      delete classifier;
      return SCIP_ERROR;
   }

   std::vector<std::vector<int>> levenshteindistances(detprobdata->getNConss(), std::vector<int>( detprobdata->getNConss(), - 1 ));

   /* read consnames */
   for( int i = 0; i < detprobdata->getNConss(); ++ i )
   {
      consnamesToCompare[i] = std::string(SCIPconsGetName(detprobdata->getCons(i)));
   }

   /* calculate levenshtein distances pairwise */
   for( int i = 0; i < detprobdata->getNConss(); ++ i )
   {
      for( int j = i + 1; j < detprobdata->getNConss(); ++ j )
      {
         levenshteindistances[i][j] = calcLevenshteinDistance(consnamesToCompare[i], consnamesToCompare[j]);
         levenshteindistances[j][i] = levenshteindistances[i][j];
      }
   }

   /* repeat doing breadth first search until every constraint is assigned to a class */
   while( nUnreachedConss > 0 )
   {
      int firstUnreached = - 1;
      currentClass ++;
      assert( helpqueue.empty() );
      for( int i = 0; i < detprobdata->getNConss(); ++ i )
      {
         if( classForCons[i] == - 1 )
         {
            firstUnreached = i;
            break;
         }
      }

      helpqueue.push(firstUnreached);
      alreadyReached[firstUnreached] = true;
      classForCons[firstUnreached] = currentClass;
      -- nUnreachedConss;

      /* consider all constraints which are connected to the current constraint by means of levenshtein distance */
      while( ! helpqueue.empty() )
      {
         int nodecons = helpqueue.front();
         helpqueue.pop();
         for( int j = 0; j < detprobdata->getNConss(); ++ j )
         {

            if( alreadyReached[j] )
               continue;

            if( j == nodecons )
               continue;

            if( levenshteindistances[j][nodecons] > connectivity )
               continue;

            alreadyReached[j] = true;
            classForCons[j] = currentClass;
            -- nUnreachedConss;
            helpqueue.push(j);
         }
      }

      /* create a new class with found constraints in ConsPartition*/
      std::stringstream text;
      text << "This class contains all constraints with a name similar to \"" << consnamesToCompare[firstUnreached] << "\".";
      classifier->addClass(consnamesToCompare[firstUnreached].c_str(), text.str().c_str(), gcg::BOTH);
   }

   /* assign constraint indices to classes */
   for( int i = 0; i < detprobdata->getNConss(); ++ i )
   {
      classifier->assignConsToClass(i, classForCons[i]);
   }

   SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, " Consclassifier levenshtein: connectivity of %d yields a classification with %d different constraint classes. \n", connectivity, currentClass + 1);

   detprobdata->addConsPartition(classifier);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

SCIP_RETCODE GCGincludeConsClassifierConsnameLevenshtein(
   GCG* gcg                /**< GCG data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata = NULL;

   SCIP_CALL(
      GCGincludeConsClassifier(gcg, CLSCONS_NAME, CLSCONS_DESC, CLSCONS_PRIORITY, CLSCONS_ENABLED, classifierdata,
                               classifierFree, classifierClassify));

   return SCIP_OKAY;
}
