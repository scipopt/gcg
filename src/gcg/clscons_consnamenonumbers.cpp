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

/**@file   clscons_consnamenonumbers.cpp
 * 
 * @brief classifies constraints according to names (without digits)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/clscons_consnamenonumbers.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include <vector>
#include <stdio.h>
#include <sstream>

#include "gcg/class_detprobdata.h"

#include "gcg/class_conspartition.h"
#include "gcg/scip_misc.h"

/* classifier properties */
#define CLSCONS_NAME                  "consnamenonumbers"       /**< name of classifier */
#define CLSCONS_DESC                  "constraint names (remove digits; check for identity)"     /**< short description of classification*/
#define CLSCONS_PRIORITY              0

#define CLSCONS_ENABLED               FALSE


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

/** removes all digits from string str */
void removeDigits(
   char *str,
   int *nremoved
)
{
   char digits[11] = "0123456789";
   * nremoved = 0;

   for( int i = 0; i < 10; ++ i )
   {
      char digit = digits[i];
      size_t j = 0;
      while( j < strlen( str ) )
      {
         if( str[j] == digit )
         {
            * nremoved = * nremoved + 1;
            for( size_t k = j; k < strlen( str ); ++ k )
            {
               str[k] = str[k + 1];
            }
         }
         else
            ++ j;
      }
   }
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

   std::vector < std::string > consnamesToCompare( detprobdata->getNConss(), "" );
   std::vector<int> nConssConstype( 0 );
   std::vector<int> classForCons = std::vector<int>( detprobdata->getNConss(), - 1 );
   std::vector < std::string > nameClasses( 0 );
   gcg::ConsPartition* classifier;

   /* firstly, remove all digits from the consnames */
   for( int i = 0; i < detprobdata->getNConss(); ++ i )
   {
      int nremoved;
      char consname[SCIP_MAXSTRLEN];
      strcpy(consname, SCIPconsGetName(detprobdata->getCons(i)));

      removeDigits(consname, &nremoved);
      consnamesToCompare[i] = std::string(consname);
   }

   for( int i = 0; i < detprobdata->getNConss(); ++ i )
   {
      /* check if string belongs to an existing name class */
      bool belongstoexistingclass = false;

      for( size_t j = 0; j < nameClasses.size(); ++ j )
      {
         if( nameClasses[j] == consnamesToCompare[i] )
         {
            belongstoexistingclass = true;
            classForCons[i] = (int) j;
            nConssConstype[j] ++;
            break;
         }
      }
      /* if not, create a new class */
      if( !belongstoexistingclass )
      {
         nameClasses.push_back(consnamesToCompare[i]);
         nConssConstype.push_back(1);
         classForCons[i] = (int) nameClasses.size() - 1;
      }
   }

   /* secondly, use these information to create a ConsPartition */
   classifier = new gcg::ConsPartition(gcg, "consnames", (int) nameClasses.size(), detprobdata->getNConss());

   /* set all class names and descriptions */
   for( int c = 0; c < classifier->getNClasses(); ++ c )
   {
      std::stringstream text;
      classifier->setClassName(c, nameClasses[c].c_str());
      text << "This class contains all constraints with name \"" << nameClasses[c] << "\".";
      classifier->setClassDescription(c, text.str().c_str());
   }

   /* copy the constraint assignment information found in first step */
   for( int i = 0; i < classifier->getNConss(); ++ i )
   {
      classifier->assignConsToClass(i, classForCons[i]);
   }

   SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, " Consclassifier \"%s\" yields a classification with %d different constraint classes \n", classifier->getName(), classifier->getNClasses());

   detprobdata->addConsPartition(classifier);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

SCIP_RETCODE GCGincludeConsClassifierForConsnamesDigitFreeIdentical(
   GCG*                 gcg                /**< GCG data structure */
) {
   GCG_CLASSIFIERDATA* classifierdata = NULL;

   SCIP_CALL(
      GCGincludeConsClassifier(gcg, CLSCONS_NAME, CLSCONS_DESC, CLSCONS_PRIORITY, CLSCONS_ENABLED, classifierdata,
                               classifierFree, classifierClassify));

   return SCIP_OKAY;
}
