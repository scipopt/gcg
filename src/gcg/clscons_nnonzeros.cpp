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

/**@file   clscons_nnonzeros.cpp
 * 
 * @brief classifies constraints according to their nonzero entries
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/clscons_nnonzeros.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include <vector>
#include <stdio.h>
#include <sstream>

#include "gcg/class_detprobdata.h"

#include "gcg/class_conspartition.h"
#include "gcg/scip_misc.h"

/* classifier properties */
#define CLSCONS_NAME                  "nnonzeros"       /**< name of classifier */
#define CLSCONS_DESC                  "nnonezero entries"     /**< short description of classification*/
#define CLSCONS_PRIORITY              0

#define CLSCONS_ENABLED               TRUE


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

static
GCG_DECL_CONSCLASSIFY(classifierClassify)
{
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

   std::vector<int> nconssforclass( 0 );
   std::vector<int> differentNNonzeros( 0 );
   std::vector<int> classForCons( detprobdata->getNConss(), - 1 );
   int counterClasses = 0;

   /* firstly, assign all constraints to classindices */
   for( int i = 0; i < detprobdata->getNConss(); ++ i )
   {
      int consnnonzeros = detprobdata->getNVarsForCons( i );
      bool nzalreadyfound = false;

      /* check if number of nonzeros belongs to an existing class index */
      for( size_t nzid = 0; nzid < differentNNonzeros.size(); ++ nzid )
      {
         if( consnnonzeros == differentNNonzeros[nzid] )
         {
            nzalreadyfound = true;
            classForCons[i] = (int) nzid;
            ++ nconssforclass[nzid];
            break;
         }
      }

      /* if not, create a new class index */
      if( ! nzalreadyfound )
      {
         classForCons[i] = counterClasses;
         ++ counterClasses;
         differentNNonzeros.push_back( consnnonzeros );
         nconssforclass.push_back( 1 );
      }
   }

   /* secondly, use these information to create a ConsPartition */
   gcg::ConsPartition* classifier = new gcg::ConsPartition(gcg, "nonzeros", (int) differentNNonzeros.size(), detprobdata->getNConss() );

   /* set class names and descriptions of every class */
   for( int c = 0; c < classifier->getNClasses(); ++ c )
   {
      std::stringstream text;
      text << differentNNonzeros[c];
      classifier->setClassName( c, text.str().c_str() );
      text.str( "" );
      text.clear();
      text << "This class contains all constraints with " << differentNNonzeros[c] << " nonzero coefficients.";
      classifier->setClassDescription( c, text.str().c_str() );
   }

   /* copy the constraint assignment information found in first step */
   for( int i = 0; i < classifier->getNConss(); ++ i )
   {
      classifier->assignConsToClass( i, classForCons[i] );
   }
   SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, " Consclassifier \"%s\" yields a classification with %d  different constraint classes \n", classifier->getName(), classifier->getNClasses() );

   detprobdata->addConsPartition(classifier);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

/** creates the handler for classifier for Nnonzeros and includes it in SCIP */
SCIP_RETCODE GCGincludeConsClassifierNNonzeros(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata = NULL;

   SCIP_CALL(
      GCGincludeConsClassifier(gcg, CLSCONS_NAME, CLSCONS_DESC, CLSCONS_PRIORITY, CLSCONS_ENABLED, classifierdata,
                               classifierFree, classifierClassify) );

   return SCIP_OKAY;
}
