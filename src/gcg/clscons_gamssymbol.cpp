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

/**@file    clscons_gamssymbol.cpp
 * 
 * @brief   gamssymbol constraint classifier (classifies by corresponding GAMS symbol)
 * @author  Stefanie Koß
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG

#include "gcg/clscons_gamssymbol.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include <vector>
#include <string>
#include <map>
#include <stdio.h>
#include <sstream>

#include "gcg/class_detprobdata.h"

#include "gcg/class_conspartition.h"
#include "gcg/scip_misc.h"

/* classifier properties */
#define CLSCONS_NAME                  "gamssymbol"                 /**< name of classifier */
#define CLSCONS_DESC                  "symbol in GAMS file"        /**< short description of classification*/
#define CLSCONS_PRIORITY              0

#define CLSCONS_ENABLED               TRUE


/*
 * Data structures
 */
struct GCG_ClassifierData
{
   std::map<std::string, int>*      constosymbol;             /**< maps constraint name to the corresponding symbol index */
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * classifier callback methods
 */

/** destructor of classifier to free user data (called when GCG is exiting) */
static
GCG_DECL_FREECONSCLASSIFIER(classifierFree)
{
   GCG_CLASSIFIERDATA* classifierdata;
   SCIP* origprob = GCGgetOrigprob(gcg);

   assert(origprob != NULL);

   classifierdata = GCGconsClassifierGetData(classifier);
   assert(classifierdata != NULL);
   assert(strcmp(GCGconsClassifierGetName(classifier), CLSCONS_NAME) == 0);

   delete classifierdata->constosymbol;

   SCIPfreeBlockMemory(origprob, &classifierdata);

   return SCIP_OKAY;
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

   int ncons = detprobdata->getNConss();
   std::vector<int> nconssForClass( 0 );        // [i] holds number of constraints for class i
   std::vector<int> symbolidxForClass( 0 );     // [i] holds symbol index for class i
   std::vector<int> classForCons( ncons, - 1 ); // [i] holds class index for constraint i -> indexing over detection internal constraint array!
   int counterClasses = 0;

   GCG_CONSCLASSIFIER* classifier = GCGfindConsClassifier(gcg, CLSCONS_NAME);
   assert(classifier != NULL);

   GCG_CLASSIFIERDATA* classdata = GCGconsClassifierGetData(classifier);
   assert(classdata != NULL);

   /* firstly, assign all constraints to classindices */
   // iterate over constraints in detection and lookup in classdata->constosymbol
   // iterating over classdata->constosymbol and lookup constraints with getIndexForCons fails with assertion if constraint is not found -> should return error value?
   for( int consid = 0; consid < detprobdata->getNConss(); ++ consid )
   {
      // int consid = detprobdata->getIndexForCons(iter.second);
      SCIP_CONS* cons = detprobdata->getCons(consid);
      std::string consname = std::string( SCIPconsGetName( cons ) );

      auto symbolidxiter = classdata->constosymbol->find(consname);
      int symbolidx;
      if( symbolidxiter != classdata->constosymbol->end() )
      {
         symbolidx = symbolidxiter->second;
      }
      else
      {
         symbolidx = -1;
      }
      
      bool classfound = false;

      /* check if class for symbol index exists */
      for( size_t classid = 0; classid < symbolidxForClass.size(); ++classid )
      {
         if( symbolidx == symbolidxForClass[classid] )
         {
            classfound = true;
            classForCons[consid] = (int) classid;
            ++nconssForClass[classid];
            break;
         }
      }

      /* if not, create a new class index */
      if( !classfound )
      {
         classForCons[consid] = counterClasses;
         ++counterClasses;
         symbolidxForClass.push_back( symbolidx );
         nconssForClass.push_back( 1 );
      }
   }
   assert( counterClasses == (int) symbolidxForClass.size() );

   /* secondly, use these information to create a ConsPartition */
   gcg::ConsPartition* partition = new gcg::ConsPartition(gcg, "gamssymbols", counterClasses, detprobdata->getNConss() );

   /* set class names and descriptions of every class */
   for( int c = 0; c < partition->getNClasses(); ++ c )
   {
      std::stringstream text;
      text << symbolidxForClass[c];
      partition->setClassName( c, text.str().c_str() );
      text.str( "" );
      text.clear();
      text << "This class contains all constraints with gams symbol index" << symbolidxForClass[c] << ".";
      partition->setClassDescription( c, text.str().c_str() );
   }

   /* copy the constraint assignment information found in first step */
   for( int i = 0; i < partition->getNConss(); ++ i )
   {
      partition->assignConsToClass( i, classForCons[i] );
   }
   SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, " Consclassifier \"%s\" yields a classification with %d  different constraint classes \n", partition->getName(), partition->getNClasses() );

   detprobdata->addConsPartition(partition);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

/** adds an entry to clsdata->constosymbol */
SCIP_RETCODE GCGconsClassifierGamssymbolAddEntry(
   GCG_CONSCLASSIFIER*   classifier,
   SCIP_CONS*            cons,
   int                   symbolIdx
   )
{
   assert(classifier != NULL);
   GCG_CLASSIFIERDATA* classdata = GCGconsClassifierGetData(classifier);
   assert(classdata != NULL);

   std::string consname = SCIPconsGetName( cons );
   classdata->constosymbol->insert({consname, symbolIdx});

   return SCIP_OKAY;
}

/** creates the handler for gamssymbol classifier and includes it in SCIP */
SCIP_RETCODE GCGincludeConsClassifierGamssymbol(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata = NULL;
   SCIP* origprob = GCGgetOrigprob(gcg);

   SCIP_CALL( SCIPallocBlockMemory(origprob, &classifierdata) );
   assert(classifierdata != NULL);
   classifierdata->constosymbol = new std::map<std::string, int>();

   SCIP_CALL( GCGincludeConsClassifier(gcg, CLSCONS_NAME, CLSCONS_DESC, CLSCONS_PRIORITY, CLSCONS_ENABLED, classifierdata, classifierFree, classifierClassify) );

   return SCIP_OKAY;
}
