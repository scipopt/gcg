/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   clscons_scipconstypes.cpp
 * 
 * @brief classifies constraints according to their scip constypes
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/clscons_scipconstypes.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include <vector>
#include <stdio.h>
#include <sstream>

#include "gcg/class_detprobdata.h"

#include "gcg/class_conspartition.h"
#include "gcg/scip_misc.h"

/* classifier properties */
#define CLSCONS_NAME                  "scipconstype"       /**< name of classifier */
#define CLSCONS_DESC                  "scip constypes"     /**< short description of classification*/
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
GCG_DECL_CONSCLASSIFY(classifierClassify) {
   gcg::DETPROBDATA *detprobdata;
   SCIP* origprob = GCGgetOrigprob(gcg);
   if( transformed )
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataPresolved(gcg);
   }
   else
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataOrig(gcg);
   }

   std::vector<consType> foundConstypes(0);
   std::vector<int> constypesIndices(0);
   std::vector<int> classForCons = std::vector<int>(detprobdata->getNConss(), -1);
   gcg::ConsPartition *classifier;

   /* firstly, assign all constraints to classindices */
   for (int i = 0; i < detprobdata->getNConss(); ++i) {
      SCIP_CONS *cons;
      bool found = false;
      cons = detprobdata->getCons(i);
      consType cT = GCGconsGetType(origprob, cons);
      size_t constype;

      /* check whether the constraint's constype is new */
      for (constype = 0; constype < foundConstypes.size(); ++constype) {
         if (foundConstypes[constype] == cT) {
            found = true;
            break;
         }
      }
      /* if it is new, create a new classindex */
      if (!found) {
         foundConstypes.push_back(GCGconsGetType(origprob, cons));
         classForCons[i] = (int) foundConstypes.size() - 1;
      } else
         classForCons[i] = (int) constype;
   }

   /* secondly, use these information to create a ConsPartition */
   classifier = new gcg::ConsPartition(gcg, "constypes", (int) foundConstypes.size(), detprobdata->getNConss());

   /* set class names and descriptions of every class */
   for (int c = 0; c < classifier->getNClasses(); ++c) {
      std::string name;
      std::stringstream text;
      switch (foundConstypes[c]) {
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
            name = "logicor";
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
      classifier->setClassName(c, name.c_str());
      text << "This class contains all constraints that are of (SCIP) constype \"" << name << "\".";
      classifier->setClassDescription(c, text.str().c_str());
   }

   /* copy the constraint assignment information found in first step */
   for (int i = 0; i < classifier->getNConss(); ++i) {
      classifier->assignConsToClass(i, classForCons[i]);
   }

   SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL,
                   " Consclassifier \"%s\" yields a classification with %d different constraint classes \n",
                   classifier->getName(), (int) foundConstypes.size());

   detprobdata->addConsPartition(classifier);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

/** creates the handler for XYZ classifier and includes it in SCIP */
SCIP_RETCODE GCGincludeConsClassifierScipConstypes(
   GCG*                 gcg                /**< GCG data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata = NULL;

   SCIP_CALL(
      GCGincludeConsClassifier(gcg, CLSCONS_NAME, CLSCONS_DESC, CLSCONS_PRIORITY, CLSCONS_ENABLED, classifierdata,
                               classifierFree, classifierClassify));

   return SCIP_OKAY;
}
