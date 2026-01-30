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

/**@file   clsvar_objvalues.cpp
 * 
 * @brief classifies variables according to their objective function values
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/clsvar_objvalues.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include <vector>
#include <stdio.h>
#include <sstream>
#include <iomanip>

#include "gcg/class_detprobdata.h"

#include "gcg/class_varpartition.h"
#include "gcg/scip_misc.h"

/* classifier properties */
#define CLSVAR_NAME        "objectivevalues"       /**< name of classifier */
#define CLSVAR_DESC                  "objective function values"     /**< short description of classification*/
#define CLSVAR_PRIORITY              0

#define CLSVAR_ENABLED               TRUE


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
GCG_DECL_VARCLASSIFY(classifierClassify)
{
   gcg::DETPROBDATA* detprobdata;
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);
   if( transformed )
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataPresolved(gcg);
   }
   else
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataOrig(gcg);
   }

   // CLASSIFICATION
   std::vector<SCIP_Real> foundobjvals; /* all found objective function values */
   std::vector<int> classforvars(detprobdata->getNVars(), -1); /* vector assigning a class index to each variable */
   int curclassindex; /* stores a var's classindex if the objective value of a var has already been found for another var */
   SCIP_Real curobjval;
   gcg::VarPartition* classifier; /* new VarPartition */

   for( int v = 0; v < detprobdata->getNVars(); ++v )
   {
      assert( detprobdata->getVar(v) != NULL );
      curobjval = SCIPvarGetObj(detprobdata->getVar(v));
      curclassindex = -1;

      /* check whether current objective funtion value already exists */
      for( size_t c = 0; c < foundobjvals.size(); ++c )
      {
         if( SCIPisEQ(origprob, curobjval, foundobjvals[c]) )
         {
            curclassindex = (int) c;
            break;
         }
      }

      /* assign var to class and save objective function value, if it is new */
      if( curclassindex == -1 )
      {
         foundobjvals.push_back(curobjval);
         classforvars[v] = (int) foundobjvals.size() - 1;
      }
      else
      {
         classforvars[v] = curclassindex;
      }
   }

   classifier = new gcg::VarPartition(gcg, "varobjvals", (int) foundobjvals.size(), detprobdata->getNVars());

   /* set up class information */
   for ( int c = 0; c < classifier->getNClasses(); ++c )
   {
      std::stringstream name;
      std::stringstream text;

      name << std::setprecision( 5 ) << foundobjvals[c];
      text << "This class contains all variables with objective function value " << name.str() << ".";

      classifier->setClassName(c, name.str().c_str());
      classifier->setClassDescription(c, text.str().c_str());
   }

   /* assign vars according to classforvars vactor */
   for ( int v = 0; v < classifier->getNVars(); ++v )
   {
      classifier->assignVarToClass(v, classforvars[v]);
   }

   SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, " Varclassifier \"%s\" yields a classification with %d different variable classes\n", classifier->getName(), classifier->getNClasses()) ;

   detprobdata->addVarPartition(classifier);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

SCIP_RETCODE GCGincludeVarClassifierObjValues(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata = NULL;

   SCIP_CALL( GCGincludeVarClassifier(gcg, CLSVAR_NAME, CLSVAR_DESC, CLSVAR_PRIORITY, CLSVAR_ENABLED, classifierdata, classifierFree, classifierClassify) );

   return SCIP_OKAY;
}
