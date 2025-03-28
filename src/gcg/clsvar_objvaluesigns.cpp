/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
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

/**@file   clsvar_objvaluesigns.cpp
 * 
 * @brief classifies variables according to their objective function value signs
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/clsvar_objvaluesigns.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include <vector>
#include <stdio.h>
#include <sstream>

#include "gcg/class_detprobdata.h"

#include "gcg/class_varpartition.h"
#include "gcg/scip_misc.h"

/* classifier properties */
#define CLSVAR_NAME        "objectivevaluesigns"       /**< name of classifier */
#define CLSVAR_DESC                  "objective function value signs"     /**< short description of classification*/
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
   gcg::VarPartition* classifier= new gcg::VarPartition(gcg, "varobjvalsigns", 3, detprobdata->getNVars() ); /* new VarPartition */
   SCIP_Real curobjval;

   /* set up class information */
   classifier->setClassName( 0, "zero" );
   classifier->setClassDescription( 0, "This class contains all variables with objective function value zero." );
   classifier->setClassDecompInfo( 0, gcg::MASTER );
   classifier->setClassName( 1, "positive" );
   classifier->setClassDescription( 1, "This class contains all variables with positive objective function value." );
   classifier->setClassDecompInfo( 1, gcg::ALL );
   classifier->setClassName( 2, "negative" );
   classifier->setClassDescription( 2, "This class contains all variables with negative objective function value." );
   classifier->setClassDecompInfo( 2, gcg::ALL );

   /* assign vars */
   for( int v = 0; v < detprobdata->getNVars(); ++v )
   {
      assert( detprobdata->getVar(v) != NULL );
      curobjval = SCIPvarGetObj(detprobdata->getVar(v));

      if( SCIPisZero(origprob, curobjval) )
      {
         classifier->assignVarToClass(v, 0);
      }
      else if ( SCIPisPositive(origprob, curobjval) )
      {
         classifier->assignVarToClass(v, 1);
      }
      else
      {
         classifier->assignVarToClass(v, 2);
      }
   }

   /* remove a class if there is no variable with the respective sign */
   classifier->removeEmptyClasses();

   SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, " Varclassifier \"%s\" yields a classification with %d different variable classes\n", classifier->getName(), classifier->getNClasses()) ;


   detprobdata->addVarPartition(classifier);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

SCIP_RETCODE GCGincludeVarClassifierObjValueSigns(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata = NULL;

   SCIP_CALL( GCGincludeVarClassifier(gcg, CLSVAR_NAME, CLSVAR_DESC, CLSVAR_PRIORITY, CLSVAR_ENABLED, classifierdata, classifierFree, classifierClassify) );

   return SCIP_OKAY;
}
