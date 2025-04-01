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

/**@file   clsvar_scipvartypes.cpp
 * 
 * @brief classifies variables according to their scip vartypes
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/clsvar_scipvartypes.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include <vector>
#include <stdio.h>
#include <sstream>

#include "gcg/class_detprobdata.h"

#include "gcg/class_varpartition.h"
#include "gcg/scip_misc.h"

/* classifier properties */
#define CLSVAR_NAME        "scipvartype"       /**< name of classifier */
#define CLSVAR_DESC                  "scipvartypes"     /**< short description of classification*/
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
   std::vector < SCIP_VARTYPE > foundVartypes( 0 );
   std::vector<int> classForVars = std::vector<int>( detprobdata->getNVars(), - 1 );
   gcg::VarPartition* classifier;

   SCIP_Bool onlycontsub;
   SCIP_Bool onlybinmaster;

   SCIPgetBoolParam(origprob, "detection/benders/onlycontsubpr", &onlycontsub);
   SCIPgetBoolParam(origprob, "detection/benders/onlybinmaster", &onlybinmaster);

   /* firstly, assign all variables to classindices */
   for( int i = 0; i < detprobdata->getNVars(); ++ i )
   {
      SCIP_VAR* var;
      bool found = false;
      var = detprobdata->getVar(i);
      SCIP_VARTYPE vT = SCIPvarGetType( var );
      size_t vartype;

      if( onlycontsub )
      {
         if ( vT == SCIP_VARTYPE_BINARY )
            vT = SCIP_VARTYPE_INTEGER;
      }

      /* check whether the variable's vartype is new */
      for( vartype = 0; vartype < foundVartypes.size(); ++ vartype )
      {
         if( foundVartypes[vartype] == vT )
         {
            found = true;
            break;
         }
      }
      /* if it is new, create a new class index */
      if( ! found )
      {
         foundVartypes.push_back( vT );
         classForVars[i] = (int) foundVartypes.size() - 1;
      }
      else
         classForVars[i] = (int) vartype;
   }

   /* secondly, use these information to create a VarPartition */
   classifier = new gcg::VarPartition(gcg, "vartypes", (int) foundVartypes.size(), detprobdata->getNVars() );

   /* set class names and descriptions of every class */
   for( int c = 0; c < classifier->getNClasses(); ++ c )
   {
      std::string name;
      std::stringstream text;
      switch( foundVartypes[c] )
      {
         case SCIP_VARTYPE_BINARY:
            name = "bin";
            if( onlybinmaster )
               classifier->setClassDecompInfo(c, gcg::LINKING);
            break;
         case SCIP_VARTYPE_INTEGER:
            name = "int";
            if( onlycontsub )
               classifier->setClassDecompInfo(c, gcg::LINKING);
            if( onlybinmaster )
               classifier->setClassDecompInfo(c, gcg::BLOCK);
            break;
         case SCIP_VARTYPE_CONTINUOUS:
            name = "cont";
            if( onlycontsub )
               classifier->setClassDecompInfo(c, gcg::BLOCK);
            if( onlybinmaster )
               classifier->setClassDecompInfo(c, gcg::BLOCK);
            break;
         default:
            name = "newVartype";
            break;
      }
      classifier->setClassName( c, name.c_str() );
      text << "This class contains all variables that are of (SCIP) vartype \"" << name << "\".";
      classifier->setClassDescription( c, text.str().c_str() );
   }

   /* copy the variable assignment information found in first step */
   for( int i = 0; i < classifier->getNVars(); ++ i )
   {
      classifier->assignVarToClass( i, classForVars[i] );
   }

   SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, " Varclassifier \"%s\" yields a classification with %d different variable classes\n", classifier->getName(), classifier->getNClasses() ) ;

   detprobdata->addVarPartition(classifier);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

SCIP_RETCODE GCGincludeVarClassifierScipVartypes(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata = NULL;

   SCIP_CALL( GCGincludeVarClassifier(gcg, CLSVAR_NAME, CLSVAR_DESC, CLSVAR_PRIORITY, CLSVAR_ENABLED, classifierdata, classifierFree, classifierClassify) );

   return SCIP_OKAY;
}
