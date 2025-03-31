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

/**@file    clsvar_gamsdomain.cpp
 * 
 * @brief   Classifies by domains from which variables are created TODO: what is together in one class?
 * @author  Stefanie Koß
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/clsvar_gamsdomain.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include <map>
#include <set>
#include <stdio.h>
#include <sstream>
#include <vector>

#include "gcg/class_detprobdata.h"

#include "gcg/class_varpartition.h"
#include "gcg/scip_misc.h"

/* classifier properties */
#define CLSVAR_NAME                  "gamsdomain"              /**< name of classifier */
#define CLSVAR_DESC                  "domain in gams file"     /**< short description of classification */
#define CLSVAR_PRIORITY              0                         /**< priority of this classifier */

#define CLSVAR_ENABLED               TRUE


/*
 * Data structures
 */

/** classifier handler data */
struct GCG_ClassifierData
{
   std::map<std::string, std::set<int>>*      vartodomain;            /**< maps variable name to the corresponding set of domain indices */
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
GCG_DECL_FREEVARCLASSIFIER(classifierFree)
{
   GCG_CLASSIFIERDATA* classifierdata = GCGvarClassifierGetData(classifier);
   assert(classifierdata != NULL);
   assert(strcmp(GCGvarClassifierGetName(classifier), CLSVAR_NAME) == 0);

   delete classifierdata->vartodomain;

   SCIPfreeMemory(GCGgetOrigprob(gcg), &classifierdata);

   return SCIP_OKAY;
}

static
GCG_DECL_VARCLASSIFY(classifierClassify)
{
   gcg::DETPROBDATA* detprobdata;

   if( transformed )
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataPresolved(gcg);
   }
   else
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataOrig(gcg);
   }

   int nvar = detprobdata->getNVars();
   std::vector<int> nvarsForClass( 0 );                // [i] holds number of variables for class i
   std::vector<std::set<int>> domainForClass( 0 );     // [i] holds domain for class i
   std::vector<int> classForVar( nvar, - 1 );          // [i] holds class index for variable i -> indexing over detection internal variable array!
   int counterClasses = 0;

   GCG_VARCLASSIFIER* classifier = GCGfindVarClassifier(gcg, CLSVAR_NAME);
   assert(classifier != NULL);

   GCG_CLASSIFIERDATA* classdata = GCGvarClassifierGetData(classifier);
   assert(classdata != NULL);

   /* firstly, assign all variables to classindices */
   // iterate over variables in detection and lookup in classdata->vartodomain
   // iterating over classdata->vartosymbol and lookup variables with getIndexForVar fails with assertion if variable is not found -> should return error value?
   for( int varid = 0; varid < detprobdata->getNVars(); ++ varid )
   {
      SCIP_VAR* var = detprobdata->getVar(varid);
      std::string varname = std::string( SCIPvarGetName( var ) );
      auto domainiter = classdata->vartodomain->find(varname);
      std::set<int> domain;
      if( domainiter != classdata->vartodomain->end() )
      {
         domain = domainiter->second;
      }
      else
      {
         domain = std::set<int>();
         domain.insert(-1);
      }
      
      bool classfound = false;

      /* check if class for domain exists */
      for( size_t classid = 0; classid < domainForClass.size(); ++classid )
      {
         if( domain == domainForClass[classid] )
         {
            classfound = true;
            classForVar[varid] = (int) classid;
            ++nvarsForClass[classid];
            break;
         }
      }

      /* if not, create a new class index */
      if( !classfound )
      {
         classForVar[varid] = counterClasses;
         ++counterClasses;
         domainForClass.push_back( domain );
         nvarsForClass.push_back( 1 );
      }
   }
   assert( counterClasses == (int) domainForClass.size() );

   /* secondly, use these information to create a ConsPartition */
   gcg::VarPartition* partition = new gcg::VarPartition(gcg, "gamsdomain", counterClasses, detprobdata->getNVars() );

   /* set class names and descriptions of every class */
   for( int c = 0; c < partition->getNClasses(); ++ c )
   {
      std::stringstream text;
      text << "{";
      for( auto iter : domainForClass[c] )
      {
         text << iter;

      }
      text << "}";
      partition->setClassName( c, text.str().c_str() );
      text.str( "" );
      text.clear();
      //text << "This class contains all variables with gams domain" << domainForClass[c] << ".";
      partition->setClassDescription( c, text.str().c_str() );
   }

   /* copy the constraint assignment information found in first step */
   for( int i = 0; i < partition->getNVars(); ++ i )
   {
      partition->assignVarToClass( i, classForVar[i] );
   }
   SCIPverbMessage(GCGgetOrigprob(gcg), SCIP_VERBLEVEL_HIGH, NULL, " Varclassifier \"%s\" yields a classification with %d  different variable classes \n", partition->getName(), partition->getNClasses() );

   detprobdata->addVarPartition(partition);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

/** adds an entry to clsdata->vartodomain */
SCIP_RETCODE GCGvarClassifierGamsdomainAddEntry(
   GCG_VARCLASSIFIER*   classifier,
   SCIP_VAR*            var,
   int                  symDomIdx[],
   int*                 symDim
   )
{
   assert(classifier != NULL);
   GCG_CLASSIFIERDATA* classdata = GCGvarClassifierGetData(classifier);
   assert(classdata != NULL);

   std::string varname = SCIPvarGetName( var );
   char varnametrans[SCIP_MAXSTRLEN];
   (void) SCIPsnprintf(varnametrans, SCIP_MAXSTRLEN, "t_%s", varname.c_str());
   std::string nametrans(varnametrans);
   std::set<int> domainset;
   for( int i = 0; i < *symDim; ++i)
   {
      domainset.insert(symDomIdx[i]);
   }
   classdata->vartodomain->insert({varname, domainset});
   classdata->vartodomain->insert({varnametrans, domainset});

   return SCIP_OKAY;
}

/** creates the handler for gamsdomain classifier and includes it in SCIP */
SCIP_RETCODE GCGincludeVarClassifierGamsdomain(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata = NULL;

   SCIP_CALL( SCIPallocMemory(GCGgetOrigprob(gcg), &classifierdata) );
   assert(classifierdata != NULL);
   classifierdata->vartodomain = new std::map<std::string, std::set<int>>();

   SCIP_CALL( GCGincludeVarClassifier(gcg, CLSVAR_NAME, CLSVAR_DESC, CLSVAR_PRIORITY, CLSVAR_ENABLED, classifierdata, classifierFree, classifierClassify) );

   return SCIP_OKAY;
}
