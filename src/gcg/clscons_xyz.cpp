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

/**@file    clscons_xyz.cpp
 * 
 * @brief   xyz constraint classifier (put your description here)
 * @author
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/clscons_xyz.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include <vector>
#include <stdio.h>
#include <sstream>

#include "gcg/class_detprobdata.h"

#include "gcg/class_conspartition.h"
#include "gcg/scip_misc.h"

/* classifier properties */
#define CLSCONS_NAME                  "xyz constraint classifier"           /**< name of classifier */
#define CLSCONS_DESC                  "constraint classifier template"      /**< short description of classification */
#define CLSCONS_PRIORITY              0                                     /**< priority of this classifier */

#define CLSCONS_ENABLED               TRUE



/*
 * Data structures
 */

/** @todo fill in the necessary classifier data */

/** classifier handler data */
struct CLSCONS_ClassifierData
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
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_FREECLASSIFIER(classifierFreeXyz)
{  /*lint --e{715}*/

   SCIPerrorMessage("Free function of classifier <%s> not implemented!\n", CLSCONS_NAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define classifierFree NULL
#endif

static
GCG_DECL_CONSCLASSIFY(classifierClassify)
{
   gcg::DETPROBDATA* detprobdata;
   if( transformed )
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataPresolved(scip);
   }
   else
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataOrig(scip);
   }

   // CLASSIFICATION
   gcg::ConsPartition* partition;
   // TODO initialize partition

   detprobdata->addConsPartition(partition);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

/** creates the handler for XYZ classifier and includes it in SCIP */
SCIP_RETCODE GCGincludeConsClassifierXyz(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata;

   /**@todo create xyz classifier data here*/
   classifierdata = NULL;

   SCIP_CALL(
      GCGincludeConsClassifier(gcg, CLSCONS_NAME, CLSCONS_DESC, CLSCONS_PRIORITY, CLSCONS_ENABLED, classifierdata,
         classifierFree, classifierClassify));

   return SCIP_OKAY;
}
