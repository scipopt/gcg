/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    clscons_xyz.cpp
 * 
 * @brief   xyz constraint classifier (put your description here)
 * @author
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "clscons_xyz.h"
#include "cons_decomp.h"
#include "cons_decomp.hpp"
#include <vector>
#include <stdio.h>
#include <sstream>

#include "class_detprobdata.h"

#include "class_conspartition.h"
#include "scip_misc.h"

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
SCIP_RETCODE SCIPincludeConsClassifierXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata;

   /**@todo create xyz classifier data here*/
   classifierdata = NULL;

   SCIP_CALL(
      GCGincludeConsClassifier(scip, CLSCONS_NAME, CLSCONS_DESC, CLSCONS_PRIORITY, CLSCONS_ENABLED, classifierdata,
         classifierFree, classifierClassify));

   return SCIP_OKAY;
}
