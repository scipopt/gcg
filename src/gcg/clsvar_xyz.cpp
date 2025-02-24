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

/**@file    clsvar_xyz.cpp
 * 
 * @brief   xyz variable classifier (put your description here)
 * @author
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "clsvar_xyz.h"
#include "cons_decomp.h"
#include "cons_decomp.hpp"
#include <vector>
#include <stdio.h>
#include <sstream>

#include "class_detprobdata.h"

#include "class_varpartition.h"
#include "scip_misc.h"

/* classifier properties */
#define CLSVAR_NAME                  "xyz variable classifier"           /**< name of classifier */
#define CLSVAR_DESC                  "variable classifier template"      /**< short description of classification */
#define CLSVAR_PRIORITY              0                                   /**< priority of this classifier */

#define CLSVAR_ENABLED               TRUE


/*
 * Data structures
 */

/** @todo fill in the necessary classifier data */

/** classifier handler data */
struct CLSVAR_ClassifierData
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
GCG_DECL_FREEVARCLASSIFIER(classifierFreeXyz)
{  /*lint --e{715}*/

   SCIPerrorMessage("Free function of classifier <%s> not implemented!\n", CLSVAR_NAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define classifierFree NULL
#endif

static
GCG_DECL_VARCLASSIFY(classifierClassify)
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
   gcg::VarPartition* partition;
   // TODO initialize partition

   detprobdata->addVarPartition(partition);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

/** creates the handler for xyz classifier and includes it in SCIP */
SCIP_RETCODE SCIPincludeVarClassifierXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata;

   /**@todo create xyz classifier data here*/
   classifierdata = NULL;

   SCIP_CALL( GCGincludeVarClassifier(scip, CLSVAR_NAME, CLSVAR_DESC, CLSVAR_PRIORITY, CLSVAR_ENABLED, classifierdata, classifierFree, classifierClassify) );

   return SCIP_OKAY;
}
