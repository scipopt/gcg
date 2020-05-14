/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2019 Operations Research, RWTH Aachen University       */
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
 * @ingroup CLASSIFIERS
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
#define DEC_CLASSIFIERNAME        "xyz variable classifier"           /**< name of classifier */
#define DEC_DESC                  "variable classifier template"      /**< short description of classification */
#define DEC_PRIORITY              0                                   /**< priority of this classifier */

#define DEC_ENABLEDORIG           TRUE                                /**< classify on original problem? */
#define DEC_ENABLEDPRESOLVED      TRUE                                /**< classify on presolved problem? */


/*
 * Data structures
 */


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * classifier callback methods
 */

/** destructor of classifier to free user data (called when GCG is exiting) */
#if 0
static
DEC_DECL_FREECLASSIFIER(classifierFreeXyz)
{  /*lint --e{715}*/

   SCIPerrorMessage("Free function of classifier <%s> not implemented!\n", DEC_CLASSIFIERNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define classifierFree NULL
#endif

/** classifier initialization method (called after problem was transformed) */
#if 0
static
DEC_DECL_INITCLASSIFIER(classifierInitXyz)
{  /*lint --e{715}*/

   SCIPerrorMessage("Init function of classifier <%s> not implemented!\n", DEC_CLASSIFIERNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define classifierInit NULL
#endif

static
DEC_DECL_VARCLASSIFY(classifierClassify)
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
   gcg::PartitionClassifier* partition;
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
   SCIP_CALL( DECincludeVarClassifier(scip, DEC_CLASSIFIERNAME, DEC_DESC, DEC_PRIORITY, DEC_ENABLEDORIG, DEC_ENABLEDPRESOLVED, classifierInit, classifierFree, classifierClassify) );

   return SCIP_OKAY;
}
