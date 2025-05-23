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

/**@file   clscons.c
 * @ingroup DECOMP
 * @brief  interface for constraint classifier
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"

#include "gcg/struct_consclassifier.h"
#include "gcg/pub_clscons.h"


/** gets user data of constraint classifier */
GCG_CLASSIFIERDATA* GCGconsClassifierGetData(
   GCG_CONSCLASSIFIER* consClassifier  /**< constraint classifier data structure */
   )
{
   assert(consClassifier != NULL);

   return consClassifier->clsdata;
}

/** sets user data of constraint classifier; user has to free old data in advance! */
void GCGconsClassifierSetData(
   GCG_CONSCLASSIFIER* consClassifier,              /**< constraint classifier data structure */
   GCG_CLASSIFIERDATA* clsdata                      /**< new constraint classifier user data */
   )
{
   assert(consClassifier != NULL);

   consClassifier->clsdata = clsdata;
}

/** gets name of constraint classifier */
const char* GCGconsClassifierGetName(
   GCG_CONSCLASSIFIER* consClassifier  /**< constraint classifier data structure */
   )
{
   assert(consClassifier != NULL);

   return consClassifier->name;
}

/** gets priority of constraint classifier */
int GCGconsClassifierGetPriority(
   GCG_CONSCLASSIFIER* consClassifier  /**< constraint classifier data structure */
   )
{
   assert(consClassifier != NULL);

   return consClassifier->priority;
}

/** returns True iff constraint classifier is enabled */
SCIP_Bool GCGconsClassifierIsEnabled(
   GCG_CONSCLASSIFIER* consClassifier  /**< constraint classifier data structure */
   )
{
   assert(consClassifier != NULL);

   return consClassifier->enabled;
}

/** gets description of constraint classifier */
const char* GCGconsClassifierGetDesc(
   GCG_CONSCLASSIFIER* consClassifier  /**< constraint classifier data structure */
   )
{
   assert(consClassifier != NULL);

   return consClassifier->description;
}
