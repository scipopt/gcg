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

/**@file   pub_clsvar.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for variable classifiers
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_PUB_CLSVAR_H__
#define __GCG_PUB_CLSVAR_H__



#include "gcg/type_varclassifier.h"
#include "gcg/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicVarClassifierMethods
 *
 * @{
 */

/** gets user data of variable classifier */
GCG_EXPORT
GCG_CLASSIFIERDATA* GCGvarClassifierGetData(
   GCG_VARCLASSIFIER* varClassifier                 /**< variable classifier data structure */
   );

/** sets user data of variable classifier; user has to free old data in advance! */
GCG_EXPORT
void GCGvarClassifierSetData(
   GCG_VARCLASSIFIER* varClassifier,                /**< variable classifier data structure */
   GCG_CLASSIFIERDATA* clsdata                      /**< new variable classifier user data */
   );

/** gets name of variable classifier */
GCG_EXPORT
const char* GCGvarClassifierGetName(
   GCG_VARCLASSIFIER* varClassifier                 /**< variable classifier data structure */
   );

/** gets priority of variable classifier */
GCG_EXPORT
int GCGvarClassifierGetPriority(
   GCG_VARCLASSIFIER* varClassifier                 /**< variable classifier data structure */
   );

/** returns True iff variable classifier is enabled */
GCG_EXPORT
SCIP_Bool GCGvarClassifierIsEnabled(
   GCG_VARCLASSIFIER* varClassifier                 /**< variable classifier data structure */
   );

/** gets description of variable classifier */
GCG_EXPORT
const char* GCGvarClassifierGetDesc(
   GCG_VARCLASSIFIER* varClassifier                 /**< variable classifier data structure */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
