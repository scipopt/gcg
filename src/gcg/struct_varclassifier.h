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

/**@file   struct_varclassifier.h
 * @ingroup DATASTRUCTURES
 * @brief  data structures for variable classifiers
 * @author William Ma
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_VARCLASSIFIER_H__
#define GCG_STRUCT_VARCLASSIFIER_H__

#include "gcg/type_varclassifier.h"


/** detector data structure */
struct GCG_VarClassifier {
   const char*           name;               /**< name of the detector */
   const char*           description;        /**< description of the detector */
   int                   priority;           /**< classifier priority */

   SCIP_Bool             enabled;            /* is enabled by default */

   GCG_CLASSIFIERDATA*   clsdata;            /**< custom data structure of the classifiers */

   GCG_DECL_FREEVARCLASSIFIER((*freeClassifier));                  /**< destructor of detector */
   GCG_DECL_VARCLASSIFY((*classify));            /**< structure detection method of detector */
};


#endif //GCG_STRUCT_VARCLASSIFIER_H__
