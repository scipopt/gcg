/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   type_varclassifier.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for variable classifier in GCG projects
 * @author William Ma
 */

#ifndef GCG_TYPE_VARCLASSIFIER_H__
#define GCG_TYPE_VARCLASSIFIER_H__

#include <scip/def.h>
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"
#include "gcg/type_classifier.h"
#include "gcg/type_gcg.h"

typedef struct GCG_VarClassifier GCG_VARCLASSIFIER;

/** destructor of classifier to free classifier data (called when GCG is exiting)
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - classifier      : classifier data structure
 */
#define GCG_DECL_FREEVARCLASSIFIER(x) SCIP_RETCODE x (GCG* gcg, GCG_VARCLASSIFIER* classifier)

/**
 * Tries to classify variables with data of the according detprobdata and store the classification in the detprobdata
 *
 * input:
 *  - gcg                  : GCG data structure
 *  - classifier           : classifier data structure
 *  - transformed          : should use data from transformed detprobdata or not
 */
#define GCG_DECL_VARCLASSIFY(x) SCIP_RETCODE x (GCG* gcg, GCG_VARCLASSIFIER* varclassifier, SCIP_Bool transformed)

#endif //GCG_TYPE_VARCLASSIFIER_H__
