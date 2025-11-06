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

/**@file   struct_mastersepacutdata.h
 * @ingroup DATASTRUCTURES
 * @brief  data structures for GCG separator cuts
 * @author Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_MASTERSEPACUTDATA_H__
#define GCG_STRUCT_MASTERSEPACUTDATA_H__

#include <scip/def.h>

#include "gcg/type_extendedmasterconsdata.h"
#include "gcg/type_mastersepacut.h"
#include "gcg/type_gcgvarhistory.h"
#include "gcg/type_sepagcg.h"

#ifdef __cplusplus
extern "C" {
#endif


/** master separator cut data structure */
struct GCG_Mastersepacut
{
   GCG_SEPA*                        sepa;                   /**< GCG Master Separator which generated the mastercut */
   GCG_VARHISTORY*                  knownvarhistory;        /**< pointer to the history of priced variables */
   GCG_MASTERSEPACUTDATA*           data;                   /**< additional data helpful to compute coefficients */
   int                              nuses;                  /**< number of times this mastersepacut is referenced */
};

#ifdef __cplusplus
}
#endif

#endif
