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

/**@file    struct_sepagcg.h
 * @ingroup DATASTRUCTURES
 * @brief   data structures for separators for master
 * @author  Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_SEPAGCG_H__
#define GCG_STRUCT_SEPAGCG_H__

#include <scip/type_sepa.h>

#include "gcg/type_sepagcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** master separator */
struct GCG_Sepa
{
   SCIP_SEPA*                             separator;                     /**< SCIP separator */
   GCG_DECL_SEPAGETCOLCOEFFICIENT         ((*sepagetcolcoefficient));    /**< compute cut coefficient using gcg column */
   GCG_DECL_SEPASETOBJECTIVE              ((*sepasetobjective));         /**< adapt pricing objectives to consider cut */
   GCG_DECL_SEPAMASTERCUTDELETE           ((*sepamastercutdelete));      /**< delete mastersepacutdata */
};

#ifdef __cplusplus
}
#endif

#endif