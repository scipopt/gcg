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

/**@file    type_sepagcg.h
 * @ingroup TYPEDEFINITIONS
 * @brief   type definitions for separators for master problem in GCG projects
 * @author  Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_TYPE_SEPAGCG_H__
#define GCG_TYPE_SEPAGCG_H__

#include <scip/def.h>
#include <scip/type_lp.h>
#include <scip/type_result.h>
#include <scip/type_scip.h>
#include <scip/type_var.h>

#include "type_gcgcol.h"
#include "type_mastersepacut.h"


#ifdef __cplusplus
extern "C" {
#endif


typedef struct GCG_Sepa GCG_SEPA;             /**< separator for master problem*/


/** method for computing the column coefficient for a cut
 *
 *  input:
 *    gcg             : GCG main data structure
 *    sepa            : the gcg separator itself
 *    cut             : cut which has to be altered
 *    gcgcol          : the column representing a new master variable which should be included in cut
 *    coef            : stores the computed coefficient
 */
#define GCG_DECL_SEPAGETCOLCOEFFICIENTS(x) SCIP_RETCODE x (GCG* gcg, GCG_SEPA* sepa, GCG_EXTENDEDMASTERCONSDATA* cut, GCG_COL* gcgcol, SCIP_Real* coeff)

/** method for computing the column coefficient for a cut
 *
 *  input:
 *    gcg             : GCG main data structure
 *    sepa            : the gcg separator itself
 *    cut             : cut which has to be altered
 *    vars            : the column representing a new master variable which should be included in cut
 *    vals            : stores the computed coefficient
 *    nvars           :
 *    coef            :
 */
#define GCG_DECL_SEPAGETVARCOEFFICIENT(x) SCIP_RETCODE x (GCG* gcg, GCG_SEPA* sepa, GCG_EXTENDEDMASTERCONSDATA* cut, SCIP_VAR** vars, SCIP_Real* vals, int nvars, int probnr, SCIP_Real* coef)


/** method for modifying the objectives of pricing problems to account for master cut
 *
 *  input:
 *    gcg             : GCG main data structure
 *    sepa            : the gcg separator itself
 *    cut             : cut which has to be altered
 *    dual            : dual for objective
 */
#define GCG_DECL_SEPASETOBJECTIVE(x) SCIP_RETCODE x (GCG* gcg, GCG_SEPA* sepa, GCG_EXTENDEDMASTERCONSDATA* cut, SCIP_Real dual)

/** method for modifying the outdated values of a gcg column
 *
 *  input:
 *    gcg             : GCG main data structure
 *    sepa            : the gcg separator itself
 *    cut             : cut which has to be altered
 *    gcgcol          : GCG column
 */
#define GCG_DECL_SEPAADJUSTCOL(x) SCIP_RETCODE x (GCG* gcg, GCG_SEPA* sepa, GCG_EXTENDEDMASTERCONSDATA* cut, GCG_COL* gcgcol)

#ifdef __cplusplus
}
#endif

#endif
