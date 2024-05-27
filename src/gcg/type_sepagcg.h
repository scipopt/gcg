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

/**@file    type_sepagcg.h
 * @ingroup TYPEDEFINITIONS
 * @brief   type definitions for separators for master problem in GCG projects
 * @author  Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_TYPE_SEPAGCG_H__
#define GCG_TYPE_SEPAGCG_H__

#include "scip/def.h"
#include "scip/type_lp.h"
#include "scip/type_result.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include "struct_mastercutdata.h"
#include "type_gcgcol.h"
#include "event_sepacuts.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct GCG_Sepa GCG_SEPA;             /**< separator for master problem*/


/** method for computing the column coefficient for a cut
 *
 *  input:
 *    scip            : SCIP main data structure (of the master problem)
 *    sepa            : the gcg separator itself
 *    cut             : cut which has to be altered
 *    gcgcol          : the column representing a new master variable which should be included in cut
 *    coef            : stores the computed coefficient
 */
#define GCG_DECL_SEPAGETCOLCOEFFICIENTS(x) SCIP_RETCODE x (SCIP* scip, GCG_SEPA* sepa, GCG_MASTERSEPACUT* cut, GCG_COL* gcgcol, SCIP_Real* coeff)

/** method for computing the column coefficient for a cut
 *
 *  input:
 *    scip            : SCIP main data structure (of the master problem)
 *    sepa            : the gcg separator itself
 *    cut             : cut which has to be altered
 *    vars            : the column representing a new master variable which should be included in cut
 *    vals            : stores the computed coefficient
 *    nvars           :
 *    coef            :
 */
#define GCG_DECL_SEPAGETVARCOEFFICIENT(x) SCIP_RETCODE x (SCIP* scip, GCG_SEPA* sepa, GCG_MASTERCUTDATA* cut, SCIP_VAR** vars, SCIP_Real* vals, int nvars, int probnr, SCIP_Real* coef)


/** method for modifying the objectives of pricing problems to account for master cut
 *
 *  input:
 *    scip            : SCIP main data structure (of the master problem)
 *    sepa            : the gcg separator itself
 *    cut             : cut which has to be altered
 *    dual            : dual for objective
 */
#define GCG_DECL_SEPASETOBJECTIVE(x) SCIP_RETCODE x (SCIP* scip, GCG_SEPA* sepa, GCG_MASTERCUTDATA* cut, SCIP_Real dual)

#ifdef __cplusplus
}
#endif

#endif
