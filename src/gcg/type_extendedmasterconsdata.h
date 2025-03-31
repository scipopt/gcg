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

/**@file    type_extendedmasterconsdata.h
 * @ingroup TYPEDEFINITIONS
 * @brief   type definitions for extended master conss in GCG projects
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_TYPE_EXTENDEDMASTERCONSDATA_H_
#define GCG_TYPE_EXTENDEDMASTERCONSDATA_H_

#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/** type of extended master constraint */
enum GCG_ExtendedMasterConsType
{
   GCG_EXTENDEDMASTERCONSTYPE_CONS,                   /**< extended master cons is represented by a constraint */
   GCG_EXTENDEDMASTERCONSTYPE_ROW                     /**< extended master cons is represented by a row */
};
typedef enum GCG_ExtendedMasterConsType GCG_EXTENDEDMASTERCONSTYPE;

typedef union GCG_ExtendedMasterCons GCG_EXTENDEDMASTERCONS;

typedef struct GCG_PricingModification GCG_PRICINGMODIFICATION;
typedef struct GCG_ExtendedMasterConsData GCG_EXTENDEDMASTERCONSDATA;

/** determine the coefficient of a column solution in the extended master cons
 *
 *  input:
 *    gcg             : GCG data structure
 *    extendedmasterconsdata   : the extended master cons data
 *    solvars         : array of column solution variables
 *    solvals         : array of column solution values
 *    nsolvars        : number of column solution variables and values
 *    probnr          : the pricing problem that the column belongs to
 *    coef            : the calculated coefficient
 */
#define GCG_DECL_EXTENDEDMASTERCONSGETCOEFF(x) SCIP_RETCODE x (GCG* gcg, GCG_EXTENDEDMASTERCONSDATA* extendedmasterconsdata, SCIP_VAR** solvars, SCIP_Real* solvals, int nsolvars, int probnr, SCIP_Real* coef)

#ifdef __cplusplus
}
#endif

#endif
