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
