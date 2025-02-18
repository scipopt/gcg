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

/**@file    type_mastercutdata.h
 * @ingroup TYPEDEFINITIONS
 * @brief   type definitions for master cuts in GCG projects
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_TYPE_MASTERCUTDATA_H_
#define GCG_TYPE_MASTERCUTDATA_H_

#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/** type of master cut */
enum GCG_MasterCutType
{
   GCG_MASTERCUTTYPE_CONS,                   /**< master cut is represented by a constraint */
   GCG_MASTERCUTTYPE_ROW                     /**< master cut is represented by a row */
};
typedef enum GCG_MasterCutType GCG_MASTERCUTTYPE;

typedef union GCG_MasterCutCut GCG_MASTERCUTCUT;

typedef struct GCG_PricingModification GCG_PRICINGMODIFICATION;
typedef struct GCG_MasterCutData GCG_MASTERCUTDATA;

/** determine the coefficient of a column solution in the mastercut
 *
 *  input:
 *    scip            : SCIP main data structure of the original problem
 *    mastercutdata   : the generic mastercut data
 *    solvars         : array of column solution variables
 *    solvals         : array of column solution values
 *    nsolvars        : number of column solution variables and values
 *    probnr          : the pricing problem that the column belongs to
 *    coef            : the calculated coefficient
 */
#define GCG_DECL_MASTERCUTGETCOEFF(x) SCIP_RETCODE x (SCIP* scip, GCG_MASTERCUTDATA* mastercutdata, SCIP_VAR** solvars, SCIP_Real* solvals, int nsolvars, int probnr, SCIP_Real* coef)

#ifdef __cplusplus
}
#endif

#endif
