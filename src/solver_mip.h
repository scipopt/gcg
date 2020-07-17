/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2020 Operations Research, RWTH Aachen University       */
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

/**@file   solver_mip.h
 * @brief  mip solver for pricing problems
 * @author Gerald Gamrath
 * @ingroup PRICINGSOLVERS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_SOLVER_MIP_H__
#define GCG_SOLVER_MIP_H__

#include "scip/scip.h"
#include "type_pricingstatus.h"
#include "type_gcgcol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the mip solver for pricing problems and includes it in GCG */
extern
SCIP_RETCODE GCGincludeSolverMip(
   SCIP*                 scip                /**< SCIP data structure */
);

extern
GCG_PRICINGSTATUS getPricingstatus(
   SCIP*                 pricingprob         /**< pricing problem SCIP data structure */
);

extern
SCIP_RETCODE createColumnFromRay(
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   SCIP*                 subproblem,
   SCIP_HASHMAP*         varmap,
   int                   probnr,             /**< problem number */
   GCG_COL**             newcol              /**< column pointer to store new column */
);

extern
SCIP_RETCODE getColumnsFromPricingprob(
   SCIP*                 scip,               /**< master problem SCIP data structure */
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   SCIP*                 subproblem,
   SCIP_HASHMAP*         varmap,
   int                   probnr,             /**< problem number */
   SCIP_Bool             checksols          /**< should solutions be checked extensively */
);

#ifdef __cplusplus
}
#endif


#endif
