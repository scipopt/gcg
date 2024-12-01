/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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
#include "def.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the mip solver for pricing problems and includes it in GCG */
GCG_EXPORT
SCIP_RETCODE GCGincludeSolverMip(
   SCIP*                 scip                /**< SCIP data structure */
);

/** get the status of the pricing problem */
GCG_EXPORT
GCG_PRICINGSTATUS getPricingstatus(
   SCIP*                 pricingprob         /**< pricing problem SCIP data structure */
);

/** extracts ray from a subproblem used to solve a pricing problem pricing problem (or directly from the pricing problem if no subproblem is specified) */
GCG_EXPORT
SCIP_RETCODE createColumnFromRay(
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   SCIP*                 subproblem,         /**< SCIP data structure that contains the actual solution (if NULL pricingprob will be used) */
   SCIP_HASHMAP*         varmap,             /**< mapping of pricingprob vars to subproblem vars (can be NULL if subproblem is NULL) */
   int                   probnr,             /**< problem number */
   GCG_COL**             newcol              /**< column pointer to store new column */
);

/** transforms feasible solutions of a subproblem used to solve a pricing problem pricing problem into columns (or directly of the pricing problem if no subproblem is specified) */
GCG_EXPORT
SCIP_RETCODE getColumnsFromPricingprob(
   SCIP*                 scip,               /**< master problem SCIP data structure */
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   SCIP*                 subproblem,         /**< SCIP data structure that contains the actual solution (if NULL pricingprob will be used) */
   SCIP_HASHMAP*         varmap,             /**< mapping of pricingprob vars to subproblem vars (can be NULL if subproblem is NULL) */
   int                   probnr,             /**< problem number */
   SCIP_Bool             checksols           /**< should solutions be checked extensively */
);

#ifdef __cplusplus
}
#endif


#endif
