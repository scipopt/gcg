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

/**@file   solver_mip.h
 * @brief  mip solver for pricing problems
 * @author Gerald Gamrath
 * @ingroup PRICINGSOLVERS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_SOLVER_MIP_H__
#define GCG_SOLVER_MIP_H__

#include "scip/scip.h"
#include "gcg/gcg.h"


#ifdef __cplusplus
extern "C" {
#endif

/** creates the mip solver for pricing problems and includes it in GCG */
GCG_EXPORT
SCIP_RETCODE GCGincludeSolverMip(
   GCG*                  gcg                 /**< GCG data structure */
);

/** get the status of the pricing problem */
GCG_EXPORT
GCG_PRICINGSTATUS getPricingstatus(
   SCIP*                 pricingprob         /**< pricing problem SCIP data structure */
);

/** extracts ray from a subproblem used to solve a pricing problem pricing problem (or directly from the pricing problem if no subproblem is specified) */
GCG_EXPORT
SCIP_RETCODE createColumnFromRay(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   SCIP*                 subproblem,         /**< SCIP data structure that contains the actual solution (if NULL pricingprob will be used) */
   SCIP_HASHMAP*         varmap,             /**< mapping of pricingprob vars to subproblem vars (can be NULL if subproblem is NULL) */
   int                   probnr,             /**< problem number */
   GCG_COL**             newcol              /**< column pointer to store new column */
);

/** transforms feasible solutions of a subproblem used to solve a pricing problem pricing problem into columns (or directly of the pricing problem if no subproblem is specified) */
GCG_EXPORT
SCIP_RETCODE getColumnsFromPricingprob(
   GCG*                  gcg,                /**< GCG data structure */
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
