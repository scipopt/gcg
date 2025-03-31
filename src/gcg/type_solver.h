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

/**@file   type_solver.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for pricing problem solvers in GCG project
 * @author Gerald Gamrath
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_TYPE_SOLVER_H__
#define GCG_TYPE_SOLVER_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "gcg/type_gcgcol.h"
#include "gcg/type_pricingstatus.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GCG_SolverData GCG_SOLVERDATA;   /**< solver data */
typedef struct GCG_Solver GCG_SOLVER;           /**< the solver */


/** destructor of pricing solver to free user data (called when GCG is exiting)
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - solver          : the pricing solver itself
 */
#define GCG_DECL_SOLVERFREE(x) SCIP_RETCODE x (GCG* gcg, GCG_SOLVER* solver)

/** initialization method of pricing solver (called after problem was transformed and solver is active)
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - solver          : the pricing solver itself
 */
#define GCG_DECL_SOLVERINIT(x) SCIP_RETCODE x (GCG* gcg, GCG_SOLVER* solver)

/** deinitialization method of pricing solver (called before transformed problem is freed and solver is active)
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - solver          : the pricing solver itself
 */
#define GCG_DECL_SOLVEREXIT(x) SCIP_RETCODE x (GCG* gcg, GCG_SOLVER* solver)

/** solving process initialization method of pricing solver (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The pricing solver may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - solver          : the pricing solver itself
 */
#define GCG_DECL_SOLVERINITSOL(x) SCIP_RETCODE x (GCG* gcg, GCG_SOLVER* solver)

/** solving process deinitialization method of pricing solver (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The pricing solver should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - solver          : the pricing solver itself
 */
#define GCG_DECL_SOLVEREXITSOL(x) SCIP_RETCODE x (GCG* gcg, GCG_SOLVER* solver)

/**
 * update method for pricing solver, used to update solver specific pricing problem data
 *
 * The pricing solver may use this method to update its own representation of the pricing problem,
 * i.e. to apply changes on variable objectives and bounds and to apply branching constraints
 */
#define GCG_DECL_SOLVERUPDATE(x) SCIP_RETCODE x (GCG* gcg, SCIP* pricingprob, GCG_SOLVER* solver, int probnr, SCIP_Bool varobjschanged, SCIP_Bool varbndschanged, SCIP_Bool consschanged)

/** solving method for pricing solver which solves the pricing problem to optimality
 *
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - pricingprob     : the pricing problem that should be solved
 *  - solver          : the pricing solver itself
 *  - probnr          : number of the pricing problem
 *  - dualsolconv     : dual solution of the corresponding convexity constraint
 *  - lowerbound      : pointer to store lower bound of pricing problem
 *  - status          : pointer to store the pricing status
 */
#define GCG_DECL_SOLVERSOLVE(x) SCIP_RETCODE x (GCG* gcg, SCIP* pricingprob, GCG_SOLVER* solver, int probnr, SCIP_Real dualsolconv, SCIP_Real* lowerbound, GCG_PRICINGSTATUS* status)

/** solving method for pricing solver using heuristic pricing only
 *
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - pricingprob     : the pricing problem that should be solved
 *  - solver          : the pricing solver itself
 *  - probnr          : number of the pricing problem
 *  - dualsolconv     : dual solution of the corresponding convexity constraint
 *  - lowerbound      : pointer to store lower bound of pricing problem
 *  - status          : pointer to store the pricing status
 */
#define GCG_DECL_SOLVERSOLVEHEUR(x) SCIP_RETCODE x (GCG* gcg, SCIP* pricingprob, GCG_SOLVER* solver, int probnr, SCIP_Real dualsolconv, SCIP_Real* lowerbound, GCG_PRICINGSTATUS* status)


#ifdef __cplusplus
}
#endif

#endif
