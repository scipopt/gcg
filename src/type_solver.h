/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   type_solver.h
 * @brief  type definitions for pricing problem solvers in gcg projects
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SOLVER_H__
#define __SCIP_TYPE_SOLVER_H__

#include "scip/def.h"
#include "scip/type_result.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GCG_SolverData GCG_SOLVERDATA;   /**< solver data */
typedef struct GCG_Solver GCG_SOLVER;           /**< the solver */


/** destructor of pricing solver to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - solver          :  the pricing solver itself
 */
#define GCG_DECL_SOLVERFREE(x) SCIP_RETCODE x (SCIP* scip, GCG_SOLVER* solver)

/** initialization method of pricing solver (called after problem was transformed and solver is active)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - solver          : the pricing solver itself
 */
#define GCG_DECL_SOLVERINIT(x) SCIP_RETCODE x (SCIP* scip, GCG_SOLVER* solver)

/** deinitialization method of pricing solver (called before transformed problem is freed and solver is active)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - solver          : the pricing solver itself
 */
#define GCG_DECL_SOLVEREXIT(x) SCIP_RETCODE x (SCIP* scip, GCG_SOLVER* solver)

/** solving process initialization method of pricing solver (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The pricing solver may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - solver          : the pricing solver itself
 */
#define GCG_DECL_SOLVERINITSOL(x) SCIP_RETCODE x (SCIP* scip, GCG_SOLVER* solver)

/** solving process deinitialization method of pricing solver (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The pricing solver should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - solver          : the pricing solver itself
 */
#define GCG_DECL_SOLVEREXITSOL(x) SCIP_RETCODE x (SCIP* scip, GCG_SOLVER* solver)


/** solving method for pricing solver
 *  
 *
 *  input:
 *  - scip            : SCIP main data structure of the master problem
 *  - solver          : the solver itself
 *  - pricingprob     : the pricing problem that should be solved
 *  - probnr          : number of the pricing problem
 *  - solvars         : pointer to store array with variables for each solution
 *  - solvals         : pointer to store array with values of the variables in the solutions
 *  - nsolvars        : pointer to store array with number of variables in the solutions
 *  - nsols           : pointer to store number of solutions
 *  - result          : the result of the solving call: 
 *                      - SCIP_SUCCESS if problem was solved to optimality
 *                      - SCIP_DIDNOTRUN if not
 */
#define GCG_DECL_SOLVERSOLVE(x) SCIP_RETCODE x (SCIP* scip, GCG_SOLVER* solver, SCIP* pricingprob, int probnr, SCIP_VAR**** solvars, SCIP_Real*** solvals, int** nsolvars, int* nsols, SCIP_STATUS* result)

/** solving method for pricing solver
 *  
 *
 *  input:
 *  - scip            : SCIP main data structure of the master problem
 *  - solver          : the solver itself
 *  - pricingprob     : the pricing problem that should be solved
 *  - probnr          : number of the pricing problem
 *  - solvars         : pointer to store array with variables for each solution
 *  - solvals         : pointer to store array with values of the variables in the solutions
 *  - nsolvars        : pointer to store array with number of variables in the solutions
 *  - nsols           : pointer to store number of solutions
 *  - result          : the result of the solving call: 
 *                      - SCIP_SUCCESS if problem was solved to optimality
 *                      - SCIP_DIDNOTRUN if not
 */
#define GCG_DECL_SOLVERSOLVEHEUR(x) SCIP_RETCODE x (SCIP* scip, GCG_SOLVER* solver, SCIP* pricingprob, int probnr, SCIP_VAR**** solvars, SCIP_Real*** solvals, int** nsolvars, int* nsols, SCIP_STATUS* result)


#ifdef __cplusplus
}
#endif

#endif
