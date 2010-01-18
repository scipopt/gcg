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


/** solving method for solver
 *  
 *
 *  input:
 *  - scip            : SCIP main data structure of the master problem
 *  - pricingprob     : the pricing problem that should be solved
 *  - probnr          : number of the pricing problem
 *  - result          : the result of the solving call: 
 *                      - SCIP_SUCCESS if problem was solved to optimality
 *                      - SCIP_DIDNOTRUN if not
 */
#define GCG_DECL_SOLVERSOLVE(x) SCIP_RETCODE x (SCIP* scip, SCIP* pricingprob, int probnr, SCIP_STATUS* result)

/** solving method for solver
 *  
 *
 *  input:
 *  - scip            : SCIP main data structure of the master problem
 *  - pricingprob     : the pricing problem that should be solved
 *  - probnr          : number of the pricing problem
 *  - result          : the result of the solving call: 
 *                      - SCIP_SUCCESS if problem was solved to optimality
 *                      - SCIP_DIDNOTRUN if not
 */
#define GCG_DECL_SOLVERSOLVEHEUR(x) SCIP_RETCODE x (SCIP* scip, SCIP* pricingprob, int probnr, SCIP_STATUS* result)


#ifdef __cplusplus
}
#endif

#endif
