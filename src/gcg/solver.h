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

/**@file    solver.h
 * @ingroup PRICING
 * @brief   private methods for GCG pricing solvers
 * @author  Henri Lotze
 * @author  Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_SOLVER_H_
#define GCG_SOLVER_H_


#include "gcg/type_solver.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @ingroup PRICING
 *
 * @{
 */

/** creates a GCG pricing solver */
GCG_EXPORT
SCIP_RETCODE GCGsolverCreate(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER**          solver,             /**< pointer to pricing solver data structure */
   const char*           name,               /**< name of solver */
   const char*           desc,               /**< description of solver */
   int                   priority,           /**< priority of solver */
   SCIP_Bool             heurenabled,        /**< flag to indicate whether heuristic solving method of the solver is enabled */
   SCIP_Bool             exactenabled,        /**< flag to indicate whether exact solving method of the solver is enabled */
   GCG_DECL_SOLVERUPDATE((*solverupdate)),   /**< update method for solver */
   GCG_DECL_SOLVERSOLVE  ((*solversolve)),   /**< solving method for solver */
   GCG_DECL_SOLVERSOLVEHEUR((*solveheur)),   /**< heuristic solving method for solver */
   GCG_DECL_SOLVERFREE   ((*solverfree)),    /**< free method of solver */
   GCG_DECL_SOLVERINIT   ((*solverinit)),    /**< init method of solver */
   GCG_DECL_SOLVEREXIT   ((*solverexit)),    /**< exit method of solver */
   GCG_DECL_SOLVERINITSOL((*solverinitsol)), /**< initsol method of solver */
   GCG_DECL_SOLVEREXITSOL((*solverexitsol)), /**< exitsol method of solver */
   GCG_SOLVERDATA*       solverdata          /**< pricing solver data */
   );

/** calls destructor and frees memory of GCG pricing solver */
GCG_EXPORT
SCIP_RETCODE GCGsolverFree(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER**          solver              /**< pointer to pricing solver data structure */
   );

/** initializes GCG pricing solver */
GCG_EXPORT
SCIP_RETCODE GCGsolverInit(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** calls exit method of GCG pricing solver */
GCG_EXPORT
SCIP_RETCODE GCGsolverExit(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** calls solving process initialization method of GCG pricing solver */
GCG_EXPORT
SCIP_RETCODE GCGsolverInitsol(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** calls solving process deinitialization method of GCG pricing solver */
GCG_EXPORT
SCIP_RETCODE GCGsolverExitsol(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** calls update method of GCG pricing solver */
GCG_EXPORT
SCIP_RETCODE GCGsolverUpdate(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP*                 pricingprob,        /**< the pricing problem that should be solved */
   GCG_SOLVER*           solver,             /**< pricing solver */
   int                   probnr,             /**< number of the pricing problem */
   SCIP_Bool             varobjschanged,     /**< have the objective coefficients changed? */
   SCIP_Bool             varbndschanged,     /**< have the lower and upper bounds changed? */
   SCIP_Bool             consschanged        /**< have the constraints changed? */
   );

/** calls heuristic or exact solving method of GCG pricing solver */
GCG_EXPORT
SCIP_RETCODE GCGsolverSolve(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP*                 pricingprob,        /**< the pricing problem that should be solved */
   GCG_SOLVER*           solver,             /**< pricing solver */
   SCIP_Bool             redcost,            /**< is reduced cost (TRUE) or Farkas (FALSE) pricing performed? */
   SCIP_Bool             heuristic,          /**< shall the pricing problem be solved heuristically? */
   int                   probnr,             /**< number of the pricing problem */
   SCIP_Real             dualsolconv,        /**< dual solution of the corresponding convexity constraint */
   SCIP_Real*            lowerbound,         /**< pointer to store lower bound of pricing problem */
   GCG_PRICINGSTATUS*    status,             /**< pointer to store the returned pricing status */
   SCIP_Bool*            solved              /**< pointer to store whether the solution method was called */
   );
/**@} */
#ifdef __cplusplus
}
#endif

#endif
