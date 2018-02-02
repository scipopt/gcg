/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file   pricer_gcg.h
 * @ingroup PUBLICMETHODS
 * @brief  GCG variable pricer
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @ingroup PRICERS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_PRICER_GCG__
#define GCG_PRICER_GCG__

#include "scip/scip.h"
#include "type_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

enum GCG_Pricetype
{
   GCG_PRICETYPE_UNKNOWN = -1,               /**< unknown pricing type */
   GCG_PRICETYPE_INIT = 0,                   /**< initial pricing */
   GCG_PRICETYPE_FARKAS = 1,                 /**< farkas pricing */
   GCG_PRICETYPE_REDCOST = 2                 /**< redcost pricing */
};
typedef enum GCG_Pricetype GCG_PRICETYPE;



/** creates the GCG variable pricer and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePricerGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origprob            /**< SCIP data structure of the original problem */
   );

/** returns the pointer to the scip instance representing the original problem */
extern
SCIP* GCGmasterGetOrigprob(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array of variables that were priced in during the solving process */
extern
SCIP_VAR** GCGmasterGetPricedvars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of variables that were priced in during the solving process */
extern
int GCGmasterGetNPricedvars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds the given constraint and the given position to the hashmap of the pricer */
extern
SCIP_RETCODE GCGmasterAddMasterconsToHashmap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint that should be added */
   int                   pos                 /**< the position of the constraint in the relaxator's masterconss array */
   );

/** includes a solver into the pricer data */
extern
SCIP_RETCODE GCGpricerIncludeSolver(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of solver */
   const char*           description,        /**< description of solver */
   int                   priority,           /**< priority of solver */
   SCIP_Bool             enabled,            /**< flag to indicate whether the solver is enabled */
   GCG_DECL_SOLVERSOLVE  ((*solversolve)),   /**< solving method for solver */
   GCG_DECL_SOLVERSOLVEHEUR((*solveheur)),   /**< heuristic solving method for solver */
   GCG_DECL_SOLVERFREE   ((*solverfree)),    /**< free method of solver */
   GCG_DECL_SOLVERINIT   ((*solverinit)),    /**< init method of solver */
   GCG_DECL_SOLVEREXIT   ((*solverexit)),    /**< exit method of solver */
   GCG_DECL_SOLVERINITSOL((*solverinitsol)), /**< initsol method of solver */
   GCG_DECL_SOLVEREXITSOL((*solverexitsol)), /**< exitsol method of solver */
   GCG_SOLVERDATA*       solverdata          /**< solverdata data structure */
   );


/** returns the solverdata of a solver */
extern
GCG_SOLVERDATA* GCGsolverGetSolverdata(
   GCG_SOLVER*           solver              /**< pointer so solver */
   );

/** sets solver data of specific solver */
extern
void GCGsolverSetSolverdata(
   GCG_SOLVER*           solver,             /**< pointer to solver  */
   GCG_SOLVERDATA*       solverdata          /**< solverdata data structure */
   );

/** writes out a list of all pricing problem solvers */
extern
void GCGpricerPrintListOfSolvers(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** prints pricing solver statistics */
extern
void GCGpricerPrintPricingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

extern
void GCGpricerPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** method to get existence of rays */
extern
SCIP_RETCODE GCGpricerExistRays(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            exist               /**< pointer to store if there exists any ray */
   );


/** transfers a primal solution of the original problem into the master variable space,
 *  i.e. creates one master variable for each block and adds the solution to the master problem  */
extern
SCIP_RETCODE GCGmasterTransOrigSolToMasterVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             origsol,            /**< the solution that should be transferred */
   SCIP_Bool*            stored              /**< pointer to store if transferred solution is feasible (or NULL) */
   );

/** create initial master variables */
extern
SCIP_RETCODE GCGmasterCreateInitialMastervars(
   SCIP*                 scip                /**< master SCIP data structure */
   );

/** get root node degeneracy */
extern
SCIP_Real GCGmasterGetDegeneracy(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get number of iterations in pricing problems */
extern
SCIP_Longint GCGmasterGetPricingSimplexIters(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** print simplex iteration statistics */
extern
SCIP_RETCODE GCGmasterPrintSimplexIters(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** set pricing objectives */
extern
SCIP_RETCODE GCGsetPricingObjs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            dualsolconv         /**< array of dual solutions corresponding to convexity constraints */
   );


#ifdef __cplusplus
}
#endif

#endif
