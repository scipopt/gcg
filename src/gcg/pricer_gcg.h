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

/**@file   pricer_gcg.h
 * @brief  GCG variable pricer
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @ingroup PRICERS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_PRICER_GCG__
#define GCG_PRICER_GCG__

#include "gcg/gcgvarhistory.h"
#include "scip/scip.h"
#include "gcg/gcg.h"
#include "gcg/type_solver.h"
#include "gcg/type_colpool.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@defgroup GCGPRICER GCG Variable Pricer
 * @ingroup PRICING_PUB
 * @{
 */

/** creates the GCG variable pricer and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludePricerGcg(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the array of variables that were priced in during the solving process */
GCG_EXPORT
SCIP_VAR** GCGmasterGetPricedvars(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the number of variables that were priced in during the solving process */
GCG_EXPORT
int GCGmasterGetNPricedvars(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** adds the given constraint and the given position to the hashmap of the pricer */
GCG_EXPORT
SCIP_RETCODE GCGmasterAddMasterconsToHashmap(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS*            cons,               /**< the constraint that should be added */
   int                   pos                 /**< the position of the constraint in the relaxator's masterconss array */
   );

/** sets the optimal LP solution in the pricerdata */
GCG_EXPORT
SCIP_RETCODE GCGmasterSetRootLPSol(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_SOL**            sol                 /**< pointer to optimal solution to root LP */
   );

#ifdef SCIP_STATISTIC
/** gets the optimal LP solution in the pricerdata */
GCG_EXPORT
SCIP_SOL* GCGmasterGetRootLPSol(
   GCG*                  gcg                 /**< GCG data structure */
   );
#endif

/** includes a solver into the pricer data */
GCG_EXPORT
SCIP_RETCODE GCGpricerIncludeSolver(
   GCG*                  gcg,                /**< GCG data structure */
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


/** returns the available pricing solvers */
GCG_EXPORT
GCG_SOLVER** GCGpricerGetSolvers(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the number of available pricing solvers */
GCG_EXPORT
int GCGpricerGetNSolvers(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** writes out a list of all pricing problem solvers */
GCG_EXPORT
void GCGpricerPrintListOfSolvers(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** prints pricing solver statistics */
GCG_EXPORT
void GCGpricerPrintPricingStatistics(
   GCG*                  gcg,                /**< GCG data structure */
   FILE*                 file                /**< output file */
   );

GCG_EXPORT
void GCGpricerPrintStatistics(
   GCG*                  gcg,                /**< GCG data structure */
   FILE*                 file                /**< output file */
   );

/** method to get existence of rays */
GCG_EXPORT
SCIP_RETCODE GCGpricerExistRays(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_Bool*            exist               /**< pointer to store if there exists any ray */
   );

/** get the number of extreme points that a pricing problem has generated so far */
GCG_EXPORT
int GCGpricerGetNPointsProb(
   GCG*                  gcg,                /**< GCG data structure */
   int                   probnr              /**< index of pricing problem */
   );

/** get the number of extreme rays that a pricing problem has generated so far */
GCG_EXPORT
int GCGpricerGetNRaysProb(
   GCG*                  gcg,                /**< GCG data structure */
   int                   probnr              /**< index of pricing problem */
   );

/** get the number of columns to be added to the master LP in the current pricing round */
GCG_EXPORT
int GCGpricerGetMaxColsRound(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** get the number of columns per pricing problem to be added to the master LP in the current pricing round */
GCG_EXPORT
int GCGpricerGetMaxColsProb(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** add a new column to the pricing storage */
GCG_EXPORT
SCIP_RETCODE GCGpricerAddCol(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_COL*              col                 /**< priced col */
   );

/** add a new column to the pricing storage and store result*/
GCG_EXPORT
SCIP_RETCODE GCGpricerAddColResult(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_COL*              col,                /**< priced col */
   SCIP_Bool*            added               /**< pointer to var that indicates whether the col was added */
   );

/** transfers a primal solution of the original problem into the master variable space,
 *  i.e. creates one master variable for each block and adds the solution to the master problem  */
GCG_EXPORT
SCIP_RETCODE GCGmasterTransOrigSolToMasterVars(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_SOL*             origsol,            /**< the solution that should be transferred */
   SCIP_Bool*            stored              /**< pointer to store if transferred solution is feasible (or NULL) */
   );

/** create initial master variables */
GCG_EXPORT
SCIP_RETCODE GCGmasterCreateInitialMastervars(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** get root node degeneracy */
GCG_EXPORT
SCIP_Real GCGmasterGetDegeneracy(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** return if artifical variables are used in current solution */
GCG_EXPORT
SCIP_Bool GCGmasterIsCurrentSolValid(
   GCG*                  gcg                 /**< GCG data structure */
   );

GCG_EXPORT
SCIP_Bool GCGmasterIsBestsolValid(
   GCG*                  gcg                 /**< GCG data structure */
   );

GCG_EXPORT
SCIP_Bool GCGmasterIsSolValid(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_SOL*             mastersol           /**< solution of the master problem, or NULL for current LP solution */
   );

/** get number of iterations in pricing problems */
GCG_EXPORT
SCIP_Longint GCGmasterGetPricingSimplexIters(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** print simplex iteration statistics */
GCG_EXPORT
SCIP_RETCODE GCGmasterPrintSimplexIters(
   GCG*                  gcg,                /**< GCG data structure */
   FILE*                 file                /**< output file */
   );

/** set pricing objectives */
GCG_EXPORT
SCIP_RETCODE GCGsetPricingObjs(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_Real*            dualsolconv         /**< array of dual solutions corresponding to convexity constraints */
   );

/** sets the dual weight for the pricing objective */
extern
void GCGsetPricingObjDualWeight(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_Real             dualweight          /**< the weighting applied to the dual variables */
   );

/** sets the Lagrangian relaxation weight for the pricing objective */
extern
void GCGsetPricingObjRelaxWeight(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_Real*            weights,            /**< the Lagrangian relaxation weight applied to the master constraints */
   int*                  weightids,          /**< the constraint ids for the relaxation weights */
   int                   nweights            /**< the number of weights */
   );

/** creates a new master variable corresponding to the given gcg column */
GCG_EXPORT
SCIP_RETCODE GCGcreateNewMasterVarFromGcgCol(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_Bool             infarkas,           /**< in Farkas pricing? */
   GCG_COL*              gcgcol,             /**< GCG column data structure */
   SCIP_Bool             force,              /**< should the given variable be added also if it has non-negative reduced cost? */
   SCIP_Bool*            added,              /**< pointer to store whether the variable was successfully added */
   SCIP_VAR**            addedvar,           /**< pointer to store the created variable */
   SCIP_Real             score               /**< score of column (or -1.0 if not specified) */
   );

/** computes the reduced cost of a column */
GCG_EXPORT
SCIP_Real GCGcomputeRedCostGcgCol(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_Bool             infarkas,           /**< in Farkas pricing? */
   GCG_COL*              gcgcol,             /**< gcg column to compute reduced cost for */
   SCIP_Real*            objvalptr           /**< pointer to store the computed objective value */
   );

/** compute master and cut coefficients of column */
GCG_EXPORT
SCIP_RETCODE GCGcomputeColMastercoefs(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_COL*              gcgcol              /**< GCG column data structure */
   );

/** includes a pricing callback plugin into the pricer data */
extern
SCIP_RETCODE GCGpricerIncludePricingcb(
   GCG*                  gcg,                /**< GCG data structure */
   const char*           name,               /**< name of pricing callback plugin */
   const char*           desc,               /**< description of pricing callback plugin */
   int                   priority,           /**< priority of pricing callback plugin */
   GCG_DECL_PRICINGCBFREE((*pricingcbfree)), /**< destructor of the pricing callback */
   GCG_DECL_PRICINGCBINIT((*pricingcbinit)), /**< initialize the pricing callback */
   GCG_DECL_PRICINGCBEXIT((*pricingcbexit)), /**< deinitialize the pricing callback */
   GCG_DECL_PRICINGCBINITSOL((*pricingcbinitsol)),/**< solving process initialization method of the pricing callback */
   GCG_DECL_PRICINGCBEXITSOL((*pricingcbexitsol)),/**< solving process deinitialization method of the pricing callback */
   GCG_DECL_PRICINGCBPREPRICING((*pricingcbprepricing)),/**< pre-pricing method of the pricing callback */
   GCG_DECL_PRICINGCBPOSTPRICING((*pricingcbpostpricing)),/**< post-pricing method of the pricing callback */
   GCG_PRICINGCBDATA*  pricingcbdata         /**< pricing callback data */
   );

/** returns the available pricing callback plugins */
extern
GCG_PRICINGCB** GCGpricerGetPricingcbs(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the number of available pricing callback plugins */
extern
int GCGpricerGetNPricingcbs(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the pricing callback plugin of the given name, or NULL if it doesn't exist */
extern
GCG_PRICINGCB* GCGpricerFindPricingcb(
   GCG*                  gcg,                /**< GCG data structure */
   const char*           name                /**< the name of the pricing callback plugin */
   );

/** get colpool */
GCG_EXPORT
GCG_COLPOOL* GCGgetColpool(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** get a weak reference to the current and latest varhistory pointer */
GCG_VARHISTORY* GCGgetCurrentVarhistoryReference(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef _OPENMP
/** get maximal number of pricing threads */
GCG_EXPORT
int GCGpricerGetMaxNThreads(
   GCG*                  gcg                 /**< GCG data structure */
   );
#endif

/** returns the GCG data structure */
GCG_EXPORT
GCG* GCGpricerGetGcg(
   SCIP*                 masterprob          /**< SCIP data structure */
   );

/**@} */
#ifdef __cplusplus
}

#endif

#endif
