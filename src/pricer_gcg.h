/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   pricer_gcg.h
 * @brief  gcg variable pricer
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRICER_GCG__
#define __SCIP_PRICER_GCG__

#include "scip/scip.h"
#include "relax_gcg.h"
#include "type_solver.h"

enum GCG_Pricetype
{
   GCG_PRICETYPE_INIT      = 0,       /**< initial pricing */
   GCG_PRICETYPE_FARKAS    = 1,       /**< farkas pricing */
   GCG_PRICETYPE_REDCOST   = 2        /**< redcost pricing */
};
typedef enum GCG_Pricetype GCG_PRICETYPE;



/** creates the gcg variable pricer and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePricerGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origprob            /**< SCIP data structure of the original problem */
   );


/* informs an original variable, that a variable in the master problem was created, that contains a part of the original variable,
 * saves this information int the original variable's data */
extern
SCIP_RETCODE GCGpricerAddMasterVarToOrigVar(
   SCIP*                 scip,
   SCIP_VAR*             origvar,
   SCIP_VAR*             var,
   SCIP_Real             val
   );

/** returns the pointer to the scip instance representing the original problem */
extern
SCIP* GCGpricerGetOrigprob(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array of variables that were priced in during the solving process */
extern
SCIP_VAR** GCGpricerGetPricedvars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of variables that were priced in during the solving process */
extern
int GCGpricerGetNPricedvars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds the given constraint and the given position to the hashmap of the pricer */
extern
SCIP_RETCODE GCGpricerAddMasterconsToHashmap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint that should be added */
   int                   pos                 /**< the position of the constraint in the relaxator's masterconss array */
   );

/** includes a solver into the pricer data */
extern
SCIP_RETCODE GCGpricerIncludeSolver(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,
   const char*           description,
   int                   priority,
   GCG_DECL_SOLVERSOLVE  ((*solversolve)),   /**<  solving method for solver */
   GCG_DECL_SOLVERSOLVEHEUR((*solveheur)),   /**<  heuristic solving method for solver */
   GCG_DECL_SOLVERFREE   ((*solverfree)),
   GCG_DECL_SOLVERINIT   ((*solverinit)),
   GCG_DECL_SOLVEREXIT   ((*solverexit)),
   GCG_DECL_SOLVERINITSOL((*solverinitsol)),
   GCG_DECL_SOLVEREXITSOL((*solverexitsol)),
   GCG_SOLVERDATA*       solverdata
   );

extern
GCG_SOLVERDATA* GCGpricerGetSolverdata(
   SCIP*                 scip,
   GCG_SOLVER*           solver
   );

extern
void GCGpricerSetSolverdata(
   SCIP*                 scip,
   GCG_SOLVER*           solver,
   GCG_SOLVERDATA*       solverdata
   );

extern
void GCGpricerPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** transfers a primal solution of the original problem into the master variable space,
 *  i.e. creates one master variable for each block and adds the solution to the master problem  */
extern
SCIP_RETCODE GCGpricerTransOrigSolToMasterVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             origsol             /**< the solution that should be transferred */
   );

#endif
