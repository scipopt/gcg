/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
   GCG_DECL_SOLVERSOLVE  ((*solversolve)),   /**<  activation method for branchrule */
   GCG_SOLVERDATA*       solverdata
   );

#endif
