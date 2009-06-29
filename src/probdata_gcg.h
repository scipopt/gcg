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

/**@file   probdata_gcg.h
 * @brief  problem data for gcg algorithm
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_GCG__
#define __SCIP_PROBDATA_GCG__

#include "scip/scip.h"
#include "tclique/tclique.h"   /* def. of clique data structures */
#include "struct_vardata.h"

/** sets up the problem data */
extern
SCIP_RETCODE SCIPcreateProbGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< problem name */     
   );


/** sets up the problem data */
extern
SCIP_RETCODE GCGprobCreateFramework(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nblocks             /**< number of blocks */
   );


/** create the convexity constraints belonging to the pricing blocks */
extern
SCIP_RETCODE GCGprobCreateConvConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a linear constraint */
extern
SCIP_RETCODE GCGcreateConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< Is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );


/** creates a variable of the original program */
extern
SCIP_RETCODE GCGcreateOrigVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable           /**< is var's column removable from the LP (due to aging or cleanup)? */
   );

/** creates a variable of a pricing problem program */
extern
SCIP_RETCODE GCGcreatePricingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar,            /**< corresponding variable in the original program */
   int                   pricingprobnr       /**< number of the pricing problem to which the variable belongs */
   );


/** adds variable to the original problem */
extern
SCIP_RETCODE GCGaddOriginalVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to add */
   );

extern
SCIP_RETCODE GCGchgOrigVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

extern
SCIP_RETCODE GCGchgOrigVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

extern
SCIP_RETCODE GCGchgOrigVarType(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_VARTYPE          vartype             /**< new type of variable */
   );

extern
SCIP_RETCODE GCGprobSetOriginalVarBlockNr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to set the block number for */
   int                   blocknr             /**< number of the block, the variable belongs to */
   );

extern
SCIP* GCGprobGetOrigprob(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
SCIP* GCGprobGetPricingprob(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   );

extern
int GCGprobGetNPricingprobs(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
void GCGprobGetMasterConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          conss,
   int*                  nconss
   );

extern
int GCGprobGetNMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
void GCGprobGetOrigMasterConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          conss,
   int*                  nconss
   );

extern
void GCGprobGetLinearOrigMasterConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          conss,
   int*                  nconss
   );

extern
SCIP_CONS* GCGprobGetConvCons(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr
   );

/** gets values of multiple original variables w.r.t. primal master solution */
extern
SCIP_RETCODE GCGgetSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables to get solution value for */
   SCIP_VAR**            vars,               /**< array with variables to get value for */
   SCIP_Real*            vals                /**< array to store solution values of variables */
   );

#endif
