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

/** creates the data for a variable of the original program */
extern
SCIP_RETCODE GCGcreateOrigVardata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< pointer to variable object */
   );

/** creates the data for all variables of the original program */
extern
SCIP_RETCODE GCGcreateOrigVarsData(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a variable of a pricing problem program */
extern
SCIP_RETCODE GCGcreatePricingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar             /**< corresponding variable in the original program */
   );

extern
SCIP* GCGprobGetMasterprob(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
SCIP_RETCODE GCGprobSetOriginalVarBlockNr(
   SCIP_VAR*             var,                /**< variable to set the block number for */
   int                   blocknr             /**< number of the block, the variable belongs to */
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
