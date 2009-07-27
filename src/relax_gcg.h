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

/**@file   relax_gcg.h
 * @brief  gcg relaxator
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RELAX_GCG_H__
#define __SCIP_RELAX_GCG_H__


#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "gcgplugins.h"
#include "struct_vardata.h"
#include "pricer_gcg.h"
#include "masterplugins.h"
#include "nodesel_master.h"




/** creates the gcg relaxator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeRelaxGcg(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a variable in a pricing problem corresponding to the given original variable */
extern
SCIP_RETCODE GCGrelaxCreatePricingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar             /**< corresponding variable in the original program */
   );

/** creates the data for a variable of the original program */
extern
SCIP_RETCODE GCGrelaxCreateOrigVardata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< pointer to variable object */
   );

/** creates the data for all variables of the original program */
extern
SCIP_RETCODE GCGrelaxCreateOrigVarsData(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* sets the number of the block, the given original variable belongs to */
extern
SCIP_RETCODE GCGrelaxSetOriginalVarBlockNr(
   SCIP_VAR*             var,                /**< variable to set the block number for */
   int                   blocknr             /**< number of the block, the variable belongs to */
   );

/* returns the master problem */
extern
SCIP* GCGrelaxGetMasterprob(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* returns the pricing problem of the given number */
extern
SCIP* GCGrelaxGetPricingprob(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   );

/* returns the number of pricing problems */
extern
int GCGrelaxGetNPricingprobs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* sets the number of pricing problems */
void GCGrelaxSetNPricingprobs(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   npricingprobs       /**< the number of pricing problems */
   );

/* returns the number of constraints in the master problem */
extern
int GCGrelaxGetNMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* returns the contraints in the master problem */
extern
SCIP_CONS** GCGrelaxGetMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* returns the contraints in the original problem that correspond to the constraints in the master problem */
extern
SCIP_CONS** GCGrelaxGetOrigMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );


/* returns the linear counterpart of the contraints in the original problem that correspond 
 * to the constraints in the master problem */
extern
SCIP_CONS** GCGrelaxGetLinearOrigMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* returns the convexity constraint for the given block */
extern
SCIP_CONS* GCGrelaxGetConvCons(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   blocknr             /**< the number of the block for which we 
                                              *   need the convexity constraint */   
   );

/* returns the current solution for the original problem */
extern
SCIP_SOL* GCGrelaxGetCurrentOrigSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transformes the current solution of the master problem into the original problem's space 
 *  and saves this solution as currentsol in the relaxator's data */
extern
SCIP_RETCODE GCGrelaxUpdateCurrentSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of fractional variables in the relaxator's current solution */
extern
int GCGrelaxGetNBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the fractional variables in the relaxator's current solution */
extern
SCIP_RETCODE GCGrelaxGetBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           branchcands,        /**< pointer to store the array of branching candidates, or NULL */
   SCIP_Real**           branchcandssol,     /**< pointer to store the array of candidate solution values, or NULL */
   SCIP_Real**           branchcandsfrac,    /**< pointer to store the array of candidate fractionalities, or NULL */
   int*                  nbranchcands,       /**< pointer to store the number of branching candidates, or NULL */
   int*                  npriobranchcands    /**< pointer to store the number of candidates with maximal priority, or NULL */
   );

/** returns solution value of the given variable in the relaxator's current solution */
extern
SCIP_Real GCGrelaxGetVarSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var
   );

/** returns solution value of the given variable in the relaxator's current solution */
extern
SCIP_RETCODE GCGrelaxLinkSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol
   );

/* transforms given values of the given original variables into values of the given master variables */
extern
void GCGrelaxTransformOrigvalsToMastervals(
   SCIP_VAR**            origvars,           /** array with (subset of the) original variables */
   SCIP_Real*            origvals,           /** array with values for the given original variables */
   int                   norigvars,          /** number of given original variables */
   SCIP_VAR**            mastervars,         /** array of (all present) master variables */
   SCIP_Real*            mastervals,         /** return value: values of the master variables */
   int                   nmastervars         /** number of master variables */
   );

#endif
