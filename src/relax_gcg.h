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
#include "type_branchgcg.h"



/** creates the gcg relaxator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeRelaxGcg(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes a branching rule into the relaxator data */
extern
SCIP_RETCODE GCGrelaxIncludeBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule for which callback methods are saved */
   GCG_DECL_BRANCHACTIVEMASTER    ((*branchactivemaster)),    /**<  activation method for branchrule */
   GCG_DECL_BRANCHDEACTIVEMASTER  ((*branchdeactivemaster)),  /**<  deactivation method for branchrule */
   GCG_DECL_BRANCHPROPMASTER      ((*branchpropmaster)),      /**<  propagation method for branchrule */
   GCG_DECL_BRANCHMASTERSOLVED    ((*branchmastersolved)),    /**<  master solved method for branchrule */
   GCG_DECL_BRANCHDATADELETE      ((*branchdatadelete))       /**<  branchdata deletion method for branchrule */
   );

/** perform activation method of the given branchrule for the given branchdata */
extern
SCIP_RETCODE GCGrelaxBranchActiveMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata          /**< data representing the branching decision */
   );

/** perform deactivation method of the given branchrule for the given branchdata */
extern
SCIP_RETCODE GCGrelaxBranchDeactiveMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata          /**< data representing the branching decision */
   );

/** perform popagation method of the given branchrule for the given branchdata */
extern
SCIP_RETCODE GCGrelaxBranchPropMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation call */
   );

/** perform method of the given branchrule that is called after the master LP is solved */
extern
SCIP_RETCODE GCGrelaxBranchMasterSolved(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   SCIP_Real             newlowerbound       /**< the new local lowerbound */
   );

/** frees branching data created by the given branchrule */
extern
SCIP_RETCODE GCGrelaxBranchDataDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA**      branchdata          /**< data representing the branching decision */
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

/* marks the constraint to be a master constraint */
extern
SCIP_RETCODE GCGrelaxMarkConsMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint that is forced to be in the master */
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

/** returns TRUE iff the pricingproblem of the given number is relevant, that means is not identical to
 *  another and represented by it */
extern
SCIP_Bool GCGrelaxIsPricingprobRelevant(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   );

/** returns the number of blocks in the original formulation, that are represented by 
 *  the pricingprob with the given number */
extern
int GCGrelaxGetNIdenticalBlocks(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
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

/* transforms given values of the given original variables into values of the given master variables */
extern
void GCGrelaxTransformOrigvalsToMastervals(
   SCIP*                 scip,               /** SCIP data structure */
   SCIP_VAR**            origvars,           /** array with (subset of the) original variables */
   SCIP_Real*            origvals,           /** array with values for the given original variables */
   int                   norigvars,          /** number of given original variables */
   SCIP_VAR**            mastervars,         /** array of (all present) master variables */
   SCIP_Real*            mastervals,         /** return value: values of the master variables */
   int                   nmastervars         /** number of master variables */
   );

/* transforms given solution of the master problem into solution of the original problem */
extern
SCIP_RETCODE GCGrelaxTransformMastersolToOrigsol(
   SCIP*                 scip,               /** SCIP data structure */
   SCIP_SOL*             mastersol,          /** solution of the master problem */
   SCIP_SOL**            origsol             /** pointer to store the new created original problem's solution */
   );
#endif
