/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   relax_gcg.h
 * @brief  gcg relaxator
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RELAX_GCG_H__
#define __SCIP_RELAX_GCG_H__

#include "scip/scip.h"
#include "type_branchgcg.h"
#include "type_decomp.h"

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

/** transformes a constraint of the original problem into the master variable space and
 *  adds it to the master problem */
extern
SCIP_RETCODE GCGrelaxTransOrigToMasterCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint that should be transformed */
   SCIP_CONS**           transcons           /**< pointer to the transformed constraint */
   );


/** marks the constraint to be a master constraint */
extern
SCIP_RETCODE GCGrelaxMarkConsMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint that is forced to be in the master */
   );

/** returns the master problem */
extern
SCIP* GCGrelaxGetMasterprob(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the pricing problem of the given number */
extern
SCIP* GCGrelaxGetPricingprob(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   );

/** returns the number of pricing problems */
extern
int GCGrelaxGetNPricingprobs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns TRUE iff the pricingproblem of the given number is relevant, that means is not identical to
 *  another and represented by it */
extern
SCIP_Bool GCGrelaxIsPricingprobRelevant(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   );

/**
 *  for a given block, return the block by which it is represented
 */
extern
int GCGrelaxGetBlockRepresentative(
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

/** returns the number of constraints in the master problem */
extern
int GCGrelaxGetNMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the contraints in the master problem */
extern
SCIP_CONS** GCGrelaxGetMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the contraints in the original problem that correspond to the constraints in the master problem */
extern
SCIP_CONS** GCGrelaxGetOrigMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** returns the linear counterpart of the contraints in the original problem that correspond
 * to the constraints in the master problem */
extern
SCIP_CONS** GCGrelaxGetLinearOrigMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the convexity constraint for the given block */
extern
SCIP_CONS* GCGrelaxGetConvCons(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   blocknr             /**< the number of the block for which we
                                              *   need the convexity constraint */
   );

/** returns the current solution for the original problem */
extern
SCIP_SOL* GCGrelaxGetCurrentOrigSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** start probing mode on master problem */
extern
SCIP_RETCODE GCGrelaxStartProbing(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** for a probing node in the original problem, create a corresponding probing node in the master problem,
 *  propagate domains and solve the LP without pricing. */
extern
SCIP_RETCODE GCGrelaxPerformProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint          maxlpiterations,    /**< maximum number of lp iterations allowed */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of performed LP iterations (or NULL) */
   SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the probing direction is infeasible */
   SCIP_Bool*            feasible            /**< pointer to store whether the probing solution is feasible */
   );

/** for a probing node in the original problem, create a corresponding probing node in the master problem,
 *  propagate domains and solve the LP with pricing. */
extern
SCIP_RETCODE GCGrelaxPerformProbingWithPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   maxpricerounds,     /**< maximum number of pricing rounds allowed */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of performed LP iterations (or NULL) */
   int*                  npricerounds,       /**< pointer to store the number of performed pricing rounds (or NULL) */
   SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the probing direction is infeasible */
   SCIP_Bool*            feasible            /**< pointer to store whether the probing solution is feasible */
   );

/** end probing mode in master problem */
extern
SCIP_RETCODE GCGrelaxEndProbing(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** transforms the current solution of the master problem into the original problem's space
 *  and saves this solution as currentsol in the relaxator's data */
extern
SCIP_RETCODE GCGrelaxUpdateCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            feasible            /**< pointer to store whether the master problem's solution is
                                              *   primal feasible*/
   );

/** transforms given values of the given original variables into values of the given master variables */
extern
void GCGrelaxTransformOrigvalsToMastervals(
   SCIP*                 scip,               /** SCIP data structure */
   SCIP_VAR**            origvars,           /** array with (subset of the) original variables */
   SCIP_Real*            origvals,           /** array with values for the given original variables */
   int                   norigvars,          /** number of given original variables */
   SCIP_VAR**            mastervars,         /** array of (all present) master variables */
   SCIP_Real*            mastervals,         /** array to store the values of the master variables */
   int                   nmastervars         /** number of master variables */
   );

/** transforms given solution of the master problem into solution of the original problem */
extern
SCIP_RETCODE GCGrelaxTransformMastersolToOrigsol(
   SCIP*                 scip,               /** SCIP data structure */
   SCIP_SOL*             mastersol,          /** solution of the master problem */
   SCIP_SOL**            origsol             /** pointer to store the new created original problem's solution */
   );

/** prints the given variable: name, type (original, master or pricing) block number,
 * and the list of all variables related to the given variable */
extern
void GCGrelaxPrintVar(
   SCIP_VAR*             var                 /**< variable that shpuld be printed */
   );

/** returns the stored primal solution of the original problem  */
extern
SCIP_SOL* GCGrelaxGetOrigPrimalSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the stored primal solution of the original problem  */
extern
void GCGrelaxSetOrigPrimalSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< solution */
   );

/** sets the structure information */
void GCGsetStructDecdecomp(
   SCIP*       scip,       /**< SCIP data structure */
   DECDECOMP*  decdecomp   /**< decomposition data structure */
   );

/** gets the structure information */
DECDECOMP* GCGgetStructDecdecomp(
   SCIP*       scip        /**< SCIP data structure */
   );

#endif
