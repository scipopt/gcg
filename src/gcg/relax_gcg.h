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

/**@file    relax_gcg.h
 * @ingroup RELAXATORS-GCG
 * @brief   GCG relaxator
 * @author  Gerald Gamrath
 * @author  Christian Puchert
 * @author  Martin Bergner
 * @author  Oliver Gaul
 * @author  Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_RELAX_GCG_H__
#define GCG_RELAX_GCG_H__

#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @ingroup RELAXATORS-GCG
 * @{
 */

/** creates the GCG relaxator and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeRelaxGcg(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** includes a branching rule into the relaxator data */
GCG_EXPORT
SCIP_RETCODE GCGrelaxIncludeBranchrule(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE**     branchrule,         /**< branching rule for which callback methods are saved */
   GCG_BRANCHRULE**      gcgbranchrule,      /**< pointer to store created GCG branch rule (can be NULL) */
   const char*           name,               /**< name of branching rule */
   const char*           desc,               /**< description of branching rule */
   int                   priority,           /**< priority of the branching rule */
   int                   maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying branching rule
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< branching rule data */
   GCG_DECL_BRANCHACTIVEMASTER((*branchactivemaster)),/**<  activation method for branchrule */
   GCG_DECL_BRANCHDEACTIVEMASTER((*branchdeactivemaster)),/**<  deactivation method for branchrule */
   GCG_DECL_BRANCHPROPMASTER((*branchpropmaster)),/**<  propagation method for branchrule */
   GCG_DECL_BRANCHMASTERSOLVED((*branchmastersolved)),/**<  master solved method for branchrule */
   GCG_DECL_BRANCHDATADELETE((*branchdatadelete)),/**<  branchdata deletion method for branchrule */
   GCG_DECL_BRANCHNEWCOL ((*branchnewcol)),  /**< new column handler method of branching rule */
   GCG_DECL_BRANCHGETEXTENDEDMASTERCONS((*branchgetextendedmastercons)), /**< extended master cons getter of branching rule */
   GCG_DECL_BRANCHGETEXTENDEDMASTERCONSCOEFF((*branchgetextendedmasterconscoeff)) /**< column coefficient calculation method for extended master conss */
   );

/** perform activation method of the given branchrule for the given branchdata */
GCG_EXPORT
SCIP_RETCODE GCGrelaxBranchActiveMaster(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE*       branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata          /**< data representing the branching decision */
   );

/** perform deactivation method of the given branchrule for the given branchdata */
GCG_EXPORT
SCIP_RETCODE GCGrelaxBranchDeactiveMaster(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE*       branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata          /**< data representing the branching decision */
   );

/** perform propagation method of the given branchrule for the given branchdata */
GCG_EXPORT
SCIP_RETCODE GCGrelaxBranchPropMaster(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE*       branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation call */
   );

/** perform method of the given branchrule that is called after the master LP is solved */
GCG_EXPORT
SCIP_RETCODE GCGrelaxBranchMasterSolved(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE*       branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   SCIP_Real             newlowerbound       /**< the new local lowerbound */
   );

/** frees branching data created by the given branchrule */
GCG_EXPORT
SCIP_RETCODE GCGrelaxBranchDataDelete(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE*       branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA**      branchdata,         /**< data representing the branching decision */
   SCIP_Bool             origbranch,         /**< true iff an origbranch triggered this call */
   SCIP_Bool             force               /**< branch data must be deleted if true */
   );

/** notifies the branching rule that a new mastervariable was created while this node was active */
GCG_EXPORT
SCIP_RETCODE GCGrelaxBranchNewCol(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE*       branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   SCIP_VAR*             mastervar           /**< new mastervariable that was created */
   );

/** gets the extendedmasterconsdata created by this branching rule, if any */
GCG_EXPORT
SCIP_RETCODE GCGrelaxBranchGetExtendedMasterCons(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE*       branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   GCG_EXTENDEDMASTERCONSDATA**   extendedmasterconsdata       /**< the extendedmasterconsdata to grab */
   );

/** get extended master conss of all active nods */
GCG_EXPORT
SCIP_RETCODE GCGrelaxBranchGetAllActiveExtendedMasterConss(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE***     branchrules,        /**< branching rules that created extended master conss (can be NULL) */
   GCG_BRANCHDATA***     branchdata,         /**< data represeting the branching decisions of the active nodes (can be NULL) */
   GCG_EXTENDEDMASTERCONSDATA***  extendedmasterconsdata,      /**< array of extended master conss generated by branching in all currently active nodes */
   int*                  nextendedmasterconss         /**< number of currently active branching rules that created extended master conss */
   );

/** transformes a constraint of the original problem into the master variable space and
 *  adds it to the master problem */
GCG_EXPORT
SCIP_RETCODE GCGrelaxTransOrigToMasterCons(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS*            cons,               /**< the constraint that should be transformed */
   SCIP_CONS**           transcons           /**< pointer to the transformed constraint */
   );

/** returns the current solution for the original problem */
GCG_EXPORT
SCIP_SOL* GCGrelaxGetCurrentOrigSol(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns whether the current solution is primal feasible in the original problem */
GCG_EXPORT
SCIP_Bool GCGrelaxIsOrigSolFeasible(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** start probing mode on both the original and master problems
 *
 *  @note This mode is intended for working on the original variables but using the master LP;
 *        it currently only supports bound changes on the original variables,
 *        but no additional rows
 */
GCG_EXPORT
SCIP_RETCODE GCGrelaxStartProbing(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_HEUR*            probingheur         /**< heuristic that started probing mode, or NULL */
   );

/** returns the  heuristic that started probing in the master problem, or NULL */
GCG_EXPORT
SCIP_HEUR* GCGrelaxGetProbingheur(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** add a new probing node the original problem together with an original branching constraint
 *
 *  @note A corresponding probing node must be added to the master problem right before solving the probing LP
 */
GCG_EXPORT
SCIP_RETCODE GCGrelaxNewProbingnodeOrig(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** add a new probing node the master problem together with a master branching constraint
 *  which ensures that bound changes are transferred to master and pricing problems
 *
 *  @note A corresponding probing node must have been added to the original problem beforehand;
 *        furthermore, this method must be called after bound changes to the original problem have been made
 */
GCG_EXPORT
SCIP_RETCODE GCGrelaxNewProbingnodeMaster(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** add a new probing node the master problem together with a master branching constraint
 *  which ensures that bound changes are transferred to master and pricing problems as well as additional
 *  constraints
 *
 *  @note A corresponding probing node must have been added to the original problem beforehand;
 *        furthermore, this method must be called after bound changes to the original problem have been made
 */
GCG_EXPORT
SCIP_RETCODE GCGrelaxNewProbingnodeMasterCons(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE*       branchrule,         /**< pointer to the branching rule */
   GCG_BRANCHDATA*       branchdata,         /**< branching data */
   SCIP_CONS**           origbranchconss,    /**< original constraints enforcing the branching decision */
   int                   norigbranchconss,   /**< number of original constraints */
   int                   maxorigbranchconss  /**< capacity of origbranchconss */
   );

/** add probing nodes to both the original and master problem;
 *  furthermore, add origbranch and masterbranch constraints to transfer branching decisions
 *  from the original to the master problem
 */
GCG_EXPORT
SCIP_RETCODE GCGrelaxBacktrackProbing(
   GCG*                  gcg,                /**< GCG data structure */
   int                   probingdepth        /**< probing depth of the node in the probing path that should be reactivated */
   );

/** solve the master probing LP without pricing */
GCG_EXPORT
SCIP_RETCODE GCGrelaxPerformProbing(
   GCG*                  gcg,                /**< GCG data structure */
   int                   maxlpiterations,    /**< maximum number of lp iterations allowed */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of performed LP iterations (or NULL) */
   SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   );

/** solve the master probing LP with pricing */
GCG_EXPORT
SCIP_RETCODE GCGrelaxPerformProbingWithPricing(
   GCG*                  gcg,                /**< GCG data structure */
   int                   maxpricerounds,     /**< maximum number of pricing rounds allowed */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of performed LP iterations (or NULL) */
   int*                  npricerounds,       /**< pointer to store the number of performed pricing rounds (or NULL) */
   SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   );

/** end probing mode in both the original and master problems */
GCG_EXPORT
SCIP_RETCODE GCGrelaxEndProbing(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** transforms the current solution of the master problem into the original problem's space
 *  and saves this solution as currentsol in the relaxator's data */
GCG_EXPORT
SCIP_RETCODE GCGrelaxUpdateCurrentSol(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** initialize master problem for solving
 */
SCIP_RETCODE GCGinitializeMasterProblemSolve(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** stash limit settings if not already stashed
 * 
 * @note This function is used to prevent that SCIP interrupts (due to limits) the solving process when the B&B trees are not synchronized.
 * Otherwise it may happen that SCIP creates a dummy node in one tree that cannot be mirrored to the other one.
 */
SCIP_RETCODE GCGstashLimitSettings(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** restore limit settings if currently stashed
 */
SCIP_RETCODE GCGrestoreLimitSettings(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef _OPENMP
/** returns OpenMP locks
 *  @returns pointer to GCG_LOCKS struct
 */
GCG_LOCKS* GCGgetLocks(
   GCG*                  gcg                 /**< GCG data structure */
   );
#endif

/** sets the pricing problem parameters */
SCIP_RETCODE GCGsetPricingProblemParameters(
   GCG_DECTYPE           dectype,            /**< the dectype of the decomp */
   SCIP*                 pricingprob,        /**< SCIP data structure of the pricing problem */
   int                   clocktype,          /**< clocktype to use in the pricing problem */
   SCIP_Real             infinity,           /**< values larger than this are considered infinity in the pricing problem */
   SCIP_Real             epsilon,            /**< absolute values smaller than this are considered zero in the pricing problem */
   SCIP_Real             sumepsilon,         /**< absolute values of sums smaller than this are considered zero in the pricing problem */
   SCIP_Real             feastol,            /**< feasibility tolerance for constraints in the pricing problem */
   SCIP_Real             lpfeastolfactor,    /**< primal feasibility tolerance factor of LP solver in the pricing problem */
   SCIP_Real             dualfeastol,        /**< feasibility tolerance for reduced costs in LP solution in the pricing problem */
   SCIP_Bool             enableppcuts        /**< should ppcuts be stored for sepa_basis */
   );

/** returns the GCG data structure */
GCG_EXPORT
GCG* GCGrelaxGetGcg(
   SCIP*                 origprob            /**< SCIP data structure */
   );

/** includes a separator into the relaxator data */
GCG_EXPORT
SCIP_RETCODE GCGrelaxIncludeSepa(
   GCG*                  gcg,                /**< SCIP data structure */
   SCIP_SEPA**           sepa,               /**< pointer to store created scip separator */
   GCG_SEPA**            gcgsepa,            /**< pointer to store created GCG separator (can be NULL) */
   const char*           name,               /**< name of separator */
   const char*           desc,               /**< description of separator */
   int                   priority,           /**< priority of separator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling separator */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   SCIP_Bool             usessubscip,        /**< does the separator use a secondary SCIP instance? */
   SCIP_Bool             delay,              /**< should separator be delayed, if other separators found cuts? */
   SCIP_DECL_SEPAEXECLP  ((*sepaexeclp)),    /**< LP solution separation method of separator */
   SCIP_DECL_SEPAEXECSOL ((*sepaexecsol)),   /**< arbitrary primal solution separation method of separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   GCG_DECL_SEPAGETCOLCOEFFICIENT((*sepagetcolcoef)),/**< method for computing the column coefficient for a cut */
   GCG_DECL_SEPAMASTERCUTDELETE((*sepamastercutdelete))/**< callback to delete the mastersepacutdata */
   );

#ifdef __cplusplus
}

#endif

#endif
