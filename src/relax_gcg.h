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

/**@file    relax_gcg.h
 * @ingroup PUBLICMETHODS
 * @brief   GCG relaxator
 * @author  Gerald Gamrath
 * @author  Christian Puchert
 * @author  Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_RELAX_GCG_H__
#define GCG_RELAX_GCG_H__

#include "scip/scip.h"
#include "type_branchgcg.h"
#include "type_decomp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the GCG relaxator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeRelaxGcg(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes a branching rule into the relaxator data */
extern
SCIP_RETCODE GCGrelaxIncludeBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule for which callback methods are saved */
   GCG_DECL_BRANCHACTIVEMASTER((*branchactivemaster)),/**<  activation method for branchrule */
   GCG_DECL_BRANCHDEACTIVEMASTER ((*branchdeactivemaster)),/**<  deactivation method for branchrule */
   GCG_DECL_BRANCHPROPMASTER((*branchpropmaster)),/**<  propagation method for branchrule */
   GCG_DECL_BRANCHMASTERSOLVED((*branchmastersolved)),/**<  master solved method for branchrule */
   GCG_DECL_BRANCHDATADELETE((*branchdatadelete))/**<  branchdata deletion method for branchrule */
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

/** perform propagation method of the given branchrule for the given branchdata */
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

/** transformes a constraint of the original problem into the master variable space and
 *  adds it to the master problem */
extern
SCIP_RETCODE GCGrelaxTransOrigToMasterCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint that should be transformed */
   SCIP_CONS**           transcons           /**< pointer to the transformed constraint */
   );

/** returns the current solution for the original problem */
extern
SCIP_SOL* GCGrelaxGetCurrentOrigSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the current solution is primal feasible in the original problem */
extern
SCIP_Bool GCGrelaxIsOrigSolFeasible(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** start probing mode on both the original and master problems
 *
 *  @note This mode is intended for working on the original variables but using the master LP;
 *        it currently only supports bound changes on the original variables,
 *        but no additional rows
 */
extern
SCIP_RETCODE GCGrelaxStartProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            probingheur         /**< heuristic that started probing mode, or NULL */
   );

/** returns the  heuristic that started probing in the master problem, or NULL */
extern
SCIP_HEUR* GCGrelaxGetProbingheur(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** add a new probing node the original problem together with an original branching constraint
 *
 *  @note A corresponding probing node must be added to the master problem right before solving the probing LP
 */
extern
SCIP_RETCODE GCGrelaxNewProbingnodeOrig(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** add a new probing node the master problem together with a master branching constraint
 *  which ensures that bound changes are transferred to master and pricing problems
 *
 *  @note A corresponding probing node must have been added to the original problem beforehand;
 *        furthermore, this method must be called after bound changes to the original problem have been made
 */
extern
SCIP_RETCODE GCGrelaxNewProbingnodeMaster(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** add probing nodes to both the original and master problem;
 *  furthermore, add origbranch and masterbranch constraints to transfer branching decisions
 *  from the original to the master problem
 */
extern
SCIP_RETCODE GCGrelaxBacktrackProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   probingdepth        /**< probing depth of the node in the probing path that should be reactivated */
   );

/** solve the master probing LP without pricing */
extern
SCIP_RETCODE GCGrelaxPerformProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   maxlpiterations,    /**< maximum number of lp iterations allowed */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of performed LP iterations (or NULL) */
   SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   );

/** solve the master probing LP with pricing */
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
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   );

/** end probing mode in both the original and master problems */
extern
SCIP_RETCODE GCGrelaxEndProbing(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transforms the current solution of the master problem into the original problem's space
 *  and saves this solution as currentsol in the relaxator's data */
extern
SCIP_RETCODE GCGrelaxUpdateCurrentSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the decomposition mode */
extern
DEC_DECMODE GCGgetDecompositionMode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the decomposition mode of the master problem. The mode is given by the existence of either the GCG pricer or
 * the GCG Benders' decomposition plugins.
 */
extern
DEC_DECMODE GCGgetMasterDecompMode(
   SCIP*                 masterprob          /**< the master problem SCIP instance */
   );

#ifdef __cplusplus
}
#endif

#endif
