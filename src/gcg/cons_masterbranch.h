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

/**@file   cons_masterbranch.h
 * @ingroup CONSHDLRS-GCG
 * @brief  constraint handler for storing the branching decisions at each node of the tree
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @author Christian Puchert
 * @author Marcel Schmickerath
 */

#ifndef GCG_CONS_MASTERBRANCH_H__
#define GCG_CONS_MASTERBRANCH_H__

#include "scip/scip.h"
#include "gcg/type_branchgcg.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif


/** creates the handler for masterbranch constraints and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeConshdlrMasterbranch(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** creates and captures a masterbranch constraint */
GCG_EXPORT
SCIP_RETCODE GCGcreateConsMasterbranch(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of the constraint */
   SCIP_NODE*            node,               /**< node at which the constraint should be created */
   SCIP_CONS*            parentcons,         /**< parent constraint */
   SCIP_BRANCHRULE*      branchrule,         /**< pointer to the branching rule */
   GCG_BRANCHDATA*       branchdata,         /**< branching data */
   SCIP_CONS**           origbranchconss,    /**< original constraints enforcing the branching decision */
   int                   norigbranchconss,   /**< number of original constraints */
   int                   maxorigbranchconss  /**< capacity origbranchconss */
   );

/** returns the name of the constraint */
GCG_EXPORT
char* GCGconsMasterbranchGetName(
   SCIP_CONS*            cons                /**< masterbranch constraint for which the data is requested */
   );

/** returns the node in the B&B tree at which the given masterbranch constraint is sticking */
GCG_EXPORT
SCIP_NODE* GCGconsMasterbranchGetNode(
   SCIP_CONS*            cons                /**< constraint pointer */
   );

/** returns the masterbranch constraint of the B&B father of the node at which the
  * given masterbranch constraint is sticking
  */
GCG_EXPORT
SCIP_CONS* GCGconsMasterbranchGetParentcons(
   SCIP_CONS*            cons                /**< constraint pointer */
   );

/** returns the number of masterbranch constraints of the children of the node at which the
  * given masterbranch constraint is sticking
  */
GCG_EXPORT
int GCGconsMasterbranchGetNChildconss(
   SCIP_CONS*            cons                /**< constraint pointer */
   );

/** returns a masterbranch constraint of a child of the node at which the
  * given masterbranch constraint is sticking
  */
GCG_EXPORT
SCIP_CONS* GCGconsMasterbranchGetChildcons(
   SCIP_CONS*            cons,                /**< constraint pointer */
   int                   childnr              /**< index of the child node */
   );

/** returns the origbranch constraint of the node in the original program corresponding to the node
  * which the given masterbranch constraint is sticking
  */
GCG_EXPORT
SCIP_CONS* GCGconsMasterbranchGetOrigcons(
   SCIP_CONS*            cons                /**< constraint pointer */
   );

/** sets the origbranch constraint of the node in the master program corresponding to the node
  * at which the given masterbranchbranch constraint is sticking
  */
GCG_EXPORT
void GCGconsMasterbranchSetOrigcons(
   SCIP_CONS*            cons,               /**< constraint pointer */
   SCIP_CONS*            origcons            /**< original branching constraint */
   );

/** returns the branching data for a given masterbranch constraint */
GCG_EXPORT
GCG_BRANCHDATA* GCGconsMasterbranchGetBranchdata(
   SCIP_CONS*            cons                /**< constraint pointer */
   );

/** sets the branching data for a given masterbranch constraint */
GCG_EXPORT
void GCGconsMasterbranchSetBranchdata(
   SCIP_CONS*            cons,               /**< masterbranch constraint for which the branching data is requested */
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   );

/** returns the branching rule of the constraint */
GCG_EXPORT
SCIP_BRANCHRULE* GCGconsMasterbranchGetBranchrule(
   SCIP_CONS*            cons                /**< masterbranch constraint for which the data is requested */
   );

/** adds a bound change on an original variable that was directly copied to the master problem */
GCG_EXPORT
SCIP_RETCODE GCGconsMasterbranchAddCopiedVarBndchg(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS*            cons,               /**< masterbranch constraint to which the bound change is added */
   SCIP_VAR*             var,                /**< variable on which the bound change was performed */
   GCG_BOUNDTYPE         boundtype,          /**< bound type of the bound change */
   SCIP_Real             newbound            /**< new bound of the variable after the bound change */
   );

/** returns the constraints in the original problem that enforce the branching decision */
GCG_EXPORT
SCIP_CONS** GCGconsMasterbranchGetOrigbranchConss(
   SCIP_CONS*            cons                /**< masterbranch constraint for which the data is requested */
   );

/** returns the number of constraints in the original problem that enforce the branching decision */
GCG_EXPORT
int GCGconsMasterbranchGetNOrigbranchConss(
   SCIP_CONS*            cons                /**< masterbranch constraint for which the data is requested */
   );

/** releases the constraints in the original problem that enforce the branching decision
 *  and frees the array holding the constraints
 */
GCG_EXPORT
SCIP_RETCODE GCGconsMasterbranchReleaseOrigbranchConss(
   GCG*                  gcg,                /**< GCG instance */
   SCIP_CONS*            cons                /**< masterbranch constraint for which the data is freed */
   );

/** returns the masterbranch constraint of the current node */
GCG_EXPORT
SCIP_CONS* GCGconsMasterbranchGetActiveCons(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the stack and the number of elements on it */
GCG_EXPORT
void GCGconsMasterbranchGetStack(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   );

/** returns the number of elements on the stack */
GCG_EXPORT
int GCGconsMasterbranchGetNStackelements(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** adds initial constraint to root node */
GCG_EXPORT
SCIP_RETCODE GCGconsMasterbranchAddRootCons(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** check whether the node was generated by generic branching */
GCG_EXPORT
SCIP_Bool GCGcurrentNodeIsGeneric(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** checks the consistency of the masterbranch constraints in the problem */
GCG_EXPORT
void GCGconsMasterbranchCheckConsistency(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
