/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   cons_origbranch.h
 * @ingroup CONSHDLRS-GCG
 * @brief  constraint handler for storing the branching decisions at each node of the tree
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CONS_ORIGBRANCH_H__
#define GCG_CONS_ORIGBRANCH_H__

#include "scip/scip.h"
#include "gcg/type_branchgcg.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for origbranch constraints and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeConshdlrOrigbranch(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** creates and captures a origbranch constraint */
GCG_EXPORT
SCIP_RETCODE GCGcreateConsOrigbranch(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_NODE*            node,               /**< the node to which this origbranch constraint belongs */
   SCIP_CONS*            parentcons,         /**< origbranch constraint associated with the father node */
   SCIP_BRANCHRULE*      branchrule,         /**< the branching rule that created the b&b node the constraint belongs to */
   GCG_BRANCHDATA*       branchdata          /**< branching data storing information about the branching restrictions at the
                                              *   corresponding node */
   );

/** returns the branch orig constraint of the current node, only needs the pointer to scip */
GCG_EXPORT
SCIP_CONS* GCGconsOrigbranchGetActiveCons(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the stack and the number of elements on it */
GCG_EXPORT
void GCGconsOrigbranchGetStack(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   );

/** sets the branching data for a given origbranch constraint */
GCG_EXPORT
void GCGconsOrigbranchSetBranchdata(
   SCIP_CONS*            cons,               /**< origbranch constraint for which the branching data is requested */
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   );

/** returns the branching data for a given origbranch constraint */
GCG_EXPORT
GCG_BRANCHDATA* GCGconsOrigbranchGetBranchdata(
   SCIP_CONS*            cons                /**< origbranch constraint for which the branching data is requested */
   );

/** returns the branching rule for a given origbranch constraint */
GCG_EXPORT
SCIP_BRANCHRULE* GCGconsOrigbranchGetBranchrule(
   SCIP_CONS*            cons                /**< origbranch constraint for which the branchrule is requested */
   );

/** returns the node in the B&B tree at which the given origbranch constraint is sticking */
GCG_EXPORT
SCIP_NODE* GCGconsOrigbranchGetNode(
   SCIP_CONS*            cons                /**< origbranch constraint for which the corresponding node is requested */
   );

/** returns the origbranch constraint of the B&B father of the node at which the
  * given origbranch constraint is sticking
  */
GCG_EXPORT
SCIP_CONS* GCGconsOrigbranchGetParentcons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the origbranch constraint of
                                              *   the father node is requested */
   );

/** returns the number of origbranch constraints of the children of the node at which the
  * given origbranch constraint is sticking
  */
GCG_EXPORT
int GCGconsOrigbranchGetNChildconss(
   SCIP_CONS*            cons                /**< constraint pointer */
   );

/** returns an origbranch constraint of a child of the node at which the
  * given origbranch constraint is sticking
  */
GCG_EXPORT
SCIP_CONS* GCGconsOrigbranchGetChildcons(
   SCIP_CONS*            cons,               /**< constraint */
   int                   childnr             /**< number of child */
   );


/** sets the masterbranch constraint of the node in the master program corresponding to the node
  * at which the given origbranchbranch constraint is sticking
  */
GCG_EXPORT
void GCGconsOrigbranchSetMastercons(
   SCIP_CONS*            cons,               /**< origbranch constraint for which the masterbranch constraint should be set */
   SCIP_CONS*            mastercons          /**< masterbranch constraint corresponding to the given origbranch constraint */
   );

/** returns the masterbranch constraint of the node in the master program corresponding to the node
  * at which the given origbranchbranch constraint is sticking
  */
GCG_EXPORT
SCIP_CONS* GCGconsOrigbranchGetMastercons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the corresponding masterbranch
                                              *   constraint is requested */
   );

/** adds initial constraint to root node */
GCG_EXPORT
SCIP_RETCODE GCGconsOrigbranchAddRootCons(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** checks the consistency of the origbranch constraints in the problem */
GCG_EXPORT
void GCGconsOrigbranchCheckConsistency(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
