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

/**@file   cons_origbranch.h
 * @brief  constraint handler for storing the graph at each node of the tree
 * @author Gerald Gamrath
 */

#ifndef CONSORIGBRANCH_H
#define CONSORIGBRANCH_H

#include "scip/scip.h"
#include "cons_masterbranch.h"
#include "struct_branchgcg.h"


/** returns the store graph constraint of the current node, needs only the pointer to scip */
extern
SCIP_CONS* GCGconsOrigbranchGetActiveCons(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** creates the handler for graph storing constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrOrigbranch(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a origbranch constraint*/
extern
SCIP_RETCODE GCGcreateConsOrigbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */          
   SCIP_CONS*            branchcons,         /**< linear constraint in the original problem */
   SCIP_VAR*             origvar,
   GCG_CONSSENSE         conssense,
   SCIP_Real             val,
   SCIP_NODE*            node,
   SCIP_CONS*            parentcons,
   SCIP_BRANCHRULE*      branchrule,
   GCG_BRANCHDATA*       branchdata
   );

/** returns the stack and the number of elements on it */
extern
void GCGconsOrigbranchGetStack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   );

/** returns the branch orig constraint of the current node, only needs the pointer to scip */
extern
SCIP_VAR* GCGconsOrigbranchGetOrigvar(
   SCIP_CONS*            cons
   );

/** returns the branch orig constraint of the current node, only needs the pointer to scip */
extern
GCG_CONSSENSE GCGconsOrigbranchGetConssense(
   SCIP_CONS*            cons
   );

/** returns the branch orig constraint of the current node, only needs the pointer to scip */
extern
SCIP_Real GCGconsOrigbranchGetVal(
   SCIP_CONS*            cons
   );

/** returns the branching data for a given origbranch constraint */
extern
GCG_BRANCHDATA* GCGconsOrigbranchGetBranchdata(
   SCIP_CONS*            cons
   );

/** returns the branchrule for a given origbranch constraint */
extern
SCIP_BRANCHRULE* GCGconsOrigbranchGetBranchrule(
   SCIP_CONS*            cons
   );

/** returns the node in the B&B tree at which the given origbranch constraint is sticking */
extern
SCIP_NODE* GCGconsOrigbranchGetNode(
   SCIP_CONS*            cons
   );

/** returns the origbranch constraint of the B&B father of the node at which the 
    given origbranch constraint is sticking */
extern
SCIP_CONS* GCGconsOrigbranchGetParentcons(
   SCIP_CONS*            cons
   );

/** returns the origbranch constraint of the first child of the node at which the 
    given origbranch constraint is sticking */
extern
SCIP_CONS* GCGconsOrigbranchGetChild1cons(
   SCIP_CONS*            cons
   );

/** returns the origbranch constraint of the second child of the node at which the 
    given origbranch constraint is sticking */
extern
SCIP_CONS* GCGconsOrigbranchGetChild2cons(
   SCIP_CONS*            cons
   );

/** sets the masterbranch constraint of the node in the master program corresponding to the node 
    at which the given origbranchbranch constraint is sticking */
extern
void GCGconsOrigbranchSetMastercons(
   SCIP_CONS*            cons,
   SCIP_CONS*            mastercons
   );

/** returns the masterbranch constraint of the node in the master program corresponding to the node 
    at which the given origbranchbranch constraint is sticking */
extern
SCIP_CONS* GCGconsOrigbranchGetMastercons(
   SCIP_CONS*            cons
   );

/** checks the consistency of the origbranch constraints in the problem */
extern
void GCGconsOrigbranchCheckConsistency(
   SCIP*                 scip
   );

#endif
