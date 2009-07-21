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

/**@file   cons_masterbranch.h
 * @brief  constraint handler for storing the graph at each node of the tree
 * @author Gerald Gamrath
 */

#ifndef CONSMASTERBRANCH_H
#define CONSMASTERBRANCH_H


/* sense of the masterbranch constraint: greater-equal or less-equal */
enum GCG_ConsSense
{
   GCG_CONSSENSE_GE   = 0,  /* greater-equal constraint */
   GCG_CONSSENSE_LE   = 1,  /* less-equal constraint */
   GCG_CONSSENSE_NONE = 2,  /* less-equal constraint */
};
typedef enum GCG_ConsSense GCG_CONSSENSE;

#include "scip/scip.h"

/** returns the store graph constraint of the current node */
extern
SCIP_CONS* GCGconsMasterbranchGetActiveCons(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates the handler for graph storing constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrMasterbranch(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a masterbranch constraint */
extern
SCIP_RETCODE GCGcreateConsMasterbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   SCIP_NODE*            node,
   SCIP_CONS*            parentcons
   );


/** returns the stack and the number of elements on it */
extern
void GCGconsMasterbranchGetStack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   );

/** returns the original variable for a given masterbranch constraint */
extern
SCIP_VAR* GCGconsMasterbranchGetOrigvar(
   SCIP_CONS*            cons
   );

/** returns the conssense for a given masterbranch constraint */
extern
GCG_CONSSENSE GCGconsMasterbranchGetConssense(
   SCIP_CONS*            cons
   );

/** returns the new bound for a given masterbranch constraint */
extern
SCIP_Real GCGconsMasterbranchGetVal(
   SCIP_CONS*            cons
   );

/** returns the node in the B&B tree at which the given masterbranch constraint is sticking */
extern
SCIP_NODE* GCGconsMasterbranchGetNode(
   SCIP_CONS*            cons
   );

/** returns the masterbranch constraint of the B&B father of the node at which the 
    given masterbranch constraint is sticking */
extern
SCIP_CONS* GCGconsMasterbranchGetParentcons(
   SCIP_CONS*            cons
   );

/** returns the masterbranch constraint of the first child of the node at which the 
    given masterbranch constraint is sticking */
extern
SCIP_CONS* GCGconsMasterbranchGetChild1cons(
   SCIP_CONS*            cons
   );

/** returns the masterbranch constraint of the second child of the node at which the 
    given masterbranch constraint is sticking */
extern
SCIP_CONS* GCGconsMasterbranchGetChild2cons(
   SCIP_CONS*            cons
   );

/** returns the origbranch constraint of the node in the original program corresponding to the node 
    at which the given masterbranch constraint is sticking */
extern
SCIP_CONS* GCGconsMasterbranchGetOrigcons(
   SCIP_CONS*            cons
   );

/** checks the consistency of the masterbranch constraints in the problem */
extern
void GCGconsMasterbranchCheckConsistency(
   SCIP*                 scip
   );

#endif
