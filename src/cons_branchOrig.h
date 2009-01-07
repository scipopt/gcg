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

/**@file   cons_branchOrig.h
 * @brief  constraint handler for storing the graph at each node of the tree
 * @author Gerald Gamrath
 */

#ifndef CONSBRANCHORIG_H
#define CONSBRANCHORIG_H


/* type of branchOrig constraint: differ, same or root */
enum GCG_ConsSense
{
   GCG_CONSSENSE_GE = 0,  /* greater-equal constraint */
   GCG_CONSSENSE_LE = 1,  /* less-equal constraint */
};
typedef enum GCG_ConsSense GCG_CONSSENSE;

#include "scip/scip.h"
#include "tclique/tclique.h"

/** returns the store graph constraint of the current node, needs the pointer to the constraint handler */
extern
SCIP_CONS* GCGconsGetActiveBranchOrigConsFromHandler(
   SCIP_CONSHDLR*        conshdlr            /**< constaint handler for store-graph constraints */
   );


/** returns the store graph constraint of the current node, needs only the pointer to scip */
extern
SCIP_CONS* GCGconsGetActiveBranchOrigCons(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** returns array of representatives of all nodes */
extern
int* GCGconsGetRepresentatives(
   SCIP*                 scip                 /**< SCIP data structure */
   );


/** creates the handler for graph storing constraints and includes it in SCIP */
SCIP_RETCODE GCGincludeConshdlrBranchOrig(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** creates and captures a branchOrig constraint, uses knowledge of the B&B-father*/
extern
SCIP_RETCODE GCGcreateConsBranchOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONS*            fatherconstraint,   /**< constraint in B&B-father */
   int                   type,               /**< type of the constraint: ROOT for root-constraint, else SAME or DIFFER */
   int                   node1,              /**< the first node of the constraint or -1 if root-constraint */
   int                   node2,              /**< the second node of the constraint or -1 if root-constraint */
   SCIP_NODE*            stickingnode        /**< the B&B-tree node at which the constraint will be sticking */     
   );


/** returns the stack and the number of elements on it */
extern
void GCGconsGetStack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   );


#endif
