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


/** returns the store graph constraint of the current node, needs only the pointer to scip */
extern
SCIP_CONS* GCGconsGetActiveOrigbranchCons(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** creates the handler for graph storing constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrOrigbranch(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a origbranch constraint*/
SCIP_RETCODE GCGcreateConsOrigbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */          
   SCIP_CONS*            branchcons          /**< linear constraint in the original problem */
   );

/** returns the stack and the number of elements on it */
extern
void GCGconsGetStack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   );


#endif
