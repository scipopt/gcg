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

/**@file   cons_integralOrig.h
 * @brief  constraint handler for storing the graph at each node of the tree
 * @author Gerald Gamrath
 */

#ifndef CONSINTEGRALORIG_H
#define CONSINTEGRALORIG_H


/* type of integralOrig constraint: differ, same or root */
enum COLOR_ConsType
{
   COLOR_CONSTYPE_DIFFER = 0,  /* constraint representing the branching decision differ(i,j) */
   COLOR_CONSTYPE_SAME   = 1,  /* constraint representing the branching decision same(i,j) */
   COLOR_CONSTYPE_ROOT   = 2   /* constraint created for the root, is created automatically */
};
typedef enum COLOR_ConsType COLOR_CONSTYPE;

#include "scip/scip.h"
#include "tclique/tclique.h"

/** creates the handler for graph storing constraints and includes it in SCIP */
SCIP_RETCODE COLORincludeConshdlrIntegralOrig(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** creates and captures a integralOrig constraint, uses knowledge of the B&B-father*/
extern
SCIP_RETCODE COLORcreateConsIntegralOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONS*            fatherconstraint,   /**< constraint in B&B-father */
   int                   type,               /**< type of the constraint: ROOT for root-constraint, else SAME or DIFFER */
   int                   node1,              /**< the first node of the constraint or -1 if root-constraint */
   int                   node2,              /**< the second node of the constraint or -1 if root-constraint */
   SCIP_NODE*            stickingnode        /**< the B&B-tree node at which the constraint will be sticking */     
   );

#endif
