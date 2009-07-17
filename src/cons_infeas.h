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

/**@file   cons_infeas.h
 * @brief  constraint handler for storing the graph at each node of the tree
 * @author Gerald Gamrath
 */

#ifndef CONSINFEAS_H
#define CONSINFEAS_H

#include "scip/scip.h"

/** creates the handler for graph storing constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrInfeas(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
