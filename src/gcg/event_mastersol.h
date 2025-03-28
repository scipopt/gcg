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

/**@file   event_mastersol.h
 * @ingroup EVENTS
 * @brief  eventhdlr to transfer solutions found in the original problem to the master problem
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_EVENT_MASTERSOL_H__
#define GCG_EVENT_MASTERSOL_H__

#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates event handler for mastersol event */
GCG_EXPORT
SCIP_RETCODE GCGincludeEventHdlrMastersol(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** return whether event has been triggered */
GCG_EXPORT
SCIP_Bool GCGeventhdlrMastersolIsTriggered(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
