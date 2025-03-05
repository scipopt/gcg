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

/**@file    gcgvarhistory.h
 * @ingroup INTERNALAPI-GCG
 * @brief   methods for managing variable history
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef GCG_GCGVARHISTORY_H__
#define GCG_GCGVARHISTORY_H__

#include "scip/scip.h"
#include "gcg/def.h"
#include "gcg/type_gcgvarhistory.h"

#ifdef __cplusplus
extern "C" {
#endif

/** get the variable behinde the pointer */
SCIP_RETCODE GCGvarhistoryGetVar(
   GCG_VARHISTORY*        pointer,           /**< pointer to the history */
   SCIP_VAR**             var                /**< pointer to store the variable */
   );

/** check if there is a next history event */
SCIP_Bool GCGvarhistoryHasNext(
   GCG_VARHISTORY*        pointer            /**< pointer to the history */
   );

/** get the next history event */
SCIP_RETCODE GCGvarhistoryNext(
   SCIP*                  scip,              /**< SCIP data structure */
   GCG_VARHISTORY**       pointer            /**< pointer to the history */
   );

/** jump to the latest history event */
SCIP_RETCODE GCGvarhistoryJumpToLatest(
   SCIP*                  scip,              /**< SCIP data structure */
   GCG_VARHISTORY**       pointer            /**< pointer to the history */
   );

/** jump to the latest history event and retrieve all new variables */
SCIP_RETCODE GCGvarhistoryJumpAndRetrieveVars(
   SCIP*                  scip,              /**< SCIP data structure */
   GCG_VARHISTORY**       pointer,           /**< pointer to the history */
   SCIP_VAR***            vars,              /**< pointer to store the variables */
   int*                   nvars              /**< pointer to store the number of variables */
   );

/** create a new history pointer to an empty existing buffer and captures it */
SCIP_RETCODE GCGvarhistoryCreate(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORY**       pointer             /**< pointer to the history */
   );

/** copy a pointer by creating a new one that points to the same buffer at the same position and capture it */
SCIP_RETCODE GCGvarhistoryCopyReference(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORY**       pointer,            /**< pointer to the history */
   GCG_VARHISTORY*        source              /**< source pointer */
   );

/** release the reference to the buffer and free the history pointer */
SCIP_RETCODE GCGvarhistoryFreeReference(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORY**       pointer             /**< pointer to the history */
   );

/** add variable to history */
SCIP_RETCODE GCGvarhistoryAddVar(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORY*        pointer,            /**< pointer to the history */
   SCIP_VAR*              var                 /**< variable */
   );

/** capture a reference to the history */
SCIP_RETCODE GCGvarhistoryCaptureBuffer(
   GCG_VARHISTORYBUFFER*  buffer              /**< buffer */
   );

/** release a reference to the history */
SCIP_RETCODE GCGvarhistoryReleaseBuffer(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORYBUFFER** buffer              /**< buffer */
   );


#ifdef __cplusplus
}
#endif

#endif
