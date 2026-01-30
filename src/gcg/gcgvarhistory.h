/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file    gcgvarhistory.h
 * @ingroup INTERNALAPI-GCG
 * @brief   methods for managing variable history
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef GCG_GCGVARHISTORY_H__
#define GCG_GCGVARHISTORY_H__

#include "scip/scip.h"

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
   GCG_VARHISTORY*        pointer            /**< pointer to the history */
   );

/** jump to the latest history event and retrieve all new variables */
SCIP_RETCODE GCGvarhistoryJumpAndRetrieveVars(
   SCIP*                  scip,              /**< SCIP data structure */
   GCG_VARHISTORY*        pointer,           /**< pointer to the history */
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
