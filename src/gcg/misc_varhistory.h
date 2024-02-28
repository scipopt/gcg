/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       */
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

/**@file    misc_varhistory.h
 * @brief   methods for managing variable history
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef GCG_HEUR_XPRINS_H__
#define GCG_HEUR_XPRINS_H__

#include "scip/scip.h"
#include "def.h"

#define GCG_VARHISTORYBUFFER_SIZE 50

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GCG_VarHistoryBuffer GCG_VARHISTORYBUFFER;

struct GCG_VarHistoryBuffer {
   SCIP_VAR*             vars[GCG_VARHISTORYBUFFER_SIZE]; /**< variables */
   int                   nvars;              /**< number of variables */
   GCG_VARHISTORYBUFFER* next;               /**< next buffer */
   int                   nuses;              /**< number of uses */
};

typedef struct GCG_VarHistoryPointer GCG_VARHISTORYPOINTER;

struct GCG_VarHistoryPointer {
   GCG_VARHISTORYBUFFER* buffer;             /**< buffer */
   int                   pos;                /**< position in the buffer */
};

/** check if there is a next history event */
GCG_EXPORT
SCIP_Bool GCGvarhistoryHasNext(
   GCG_VARHISTORYPOINTER* pointer            /**< pointer to the history */
   );

/** get the next history event */
GCG_EXPORT
SCIP_RETCODE GCGvarhistoryNext(
   SCIP*                  scip,              /**< SCIP data structure */
   GCG_VARHISTORYPOINTER**pointer            /**< pointer to the history */
   );

/** add variable to history */
GCG_EXPORT
SCIP_RETCODE GCGvarhistoryAddVar(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORYPOINTER* pointer,            /**< pointer to the history */
   SCIP_VAR*              var                 /**< variable */
   );

/** capture a reference to the history */
GCG_EXPORT
SCIP_RETCODE GCGvarhistoryCaptureBuffer(
   GCG_VARHISTORYBUFFER*  buffer              /**< buffer */
   );

/** release a reference to the history */
GCG_EXPORT
SCIP_RETCODE GCGvarhistoryReleaseBuffer(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORYBUFFER** buffer              /**< buffer */
   );


#ifdef __cplusplus
}
#endif

#endif
