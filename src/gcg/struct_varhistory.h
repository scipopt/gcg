/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file    struct_varhistory.h
 * @ingroup DATASTRUCTURES
 * @brief   data structures for managing variable history
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_VARHISTORY_H_
#define GCG_STRUCT_VARHISTORY_H_

#include <scip/type_cons.h>
#include <scip/type_lp.h>
#include <scip/type_var.h>
#include "type_varhistory.h"

#ifdef __cplusplus
extern "C" {
#endif

#define GCG_VARHISTORYBUFFER_SIZE 50

struct GCG_VarHistoryBuffer {
   int                   nvars;              /**< number of variables */
   GCG_VARHISTORYBUFFER* next;               /**< next buffer */
   int                   nuses;              /**< number of uses */
   SCIP_VAR*             vars[GCG_VARHISTORYBUFFER_SIZE]; /**< variables */
};

struct GCG_VarHistory {
   GCG_VARHISTORYBUFFER* buffer;             /**< buffer */
   int                   pos;                /**< position in the buffer */
};

#ifdef __cplusplus
}
#endif

#endif
