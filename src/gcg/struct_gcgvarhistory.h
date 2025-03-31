/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file    struct_gcgvarhistory.h
 * @ingroup DATASTRUCTURES
 * @brief   data structures for managing variable history
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_GCGVARHISTORY_H_
#define GCG_STRUCT_GCGVARHISTORY_H_

#include <scip/type_cons.h>
#include <scip/type_lp.h>
#include <scip/type_var.h>
#include "gcg/type_gcgvarhistory.h"

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
