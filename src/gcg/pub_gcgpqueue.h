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

/**@file   pub_decomp.h
 * @ingroup PUBLICCOREAPI
 * @ingroup DATASTRUCTURES
 * @brief  public methods for working with priority queues
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PUB_GCGPQUEUE_H__
#define GCG_PUB_GCGPQUEUE_H__


#include "scip/type_scip.h"
#include "scip/type_retcode.h"
#include "scip/type_var.h"
#include "scip/type_cons.h"
#include "scip/type_misc.h"
#include "gcg/type_locks.h"
#include "gcg/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Priority Queue
 */

/**@defgroup PriorityQueue Priority Queue
 * @ingroup DATASTRUCTURES
 * @{
 */

/** creates priority queue */
GCG_EXPORT
SCIP_RETCODE GCGpqueueCreate(
   SCIP*                scip,                /** SCIP data structure */
   GCG_PQUEUE**         pqueue,              /**< pointer to a priority queue */
   int                   initsize,           /**< initial number of available element slots */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   GCG_LOCK*             memorylock          /**< memory lock */
   );

/** frees priority queue, but not the data elements themselves */
GCG_EXPORT
void GCGpqueueFree(
   GCG_PQUEUE**         pqueue              /**< pointer to a priority queue */
   );

/** clears the priority queue, but doesn't free the data elements themselves */
GCG_EXPORT
void GCGpqueueClear(
   GCG_PQUEUE*          pqueue              /**< priority queue */
   );

/** inserts element into priority queue */
GCG_EXPORT
SCIP_RETCODE GCGpqueueInsert(
   GCG_PQUEUE*          pqueue,             /**< priority queue */
   void*                 elem                /**< element to be inserted */
   );

/** removes and returns best element from the priority queue */
GCG_EXPORT
void* GCGpqueueRemove(
   GCG_PQUEUE*          pqueue              /**< priority queue */
   );

/** resorts priority queue after changing the key values */
GCG_EXPORT
SCIP_RETCODE GCGpqueueResort(
   GCG_PQUEUE*           pqueue              /**< priority queue */
   );

/** set the comperator of the priority queue */
GCG_EXPORT
SCIP_RETCODE GCGpqueueSetComperator(
   GCG_PQUEUE*           pqueue,             /**< priority queue */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/**< delete item at position pos and insert last item at this position and resort pqueue */
GCG_EXPORT
SCIP_RETCODE GCGpqueueDelete(
   GCG_PQUEUE*          pqueue,             /**< priority queue */
   int                  pos,                /**< position of item that should be deleted */
   void**               elem                /**< pointer to store element that was deleted from pqueue */
);

/** returns the best element of the queue without removing it */
GCG_EXPORT
void* GCGpqueueFirst(
   GCG_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the number of elements in the queue */
GCG_EXPORT
int GCGpqueueNElems(
   GCG_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the elements of the queue; changing the returned array may destroy the queue's ordering! */
GCG_EXPORT
void** GCGpqueueElems(
   GCG_PQUEUE*          pqueue              /**< priority queue */
   );

/**@} */


#ifdef __cplusplus
}
#endif
#endif
