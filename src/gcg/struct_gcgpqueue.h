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

/**@file   struct_gcgpqueue.h
 * @ingroup DATASTRUCTURES
 * @brief  data structure for priority queue
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_GCGPQUEUE_H_
#define GCG_STRUCT_GCGPQUEUE_H_

#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/scip.h"

#include "gcg/type_gcgpqueue.h"
#include "gcg/type_locks.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @brief priority queue data structure
 *
 *  Elements are stored in an array, which grows dynamically in size as new elements are added to the queue.
 *  The ordering is done through a pointer comparison function.
 *  The array is organized as follows. The root element (that is the "best" element $r$ with $r <= x$ for all $x$)
 *  is stored in position 0. The children of an element at position $p$ are stored at positions $q_1 = 2*p+1$ and
 *  $q_2 = 2*p+2$. That means, the parent of the element at position $q$ is at position $p = (q-1)/2$.
 *  At any time, the condition holds that $p <= q$ for each parent $p$ and its children $q$.
 *  Insertion and removal of single elements needs time $O(log n)$.
 */
struct GCG_PQueue
{
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp));        /**< compares two data elements */
   void**                slots;              /**< array of element slots */
   int                   len;                /**< number of used element slots */
   int                   size;               /**< total number of available element slots */
#ifdef _OPENMP
   GCG_LOCK*             memorylock;         /**< memory lock */
#endif
};

#ifdef __cplusplus
}
#endif

#endif /* STRUCT_GCGPQUEUE_H_ */
