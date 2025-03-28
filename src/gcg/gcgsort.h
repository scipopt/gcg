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

/**@file   gcgsort.h
 * @brief  sorting functions, adapted from SCIP's sorttpl to include userdata
 * @author Tobias Oelschlegel
 **/

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_SORT_H__
#define GCG_SORT_H__

#include "gcg/def.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** compares two element indices
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 **/
#define GCG_DECL_SORTINDCOMP(x) int x (void* userdata, void* dataptr, int ind1, int ind2)

/** compares two data element pointers
 *  result:
 *    < 0: elem1 comes before (is better than) elem2
 *    = 0: both elements have the same value
 *    > 0: elem2 comes after (is worse than) elem2
 **/
#define GCG_DECL_SORTPTRCOMP(x) int x (void* userdata, void* elem1, void* elem2)

GCG_EXPORT
void GCGsortPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   GCG_DECL_SORTPTRCOMP((*ptrcomp)),         /**< data element comparator */
   void*                 userdata,           /**< userdata that is supplied to the comparator function */
   int                   len                 /**< length of array */
   );

/** sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order */
GCG_EXPORT
void GCGsortPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   GCG_DECL_SORTPTRCOMP((*ptrcomp)),         /**< data element comparator */
   void*                 userdata,           /**< userdata that is supplied to the comparator function */
   int                   len                 /**< length of arrays */
   );

#ifdef __cplusplus
}
#endif

#endif
