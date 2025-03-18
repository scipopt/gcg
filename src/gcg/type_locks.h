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

/**@file    type_locks.h
 * @ingroup TYPEDEFINITIONS
 * @brief   type definitions for locks data structure
 * @author  Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_TYPE_GCGLOCKS_H__
#define GCG_TYPE_GCGLOCKS_H__

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _OPENMP

#define GCG_LOCK omp_nest_lock_t

typedef struct GCG_Locks GCG_LOCKS; /**< data structure to store OpenMP locks */

#define GCG_SET_LOCK(lockptr) omp_set_nest_lock(lockptr)

#define GCG_UNSET_LOCK(lockptr) omp_unset_nest_lock(lockptr)

#define GCG_INIT_LOCK(lockptr) omp_init_nest_lock(lockptr)

#define GCG_DESTROY_LOCK(lockptr) omp_destroy_nest_lock(lockptr)

#else

#define GCG_LOCK void

typedef void GCG_LOCKS;

#define GCG_SET_LOCK(lockname) 0

#define GCG_UNSET_LOCK(lockname) 0

#define GCG_INIT_LOCK(lockptr) 0

#define GCG_DESTROY_LOCK(lockptr) 0

#endif

#ifdef __cplusplus
}
#endif

#endif
