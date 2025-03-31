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

/**@file   type_origdiving.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for GCG diving heuristics on the original variables
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_ORIGDIVING_H__
#define __SCIP_TYPE_ORIGDIVING_H__

#include "scip/type_scip.h"
#include "gcg/type_gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GCG_DivingData GCG_DIVINGDATA;   /**< locally defined diving data */


/** destructor of diving heuristic to free user data (called when GCG is exiting)
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - heur            : the diving heuristic itself
 */
#define GCG_DECL_DIVINGFREE(x) SCIP_RETCODE x (GCG* gcg, SCIP_HEUR* heur)

/** initialization method of diving heuristic (called after problem was transformed)
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - heur            : the diving heuristic itself
 */
#define GCG_DECL_DIVINGINIT(x) SCIP_RETCODE x (GCG* gcg, SCIP_HEUR* heur)

/** deinitialization method of diving heuristic (called before transformed problem is freed)
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - heur            : the diving heuristic itself
 */
#define GCG_DECL_DIVINGEXIT(x) SCIP_RETCODE x (GCG* gcg, SCIP_HEUR* heur)

/** solving process initialization method of diving heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The diving heuristic may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - heur            : the diving heuristic itself
 */
#define GCG_DECL_DIVINGINITSOL(x) SCIP_RETCODE x (GCG* gcg, SCIP_HEUR* heur)

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The diving heuristic should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - heur            : the diving heuristic itself
 */
#define GCG_DECL_DIVINGEXITSOL(x) SCIP_RETCODE x (GCG* gcg, SCIP_HEUR* heur)

/** execution initialization method of diving heuristic (called when execution of diving heuristic is about to begin)
 *
 *  This method is called when the execution of the diving heuristic starts, before the diving loop.
 *  The diving heuristic may use this call to collect data which is specific to this call of the heuristic.
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - heur            : the diving heuristic itself
 */
#define GCG_DECL_DIVINGINITEXEC(x) SCIP_RETCODE x (GCG* gcg, SCIP_HEUR* heur)

/** execution deinitialization method of diving heuristic (called when execution data is freed)
 *
 *  This method is called before the execution of the heuristic stops.
 *  The diving heuristic should use this call to clean up its execution specific data.
 *
 *  input:
 *  - gcg             : GCG data structure
 *  - heur            : the diving heuristic itself
 */
#define GCG_DECL_DIVINGEXITEXEC(x) SCIP_RETCODE x (GCG* gcg, SCIP_HEUR* heur)

/** variable selection method of diving heuristic
 *
 *  Selects a master variable to dive on
 *
 *  input:
 *  - gcg              : GCG data structure
 *  - heur             : the diving heuristic itself
 *  - tabulist         : an array containing variables that must not be chosen
 *  - tabulistsize     : the size of the array
 *  - bestcand         : pointer to store the SCIP_VAR* returned by the selection rule
 *  - bestcandmayround : pointer to store whether the variable may be rounded without losing LP feasibility
 *  - bestcandroundup  : pointer to store whether the variable is to be rounded up
 */
#define GCG_DECL_DIVINGSELECTVAR(x) SCIP_RETCODE x (GCG* gcg, SCIP_HEUR* heur, SCIP_VAR** tabulist, int tabulistsize, SCIP_VAR** bestcand, SCIP_Bool* bestcandmayround, SCIP_Bool* bestcandroundup)

#ifdef __cplusplus
}
#endif

#endif
