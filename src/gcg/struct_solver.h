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

/**@file   struct_solver.h
 * @ingroup DATASTRUCTURES
 * @brief  data structures for solvers
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_SOLVER_H__
#define GCG_STRUCT_SOLVER_H__

#include "gcg/type_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

/** pricing problem solver data structure */
struct GCG_Solver
{
   char*                 name;               /**< solver name */
   char*                 desc;               /**< solver description */
   int                   priority;           /**< solver priority */
   SCIP_Bool             heurenabled;        /**< switch for heuristic solving method */
   SCIP_Bool             exactenabled;       /**< switch for exact solving method */
   GCG_SOLVERDATA*       solverdata;         /**< private solver data structure */

   GCG_DECL_SOLVERFREE((*solverfree));       /**< destruction method */
   GCG_DECL_SOLVERINIT((*solverinit));       /**< initialization method */
   GCG_DECL_SOLVEREXIT((*solverexit));       /**< deinitialization method */
   GCG_DECL_SOLVERINITSOL((*solverinitsol)); /**< solving process initialization method */
   GCG_DECL_SOLVEREXITSOL((*solverexitsol)); /**< solving process deinitialization method */
   GCG_DECL_SOLVERUPDATE((*solverupdate));   /**< update method */
   GCG_DECL_SOLVERSOLVE((*solversolve));     /**< solving callback method */
   GCG_DECL_SOLVERSOLVEHEUR((*solversolveheur)); /**< heuristic solving callback method */

   SCIP_CLOCK*           optfarkasclock;     /**< optimal farkas pricing time */
   SCIP_CLOCK*           optredcostclock;    /**< optimal reduced cost pricing time */
   SCIP_CLOCK*           heurfarkasclock;    /**< heuristic farkas pricing time */
   SCIP_CLOCK*           heurredcostclock;   /**< heuristic reduced cost pricing time */
   int                   optfarkascalls;     /**< optimal farkas pricing calls */
   int                   optredcostcalls;    /**< optimal reduced cost pricing calls */
   int                   heurfarkascalls;    /**< heuristic farkas pricing calls */
   int                   heurredcostcalls;   /**< heuristic reduced cost pricing calls */
};


#ifdef __cplusplus
}
#endif

#endif
