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

/**@file    solver.h
 * @ingroup PUBLICCOREAPI
 * @brief   public methods for GCG pricing solvers
 * @author  Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_PUB_SOLVER_H_
#define GCG_PUB_SOLVER_H_


#include "gcg/type_solver.h"
#include "scip/scip.h"
#include "gcg/def.h"
#include "gcg/type_gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @ingroup PRICING_PUB
 * @{
 */

/** compares two solvers w. r. t. their priorities */
GCG_EXPORT
SCIP_DECL_SORTPTRCOMP(GCGsolverComp);

/** gets user data of GCG pricing solver */
GCG_EXPORT
GCG_SOLVERDATA* GCGsolverGetData(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** sets user data of GCG pricing solver */
GCG_EXPORT
void GCGsolverSetData(
   GCG_SOLVER*           solver,             /**< pricing solver */
   GCG_SOLVERDATA*       solverdata          /**< pricing solver data */
   );

/** gets name of GCG pricing solver */
GCG_EXPORT
const char* GCGsolverGetName(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets description of GCG pricing solver */
GCG_EXPORT
const char* GCGsolverGetDesc(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets priority of GCG pricing solver */
GCG_EXPORT
int GCGsolverGetPriority(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets whether heuristic solving method of GCG pricing solver is enabled */
GCG_EXPORT
SCIP_Bool GCGsolverIsHeurEnabled(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets whether exact solving method of GCG pricing solver is enabled */
GCG_EXPORT
SCIP_Bool GCGsolverIsExactEnabled(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets number of exact Farkas pricing calls of pricing solver */
GCG_EXPORT
int GCGsolverGetOptFarkasCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets number of exact reduced cost pricing calls of pricing solver */
GCG_EXPORT
int GCGsolverGetOptRedcostCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets number of heuristic Farkas pricing calls of pricing solver */
GCG_EXPORT
int GCGsolverGetHeurFarkasCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets number of heuristic reduced cost pricing calls of pricing solver */
GCG_EXPORT
int GCGsolverGetHeurRedcostCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets exact Farkas pricing time of pricing solver */
GCG_EXPORT
SCIP_Real GCGsolverGetOptFarkasTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets exact reduced cost pricing time of pricing solver */
GCG_EXPORT
SCIP_Real GCGsolverGetOptRedcostTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets heuristic Farkas pricing time of pricing solver */
GCG_EXPORT
SCIP_Real GCGsolverGetHeurFarkasTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets heuristic reduced cost pricing time of pricing solver */
GCG_EXPORT
SCIP_Real GCGsolverGetHeurRedcostTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** @} */
#ifdef __cplusplus
}

#endif

#endif
