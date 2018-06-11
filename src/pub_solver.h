/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file    solver.h
 * @ingroup PUBLICMETHODS
 * @brief   public methods for GCG pricing solvers
 * @author  Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_PUB_SOLVER_H_
#define GCG_PUB_SOLVER_H_

#include "type_solver.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** compares two solvers w. r. t. their priorities */
EXTERN
SCIP_DECL_SORTPTRCOMP(GCGsolverComp);

/** gets user data of GCG pricing solver */
EXTERN
GCG_SOLVERDATA* GCGsolverGetData(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** sets user data of GCG pricing solver */
EXTERN
void GCGsolverSetData(
   GCG_SOLVER*           solver,             /**< pricing solver */
   GCG_SOLVERDATA*       solverdata          /**< pricing solver data */
   );

/** gets name of GCG pricing solver */
EXTERN
const char* GCGsolverGetName(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets description of GCG pricing solver */
EXTERN
const char* GCGsolverGetDesc(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets priority of GCG pricing solver */
EXTERN
int GCGsolverGetPriority(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets whether GCG pricing solver is enabled */
EXTERN
SCIP_Bool GCGsolverIsEnabled(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets number of exact Farkas pricing calls of pricing solver */
EXTERN
int GCGsolverGetOptFarkasCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets number of exact reduced cost pricing calls of pricing solver */
EXTERN
int GCGsolverGetOptRedcostCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets number of heuristic Farkas pricing calls of pricing solver */
EXTERN
int GCGsolverGetHeurFarkasCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets number of heuristic reduced cost pricing calls of pricing solver */
EXTERN
int GCGsolverGetHeurRedcostCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets exact Farkas pricing time of pricing solver */
EXTERN
SCIP_Real GCGsolverGetOptFarkasTime(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets exact reduced cost pricing time of pricing solver */
EXTERN
SCIP_Real GCGsolverGetOptRedcostTime(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets heuristic Farkas pricing time of pricing solver */
EXTERN
SCIP_Real GCGsolverGetHeurFarkasTime(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

/** gets heuristic reduced cost pricing time of pricing solver */
EXTERN
SCIP_Real GCGsolverGetHeurRedcostTime(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_SOLVER*           solver              /**< pricing solver */
   );

#ifdef __cplusplus
}
#endif

#endif
