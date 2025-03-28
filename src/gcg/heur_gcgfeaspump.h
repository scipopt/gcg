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

/**@file   heur_gcgfeaspump.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Objective Feasibility Pump 2.0
 * @author Timo Berthold
 * @author Domenico Salvagnin
 *
 * The fundamental idea of the Feasibility Pump is to construct two sequences of points which hopefully converge to a
 * feasible solution. One sequence consists of LP-feasiblepoints, the other one of integer feasible points.  They are
 * produced by alternately rounding an LP-feasible point and solvng an LP that finds a point on the LP polyhedron which
 * is closest to the rounded, integral point (w.r.t. Manhattan distance).
 *
 * The version implemented in SCIP supports using an Objective Feasibility Pump that uses a convex combination of the
 * Manhattan distance and the original LP objective for reoptimization. It further features Feasibility Pump 2.0
 * capabilities, hence propagating the fixings for a faster convergence.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_GCGFEASPUMP_H__
#define __SCIP_HEUR_GCGFEASPUMP_H__


#include "scip/scip.h"
#include "gcg/gcg.h"


#ifdef __cplusplus
extern "C" {
#endif

/** creates the feaspump primal heuristic and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeHeurGcgfeaspump(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
