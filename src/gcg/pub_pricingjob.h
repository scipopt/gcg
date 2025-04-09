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

/**@file   pub_pricingjob.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for working with pricing jobs
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PUB_PRICINGJOB_H__
#define GCG_PUB_PRICINGJOB_H__


#include "gcg/type_pricingjob.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * GCG Pricing Job
 */

/**
 * @ingroup PRICINGJOB
 *
 * @{
 */


/** get the pricing problem structure associated with a pricing job */
GCG_EXPORT
GCG_PRICINGPROB* GCGpricingjobGetPricingprob(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** get the pricing solver with which the pricing job is to be performed */
GCG_EXPORT
GCG_SOLVER* GCGpricingjobGetSolver(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** get the chunk of a pricing job */
GCG_EXPORT
SCIP_Real GCGpricingjobGetChunk(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** get the score of a pricing job */
GCG_EXPORT
SCIP_Real GCGpricingjobGetScore(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** return whether the pricing job is to be performed heuristically */
GCG_EXPORT
SCIP_Bool GCGpricingjobIsHeuristic(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** get the number of heuristic pricing iterations of the pricing job */
GCG_EXPORT
int GCGpricingjobGetNHeurIters(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** returns TRUE iff the solver was changed after the last solver call */
GCG_EXPORT
SCIP_Bool GCGpricingjobSolverChanged(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** inform the pricing job that the current solver was called */
GCG_EXPORT
void GCGpricingjobSolverCalled(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/**@} */
#ifdef __cplusplus
}
#endif
#endif
