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

/**@file   pricingjob.h
 * @ingroup PRICING_PRIV
 * @brief  private methods for working with pricing jobs, to be used by the pricing controller only
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PRICINGJOB_H__
#define GCG_PRICINGJOB_H__


#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create a pricing job */
SCIP_RETCODE GCGpricingjobCreate(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGJOB**      pricingjob,         /**< pricing job to be created */
   GCG_PRICINGPROB*      pricingprob,        /**< data structure of the corresponding pricing problem */
   GCG_SOLVER*           solver,             /**< pricing solver responsible for the pricing job */
   int                   chunk               /**< chunk that the pricing problem should belong to */
);

/** free a pricing job */
void GCGpricingjobFree(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGJOB**      pricingjob          /**< pricing job to be freed */
);

/** setup a pricing job at the beginning of the pricing loop */
SCIP_RETCODE GCGpricingjobSetup(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   SCIP_Bool             heuristic,          /**< shall the pricing job be performed heuristically? */
   int                   scoring,            /**< scoring parameter */
   int                   nroundscol,         /**< number of previous pricing rounds for which the number of improving columns should be counted */
   SCIP_Real             dualsolconv,        /**< dual solution value of corresponding convexity constraint */
   int                   npointsprob,        /**< total number of extreme points generated so far by the pricing problem */
   int                   nraysprob           /**< total number of extreme rays generated so far by the pricing problem */
   );

/** reset the pricing solver to be used to the one with the highest priority */
void GCGpricingjobResetSolver(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** get the next pricing solver to be used, or NULL of there is none */
void GCGpricingjobNextSolver(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** set the pricing job to be performed exactly */
void GCGpricingjobSetExact(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** reset number of heuristic pricing iterations of a pricing job */
void GCGpricingjobResetHeuristic(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** update number of heuristic pricing iterations of a pricing job */
void GCGpricingjobIncreaseNHeurIters(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

#ifdef __cplusplus
}
#endif

#endif
