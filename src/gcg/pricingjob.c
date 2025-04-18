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

/**@file   pricingjob.c
 * @brief  methods for working with pricing jobs
 * @author Christian Puchert
 *
 * Various methods to work with pricing jobs
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/pricingjob.h"
#include "gcg/pub_gcgcol.h"
#include "gcg/pub_pricingjob.h"

#include "gcg/gcg.h"
#include "gcg/pricer_gcg.h"
#include "gcg/pub_pricingprob.h"
#include "gcg/pub_solver.h"

#include "scip/scip.h"


/** create a pricing job */
SCIP_RETCODE GCGpricingjobCreate(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGJOB**      pricingjob,         /**< pricing job to be created */
   GCG_PRICINGPROB*      pricingprob,        /**< data structure of the corresponding pricing problem */
   GCG_SOLVER*           solver,             /**< pricing solver responsible for the pricing job */
   int                   chunk               /**< chunk that the pricing problem should belong to */
)
{
   SCIP* masterprob = GCGgetMasterprob(gcg);
   SCIP_CALL( SCIPallocBlockMemory(masterprob, pricingjob) );

   (*pricingjob)->pricingprob = pricingprob;
   (*pricingjob)->solver = solver;
   (*pricingjob)->chunk = chunk;
   (*pricingjob)->score = 0.0;
   (*pricingjob)->heuristic = FALSE;
   (*pricingjob)->nheuriters = 0;
   (*pricingjob)->solverchanged = TRUE;

   return SCIP_OKAY;
}

/** free a pricing job */
void GCGpricingjobFree(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGJOB**      pricingjob          /**< pricing job to be freed */
)
{
   SCIP* masterprob = GCGgetMasterprob(gcg);
   SCIPfreeBlockMemory(masterprob, pricingjob);
   *pricingjob = NULL;
}

/** setup a pricing job at the beginning of the pricing loop
 *
 *  The solver in the pricing job is reset during this call. This is to ensure that a valid solver is available for the
 *  pricing job. In some cases, such as when columns could not be found from a particular solver, the solver could be
 *  set to NULL. As such, the highest priority solver is assigned to the pricing job in this function.
 */
SCIP_RETCODE GCGpricingjobSetup(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   SCIP_Bool             heuristic,          /**< is heuristic pricing still possible at the current solving stage */
   int                   scoring,            /**< scoring parameter */
   int                   nroundscol,         /**< number of previous pricing rounds for which the number of improving columns should be counted */
   SCIP_Real             dualsolconv,        /**< dual solution value of corresponding convexity constraint */
   int                   npointsprob,        /**< total number of extreme points generated so far by the pricing problem */
   int                   nraysprob           /**< total number of extreme rays generated so far by the pricing problem */
   )
{
   GCG_PRICINGPROB* pricingprob = GCGpricingjobGetPricingprob(pricingjob);
   assert(pricingprob != NULL);

   /* set the score; the larger, the better */
   switch( scoring )
   {
   case 'i':
      pricingjob->score = - (SCIP_Real) GCGpricingprobGetProbnr(pricingprob);
      break;
   case 'd':
      pricingjob->score = dualsolconv;
      break;
   case 'r':
      pricingjob->score = -(0.2 * npointsprob + nraysprob);
      break;
   case 'l':
      pricingjob->score = (SCIP_Real) GCGpricingprobGetNColsLastRounds(pricingprob, nroundscol);
      break;
   default:
      pricingjob->score = 0.0;
      break;
   }

   GCGpricingjobResetSolver(gcg, pricingjob);
   if( heuristic && pricingjob->solver != NULL && GCGsolverIsHeurEnabled(pricingjob->solver) )
      GCGpricingjobResetHeuristic(pricingjob);
   else
      GCGpricingjobSetExact(pricingjob);

   return SCIP_OKAY;
}

/** get the pricing problem structure associated with a pricing job */
GCG_PRICINGPROB* GCGpricingjobGetPricingprob(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   return pricingjob->pricingprob;
}
/** get the pricing solver with which the pricing job is to be performed */
GCG_SOLVER* GCGpricingjobGetSolver(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   return pricingjob->solver;
}

/** reset the pricing solver to be used to the one with the highest priority */
void GCGpricingjobResetSolver(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   GCG_SOLVER** solvers;
   int nsolvers;
   int i;

   solvers = GCGpricerGetSolvers(gcg);
   nsolvers = GCGpricerGetNSolvers(gcg);

   /* get first available solver;
    * assumption: solvers are sorted by priority
    */
   pricingjob->solver = NULL;
   for( i = 0; i < nsolvers; ++i )
   {
      if( GCGsolverIsHeurEnabled(solvers[i]) || GCGsolverIsExactEnabled(solvers[i]) )
      {
         pricingjob->solver = solvers[i];
         break;
      }
   }
   pricingjob->solverchanged = TRUE;

   assert(pricingjob->solver != NULL);
}

/** get the next pricing solver to be used, or NULL of there is none */
void GCGpricingjobNextSolver(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   GCG_SOLVER** solvers;
   int nsolvers;

   int pos;
   int i;

   solvers = GCGpricerGetSolvers(gcg);
   nsolvers = GCGpricerGetNSolvers(gcg);

   /* get position of current solver */
   for( pos = 0; pos < nsolvers; ++pos )
      if( solvers[pos] == pricingjob->solver )
         break;
   assert(pos < nsolvers);

   /* get next available solver;
    * assumption: solvers are sorted by priority
    */
   pricingjob->solver = NULL;
   for( i = pos + 1; i < nsolvers; ++i )
   {
      if( GCGsolverIsHeurEnabled(solvers[i]) || GCGsolverIsExactEnabled(solvers[i]) )
      {
         pricingjob->solver = solvers[i];
         break;
      }
   }
   pricingjob->solverchanged = TRUE;
}

/** returns TRUE iff the solver was changed after the last solver call */
SCIP_Bool GCGpricingjobSolverChanged(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   return pricingjob->solverchanged;
}

/** inform the pricing job that the current solver was called */
void GCGpricingjobSolverCalled(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   pricingjob->solverchanged = FALSE;
}

/** get the chunk of a pricing job */
SCIP_Real GCGpricingjobGetChunk(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   return pricingjob->chunk;
}

/** get the score of a pricing job */
SCIP_Real GCGpricingjobGetScore(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   return pricingjob->score;
}

/** return whether the pricing job is to be performed heuristically */
SCIP_Bool GCGpricingjobIsHeuristic(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   return pricingjob->heuristic;
}

/** set the pricing job to be performed exactly */
void GCGpricingjobSetExact(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   pricingjob->heuristic = FALSE;
}

/** reset number of heuristic pricing iterations of a pricing job */
void GCGpricingjobResetHeuristic(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   if( pricingjob->solver != NULL && GCGsolverIsHeurEnabled(pricingjob->solver) )
      pricingjob->heuristic = TRUE;
   else
      pricingjob->heuristic = FALSE;
   pricingjob->nheuriters = 0;
}

/** update number of heuristic pricing iterations of a pricing job */
void GCGpricingjobIncreaseNHeurIters(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   if( pricingjob->heuristic )
      ++pricingjob->nheuriters;
}

/** get the number of heuristic pricing iterations of the pricing job */
int GCGpricingjobGetNHeurIters(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   return pricingjob->nheuriters;
}
