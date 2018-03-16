/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       */
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

/**@file   pricingjob.c
 * @brief  methods for working with pricing jobs
 * @author Christian Puchert
 *
 * Various methods to work with pricing jobs
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "pricingjob.h"
#include "pub_pricingjob.h"
#include "pub_pricingprob.h"

#include "gcg.h"

#include "scip/scip.h"

#include <assert.h>

/** create a pricing job */
SCIP_RETCODE GCGpricingjobCreate(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGJOB**      pricingjob,         /**< pricing job to be created */
   GCG_PRICINGPROB*      pricingprob,        /**< data structure of the corresponding pricing problem */
   GCG_SOLVER*           solver,             /**< pricing solver responsible for the pricing job */
   int                   chunk               /**< chunk that the pricing problem should belong to */
)
{
   SCIP_CALL( SCIPallocMemory(scip, pricingjob) );

   (*pricingjob)->pricingprob = pricingprob;
   (*pricingjob)->solver = solver;
   (*pricingjob)->chunk = chunk;
   (*pricingjob)->score = 0.0;
   (*pricingjob)->heuristic = FALSE;
   (*pricingjob)->pricingstatus = SCIP_STATUS_UNKNOWN;
   (*pricingjob)->nheuriters = 0;

   return SCIP_OKAY;
}

/** free a pricing job */
void GCGpricingjobFree(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGJOB**      pricingjob          /**< pricing job to be freed */
)
{
   SCIPfreeMemory(scip, pricingjob);
   *pricingjob = NULL;
}

/** setup a pricing job at the beginning of the pricing loop */
SCIP_RETCODE GCGpricingjobSetup(
   SCIP*                 scip,               /**< master SCIP instance */
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   SCIP_Bool             heuristic,          /**< shall the pricing job be performed heuristically? */
   int                   scoring,            /**< scoring parameter */
   int                   nroundscol,         /**< number of previous pricing rounds for which the number of improving columns should be counted */
   SCIP_Real             dualsolconv,        /**< dual solution value of corresponding convexity constraint */
   int                   npointsprob,        /**< total number of extreme points generated so far by the pricing problem */
   int                   nraysprob           /**< total number of extreme rays generated so far by the pricing problem */
   )
{
   int i;

   /* There should be no remaining columns from the previous iteration */
   assert(pricingjob->ncols == 0);
   assert(pricingjob->nimpcols == 0);

   pricingjob->heuristic = heuristic;

   /* set the score; the larger, the better */
   switch( scoring )
   {
   case 'i':
      pricingjob->score = (SCIP_Real) pricingjob->probnr;
      break;
   case 'd':
      pricingjob->score = dualsolconv;
      break;
   case 'r':
      pricingjob->score = 0.2 * npointsprob + nraysprob;
      break;
   case 'l':
      pricingjob->score = 0.0;
      for( i = 0; i < nroundscol; ++i )
         pricingjob->score += (SCIP_Real) pricingjob->ncolsround[i];
      break;
   default:
      pricingjob->score = 0.0;
      break;
   }

   /* initialize result variables */
   pricingjob->pricingstatus = SCIP_STATUS_UNKNOWN;
   pricingjob->nheuriters = 0;

   return SCIP_OKAY;
}

/** update a pricing job after the pricing problem has been solved */
SCIP_RETCODE GCGpricingjobUpdate(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   SCIP_STATUS           status,             /**< status after solving the pricing problem */
   SCIP_Real             lowerbound,         /**< lower bound returned by the pricing problem */
   GCG_COL**             cols,               /**< columns found by the last solving of the pricing problem */
   int                   ncols               /**< number of columns found */
   )
{
   int i;
   int j;
   int k;

   pricingjob->pricingstatus = status;
   pricingjob->lowerbound = lowerbound;

   if( pricingjob->colssize < pricingjob->ncols + ncols )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &pricingjob->cols, pricingjob->ncols + ncols) );
      pricingjob->colssize = pricingjob->ncols + ncols;
      for( i = pricingjob->ncols; i < pricingjob->colssize; ++i )
         pricingjob->cols[i] = NULL;
   }

   /* add new columns; ensure that the column array remains sorted by reduced costs */
   for( i = pricingjob->ncols + ncols - 1, j = pricingjob->ncols-1, k = ncols-1; k >= 0; --i )
   {
      if( j >= 0 && SCIPisDualfeasGT(scip, GCGcolGetRedcost(pricingjob->cols[j]), GCGcolGetRedcost(cols[k])) )
      {
         pricingjob->cols[i] = pricingjob->cols[j];
         --j;
      }
      else
      {
         if( SCIPisDualfeasNegative(scip, GCGcolGetRedcost(cols[k])) )
            ++pricingjob->nimpcols;

         pricingjob->cols[i] = cols[k];
         --k;
      }
   }

   pricingjob->ncols += ncols;

   return SCIP_OKAY;
}

/** update solving statistics of a pricing job */
void GCGpricingjobUpdateSolvingStats(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   ++pricingjob->nsolves;
   if( pricingjob->heuristic )
      ++pricingjob->nheuriters;
}

/** get the pricing solver with which the pricing job is to be performed */
GCG_SOLVER* GCGpricingjobGetSolver(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   return pricingjob->solver;
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

/** set the pricing job to be performed heuristically */
void GCGpricingjobSetHeuristic(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   pricingjob->heuristic = TRUE;
}

/** set the pricing job to be performed exactly */
void GCGpricingjobSetExact(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   pricingjob->heuristic = FALSE;
}

/* get the status of a pricing job */
SCIP_STATUS GCGpricingjobGetStatus(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   return pricingjob->pricingstatus;
}

/* set the status of a pricing job */
void GCGpricingjobSetStatus(
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   SCIP_STATUS           status              /**< new pricing status */
   )
{
   pricingjob->pricingstatus = status;
}

/* get the number of heuristic pricing iterations of the pricing job */
int GCGpricingjobGetNHeurIters(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   return pricingjob->nheuriters;
}
