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

/**@file   pricingprob.c
 * @brief  methods for working with pricing problems
 * @author Christian Puchert
 *
 * Various methods to work with pricing problems
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "pricingprob.h"
#include "pub_pricingprob.h"

#include "gcg.h"

#include "scip/scip.h"

#include <assert.h>

/** create a pricing problem */
SCIP_RETCODE GCGpricingprobCreate(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGPROB**     pricingprob,        /**< pricing problem to be created */
   SCIP*                 pricingscip,        /**< SCIP data structure of the corresponding pricing problem */
   int                   probnr,             /**< index of the corresponding pricing problem */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
)
{
   SCIP_CALL( SCIPallocMemory(scip, pricingprob) );

   (*pricingprob)->pricingscip = pricingscip;
   (*pricingprob)->probnr = probnr;
   (*pricingprob)->nsolves = 0;
   (*pricingprob)->nheuriters = 0;
   (*pricingprob)->status = SCIP_STATUS_UNKNOWN;
   (*pricingprob)->nimpcols = 0;

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &(*pricingprob)->ncolsround, nroundscol) );


   return SCIP_OKAY;
}

/** free a pricing problem */
void GCGpricingprobFree(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGPROB**     pricingprob         /**< pricing problem to be freed */
)
{
   SCIPfreeMemoryArray(scip, &(*pricingprob)->ncolsround);
   SCIPfreeMemoryArray(scip, &(*pricingprob)->nsolvers);
   SCIPfreeMemory(scip, pricingprob);
   *pricingprob = NULL;
}

/** reset the pricing problem statistics for the current pricing round */
void GCGpricingprobReset(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   )
{
   pricingprob->nsolves = 0;
   pricingprob->nheuriters = 0;
   pricingprob->status = SCIP_STATUS_UNKNOWN;
   pricingprob->lowerbound = -SCIPinfinity(scip);
   pricingprob->nimpcols = 0;
}

/** update solution information of a pricing problem */
void GCGpricingprobUpdate(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   int                   nsolves,            /**< additional number of times the pricing problem was solved */
   SCIP_STATUS           status,             /**< new pricing status */
   SCIP_Real             lowerbound,         /**< new lower bound */
   int                   nimpcols            /**< additional number of found improving columns */
   )
{
   pricingprob->nsolves += nsolves;
   pricingprob->status = status;
   if( SCIPisDualfeasGT(scip, lowerbound, pricingprob->lowerbound) )
      pricingprob->lowerbound = lowerbound;
   pricingprob->nimpcols += nimpcols;
}

/** get the SCIP instance corresponding to the pricing problem */
SCIP* GCGpricingprobGetPricingscip(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->pricingscip;
}

/** get the index of the corresponding pricing problem */
int GCGpricingprobGetProbnr(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->probnr;
}

/** get the number of times the pricing problem was solved during the loop */
int GCGpricingprobGetNSolves(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->nsolves;
}

/* get the status of a pricing problem */
SCIP_STATUS GCGpricingprobGetStatus(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->pricingstatus;
}

/* get the lower bound of a pricing problem */
SCIP_Real GCGpricingprobGetLowerbound(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->lowerbound;
}

/* get the number of improving columns found for this pricing problem */
int GCGpricingprobGetNImpCols(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->nimpcols;
}

/* get the total number of improving colums found in the last pricing rounds */
int GCGpricingprobGetNColsLastRounds(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
   )
{
   int ncols;
   int i;

   ncols = 0;
   for( i = 0; i < nroundcols; ++i )
      ncols += pricingprob->ncolsround[i];

   return ncols;
}

/* update numbers of improving columns over the last pricing rounds */
void GCGpricingprobUpdateNColsround(
   GCG_PRICINGPROB*      pricingjob,         /**< pricing problem structure */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
   )
{
   int i;

   for( i = nroundscol-1; i > 0; --i )
      pricingprob->ncolsround[i] = pricingprob->ncolsround[i-1];
   pricingprob->ncolsround[0] = pricingprob->nimpcols;
}
