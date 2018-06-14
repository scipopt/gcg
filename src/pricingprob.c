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
#include "pricestore_gcg.h"
#include "pub_colpool.h"
#include "pub_gcgcol.h"
#include "pub_pricingjob.h"

#include "scip/scip.h"


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
   (*pricingprob)->branchconss = NULL;
   (*pricingprob)->branchduals = NULL;
   (*pricingprob)->nbranchconss = 0;
   (*pricingprob)->branchconsssize = 0;
   (*pricingprob)->branchconsidx = 0;
   (*pricingprob)->consisadded = TRUE;
   (*pricingprob)->status = GCG_PRICINGSTATUS_UNKNOWN;
   (*pricingprob)->lowerbound = -SCIPinfinity(scip);
   (*pricingprob)->nimpcols = 0;
   (*pricingprob)->nsolves  = 0;

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &(*pricingprob)->ncolsround, nroundscol) );


   return SCIP_OKAY;
}

/** free a pricing problem */
void GCGpricingprobFree(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGPROB**     pricingprob         /**< pricing problem to be freed */
)
{
   SCIPfreeMemoryArrayNull(scip, &(*pricingprob)->branchduals);
   SCIPfreeMemoryArrayNull(scip, &(*pricingprob)->branchconss);
   SCIPfreeMemoryArray(scip, &(*pricingprob)->ncolsround);
   SCIPfreeMemory(scip, pricingprob);
   *pricingprob = NULL;
}

/** initialize pricing problem at the beginning of the pricing round */
void GCGpricingprobInitPricing(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   assert(pricingprob->nimpcols == 0);

   pricingprob->nbranchconss = 0;
   pricingprob->branchconsidx = 0;
   pricingprob->consisadded = TRUE;
}

/** uninitialize pricing problem at the beginning of the pricing round */
void GCGpricingprobExitPricing(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
   )
{
   int i;

   for( i = nroundscol-1; i > 0; --i )
      pricingprob->ncolsround[i] = pricingprob->ncolsround[i-1];
   pricingprob->ncolsround[0] = pricingprob->nimpcols;

   pricingprob->nimpcols = 0;
}

/** add generic branching data (constraint and dual value) to the current pricing problem */
SCIP_RETCODE GCGpricingprobAddGenericBranchData(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   SCIP_CONS*            branchcons,         /**< generic branching constraint */
   SCIP_Real             branchdual          /**< corresponding dual solution value */
   )
{
   /* allocate memory, if necessary */
   if( pricingprob->branchconsssize == pricingprob->nbranchconss )
   {
      assert((pricingprob->branchconsssize == 0) == (pricingprob->branchconss == NULL));
      assert((pricingprob->branchconsssize == 0) == (pricingprob->branchduals == NULL));

      pricingprob->branchconsssize = SCIPcalcMemGrowSize(scip, pricingprob->branchconsssize+1);

      if( pricingprob->branchconss == NULL )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &pricingprob->branchconss, pricingprob->branchconsssize) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &pricingprob->branchduals, pricingprob->branchconsssize) );
      }
      else
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, &pricingprob->branchconss, pricingprob->branchconsssize) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &pricingprob->branchduals, pricingprob->branchconsssize) );
      }

   }

   /* add constraint and dual solution value */
   pricingprob->branchconss[pricingprob->nbranchconss] = branchcons;
   pricingprob->branchduals[pricingprob->nbranchconss] = branchdual;
   ++pricingprob->nbranchconss;
   ++pricingprob->branchconsidx;

   return SCIP_OKAY;
}

/** reset the pricing problem statistics for the current pricing round */
void GCGpricingprobReset(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   assert(pricingprob->nimpcols == 0);

   pricingprob->branchconsidx = pricingprob->nbranchconss;
   pricingprob->status = GCG_PRICINGSTATUS_UNKNOWN;
   pricingprob->lowerbound = -SCIPinfinity(scip);
   pricingprob->nsolves = 0;
}

/** update solution information of a pricing problem */
void GCGpricingprobUpdate(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   GCG_PRICINGSTATUS     status,             /**< status of last pricing job */
   SCIP_Real             lowerbound,         /**< new lower bound */
   int                   nimpcols            /**< number of new improving columns */
   )
{
   /* if the solver was not applicable to the problem, there is nothing to be done */
   if( status == GCG_PRICINGSTATUS_NOTAPPLICABLE )
      return;

   /* update status, lower bound and number of improving columns */
   pricingprob->status = status;
   if( SCIPisDualfeasGT(scip, lowerbound, pricingprob->lowerbound) )
      pricingprob->lowerbound = lowerbound;
   pricingprob->nimpcols += nimpcols;

   ++pricingprob->nsolves;
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

/** get generic branching data corresponding to the pricing problem */
void GCGpricingprobGetGenericBranchData(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   SCIP_CONS***          branchconss,        /**< pointer to store branching constraints array, or NULL */
   SCIP_Real**           branchduals,        /**< pointer to store array of corresponding dual values, or NULL */
   int*                  nbranchconss        /**< pointer to store number of generic branching constraints, or NULL */
   )
{
   if( branchconss != NULL )
      *branchconss = pricingprob->branchconss;
   if( branchduals != NULL )
      *branchduals = pricingprob->branchduals;
   if( nbranchconss != NULL )
      *nbranchconss = pricingprob->nbranchconss;
}

/** get the number of generic branching constraints corresponding to the pricing problem */
int GCGpricingprobGetNGenericBranchconss(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->nbranchconss;
}

/** get index of current generic branching constraint considered the pricing problem */
int GCGpricingprobGetBranchconsIdx(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->branchconsidx;
}

/** check if the current generic branching constraint has already been added */
SCIP_Bool GCGpricingprobBranchconsIsAdded(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->consisadded;
}

/** mark the current generic branching constraint to be added */
void GCGpricingprobMarkBranchconsAdded(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   pricingprob->consisadded = TRUE;
}

/** add the information that the next branching constraint must be added */
void GCGpricingprobNextBranchcons(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   assert(pricingprob->branchconsidx >= 1);
   --pricingprob->branchconsidx;
   pricingprob->consisadded = FALSE;
   pricingprob->status = GCG_PRICINGSTATUS_UNKNOWN;
}

/** get the status of a pricing problem */
GCG_PRICINGSTATUS GCGpricingprobGetStatus(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->status;
}

/** get the lower bound of a pricing problem */
SCIP_Real GCGpricingprobGetLowerbound(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->lowerbound;
}

/** get the number of improving columns found for this pricing problem */
int GCGpricingprobGetNImpCols(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->nimpcols;
}

/** get the number of times the pricing problem was solved during the loop */
int GCGpricingprobGetNSolves(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->nsolves;
}

/** get the total number of improving colums found in the last pricing rounds */
int GCGpricingprobGetNColsLastRounds(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
   )
{
   int ncols;
   int i;

   ncols = 0;
   for( i = 0; i < nroundscol; ++i )
      ncols += pricingprob->ncolsround[i];

   return ncols;
}
