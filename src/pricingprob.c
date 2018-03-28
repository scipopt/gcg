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
   int                   colssize,           /**< size of column array */
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
   (*pricingprob)->nextconsidx = 0;
   (*pricingprob)->pricingstatus = SCIP_STATUS_UNKNOWN;
   (*pricingprob)->lowerbound = -SCIPinfinity(scip);
   (*pricingprob)->colssize = colssize;
   (*pricingprob)->ncols = 0;
   (*pricingprob)->nimpcols = 0;
   (*pricingprob)->nsolves  = 0;

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &(*pricingprob)->cols, colssize) );
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
   SCIPfreeMemoryArray(scip, &(*pricingprob)->cols);
   SCIPfreeMemory(scip, pricingprob);
   *pricingprob = NULL;
}

/** initialize pricing problem at the beginning of the pricing round */
void GCGpricingprobInitPricing(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   assert(pricingprob->ncols == 0);
   assert(pricingprob->nimpcols == 0);

   pricingprob->nbranchconss = 0;
   pricingprob->nextconsidx = 0;
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
   ++pricingprob->nextconsidx;

   assert(pricingprob->nbranchconss == pricingprob->lastconsidx);

   return SCIP_OKAY;
}

/** reset the pricing problem statistics for the current pricing round */
void GCGpricingprobReset(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   assert(pricingprob->ncols == 0);
   assert(pricingprob->nimpcols == 0);

   pricingprob->nextconsidx = pricingprob->nbranchconss;
   pricingprob->pricingstatus = SCIP_STATUS_UNKNOWN;
   pricingprob->lowerbound = -SCIPinfinity(scip);
   pricingprob->nsolves = 0;
}

/** update solution information of a pricing problem */
void GCGpricingprobUpdate(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   SCIP_STATUS           status,             /**< new pricing status */
   SCIP_Real             lowerbound,         /**< new lower bound */
   GCG_COL**             cols,               /**< sorted array of columns found by the last solver call */
   int                   ncols               /**< number of found columns */
   )
{
   int i;
   int j;
   int pos;
   int nnewcols;

   /* update status and lower bound */
   pricingprob->pricingstatus = status;
   if( SCIPisDualfeasGT(scip, lowerbound, pricingprob->lowerbound) )
      pricingprob->lowerbound = lowerbound;

   /*
    * ensure that no columns are stored double:
    * for each new column, check whether it is equal to a column in 'pricingprob->cols';
    * if this is the case, free it
    */
   nnewcols = ncols;
   for( j = ncols-1, pos = pricingprob->ncols-1; j >= 0; --j )
   {
      SCIP_Bool equal = FALSE;

      /* columns with higher reduced cost cannot be equal */
      for( ; pos >= 0 && SCIPisDualfeasGT(scip, GCGcolGetRedcost(pricingprob->cols[pos]), GCGcolGetRedcost(cols[j])); --pos );

      /* among columns with equal reduced cost, check which are equal to the given one */
      for( i = pos; i >= 0 && SCIPisDualfeasEQ(scip, GCGcolGetRedcost(pricingprob->cols[i]), GCGcolGetRedcost(cols[j])) & !equal; --i )
         if( GCGcolIsEq(pricingprob->cols[i], cols[j]) )
            equal = TRUE;

      if( equal )
      {
         GCGfreeGcgCol(&cols[j]);
         cols[j] = NULL;
         --nnewcols;
      }
   }

   /* add new columns, i.e. sort columns given in 'cols' into 'pricingprob->cols';
    * ensure that the array 'pricingprob->cols' remains sorted by reduced costs;
    * only the best 'pricingprob->ncols' are stored, the rest is freed
    */
   for( pos = pricingprob->ncols + nnewcols - 1, i = pricingprob->ncols - 1, j = ncols - 1; j >= 0; )
   {
      if( cols[j] == NULL )
      {
         --j;
         continue;
      }

      /* case 1: keep current column from 'pricingprob->cols' */
      if( i >= 0 && SCIPisDualfeasGT(scip, GCGcolGetRedcost(pricingprob->cols[i]), GCGcolGetRedcost(cols[j])) )
      {
         if( pos < pricingprob->colssize )
            pricingprob->cols[pos] = pricingprob->cols[i];
         else
         {
            if( SCIPisDualfeasNegative(scip, GCGcolGetRedcost(pricingprob->cols[i])) )
               --pricingprob->nimpcols;
            assert(pricingprob->nimpcols >= 0);
            GCGfreeGcgCol(&pricingprob->cols[i]);
         }

         --i;
         --pos;
      }
      /* case 2: add current column from 'cols' to 'pricingprob->cols' */
      else
      {
         if( pos < pricingprob->colssize )
         {
            if( SCIPisDualfeasNegative(scip, GCGcolGetRedcost(cols[j])) )
               ++pricingprob->nimpcols;
            pricingprob->cols[pos] = cols[j];
         }
         else
            GCGfreeGcgCol(&cols[j]);

         --j;
         --pos;
      }
   }

   pricingprob->ncols = MIN(pricingprob->ncols + nnewcols, pricingprob->colssize);

   ++pricingprob->nsolves;
}

/** for a pricing problem, move its columns to the pricing store or column pool */
SCIP_RETCODE GCGpricingprobMoveCols(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   GCG_COLPOOL*          colpool,            /**< GCG column pool */
   GCG_PRICESTORE*       pricestore,         /**< GCG pricing store */
   SCIP_Bool             usecolpool,         /**< use column pool? */
   SCIP_Bool             usepricestore       /**< use price store? */
   )
{
   SCIP_Bool added;

   assert(pricingprob->cols != NULL || pricingprob->ncols == 0);

   for( int i = 0; i < pricingprob->ncols; ++i )
   {
      SCIPdebugMessage("  (prob %d) column %d/%d <%p>: ", pricingprob->probnr, i+1, pricingprob->ncols, (void*) pricingprob->cols[i]);

      added = FALSE;

      if( usepricestore && SCIPisDualfeasNegative(scip, GCGcolGetRedcost(pricingprob->cols[i])) )
      {
         SCIP_CALL( GCGcomputeColMastercoefs(scip, pricingprob->cols[i]) );

         SCIP_CALL( GCGpricestoreAddCol(scip, pricestore, pricingprob->cols[i], FALSE) );
         added = TRUE;
      }
      else if( usecolpool )
      {
         SCIP_CALL( GCGcolpoolAddCol(colpool, pricingprob->cols[i], &added) );
      }

      if( !added )
      {
         GCGfreeGcgCol(&pricingprob->cols[i]);
         SCIPdebugPrintf("freed.\n");
      }
      else
      {
         SCIPdebugPrintf("added to column pool or price store.\n");
      }

      pricingprob->cols[i] = NULL;
   }

   pricingprob->ncols = 0;
   pricingprob->nimpcols = 0;

   return SCIP_OKAY;
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

/** get index of next generic branching constraint added to the pricing problem */
int GCGpricingprobGetNextConsIdx(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->nextconsidx;
}

/** decrease index of next generic branching constraint added to the pricing problem */
void GCGpricingprobDecreaseNextConsIdx(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   --pricingprob->nextconsidx;
}

/** get the status of a pricing problem */
SCIP_STATUS GCGpricingprobGetStatus(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->pricingstatus;
}

/** get the lower bound of a pricing problem */
SCIP_Real GCGpricingprobGetLowerbound(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->lowerbound;
}

/** get the best column found by solving a particular pricing problem */
GCG_COL* GCGpricingprobGetBestCol(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   if( pricingprob->ncols == 0 )
      return NULL;
   else
   {
      assert(pricingprob->cols[0] != NULL);
      return pricingprob->cols[0];
   }
}

/** get the best reduced cost of a pricing problem */
SCIP_Real GCGpricingprobGetBestRedcost(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   if( pricingprob->ncols == 0 )
      return 0.0;
   else
   {
      assert(pricingprob->cols[0] != NULL);
      return GCGcolGetRedcost(pricingprob->cols[0]);
   }
}

/** get the number of columns found for this pricing problem */
int GCGpricingprobGetNCols(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   return pricingprob->ncols;
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

/** update numbers of improving columns over the last pricing rounds */
void GCGpricingprobUpdateNColsround(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
   )
{
   int i;

   for( i = nroundscol-1; i > 0; --i )
      pricingprob->ncolsround[i] = pricingprob->ncolsround[i-1];
   pricingprob->ncolsround[0] = pricingprob->nimpcols;
}
