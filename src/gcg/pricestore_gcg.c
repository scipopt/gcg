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

/**@file   pricestore_gcg.c
 * @brief  methods for storing priced cols (based on SCIP's separation storage)
 * @author Jonas Witt
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/tree.h"
#include "scip/reopt.h"
#include "scip/event.h"
#include "scip/cons.h"
#include "scip/debug.h"

#include "gcg.h"
#include "pricestore_gcg.h"
#include "struct_pricestore_gcg.h"
#include "pricer_gcg.h"

#ifdef _OPENMP
#include "struct_locks.h"
#endif

/*
 * dynamic memory arrays
 */

/** resizes cols and score arrays to be able to store at least num entries */
static
SCIP_RETCODE pricestoreEnsureColsMem(
   GCG_PRICESTORE*       pricestore,          /**< price storage */
   int                   arrayindex,          /**< index of the arrays */
   int                   num                  /**< minimal number of slots in array */
   )
{
   int retcode;
   assert(pricestore != NULL);
   assert(pricestore->scip != NULL);

   if( num > pricestore->colssize[arrayindex] )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(pricestore->scip, num);

      GCG_SET_LOCK(&pricestore->locks->memorylock);
      retcode = SCIPreallocBlockMemoryArray(pricestore->scip, &pricestore->cols[arrayindex], pricestore->colssize[arrayindex], newsize);
      if( retcode == SCIP_OKAY )
         retcode = SCIPreallocBlockMemoryArray(pricestore->scip, &pricestore->objparallelisms[arrayindex], pricestore->colssize[arrayindex], newsize);
      if( retcode == SCIP_OKAY )
         retcode = SCIPreallocBlockMemoryArray(pricestore->scip, &pricestore->orthogonalities[arrayindex], pricestore->colssize[arrayindex], newsize);
      if( retcode == SCIP_OKAY )
         retcode = SCIPreallocBlockMemoryArray(pricestore->scip, &pricestore->scores[arrayindex], pricestore->colssize[arrayindex], newsize);
      GCG_UNSET_LOCK(&pricestore->locks->memorylock);

      SCIP_CALL(retcode);
      pricestore->colssize[arrayindex] = newsize;
   }
   assert(num <= pricestore->colssize[arrayindex]);

   return SCIP_OKAY;
}

/** returns the index of the arrays that are used to store data of a column */
static
int pricestoreGetArrayIndex(
   GCG_PRICESTORE*       pricestore,          /**< price storage */
   GCG_COL*              col                  /**< column */
   )
{
   assert(col != NULL);
   assert(col->probnr < pricestore->narrays);
   return col->probnr;
}

/** creates price storage */
SCIP_RETCODE GCGpricestoreCreate(
   SCIP*                 scip,                /**< SCIP data structure (master problem) */
   SCIP*                 origscip,            /**< SCIP data structure (original problem) */
   GCG_PRICESTORE**      pricestore,          /**< pointer to store price storage */
   SCIP_Real             efficiacyfac,          /**< factor of -redcost/norm in score function */
   SCIP_Real             objparalfac,         /**< factor of objective parallelism in score function */
   SCIP_Real             orthofac,            /**< factor of orthogonalities in score function */
   SCIP_Real             mincolorth,          /**< minimal orthogonality of columns to add
                                                  (with respect to columns added in the current round) */
   GCG_EFFICIACYCHOICE   efficiacychoice,     /**< choice to base efficiacy on */
   int                   hashtablesize        /**< size of hashtable */
   )
{
   int i;
   assert(pricestore != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, pricestore) );

   (*pricestore)->scip = scip;
   (*pricestore)->ncolstotal = 0;
   (*pricestore)->ncolsfound = 0;
   (*pricestore)->ncolsfoundround = 0;
   (*pricestore)->ncolsapplied = 0;
   (*pricestore)->narrays = GCGgetNPricingprobs(GCGgetOriginalprob(scip));
   (*pricestore)->infarkas = FALSE;
   (*pricestore)->forcecols = FALSE;
   (*pricestore)->efficiacyfac = efficiacyfac;   /* factor of efficiacies in score function */
   (*pricestore)->objparalfac = objparalfac;     /* factor of objective parallelism in score function */
   (*pricestore)->orthofac = orthofac;           /* factor of orthogonalities in score function */
   (*pricestore)->mincolorth = mincolorth;       /* minimal orthogonality of columns to add
                                                      (with respect to columns added in the current round) */
   (*pricestore)->efficiacychoice = efficiacychoice;

#ifdef _OPENMP
   (*pricestore)->locks = GCGgetLocks(origscip);
#endif

   SCIP_CALL( SCIPhashtableCreate(&(*pricestore)->hashtable, SCIPblkmem(scip),
         hashtablesize, GCGhashGetKeyCol, GCGhashKeyEqCol, GCGhashKeyValCol, (void*)pricestore) );
   
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*pricestore)->cols, (*pricestore)->narrays) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*pricestore)->objparallelisms, (*pricestore)->narrays) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*pricestore)->orthogonalities, (*pricestore)->narrays) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*pricestore)->scores, (*pricestore)->narrays) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*pricestore)->colssize, (*pricestore)->narrays) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*pricestore)->ncols, (*pricestore)->narrays) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*pricestore)->nforcedcols, (*pricestore)->narrays) );

   for( i = 0; i < (*pricestore)->narrays; ++i )
   {
      (*pricestore)->cols[i] = NULL;
      (*pricestore)->objparallelisms[i] = NULL;
      (*pricestore)->orthogonalities[i] = NULL;
      (*pricestore)->scores[i] = NULL;
      (*pricestore)->colssize[i] = 0;
      (*pricestore)->ncols[i] = 0;
      (*pricestore)->nforcedcols[i] = 0;
   }

   return SCIP_OKAY;
}

/** frees price storage */
SCIP_RETCODE GCGpricestoreFree(
   SCIP*                 scip,                /**< SCIP data structure */
   GCG_PRICESTORE**      pricestore           /**< pointer to store price storage */
   )
{
   int i;
   assert(scip == (*pricestore)->scip);
   assert(pricestore != NULL);
   assert(*pricestore != NULL);
   assert((*pricestore)->ncolstotal == 0);

   SCIPhashtableFree(&(*pricestore)->hashtable);

   for( i = 0; i < (*pricestore)->narrays; ++i )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*pricestore)->cols[i], (*pricestore)->colssize[i]);
      SCIPfreeBlockMemoryArrayNull(scip, &(*pricestore)->objparallelisms[i], (*pricestore)->colssize[i]);
      SCIPfreeBlockMemoryArrayNull(scip, &(*pricestore)->orthogonalities[i], (*pricestore)->colssize[i]);
      SCIPfreeBlockMemoryArrayNull(scip, &(*pricestore)->scores[i], (*pricestore)->colssize[i]);
   }
   SCIPfreeBlockMemoryArray(scip, &(*pricestore)->cols, (*pricestore)->narrays);
   SCIPfreeBlockMemoryArray(scip, &(*pricestore)->objparallelisms, (*pricestore)->narrays);
   SCIPfreeBlockMemoryArray(scip, &(*pricestore)->orthogonalities, (*pricestore)->narrays);
   SCIPfreeBlockMemoryArray(scip, &(*pricestore)->scores, (*pricestore)->narrays);
   SCIPfreeBlockMemoryArray(scip, &(*pricestore)->colssize, (*pricestore)->narrays);
   SCIPfreeBlockMemoryArray(scip, &(*pricestore)->ncols, (*pricestore)->narrays);
   SCIPfreeBlockMemoryArray(scip, &(*pricestore)->nforcedcols, (*pricestore)->narrays);
   SCIPfreeBlockMemory(scip, pricestore);

   return SCIP_OKAY;
}

/** informs price storage, that Farkas pricing starts now */
void GCGpricestoreStartFarkas(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   )
{
   assert(pricestore != NULL);
   assert(pricestore->ncolstotal == 0);

   pricestore->infarkas = TRUE;
}

/** informs price storage, that Farkas pricing is now finished */
void GCGpricestoreEndFarkas(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   )
{
   assert(pricestore != NULL);
   assert(pricestore->ncolstotal == 0);

   pricestore->infarkas = FALSE;
}

/** informs price storage, that the following cols should be used in any case */
void GCGpricestoreStartForceCols(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   )
{
   assert(pricestore != NULL);
   assert(!pricestore->forcecols);

   pricestore->forcecols = TRUE;
}

/** informs price storage, that the following cols should no longer be used in any case */
void GCGpricestoreEndForceCols(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   )
{
   assert(pricestore != NULL);
   assert(pricestore->forcecols);

   pricestore->forcecols = FALSE;
}

/** removes a non-forced col from the price storage */
static
void pricestoreDelCol(
   GCG_PRICESTORE*       pricestore,          /**< price storage */
   int                   arrayindex,          /**< index of the arrays */
   int                   pos,                 /**< position of col to delete */
   SCIP_Bool             freecol              /**< should col be freed */
   )
{
   assert(pricestore != NULL);
   assert(pricestore->cols != NULL);
   assert(pricestore->nforcedcols[arrayindex] <= pos && pos < pricestore->ncols[arrayindex]);
   
   SCIPhashtableRemove(pricestore->hashtable, pricestore->cols[arrayindex][pos]);
   pricestore->cols[arrayindex][pos]->pos = -1;

   /* free the column */
   if( freecol )
      GCGfreeGcgCol(&(pricestore->cols[arrayindex][pos]));

   pricestore->ncols[arrayindex]--;
   if( pos != pricestore->ncols[arrayindex] )
   {
      /* move last col to the empty position */
      pricestore->cols[arrayindex][pos] = pricestore->cols[arrayindex][pricestore->ncols[arrayindex]];
      pricestore->objparallelisms[arrayindex][pos] = pricestore->objparallelisms[arrayindex][pricestore->ncols[arrayindex]];
      pricestore->orthogonalities[arrayindex][pos] = pricestore->orthogonalities[arrayindex][pricestore->ncols[arrayindex]];
      pricestore->scores[arrayindex][pos] = pricestore->scores[arrayindex][pricestore->ncols[arrayindex]];
      pricestore->cols[arrayindex][pos]->pos = pos;
   }
}

/** for a given column, check if an identical column already exists in the price storage;
 *  if one exists, return its position, otherwise, return -1
 */
static
int pricestoreFindEqualCol(
   GCG_PRICESTORE*       pricestore,         /**< price storage */
   GCG_COL*              col                 /**< column to be checked */
   )
{
   GCG_COL* othercol = (GCG_COL*)SCIPhashtableRetrieve(pricestore->hashtable, (void*)col);

   if( othercol != NULL )
      return othercol->pos;
   else
      return -1;
}

/** adds col to price storage;
 *  if the col should be forced to enter the LP, an infinite score will be used
 */
SCIP_RETCODE GCGpricestoreAddCol(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_PRICESTORE*       pricestore,         /**< price storage */
   GCG_COL*              col,                /**< priced col */
   SCIP_Bool             forcecol,           /**< should the col be forced to enter the LP? */
   SCIP_Bool*            added
   )
{
   SCIP_Real colobjparallelism;
   SCIP_Real colscore;

   int oldpos;
   int pos;
   int arrayindex;

   assert(pricestore != NULL);
   assert(col != NULL);
   assert(col->pos == -1);
   assert(added != NULL);

   arrayindex = pricestoreGetArrayIndex(pricestore, col);
   assert(pricestore->nforcedcols[arrayindex] <= pricestore->ncols[arrayindex]);

   /* a col is forced to enter the LP if
    *  - we construct the initial LP, or
    *  - it has infinite score factor, or
    * if it is a non-forced col and no cols should be added, abort
    */
   forcecol = forcecol || pricestore->forcecols;

   GCGcolComputeNorm(scip, col);

   if( forcecol )
   {
      colscore = SCIPinfinity(scip);
      colobjparallelism = 1.0;
   }
   else
   {
      /* initialize values to invalid (will be initialized during col filtering) */
      colscore = SCIP_INVALID;

      if( SCIPisPositive(scip, pricestore->objparalfac) )
         colobjparallelism = GCGcolComputeDualObjPara(scip, col);
      else
         colobjparallelism = 0.0; /* no need to calculate it */
   }

   GCG_SET_LOCK(&pricestore->locks->pricestorelock);
   oldpos = pricestoreFindEqualCol(pricestore, col);
   GCG_UNSET_LOCK(&pricestore->locks->pricestorelock);

   pos = -1;

   /* If the column is no duplicate of an existing one, add it */
   if( oldpos == -1 )
   {
#ifndef NDEBUG
      int i;
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_SOL* sol;
      SCIP_Bool feasible;
      if( SCIPgetStage(col->pricingprob) < SCIP_STAGE_PRESOLVING )
      {
         SCIPcreateSol(col->pricingprob, &sol, NULL);
         SCIPsetSolVals(col->pricingprob, sol, col->nvars, col->vars, col->vals);
         SCIPcheckSolOrig(col->pricingprob, sol, &feasible, TRUE, TRUE);
         if( !feasible )
            SCIPprintSol(col->pricingprob, sol, NULL, FALSE);
         assert(feasible);
         SCIPfreeSol(col->pricingprob, &sol);
      }
      for( i = 0; i < col->nvars; ++i )
      {
         var = col->vars[i];
         val = col->vals[i];

         assert(SCIPisFeasGE(col->pricingprob, val, SCIPvarGetLbGlobal(GCGpricingVarGetOrigvars(var)[0])) &&
                SCIPisFeasLE(col->pricingprob, val, SCIPvarGetUbGlobal(GCGpricingVarGetOrigvars(var)[0])));
      }
#endif
      /* get enough memory to store the col */
      SCIP_CALL( pricestoreEnsureColsMem(pricestore, arrayindex, pricestore->ncols[arrayindex]+1) );
      assert(pricestore->ncols[arrayindex] < pricestore->colssize[arrayindex]);

      if( forcecol )
      {
         /* make room at the beginning of the array for forced col */
         pos = pricestore->nforcedcols[arrayindex];
         pricestore->cols[arrayindex][pricestore->ncols[arrayindex]] = pricestore->cols[arrayindex][pos];
         pricestore->objparallelisms[arrayindex][pricestore->ncols[arrayindex]] = pricestore->objparallelisms[arrayindex][pos];
         pricestore->orthogonalities[arrayindex][pricestore->ncols[arrayindex]] = pricestore->orthogonalities[arrayindex][pos];
         pricestore->scores[arrayindex][pricestore->ncols[arrayindex]] = pricestore->scores[arrayindex][pos];
         pricestore->cols[arrayindex][pricestore->ncols[arrayindex]]->pos = pricestore->ncols[arrayindex];
         pricestore->nforcedcols[arrayindex]++;
      }
      else
         pos = pricestore->ncols[arrayindex];

      pricestore->ncols[arrayindex]++;
      #pragma omp atomic update
      pricestore->ncolstotal++;

      /* update statistics of total number of found cols */
      #pragma omp atomic update
      pricestore->ncolsfound++;
      #pragma omp atomic update
      pricestore->ncolsfoundround++;
   }
   /* Otherwise, if the new column is forced and the duplicate one is not,
    * remove the duplicate and replace it by the new column
    */
   else if( forcecol && oldpos >= pricestore->nforcedcols[arrayindex] )
   {
      assert(GCGcolIsEq(pricestore->cols[arrayindex][oldpos], col));
      GCGfreeGcgCol(&pricestore->cols[arrayindex][oldpos]);
      pricestore->cols[arrayindex][oldpos] = pricestore->cols[arrayindex][pricestore->nforcedcols[arrayindex]];
      pricestore->objparallelisms[arrayindex][oldpos] = pricestore->objparallelisms[arrayindex][pricestore->nforcedcols[arrayindex]];
      pricestore->orthogonalities[arrayindex][oldpos] = pricestore->orthogonalities[arrayindex][pricestore->nforcedcols[arrayindex]];
      pricestore->scores[arrayindex][oldpos] = pricestore->scores[arrayindex][pricestore->nforcedcols[arrayindex]];

      pos = pricestore->nforcedcols[arrayindex];
      pricestore->nforcedcols[arrayindex]++;
   }
   /* The column already exists and is not forced, free it */
   else
   {
      assert(GCGcolIsEq(pricestore->cols[arrayindex][oldpos], col));
      /* @todo: This is a little dangerous */
      assert(col->pos == -1);
      GCGfreeGcgCol(&col);
   }

   if( pos > -1 )
   {
      SCIPdebugMessage("adding col %p to price storage of size %d (forcecol=%u)\n",
         (void*)col, pricestore->ncolstotal, forcecol);

      /* add col to arrays */
      pricestore->cols[arrayindex][pos] = col;
      pricestore->objparallelisms[arrayindex][pos] = colobjparallelism;
      pricestore->orthogonalities[arrayindex][pos] = 1.0;
      pricestore->scores[arrayindex][pos] = colscore;
      assert(col->pos == -1);
      col->pos = pos;

      GCG_SET_LOCK(&pricestore->locks->pricestorelock);
      GCG_SET_LOCK(&pricestore->locks->memorylock);
      SCIPhashtableInsert(pricestore->hashtable, col);
      GCG_UNSET_LOCK(&pricestore->locks->memorylock);
      GCG_UNSET_LOCK(&pricestore->locks->pricestorelock);

      *added = TRUE;
   }
   else
   {
      *added = FALSE;
   }

   return SCIP_OKAY;
}

/** updates the orthogonalities and scores of the non-forced cols after the given col was added to the LP */
static
SCIP_RETCODE pricestoreUpdateOrthogonalities(
   GCG_PRICESTORE*       pricestore,          /**< price storage */
   GCG_COL*              col,                /**< col that was applied */
   SCIP_Real             mincolorthogonality /**< minimal orthogonality of cols to apply to LP */
   )
{
   int pos;
   int i;

   assert(pricestore != NULL);

   for( i = 0; i < pricestore->narrays; ++i)
   {
      pos = pricestore->nforcedcols[i];
      while( pos < pricestore->ncols[i] )
      {
         SCIP_Real thisortho;

         /* update orthogonality */
         thisortho = GCGcolComputeOrth(pricestore->scip, col, pricestore->cols[i][pos]);

         if( thisortho < pricestore->orthogonalities[i][pos] )
         {
            if( thisortho < mincolorthogonality )
            {
               /* col is too parallel: delete the col */
               SCIPdebugMessage("    -> deleting parallel col %p after adding %p (pos=%d, orthogonality=%g, score=%g)\n",
                  (void*) pricestore->cols[i][pos], (void*) col, pos, thisortho, pricestore->scores[i][pos]);
               pricestoreDelCol(pricestore, i, pos, TRUE);
               continue;
            }
            else
            {
               SCIP_Real colefficiacy;

               /* calculate col's efficacy */
               switch ( pricestore->efficiacychoice )
               {
               case GCG_EFFICIACYCHOICE_DANTZIG:
                  colefficiacy = -1.0 * GCGcolGetRedcost(pricestore->cols[i][pos]);
                  break;
               case GCG_EFFICIACYCHOICE_STEEPESTEDGE:
                  colefficiacy = -1.0 * GCGcolGetRedcost(pricestore->cols[i][pos])/ GCGcolGetNorm(col);
                  break;
               case GCG_EFFICIACYCHOICE_LAMBDA:
                  SCIPerrorMessage("Lambda pricing not yet implemented.\n");
                  return SCIP_INVALIDCALL;
               default:
                  SCIPerrorMessage("Invalid efficiacy choice.\n");
                  return SCIP_INVALIDCALL;
               }

               /* recalculate score */
               pricestore->orthogonalities[i][pos] = thisortho;
               assert( pricestore->objparallelisms[i][pos] != SCIP_INVALID ); /*lint !e777*/
               assert( pricestore->scores[i][pos] != SCIP_INVALID ); /*lint !e777*/


               pricestore->scores[i][pos] = pricestore->efficiacyfac * colefficiacy
                  + pricestore->objparalfac * pricestore->objparallelisms[i][pos]
                  + pricestore->orthofac * thisortho;
            }
         }
         pos++;
      }
   }
   return SCIP_OKAY;
}

/** adds the given col to priced vars, updates the orthogonalities and scores of remaining cols, and frees non-forced cols */
static
SCIP_RETCODE pricestoreApplyCol(
   GCG_PRICESTORE*       pricestore,         /**< price storage */
   GCG_COL*              col,                /**< col to apply to the LP */
   SCIP_Bool             force,              /**< force column */
   SCIP_Real             mincolorthogonality,/**< minimal orthogonality of cols to apply to LP */
   SCIP_Real             score,              /**< score of column (or -1.0 if not specified) */
   SCIP_Bool*            added               /**< pointer to store whether the column was added */
   )
{
   int arrayindex;
   assert(pricestore != NULL);
   assert(added != NULL);

   arrayindex = pricestoreGetArrayIndex(pricestore, col);

   SCIP_CALL( GCGcreateNewMasterVarFromGcgCol(pricestore->scip, pricestore->infarkas, col, force, added, NULL, score) );
   assert(*added);

   if( !force )
   {
      assert(pricestore->cols[arrayindex][col->pos] == col);
      pricestoreDelCol(pricestore, arrayindex, col->pos, FALSE);
   }

   /* update the orthogonalities if needed */
   if( SCIPisGT(pricestore->scip, mincolorthogonality, SCIPepsilon(pricestore->scip)) || SCIPisPositive(pricestore->scip, pricestore->orthofac))
      SCIP_CALL( pricestoreUpdateOrthogonalities(pricestore, col, mincolorthogonality) );
   
   if( !force )
   {
      GCGfreeGcgCol(&col);
   }

   return SCIP_OKAY;
}

/** returns the position of the best non-forced col in the cols array */
static
GCG_COL* pricestoreGetBestCol(
   GCG_PRICESTORE*       pricestore,          /**< price storage */
   int*                  ncolsapplied,        /**< number of cols applied per array index */
   int                   maxpricecolsprob     /**< maximum number of cols allowed per array index */
   )
{
   SCIP_Real bestscore;
   GCG_COL* bestcol;
   int pos;
   int i;

   assert(pricestore != NULL);

   bestscore = SCIP_REAL_MIN;
   bestcol = NULL;
   for( i = 0; i < pricestore->narrays; ++i )
   {
      if( ncolsapplied[i] < maxpricecolsprob )
      {
         for( pos = pricestore->nforcedcols[i]; pos < pricestore->ncols[i]; pos++ )
         {
            /* check if col is current best col */
            assert( pricestore->scores[i][pos] != SCIP_INVALID ); /*lint !e777*/
            if( pricestore->scores[i][pos] > bestscore )
            {
               bestscore = pricestore->scores[i][pos];
               bestcol = pricestore->cols[i][pos];
            }
         }
      }
   }

   return bestcol;
}

/** computes score for dual solution and initialized orthogonalities */
static
SCIP_RETCODE computeScore(
   GCG_PRICESTORE*       pricestore,          /**< price storage */
   int                   arrayindex,          /**< index of the arrays */
   int                   pos                  /**< position of col to handle */
   )
{
   GCG_COL* col;
   SCIP_Real colefficiacy;
   SCIP_Real colscore;

   col = pricestore->cols[arrayindex][pos];

   /* calculate cut's efficacy */
   switch ( pricestore->efficiacychoice )
   {
   case GCG_EFFICIACYCHOICE_DANTZIG:
      colefficiacy = -1.0 * GCGcolGetRedcost(pricestore->cols[arrayindex][pos]);
      break;
   case GCG_EFFICIACYCHOICE_STEEPESTEDGE:
      colefficiacy = -1.0 * GCGcolGetRedcost(pricestore->cols[arrayindex][pos])/ GCGcolGetNorm(col);
      break;
   case GCG_EFFICIACYCHOICE_LAMBDA:
      SCIPerrorMessage("Lambda pricing not yet implemented.\n");
      return SCIP_INVALIDCALL;
   default:
      SCIPerrorMessage("Invalid efficiacy choice.\n");
      return SCIP_INVALIDCALL;
   }

   assert( pricestore->objparallelisms[arrayindex][pos] != SCIP_INVALID ); /*lint !e777*/
   colscore = pricestore->efficiacyfac * colefficiacy
            + pricestore->objparalfac * pricestore->objparallelisms[arrayindex][pos]
            + pricestore->orthofac * 1.0;
   assert( !SCIPisInfinity(pricestore->scip, colscore) );

   pricestore->scores[arrayindex][pos] = colscore;

   /* make sure that the orthogonalities are initialized to 1.0 */
   pricestore->orthogonalities[arrayindex][pos] = 1.0;

   return SCIP_OKAY;
}

/** adds cols to priced vars and clears price storage */
SCIP_RETCODE GCGpricestoreApplyCols(
   GCG_PRICESTORE*       pricestore,         /**< price storage */
   GCG_COLPOOL*          colpool,            /**< GCG column pool */
   SCIP_Bool             usecolpool,         /**< use column pool? */
   int*                  nfoundvars          /**< pointer to store number of variables that were added to the problem */
   )
{
   SCIP* scip;
   SCIP_Bool added;
   int* ncolsappliedprob;
   SCIP_Real mincolorthogonality;
   int maxpricecols;
   int maxpricecolsprob;
   int ncolsapplied;
   int pos;
   int i;

   assert(pricestore != NULL);

   scip = pricestore->scip;

   SCIPdebugMessage("applying %d cols\n", pricestore->ncolstotal);

   /* get maximal number of cols to add to the LP */
   maxpricecols = GCGpricerGetMaxColsRound(scip);
   maxpricecolsprob = GCGpricerGetMaxColsProb(scip);

   ncolsapplied = 0;
   SCIP_CALL( SCIPallocClearBufferArray(scip, &ncolsappliedprob, GCGgetNPricingprobs(GCGmasterGetOrigprob(scip))) );

   /* set minimal col orthogonality */
   mincolorthogonality = pricestore->mincolorth;
   mincolorthogonality = MAX(mincolorthogonality, SCIPepsilon(scip)); /*lint !e666 */

   /* Compute scores for all non-forced cols and initialize orthogonalities - make sure all cols are initialized again for the current dual solution */
   for( i = 0; i < pricestore->narrays; ++i )
   {
      for( pos = pricestore->nforcedcols[i]; pos < pricestore->ncols[i]; pos++ )
      {
         SCIP_CALL( computeScore(pricestore, i, pos) );
      }
   }

   /* apply all forced cols */
   for( i = 0; i < pricestore->narrays; ++i )
   {
      for( pos = 0; pos < pricestore->nforcedcols[i]; pos++ )
      {
         GCG_COL* col;
         int probnr;

         col = pricestore->cols[i][pos];
         assert(SCIPisInfinity(scip, pricestore->scores[i][pos]));

         probnr = GCGcolGetProbNr(col);

         /* add col to the priced vars and update orthogonalities */
         SCIPdebugMessage(" -> applying forced col %p (probnr = %d)\n", (void*) col, probnr);

         SCIP_CALL( pricestoreApplyCol(pricestore, col, TRUE, mincolorthogonality, pricestore->scores[i][pos], &added) );

         if( added )
         {
            ++ncolsapplied;
            ++ncolsappliedprob[probnr];
         }
      }
   }

   /* apply non-forced cols */
   while( TRUE )
   {
      GCG_COL* col;
      int bestpos;
      SCIP_Real score;
      int probnr;

      /* get best non-forced col */
      col = NULL;
      if( ncolsapplied < maxpricecols )
      {
         col = pricestoreGetBestCol(pricestore, ncolsappliedprob, maxpricecolsprob);
      }
      else
      {
         for( i = 0; i < pricestore->narrays; ++i )
         {
            if( pricestore->nforcedcols[i] < pricestore->ncols[i] )
            {
               col = pricestore->cols[i][pricestore->nforcedcols[i]];
               assert(col != NULL);
            }
         }
      }
      if( col == NULL )
         break;
      i = pricestoreGetArrayIndex(pricestore, col);
      bestpos = col->pos;
      assert(pricestore->nforcedcols[i] <= bestpos && bestpos < pricestore->ncols[i]);
      assert(pricestore->scores[i][bestpos] != SCIP_INVALID ); /*lint !e777*/
      score = pricestore->scores[i][bestpos];
      assert(!SCIPisInfinity(scip, pricestore->scores[i][bestpos]));
      probnr = GCGcolGetProbNr(col);

      /* Do not add (non-forced) non-violated cols.
       * Note: do not take SCIPsetIsEfficacious(), because constraint handlers often add cols w.r.t. SCIPsetIsFeasPositive().
       * Note2: if pricerating/feastolfac != -1, constraint handlers may even add cols w.r.t. SCIPsetIsPositive(); those are currently rejected here
       */
      if( SCIPisDualfeasNegative(scip, GCGcolGetRedcost(col)) && ncolsapplied < maxpricecols )
      {
         assert(ncolsappliedprob[probnr] < maxpricecolsprob);
         /* add col to the LP and update orthogonalities */
         SCIP_CALL( pricestoreApplyCol(pricestore, col, FALSE, mincolorthogonality, score, &added) );
         if( added )
         {
            SCIPdebugMessage(" -> applying col %p (pos=%d/%d, probnr=%d, efficacy=%g, objparallelism=%g, orthogonality=%g, score=%g)\n",
                  (void*) col, bestpos+1, pricestore->ncols[i], probnr, GCGcolGetRedcost(col), pricestore->objparallelisms[i][bestpos],
                  pricestore->orthogonalities[i][bestpos], pricestore->scores[i][bestpos]);

            ++ncolsapplied;
            ++ncolsappliedprob[probnr];
         }
      }
      else if( usecolpool )
      {
         pricestoreDelCol(pricestore, i, bestpos, FALSE);
         SCIP_CALL( GCGcolpoolAddCol(colpool, col, TRUE) );
      }
      else
      {
         pricestoreDelCol(pricestore, i, bestpos, TRUE);
      }
   }

   *nfoundvars = ncolsapplied;

   /* clear the price storage and reset statistics for price round */
   GCGpricestoreClearCols(pricestore);

   SCIPfreeBufferArray(scip, &ncolsappliedprob);

   return SCIP_OKAY;
}

/** clears the price storage without adding the cols to the LP */
void GCGpricestoreClearCols(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   )
{
   int c;
   int i;

   assert(pricestore != NULL);

   SCIPdebugMessage("clearing %d cols\n", pricestore->ncolstotal);

   /* release cols */
   for( i = 0; i < pricestore->narrays; ++i )
   {
      for( c = 0; c < pricestore->ncols[i]; ++c )
      {
         GCGfreeGcgCol(&(pricestore->cols[i][c]));
      }
      /* reset counters */
      pricestore->ncols[i] = 0;
      pricestore->nforcedcols[i] = 0;
   }

   SCIPhashtableRemoveAll(pricestore->hashtable);

   /* reset counters */
   pricestore->ncolstotal = 0;
   pricestore->ncolsfoundround = 0;

   /* if we have just finished the initial LP construction, free the (potentially large) cols array */
   if( pricestore->infarkas )
   {
      for( i = 0; i < pricestore->narrays; ++i )
      {
         SCIPfreeBlockMemoryArrayNull(pricestore->scip, &pricestore->cols[i], pricestore->colssize[i]);
         SCIPfreeBlockMemoryArrayNull(pricestore->scip, &pricestore->objparallelisms[i], pricestore->colssize[i]);
         SCIPfreeBlockMemoryArrayNull(pricestore->scip, &pricestore->orthogonalities[i], pricestore->colssize[i]);
         SCIPfreeBlockMemoryArrayNull(pricestore->scip, &pricestore->scores[i], pricestore->colssize[i]);

         pricestore->colssize[i] = 0;
      }
   }
}

/** get cols in the price storage */
GCG_COL** GCGpricestoreGetCols(
   GCG_PRICESTORE*       pricestore,          /**< price storage */
   int                   arrayindex           /**< index of the arrays */
   )
{
   assert(pricestore != NULL);

   return pricestore->cols[arrayindex];
}

/** get number of cols in the price storage */
int GCGpricestoreGetNCols(
   GCG_PRICESTORE*       pricestore,          /**< price storage */
   int                   arrayindex           /**< index of the arrays */
   )
{
   assert(pricestore != NULL);

   return pricestore->ncols[arrayindex];
}

/** get number of cols in the price storage */
int GCGpricestoreGetNColsTotal(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->ncolstotal;
}

/** get total number of cols found so far */
int GCGpricestoreGetNColsFound(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->ncolsfound;
}

/** get number of cols found so far in current price round */
int GCGpricestoreGetNColsFoundRound(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->ncolsfoundround;
}

/** get total number of cols applied to the LPs */
int GCGpricestoreGetNColsApplied(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->ncolsapplied;
}
