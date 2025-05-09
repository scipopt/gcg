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

/**@file   colpool.c
 * @brief  methods for storing cols in a col pool (based on SCIP's cut pool)
 * @author Jonas Witt
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/clock.h"

#include "gcg/gcg.h"
#include "gcg/pub_gcgcol.h"
#include "gcg/colpool.h"
#include "gcg/struct_colpool.h"
#include "gcg/pricestore_gcg.h"
#include "gcg/struct_pricestore_gcg.h"
#include "gcg/pricer_gcg.h"
#include "gcg/struct_gcgcol.h"

#define GCG_USESMALLTABLES FALSE
#define GCG_HASHSIZE_COLPOOLS_SMALL 100 /**< size of hash table in col pools for small problems */
#define GCG_HASHSIZE_COLPOOLS       500 /**< size of hash table in col pools */

/*
 * dynamic memory arrays
 */

/** resizes cols array to be able to store at least num entries */
static
SCIP_RETCODE colpoolEnsureColsMem(
   GCG_COLPOOL*          colpool,            /**< col pool */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(colpool != NULL);

   if( num > colpool->colssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(colpool->scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(colpool->scip, &colpool->cols, colpool->colssize, newsize) );
      colpool->colssize = newsize;
   }
   assert(num <= colpool->colssize);

   return SCIP_OKAY;
}

/*
 * Col methods
 */

/*
 * Colpool methods
 */

/** creates col pool */
SCIP_RETCODE GCGcolpoolCreate(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_COLPOOL**         colpool,            /**< pointer to store col pool */
   int                   agelimit            /**< maximum age a col can reach before it is deleted from the pool (-1 fpr no limit) */
   )
{
   SCIP* scip = GCGgetMasterprob(gcg);
   assert(colpool != NULL);
   assert(agelimit >= -1);

   SCIP_CALL( SCIPallocMemory(scip, colpool) );

   SCIP_CALL( SCIPcreateClock(scip, &(*colpool)->poolclock) );

   SCIP_CALL( SCIPhashtableCreate(&(*colpool)->hashtable, SCIPblkmem(scip),
         (GCG_USESMALLTABLES ? GCG_HASHSIZE_COLPOOLS_SMALL :  GCG_HASHSIZE_COLPOOLS),
         GCGhashGetKeyCol, GCGhashKeyEqCol, GCGhashKeyValCol, (void*) scip) );

   (*colpool)->gcg = gcg;
   (*colpool)->scip = scip;
   (*colpool)->nodenr = -1;
   (*colpool)->infarkas = FALSE;
   (*colpool)->cols = NULL;
   (*colpool)->colssize = 0;
   (*colpool)->ncols = 0;
   (*colpool)->agelimit = agelimit;
   (*colpool)->processedlp = -1;
   (*colpool)->processedlpsol = -1;
   (*colpool)->firstunprocessed = 0;
   (*colpool)->firstunprocessedsol = 0;
   (*colpool)->maxncols = 0;
   (*colpool)->ncalls = 0;
   (*colpool)->ncolsfound = 0;

   return SCIP_OKAY;
}

/** frees col pool */
SCIP_RETCODE GCGcolpoolFree(
   GCG_COLPOOL**        colpool             /**< pointer to store col pool */
   )
{
   SCIP* scip;
   assert(colpool != NULL);
   assert(*colpool != NULL);

   scip = (*colpool)->scip;

   /* remove all cols from the pool */
   SCIP_CALL( GCGcolpoolClear(*colpool) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Pricing time in colpool = %f sec\n", GCGcolpoolGetTime(*colpool));

   /* free clock */
   SCIP_CALL( SCIPfreeClock(scip, &(*colpool)->poolclock) );

   /* free hash table */
   SCIPhashtableFree(&(*colpool)->hashtable);

   SCIPfreeBlockMemoryArrayNull(scip, &(*colpool)->cols, (*colpool)->colssize);
   SCIPfreeMemory(scip, colpool);

   return SCIP_OKAY;
}

/** removes the col from the col pool */
static
SCIP_RETCODE colpoolDelCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col,                /**< col to remove */
   SCIP_Bool             freecol             /**< should the col be freed? */
   )
{
   int pos;

   assert(colpool != NULL);
   assert(colpool->firstunprocessed <= colpool->ncols);
   assert(colpool->firstunprocessedsol <= colpool->ncols);
   assert(col != NULL);

   pos = col->pos;
   col->pos = -1;
   assert(0 <= pos && pos < colpool->ncols);
   assert(colpool->cols[pos] == col);

   /* remove the col from the hash table */
   assert(SCIPhashtableExists(colpool->hashtable, (void*)col));
   SCIP_CALL( SCIPhashtableRemove(colpool->hashtable, (void*)col) );

   /* free the col */
   if( freecol )
      SCIP_CALL( GCGfreeGcgCol(&colpool->cols[pos]) );

   /* move the last col of the pool to the free position */
   if( pos < colpool->ncols-1 )
   {
      colpool->cols[pos] = colpool->cols[colpool->ncols-1];
      colpool->cols[pos]->pos = pos;
   }

   colpool->ncols--;

   return SCIP_OKAY;
}


/** removes all rows from the col pool */
SCIP_RETCODE GCGcolpoolClear(
   GCG_COLPOOL*          colpool             /**< col pool */
   )
{
   int i;

   assert(colpool != NULL);

   /* free cols (in reverse order!) */
   for( i = colpool->ncols - 1; i >= 0; --i )
   {
      SCIP_CALL( colpoolDelCol(colpool, colpool->cols[i], TRUE) );
   }
   colpool->ncols = 0;

   return SCIP_OKAY;
}

/** if not already existing, adds col to col pool and captures it */
SCIP_RETCODE GCGcolpoolAddCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col,                /**< column to add */
   SCIP_Bool             freeduplicate       /**< shouldl the col be freed if it is a duplicate? */
   )
{
   assert(colpool != NULL);
   assert(col != NULL);

   /* check in hash table, if col already exists in the pool */
   if( SCIPhashtableRetrieve(colpool->hashtable, (void*)col) == NULL )
   {
      SCIP_CALL( GCGcolpoolAddNewCol(colpool, col) );
   }
   else if( freeduplicate )
   {
      assert(col->pos == -1);
      SCIP_CALL( GCGfreeGcgCol(&col) );
   }

   return SCIP_OKAY;
}

/** adds row to col pool and captures it; doesn't check for multiple cols */
SCIP_RETCODE GCGcolpoolAddNewCol(
   GCG_COLPOOL*         colpool,            /**< col pool */
   GCG_COL*             col                 /**< column to add */
   )
{

   assert(colpool != NULL);
   assert(col != NULL);
   assert(col->pos == -1);

   col->pos = colpool->ncols;

   /* add col to the pool */
   SCIP_CALL( colpoolEnsureColsMem(colpool, colpool->ncols+1) );
   colpool->cols[colpool->ncols] = col;
   colpool->ncols++;
   colpool->maxncols = MAX(colpool->maxncols, colpool->ncols);

   /* insert col in the hash table */
   SCIP_CALL( SCIPhashtableInsert(colpool->hashtable, (void*)col) );

   return SCIP_OKAY;
}

/** removes the LP row from the col pool */
SCIP_RETCODE GCGcolpoolDelCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col,                /**< col to remove */
   SCIP_Bool             freecol             /**< should the col be freed? */
   )
{
   assert(colpool != NULL);
   assert(col != NULL);

   /* find the col in hash table */
   col = (GCG_COL*)SCIPhashtableRetrieve(colpool->hashtable, (void*)col);
   if( col == NULL )
   {
      SCIPerrorMessage("col %p is not existing in colpool %p\n", (void*)col, (void*)colpool);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( colpoolDelCol(colpool, col, freecol) );

   return SCIP_OKAY;
}


/** prices cols of the col pool */
SCIP_RETCODE GCGcolpoolPrice(
   GCG_COLPOOL*          colpool,            /**< col pool */
   SCIP_SOL*             sol,                /**< solution to be separated (or NULL for LP-solution) */
   int*                  nfoundvars          /**< pointer to store the result of the separation call */
   )
{
   GCG_COL* col;
   int firstunproc;
   int c;

   assert(colpool != NULL);
   assert(colpool->firstunprocessed <= colpool->ncols);
   assert(colpool->firstunprocessedsol <= colpool->ncols);
   assert(nfoundvars != NULL);
   assert(SCIPnodeGetType(SCIPgetCurrentNode(colpool->scip)) != SCIP_NODETYPE_PROBINGNODE);

   colpool->ncalls++;

   SCIPdebugMessage("separating%s col pool %p with %d cols, beginning with col %d\n", ( sol == NULL ) ? "" : " solution from", (void*)colpool, colpool->ncols, firstunproc);

   /* start timing */
   SCIP_CALL( SCIPstartClock(colpool->scip, colpool->poolclock) );

   /* process all unprocessed cols in the pool */
   *nfoundvars = 0;

   for( c = colpool->ncols - 1; c >= 0; --c )
   {
      SCIP_Real redcost;

      col = colpool->cols[c];
      assert(col != NULL);
      assert(col->pos == c);

      redcost = GCGcolGetRedcost(col);

      if( SCIPisDualfeasNegative(colpool->scip, redcost) )
      {
         SCIP_Bool added;
         /* insert col in separation storage */
         SCIPdebugMessage(" -> col %p from the col pool (redcost: %g)\n",
            (void*)col, redcost );

         SCIP_CALL( colpoolDelCol(colpool, col, FALSE) );
         SCIP_CALL( GCGpricerAddColResult(colpool->gcg, col, &added) );

         if( added )
            *nfoundvars = *nfoundvars + 1;

         col->age = 0;
      }
      else
      {
         col->age++;
         if( GCGcolIsAged(col, colpool->agelimit) )
         {
            SCIP_CALL( colpoolDelCol(colpool, col, TRUE) );
         }
      }
   }

   /* update the number of found cols */
   colpool->ncolsfound += *nfoundvars; /*lint !e776*/

   /* stop timing */
   SCIP_CALL( SCIPstopClock(colpool->scip, colpool->poolclock) );

   return SCIP_OKAY;
}

/** update node at which columns of column pool are feasible */
SCIP_RETCODE GCGcolpoolUpdateNode(
   GCG_COLPOOL*         colpool             /**< col pool */
   )
{
   assert(colpool != NULL);
   assert(SCIPnodeGetType(SCIPgetCurrentNode(colpool->scip)) != SCIP_NODETYPE_PROBINGNODE);

   if( colpool->nodenr < 0 )
   {
      colpool->nodenr = SCIPnodeGetNumber(SCIPgetCurrentNode(colpool->scip));
   }
   else if( colpool->nodenr != SCIPnodeGetNumber(SCIPgetCurrentNode(colpool->scip)) )
   {
      SCIP_CALL( GCGcolpoolClear(colpool) );

      colpool->nodenr = SCIPnodeGetNumber(SCIPgetCurrentNode(colpool->scip));
   }

   return SCIP_OKAY;
}

/** update reduced cost and compute master coefs of columns in column pool */
SCIP_RETCODE GCGcolpoolUpdateRedcost(
   GCG_COLPOOL*         colpool             /**< col pool */
   )
{
   GCG_COL** cols;
   int ncols;

   int i;

   ncols = GCGcolpoolGetNCols(colpool);
   cols = GCGcolpoolGetCols(colpool);

   for( i = 0; i < ncols; ++i )
   {
      GCG_COL* col;
      SCIP_Real redcost;

      col = cols[i];

      SCIP_CALL( GCGcomputeColMastercoefs(colpool->gcg, col) );

      redcost = GCGcomputeRedCostGcgCol(colpool->gcg, colpool->infarkas, col, NULL);

      GCGcolUpdateRedcost(col, redcost, FALSE);
   }

   return SCIP_OKAY;
}

/** gets number of cols in the col pool */
void GCGcolpoolStartFarkas(
   GCG_COLPOOL*         colpool             /**< col pool */
   )
{
   colpool->infarkas = TRUE;
}

/** gets number of cols in the col pool */
void GCGcolpoolEndFarkas(
   GCG_COLPOOL*         colpool             /**< col pool */
   )
{
   colpool->infarkas = FALSE;
}


/** gets array of cols in the col pool */
GCG_COL** GCGcolpoolGetCols(
   GCG_COLPOOL*         colpool             /**< col pool */
   )
{
   assert(colpool != NULL);

   return colpool->cols;
}

/** gets number of cols in the col pool */
int GCGcolpoolGetNCols(
   GCG_COLPOOL*         colpool             /**< col pool */
   )
{
   assert(colpool != NULL);

   return colpool->ncols;
}

/** gets maximum number of cols that were stored in the col pool at the same time */
int GCGcolpoolGetMaxNCols(
   GCG_COLPOOL*         colpool             /**< col pool */
   )
{
   assert(colpool != NULL);

   return colpool->maxncols;
}

/** gets time in seconds used for separating cols from the pool */
SCIP_Real GCGcolpoolGetTime(
   GCG_COLPOOL*         colpool             /**< col pool */
   )
{
   assert(colpool != NULL);

   return SCIPgetClockTime(colpool->scip, colpool->poolclock);
}

/** get number of times, the col pool was separated */
SCIP_Longint GCGcolpoolGetNCalls(
   GCG_COLPOOL*         colpool             /**< col pool */
   )
{
   assert(colpool != NULL);

   return colpool->ncalls;
}

/** get total number of cols that were separated from the col pool */
SCIP_Longint GCGcolpoolGetNColsFound(
   GCG_COLPOOL*         colpool             /**< col pool */
   )
{
   assert(colpool != NULL);

   return colpool->ncolsfound;
}

SCIP_RETCODE GCGcolpoolPropagateGlobalBounds(
   GCG_COLPOOL*          colpool
   )
{
   GCG_COL* col;
   int c;
   int i;

   for( c = colpool->ncols - 1; c >= 0; --c )
   {
      col = colpool->cols[c];
      assert(col != NULL);
      assert(col->pricingprob != NULL);

      for( i = 0; i < col->nvars; ++i )
      {
         assert(GCGvarIsPricing(col->vars[i]) && GCGpricingVarGetNOrigvars(col->vars[i]) > 0 && GCGpricingVarGetOrigvars(col->vars[i])[0] != NULL);
         if( SCIPisFeasLT(col->pricingprob, col->vals[i], SCIPvarGetLbGlobal(GCGpricingVarGetOrigvars(col->vars[i])[0])) ||
             SCIPisFeasGT(col->pricingprob, col->vals[i], SCIPvarGetUbGlobal(GCGpricingVarGetOrigvars(col->vars[i])[0])) )
         {
            colpoolDelCol(colpool, col, TRUE);
            break;
         }
      }
   }

   return SCIP_OKAY;
}
