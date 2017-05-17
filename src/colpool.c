/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   colpool.c
 * @brief  methods for storing cols in a col pool (based on SCIP's cut pool)
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/lp.h"
#include "scip/cons.h"
#include "scip/sepa.h"
#include "scip/sepastore.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"

#include "pub_gcgcol.h"
#include "colpool.h"
#include "struct_colpool.h"
#include "pricestore_gcg.h"
#include "struct_pricestore_gcg.h"

#define SCIP_HASHTABLE_USESMALL FALSE /**< size of hash table in col pools for small problems */
#define SCIP_HASHSIZE_COLPOOLS_SMALL 100 /**< size of hash table in col pools for small problems */
#define SCIP_HASHSIZE_COLPOOLS       500 /**< size of hash table in col pools */

/*
 * Hash functions
 */

/** gets the hash key of a col */
static
SCIP_DECL_HASHGETKEY(hashGetKeyCol)
{  /*lint --e{715}*/
   GCG_COL* col;

   col = (GCG_COL*)elem;
   assert(col != NULL);

   /* the key of a col is the col itself */
   return col;
}

/** returns TRUE iff both cols are identical */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqCol)
{  /*lint --e{715}*/
   /* Warning: The comparison of real values is made against default epsilon.
    *          This is ugly, but we have no settings at hand.
    */
   SCIP* scip;
   GCG_COL* col1;
   GCG_COL* col2;
   int i;

   scip = (SCIP*) userptr;
   col1 = (GCG_COL*)key1;
   col2 = (GCG_COL*)key2;
   assert(col1 != NULL);
   assert(col2 != NULL);

   assert(col1->vars != NULL || col1->nvars == 0);
   assert(col2->vars != NULL || col2->nvars == 0);

   /* compare the trivial characteristics of the cols */
   if( col1->probnr != col2->probnr
      || col1->isray != col2->isray
      || col1->nvars != col2->nvars
      //|| !SCIPisEQ(scip, col1->redcost, col2->redcost)
       )
      return FALSE;

   /* compare variables and coresponding values in sorted arrays */
   for( i = 0; i < col1->nvars; ++i )
   {
      if( col1->vars[i] != col2->vars[i]
         || !SCIPisEQ(scip, col1->vals[i], col2->vals[i]))
         return FALSE;
   }

   return TRUE;
}

static
SCIP_DECL_HASHKEYVAL(hashKeyValCol)
{  /*lint --e{715}*/
   GCG_COL* col;
   unsigned int keyval;
   SCIP_Real maxval;
   SCIP_Real minval;
   SCIP* set;

   set = (SCIP*) userptr;
   col = (GCG_COL*)key;
   assert(col != NULL);

   keyval = SCIPhashTwo(SCIPpositiveRealHashCode(0.0, 8),
                        SCIPcombineThreeInt(col->probnr, col->nvars, col->isray));
   return keyval;
}



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
      SCIP_CALL( SCIPreallocMemoryArray(colpool->scip, &colpool->cols, newsize) );
      colpool->colssize = newsize;
   }
   assert(num <= colpool->colssize);

   return SCIP_OKAY;
}



/*
 * Col methods
 */

/* TODO: is this needed? */
///** returns the ratio of LPs where the row belonging to this col was active in an LP solution, i.e.
// *  where the age of its row has not been increased
// *
// *  @see SCIPcolGetAge() to get the age of a col
// */
//SCIP_Real SCIPcolGetLPActivityQuot(
//   GCG_COL*             col                 /**< col */
//   )
//{
//   SCIP_Longint nlpsaftercreation;
//   SCIP_Longint activeinlpcounter;
//
//   assert(col != NULL);
//   assert(col->row != NULL);
//
//   nlpsaftercreation = SCIProwGetNLPsAfterCreation(col->row);
//   activeinlpcounter = SCIProwGetActiveLPCount(col->row);
//
//   return (nlpsaftercreation > 0 ? activeinlpcounter / (SCIP_Real)nlpsaftercreation : 0.0);
//}

/*
 * Colpool methods
 */

/** creates col pool */
SCIP_RETCODE GCGcolpoolCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_COLPOOL**         colpool,            /**< pointer to store col pool */
   int                   agelimit,           /**< maximum age a col can reach before it is deleted from the pool */
   SCIP_Bool             globalcolpool       /**< is this the global col pool of SCIP? */
   )
{
   assert(colpool != NULL);
   assert(agelimit >= -1);

   SCIP_CALL( SCIPallocMemory(scip, colpool) );

   SCIP_CALL( SCIPcreateClock(scip, &(*colpool)->poolclock) );

   SCIP_CALL( SCIPhashtableCreate(&(*colpool)->hashtable, SCIPblkmem(scip),
         (SCIP_HASHTABLE_USESMALL ? SCIP_HASHSIZE_COLPOOLS_SMALL : SCIP_HASHSIZE_COLPOOLS),
         hashGetKeyCol, hashKeyEqCol, hashKeyValCol, (void*) scip) );

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
   (*colpool)->globalcolpool = globalcolpool;

   return SCIP_OKAY;
}

/** frees col pool */
SCIP_RETCODE GCGcolpoolFree(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_COLPOOL**        colpool             /**< pointer to store col pool */
   )
{
   assert(scip == (*colpool)->scip);
   assert(colpool != NULL);
   assert(*colpool != NULL);

   /* remove all cols from the pool */
   SCIP_CALL( GCGcolpoolClear(*colpool) );

   /* free clock */
   SCIPfreeClock(scip, &(*colpool)->poolclock);

   /* free hash table */
   SCIPhashtableFree(&(*colpool)->hashtable);

   SCIPfreeMemoryArrayNull(scip, &(*colpool)->cols);
   SCIPfreeMemory(scip, colpool);

   return SCIP_OKAY;
}

/** removes the col from the col pool */
static
SCIP_RETCODE colpoolDelCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col,                /**< col to remove */
   SCIP_Bool             free                /**< should the col be freed? */
   )
{
   int pos;

   assert(colpool != NULL);
   assert(colpool->firstunprocessed <= colpool->ncols);
   assert(colpool->firstunprocessedsol <= colpool->ncols);
   assert(col != NULL);

   pos = col->pos;
   assert(0 <= pos && pos < colpool->ncols);
   assert(colpool->cols[pos] == col);

   /* remove the col from the hash table */
   assert(SCIPhashtableExists(colpool->hashtable, (void*)col));
   SCIP_CALL( SCIPhashtableRemove(colpool->hashtable, (void*)col) );

   /* free the col */
   if( free )
      GCGfreeGcgCol(&colpool->cols[pos]);

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

   SCIPinfoMessage(colpool->scip, NULL, "clear colpool \n");

   /* free cols (in reverse order!) */
   for( i = colpool->ncols - 1; i >= 0; --i )
   {
      colpoolDelCol(colpool, colpool->cols[i], TRUE);
   }
   colpool->ncols = 0;

   return SCIP_OKAY;
}

/** if not already existing, adds row to col pool and captures it */
SCIP_RETCODE GCGcolpoolAddCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col,                /**< column to add */
   SCIP_Bool*            success             /**< pointer to store if col was added */
   )
{
   assert(colpool != NULL);
   assert(col != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check in hash table, if col already exists in the pool */
   if( SCIPhashtableRetrieve(colpool->hashtable, (void*)col) == NULL )
   {
      SCIP_CALL( GCGcolpoolAddNewCol(colpool, col) );
      *success = TRUE;
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
   SCIP_Bool             free                /**< should the col be freed? */
   )
{
   assert(colpool != NULL);
   assert(col != NULL);

   /* find the col in hash table */
   col = (GCG_COL*)SCIPhashtableRetrieve(colpool->hashtable, (void*)col);
   if( col == NULL )
   {
      SCIPerrorMessage("col %p is not existing in colpool %p\n", col, colpool);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( colpoolDelCol(colpool, col, free) );

   return SCIP_OKAY;
}


/** prices cols of the col pool */
SCIP_RETCODE GCGcolpoolPrice(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_PRICESTORE*       pricestore,         /**< GCG price storage */
   SCIP_SOL*             sol,                /**< solution to be separated (or NULL for LP-solution) */
   SCIP_Bool             colpoolisdelayed,   /**< is the colpool delayed (count cols found)? */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_Bool*            foundvars           /**< pointer to store the result of the separation call */
   )
{
   GCG_COL* col;
   SCIP_Bool found;
   int firstunproc;
   int oldncols;
   int c;

   assert(colpool != NULL);
   assert(colpool->firstunprocessed <= colpool->ncols);
   assert(colpool->firstunprocessedsol <= colpool->ncols);
   assert(foundvars != NULL);

   colpool->ncalls++;
   found = FALSE;

   SCIPdebugMessage("separating%s col pool %p with %d cols, beginning with col %d\n", ( sol == NULL ) ? "" : " solution from", (void*)colpool, colpool->ncols, firstunproc);

   /* start timing */
   SCIPstartClock(colpool->scip, colpool->poolclock);

   /* remember the current total number of found cols */
   //oldncols = getNCols();
   oldncols = 0;

   /* process all unprocessed cols in the pool */
   *foundvars = FALSE;

   for( c = colpool->ncols - 1; c >= 0; --c )
   {
      SCIP_Real redcost;

      col = colpool->cols[c];
      assert(col != NULL);
      assert(col->pos == c);

      redcost = GCGcolGetRedcost(col);

      /* TODO: use reduced costs? */
      if( SCIPisDualfeasNegative(scip, redcost) )
      {
         /* insert col in separation storage */
         SCIPdebugMessage(" -> col %p from the col pool (redcost: %g)\n",
            (void*)col, redcost );

         SCIP_CALL( GCGpricestoreAddCol(scip, pricestore, NULL, col, FALSE) );

         SCIP_CALL( colpoolDelCol(colpool, col, FALSE) );

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
   colpool->ncolsfound = 0; /*lint !e776*/
//   SCIPsepastoreGetNCols(sepastore) - oldncols
   /* stop timing */
   SCIPstopClock(colpool->scip, colpool->poolclock);

   return SCIP_OKAY;
}

/** gets array of cols in the col pool */
SCIP_RETCODE GCGcolpoolUpdateNode(
   GCG_COLPOOL*         colpool             /**< col pool */
   )
{
   assert(colpool != NULL);

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

      SCIP_CALL( GCGcomputeColMastercoefs(colpool->scip, col) );

      redcost = GCGcomputeRedCostGcgCol(colpool->scip, colpool->infarkas, col, NULL);

      SCIP_CALL( GCGcolUpdateRedcost(col, redcost, FALSE) );
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

   return SCIPclockGetTime(colpool->poolclock);
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

