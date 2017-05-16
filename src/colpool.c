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
 * @brief  methods for storing cols in a col pool
 * @author Tobias Achterberg
 * @author Stefan Heinz
 * @author Gerald Gamrath
 * @author Marc Pfetsch
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

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
   SCIP_ROW* row1;
   SCIP_ROW* row2;

   row1 = (SCIP_ROW*)key1;
   row2 = (SCIP_ROW*)key2;
   assert(row1 != NULL);
   assert(row2 != NULL);

   /* Sort the column indices of both rows.
    *
    * The columns in a row are divided into two parts: LP columns, which are currently in the LP and non-LP columns;
    * we sort the rows, but that only ensures that within these two parts, columns are sorted w.r.t. their index.
    * Normally, this should be suficient, because a column contained in both rows should either be one of the LP columns
    * for both or one of the non-LP columns for both.
    * However, directly after a row was created, before it is added to the LP, the row is not linked to all its
    * columns and all columns are treated as non-LP columns.
    * Therefore, if exactly one of the rows has no LP columns, we cannot rely on the partition, because this row might
    * just have been created and also columns that are in the LP might be in the non-LP columns part.
    */
   SCIProwSort(row1);
   SCIProwSort(row2);
   assert(row1->lpcolssorted);
   assert(row1->nonlpcolssorted);
   assert(row1->validminmaxidx);
   assert(row2->lpcolssorted);
   assert(row2->nonlpcolssorted);
   assert(row2->validminmaxidx);

   /* currently we are only handling rows which are completely linked or not linked at all */
   assert(row1->nunlinked == 0 || row1->nlpcols == 0);
   assert(row2->nunlinked == 0 || row2->nlpcols == 0);

   /* compare the trivial characteristics of the rows */
   if( row1->len != row2->len
      || row1->minidx != row2->minidx
      || row1->maxidx != row2->maxidx
      || row1->nummaxval != row2->nummaxval
      || row1->numminval != row2->numminval
      || REALABS(row1->lhs - row2->lhs) > SCIP_DEFAULT_EPSILON
      || REALABS(row1->rhs - row2->rhs) > SCIP_DEFAULT_EPSILON
      || REALABS(row1->maxval - row2->maxval) > SCIP_DEFAULT_EPSILON
      || REALABS(row1->minval - row2->minval) > SCIP_DEFAULT_EPSILON
       )
      return FALSE;

   /* both rows have LP columns, or none of them has, or one has only LP colums and the other only non-LP columns,
    * so we can rely on the sorting of the columns
    */
   if( (row1->nlpcols == 0) == (row2->nlpcols == 0)
      || (row1->nlpcols == 0 && row2->nlpcols == row2->len)
      || (row1->nlpcols == row1->len && row2->nlpcols == 0) )
   {
      int i;

      if( (row1->nlpcols == 0) == (row2->nlpcols == 0) )
      {
#ifndef NDEBUG
         /* in debug mode, we check that we can rely on the partition into LP columns and non-LP columns */
         int i2;

         i = 0;
         i2 = row2->nlpcols;
         while( i < row1->nlpcols && i2 < row2->len )
         {
            assert(row1->cols[i] != row2->cols[i2]);
            if( row1->cols[i]->index < row2->cols[i2]->index )
               ++i;
            else
            {
               assert(row1->cols[i]->index > row2->cols[i2]->index);
               ++i2;
            }
         }
         assert(i == row1->nlpcols || i2 == row2->len);

         i = row1->nlpcols;
         i2 = 0;
         while( i < row1->len && i2 < row2->nlpcols )
         {
            assert(row1->cols[i] != row2->cols[i2]);
            if( row1->cols[i]->index < row2->cols[i2]->index )
               ++i;
            else
            {
               assert(row1->cols[i]->index > row2->cols[i2]->index);
               ++i2;
            }
         }
         assert(i == row1->len || i2 == row2->nlpcols);
#endif

         /* both rows are linked and the number of lpcolumns is not equal so they cannot be equal */
         if( row1->nlpcols != row2->nlpcols )
            return FALSE;
      }

      /* compare the columns of the rows */
      for( i = 0; i < row1->len; ++i )
      {
         if( row1->cols[i] != row2->cols[i] )
            return FALSE;
      }

      /* compare the coefficients of the rows */
      for( i = 0; i < row1->len; ++i )
      {
         if( REALABS(row1->vals[i] - row2->vals[i]) > SCIP_DEFAULT_EPSILON )
            return FALSE;
      }
   }
   /* one row has LP columns, but the other not, that could be because the one without was just created and isn't
    * linked yet; in this case, one column could be an LP column in one row and a non-LP column in the other row, so we
    * cannot rely on the partition; thus, we iteratively check whether the next column of row1 is either the next LP
    * column of row2 or the next non-LP column of row2 and the coefficients are equal
    */
   else
   {
      int i1;
      int ilp;
      int inlp;

      /* ensure that row1 is the row without LP columns, switch the rows, if neccessary */
      if( row2->nlpcols == 0 )
      {
         SCIP_ROW* tmprow;
         tmprow = row2;
         row2 = row1;
         row1 = tmprow;
      }
      assert(row1->nlpcols == 0 && row2->nlpcols > 0);

      ilp = 0;
      inlp = row2->nlpcols;

      /* compare the columns and coefficients of the rows */
      for( i1 = 0; i1 < row1->len; ++i1 )
      {
         /* current column of row1 is the current LP column of row2, check the coefficient */
         if( ilp < row2->nlpcols && row1->cols[i1] == row2->cols[ilp] )
         {
            if( REALABS(row1->vals[i1] - row2->vals[ilp]) > SCIP_DEFAULT_EPSILON )
               return FALSE;
            else
               ++ilp;
         }
         /* current column of row1 is the current non-LP column of row2, check the coefficient */
         else if( inlp < row2->len && row1->cols[i1] == row2->cols[inlp] )
         {
            if( REALABS(row1->vals[i1] - row2->vals[inlp]) > SCIP_DEFAULT_EPSILON )
               return FALSE;
            else
               ++inlp;
         }
         /* current column of row1 is neither the current LP column of row2, nor the current non-LP column of row 2 */
         else
            return FALSE;
      }
   }

   return TRUE;
}

static
SCIP_DECL_HASHKEYVAL(hashKeyValCol)
{  /*lint --e{715}*/
   SCIP_ROW* row;
   unsigned int keyval;
   SCIP_Real maxval;
   SCIP_Real minval;
   SCIP_SET* set;

   set = (SCIP_SET*) userptr;
   row = (SCIP_ROW*)key;
   assert(row != NULL);

   maxval = SCIProwGetMaxval(row, set);
   minval = SCIProwGetMinval(row, set);
   assert(row->nummaxval > 0);
   assert(row->numminval > 0);
   assert(row->validminmaxidx);

   keyval = SCIPhashFour(SCIPpositiveRealHashCode(maxval, 8),
                         SCIPpositiveRealHashCode(minval, 12),
                         SCIPcombineThreeInt(row->maxidx, row->len, row->minidx),
                         SCIPcombineTwoInt(row->nummaxval, row->numminval));

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

/** removes all rows from the col pool */
SCIP_RETCODE GCGcolpoolClear(
   GCG_COLPOOL*          colpool             /**< col pool */
   )
{
   int i;

   assert(colpool != NULL);

   /* free cols */
   for( i = 0; i < colpool->ncols; ++i )
   {
      GCGfreeGcgCol(&colpool->cols[i]);
   }
   colpool->ncols = 0;

   return SCIP_OKAY;
}

/** if not already existing, adds row to col pool and captures it */
SCIP_RETCODE GCGcolpoolAddCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col                 /**< column to add */
   )
{
   assert(colpool != NULL);
   assert(col != NULL);

   /* check in hash table, if col already exists in the pool */
   if( SCIPhashtableRetrieve(colpool->hashtable, (void*)col) == NULL )
   {
      SCIP_CALL( GCGcolpoolAddNewCol(colpool, col) );
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

/** removes the col from the col pool */
static
SCIP_RETCODE colpoolDelCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col                 /**< col to remove */
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

/** removes the LP row from the col pool */
SCIP_RETCODE GCGcolpoolDelCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col                 /**< col to remove */
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

   SCIP_CALL( colpoolDelCol(colpool, col) );

   return SCIP_OKAY;
}


/** prices cols of the col pool */
SCIP_RETCODE GCGcolpoolPrice(
   GCG_COLPOOL*          colpool,            /**< col pool */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
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
//   oldncols = SCIPsepastoreGetNCols(sepastore);
   oldncols = 0;

   /* process all unprocessed cols in the pool */
   *foundvars = FALSE;
   for( c = firstunproc; c < colpool->ncols; ++c )
   {
      SCIP_Longint proclp;

      col = colpool->cols[c];
      assert(col != NULL);
      assert(col->pos == c);

//      if( proclp < stat->lpcount )
//      {
//         SCIP_ROW* row;
//
//         if ( sol == NULL )
//            col->processedlp = stat->lpcount;
//         else
//            col->processedlpsol = stat->lpcount;
//
//         row = col->row;
//         if( !SCIProwIsInLP(row) )
//         {
//            /* TODO: use reduced costs? */
//            if( (sol == NULL && SCIProwIsLPEfficacious(row, set, stat, lp, root)) || (sol != NULL && SCIProwIsSolEfficacious(row, set, stat, sol, root)) )
//            {
//               /* insert col in separation storage */
//               SCIPdebugMessage(" -> separated col <%s> from the col pool (feasibility: %g)\n",
//                  SCIProwGetName(row), ( sol == NULL ) ? SCIProwGetLPFeasibility(row, set, stat, lp) : SCIProwGetSolFeasibility(row, set, stat, sol) );
//               SCIP_CALL( SCIPsepastoreAddCol(sepastore, blkmem, set, stat, eventqueue, eventfilter, lp, sol, row, FALSE, root, &foundvars) );
//
//               col->age = 0;
//            }
//            else
//            {
//               col->age++;
//               if( colIsAged(col, colpool->agelimit) )
//               {
//                  SCIP_CALL( colpoolDelCol(colpool, col) );
//               }
//            }
//         }
//      }
   }

   /* update the number of found cols */
   colpool->ncolsfound += SCIPsepastoreGetNCols(sepastore) - oldncols; /*lint !e776*/

   /* stop timing */
   SCIPstopClock(colpool->scip, colpool->poolclock);

   return SCIP_OKAY;
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

