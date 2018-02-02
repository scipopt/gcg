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

/**@file   class_colpool.cpp
 * @brief  class with functions for colpool
 * @author Jonas Witt
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_colpool.h"
#include "pricer_gcg.h"
#include "gcg.h"
#include "sepa_master.h"
#include "objscip/objscip.h"
#include "scip/cons_linear.h"
#include "gcgpqueue.h"
#include "pub_gcgpqueue.h"
#include "pub_gcgcol.h"

#include <exception>

#define SCIP_CALL_EXC(x)   do                                                                                 \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) !=  SCIP_OKAY )                                                \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             throw std::exception();                                                          \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

namespace gcg {

   Colpool::Colpool(
      SCIP*             scip_,               /**< SCIP data structure */
      int               agelimit_,           /**< maximum age a column can reach before it is deleted from the pool */
      int               maxncolssoft_,       /**< soft maximal number of columns stored in the pool at the same time */
      int               maxncolshard_        /**< hard maximal number of columns stored in the pool at the same time */
      ):scip(scip_), pqueue((GCG_PQueue*) NULL), agelimit(agelimit_), maxncolssoft(maxncolssoft_), maxncolshard(maxncolshard_), nodenr(-1)
   {
      SCIP_CALL_EXC( GCGpqueueCreate(&pqueue, maxncolshard, (SCIP_Real) sizeof(GCG_COL*), GCGcolCompRedcost) );
   }
   Colpool::~Colpool()
   {
      GCG_COL** gcgcols;
      int ngcgcols;

      int i;

      ngcgcols = getNCols();
      gcgcols = getCols();

      for( i = 0; i < ngcgcols; ++i )
      {
         GCGfreeGcgCol(&(gcgcols[i]));
      }
      GCGpqueueFree(&pqueue);
      pqueue = NULL;
   }

   /** add gcg column to column pool */
   SCIP_RETCODE Colpool::addCol(
      GCG_COL*          gcgcol,             /**< gcg column to add */
      SCIP_Bool*        success             /**< bool returns if colum was succesfully added
                                                 (number of columns is not bigger than maxncols or column already exists) */
   )
   {
      *success = FALSE;

      if( GCGpqueueNElems(pqueue) >= maxncolshard )
      {
         return SCIP_OKAY;
      }


      if( !existsCol(gcgcol) )
      {
         SCIP_CALL( GCGpqueueInsert(pqueue, (void*) gcgcol) );
         *success = TRUE;
      }

      return SCIP_OKAY;
   }

   /** return if column already exists in column pool */
   SCIP_Bool Colpool::existsCol(
      GCG_COL*          newcol
      )
   {
      GCG_COL** cols;
      int ncols;
      int i;

      cols = getCols();
      ncols = getNCols();

      for( i = 0; i < ncols; ++i)
      {
         GCG_COL* col;

         col = cols[i];

         if( GCGcolIsEq(newcol, col) )
         {
            return TRUE;
         }
      }

      return FALSE;
   }

   /**< get best column in column pool and remove it from column pool */
   SCIP_RETCODE Colpool::getBestCol(
      GCG_COL**         gcgcol              /**< pointer to store gcg column */
   )
   {
      *gcgcol = (GCG_COL*) GCGpqueueRemove(pqueue);

      return SCIP_OKAY;
   }


   /**< get best column's reduced cost */
   SCIP_Real Colpool::getBestColRedcost()
   {
      GCG_COL* gcgcol;
      SCIP_Real redcost;

      gcgcol = (GCG_COL*) GCGpqueueFirst(pqueue);

      if( gcgcol != NULL )
      {
         redcost = GCGcolGetRedcost(gcgcol);
      }
      else
      {
         redcost = SCIPinfinity(scip);
      }

      return redcost;
   }

   /**< get best column's probnr (or -1 if colpool is empty) */
   int Colpool::getBestColProbNr()
   {
      GCG_COL* gcgcol;
      int probnr;

      gcgcol = (GCG_COL*) GCGpqueueFirst(pqueue);

      if( gcgcol != NULL)
      {
         probnr = GCGcolGetProbNr(gcgcol);
      }
      else
      {
         probnr = -1;
      }
      return probnr;
   }

   /**< get reduced cost of column at specific postition */
   SCIP_Real Colpool::getColRedcost(
      int               pos                 /**< position of column */
   )
   {
      GCG_COL** gcgcols;

      gcgcols = getCols();

#ifdef SCIP_DEBUG
      int ngcgcols;


      ngcgcols = getNCols();

      assert(0 <= pos && pos < ngcgcols);
#endif

      return GCGcolGetRedcost(gcgcols[pos]);
   }

   /**< get age of column at specific postition */
   int Colpool::getColAge(
      int               pos                 /**< position of column */
   )
   {
      GCG_COL** gcgcols;

      gcgcols = getCols();

#ifdef SCIP_DEBUG
      int ngcgcols;


      ngcgcols = getNCols();

      assert(0 <= pos && pos < ngcgcols);
#endif

      return GCGcolGetAge(gcgcols[pos]);
   }

   /**< get columns in column pool */
   GCG_COL** Colpool::getCols()
   {
      return (GCG_COL**) GCGpqueueElems(pqueue);
   }

   /**< get number of columns in column pool */
   int Colpool::getNCols()
   {
      return GCGpqueueNElems(pqueue);
   }

   /**< delete all columns that are older than agelimit
    * WARNING: This method changes the order in which the colums are stored.
    * Use GCGpqueueResort() to resort the columns by reduced cost again */
   SCIP_RETCODE Colpool::deleteOldColumns()
   {
      /* todo: get comperator of pqueue */

      SCIP_CALL( GCGpqueueSetComperator(pqueue, GCGcolCompAge) );

      SCIP_CALL( GCGpqueueResort(pqueue) );

      while( GCGpqueueNElems(pqueue) > 0 )
      {
         GCG_COL* gcgcol;

         gcgcol = (GCG_COL*) GCGpqueueFirst(pqueue);

         if( GCGcolGetAge(gcgcol) > agelimit)
         {
            gcgcol = (GCG_COL*) GCGpqueueRemove(pqueue);

            GCGfreeGcgCol(&gcgcol);
         }
         else
            break;
      }

      /* todo: use previous comperator of pqueue */
      SCIP_CALL( GCGpqueueSetComperator(pqueue, GCGcolCompRedcost) );

      return SCIP_OKAY;
   }

   /**< delete the oldest columns such that number of columns in colpool is
    *   lower than or equal to maxncolssoft
    * WARNING: This method changes the order in which the colums are stored.
    * Use GCGpqueueResort() to resort the columns by reduced cost again  */
   SCIP_RETCODE Colpool::deleteOldestColumns()
   {
      if( GCGpqueueNElems(pqueue) <= maxncolssoft )
      {
         return SCIP_OKAY;
      }

      if( maxncolssoft == 0 )
      {
         SCIP_CALL( deleteAllColumns() );

         return SCIP_OKAY;
      }

      /* todo: get comperator of pqueue */

      SCIP_CALL( GCGpqueueSetComperator(pqueue, GCGcolCompAge) );

      SCIP_CALL( GCGpqueueResort(pqueue) );

      while( GCGpqueueNElems(pqueue) > maxncolssoft )
      {
         GCG_COL* gcgcol;

         gcgcol = (GCG_COL*) GCGpqueueRemove(pqueue);

         GCGfreeGcgCol(&gcgcol);
      }

      /* todo: use previous comperator of pqueue */
      SCIP_CALL( GCGpqueueSetComperator(pqueue, GCGcolCompRedcost) );

      return SCIP_OKAY;
   }

   /**< delete all columns in colpool */
   SCIP_RETCODE Colpool::deleteAllColumns()
   {
      GCG_COL** cols;

      int ncols;
      int i;

      ncols = GCGpqueueNElems(pqueue);
      cols = (GCG_COL**) GCGpqueueElems(pqueue);

      for(i = 0; i < ncols; ++i)
      {
         GCGfreeGcgCol(&(cols[i]));
      }

      GCGpqueueClear(pqueue);

      return SCIP_OKAY;
   }

   /**< resort columns (after reduce cost have changed) */
   SCIP_RETCODE Colpool::resortColumns()
   {
      SCIP_CALL( GCGpqueueResort(pqueue) );

      return SCIP_OKAY;
   }

   SCIP_RETCODE Colpool::setSoftlimit(
      int               newsoftlimit
   )
   {
      maxncolssoft = newsoftlimit;

      return SCIP_OKAY;
   }

   SCIP_RETCODE Colpool::updateNode(
   )
   {
      if( nodenr < 0 )
      {
         nodenr = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
      }
      else if( nodenr != SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) )
      {
         SCIP_CALL( deleteAllColumns() );

         nodenr = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
      }

      return SCIP_OKAY;
   }

} /* namespace gcg */
