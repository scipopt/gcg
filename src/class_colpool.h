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

/**@file   class_colpool.h
 * @brief  class with functions for column pool
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CLASS_COLPOOL_H__
#define GCG_CLASS_COLPOOL_H__

#include "objscip/objscip.h"
#include "class_pricingtype.h"
#include "type_gcgcol.h"
#include "gcgpqueue.h"
#include "pub_gcgpqueue.h"

namespace gcg {

class Colpool
{   /*lint -esym(1712,Colpool)*/

private:
   SCIP*                scip;               /**< SCIP data structure */
   GCG_PQUEUE*          pqueue;             /**< priority queue for storing columns */
   int                  agelimit;           /**< maximum age a column can reach before it is deleted from the pool */
   int                  maxncolssoft;       /**< soft maximal number of columns stored in the pool at the same time */
   int                  maxncolshard;       /**< hard maximal number of columns stored in the pool at the same time */
   int                  nodenr;             /**< node at which columns in colpool respect branching decisions */

public:

   /** constructor */
   Colpool(
      SCIP*             scip,               /**< SCIP data structure */
      int               agelimit,           /**< maximum age a column can reach before it is deleted from the pool */
      int               maxncolssoft,       /**< soft maximal number of columns stored in the pool at the same time */
      int               maxncolshard        /**< hard maximal number of columns stored in the pool at the same time */
      );

   ~Colpool();

   /** add gcg column to column pool */
   SCIP_RETCODE addCol(
      GCG_COL*          gcgcol,             /**< gcg column to add */
      SCIP_Bool*        success             /**< bool returns if colum was succesfully added (number of columns is not bigger than maxncols) */
   );

   /** return if column already exists in column pool */
   SCIP_Bool existsCol(
      GCG_COL*          gcgcol
      );

   /**< get best column in column pool and remove it from column pool */
   SCIP_RETCODE getBestCol(
      GCG_COL**         gcgcol              /**< pointer to store gcg column */
      );

   /**< get best column's reduced cost */
   SCIP_Real getBestColRedcost();

   /**< get best column's probnr */
   int getBestColProbNr();

   /**< get age of column at specific postition */
   int getColAge(
      int               pos                 /**< position of column */
   );

   /**< get reduced cost of column at specific postition */
   SCIP_Real getColRedcost(
      int               pos                 /**< position of column */
      );

   /**< get number of columns in column pool */
   int getNCols();

   /**< delete all columns that are older than agelimit
    * WARNING: This method changes the order in which the colums are stored.
    * Use GCGpqueueResort() to resort the columns by reduced cost again */
   SCIP_RETCODE deleteOldColumns();

   /**< delete the oldest columns such that number of columns in colpool is
    *   lower than or equal to maxncolssoft
    * WARNING: This method changes the order in which the colums are stored.
    * Use GCGpqueueResort() to resort the columns by reduced cost again  */
   SCIP_RETCODE deleteOldestColumns();

   /**< delete all columns in colpool */
   SCIP_RETCODE deleteAllColumns();

   /**< resort columns (after reduce cost have changed) */
   SCIP_RETCODE resortColumns();

   /**< get columns in column pool */
   GCG_COL** getCols();

   SCIP_RETCODE setSoftlimit(
      int               newsoftlimit
      );

   SCIP_RETCODE updateNode(
      );

private:



};

} /* namespace gcg */
#endif /* GCG_CLASS_COLPOOL_H__ */
