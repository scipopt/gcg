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

/**@file   pricingjob.h
 * @brief  private methods for working with pricing jobs, to be used by the pricing controller only
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PRICINGJOB_H__
#define GCG_PRICINGJOB_H__

#include "struct_pricingjob.h"
#include "type_pricingjob.h"
#include "type_gcgcol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create a pricing job */
EXTERN
SCIP_RETCODE GCGpricingjobCreate(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGJOB**      pricingjob,         /**< pricing job to be created */
   SCIP*                 pricingscip,        /**< SCIP data structure of the corresponding pricing problem */
   int                   probnr,             /**< index of the corresponding pricing problem */
   int                   chunk,              /**< chunk that the pricing problem should belong to */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
);

/** free a pricing job */
EXTERN
void GCGpricingjobFree(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGJOB**      pricingjob          /**< pricing job to be freed */
);

/** setup a pricing job at the beginning of the pricing loop */
EXTERN
SCIP_RETCODE GCGpricingjobSetup(
   SCIP*                 scip,               /**< master SCIP instance */
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   SCIP_Bool             heuristic,          /**< shall the pricing job be performed heuristically? */
   int                   maxcolsprob,        /**< maximum number of columns that the problem should be looking for */
   int                   scoring,            /**< scoring parameter */
   int                   nroundscol,         /**< number of previous pricing rounds for which the number of improving columns should be counted */
   SCIP_Real             dualsolconv,        /**< dual solution value of corresponding convexity constraint */
   int                   npointsprob,        /**< total number of extreme points generated so far by the pricing problem */
   int                   nraysprob,          /**< total number of extreme rays generated so far by the pricing problem */
   int                   maxcols             /**< maximum number of columns to be generated in total */
   );

/** update a pricing job after the pricing problem has been solved */
EXTERN
SCIP_RETCODE GCGpricingjobUpdate(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   SCIP_STATUS           status,             /**< status after solving the pricing problem */
   SCIP_Real             lowerbound,         /**< lower bound returned by the pricing problem */
   GCG_COL**             cols,               /**< columns found by the last solving of the pricing problem */
   int                   ncols               /**< number of columns found */
   );

/** update solving statistics of a pricing job */
EXTERN
void GCGpricingjobUpdateSolvingStats(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** increase the solution limit of a pricing job */
EXTERN
SCIP_RETCODE GCGpricingjobIncreaseSollimit(
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   int                   maxcolsprob         /**< maximum number of columns that the problem should be looking for */
   );

/** set the pricing job to be performed heuristically */
EXTERN
void GCGpricingjobSetHeuristic(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** set the pricing job to be performed exactly */
EXTERN
void GCGpricingjobSetExact(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** set the lower bound of a pricing job */
EXTERN
void GCGpricingjobSetLowerbound(
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   SCIP_Real             lowerbound          /**< new lower bound */
   );

/* get a column array of a pricing job */
EXTERN
GCG_COL** GCGpricingjobGetCols(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/* set the number of columns found by a pricing job */
EXTERN
void GCGpricingjobSetNCols(
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   int                   ncols               /**< number of columns */
   );

/* set the number of improving columns found by a pricing job */
EXTERN
void GCGpricingjobSetNImpCols(
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   int                   nimpcols            /**< number of improving columns */
   );

/* update numbers of improving columns over the last pricing rounds */
EXTERN
void GCGpricingjobUpdateNColsround(
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
   );

/* get the number of heuristic pricing iterations of the pricing job */
EXTERN
int GCGpricingjobGetNHeurIters(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

#ifdef __cplusplus
}
#endif

#endif
