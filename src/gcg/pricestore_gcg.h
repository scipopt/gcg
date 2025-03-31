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

/**@file   pricestore_gcg.h
 * @brief  methods for storing priced cols (based on SCIP's separation storage)
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_PRICESTORE_H__
#define __GCG_PRICESTORE_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_implics.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"
#include "scip/type_reopt.h"
#include "scip/type_branch.h"

#include "gcg/pub_colpool.h"
#include "gcg/pub_gcgcol.h"
#include "gcg/type_pricestore_gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates price storage */
GCG_EXPORT
SCIP_RETCODE GCGpricestoreCreate(
   GCG*                  gcg,                 /**< GCG data structure */
   GCG_PRICESTORE**      pricestore,          /**< pointer to store price storage */
   SCIP_Real             redcostfac,          /**< factor of -redcost/norm in score function */
   SCIP_Real             objparalfac,         /**< factor of objective parallelism in score function */
   SCIP_Real             orthofac,            /**< factor of orthogonalities in score function */
   SCIP_Real             mincolorth,          /**< minimal orthogonality of columns to add
                                                  (with respect to columns added in the current round) */
   GCG_EFFICIACYCHOICE   efficiacychoice,     /**< choice to base efficiacy on */
   int                   hashtablesize        /**< size of hashtable */
   );

/** frees price storage */
GCG_EXPORT
SCIP_RETCODE GCGpricestoreFree(
   GCG_PRICESTORE**      pricestore           /**< pointer to store price storage */
   );

/** informs price storage, that Farkas pricing starts now */
GCG_EXPORT
void GCGpricestoreStartFarkas(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** informs price storage, that Farkas pricing is now finished */
GCG_EXPORT
void GCGpricestoreEndFarkas(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** informs price storage, that the following cols should be used in any case */
GCG_EXPORT
void GCGpricestoreStartForceCols(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** informs price storage, that the following cols should no longer be used in any case */
GCG_EXPORT
void GCGpricestoreEndForceCols(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** adds col to price storage;
 *  if the col should be forced to enter the LP, an infinite score has to be used
 */
GCG_EXPORT
SCIP_RETCODE GCGpricestoreAddCol(
   GCG_PRICESTORE*       pricestore,         /**< price storage */
   GCG_COL*              col,                /**< priced col */
   SCIP_Bool             forcecol,           /**< should the col be forced to enter the LP? */
   SCIP_Bool             fromcolpool,        /**< is column from colpool */
   SCIP_Bool*            added               /**< pointer to var that indicates whether the col was added */
   );

/** adds cols to priced vars and clears price storage */
GCG_EXPORT
SCIP_RETCODE GCGpricestoreApplyCols(
   GCG_PRICESTORE*       pricestore,         /**< price storage */
   GCG_COLPOOL*          colpool,            /**< GCG column pool */
   SCIP_Bool             usecolpool,         /**< use column pool? */
   int*                  nfoundvars          /**< pointer to store number of variables that were added to the problem */
   );

/** clears the price storage without adding the cols to priced vars */
GCG_EXPORT
void GCGpricestoreClearCols(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** get cols in the price storage */
GCG_EXPORT
GCG_COL** GCGpricestoreGetCols(
   GCG_PRICESTORE*       pricestore,          /**< price storage */
   int                   arrayindex           /**< index of the arrays */
   );

/** get number of cols in the price storage */
GCG_EXPORT
int GCGpricestoreGetNCols(
   GCG_PRICESTORE*       pricestore,          /**< price storage */
   int                   arrayindex           /**< index of the arrays */
   );

/** get number of cols in the price storage */
GCG_EXPORT
int GCGpricestoreGetNColsTotal(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** get total number of cols found so far */
GCG_EXPORT
int GCGpricestoreGetNColsFound(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** get number of cols found so far in current price round */
GCG_EXPORT
int GCGpricestoreGetNColsFoundRound(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** get total number of cols applied to the LPs */
GCG_EXPORT
int GCGpricestoreGetNColsApplied(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

#ifdef __cplusplus
}
#endif

#endif
