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

/**@file   struct_pricestore_gcg.h
 * @ingroup DATASTRUCTURES
 * @brief  datastructures for storing priced cols
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_STRUCT_PRICESTORE_H__
#define __GCG_STRUCT_PRICESTORE_H__


#include "scip/def.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "gcg/type_pricestore_gcg.h"
#include "gcg/pub_gcgcol.h"
#include "gcg/type_gcgcol.h"
#include "gcg/type_locks.h"

#ifdef __cplusplus
extern "C" {
#endif

/** storage for priced cols */
struct GCG_PriceStore
{
   GCG*                  gcg;                /**< GCG data structure */
   SCIP*                 masterprob;         /**< SCIP data structure */
   GCG_COL***            cols;               /**< array with priced cols sorted by score */
   SCIP_HASHTABLE*       hashtable;          /**< hashtable that maps the cols to their indices in the cols array */
   SCIP_Real**           objparallelisms;    /**< parallelism of col to the objective function */
   SCIP_Real**           orthogonalities;    /**< minimal orthogonality of col with all other cols of larger score */
   SCIP_Real**           scores;             /**< score for each priced col: weighted sum of efficacy and orthogonality */
   int*                  colssize;           /**< size of cols and score arrays */
   int*                  ncols;              /**< number of priced cols per problem */
   int*                  nforcedcols;        /**< number of forced priced cols (first positions in cols array) */
   int                   ncolstotal;         /**< number of priced cols (max. is set->price_maxcols) */
   int                   ncolsfound;         /**< total number of cols found so far */
   int                   ncolsfoundround;    /**< number of cols found so far in this pricing round */
   int                   ncolsapplied;       /**< total number of cols applied to the LPs */
   int                   narrays;            /**< number of allocated arrays (i.e., size of cols, scores, etc.) */
   SCIP_Bool             infarkas;           /**< is the price storage currently being filled with the columns from farkas pricing? */
   SCIP_Bool             forcecols;          /**< should the cols be used despite the number of cols parameter limit? */
   SCIP_Real             efficiacyfac;       /**< factor of efficiacy in score function */
   SCIP_Real             objparalfac;        /**< factor of objective parallelism in score function */
   SCIP_Real             orthofac;           /**< factor of orthogonalities in score function */
   SCIP_Real             mincolorth;         /**< minimal orthogonality of columns to add
                                                  (with respect to columns added in the current round) */
   GCG_EFFICIACYCHOICE   efficiacychoice;    /**< choice to base efficiacy on */

#ifdef _OPENMP
   GCG_LOCKS*            locks;              /**< OpenMP locks */
#endif
};

#ifdef __cplusplus
}
#endif

#endif
