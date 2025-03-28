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

/**@file   pub_colpool.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for storing cols in a col pool
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_COLPOOL_H__
#define __SCIP_PUB_COLPOOL_H__



#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_retcode.h"
#include "scip/type_sol.h"
#include "gcg/def.h"
#include "gcg/type_gcgcol.h"
#include "gcg/type_colpool.h"
#include "gcg/type_pricestore_gcg.h"
#include "gcg/type_gcg.h"


#ifdef __cplusplus
extern "C" {
#endif

/**@defgroup GCG_COLPOOL GCG Column Pool
 * @ingroup DATASTRUCTURES
 * @{
 */

/** gets array of cols in the col pool */
GCG_EXPORT
GCG_COL** GCGcolpoolGetCols(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** get number of cols in the col pool */
GCG_EXPORT
int GCGcolpoolGetNCols(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** get maximum number of cols that were stored in the col pool at the same time */
GCG_EXPORT
int GCGcolpoolGetMaxNCols(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** gets time in seconds used for pricing cols from the pool */
GCG_EXPORT
SCIP_Real GCGcolpoolGetTime(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** get number of times, the col pool was separated */
GCG_EXPORT
SCIP_Longint GCGcolpoolGetNCalls(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** get total number of cols that were priced from the col pool */
GCG_EXPORT
SCIP_Longint GCGcolpoolGetNColsFound(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** creates col pool */
GCG_EXPORT
SCIP_RETCODE GCGcolpoolCreate(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_COLPOOL**         colpool,            /**< pointer to store col pool */
   int                   agelimit            /**< maximum age a col can reach before it is deleted from the pool */
   );

/** frees col pool */
GCG_EXPORT
SCIP_RETCODE GCGcolpoolFree(
   GCG_COLPOOL**        colpool             /**< pointer to store col pool */
   );

/** removes all cols from the col pool */
GCG_EXPORT
SCIP_RETCODE GCGcolpoolClear(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** if not already existing, adds col to col pool and captures it */
GCG_EXPORT
SCIP_RETCODE GCGcolpoolAddCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col,                /**< column to add */
   SCIP_Bool             freeduplicate       /**< shouldl the col be freed if it is a duplicate? */
   );

/** adds col to col pool and captures it; doesn't check for multiple cols */
GCG_EXPORT
SCIP_RETCODE GCGcolpoolAddNewCol(
   GCG_COLPOOL*         colpool,            /**< col pool */
   GCG_COL*             col                 /**< column to add */
   );

/** removes the col from the col pool */
GCG_EXPORT
SCIP_RETCODE GCGcolpoolDelCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col,                /**< col to remove */
   SCIP_Bool             freecol             /**< should the col be freed? */
   );

/** update node at which columns of column pool are feasible */
GCG_EXPORT
SCIP_RETCODE GCGcolpoolUpdateNode(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** update reduced cost of columns in column pool */
GCG_EXPORT
SCIP_RETCODE GCGcolpoolUpdateRedcost(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** gets number of cols in the col pool */
GCG_EXPORT
void GCGcolpoolStartFarkas(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** gets number of cols in the col pool */
GCG_EXPORT
void GCGcolpoolEndFarkas(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** prices cols of the col pool */
GCG_EXPORT
SCIP_RETCODE GCGcolpoolPrice(
   GCG_COLPOOL*          colpool,            /**< col pool */
   SCIP_SOL*             sol,                /**< solution to be separated (or NULL for LP-solution) */
   int*                  nfoundvars          /**< pointer to store the result of the separation call */
   );

/** removes cols that violate global bounds */
GCG_EXPORT
SCIP_RETCODE GCGcolpoolPropagateGlobalBounds(
   GCG_COLPOOL*          colpool             /**< col pool */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
