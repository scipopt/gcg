/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
#include "type_colpool.h"
#include "type_gcgcol.h"
#include "type_pricestore_gcg.h"


#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicColMethods
 *
 * @{
 */

/** gets array of cols in the col pool */
EXTERN
GCG_COL** GCGcolpoolGetCols(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** get number of cols in the col pool */
EXTERN
int GCGcolpoolGetNCols(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** get maximum number of cols that were stored in the col pool at the same time */
EXTERN
int GCGcolpoolGetMaxNCols(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** gets time in seconds used for pricing cols from the pool */
EXTERN
SCIP_Real GCGcolpoolGetTime(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** get number of times, the col pool was separated */
EXTERN
SCIP_Longint GCGcolpoolGetNCalls(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** get total number of cols that were priced from the col pool */
EXTERN
SCIP_Longint GCGcolpoolGetNColsFound(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** creates col pool */
EXTERN
SCIP_RETCODE GCGcolpoolCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_COLPOOL**         colpool,            /**< pointer to store col pool */
   int                   agelimit            /**< maximum age a col can reach before it is deleted from the pool */
   );

/** frees col pool */
EXTERN
SCIP_RETCODE GCGcolpoolFree(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_COLPOOL**        colpool             /**< pointer to store col pool */
   );

/** removes all cols from the col pool */
EXTERN
SCIP_RETCODE GCGcolpoolClear(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** if not already existing, adds col to col pool and captures it */
EXTERN
SCIP_RETCODE GCGcolpoolAddCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col,                /**< column to add */
   SCIP_Bool*            success             /**< pointer to store if col was added */
   );

/** adds col to col pool and captures it; doesn't check for multiple cols */
EXTERN
SCIP_RETCODE GCGcolpoolAddNewCol(
   GCG_COLPOOL*         colpool,            /**< col pool */
   GCG_COL*             col                 /**< column to add */
   );

/** removes the col from the col pool */
EXTERN
SCIP_RETCODE GCGcolpoolDelCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col,                /**< col to remove */
   SCIP_Bool             freecol             /**< should the col be freed? */
   );

/** gets array of cols in the col pool */
EXTERN
SCIP_RETCODE GCGcolpoolUpdateNode(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** update reduced cost of columns in column pool */
EXTERN
SCIP_RETCODE GCGcolpoolUpdateRedcost(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** gets number of cols in the col pool */
EXTERN
void GCGcolpoolStartFarkas(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** gets number of cols in the col pool */
EXTERN
void GCGcolpoolEndFarkas(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** prices cols of the col pool */
EXTERN
SCIP_RETCODE GCGcolpoolPrice(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_PRICESTORE*       pricestore,         /**< GCG price storage */
   SCIP_SOL*             sol,                /**< solution to be separated (or NULL for LP-solution) */
   SCIP_Bool*            foundvars           /**< pointer to store the result of the separation call */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
