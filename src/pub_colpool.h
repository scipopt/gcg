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

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicColMethods
 *
 * @{
 */

///** returns the ratio of LPs where the row belonging to this col was active in an LP solution, i.e.
// *  where the age of its row has not been increased
// *
// *  @see SCIPcolGetAge() to get the age of a col
// */
//EXTERN
//SCIP_Real SCIPcolGetLPActivityQuot(
//   GCG_COL*             col                 /**< col */
//   );

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

/** gets time in seconds used for separating cols from the pool */
EXTERN
SCIP_Real GCGcolpoolGetTime(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** get number of times, the col pool was separated */
EXTERN
SCIP_Longint GCGcolpoolGetNCalls(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** get total number of cols that were separated from the col pool */
EXTERN
SCIP_Longint GCGcolpoolGetNColsFound(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** creates col pool */
EXTERN
SCIP_RETCODE GCGcolpoolCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_COLPOOL**         colpool,            /**< pointer to store col pool */
   int                   agelimit,           /**< maximum age a col can reach before it is deleted from the pool */
   SCIP_Bool             globalcolpool       /**< is this the global col pool of SCIP? */
   );

/** frees col pool */
EXTERN
SCIP_RETCODE GCGcolpoolFree(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_COLPOOL**        colpool             /**< pointer to store col pool */
   );

/** removes all rows from the col pool */
EXTERN
SCIP_RETCODE GCGcolpoolClear(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** if not already existing, adds row to col pool and captures it */
EXTERN
SCIP_RETCODE GCGcolpoolAddCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col,                /**< column to add */
   SCIP_Bool*            success             /**< pointer to store if col was added */
   );

/** adds row to col pool and captures it; doesn't check for multiple cols */
EXTERN
SCIP_RETCODE GCGcolpoolAddNewCol(
   GCG_COLPOOL*         colpool,            /**< col pool */
   GCG_COL*             col                 /**< column to add */
   );

/** removes the LP row from the col pool */
EXTERN
SCIP_RETCODE GCGcolpoolDelCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col,                /**< col to remove */
   SCIP_Bool             free                /**< should the col be freed? */
   );
/* @} */

#ifdef __cplusplus
}
#endif

#endif
