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

/**@file   colpool.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for storing cols in a col pool
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_COLPOOL_H__
#define __GCG_COLPOOL_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_sepastore.h"
#include "type_colpool.h"
#include "pub_colpool.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates col pool */
extern
SCIP_RETCODE GCGcolpoolCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_COLPOOL**         colpool,            /**< pointer to store col pool */
   int                   agelimit,           /**< maximum age a col can reach before it is deleted from the pool */
   SCIP_Bool             globalcolpool       /**< is this the global col pool of SCIP? */
   );

/** frees col pool */
extern
SCIP_RETCODE GCGcolpoolFree(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_COLPOOL**        colpool             /**< pointer to store col pool */
   );

/** removes all rows from the col pool */
extern
SCIP_RETCODE GCGcolpoolClear(
   GCG_COLPOOL*         colpool             /**< col pool */
   );

/** if not already existing, adds row to col pool and captures it */
extern
SCIP_RETCODE GCGcolpoolAddCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col                 /**< column to add */
   );

/** adds row to col pool and captures it; doesn't check for multiple cols */
extern
SCIP_RETCODE GCGcolpoolAddNewCol(
   GCG_COLPOOL*         colpool,            /**< col pool */
   GCG_COL*             col                 /**< column to add */
   );

/** removes the LP row from the col pool */
extern
SCIP_RETCODE GCGcolpoolDelCol(
   GCG_COLPOOL*          colpool,            /**< col pool */
   GCG_COL*              col                 /**< col to remove */
   );

/** prices cols of the col pool */
extern
SCIP_RETCODE GCGcolpoolPrice(
   GCG_COLPOOL*          colpool,            /**< col pool */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< solution to be separated (or NULL for LP-solution) */
   SCIP_Bool             colpoolisdelayed,   /**< is the colpool delayed (count cols found)? */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_Bool*            foundvars           /**< pointer to store the result of the separation call */
   );

#ifdef __cplusplus
}
#endif

#endif
