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

/**@file   pricestore.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for storing pricerated cols
 * @author Tobias Achterberg
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

#include "pub_gcgcol.h"
#include "type_pricestore_gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates price storage */
extern
SCIP_RETCODE GCGpricestoreCreate(
   SCIP*                 scip,                /**< SCIP data structure */
   GCG_PRICESTORE**      pricestore           /**< pointer to store price storage */
   );

/** frees priceration storage */
extern
SCIP_RETCODE GCGpricestoreFree(
   SCIP*                 scip,                /**< SCIP data structure */
   GCG_PRICESTORE**      pricestore           /**< pointer to store price storage */
   );

/** informs price storage, that the setup in Farkas pricing starts now */
extern
void GCGpricestoreStartFarkas(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** informs price storage, that the setup in Farkas pricing is now finished */
extern
void GCGpricestoreEndFarkas(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** informs priceration storage, that the following cols should be used in any case */
extern
void GCGpricestoreStartForceCols(
   GCG_PRICESTORE*       pricestore           /**< priceration storage */
   );

/** informs priceration storage, that the following cols should no longer be used in any case */
extern
void GCGpricestoreEndForceCols(
   GCG_PRICESTORE*       pricestore           /**< priceration storage */
   );

/** adds col to priceration storage and captures it;
 *  if the col should be forced to enter the LP, an infinite score has to be used
 */
extern
SCIP_RETCODE GCGpricestoreAddCol(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_PRICESTORE*       pricestore,         /**< price storage */
   SCIP_SOL*             sol,                /**< primal solution that was pricerated, or NULL for LP solution */
   GCG_COL*              col,                /**< pricerated col */
   SCIP_Bool             forcecol            /**< should the col be forced to enter the LP? */
   );

/** adds cols to the LP and clears priceration storage */
extern
SCIP_RETCODE GCGpricestoreApplyCols(
   GCG_PRICESTORE*       pricestore,          /**< price storage */
   SCIP_Bool             root                 /**< are we at the root node? */
   );

/** clears the priceration storage without adding the cols to the LP */
extern
SCIP_RETCODE GCGpricestoreClearCols(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** removes cols that are inefficacious w.r.t. the current LP solution from priceration storage without adding the cols to the LP */
extern
SCIP_RETCODE GCGpricestoreRemoveInefficaciousCols(
   GCG_PRICESTORE*       pricestore,         /**< price storage */
   SCIP_Bool             root                /**< are we at the root node? */
   );

/** get cols in the priceration storage */
extern
GCG_COL** GCGpricestoreGetCols(
   GCG_PRICESTORE*       pricestore           /**< price storage */
   );

/** get number of cols in the priceration storage */
extern
int GCGpricestoreGetNCols(
   GCG_PRICESTORE*       pricestore           /**< priceration storage */
   );

/** get total number of cols found so far */
extern
int GCGpricestoreGetNColsFound(
   GCG_PRICESTORE*       pricestore           /**< priceration storage */
   );

/** get number of cols found so far in current priceration round */
extern
int GCGpricestoreGetNColsFoundRound(
   GCG_PRICESTORE*       pricestore           /**< priceration storage */
   );

/** get total number of cols applied to the LPs */
extern
int GCGpricestoreGetNColsApplied(
   GCG_PRICESTORE*       pricestore           /**< priceration storage */
   );

#ifdef __cplusplus
}
#endif

#endif
