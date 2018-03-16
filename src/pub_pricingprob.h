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

/**@file   pub_pricingprob.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for working with pricing problems
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PUB_PRICINGPROB_H__
#define GCG_PUB_PRICINGPROB_H__

#include "type_pricingprob.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * GCG Pricing Problem
 */

/**@defgroup GCG_PricingProb gcg pricingprob
 *
 * @{
 */


/** free all columns of a pricing job */
/* @todo: Move this to pricingjob.h ? */
EXTERN
void GCGpricingjobFreeCols(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** get the SCIP instance corresponding to the pricing job */
EXTERN
SCIP* GCGpricingjobGetPricingscip(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** get the index of the corresponding pricing problem */
EXTERN
int GCGpricingjobGetProbnr(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** return whether the pricing job is to be performed heuristically */
EXTERN
SCIP_Bool GCGpricingjobIsHeuristic(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** get the chunk of a pricing job */
EXTERN
SCIP_Real GCGpricingjobGetChunk(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** get the score of a pricing job */
EXTERN
SCIP_Real GCGpricingjobGetScore(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/** get the number of times the pricing job was performed during the loop */
EXTERN
int GCGpricingjobGetNSolves(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/* get the status of a pricing job */
EXTERN
SCIP_STATUS GCGpricingjobGetStatus(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/* get the lower bound of a pricing job */
EXTERN
SCIP_Real GCGpricingjobGetLowerbound(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/* get a column found by a pricing job */
EXTERN
GCG_COL* GCGpricingjobGetCol(
   GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
   int                   idx                 /**< index of a column */
   );

/* get the number of columns found by a pricing job */
EXTERN
int GCGpricingjobGetNCols(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/* get the number of improving columns found by a pricing job */
EXTERN
int GCGpricingjobGetNImpCols(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   );

/**@} */


#ifdef __cplusplus
}
#endif
#endif
