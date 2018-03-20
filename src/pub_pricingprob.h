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


/** get the SCIP instance corresponding to the pricing problem */
EXTERN
SCIP* GCGpricingprobGetPricingscip(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** get the index of the corresponding pricing problem */
EXTERN
int GCGpricingprobGetProbnr(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** get the number of times the pricing problem was solved during the loop */
EXTERN
int GCGpricingprobGetNSolves(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/* get the status of a pricing problem */
EXTERN
SCIP_STATUS GCGpricingprobGetStatus(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/* get the lower bound of a pricing problem */
EXTERN
SCIP_Real GCGpricingprobGetLowerbound(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/* get the best column found by solving a particular pricing problem */
EXTERN
GCG_COL* GCGpricingprobGetBestCol(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/* get the best reduced cost of a pricing problem */
EXTERN
SCIP_Real GCGpricingprobGetBestRedcost(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/* get the number of columns found for this pricing problem */
EXTERN
int GCGpricingprobGetNCols(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/* get the number of improving columns found for this pricing problem */
EXTERN
int GCGpricingprobGetNImpCols(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/* get the total number of improving colums found in the last pricing rounds */
EXTERN
int GCGpricingprobGetNColsLastRounds(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
   );

/**@} */


#ifdef __cplusplus
}
#endif
#endif
