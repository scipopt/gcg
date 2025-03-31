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

/**@file   pub_pricingprob.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for working with pricing problems
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PUB_PRICINGPROB_H__
#define GCG_PUB_PRICINGPROB_H__


#include "gcg/type_pricingprob.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * GCG Pricing Problem
 */

/**
 * @ingroup PRICINGPROB
 * @{
 */


/** get the SCIP instance corresponding to the pricing problem */
GCG_EXPORT
SCIP* GCGpricingprobGetPricingscip(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** get the index of the corresponding pricing problem */
GCG_EXPORT
int GCGpricingprobGetProbnr(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** get generic branching data corresponding to the pricing problem */
GCG_EXPORT
void GCGpricingprobGetGenericBranchData(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   SCIP_CONS***          branchconss,        /**< pointer to store branching constraints array, or NULL */
   SCIP_Real**           branchduals,        /**< pointer to store array of corresponding dual values, or NULL */
   int*                  nbranchconss        /**< pointer to store number of generic branching constraints, or NULL */
   );

/** get the number of generic branching constraints corresponding to the pricing problem */
GCG_EXPORT
int GCGpricingprobGetNGenericBranchconss(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** get index of current generic branching constraint considered the pricing problem */
GCG_EXPORT
int GCGpricingprobGetBranchconsIdx(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** check if the current generic branching constraint has already been added */
GCG_EXPORT
SCIP_Bool GCGpricingprobBranchconsIsAdded(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** mark the current generic branching constraint to be added */
GCG_EXPORT
void GCGpricingprobMarkBranchconsAdded(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** get the status of a pricing problem */
GCG_EXPORT
GCG_PRICINGSTATUS GCGpricingprobGetStatus(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** get the lower bound of a pricing problem */
GCG_EXPORT
SCIP_Real GCGpricingprobGetLowerbound(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** get the number of improving columns found for this pricing problem */
GCG_EXPORT
int GCGpricingprobGetNImpCols(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** get the number of times the pricing problem was solved during the loop */
GCG_EXPORT
int GCGpricingprobGetNSolves(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   );

/** get the total number of improving colums found in the last pricing rounds */
GCG_EXPORT
int GCGpricingprobGetNColsLastRounds(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   int                   nroundscol          /**< number of previous pricing rounds for which the number of improving columns should be counted */
   );

/**@} */

#ifdef __cplusplus
}
#endif
#endif
