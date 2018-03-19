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

/**@file   struct_pricingprob.h
 * @brief  data structure to store pricing problem information
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_PRICINGPROB_H_
#define GCG_STRUCT_PRICINGPROB_H_

#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/scip.h"

#include "type_pricingprob.h"
#include "type_gcgcol.h"
#include "type_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

struct GCG_PricingProb
{
   /* problem data */
   SCIP*                pricingscip;        /**< SCIP data structure */
   int                  probnr;             /**< (block) index of the corresponding pricing problem */

   /* result values */
   int                  nsolves;            /**< number of times the pricing problem was solved during the loop */
   SCIP_STATUS          pricingstatus;      /**< current solution status of the pricing problem */
   SCIP_Real            lowerbound;         /**< lower bound obtained by solving the pricing problem */
   GCG_COL**            cols;               /**< array of columns found in the current pricing round */
   int                  colssize;           /**< size of column array */
   int                  ncols;              /**< number of columns found in the current pricing round */
   int                  nimpcols;           /**< number of improving columns found in the current pricing round */

   /* statistics */
   int*                 ncolsround;         /**< number of improving columns found in the last rounds */
};

#ifdef __cplusplus
}
#endif

#endif /* STRUCT_PRICINGPROB_H_ */
