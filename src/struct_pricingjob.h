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

/**@file   struct_gcgcol.h
 * @brief  struct to store pricing jobs
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_PRICINGJOB_H_
#define GCG_STRUCT_PRICINGJOB_H_

#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/scip.h"

#include "type_pricingjob.h"
#include "type_gcgcol.h"

#ifdef __cplusplus
extern "C" {
#endif

struct GCG_PricingJob
{
   /* problem data */
   SCIP*                pricingscip;        /**< SCIP data structure of the corresponding pricing problem */
   int                  probnr;             /**< index of the corresponding pricing problem */

   /* strategic parameters */
   SCIP_Real            score;              /**< current score of the pricing job */
   SCIP_Bool            heuristic;          /**< shall the pricing problem be solved heuristically? */

   /* result values */
   SCIP_STATUS          pricingstatus;      /**< current solution status of the pricing problem */
   SCIP_Real            lowerbound;         /**< lower bound obtained by solving the pricing problem */
   GCG_COL**            cols;               /**< array of columns found by the pricing problem */
   int                  ncols;              /**< number of columns found by the pricing problem */
};

#ifdef __cplusplus
}
#endif

#endif /* STRUCT_PRICINGJOB_H_ */
