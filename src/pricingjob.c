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

/**@file   pricingjob.c
 * @brief  methods for working with pricing jobs
 * @author Christian Puchert
 *
 * Various methods to work with pricing jobs
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "pub_pricingjob.h"

#include "gcg.h"

#include <assert.h>

/** create a pricing job */
SCIP_RETCODE GCGcreatePricingjob(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGJOB**      pricingjob,         /**< pricing job to be created */
   SCIP*                 pricingscip,        /**< SCIP data structure of the corresponding pricing problem */
   int                   probnr              /**< index of the corresponding pricing problem */
)
{
   SCIP_CALL( SCIPallocMemory(scip, pricingjob) );

   (*pricingjob)->pricingscip = pricingscip;
   (*pricingjob)->probnr = probnr;
   (*pricingjob)->heuristic = FALSE;
   (*pricingjob)->cols = NULL;
   (*pricingjob)->ncols = 0;

   return SCIP_OKAY;
}

/** free a pricing job */
void GCGfreePricingjob(
   SCIP*                 scip,               /**< SCIP data structure (master problem) */
   GCG_PRICINGJOB**      pricingjob          /**< pricing job to be freed */
)
{
   SCIPfreeMemory(scip, pricingjob);
   *pricingjob = NULL;
}
