/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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
 * @brief  struct to store columns (solutions from a pricing problem)
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_GCGCOL_H_
#define GCG_STRUCT_GCGCOL_H_

#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/scip.h"

#include "type_gcgcol.h"

#ifdef __cplusplus
extern "C" {
#endif

struct GCG_Col
{
   SCIP*                pricingprob;        /**< SCIP data structure (pricing problem)*/
   int                  probnr;             /**< number of corresponding pricing problem */
   SCIP_VAR**           vars;               /**< (sorted) array of variables of corresponding pricing problem */
   SCIP_Real*           vals;               /**< array of solution values (belonging to vars) */
   int                  nvars;              /**< number of variables */
   SCIP_Bool            isray;              /**< is the column a ray? */
   SCIP_Real            redcost;            /**< last known reduced cost */
   int                  age;                /**< age of column (number of iterations since it was created;
                                                 each time reduced cost are calculated counts as an interation) */
};

#ifdef __cplusplus
}
#endif

#endif /* STRUCT_GCGCOL_H_ */
