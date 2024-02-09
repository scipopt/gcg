/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file    struct_mastercutdata.h
 * @ingroup DATASTRUCTURES
 * @brief   data structures for GCG mastercut data
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_MASTERCUTDATA_H_
#define GCG_STRUCT_MASTERCUTDATA_H_

#include <scip/type_cons.h>
#include <scip/type_lp.h>
#include <scip/type_var.h>

#ifdef __cplusplus
extern "C" {
#endif

/** data for master variables */
struct GCG_MasterCutData
{
   SCIP_ROW*             mastercons;         /**< row in the master problem that represents the master cut */
   SCIP_VAR**            pricingvars;        /**< array of additional variables in the pricing programs inferred from the master cut */
   int                   npricingvars;       /**< number of additional variables in the pricing programs */
   SCIP_CONS**           npricingconss;      /**< array of additional constraints in the pricing programs inferred from the master cut */
   int                   nnpricingconss;     /**< number of additional constraints in the pricing programs */
   int                   blocknr;            /**< block number of the master cut */
};

#ifdef __cplusplus
}
#endif

#endif
