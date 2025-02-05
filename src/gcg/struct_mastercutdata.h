/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

#include <scip/def.h>
#include <scip/type_cons.h>
#include <scip/type_lp.h>
#include <scip/type_var.h>
#include "type_mastercutdata.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data for a pricing problem modification */
struct GCG_PricingModification
{
   int                   blocknr;            /**< block number of the master cut */
   SCIP_VAR*             coefvar;            /**< variable in the pricing problem inferred from the master cut
                                                * always has the objective coefficient of the negated dual value of the master cut
                                                * its solution value corresponds to the coefficient of the new mastervariable in the master cut */
   SCIP_VAR**            additionalvars;     /**< array of additional variables with no objective coefficient in the pricing programs inferred from the master cut */
   int                   nadditionalvars;    /**< number of additional variables in the pricing programs */
   SCIP_CONS**           additionalconss;    /**< array of additional constraints in the pricing programs inferred from the master cut */
   int                   nadditionalconss;   /**< number of additional constraints in the pricing programs */
};

/** type of master cut */
enum GCG_MasterCutType
{
   GCG_MASTERCUTTYPE_CONS,                   /**< master cut is represented by a constraint */
   GCG_MASTERCUTTYPE_ROW                     /**< master cut is represented by a row */
};
typedef enum GCG_MasterCutType GCG_MASTERCUTTYPE;

/** cut of the master cut */
union GCG_MasterCutCut
{
   SCIP_CONS*           cons;                /**< constraint in the master problem that represents the master cut, iff type == Cons */
   SCIP_ROW*            row;                 /**< row in the master problem that represents the master cut, iff type == Row */
};
typedef union GCG_MasterCutCut GCG_MASTERCUTCUT;

/** data for master cuts */
struct GCG_MasterCutData
{
   GCG_MASTERCUTTYPE     type;               /**< type of the master cut */
   GCG_MASTERCUTCUT      cut;                /**< constraint or row in the master problem that represents the master cut */
   GCG_PRICINGMODIFICATION* pricingmodifications; /**< array of pricing modifications for the master cut */
   int                   npricingmodifications; /**< number of pricing modifications for the master cut */
   void*                 data;               /**< any data that might be required to calculate the coefficient of a column solution */
   GCG_DECL_MASTERCUTGETCOEFF ((*mastercutGetCoeff)); /**< callback to calculate the coefficient of a column solution */
};

#ifdef __cplusplus
}
#endif

#endif
