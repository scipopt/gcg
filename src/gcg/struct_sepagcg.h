/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       */
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

/**@file    struct_sepagcg.h
 * @ingroup DATASTRUCTURES
 * @brief   data structures for separators for master
 * @author  Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_SEPAGCG_H__
#define GCG_STRUCT_SEPAGCG_H__

#include <scip/type_sepa.h>

#include "type_sepagcg.h"
#include "event_sepacuts.h"

#ifdef __cplusplus
extern "C" {
#endif

/** master separator */
struct GCG_Sepa
{
   SCIP_SEPA*                             separator;                     /**< SCIP separator */
   GCG_DECL_SEPAGETCOLCOEFFICIENTS        ((*gcgsepagetcolcoefficient)); /**< compute cut coefficient using gcg column */
   GCG_DECL_SEPAGETVARCOEFFICIENT         ((*gcgsepagetvarcoefficient)); /**< compute cut coefficient using variable values */
   GCG_DECL_SEPASETOBJECTIVE              ((*gcgsepasetobjective));      /**< adapt pricing objectives to consider cut */
   GCG_DECL_SEPAADJUSTCOL                 ((*gcgsepaadjustcol));         /**< modify outdated column to consider cut */
};

#ifdef __cplusplus
}
#endif

#endif