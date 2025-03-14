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

/**@file   struct_gcg.h
 * @ingroup DATASTRUCTURES
 * @brief  gcg data structure
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_GCG_H_
#define GCG_STRUCT_GCG_H_

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_relax.h"

#include "gcg/type_gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Column data structure */
struct Gcg
{
   SCIP*                origprob;         /**< SCIP data structure of origprob */
   SCIP*                masterprob;       /**< SCIP data structure of masterprob */
   SCIP*                dwmasterprob;     /**< SCIP data structure of DW masterprob */
   SCIP*                bendersmasterprob;/**< SCIP data structure of Benders masterprob */
   SCIP_RELAX*          relax;            /**< GCG's relaxation handler */
};

#ifdef __cplusplus
}
#endif

#endif /* STRUCT_GCG_H_ */
