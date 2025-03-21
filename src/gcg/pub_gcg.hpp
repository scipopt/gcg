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

/**@file   pub_gcg.hpp
 * @ingroup PUBLICCOREAPI
 * @brief  public C++ methods for working with gcg structure
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PUB_GCG_HPP__
#define GCG_PUB_GCG_HPP__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_retcode.h"
#include "gcg/def.h"
#include "gcg/type_gcg.h"
#include "gcg/objpricer_gcg.h"


#ifdef NDEBUG
#include "gcg/struct_gcg.h"
#endif

/*
 * GCG structure
 */

/**@defgroup GCG GCG Struct
 * @ingroup DATASTRUCTURES
 * @{
 */

#ifdef NDEBUG
#define GCGgetObjPricer(gcg)               (reinterpret_cast<ObjPricerGcg*>(gcg->pricer))
#else
/** gets the GCG pricer */
GCG_EXPORT
ObjPricerGcg* GCGgetObjPricer(
   GCG*                 gcg               /**< GCG data structure */
   );
#endif

/**@} */

#endif
