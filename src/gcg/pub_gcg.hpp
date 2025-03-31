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
