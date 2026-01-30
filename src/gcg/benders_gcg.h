/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   benders_gcg.h
 * @ingroup BENDERS-GCG
 * @brief  GCG Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERS_GCG_H__
#define __SCIP_BENDERS_GCG_H__


#include "gcg/gcg.h"
#include "scip/scip.h"
#include "scip/bendersdefcuts.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the GCG Benders' decomposition and includes it in SCIP
 *
 *  @ingroup BendersIncludes
 */
GCG_EXPORT
SCIP_RETCODE GCGincludeBendersGcg(
   GCG*                  gcg                 /**< GCG data structure */
   );

/**@addtogroup BENDERS-GCG
 *
 * @{
 */

/** returns the last relaxation solution */
GCG_EXPORT
SCIP_SOL* GCGbendersGetRelaxSol(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   );

/** returns the GCG data structure */
GCG_EXPORT
GCG* GCGbendersGetGcg(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
