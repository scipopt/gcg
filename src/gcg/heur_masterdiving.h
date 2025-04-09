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

/**@file   heur_masterdiving.h
 * @ingroup DIVINGHEURISTICS
 * @brief  primal heuristic interface for LP diving heuristics on the master variables
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_HEUR_MASTERDIVING_H__
#define GCG_HEUR_MASTERDIVING_H__



#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** gets diving rule specific data of a diving heuristic */
GCG_EXPORT
GCG_MASTER_DIVINGDATA* GCGheurGetDivingDataMaster(
   SCIP_HEUR*               heur                    /**< primal heuristic */
   );

/** sets diving rule specific data of a diving heuristic */
GCG_EXPORT
void GCGheurSetDivingDataMaster(
   SCIP_HEUR*               heur,                   /**< primal heuristic */
   GCG_MASTER_DIVINGDATA*   divingdata              /**< diving rule specific data */
   );

/** creates a master diving heuristic and includes it in GCG */
GCG_EXPORT
SCIP_RETCODE GCGincludeDivingHeurMaster(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_HEUR**           heur,               /**< pointer to diving heuristic */
   const char*           name,               /**< name of primal heuristic */
   const char*           desc,               /**< description of primal heuristic */
   char                  dispchar,           /**< display character of primal heuristic */
   int                   priority,           /**< priority of the primal heuristic */
   int                   freq,               /**< frequency for calling primal heuristic */
   int                   freqofs,            /**< frequency offset for calling primal heuristic */
   int                   maxdepth,           /**< maximal depth level to call heuristic at (-1: no limit) */
   GCG_DECL_MASTER_DIVINGFREE   ((*divingfree)),         /**< destructor of diving heuristic */
   GCG_DECL_MASTER_DIVINGINIT   ((*divinginit)),         /**< initialize diving heuristic */
   GCG_DECL_MASTER_DIVINGEXIT   ((*divingexit)),         /**< deinitialize diving heuristic */
   GCG_DECL_MASTER_DIVINGINITSOL ((*divinginitsol)),     /**< solving process initialization method of diving heuristic */
   GCG_DECL_MASTER_DIVINGEXITSOL ((*divingexitsol)),     /**< solving process deinitialization method of diving heuristic */
   GCG_DECL_MASTER_DIVINGINITEXEC ((*divinginitexec)),   /**< execution initialization method of diving heuristic */
   GCG_DECL_MASTER_DIVINGEXITEXEC ((*divingexitexec)),   /**< execution deinitialization method of diving heuristic */
   GCG_DECL_MASTER_DIVINGSELECTVAR ((*divingselectvar)), /**< variable selection method of diving heuristic */
   GCG_MASTER_DIVINGDATA*       divingdata          /**< diving rule specific data (or NULL) */
   );

/** creates event handler for masterdiving event */
GCG_EXPORT
SCIP_RETCODE GCGincludeEventHdlrMasterdiving(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
