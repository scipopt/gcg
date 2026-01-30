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

/**@file    branch_compbnd.h
 * @ingroup BRANCHINGRULES-GCG
 * @brief   component bound branching rule
 * @author  Til Mohr
 *
 * This is an implementation of the component bound branching rule based on the papers:
 *
 * J. Desrosiers, M. L¨ubbecke, G. Desaulniers,
 * J. B. Gauthier (Juin 2024). Branch-and-Price, Technical report,
 * Les Cahiers du GERAD G–2024–36, GERAD, HEC Montr´eal, Canada.
 *
 * Vanderbeck, François, and Laurence A. Wolsey. "An exact algorithm for IP column generation."
 * Operations research letters 19.4 (1996): 151-159.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_BRANCH_COMPBND_H__
#define GCG_BRANCH_COMPBND_H__


#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
   GCG_BRANCH_DOWN = 0,
   GCG_BRANCH_UP = 1
} GCG_BRANCH_TYPE;

typedef enum {
   GCG_COMPBND_SENSE_GE = 0,
   GCG_COMPBND_SENSE_LE = 1
} GCG_COMPBND_SENSE;

/** component bound structure */
struct ComponentBound
{
   SCIP_VAR*             component;          /**< variable to which this bound belongs */
   GCG_COMPBND_SENSE     sense;              /**< sense of the bound */
   int                   bound;              /**< bound value */
};
typedef struct ComponentBound GCG_COMPBND;

/** creates the component bound branching rule and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeBranchruleCompBnd(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns true when the branch rule is the generic branchrule */
SCIP_Bool GCGisBranchruleCompBnd(
   SCIP_BRANCHRULE*      branchrule          /**< branchrule to check */
);

#ifdef __cplusplus
}
#endif

#endif
