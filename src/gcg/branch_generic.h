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

/**@file   branch_generic.h
 * @ingroup BRANCHINGRULES-GCG
 * @brief  branching rule based on vanderbeck's generic branching scheme
 * @author Marcel Schmickerath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_GENERIC_H__
#define __SCIP_BRANCH_GENERIC_H__


#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
   GCG_COMPSENSE_GE = 1,
   GCG_COMPSENSE_LT = 0
} GCG_COMPSENSE;

/** component bound structure */
struct ComponentBoundSequence
{
   SCIP_VAR*             component;          /**< variable to which this bound belongs */
   GCG_COMPSENSE         sense;              /**< sense of the bound */
   SCIP_Real             bound;              /**< bound value */
};
typedef struct ComponentBoundSequence GCG_COMPSEQUENCE;

/** strip structure */
struct GCG_Strip
{
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_VAR*             mastervar;          /**< master variable */
   GCG_COMPSEQUENCE**    C;                  /**< current set of comp bound sequences */
   int                   Csize;              /**< number of component bound sequences */
   int*                  sequencesizes;      /**< array of sizes of component bound sequences */
};
typedef struct GCG_Strip GCG_STRIP;

/** creates the generic branching rule and includes it in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeBranchruleGeneric(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** get component bound sequence */
GCG_COMPSEQUENCE* GCGbranchGenericBranchdataGetConsS(
   GCG_BRANCHDATA*       branchdata          /**< branching data to initialize */
   );

/** get size of component bound sequence */
int GCGbranchGenericBranchdataGetConsSsize(
   GCG_BRANCHDATA*       branchdata          /**< branching data to initialize */
   );

/** get id of pricing problem (or block) to which the constraint belongs */
int GCGbranchGenericBranchdataGetConsblocknr(
   GCG_BRANCHDATA*       branchdata          /**< branching data to initialize */
   );

/** get master constraint */
SCIP_CONS* GCGbranchGenericBranchdataGetMastercons(
   GCG_BRANCHDATA*       branchdata          /**< branching data to initialize */
   );

/** returns true when the branch rule is the generic branchrule */
SCIP_Bool GCGisBranchruleGeneric(
   SCIP_BRANCHRULE*      branchrule          /**< branchrule to check */
);

#ifdef __cplusplus
}
#endif

#endif
