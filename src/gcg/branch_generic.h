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
