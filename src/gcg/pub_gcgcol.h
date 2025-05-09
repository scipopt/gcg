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

/**@file   pub_gcgcol.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for working with gcg columns
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PUB_GCGCOL_H__
#define GCG_PUB_GCGCOL_H__

#include "gcg/type_gcgcol.h"

#include "scip/type_scip.h"
#include "scip/type_retcode.h"
#include "scip/type_var.h"
#include "scip/type_cons.h"
#include "scip/type_misc.h"
#include "gcg/def.h"
#include "gcg/type_gcg.h"
#include "gcg/type_extendedmasterconsdata.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * GCG Column
 */

/**@defgroup GCG_COLUMN GCG Column
 * @ingroup DATASTRUCTURES
 * @{
 */

/** create a gcg column */
GCG_EXPORT
SCIP_RETCODE GCGcreateGcgCol(
   GCG*                 gcg,                /**< GCG data structure */
   SCIP*                pricingprob,        /**< SCIP data structure */
   GCG_COL**            gcgcol,             /**< pointer to store gcg column */
   int                  prob,               /**< number of corresponding pricing problem */
   SCIP_VAR**           vars,               /**< (sorted) array of variables of corresponding pricing problem */
   SCIP_Real*           vals,               /**< array of solution values (belonging to vars) */
   int                  nvars,              /**< number of variables */
   SCIP_Bool            isray,              /**< is the column a ray? */
   SCIP_Real            redcost             /**< last known reduced cost */
);

/** free a gcg column */
GCG_EXPORT
SCIP_RETCODE GCGfreeGcgCol(
   GCG_COL**            gcgcol              /**< pointer to store gcg column */
);

/** create a gcg column from a solution to a pricing problem */
GCG_EXPORT
SCIP_RETCODE GCGcreateGcgColFromSol(
   GCG*                 gcg,                /**< GCG data structure */
   SCIP*                pricingprob,        /**< SCIP data structure (pricing problem) */
   SCIP*                subproblem,         /**< SCIP data structure that contains the actual solution (if NULL pricingprob will be used) */
   SCIP_HASHMAP*        varmap,             /**< mapping of pricingprob vars to subproblem vars (can be NULL if subproblem is NULL) */
   GCG_COL**            gcgcol,             /**< pointer to store gcg column */
   int                  prob,               /**< number of corresponding pricing problem */
   SCIP_SOL*            sol,                /**< solution of pricing problem with index prob */
   SCIP_Bool            isray,              /**< is column a ray? */
   SCIP_Real            redcost             /**< last known reduced cost */
);

/** get pricing problem index of gcg column */
GCG_EXPORT
int GCGcolGetProbNr(
   GCG_COL*             gcgcol              /**< gcg column structure */
);

/** get pricing problem of gcg column */
GCG_EXPORT
SCIP* GCGcolGetPricingProb(
   GCG_COL*             gcgcol              /**< gcg column structure */
);

/** get variables of gcg column */
GCG_EXPORT
SCIP_VAR** GCGcolGetVars(
   GCG_COL*             gcgcol              /**< gcg column structure */
);

/** get values of gcg column */
GCG_EXPORT
SCIP_Real* GCGcolGetVals(
   GCG_COL*             gcgcol              /**< gcg column structure */
);

/** get number of variables of gcg column */
GCG_EXPORT
int GCGcolGetNVars(
   GCG_COL*             gcgcol              /**< gcg column structure */
);

/** is gcg column a ray? */
GCG_EXPORT
SCIP_Bool GCGcolIsRay(
   GCG_COL*             gcgcol              /**< gcg column structure */
);

/** get reduced cost of gcg column */
GCG_EXPORT
SCIP_Real GCGcolGetRedcost(
   GCG_COL*             gcgcol              /**< gcg column structure */
);

/** comparison method for sorting gcg columns by non-decreasing reduced cost */
GCG_EXPORT
SCIP_DECL_SORTPTRCOMP(GCGcolCompRedcost);

/** comparison method for sorting gcg columns by non-increasing age */
GCG_EXPORT
SCIP_DECL_SORTPTRCOMP(GCGcolCompAge);

/** comparison method for gcg columns. Returns TRUE iff columns are equal */
GCG_EXPORT
SCIP_Bool GCGcolIsEq(
   GCG_COL*             gcgcol1,            /**< first gcg column structure */
   GCG_COL*             gcgcol2             /**< second gcg column structure */
);

/** update reduced cost of variable and increase age */
GCG_EXPORT
void GCGcolUpdateRedcost(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   SCIP_Real            redcost,            /**< new reduced cost */
   SCIP_Bool            growold             /**< increase age counter? */
   );

/** return solution value of variable in gcg column */
GCG_EXPORT
SCIP_Real GCGcolGetSolVal(
   GCG_COL*             gcgcol,             /**< gcg column */
   SCIP_VAR*            var                 /**< variable */
   );

/** get master coefficients of column */
GCG_EXPORT
SCIP_Real* GCGcolGetMastercoefs(
   GCG_COL*             gcgcol              /**< gcg column structure */
   );

/** get number of master coefficients of column */
GCG_EXPORT
int GCGcolGetNMastercoefs(
   GCG_COL*             gcgcol              /**< gcg column structure */
   );

/** set master coefficients of column */
GCG_EXPORT
SCIP_RETCODE GCGcolSetMastercoefs(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   SCIP_Real*           mastercoefs,        /**< array of master coefficients */
   int                  nmastercoefs        /**< number of master coefficients */
   );

/** set norm of column */
GCG_EXPORT
void GCGcolSetNorm(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   SCIP_Real            norm                /**< norm of column */
   );

/** get norm of column */
GCG_EXPORT
SCIP_RETCODE GCGcolComputeNorm(
   GCG*                 gcg,                /**< GCG data structure */
   GCG_COL*             gcgcol              /**< gcg column structure */
   );

/** set master coefficients of column as initialized */
GCG_EXPORT
SCIP_RETCODE GCGcolSetInitializedCoefs(
   GCG_COL*             gcgcol              /**< gcg column structure */
   );

/** return if master coefficients of column have been initialized */
GCG_EXPORT
SCIP_Bool GCGcolGetInitializedCoefs(
   GCG_COL*             gcgcol              /**< gcg column structure */
   );

/** get master coefficients of column */
GCG_EXPORT
int* GCGcolGetLinkvars(
   GCG_COL*             gcgcol              /**< gcg column structure */
   );

/** get number of master coefficients of column */
GCG_EXPORT
int GCGcolGetNLinkvars(
   GCG_COL*             gcgcol              /**< gcg column structure */
   );

/** set master coefficients information of column */
GCG_EXPORT
SCIP_RETCODE GCGcolSetLinkvars(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   int*                 linkvars,           /**< array of linking variable indices for gcgcol->var */
   int                  nlinkvars           /**< number of linking variables in gcgcol->var */
   );

/** get original separator cut coefficients of column in master problem */
GCG_EXPORT
SCIP_Real* GCGcolGetOriginalSepaMastercuts(
   GCG_COL*             gcgcol              /**< gcg column structure */
   );

/** get number of original separator cut coefficients of column in master problem */
GCG_EXPORT
int GCGcolGetNOriginalSepaMastercuts(
   GCG_COL*             gcgcol              /**< gcg column structure */
   );

/** get norm of column */
GCG_EXPORT
SCIP_Real GCGcolGetNorm(
   GCG_COL*             gcgcol              /**< gcg column structure */
   );

/** update original separator cut coefficients information of column in the amster problem */
GCG_EXPORT
SCIP_RETCODE GCGcolUpdateOriginalSepaMastercuts(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   SCIP_Real*           neworiginalsepamastercuts,/**< pointer to new array of master cut coefficients */
   int                  nneworiginalsepamastercuts/**< new number of master cut coefficients */
   );

/** set extended master cons coefficients information of column in the master problem
 * @note the arrays will be freed by the column, they must be allocated using the pricingscip the column belongs to
 */
GCG_EXPORT
SCIP_RETCODE GCGcolSetExtendedmasterconss(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   SCIP_Real*           extendedmasterconss,  /**< pointer to array of extended master cons coefficients */
   SCIP_Real*           extendedmasterconsbounds,/**< pointer to array of extended master cons bounds */
   int                  nextendedmasterconss  /**< number of extended master cons coefficients */
   );

/** gets the age of the col */
GCG_EXPORT
int GCGcolGetAge(
   GCG_COL*             col                 /**< col */
   );

/** returns whether the col's age exceeds the age limit */
GCG_EXPORT
SCIP_Bool GCGcolIsAged(
   GCG_COL*              col,                /**< col to check */
   int                   agelimit            /**< maximum age a col can reach before it is deleted from the pool, or -1 */
   );

/** compute parallelism of column to dual objective */
GCG_EXPORT
SCIP_RETCODE GCGcolComputeDualObjPara(
   GCG*                 gcg,                /**< GCG data structure */
   GCG_COL*             gcgcol,             /**< gcg column */
   SCIP_Real*           para                /**< pointer to store the parallelism of column to dual objective */
   );

/** compute orthogonality of two gcg columns */
GCG_EXPORT
SCIP_RETCODE GCGcolComputeOrth(
   GCG*                 gcg,                /**< GCG data structure */
   GCG_COL*             gcgcol1,            /**< first gcg column */
   GCG_COL*             gcgcol2,            /**< second gcg column */
   SCIP_Real*           orth                /**< pointer to store the orthogonality of two gcg columns */
   );

/** returns the inferred (coefficient) pricing variables solution values */
GCG_EXPORT
SCIP_Real* GCGcolGetInferredPricingVals(
   GCG_COL*              gcgcol             /**< gcgcol */
   );

/** returns the inferred (coefficient) pricing variables */
GCG_EXPORT
SCIP_VAR** GCGcolGetInferredPricingVars(
   GCG_COL*              gcgcol             /**< gcgcol */
   );

/** returns the number of inferred (coefficient) pricing variables */
int GCGcolGetNInferredPricingVars(
   GCG_COL*              gcgcol             /**< gcgcol */
   );

/** gets the hash key of a col */
GCG_EXPORT
SCIP_DECL_HASHGETKEY(GCGhashGetKeyCol);

/** returns TRUE iff both cols are identical */
GCG_EXPORT
SCIP_DECL_HASHKEYEQ(GCGhashKeyEqCol);

/** calculates the hash key value of a col */
GCG_EXPORT
SCIP_DECL_HASHKEYVAL(GCGhashKeyValCol);

/**@} */

#ifdef __cplusplus
}
#endif
#endif
