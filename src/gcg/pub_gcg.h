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

/**@file   pub_gcg.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for working with gcg structure
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PUB_GCG_H__
#define GCG_PUB_GCG_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_retcode.h"
#include "gcg/def.h"
#include "gcg/type_gcg.h"


#ifdef NDEBUG
#include "gcg/struct_gcg.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * GCG structure
 */

/**@defgroup GCG GCG Struct
 * @ingroup DATASTRUCTURES
 * @{
 */

/** create a gcg struct */
GCG_EXPORT
SCIP_RETCODE GCGcreate(
   GCG**                gcg                 /**< pointer to store GCG data structure */
   );

/** free a gcg column */
GCG_EXPORT
SCIP_RETCODE GCGfree(
   GCG**                gcgcol              /**< pointer to gcg structure */
   );

/** checks whether the scip is the original scip instance
 * @returns whether the scip is the original scip instance */
GCG_EXPORT
SCIP_Bool GCGisOriginal(
   SCIP*                scip                /**< SCIP data structure */
   );

/** checks whether the scip is the master problem scip
 * @returns whether the scip is the master problem scip */
GCG_EXPORT
SCIP_Bool GCGisMaster(
   SCIP*                scip                /**< SCIP data structure */
   );

/** print out GCG statistics
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGprintStatistics(
   GCG*                 gcg,                /**< GCG data structure */
   FILE*                file                /**< output file or NULL for standard output */
   );

/** print out complete detection statistics
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGprintCompleteDetectionStatistics(
   GCG*                 gcg,                /**< GCG data structure */
   FILE*                file                /**< output file or NULL for standard output */
   );

/** print name of current instance to given output
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGprintInstanceName(
   GCG*                 gcg,                /**< GCG data structure */
   FILE*                file                /**< output file or NULL for standard output */
   );

GCG_EXPORT
SCIP_RETCODE GCGprintBlockcandidateInformation(
   GCG*                 gcg,                /**< GCG data structure */
   FILE*                file                /**< output file or NULL for standard output */
   );

GCG_EXPORT
SCIP_RETCODE GCGprintCompleteDetectionTime(
   GCG*                 gcg,                /**< GCG data structure */
   FILE*                file                /**< output file or NULL for standard output */
   );

GCG_EXPORT
SCIP_RETCODE GCGprintPartitionInformation(
   GCG*                 gcg,                /**< GCG data structure */
   FILE*                file                /**< output file or NULL for standard output */
   );

GCG_EXPORT
SCIP_RETCODE GCGprintDecompInformation(
   GCG*                 gcg,                /**< GCG data structure */
   FILE*                file                /**< output file or NULL for standard output */
   );

/** gets the total memory used after problem creation stage for all pricingproblems */
GCG_EXPORT
SCIP_Real GCGgetPricingprobsMemUsed(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** returns the average degeneracy */
GCG_EXPORT
SCIP_Real GCGgetDegeneracy(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** transforms given values of the given original variables into values of the given master variables
 * @returns the sum of the values of the corresponding master variables that are fixed */
GCG_EXPORT
SCIP_Real GCGtransformOrigvalsToMastervals(
   GCG*                 gcg,                /**< GCG data structure */
   SCIP_VAR**           origvars,           /**< array with (subset of the) original variables */
   SCIP_Real*           origvals,           /**< array with values for the given original variables */
   int                  norigvars,          /**< number of given original variables */
   SCIP_VAR**           mastervars,         /**< array of (all present) master variables */
   SCIP_Real*           mastervals,         /**< array to store the values of the master variables */
   int                  nmastervars         /**< number of master variables */
   );

/** transforms given solution of the master problem into solution of the original problem
 *  @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGtransformMastersolToOrigsol(
   GCG*                 gcg,                /**< GCG data structure */
   SCIP_SOL*            mastersol,          /**< solution of the master problem */
   SCIP_SOL**           origsol,            /**< pointer to store the new created original problem's solution */
   SCIP_Bool            ignorelocalvarbnds, /**< check global or local varbounds */
   SCIP_Bool*           violatesvarbnds     /**< pointer to variable to store whether the solution violates orig varbnds (can be NULL) */
   );

/** Checks whether the constraint belongs to GCG or not
 *  @returns whether the constraint belongs to GCG or not */
GCG_EXPORT
SCIP_Bool GCGisConsGCGCons(
   SCIP_CONS*           cons                /**< constraint to check */
   );


/** returns the pricing problem of the given number */
GCG_EXPORT
SCIP* GCGgetPricingprob(
   GCG*                 gcg,                /**< GCG data structure */
   int                  pricingprobnr       /**< number of the pricing problem */
   );

/** returns the number of pricing problems */
GCG_EXPORT
int GCGgetNPricingprobs(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** returns TRUE iff the pricingproblem of the given number is relevant, that means is not identical to
 *  another and represented by it */
GCG_EXPORT
SCIP_Bool GCGisPricingprobRelevant(
   GCG*                 gcg,                /**< GCG data structure */
   int                  pricingprobnr       /**< number of the pricing problem */
   );

/**
 *  for a given block, return the block by which it is represented
 */
GCG_EXPORT
int GCGgetBlockRepresentative(
   GCG*                 gcg,                /**< GCG data structure */
   int                  pricingprobnr       /**< number of the pricing problem */
   );

/** returns the number of relevant pricing problems */
GCG_EXPORT
int GCGgetNRelPricingprobs(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** returns the number of blocks in the original formulation, that are represented by
 *  the pricingprob with the given number */
GCG_EXPORT
int GCGgetNIdenticalBlocks(
   GCG*                 gcg,                /**< GCG data structure */
   int                  pricingprobnr       /**< number of the pricing problem */
   );

/** returns the number of constraints in the master problem */
GCG_EXPORT
int GCGgetNMasterConss(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** returns the contraints in the master problem */
GCG_EXPORT
SCIP_CONS** GCGgetMasterConss(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** returns the contraints in the original problem that correspond to the constraints in the master problem */
GCG_EXPORT
SCIP_CONS** GCGgetOrigMasterConss(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** returns the convexity constraint for the given block */
GCG_EXPORT
SCIP_CONS* GCGgetConvCons(
   GCG*                 gcg,                /**< GCG data structure */
   int                  blocknr             /**< the number of the block for which we
                                              *   need the convexity constraint */
   );

/** returns whether the master problem is a set covering problem */
GCG_EXPORT
SCIP_Bool GCGisMasterSetCovering(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** returns whether the master problem is a set partitioning problem */
GCG_EXPORT
SCIP_Bool GCGisMasterSetPartitioning(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** returns whether the relaxator has been initialized */
GCG_EXPORT
SCIP_Bool GCGrelaxIsInitialized(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** return linking constraints for variables */
GCG_EXPORT
SCIP_CONS** GCGgetVarLinkingconss(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** return blocks of linking constraints for variables */
GCG_EXPORT
int* GCGgetVarLinkingconssBlock(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** return number of linking constraints for variables */
GCG_EXPORT
int GCGgetNVarLinkingconss(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** return number of linking variables */
GCG_EXPORT
int GCGgetNLinkingvars(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** return number of variables directly transferred to the master problem */
GCG_EXPORT
int GCGgetNTransvars(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** returns the relaxation solution from the Benders' decomposition */
GCG_EXPORT
SCIP_SOL* GCGgetBendersRelaxationSol(
   GCG*                 gcg                 /**< GCG data structure */
   );

#ifdef NDEBUG
#define GCGgetOrigprob(gcg)               (gcg->origprob)
#else
/** returns the original problem */
GCG_EXPORT
SCIP* GCGgetOrigprob(
   GCG*                 gcg                  /**< GCG data structure */
   );
#endif

#ifdef NDEBUG
#define GCGgetMasterprob(gcg)               (gcg->masterprob)
#else
/** returns the active master problem (may change until solving is initiated) */
GCG_EXPORT
SCIP* GCGgetMasterprob(
   GCG*                 gcg                  /**< GCG data structure */
   );
#endif

#ifdef NDEBUG
#define GCGgetBendersMasterprob(gcg)               (gcg->bendersmasterprob)
#else
/** returns the benders master problem (also used to solve the original problem directly) */
GCG_EXPORT
SCIP* GCGgetBendersMasterprob(
   GCG*                 gcg                  /**< GCG data structure */
   );
#endif

#ifdef NDEBUG
#define GCGgetDwMasterprob(gcg)               (gcg->dwmasterprob)
#else
/** returns the dw master problem */
GCG_EXPORT
SCIP* GCGgetDwMasterprob(
   GCG*                 gcg                  /**< GCG data structure */
   );
#endif

/** returns the GCG data structure */
GCG_EXPORT
GCG* GCGmasterGetGcg(
   SCIP*                masterprob           /**< SCIP data structure */
   );

/** returns the GCG data structure */
GCG_EXPORT
GCG* GCGorigGetGcg(
   SCIP*                origprob             /**< SCIP data structure */
   );

/** returns the decomposition mode */
GCG_EXPORT
GCG_DECMODE GCGgetDecompositionMode(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** gets the structure information */
GCG_EXPORT
GCG_DECOMP* GCGgetStructDecomp(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** returns the visualization parameters */
GCG_EXPORT
GCG_PARAMDATA* GCGgetParamsVisu(
   GCG*                  gcg                /**< GCG data structure */
   );


/** return root node clock */
GCG_EXPORT
SCIP_CLOCK* GCGgetRootNodeTime(
   GCG*                  gcg                 /**< GCG data structure */
   );


/** initializes solving data structures and transforms problem for solving with GCG
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post When calling this method in the \ref SCIP_STAGE_PROBLEM stage, the \SCIP stage is changed to \ref
 *        SCIP_STAGE_TRANSFORMED; otherwise, the stage is not changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
GCG_EXPORT
SCIP_RETCODE GCGtransformProb(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** transforms and presolves the problem suitable for GCG
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *
 *  @post After calling this method \SCIP reaches one of the following stages:
 *        - \ref SCIP_STAGE_PRESOLVING if the presolving process was interrupted
 *        - \ref SCIP_STAGE_PRESOLVED if the presolving process was finished and did not solve the problem
 *        - \ref SCIP_STAGE_SOLVED if the problem was solved during presolving
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
GCG_EXPORT
SCIP_RETCODE GCGpresolve(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** transforms and detects the problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *
 *  @post When calling this method in the \ref SCIP_STAGE_PROBLEM stage, the \SCIP stage is changed to \ref
 *        SCIP_STAGE_TRANSFORMED; otherwise, the stage is not changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
GCG_EXPORT
SCIP_RETCODE GCGdetect(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** transforms, resolves, detects, and solves the problem using GCG
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 * @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 * @post After calling this method \SCIP reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *        - \ref SCIP_STAGE_PRESOLVING if the solution process was interrupted during presolving
 *        - \ref SCIP_STAGE_SOLVING if the solution process was interrupted during the tree search
 *        - \ref SCIP_STAGE_SOLVED if the solving process was not interrupted
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
GCG_EXPORT
SCIP_RETCODE GCGsolve(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** @} */

/** gets GCG's global dual bound
 *
 *  Computes the global dual bound while considering the original
 *  problem SCIP instance and the master problem SCIP instance.
 *
 *  @return the global dual bound
 */
GCG_EXPORT
SCIP_Real GCGgetDualbound(
   GCG*                 gcg               /**< GCG data structure */
   );

/** gets GCG's global primal bound
 *
 *  Computes the global primal bound while considering the original
 *  problem SCIP instance and the master problem SCIP instance.
 *
 *  @return the global dual bound
 */
GCG_EXPORT
SCIP_Real GCGgetPrimalbound(
   GCG*                 gcg               /**< GCG data structure */
   );

/** gets GCG's global gap
 *
 *  Computes the global gap based on the gloal dual bound and the
 *  global primal bound.
 *
 *  @return the global dual bound
 *
 *  @see GCGgetDualbound()
 *  @see GCGgetPrimalbound()
 */
GCG_EXPORT
SCIP_Real GCGgetGap(
   GCG*                 gcg               /**< GCG data structure */
   );

#ifdef NDEBUG
#define GCGgetRelax(gcg)               (gcg->relax)
#else
/** gets GCG's relaxator */
GCG_EXPORT
SCIP_RELAX* GCGgetRelax(
   GCG*                 gcg               /**< GCG data structure */
   );
#endif

#ifdef NDEBUG
#define GCGgetSepaorig(gcg)               (gcg->sepaorig)
#else
/** gets the orig separator */
GCG_EXPORT
SCIP_SEPA* GCGgetSepaorig(
   GCG*                 gcg               /**< GCG data structure */
   );
#endif

#ifdef NDEBUG
#define GCGgetMastersepacutEventhdlr(gcg)               (gcg->mastersepacuthdlr)
#else
/** gets the master separator cut event handler */
GCG_EXPORT
SCIP_EVENTHDLR* GCGgetMastersepacutEventhdlr(
   GCG*                 gcg               /**< GCG data structure */
   );
#endif


/**@} */

#ifdef __cplusplus
}
#endif
#endif
