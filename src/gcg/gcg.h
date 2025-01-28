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

/**@file   gcg.h
 * @ingroup PUBLICCOREAPI
 * @brief  GCG interface methods
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_STATISTIC */

#ifndef GCG_H_
#define GCG_H_

#include "scip/scip.h"
#include "def.h"

#include "type_branchgcg.h"
#include "type_decomp.h"
#include "type_detector.h"
#include "type_solver.h"
#include "type_score.h"

/* include public interfaces, s.t. the user only needs to include gcg.h */
#include "pub_clscons.h"
#include "pub_clsvar.h"
#include "pub_decomp.h"
#include "pub_gcgvar.h"
#include "pub_score.h"

#include "relax_gcg.h"
#include "gcg_general.h"

#ifdef __cplusplus
extern "C" {
#endif

/** checks whether the scip is the original scip instance
 * @returns whether the scip is the original scip instance */
GCG_EXPORT
SCIP_Bool GCGisOriginal(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** checks whether the scip is the master problem scip
 * @returns whether the scip is the master problem scip */
GCG_EXPORT
SCIP_Bool GCGisMaster(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** print out GCG statistics
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGprintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
);

/** print out complete detection statistics
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGprintCompleteDetectionStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
);

/** print name of current instance to given output
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGprintInstanceName(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
);

GCG_EXPORT
SCIP_RETCODE GCGprintMiplibStructureInformation(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_DIALOGHDLR*      dialoghdlr         /**< dialog handler */
   );

GCG_EXPORT
SCIP_RETCODE GCGprintBlockcandidateInformation(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
);

GCG_EXPORT
SCIP_RETCODE GCGprintCompleteDetectionTime(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
);

GCG_EXPORT
SCIP_RETCODE GCGprintPartitionInformation(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
);

GCG_EXPORT
SCIP_RETCODE GCGprintDecompInformation(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
);

/** gets the total memory used after problem creation stage for all pricingproblems */
GCG_EXPORT
SCIP_Real GCGgetPricingprobsMemUsed(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** prints out the degeneracy of the problem */
GCG_EXPORT
void GCGprintDegeneracy(
   SCIP*                 scip,               /**< SCIP data structure */
   double                degeneracy          /**< degeneracy to print*/
   );

/** returns the average degeneracy */
GCG_EXPORT
SCIP_Real GCGgetDegeneracy(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transforms given values of the given original variables into values of the given master variables
 * @returns the sum of the values of the corresponding master variables that are fixed */
GCG_EXPORT
SCIP_Real GCGtransformOrigvalsToMastervals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            origvars,           /**< array with (subset of the) original variables */
   SCIP_Real*            origvals,           /**< array with values for the given original variables */
   int                   norigvars,          /**< number of given original variables */
   SCIP_VAR**            mastervars,         /**< array of (all present) master variables */
   SCIP_Real*            mastervals,         /**< array to store the values of the master variables */
   int                   nmastervars         /**< number of master variables */
   );

/** transforms given solution of the master problem into solution of the original problem
 *  @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGtransformMastersolToOrigsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             mastersol,          /**< solution of the master problem */
   SCIP_SOL**            origsol             /**< pointer to store the new created original problem's solution */
   );

/** Checks whether the constraint belongs to GCG or not
 *  @returns whether the constraint belongs to GCG or not */
GCG_EXPORT
SCIP_Bool GCGisConsGCGCons(
   SCIP_CONS*            cons                /**< constraint to check */
   );


/** returns the original problem for the given master problem */
GCG_EXPORT
SCIP* GCGgetOriginalprob(
   SCIP*                 masterprob          /**< the SCIP data structure for the master problem */
   );

/** returns the master problem */
GCG_EXPORT
SCIP* GCGgetMasterprob(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the pricing problem of the given number */
GCG_EXPORT
SCIP* GCGgetPricingprob(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   );

/** returns the number of pricing problems */
GCG_EXPORT
int GCGgetNPricingprobs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns TRUE iff the pricingproblem of the given number is relevant, that means is not identical to
 *  another and represented by it */
GCG_EXPORT
SCIP_Bool GCGisPricingprobRelevant(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   );

/**
 *  for a given block, return the block by which it is represented
 */
GCG_EXPORT
int GCGgetBlockRepresentative(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   );

/** returns the number of relevant pricing problems */
GCG_EXPORT
int GCGgetNRelPricingprobs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of blocks in the original formulation, that are represented by
 *  the pricingprob with the given number */
GCG_EXPORT
int GCGgetNIdenticalBlocks(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   );

/** returns the number of constraints in the master problem */
GCG_EXPORT
int GCGgetNMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the contraints in the master problem */
GCG_EXPORT
SCIP_CONS** GCGgetMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the contraints in the original problem that correspond to the constraints in the master problem */
GCG_EXPORT
SCIP_CONS** GCGgetOrigMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the convexity constraint for the given block */
GCG_EXPORT
SCIP_CONS* GCGgetConvCons(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   blocknr             /**< the number of the block for which we
                                              *   need the convexity constraint */
   );

/** returns whether the master problem is a set covering problem */
GCG_EXPORT
SCIP_Bool GCGisMasterSetCovering(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the master problem is a set partitioning problem */
GCG_EXPORT
SCIP_Bool GCGisMasterSetPartitioning(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the relaxator has been initialized */
GCG_EXPORT
SCIP_Bool GCGrelaxIsInitialized(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return linking constraints for variables */
GCG_EXPORT
SCIP_CONS** GCGgetVarLinkingconss(
   SCIP*                 scip                /**< SCIP data structure */
  );

/** return blocks of linking constraints for variables */
GCG_EXPORT
int* GCGgetVarLinkingconssBlock(
   SCIP*                 scip                /**< SCIP data structure */
  );

/** return number of linking constraints for variables */
GCG_EXPORT
int GCGgetNVarLinkingconss(
   SCIP*                 scip                /**< SCIP data structure */
  );

/** return number of linking variables */
GCG_EXPORT
int GCGgetNLinkingvars(
   SCIP*                 scip                /**< SCIP data structure */
  );

/** return number of variables directly transferred to the master problem */
GCG_EXPORT
int GCGgetNTransvars(
   SCIP*                 scip                /**< SCIP data structure */
  );

/** returns the auxiliary variable for the given pricing probblem */
GCG_EXPORT
SCIP_VAR* GCGgetAuxiliaryVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   );

/** returns the relaxation solution from the Benders' decomposition */
GCG_EXPORT
SCIP_SOL* GCGgetBendersRelaxationSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif
#endif /* GCG_H_ */
