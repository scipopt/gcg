/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2014 Operations Research, RWTH Aachen University       */
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
 * @brief  GCG interface methods
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef GCG_H_
#define GCG_H_

#include "scip/scip.h"

#include "type_branchgcg.h"
#include "type_decomp.h"
#include "type_detector.h"
#include "type_solver.h"

#include "pub_gcgvar.h"
#include "pub_decomp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns whether the scip is the original scip instance */
extern
SCIP_Bool GCGisOriginal(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the scip is the master problem scip */
extern
SCIP_Bool GCGisMaster(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** print out GCG statistics */
SCIP_RETCODE GCGprintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
);

/** gets the total memory used after problem creation stage for all pricingproblems */
extern
SCIP_Real GCGgetPricingprobsMemUsed(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** prints out the degeneracy of the problem */
extern
void GCGprintDegeneracy(
   SCIP*                 scip,               /**< SCIP data structure */
   double                degeneracy          /**< degeneracy to print*/
   );

/** returns the average degeneracy */
extern
SCIP_Real GCGgetDegeneracy(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transforms given values of the given original variables into values of the given master variables */
extern
void GCGtransformOrigvalsToMastervals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            origvars,           /**< array with (subset of the) original variables */
   SCIP_Real*            origvals,           /**< array with values for the given original variables */
   int                   norigvars,          /**< number of given original variables */
   SCIP_VAR**            mastervars,         /**< array of (all present) master variables */
   SCIP_Real*            mastervals,         /**< array to store the values of the master variables */
   int                   nmastervars         /**< number of master variables */
   );

/** transforms given solution of the master problem into solution of the original problem */
extern
SCIP_RETCODE GCGtransformMastersolToOrigsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             mastersol,          /**< solution of the master problem */
   SCIP_SOL**            origsol             /**< pointer to store the new created original problem's solution */
   );

/** returns whether the constraint belongs to GCG or not */
extern
SCIP_Bool GCGisConsGCGCons(
   SCIP_CONS*            cons                /**< constraint to check */
   );

#ifdef __cplusplus
}
#endif
#endif /* GCG_H_ */
