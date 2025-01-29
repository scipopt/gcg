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

/**@file   pub_decomp.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for working with decomposition structures
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PUB_DECOMP_H__
#define GCG_PUB_DECOMP_H__

#include "type_decomp.h"
#include "scip/type_scip.h"
#include "scip/type_retcode.h"
#include "scip/type_var.h"
#include "scip/type_cons.h"
#include "scip/type_misc.h"
#include "def.h"
#include "type_detector.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @ingroup DECOMP
 * @{
 */

/** score data structure
 * @ingroup DATASTRUCTURES
**/
struct Dec_Scores
{
   SCIP_Real             borderscore;        /**< score of the border */
   SCIP_Real             densityscore;       /**< score of block densities */
   SCIP_Real             linkingscore;       /**< score related to interlinking blocks */
   SCIP_Real             totalscore;         /**< accumulated score */
   SCIP_Real             maxwhitescore;      /**< score related to max white measure (i.e. fraction of white (nonblock and nonborder) matrix area ) */
};
typedef struct Dec_Scores DEC_SCORES;

/** converts the GCG_DECTYPE enum to a string */
const char *GCGdecompGetStrType(
   GCG_DECTYPE           type                /**< decomposition type */
   );

/** initializes the decomposition to absolutely nothing */
GCG_EXPORT
SCIP_RETCODE GCGdecompCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP**          decomp              /**< pointer to the decomposition data structure */
   );

/** frees the decomposition */
GCG_EXPORT
SCIP_RETCODE GCGdecompFree(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP**          decomp              /**< pointer to the decomposition data structure */
   );

/** sets the type of the decomposition */
GCG_EXPORT
SCIP_RETCODE GCGdecompSetType(
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   GCG_DECTYPE           type                /**< type of the decomposition */
   );

/** gets the type of the decomposition */
GCG_EXPORT
GCG_DECTYPE GCGdecompGetType(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

GCG_EXPORT
SCIP_Real GCGdecompGetMaxwhiteScore(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

GCG_EXPORT
void GCGdecompSetMaxWhiteScore(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_Real             maxwhitescore       /**< score related to max white measure (i.e. fraction of white (nonblock and nonborder) matrix area ) */
   );

/** sets the presolved flag for decomposition */
GCG_EXPORT
void GCGdecompSetPresolved(
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_Bool             presolved           /**< presolved flag for decomposition */
   );

/** gets the presolved flag for decomposition */
GCG_EXPORT
SCIP_Bool GCGdecompGetPresolved(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** sets the number of blocks for decomposition */
GCG_EXPORT
void GCGdecompSetNBlocks(
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   int                   nblocks             /**< number of blocks for decomposition */
   );

/** gets the number of blocks for decomposition */
GCG_EXPORT
int GCGdecompGetNBlocks(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** copies the input subscipvars array to the given decomposition */
GCG_EXPORT
SCIP_RETCODE GCGdecompSetSubscipvars(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_VAR***           subscipvars,        /**< subscipvars array  */
   int*                  nsubscipvars        /**< number of subscipvars per block */
   );

/** returns the subscipvars array of the given decomposition */
GCG_EXPORT
SCIP_VAR*** GCGdecompGetSubscipvars(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** returns the nsubscipvars array of the given decomposition */
GCG_EXPORT
int* GCGdecompGetNSubscipvars(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** copies the input subscipconss array to the given decomposition */
GCG_EXPORT
SCIP_RETCODE GCGdecompSetSubscipconss(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_CONS***          subscipconss,       /**< subscipconss array  */
   int*                  nsubscipconss       /**< number of subscipconss per block */
   );

/** returns the subscipconss array of the given decomposition */
GCG_EXPORT
SCIP_CONS*** GCGdecompGetSubscipconss(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** returns the nsubscipconss array of the given decomposition */
GCG_EXPORT
int* GCGdecompGetNSubscipconss(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** copies the input linkingconss array to the given decomposition */
GCG_EXPORT
SCIP_RETCODE GCGdecompSetLinkingconss(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_CONS**           linkingconss,       /**< linkingconss array  */
   int                   nlinkingconss       /**< number of linkingconss per block */
   );

/** returns the linkingconss array of the given decomposition */
GCG_EXPORT
SCIP_CONS** GCGdecompGetLinkingconss(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** returns the nlinkingconss array of the given decomposition */
GCG_EXPORT
int GCGdecompGetNLinkingconss(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** copies the input linkingvars array to the given decomposition */
GCG_EXPORT
SCIP_RETCODE GCGdecompSetLinkingvars(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_VAR**            linkingvars,        /**< linkingvars array  */
   int                   nlinkingvars,       /**< number of linkingvars per block */
   int                   nfixedlinkingvars,  /**< number of linking variables that are fixed */
   int                   nmastervars         /**< number of linkingvars that are purely master variables */
   );

/** returns the linkingvars array of the given decomposition */
GCG_EXPORT
SCIP_VAR** GCGdecompGetLinkingvars(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** returns the number of master variables of the given decomposition */
GCG_EXPORT
int GCGdecompGetNMastervars(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );


/** returns the nlinkingvars array of the given decomposition */
GCG_EXPORT
int GCGdecompGetNLinkingvars(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** returns the nlinkingvars array of the given decomposition */
GCG_EXPORT
int GCGdecompGetNFixedLinkingvars(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );


/** copies the input stairlinkingvars array to the given decomposition */
GCG_EXPORT
SCIP_RETCODE GCGdecompSetStairlinkingvars(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_VAR***           stairlinkingvars,   /**< stairlinkingvars array  */
   int*                  nstairlinkingvars   /**< number of linkingvars per block */
   );

/** returns the stairlinkingvars array of the given decomposition */
GCG_EXPORT
SCIP_VAR*** GCGdecompGetStairlinkingvars(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** returns the nstairlinkingvars array of the given decomposition */
GCG_EXPORT
int* GCGdecompGetNStairlinkingvars(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** returns the total number of stairlinkingvars array of the given decomposition */
GCG_EXPORT
int GCGdecompGetNTotalStairlinkingvars(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );


/** sets the vartoblock hashmap of the given decomposition */
GCG_EXPORT
void GCGdecompSetVartoblock(
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_HASHMAP*         vartoblock          /**< Vartoblock hashmap */
   );

/** returns the vartoblock hashmap of the given decomposition */
GCG_EXPORT
SCIP_HASHMAP* GCGdecompGetVartoblock(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** sets the constoblock hashmap of the given decomposition */
GCG_EXPORT
void GCGdecompSetConstoblock(
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_HASHMAP*         constoblock         /**< Constoblock hashmap */
   );

/** returns the constoblock hashmap of the given decomposition */
GCG_EXPORT
SCIP_HASHMAP* GCGdecompGetConstoblock(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** sets the varindex hashmap of the given decomposition */
GCG_EXPORT
void GCGdecompSetVarindex(
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_HASHMAP*         varindex            /**< Varindex hashmap */
   );

/** returns the varindex hashmap of the given decomposition */
GCG_EXPORT
SCIP_HASHMAP* GCGdecompGetVarindex(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** sets the consindex hashmap of the given decomposition */
GCG_EXPORT
void GCGdecompSetConsindex(
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_HASHMAP*         consindex           /**< Consindex hashmap */
   );

/** returns the consindex hashmap of the given decomposition */
GCG_EXPORT
SCIP_HASHMAP* GCGdecompGetConsindex(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** completely initializes decomposition structure from the values of the hashmaps */
GCG_EXPORT
SCIP_RETCODE GCGfilloutDecompFromHashmaps(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_HASHMAP*         vartoblock,         /**< variable to block hashmap */
   SCIP_HASHMAP*         constoblock,        /**< constraint to block hashmap */
   int                   nblocks,            /**< number of blocks */
   SCIP_Bool             staircase           /**< should the decomposition be a staircase structure */
   );

/** completely fills out decomposition structure from only the constraint partition */
GCG_EXPORT
SCIP_RETCODE GCGfilloutDecompFromConstoblock(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_HASHMAP*         constoblock,        /**< constraint to block hashmap, start with 1 for first block and nblocks+1 for linking constraints */
   int                   nblocks,            /**< number of blocks */
   SCIP_Bool             staircase           /**< should the decomposition be a staircase structure */
   );

/** sets the detector for the given decomposition */
GCG_EXPORT
void GCGdecompSetDetector(
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   GCG_DETECTOR*         detector            /**< detector data structure */
   );

/** gets the detector for the given decomposition */
GCG_EXPORT
GCG_DETECTOR* GCGdecompGetDetector(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** gets the detectors for the given decomposition */
GCG_EXPORT
GCG_DETECTOR** GCGdecompGetDetectorChain(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** gets the number of detectors for the given decomposition */
GCG_EXPORT
int GCGdecompGetDetectorChainSize(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** gets the id of the original partialdec */
GCG_EXPORT
int GCGdecompGetPartialdecID(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** sets the detector clock times of the detectors of the detector chain */
GCG_EXPORT
void GCGdecompSetDetectorClockTimes(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_Real*            detectorClockTimes  /**< time used by the detectors */
   );

/** gets the detector clock times of the detectors of the detector chain */
GCG_EXPORT
SCIP_Real* GCGdecompGetDetectorClockTimes(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** sets the detector clock times of the detectors of the detector chain */
GCG_EXPORT
SCIP_RETCODE GCGdecompSetDetectorChainString(
   SCIP*                 scip,                /**< SCIP data structure */
   GCG_DECOMP*           decomp,              /**< decomposition data structure */
   const char*           detectorchainstring  /**< string for the detector information working on that decomposition */
   );


GCG_EXPORT
char* GCGdecompGetDetectorChainString(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );


/** sets the percentages of variables assigned to the border of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
void GCGdecompSetDetectorPctVarsToBorder(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,              /**< decomposition data structure */
   SCIP_Real*            pctVarsToBorder
   );

/** gets the percentages of variables assigned to the border of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
SCIP_Real* GCGdecompGetDetectorPctVarsToBorder(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** sets the percentages of constraints assigned to the border of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
void GCGdecompSetDetectorPctConssToBorder(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,              /**< decomposition data structure */
   SCIP_Real*            pctConssToBorder
   );

/** gets the percentages of constraints assigned to the border of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
SCIP_Real* GCGdecompGetDetectorPctConssToBorder(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** sets the percentages of variables assigned to some block of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
void GCGdecompSetDetectorPctVarsToBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,              /**< decomposition data structure */
   SCIP_Real*            pctVarsToBlock
   );

/** gets the percentages of variables assigned to some block of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
SCIP_Real* GCGdecompGetDetectorPctVarsToBlock(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** sets the percentages of constraints assigned to some block of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
void GCGdecompSetDetectorPctConssToBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,              /**< decomposition data structure */
   SCIP_Real*            pctConssToBlock
   );

/** gets the percentages of constraints assigned to some block of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
SCIP_Real* GCGdecompGetDetectorPctConssToBlock(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );


/** sets the percentages of variables assigned to some block of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
void GCGdecompSetDetectorPctVarsFromOpen(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,              /**< decomposition data structure */
   SCIP_Real*            pctVarsFromOpen
   );

/** gets the percentages of variables assigned to some block of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
SCIP_Real* GCGdecompGetDetectorPctVarsFromOpen(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** sets the percentages of constraints assigned to some block of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
void GCGdecompSetDetectorPctConssFromOpen(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,              /**< decomposition data structure */
   SCIP_Real*            pctConssToBorder
   );

/** gets the percentages of constraints assigned to some block of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
SCIP_Real* GCGdecompGetDetectorPctConssFromOpen(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** sets the number of new blocks of the corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
void GCGdecompSetNNewBlocks(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   int*                  nNewBlocks          /**< number of new blocks on this decomposition */
   );

/** gets the number of new blocks corresponding detectors (of the detector chain) on this decomposition */
GCG_EXPORT
int* GCGdecompGetNNewBlocks(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** transforms all constraints and variables, updating the arrays */
GCG_EXPORT
SCIP_RETCODE GCGdecompTransform(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/**
 * Remove all those constraints that were removed from the problem after the decomposition had been created
 */
GCG_EXPORT
SCIP_RETCODE GCGdecompRemoveDeletedConss(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/**
 *  Adds all those constraints that were added to the problem after the decomposition as created
 */
GCG_EXPORT
SCIP_RETCODE GCGdecompAddRemainingConss(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** checks the validity of the decomposition data structure */
GCG_EXPORT
SCIP_RETCODE GCGdecompCheckConsistency(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** creates a decomposition with all constraints in the master */
GCG_EXPORT
SCIP_RETCODE GCGcreateBasicDecomp(
   SCIP*                 scip,                /**< SCIP data structure */
   GCG_DECOMP**          decomp,              /**< decomposition data structure */
   SCIP_Bool             solveorigprob        /**< is the original problem being solved? */
   );

/** creates a decomposition with provided constraints in the master
 * The function will put the remaining constraints in one or more pricing problems
 * depending on whether the subproblems decompose with no variables in common.
 */
GCG_EXPORT
SCIP_RETCODE GCGcreateDecompFromMasterconss(
   SCIP*                 scip,                /**< SCIP data structure */
   GCG_DECOMP**          decomp,              /**< decomposition data structure */
   SCIP_CONS**           conss,               /**< constraints to be put in the master */
   int                   nconss               /**< number of constraints in the master */
   );

/** return the number of variables and binary, integer, implied integer, continuous variables of all subproblems */
GCG_EXPORT
void GCGgetSubproblemVarsData(
   SCIP*                 scip,                /**< SCIP data structure */
   GCG_DECOMP*           decomp,              /**< decomposition data structure */
   int*                  nvars,               /**< pointer to array of size nproblems to store number of subproblem vars or NULL */
   int*                  nbinvars,            /**< pointer to array of size nproblems to store number of binary subproblem vars or NULL */
   int*                  nintvars,            /**< pointer to array of size nproblems to store number of integer subproblem vars or NULL */
   int*                  nimplvars,           /**< pointer to array of size nproblems to store number of implied subproblem vars or NULL */
   int*                  ncontvars,           /**< pointer to array of size nproblems to store number of continuous subproblem vars or NULL */
   int                   nproblems            /**< size of the arrays*/
   );

/** return the number of variables and binary, integer, implied integer, continuous variables of the master */
GCG_EXPORT
void GCGgetLinkingVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   int*                  nvars,              /**< pointer to store number of linking vars or NULL */
   int*                  nbinvars,           /**< pointer to store number of binary linking vars or NULL */
   int*                  nintvars,           /**< pointer to store number of integer linking vars or NULL */
   int*                  nimplvars,          /**< pointer to store number of implied linking vars or NULL */
   int*                  ncontvars           /**< pointer to store number of continuous linking vars or NULL */
   );

/**
 * returns the number of nonzeros of each column of the constraint matrix both in the subproblem and in the master
 * @note For linking variables, the number of nonzeros in the subproblems corresponds to the number on nonzeros
 * in the border
 *
 * @note The arrays have to be allocated by the caller
 *
 * @pre This function assumes that constraints are partitioned in the decomp structure, no constraint is present in more than one block
 *
 */
GCG_EXPORT
SCIP_RETCODE GCGgetDensityData(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< pointer to array store variables belonging to density */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< pointer to array to store constraints belonging to the density */
   int                   nconss,             /**< number of constraints */
   int*                  varsubproblemdensity, /**< pointer to array to store the nonzeros for the subproblems */
   int*                  varmasterdensity,   /**< pointer to array to store the nonzeros for the master */
   int*                  conssubproblemdensity, /**< pointer to array to store the nonzeros for the subproblems */
   int*                  consmasterdensity   /**< pointer to array to store the nonzeros for the master */
);

/**
 *  calculates the number of up and down locks of variables for a given decomposition in both the original problem and the pricingproblems
 *
 *  @note All arrays need to be allocated by the caller
 *
 *  @warning This function needs a lot of memory (nvars*nblocks+1) array entries
 */
GCG_EXPORT
SCIP_RETCODE GCGgetVarLockData(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< pointer to array store variables belonging to density */
   int                   nvars,              /**< number of variables */
   int                   nsubproblems,       /**< number of sub problems */
   int**                 subsciplocksdown,   /**< pointer to two dimensional array to store the down locks for the subproblems */
   int**                 subsciplocksup,     /**< pointer to two dimensional array to store the down locks for the subproblems */
   int*                  masterlocksdown,    /**< pointer to array to store the down locks for the master */
   int*                  masterlocksup       /**< pointer to array to store the down locks for the master */
   );

/**
 * returns the maximum white score ( if it is not calculated yet is decomp is evaluated)
 */
GCG_EXPORT
SCIP_Real GCGdecompGetMaxWhiteScore(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );


/** computes the score of the given decomposition based on the border, the average density score and the ratio of
 * linking variables
 */
GCG_EXPORT
SCIP_RETCODE GCGdecompEvaluateDecomposition(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   DEC_SCORES*           score               /**< returns the score of the decomposition */
   );

/** returns the number of constraints saved in the decomposition */
GCG_EXPORT
int GCGdecompGetNConss(
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** display statistics about the decomposition */
GCG_EXPORT
SCIP_RETCODE GCGprintDecompStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file or NULL for standard output */
   GCG_DECOMP*           decomp              /**< decomp that should be evaluated */
   );

/** returns whether both structures lead to the same decomposition */
GCG_EXPORT
SCIP_Bool GCGdecompositionsAreEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp1,            /**< first decomp data structure */
   GCG_DECOMP*           decomp2             /**< second decomp data structure */
);


/** filters similar decompositions from a given list and moves them to the end
 * @return the number of unique decompositions
 */
GCG_EXPORT
int GCGfilterSimilarDecompositions(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP**          decs,               /**< array of decompositions */
   int                   ndecs               /**< number of decompositions */
);

/** returns the number of the block that the constraint is with respect to the decomposition */
GCG_EXPORT
SCIP_RETCODE GCGdetermineConsBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   int                   *block              /**< block of the constraint (or nblocks for master) */
);

/** move a master constraint to pricing problem */
GCG_EXPORT
SCIP_RETCODE GCGdecompMoveLinkingConsToPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   int                   consindex,          /**< index of constraint to move */
   int                   block               /**< block of the pricing problem where to move */
   );

/** tries to assign masterconss to pricing problem */
GCG_EXPORT
SCIP_RETCODE GCGdecompTryAssignMasterconssToExistingPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   int*                  transferred         /**< number of master constraints reassigned */
   );

/** tries to assign masterconss to new pricing problem */
GCG_EXPORT
SCIP_RETCODE GCGdecompTryAssignMasterconssToNewPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   GCG_DECOMP**          newdecomp,          /**< new decomposition, if successful */
   int*                  transferred         /**< number of master constraints reassigned */
   );

/** polish the decomposition and try to greedily assign master constraints to pricing problem where useful */
GCG_EXPORT
SCIP_RETCODE GCGcreatePolishedDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   GCG_DECOMP**          newdecomp           /**< new decomposition, if successful */
   );

/** permutes the decomposition according to the permutation seed */
GCG_EXPORT
SCIP_RETCODE GCGpermuteDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< decomposition data structure */
   SCIP_RANDNUMGEN*      randnumgen          /**< random number generator */
   );

/** @} */
#ifdef __cplusplus
}

#endif
#endif
