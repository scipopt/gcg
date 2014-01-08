/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   branch_generic.h
 * @brief  branching rule based on vanderbeck's generic branching scheme
 * @author Marcel Schmickerath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_GENERIC_H__
#define __SCIP_BRANCH_GENERIC_H__


#include "scip/scip.h"
#include "type_branchgcg.h"

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

/** creates the most infeasible LP branching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchruleGeneric(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** initializes branchdata */
extern
SCIP_RETCODE GCGbranchGenericCreateBranchdata(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_BRANCHDATA**      branchdata          /**< branching data to initialize */
   );

/** computes the generator of mastervar for the entry in origvar */
extern
SCIP_Real getGeneratorEntry(
   SCIP_VAR*             mastervar,          /**< current mastervariable */
   SCIP_VAR*             origvar             /**< corresponding origvar */
   );

/** get component bound sequence */
extern
GCG_COMPSEQUENCE* GCGbranchGenericBranchdataGetConsS(
   GCG_BRANCHDATA*       branchdata          /**< branching data to initialize */
   );

/** get size of component bound sequence */
extern
int GCGbranchGenericBranchdataGetConsSsize(
   GCG_BRANCHDATA*       branchdata          /**< branching data to initialize */
   );

/** get id of pricing problem (or block) to which the constraint belongs */
extern
int GCGbranchGenericBranchdataGetConsblocknr(
   GCG_BRANCHDATA*       branchdata          /**< branching data to initialize */
   );

/** get master constraint */
extern
SCIP_CONS* GCGbranchGenericBranchdataGetMastercons(
   GCG_BRANCHDATA*       branchdata          /**< branching data to initialize */
   );

/** prepares informations for using the generic branching scheme */
extern
SCIP_RETCODE GCGbranchGenericInitbranch(
   SCIP*                 masterscip,         /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_RESULT*          result,             /**< pointer to store the result of the branching call */
   int*                  checkedblocks,
   int                   ncheckedblocks,
   GCG_STRIP***          checkedblockssortstrips,
   int*                  checkedblocksnsortstrips
   );

/** returns true when the branch rule is the generic branchrule */
SCIP_Bool GCGisBranchruleGeneric(
   SCIP_BRANCHRULE*      branchrule          /**< branchrule to check */
);

#ifdef __cplusplus
}
#endif

#endif
