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

struct ComponentBoundSequence
{
   SCIP_VAR* component;
   GCG_COMPSENSE sense;
   SCIP_Real bound;
};

typedef struct ComponentBoundSequence GCG_COMPSEQUENCE;

struct GCG_Strip
{
   SCIP*                scip;
   SCIP_VAR*            mastervar;
   GCG_COMPSEQUENCE**   C;             /**< current set of comp bound sequences */
   int                  Csize;
   int*                 sequencesizes;
};
typedef struct GCG_Strip GCG_STRIP;

/** creates the most infeasible LP branching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchruleGeneric(
   SCIP*                scip                /**< SCIP data structure */
   );

/** initializes branchdata */
extern
SCIP_RETCODE GCGbranchGenericCreateBranchdata(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_BRANCHDATA**     branchdata          /**< branching data to initialize */
   );

/** computes the generator of mastervar for the entry in origvar */
extern
SCIP_Real getGeneratorEntry(
   SCIP_VAR*            mastervar,          /**< current mastervariable */
   SCIP_VAR*            origvar             /**< corresponding origvar */
   );

extern
GCG_COMPSEQUENCE* GCGbranchGenericBranchdataGetConsS(
   GCG_BRANCHDATA*      branchdata          /**< branching data to initialize */
   );

extern
int GCGbranchGenericBranchdataGetConsSsize(
   GCG_BRANCHDATA*      branchdata          /**< branching data to initialize */
   );

extern
int GCGbranchGenericBranchdataGetConsblocknr(
   GCG_BRANCHDATA*      branchdata          /**< branching data to initialize */
   );

extern
SCIP_CONS* GCGbranchGenericBranchdataGetMastercons(
   GCG_BRANCHDATA*      branchdata          /**< branching data to initialize */
   );

/** prepares informations for using the generic branching scheme
 * @return SCIP_RETCODE */
extern
SCIP_RETCODE GCGbranchGenericInitbranch(
   SCIP*                masterscip,         /**< */
   SCIP_BRANCHRULE*     branchrule,
   SCIP_RESULT*         result,
   int*                 checkedblocks,
   int                  ncheckedblocks,
   GCG_STRIP***         checkedblockssortstrips,
   int*                 checkedblocksnsortstrips
   );

#ifdef __cplusplus
}
#endif

#endif