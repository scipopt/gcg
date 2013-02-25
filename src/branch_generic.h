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
   int component;
   GCG_COMPSENSE sense;
   SCIP_Real bound;
};

typedef struct ComponentBoundSequence GCG_COMPSEQUENCE;

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

/** computes the generator of mastervar */
extern
SCIP_RETCODE getGenerators(
   SCIP*                scip,               /**< */
   SCIP_Real**          generator,          /**< */
   int*                 generatorsize,      /**< */
   SCIP_Bool**          compisinteger,      /**< */
   int                  blocknr,            /**< */
   SCIP_VAR**           mastervars,         /**< */
   int                  nmastervars,        /**< */
   SCIP_VAR*            mastervar           /**< */
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
#ifdef __cplusplus
}
#endif

#endif
