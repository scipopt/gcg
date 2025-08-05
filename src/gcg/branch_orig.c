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

/**@file   branch_orig.c
 * @brief  branching rule for the original problem in GCG
 * @author Gerald Gamrath
 * @author Marcel Schmickerath
 * @author Christian Puchert
 * @author Jonas Witt
 * @author Oliver Gaul
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/
#include <assert.h>
#include <string.h>

#include "gcg/branch_orig.h"
#include "gcg/gcg.h"
#include "gcg/branch_relpsprob.h"
#include "gcg/branch_bpstrong.h"
#include "gcg/cons_integralorig.h"
#include "gcg/cons_masterbranch.h"
#include "gcg/cons_origbranch.h"
#include "gcg/relax_gcg.h"
#include "gcg/pricer_gcg.h"
#include "gcg/type_branchgcg.h"

#include "scip/cons_linear.h"


#define BRANCHRULE_NAME          "orig"
#define BRANCHRULE_DESC          "branching for the original program in generic column generation"
#define BRANCHRULE_PRIORITY      100
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define DEFAULT_ENFORCEBYCONS  FALSE
#define DEFAULT_USEPSEUDO      TRUE
#define DEFAULT_MOSTFRAC       FALSE
#define DEFAULT_USERANDOM      FALSE
#define DEFAULT_USEPSSTRONG    FALSE

/* strong branching */
#define DEFAULT_USESTRONG        FALSE

#define DEFAULT_MINPHASE0OUTCANDS      10
#define DEFAULT_MAXPHASE0OUTCANDS      50
#define DEFAULT_MAXPHASE0OUTCANDSFRAC  0.7
#define DEFAULT_PHASE1GAPWEIGHT        0.25

#define DEFAULT_MINPHASE1OUTCANDS      3
#define DEFAULT_MAXPHASE1OUTCANDS      20
#define DEFAULT_MAXPHASE1OUTCANDSFRAC  0.7
#define DEFAULT_PHASE2GAPWEIGHT        1
/**/


/** branching rule data */
struct SCIP_BranchruleData
{
   GCG*                  gcg;                   /**< GCG data structure */
   SCIP_Bool             enforcebycons;         /**< should bounds on variables be enforced by constraints(TRUE) or by
                                                   * bounds(FALSE) */
   SCIP_Bool             usepseudocosts;        /**< should pseudocosts be used to determine the variable on which the
                                                   * branching is performed? */
   SCIP_Bool             mostfrac;              /**< should branching be performed on the most fractional variable?
                                                   * (only if usepseudocosts = FALSE) */
   SCIP_Bool             userandom;             /**< should the variable on which the branching is performed be
                                                   * selected randomly? (only if usepseudocosts = mostfrac = FALSE) */
   SCIP_Bool             usepsstrong;           /**< should strong branching with propagation be used to determine the
                                                   * variable on which the branching is performed?
                                                   * (only if usepseudocosts = mostfrac = random = FALSE)*/
   SCIP_Bool             usestrong;             /**< should strong branching be used to determine the variable on which
                                                   * the branching is performed? */
   GCG_BRANCHRULE*       gcgbranchrule;         /**< GCG branchrule structure */
};

/** branching data for branching decisions */
struct GCG_BranchData
{
   SCIP_VAR*             origvar;            /**< original variable on which the branching is done */
   GCG_BOUNDTYPE         boundtype;          /**< type of the new bound of original variable */
   SCIP_Real             newbound;           /**< new lower/upper bound of the original variable */
   SCIP_Real             oldbound;           /**< old lower/upper bound of the pricing variable */
   SCIP_Real             oldvalue;           /**< old value of the original variable */
   SCIP_Real             olddualbound;       /**< dual bound before the branching was performed */
   SCIP_CONS*            cons;               /**< constraint that enforces the branching restriction in the original
                                              *   problem, or NULL if this is done by variable bounds */
};

/* returns TRUE iff
 * iter = 0 and branchcand is an integer variable belonging to a unique block with fractional value, or
 * iter = 1 and branchcand is an integer variable that belongs to no block but was directly transferred to the
 *            master problem and which has a fractional value in the current solution
 */
static
SCIP_Bool getUniqueBlockFlagForIter(
   GCG* gcg,
   SCIP_VAR* branchcand,
   int iter
   )
{
   SCIP* scip = GCGgetOrigprob(gcg);
   assert(GCGvarIsOriginal(branchcand));
   /* continue if variable belongs to a block in second iteration*/
   if (iter == 0)
   {
      /* variable belongs to no block */
      if (GCGvarGetBlock(branchcand) == -1)
         return FALSE;

      /* block is not unique (non-linking variables) */
      if (!GCGoriginalVarIsLinking(branchcand) && GCGgetNIdenticalBlocks(gcg, GCGvarGetBlock(branchcand)) != 1)
         return FALSE;

      /* check that blocks of linking variable are unique */
      if (GCGoriginalVarIsLinking(branchcand))
      {
         int nvarblocks;
         int *varblocks;
         SCIP_Bool unique;
         int j;

         nvarblocks = GCGlinkingVarGetNBlocks(branchcand);
         SCIP_CALL(SCIPallocBufferArray(scip, &varblocks, nvarblocks));
         SCIP_CALL(GCGlinkingVarGetBlocks(branchcand, nvarblocks, varblocks));

         unique = TRUE;
         for (j = 0; j < nvarblocks; ++j)
            if (GCGgetNIdenticalBlocks(gcg, varblocks[j]) != 1)
               unique = FALSE;

         SCIPfreeBufferArray(scip, &varblocks);

         if (!unique)
            return FALSE;
      }
      /* candidate is valid in first iteration */
      return TRUE;

   }
   else /* iter == 1 */
   {
      if (GCGvarGetBlock(branchcand) != -1)
         return FALSE;

      /* candidate is valid in second iteration */
      return TRUE;
   }
   return FALSE;
}

/** branches on a integer variable x
 *  if solution value x' is fractional, two child nodes will be created
 *  (x <= floor(x'), x >= ceil(x')),
 *  if solution value is integral, the bounds of x are finite, then two child nodes will be created
 *  (x <= x", x >= x"+1 with x" = floor((lb + ub)/2)),
 *  otherwise (up to) three child nodes will be created
 *  (x <= x'-1, x == x', x >= x'+1)
 *  if solution value is equal to one of the bounds and the other bound is infinite, only two child nodes
 *  will be created (the third one would be infeasible anyway)
 */
static
SCIP_RETCODE branchVar(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE*       branchrule,         /**< pointer to the original variable branching rule */
   SCIP_VAR*             branchvar,          /**< variable to branch on */
   SCIP_Real             solval,             /**< value of the variable in the current solution */
   SCIP_Bool             upinf,              /**< have we seen during strong branching that the upbranch is
                                                * infeasible? */
   SCIP_Bool             downinf             /**< have we seen during strong branching that the downbranch is
                                                * infeasible? */
   )
{
   SCIP* origprob;
   SCIP* masterprob;
   /* data for b&b child creation */
   SCIP_Real downub;
   SCIP_Real fixval;
   SCIP_Real uplb;
   SCIP_BRANCHRULEDATA* branchruledata;

   origprob = GCGgetOrigprob(gcg);
   masterprob = GCGgetMasterprob(gcg);

   branchruledata = GCGbranchGetScipBranchruledata(branchrule);
   assert(branchruledata != NULL);

   assert(strcmp(SCIPbranchruleGetName(GCGbranchGetScipBranchrule(branchrule)), BRANCHRULE_NAME) == 0);

   assert(origprob != NULL);
   assert(branchrule != NULL);
   assert(branchvar != NULL);

   downub = SCIP_INVALID;
   fixval = SCIP_INVALID;
   uplb = SCIP_INVALID;

   if( SCIPisFeasIntegral(origprob, solval) )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbLocal(branchvar);
      ub = SCIPvarGetUbLocal(branchvar);

      /* if there was no explicit value given for branching, the variable has a finite domain and the current LP/pseudo
       * solution is one of the bounds, we branch in the center of the domain */
      if( !SCIPisInfinity(origprob, -lb) && !SCIPisInfinity(origprob, ub) )
      {
         SCIP_Real center;

         /* create child nodes with x <= x", and x >= x"+1 with x" = floor((lb + ub)/2);
          * if x" is integral, make the interval smaller in the child in which the current solution x'
          * is still feasible
          */
         center = (ub + lb) / 2.0;
         if( solval <= center )
         {
            downub = SCIPfeasFloor(origprob, center);
            uplb = downub + 1.0;
         }
         else
         {
            uplb = SCIPfeasCeil(origprob, center);
            downub = uplb - 1.0;
         }
      }
      else
      {
         /* create child nodes with x <= x'-1, x = x', and x >= x'+1 */
         assert(SCIPisEQ(origprob, SCIPfeasCeil(origprob, solval), SCIPfeasFloor(origprob, solval)));

         fixval = solval;

         /* create child node with x <= x'-1, if this would be feasible */
         if( SCIPisFeasGE(origprob, fixval-1.0, lb) )
            downub = fixval - 1.0;

         /* create child node with x >= x'+1, if this would be feasible */
         if( SCIPisFeasLE(origprob, fixval+1.0, ub) )
            uplb = fixval + 1.0;
      }
      SCIPdebugMessage("integral branch on variable <%s> with value %g, priority %d (current lower bound: %g)\n",
                       SCIPvarGetName(branchvar), solval, SCIPvarGetBranchPriority(branchvar),
                       SCIPgetLocalLowerbound(masterprob));
   }
   else
   {
      /* create child nodes with x <= floor(x'), and x >= ceil(x') */
      downub = SCIPfeasFloor(origprob, solval);
      uplb = downub + 1.0;
      assert( SCIPisEQ(origprob, SCIPfeasCeil(origprob, solval), uplb) );
   }

   if( uplb != SCIP_INVALID && !upinf ) /*lint !e777*/
   {
      SCIP_CONS* cons;
      SCIP_NODE* child;
      SCIP_CONS** origbranchconss;
      GCG_BRANCHDATA* branchdata;
      char name[SCIP_MAXSTRLEN];

      int norigbranchconss;
      int maxorigbranchconss;

      origbranchconss = NULL;
      norigbranchconss = 0;
      maxorigbranchconss = 0;

      /* create child node x >= uplb */
      SCIP_CALL( SCIPcreateChild(masterprob, &child, 0.0, SCIPgetLocalTransEstimate(masterprob)) );

      SCIP_CALL( SCIPallocBlockMemory(origprob, &branchdata) );

      branchdata->origvar = branchvar;
      branchdata->oldvalue = solval;
      branchdata->olddualbound = SCIPgetLocalLowerbound(masterprob);
      branchdata->boundtype = GCG_BOUNDTYPE_LOWER;
      branchdata->newbound = uplb;
      branchdata->oldbound = SCIPvarGetLbLocal(branchvar);

      SCIPdebugMessage(" -> creating child: <%s> >= %g\n",
         SCIPvarGetName(branchvar), uplb);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchdata->origvar),
         ">=", branchdata->newbound);

      if( branchruledata->enforcebycons )
      {
         /* enforce new bounds by linear constraints */
         SCIP_CONS* consup;

         SCIPdebugMessage("enforced by cons\n");

         norigbranchconss = 1;
         maxorigbranchconss = SCIPcalcMemGrowSize(origprob, 1);
         SCIP_CALL( SCIPallocBlockMemoryArray(origprob, &origbranchconss, maxorigbranchconss) );

         /* create corresponding constraints */
         SCIP_CALL( SCIPcreateConsLinear(origprob, &consup, name, 0, NULL, NULL,
            SCIPceil(origprob, solval), SCIPinfinity(origprob),
            TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );

         SCIP_CALL( SCIPaddCoefLinear(origprob, consup, branchvar, 1.0) );

         origbranchconss[0] = consup;
         branchdata->cons = consup;
      }
      else
         branchdata->cons = NULL;

      /* create and add the masterbranch constraint */
      SCIP_CALL( GCGcreateConsMasterbranch(gcg, &cons, name, child,
         GCGconsMasterbranchGetActiveCons(gcg), branchrule, branchdata, origbranchconss, norigbranchconss,
         maxorigbranchconss) );
      SCIP_CALL( SCIPaddConsNode(masterprob, child, cons, NULL) );
   }

   if( downub != SCIP_INVALID && !downinf ) /*lint !e777*/
   {
      SCIP_CONS* cons;
      SCIP_NODE* child;
      SCIP_CONS** origbranchconss;
      GCG_BRANCHDATA* branchdata;
      char name[SCIP_MAXSTRLEN];

      int norigbranchconss;
      int maxorigbranchconss;

      origbranchconss = NULL;
      norigbranchconss = 0;
      maxorigbranchconss = 0;

      /* create child node x <= downub */
      SCIP_CALL( SCIPcreateChild(masterprob, &child, 0.0, SCIPgetLocalTransEstimate(masterprob)) );

      SCIP_CALL( SCIPallocBlockMemory(origprob, &branchdata) );

      branchdata->origvar = branchvar;
      branchdata->oldvalue = solval;
      branchdata->olddualbound = SCIPgetLocalLowerbound(masterprob);
      branchdata->boundtype = GCG_BOUNDTYPE_UPPER;
      branchdata->newbound = downub;
      branchdata->oldbound = SCIPvarGetUbLocal(branchvar);

      SCIPdebugMessage(" -> creating child: <%s> <= %g\n",
         SCIPvarGetName(branchvar), downub);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchdata->origvar),
         "<=", branchdata->newbound);

      /* enforce branching decision by a constraint rather than by bound changes */
      if( branchruledata->enforcebycons )
      {
         /* enforce new bounds by linear constraints */
         SCIP_CONS* consdown;

         norigbranchconss = 1;
         maxorigbranchconss = SCIPcalcMemGrowSize(origprob, 1);
         SCIP_CALL( SCIPallocBlockMemoryArray(origprob, &origbranchconss, maxorigbranchconss) );

         /* create corresponding constraints */
         SCIP_CALL( SCIPcreateConsLinear(origprob, &consdown, name, 0, NULL, NULL,
            -1.0 * SCIPinfinity(origprob), SCIPfloor(origprob, solval),
            TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPaddCoefLinear(origprob, consdown, branchvar, 1.0) );

         origbranchconss[0] = consdown;
         branchdata->cons = consdown;
      }
      else
         branchdata->cons = NULL;

      /* create and add the masterbranch constraint */
      SCIP_CALL( GCGcreateConsMasterbranch(gcg, &cons, name, child,
         GCGconsMasterbranchGetActiveCons(gcg), branchrule, branchdata, origbranchconss, norigbranchconss,
         maxorigbranchconss) );
      SCIP_CALL( SCIPaddConsNode(masterprob, child, cons, NULL) );
   }

   if( fixval != SCIP_INVALID ) /*lint !e777*/
   {
      SCIP_CONS* cons;
      SCIP_NODE* child;
      SCIP_CONS** origbranchconss;
      GCG_BRANCHDATA* branchdata;
      char name[SCIP_MAXSTRLEN];

      int norigbranchconss;
      int maxorigbranchconss;

      origbranchconss = NULL;
      norigbranchconss = 0;
      maxorigbranchconss = 0;

      /* create child node x = fixval */
      SCIP_CALL( SCIPcreateChild(masterprob, &child, 0.0, SCIPgetLocalTransEstimate(masterprob)) );

      SCIP_CALL( SCIPallocBlockMemory(origprob, &branchdata) );

      branchdata->origvar = branchvar;
      branchdata->oldvalue = solval;
      branchdata->olddualbound = SCIPgetLocalLowerbound(masterprob);
      branchdata->boundtype = GCG_BOUNDTYPE_FIXED;
      branchdata->newbound = fixval;
      branchdata->oldbound = SCIPvarGetUbLocal(branchvar);

      SCIPdebugMessage(" -> creating child: <%s> == %g\n",
         SCIPvarGetName(branchvar), fixval);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchdata->origvar),
         "==", branchdata->newbound);

      /* enforce branching decision by a constraint rather than by bound changes */
      if( branchruledata->enforcebycons )
      {
         /* enforce new bounds by linear constraints */
         SCIP_CONS* consfix;

         norigbranchconss = 1;
         maxorigbranchconss = SCIPcalcMemGrowSize(origprob, 1);
         SCIP_CALL( SCIPallocBlockMemoryArray(origprob, &origbranchconss, maxorigbranchconss) );

         /* create corresponding constraints */
         SCIP_CALL( SCIPcreateConsLinear(origprob, &consfix, name, 0, NULL, NULL,
            fixval, fixval, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPaddCoefLinear(origprob, consfix, branchvar, 1.0) );

         origbranchconss[0] = consfix;
         branchdata->cons = consfix;
      }
      else
         branchdata->cons = NULL;

      /* create and add the masterbranch constraint */
      SCIP_CALL( GCGcreateConsMasterbranch(gcg, &cons, name, child,
         GCGconsMasterbranchGetActiveCons(gcg), branchrule, branchdata, origbranchconss, norigbranchconss,
         maxorigbranchconss) );
      SCIP_CALL( SCIPaddConsNode(masterprob, child, cons, NULL) );
   }

   return SCIP_OKAY;
}

/* Evaluates the given variable based on a score function of choice. Higher scores are given to better
 * variables.
 */
static SCIP_Real score_function(
    SCIP *scip,
    SCIP_BRANCHRULEDATA*      branchruledata,         /* pointer to branching rule data */
    SCIP_VAR *var,            /* var to be scored */
    SCIP_Real solval,         /* the var's current solution value */
    SCIP_Real *score         /* stores the computed score */
   )
{
   assert(branchruledata != NULL);

   /* define score functions and calculate score for all variables for sorting dependent on used heuristic */
   if( branchruledata->usepseudocosts )
   {
      *score = SCIPgetVarPseudocostScore(scip, var, solval);
   }
   else
   {
      if( !branchruledata->mostfrac )
         return 1;

      *score = solval - SCIPfloor(scip, solval);
      *score = MIN(*score, 1.0 - *score);
   }

   return SCIP_OKAY;
}

/** branching method for relaxation solutions */
static
SCIP_RETCODE branchExtern(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE*       branchrule,         /**< pointer to the original variable branching rule */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching call */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP* origprob;

   /* branching candidates */
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   int nbranchcands;
   int npriobranchcands;

   SCIP_Bool upinf;
   SCIP_Bool downinf;

   /* values for choosing the variable to branch on */
   SCIP_VAR* branchvar;
   SCIP_Real solval;

   SCIP_Real maxscore;
   SCIP_Real score;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(GCGbranchGetScipBranchrule(branchrule)), BRANCHRULE_NAME) == 0);
   assert(gcg != NULL);
   assert(result != NULL);

   origprob = GCGgetOrigprob(gcg);

   assert(SCIPisRelaxSolValid(origprob));

   branchruledata = GCGbranchGetScipBranchruledata(branchrule);
   assert(branchruledata != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get the branching candidates */
   SCIP_CALL( SCIPgetExternBranchCands(origprob, &branchcands, &branchcandssol, NULL, &nbranchcands,
         &npriobranchcands, NULL, NULL, NULL) );

   branchvar = NULL;
   solval = 0.0;

   upinf = FALSE;
   downinf = FALSE;

   maxscore = -1.0;

   if( !branchruledata->usestrong )
   {
      if( branchruledata->usepseudocosts || branchruledata->mostfrac || branchruledata->userandom )
      {
         /* iter = 0: integer variables belonging to a unique block with fractional value,
          * iter = 1: we did not find enough variables to branch on so far, so we look for integer variables that
          * belong to no block
          * but were directly transferred to the master problem and which have a fractional value in the current
          * solution
          */
         for( int iter = 0; iter <= 1 && branchvar == NULL; iter++ )
         {
            for( int i = 0; i < npriobranchcands; i++ )
            {
               if( !getUniqueBlockFlagForIter(gcg, branchcands[i], iter) )
                  continue;

               if( !branchruledata->userandom )
               {
                  SCIP_CALL( score_function(origprob, branchruledata, branchcands[i], branchcandssol[i], &score) );

                  if( score > maxscore )
                  {
                     maxscore = score;
                     branchvar = branchcands[i];
                     solval = branchcandssol[i];
                  }
               }
               else
               {
                  branchvar = branchcands[i];
                  solval = branchcandssol[i];
                  break;
               }
            }
         }
      }
      else if( branchruledata->usepsstrong )
      {
         SCIP_CALL( GCGgetRelpsprobBranchVar(gcg, branchcands, branchcandssol, npriobranchcands,
               npriobranchcands, result, &branchvar) );
         assert(branchvar != NULL || *result == SCIP_CUTOFF);
         assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF);

         if( *result == SCIP_CUTOFF )
            return SCIP_OKAY;

         solval = SCIPgetRelaxSolVal(origprob, branchvar);
      }


   }
   else
   {
      SCIP_CALL( GCGbranchSelectCandidateStrongBranchingOrig(gcg, branchrule, &branchvar, &upinf, &downinf, result,
                                                             &branchruledata->usestrong) );
   }

   if( upinf && downinf )
      return SCIP_OKAY;

   if( branchvar == NULL )
   {
      SCIPdebugMessage("Original branching rule could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(!(upinf && downinf));

   SCIPdebugMessage("Original branching rule selected variable %s%s\n",
                    SCIPvarGetName(branchvar), (upinf || downinf)? ", which is infeasible in one direction" : "");

   SCIP_CALL( branchVar(gcg, branchrule, branchvar, solval, upinf, downinf) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/*
 * Callback methods for enforcing branching constraints
 */

#define branchDeactiveMasterOrig NULL
#define branchPropMasterOrig NULL

/** callback activation method */
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterOrig)
{
   SCIP* masterprob;
   SCIP_CONS* mastercons;

   masterprob = GCGgetMasterprob(gcg);

   assert(masterprob != NULL);
   assert(branchdata != NULL);

   /* branching restrictions are enforced by variable bounds, this is done automatically, so we can abort here */
   if( branchdata->cons == NULL )
      return SCIP_OKAY;

   assert(branchdata->origvar != NULL);   

   SCIPdebugMessage("branchActiveMasterOrig: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == GCG_BOUNDTYPE_LOWER ?
                            ">=" : branchdata->boundtype == GCG_BOUNDTYPE_UPPER ? "<=" : "==" ), branchdata->newbound);

   /* transform constraint to the master variable space */
   SCIP_CALL( GCGrelaxTransOrigToMasterCons(gcg, branchdata->cons, &mastercons) );
   assert(mastercons != NULL);

   /* add constraint to the master problem */
   SCIP_CALL( SCIPaddConsNode(masterprob, SCIPgetCurrentNode(masterprob), mastercons, NULL) );

   /* constraint was added locally to the node where it is needed, so we do not need to care about this
    * at the next activation of the node and can set the constraint pointer to NULL */
   SCIP_CALL( SCIPreleaseCons(masterprob, &branchdata->cons) );
   branchdata->cons = NULL;

   return SCIP_OKAY;
}

/** callback solved method */
static
GCG_DECL_BRANCHMASTERSOLVED(branchMasterSolvedOrig)
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   SCIP* masterprob = GCGgetMasterprob(gcg);
   assert(origprob != NULL);
   assert(GCGisOriginal(origprob));
   assert(branchdata != NULL);
   assert(branchdata->origvar != NULL);

   SCIPdebugMessage("branchMasterSolvedOrig: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == GCG_BOUNDTYPE_LOWER ?
                            ">=" : branchdata->boundtype == GCG_BOUNDTYPE_UPPER ? "<=" : "==" ), branchdata->newbound);

   if( !SCIPisInfinity(origprob, newlowerbound) && SCIPgetStage(masterprob) == SCIP_STAGE_SOLVING
      && SCIPisRelaxSolValid(masterprob) )
   {
      SCIP_CALL( SCIPupdateVarPseudocost(origprob, branchdata->origvar,
            SCIPgetRelaxSolVal(origprob, branchdata->origvar) - branchdata->oldvalue,
            newlowerbound - branchdata->olddualbound, 1.0) );
   }

   return SCIP_OKAY;
}

/** callback deletion method for branching data */
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteOrig)
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);
   assert(branchdata != NULL);

   if( !origbranch && !force )
      return SCIP_OKAY;

   if( *branchdata == NULL )
      return SCIP_OKAY;

   SCIPdebugMessage("branchDataDeleteOrig: %s %s %f\n", SCIPvarGetName((*branchdata)->origvar),
      ( (*branchdata)->boundtype == GCG_BOUNDTYPE_LOWER ?
        ">=" : (*branchdata)->boundtype == GCG_BOUNDTYPE_UPPER ? "<=" : "==" ),
      (*branchdata)->newbound);

   /* release constraint */
   if( (*branchdata)->cons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(origprob, &(*branchdata)->cons) );
      (*branchdata)->cons = NULL;
   }

   SCIPfreeBlockMemoryNull(origprob, branchdata);
   *branchdata = NULL;

   return SCIP_OKAY;
}

static
SCIP_DECL_BRANCHFREE(branchFreeOrig)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);

   SCIPfreeBlockMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpOrig)
{  /*lint --e{715}*/
   SCIP* origprob;
   SCIP_BRANCHRULEDATA* branchruledata;

   //SCIPdebugMessage("Execlp method of orig branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get original problem */
   origprob = GCGgetOrigprob(branchruledata->gcg);
   assert(origprob != NULL);

   if( GCGcurrentNodeIsGeneric(branchruledata->gcg) )
   {
      SCIPdebugMessage("Not executing orig branching, node was branched by generic branchrule\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* if the transferred master solution is feasible, the current node is solved to optimality and can be pruned */
   if( GCGrelaxIsOrigSolFeasible(branchruledata->gcg) )
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMessage("solution was feasible, node can be cut off!");
   }

   if( SCIPgetNExternBranchCands(origprob) > 0 )
   {
      assert(SCIPisRelaxSolValid(origprob));
      SCIP_CALL( branchExtern(branchruledata->gcg, branchruledata->gcgbranchrule, result) );
   }

   return SCIP_OKAY;
}

/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextOrig)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   SCIPdebugMessage("Execext method of orig branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( GCGcurrentNodeIsGeneric(branchruledata->gcg) )
   {
      SCIPdebugMessage("Not executing orig branching, node was branched by generic branchrule\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* if the transferred master solution is feasible, the current node is solved to optimality and can be pruned */
   if( GCGrelaxIsOrigSolFeasible(branchruledata->gcg) )
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMessage("solution was feasible, node can be cut off!");
   }
   SCIP_CALL( branchExtern(branchruledata->gcg, branchruledata->gcgbranchrule, result) );

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitOrig)
{
   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsOrig)
{  /*lint --e{715}*/
   int i;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP* origprob;

   /* branching candidates */
   SCIP_VAR** branchcands;
   int nbranchcands;
   int npriobranchcands;

   /* values for choosing the variable to branch on */
   SCIP_VAR* branchvar;
   SCIP_Real solval;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPdebugMessage("EXECPS orig branching rule\n");

   /* get original problem */
   origprob = GCGgetOrigprob(branchruledata->gcg);
   assert(origprob != NULL);

   if( GCGcurrentNodeIsGeneric(branchruledata->gcg) )
   {
      SCIPdebugMessage("Not executing orig branching, node was branched by generic branchrule\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Execps method of orig branching\n");

   *result = SCIP_DIDNOTRUN;
   if( SCIPgetStage(scip) > SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* get the branching candidates */
   SCIP_CALL( SCIPgetPseudoBranchCands(origprob, &branchcands, &nbranchcands, &npriobranchcands) );

   branchvar = NULL;
   solval = 0.0;

   /* branch on an integer variable belonging to a unique block with fractional value */
   for( i = 0; i < npriobranchcands; i++ )
   {
      assert(GCGvarIsOriginal(branchcands[i]));

      /* variable belongs to no block or the block is not unique */
      if( GCGvarGetBlock(branchcands[i]) <= -1 || GCGgetNIdenticalBlocks(branchruledata->gcg, GCGvarGetBlock(branchcands[i])) != 1 )
         continue;

      branchvar = branchcands[i];
      lb = SCIPvarGetLbLocal(branchvar);
      ub = SCIPvarGetUbLocal(branchvar);

      assert(ub - lb > 0.8);

      /*  if the bounds of the branching variable x are finite, then the solution value
       *  is floor((lb + ub)/2)) + 0.5,
       *  otherwise the solution value is set to a finite bound
       *  if no finite bound exists, the solution value is set to 0.
       */
      if( !SCIPisInfinity(origprob, ub) && !SCIPisInfinity(origprob, -lb) )
         solval =  SCIPfeasFloor(scip, (ub + lb) / 2.0) + 0.5;
      else if( !SCIPisInfinity(origprob, -lb) )
         solval = lb;
      else if( !SCIPisInfinity(origprob, ub) )
         solval = ub;
      else
         solval = 0.0;

      break;
   }

   /* we did not find a variable to branch on so far, so we look for an unfixed linking variable or an integer variable
    * that belongs to no block but was directly transferred to the master problem
    */
   if( branchvar == NULL )
   {
      for( i = 0; i < npriobranchcands; i++ )
      {
         assert(GCGvarIsOriginal(branchcands[i]));

         /* continue if variable belongs to a block */
         if( GCGvarGetBlock(branchcands[i]) > -1 )
            continue;

         /* check that blocks of linking variable are unique */
         if( GCGoriginalVarIsLinking(branchcands[i]) )
         {
            int nvarblocks;
            int* varblocks;
            SCIP_Bool unique;
            int j;

            nvarblocks = GCGlinkingVarGetNBlocks(branchcands[i]);
            SCIP_CALL( SCIPallocBufferArray(origprob, &varblocks, nvarblocks) );
            SCIP_CALL( GCGlinkingVarGetBlocks(branchcands[i], nvarblocks, varblocks) );

            unique = TRUE;
            for( j = 0; j < nvarblocks; ++j )
               if( GCGgetNIdenticalBlocks(branchruledata->gcg, varblocks[j]) != 1 )
                  unique = FALSE;

            SCIPfreeBufferArray(origprob, &varblocks);

            if( !unique )
               continue;
         }

         branchvar = branchcands[i];
         lb = SCIPvarGetLbLocal(branchvar);
         ub = SCIPvarGetUbLocal(branchvar);

         assert(ub - lb > 0.8);

         /*  if the bounds of the branching variable x are finite, then the solution value
          *  is floor((lb + ub)/2)) + 0.5,
          *  otherwise the solution value is set to a finite bound
          *  if no finite bound exists, the solution value is set to 0.
          */
         if( !SCIPisInfinity(origprob, ub) && !SCIPisInfinity(origprob, -lb) )
            solval =  SCIPfeasFloor(scip, (ub + lb) / 2.0) + 0.5;
         else if( !SCIPisInfinity(origprob, -lb) )
            solval = lb;
         else if( !SCIPisInfinity(origprob, ub) )
            solval = ub;
         else
            solval = 0.0;

         break;
      }
   }

   if( branchvar == NULL )
   {
      SCIPdebugMessage("Original branching rule could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(branchvar != NULL);

   SCIP_CALL( branchVar(branchruledata->gcg, branchruledata->gcgbranchrule, branchvar, solval, FALSE, FALSE) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}



/*
 * branching specific interface methods
 */

/** creates the branching on original variable branching rule and includes it in SCIP */
SCIP_RETCODE GCGincludeBranchruleOrig(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* origprob;
   SCIP* masterprob;
   SCIP_BRANCHRULE* branchrule;
   GCG_BRANCHRULE* gcgbranchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   SCIPdebugMessage("Include orig branching rule\n");
   assert(gcg != NULL);

   /* get original problem */
   origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);
   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);

   /* alloc branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(masterprob, &branchruledata) );

   branchruledata->gcg = gcg;

   /* include branching rule */
   SCIP_CALL( GCGrelaxIncludeBranchrule(branchruledata->gcg, &branchrule, &gcgbranchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata, branchActiveMasterOrig,
         branchDeactiveMasterOrig, branchPropMasterOrig, branchMasterSolvedOrig, branchDataDeleteOrig, NULL, NULL, NULL) );
   assert(branchrule != NULL);
   branchruledata->gcgbranchrule = gcgbranchrule;

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleInit(masterprob, branchrule, branchInitOrig) );
   SCIP_CALL( SCIPsetBranchruleExecLp(masterprob, branchrule, branchExeclpOrig) );
   SCIP_CALL( SCIPsetBranchruleExecExt(masterprob, branchrule, branchExecextOrig) );
   SCIP_CALL( SCIPsetBranchruleExecPs(masterprob, branchrule, branchExecpsOrig) );
   SCIP_CALL( SCIPsetBranchruleFree(masterprob, branchrule, branchFreeOrig) );

   /* add original variable branching rule parameters */

   SCIP_CALL( SCIPaddBoolParam(origprob, "branching/orig/enforcebycons",
         "should bounds on variables be enforced by constraints(TRUE) or by bounds(FALSE)",
         &branchruledata->enforcebycons, FALSE, DEFAULT_ENFORCEBYCONS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "branching/orig/usepseudocosts",
         "should pseudocosts be used to determine the variable on which the branching is performed?",
         &branchruledata->usepseudocosts, FALSE, DEFAULT_USEPSEUDO, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "branching/orig/mostfrac",
         "should branching be performed on the most fractional variable? (only if usepseudocosts = FALSE)",
         &branchruledata->mostfrac, FALSE, DEFAULT_MOSTFRAC, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "branching/orig/userandom",
         "should the variable on which the branching is performed be selected randomly? (only if usepseudocosts = mostfrac = FALSE)",
         &branchruledata->userandom, FALSE, DEFAULT_USERANDOM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "branching/orig/usepsstrong",
         "should strong branching with propagation be used to determine the variable on which the branching is performed? (only if usepseudocosts = mostfrac = random = FALSE)",
         &branchruledata->usepsstrong, FALSE, DEFAULT_USEPSSTRONG, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "branching/orig/usestrong",
         "should strong branching be used to determine the variable on which the branching is performed?",
         &branchruledata->usestrong, FALSE, DEFAULT_USESTRONG, NULL, NULL) );

   /* strong branching */
   SCIP_CALL( SCIPaddIntParam(origprob, "branching/orig/minphase0outcands",
         "minimum number of output candidates from phase 0 during strong branching",
         NULL, FALSE, DEFAULT_MINPHASE0OUTCANDS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "branching/orig/maxphase0outcands",
         "maximum number of output candidates from phase 0 during strong branching",
         NULL, FALSE, DEFAULT_MAXPHASE0OUTCANDS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "branching/orig/maxphase0outcandsfrac",
         "maximum number of output candidates from phase 0 as fraction of total cands during strong branching",
         NULL, FALSE, DEFAULT_MAXPHASE0OUTCANDSFRAC, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "branching/orig/phase1gapweight",
         "how much impact should the node gap have on the number of precisely evaluated candidates in phase 1 during strong branching?",
         NULL, FALSE, DEFAULT_PHASE1GAPWEIGHT, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "branching/orig/minphase1outcands",
         "minimum number of output candidates from phase 1 during strong branching",
         NULL, FALSE, DEFAULT_MINPHASE1OUTCANDS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "branching/orig/maxphase1outcands",
         "maximum number of output candidates from phase 1 during strong branching",
         NULL, FALSE, DEFAULT_MAXPHASE1OUTCANDS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "branching/orig/maxphase1outcandsfrac",
         "maximum number of output candidates from phase 1 as fraction of phase 1 cands during strong branching",
         NULL, FALSE, DEFAULT_MAXPHASE1OUTCANDSFRAC, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "branching/orig/phase2gapweight",
         "how much impact should the node gap have on the number of precisely evaluated candidates in phase 2 during strong branching?",
         NULL, FALSE, DEFAULT_PHASE2GAPWEIGHT, 0, 1, NULL, NULL) );

   return SCIP_OKAY;
}

/** get the original variable on which the branching was performed */
SCIP_VAR* GCGbranchOrigGetOrigvar(
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   )
{
   assert(branchdata != NULL);

   return branchdata->origvar;
}

/** get the type of the new bound which resulted of the performed branching */
GCG_BOUNDTYPE GCGbranchOrigGetBoundtype(
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   )
{
   assert(branchdata != NULL);

   return branchdata->boundtype;
}

/** get the new bound which resulted of the performed branching */
SCIP_Real GCGbranchOrigGetNewbound(
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   )
{
   assert(branchdata != NULL);

   return branchdata->newbound;
}

/** updates extern branching candidates before branching */
SCIP_RETCODE GCGbranchOrigUpdateExternBranchcands(
   GCG*                  gcg                /**< GCG data structure */
)
{
   SCIP* scip;
   SCIP_VAR** origvars;
   int norigvars;
   int i;

   scip = GCGgetOrigprob(gcg);
   assert(GCGisOriginal(scip));

   origvars = SCIPgetVars(scip);
   norigvars = SCIPgetNVars(scip);
   assert(origvars != NULL);

   SCIPclearExternBranchCands(scip);

   /* store branching candidates */
   for( i = 0; i < norigvars; i++ )
   {
      if( SCIPvarGetType(origvars[i]) <= SCIP_VARTYPE_INTEGER && !SCIPisFeasIntegral(scip,
                         SCIPgetRelaxSolVal(scip, origvars[i])) )
      {
         assert(!SCIPisEQ(scip, SCIPvarGetLbLocal(origvars[i]), SCIPvarGetUbLocal(origvars[i])));

         SCIP_CALL( SCIPaddExternBranchCand(scip, origvars[i], SCIPgetRelaxSolVal(scip,
            origvars[i]) - SCIPfloor(scip, SCIPgetRelaxSolVal(scip, origvars[i])),
            SCIPgetRelaxSolVal(scip, origvars[i])) );
      }
   }
   SCIPdebugMessage("updated relaxation branching candidates\n");

   return SCIP_OKAY;
}
