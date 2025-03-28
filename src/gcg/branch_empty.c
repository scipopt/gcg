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

/**@file   branch_empty.c
 * @brief  branching rule for the original problem while real branching is applied in the master
 * @author Marcel Schmickerath
 * @author Martin Bergner
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */
#include <assert.h>
#include <string.h>

#include "gcg/branch_empty.h"
#include "gcg/relax_gcg.h"
#include "gcg/gcg.h"
#include "gcg/branch_orig.h"
#include "gcg/cons_masterbranch.h"
#include "gcg/cons_origbranch.h"
#include "scip/branch_allfullstrong.h"
#include "scip/branch_fullstrong.h"
#include "scip/branch_inference.h"
#include "scip/branch_mostinf.h"
#include "scip/branch_leastinf.h"
#include "scip/branch_pscost.h"
#include "scip/branch_random.h"
#include "scip/branch_relpscost.h"
#include "scip/cons_varbound.h"
#include "gcg/type_branchgcg.h"

#define BRANCHRULE_NAME          "empty"
#define BRANCHRULE_DESC          "branching rule for the original problem while real branching is applied in the master"
#define BRANCHRULE_PRIORITY      1000000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

struct SCIP_BranchruleData
{
   GCG*                 gcg;        /**< GCG data structure */
};

/*
 * Callback methods for enforcing branching constraints
 */

/** copy default SCIP branching rules in order to solve restrictions of the original problem as a subSCIP without
 *  Dantzig-Wolfe decomposition
 */
static
SCIP_RETCODE includeSCIPBranchingRules(
   SCIP*                 scip
)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeBranchruleAllfullstrong(scip) );
   SCIP_CALL( SCIPincludeBranchruleFullstrong(scip) );
   SCIP_CALL( SCIPincludeBranchruleInference(scip) );
   SCIP_CALL( SCIPincludeBranchruleMostinf(scip) );
   SCIP_CALL( SCIPincludeBranchruleLeastinf(scip) );
   SCIP_CALL( SCIPincludeBranchrulePscost(scip) );
   SCIP_CALL( SCIPincludeBranchruleRandom(scip) );
   SCIP_CALL( SCIPincludeBranchruleRelpscost(scip) );

   return SCIP_OKAY;
}

/** for a new branch-and-bound node on the master problem
 *  add an original branching constraint that holds the branching decision to the corresponding node in the original problem
 */
static
SCIP_RETCODE createOrigbranchConstraint(
   GCG*                  gcg,
   SCIP_NODE*            childnode,
   SCIP_CONS*            masterbranchchildcons
)
{
   SCIP* origprob;
   char* consname;
   SCIP_BRANCHRULE* branchrule;
   GCG_BRANCHDATA* branchdata;
   SCIP_CONS* origcons;
   SCIP_CONS** origbranchconss;
   int norigbranchconss;

   int i;

   origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);
   assert(masterbranchchildcons != NULL);

   /* get name and branching information from the corresponding masterbranch constraint */
   consname = GCGconsMasterbranchGetName(masterbranchchildcons);
   branchrule = GCGconsMasterbranchGetBranchrule(masterbranchchildcons);
   branchdata = GCGconsMasterbranchGetBranchdata(masterbranchchildcons);

   /* create an origbranch constraint and add it to the node */
   SCIPdebugMessage("Create original branching constraint %s\n", consname);
   SCIP_CALL( GCGcreateConsOrigbranch(gcg, &origcons, consname, childnode,
            GCGconsOrigbranchGetActiveCons(gcg), branchrule, branchdata) );
   if( branchdata == NULL )
   {
      SCIPdebugMessage("  origbranch with no branchdata created\n");
   }
   SCIP_CALL( SCIPaddConsNode(origprob, childnode, origcons, NULL) );

   /* add those constraints to the node that enforce the branching decision in the original problem */
   origbranchconss = GCGconsMasterbranchGetOrigbranchConss(masterbranchchildcons);
   norigbranchconss = GCGconsMasterbranchGetNOrigbranchConss(masterbranchchildcons);
   for( i = 0; i < norigbranchconss; ++i )
   {
      SCIP_CALL( SCIPaddConsNode(origprob, childnode, origbranchconss[i], NULL) );
      SCIPdebugMessage("  add cons %s to node\n", SCIPconsGetName(origbranchconss[i]));
   }

   /* notify the original and master branching constraint about each other */
   GCGconsOrigbranchSetMastercons(origcons, masterbranchchildcons);
   GCGconsMasterbranchSetOrigcons(masterbranchchildcons, origcons);
   SCIPdebugMessage("  link branching conss %s <-> %s\n", SCIPconsGetName(masterbranchchildcons), SCIPconsGetName(origcons));

   SCIP_CALL( SCIPreleaseCons(origprob, &origcons) );

   /* release array of original branching constraints */
   SCIP_CALL( GCGconsMasterbranchReleaseOrigbranchConss(gcg, masterbranchchildcons) );

   return SCIP_OKAY;
}

/* apply a branching decision on the original variables to the corresponding node */
static
SCIP_RETCODE applyOriginalBranching(
   GCG*                  gcg,
   SCIP_NODE*            childnode,
   SCIP_CONS*            masterbranchchildcons
   )
{
   GCG_BRANCHDATA* branchdata;
   SCIP_VAR* boundvar;
   GCG_BOUNDTYPE boundtype;
   SCIP_Real newbound;
   SCIP* origprob;
   SCIP* masterprob;

   /* get branching decision */
   branchdata = GCGconsMasterbranchGetBranchdata(masterbranchchildcons);
   assert(branchdata != NULL);
   boundvar = GCGbranchOrigGetOrigvar(branchdata);
   boundtype = GCGbranchOrigGetBoundtype(branchdata);
   newbound = GCGbranchOrigGetNewbound(branchdata);

   masterprob = GCGgetMasterprob(gcg);
   origprob = GCGgetOrigprob(gcg);

   assert(boundvar != NULL);
   assert(boundtype == GCG_BOUNDTYPE_LOWER || boundtype == GCG_BOUNDTYPE_UPPER || boundtype == GCG_BOUNDTYPE_FIXED);
   assert(SCIPgetStage(masterprob) <= SCIP_STAGE_SOLVING);

   if( boundtype == GCG_BOUNDTYPE_LOWER || boundtype == GCG_BOUNDTYPE_FIXED )
   {
      if( SCIPisLE(origprob, newbound, SCIPvarGetUbLocal(boundvar)) )
      {
         if( SCIPisGT(origprob, newbound, SCIPvarGetLbLocal(boundvar)) )
            SCIP_CALL( SCIPchgVarLbNode(origprob, childnode, boundvar, newbound) );
      }
      else
      {
         // cut off child nodes
         SCIP_NODE* masterchildnode = GCGconsMasterbranchGetNode(masterbranchchildcons);
         SCIPupdateNodeLowerbound(masterprob, masterchildnode, SCIPinfinity(masterprob));
         SCIPupdateNodeLowerbound(origprob, childnode, SCIPinfinity(origprob));
      }
   }

   if( boundtype == GCG_BOUNDTYPE_UPPER || boundtype == GCG_BOUNDTYPE_FIXED )
   {
      if( SCIPisGE(origprob, newbound, SCIPvarGetLbLocal(boundvar)) )
      {
         if( SCIPisLT(origprob, newbound, SCIPvarGetUbLocal(boundvar)) )
            SCIP_CALL( SCIPchgVarUbNode(origprob, childnode, boundvar, newbound) );
      }
      else
      {
         // cut off child nodes
         SCIP_NODE* masterchildnode = GCGconsMasterbranchGetNode(masterbranchchildcons);
         SCIPupdateNodeLowerbound(masterprob, masterchildnode, SCIPinfinity(masterprob));
         SCIPupdateNodeLowerbound(origprob, childnode, SCIPinfinity(origprob));
      }
   }

   if( GCGvarGetBlock(boundvar) == -1 )
   {
      SCIP_CALL( GCGconsMasterbranchAddCopiedVarBndchg(gcg, masterbranchchildcons, boundvar, boundtype, newbound) );
   }

   return SCIP_OKAY;
}

/** creates branch-and-bound nodes in the original problem corresponding to those in the master problem */
static
SCIP_RETCODE createBranchNodesInOrigprob(
   GCG*                  gcg,                /**< SCIP data structure */
   SCIP_RESULT*          result              /**< result pointer */
)
{
   SCIP* origprob;
   SCIP* masterprob;
   SCIP_BRANCHRULE* branchrule;
   SCIP_CONS* masterbranchcons;
   int nchildnodes;
   SCIP_Bool enforcebycons;

   int i;

   assert(gcg != NULL);
   origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get master problem */
   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);

   if( SCIPgetStage(masterprob) > SCIP_STAGE_SOLVING)
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   GCGrestoreLimitSettings(gcg);

   /* get masterbranch constraint at the current node */
   masterbranchcons = GCGconsMasterbranchGetActiveCons(gcg);

   /* @todo: Why should this happen? */
   if( masterbranchcons == NULL )
      return SCIP_OKAY;

   /* get the children of the current node */
   nchildnodes = GCGconsMasterbranchGetNChildconss(masterbranchcons);

   /* check if the focus node of the master problem has children */
   if( nchildnodes <= 0 && SCIPgetStage(masterprob) != SCIP_STAGE_SOLVED && SCIPgetNChildren(masterprob) >= 1 )
   {
      SCIP_NODE* child;

      SCIPdebugMessage("create dummy child in origprob, because there is also a child in the master\n");

      /* create dummy child */
      SCIP_CALL( SCIPcreateChild(origprob, &child, 0.0, SCIPgetLocalTransEstimate(origprob)) );

      *result = SCIP_BRANCHED;
      return SCIP_OKAY;
   }

   if( nchildnodes <= 0 )
   {
      SCIPdebugMessage("node cut off, since there is no successor node\n");

      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPgetBoolParam(origprob, "branching/orig/enforcebycons", &enforcebycons) );

   /* for each child, create a corresponding node in the original problem as well as an origbranch constraint */
   for( i = 0; i < nchildnodes; ++i )
   {
      SCIP_NODE* childnode;
      SCIP_CONS* masterbranchchildcons = GCGconsMasterbranchGetChildcons(masterbranchcons, i);
      assert(masterbranchchildcons != NULL);

      /* create a child node and an origbranch constraint holding the branching decision */
      SCIP_CALL( SCIPcreateChild(origprob, &childnode, 0.0, SCIPgetLocalTransEstimate(origprob)) );
      SCIP_CALL( createOrigbranchConstraint(gcg, childnode, masterbranchchildcons) );

      /* get branching rule */
      branchrule = GCGconsMasterbranchGetBranchrule(masterbranchchildcons);

      /* If a branching decision on an original variable was made, apply it */
      if( !enforcebycons && branchrule != NULL && strcmp(SCIPbranchruleGetName(branchrule), "orig") == 0 )
      {
         SCIP_CALL( applyOriginalBranching(gcg, childnode, masterbranchchildcons) );
      }

      /* @fixme: this should actually be an assertion */
      if( SCIPnodeGetNumber(GCGconsOrigbranchGetNode(GCGconsOrigbranchGetActiveCons(gcg))) != SCIPnodeGetNumber(GCGconsMasterbranchGetNode(GCGconsMasterbranchGetActiveCons(gcg))) )
      {
   #ifdef SCIP_DEBUG
         SCIPwarningMessage(origprob, "norignodes = %lld; nmasternodes = %lld\n",
            SCIPnodeGetNumber(GCGconsOrigbranchGetNode(GCGconsOrigbranchGetActiveCons(gcg))),
            SCIPnodeGetNumber(GCGconsMasterbranchGetNode(GCGconsMasterbranchGetActiveCons(gcg))));
   #endif
      }
   }

   *result = SCIP_BRANCHED;

   assert(nchildnodes > 0);

   return SCIP_OKAY;
}


/** copy method for empty branching rule */
static
SCIP_DECL_BRANCHCOPY(branchCopyEmpty)
{
   assert(branchrule != NULL);
   assert(scip != NULL);

   /* SubSCIPs are solved with SCIP rather than GCG;
    * therefore, only the default SCIP branching rules are included into the subSCIP.
    */
   SCIP_CALL( includeSCIPBranchingRules(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeEmpty)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPfreeMemory(scip, &branchruledata);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
#define branchInitEmpty NULL

/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitEmpty NULL

/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolEmpty NULL

/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolEmpty NULL

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpEmpty)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->gcg != NULL);

   SCIP_CALL( createBranchNodesInOrigprob(branchruledata->gcg, result) );

   return SCIP_OKAY;
}

/** branching execution method relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextEmpty)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->gcg != NULL);

   SCIP_CALL( createBranchNodesInOrigprob(branchruledata->gcg, result) );

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsEmpty)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->gcg != NULL);

   SCIP_CALL( createBranchNodesInOrigprob(branchruledata->gcg, result) );

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the empty branching rule and includes it in SCIP */
SCIP_RETCODE GCGincludeBranchruleEmpty(
   GCG*                 gcg                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   /* create inference branching rule data */
   SCIP_CALL( SCIPallocMemory(origprob, &branchruledata) );
   branchruledata->gcg = gcg;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(origprob, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
      BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
      branchCopyEmpty, branchFreeEmpty, branchInitEmpty, branchExitEmpty, branchInitsolEmpty,
      branchExitsolEmpty, branchExeclpEmpty, branchExecextEmpty, branchExecpsEmpty,
      branchruledata) );

   return SCIP_OKAY;
}
