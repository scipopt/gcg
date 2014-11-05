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

/**@file   branch_empty.c
 * @brief  branching rule for original problem in GCG while real branching is applied in the master
 * @author Marcel Schmickerath
 * @author Martin Bergner
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/
#include <assert.h>
#include <string.h>

#include "branch_empty.h"
#include "relax_gcg.h"
#include "gcg.h"
#include "cons_masterbranch.h"
#include "cons_origbranch.h"
#include "scip/branch_allfullstrong.h"
#include "scip/branch_fullstrong.h"
#include "scip/branch_inference.h"
#include "scip/branch_mostinf.h"
#include "scip/branch_leastinf.h"
#include "scip/branch_pscost.h"
#include "scip/branch_random.h"
#include "scip/branch_relpscost.h"
#include "pricer_gcg.h"
#include "scip/cons_varbound.h"
#include "type_branchgcg.h"

#define BRANCHRULE_NAME          "empty"
#define BRANCHRULE_DESC          "empty branching in generic column generation"
#define BRANCHRULE_PRIORITY      1000000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0


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

/** for a new branch-and-bound node on the master problem, create a corresponding node in the original problem
 *  as well as an origbranch constraint that holds the branching decision
 */
SCIP_RETCODE GCGcreateConsOrigbranchNode(
   SCIP*                 scip,
   SCIP_CONS*            masterbranchchildcons
)
{
   SCIP_NODE* child;
   SCIP_CONS*  origbranch;
   SCIP_CONS** origbranchconss;
   int norigbranchconss;
   SCIP_Bool enforcebycons;

   int i;

   assert(scip != NULL);
   assert(masterbranchchildcons != NULL);

   /* create a child node and an origbranch constraint holding the branching decision */
   SCIP_CALL( SCIPcreateChild(scip, &child, 0.0, SCIPgetLocalTransEstimate(scip)) );

   SCIPdebugMessage("Create original branching constraint %s\n", GCGconsMasterbranchGetOrigbranchConsName(masterbranchchildcons));

   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranch, GCGconsMasterbranchGetOrigbranchConsName(masterbranchchildcons), child,
            GCGconsOrigbranchGetActiveCons(scip), GCGconsMasterbranchGetOrigbranchrule(masterbranchchildcons), GCGconsMasterbranchGetOrigbranchdata(masterbranchchildcons)) );

   if( GCGconsMasterbranchGetOrigbranchdata(masterbranchchildcons) == NULL )
   {
      SCIPdebugMessage("origbranch with no branchdata created\n");
   }

   SCIP_CALL( SCIPaddConsNode(scip, child, origbranch, NULL) );

   /* get previous original branching constraints and add them to the new node */
   origbranchconss = GCGconsMasterbranchGetOrigbranchConss(masterbranchchildcons);
   norigbranchconss = GCGconsMasterbranchGetNOrigbranchConss(masterbranchchildcons);
   for( i = 0; i < norigbranchconss; ++i )
   {
      SCIP_CALL( SCIPaddConsNode(scip, child, origbranchconss[i], NULL) );
   }

   /* If a branching decision on an original variable was made, apply it */
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/enforcebycons", &enforcebycons) );
   if( !enforcebycons )
   {
      if( GCGconsMasterbranchGetOrigboundtype(masterbranchchildcons) == GCG_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPchgVarLbNode(scip, child,
            GCGconsMasterbranchGetOrigboundvar(masterbranchchildcons),
            GCGconsMasterbranchGetOrigbound(masterbranchchildcons)) );
      }
      else if( GCGconsMasterbranchGetOrigboundtype(masterbranchchildcons) == GCG_BOUNDTYPE_UPPER )
      {
         SCIP_CALL( SCIPchgVarUbNode(scip, child,
            GCGconsMasterbranchGetOrigboundvar(masterbranchchildcons),
            GCGconsMasterbranchGetOrigbound(masterbranchchildcons)) );
      }
   }

   if( GCGconsMasterbranchGetPropagatebndchg(masterbranchchildcons) )
   {
      assert(GCGconsMasterbranchGetOrigboundtype(masterbranchchildcons) != GCG_BOUNDTYPE_NONE);

      SCIP_CALL( GCGconsOrigbranchAddPropBoundChg(scip, origbranch,
         GCGconsMasterbranchGetOrigboundvar(masterbranchchildcons),
         GCGconsMasterbranchGetOrigboundtype(masterbranchchildcons),
         GCGconsMasterbranchGetOrigbound(masterbranchchildcons)) );
   }

   GCGconsOrigbranchSetMastercons(origbranch, masterbranchchildcons);
   GCGconsMasterbranchSetOrigcons(masterbranchchildcons, origbranch);

   SCIP_CALL( SCIPreleaseCons(scip, &origbranch) );

   /* release array of original branching constraints */
   SCIP_CALL( GCGconsMasterbranchReleaseOrigbranchConss(GCGgetMasterprob(scip), scip, masterbranchchildcons) );

   /* @fixme: this should actually be an assertion */
   if( SCIPnodeGetNumber(GCGconsOrigbranchGetNode(GCGconsOrigbranchGetActiveCons(scip))) != SCIPnodeGetNumber(GCGconsMasterbranchGetNode(GCGconsMasterbranchGetActiveCons(GCGgetMasterprob(scip)))) )
   {
#ifdef SCIP_DEBUG
      SCIPwarningMessage(scip, "norignodes = %lld; nmasternodes = %lld\n",
         SCIPnodeGetNumber(GCGconsOrigbranchGetNode(GCGconsOrigbranchGetActiveCons(scip))),
         SCIPnodeGetNumber(GCGconsMasterbranchGetNode(GCGconsMasterbranchGetActiveCons(GCGgetMasterprob(scip)))));
#endif
   }

   assert(SCIPgetNNodes(scip) == SCIPgetNNodes(GCGgetMasterprob(scip)));
   /*assert(SCIPnodeGetNumber(GCGconsOrigbranchGetNode(GCGconsOrigbranchGetActiveCons(scip))) == SCIPnodeGetNumber(GCGconsMasterbranchGetNode(GCGconsMasterbranchGetActiveCons(GCGrelaxGetMasterprob(scip)))));*/

   return SCIP_OKAY;
}

/** creates branch-and-bound nodes in the original problem corresponding to those in the master problem */
static
SCIP_RETCODE createBranchNodesInOrigprob(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< result pointer */
)
{
   SCIP* masterscip;
   SCIP_Bool feasible;
   SCIP_CONS* masterbranchcons;
   int nchildnodes;

   int i;

   feasible = TRUE;

   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* @fixme: This looks like a hack */
   if( GCGrelaxGetCurrentOrigSol(scip) == NULL )
   {
      SCIP_CALL( GCGrelaxUpdateCurrentSol(scip, &feasible) );
   }
   else
   {
      /* check whether the current original solution is integral */
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), TRUE, TRUE, TRUE, TRUE, &feasible) );
#else
      SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), FALSE, TRUE, TRUE, TRUE, &feasible) );
#endif
   }

   if( feasible )
   {
      SCIPdebugMessage("node cut off, since origsol was feasible, solval = %f\n",
            SCIPgetSolOrigObj(scip, GCGrelaxGetCurrentOrigSol(scip)));

      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* get master problem */
   masterscip = GCGgetMasterprob(scip);
   assert(masterscip != NULL);

   /* get masterbranch constraint at the current node */
   masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);
   /* @todo: Why should this happen? */
   if( masterbranchcons == NULL )
      return SCIP_OKAY;

   /* get the children of the current node */
   nchildnodes = GCGconsMasterbranchGetNChildcons(masterbranchcons);
   if( nchildnodes <= 0 )
   {
      SCIPdebugMessage("node cut off, since there is no successor node\n");

      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* for each child, create a corresponding node in the original problem as well as an origbranch constraint */
   for( i = 0; i < nchildnodes; ++i )
   {
      SCIP_CONS* masterbranchchildcons = GCGconsMasterbranchGetChildcons(masterbranchcons, i);
      assert(masterbranchchildcons != NULL);

      SCIP_CALL( GCGcreateConsOrigbranchNode(scip, masterbranchchildcons) );
   }

   *result = SCIP_BRANCHED;

   assert(nchildnodes > 0);

   return SCIP_OKAY;
}


/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeEmpty NULL

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
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIP_CALL( createBranchNodesInOrigprob(scip, result) );

   return SCIP_OKAY;
}

/** branching execution method relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextEmpty)
{  /*lint --e{715}*/
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIP_CALL( createBranchNodesInOrigprob(scip, result) );

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsEmpty)
{  /*lint --e{715}*/
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIP_CALL( createBranchNodesInOrigprob(scip, result) );

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the empty branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleEmpty(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create inference branching rule data */
   branchruledata = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
      BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
      branchCopyEmpty, branchFreeEmpty, branchInitEmpty, branchExitEmpty, branchInitsolEmpty,
      branchExitsolEmpty, branchExeclpEmpty, branchExecextEmpty, branchExecpsEmpty,
      branchruledata) );

   return SCIP_OKAY;
}
