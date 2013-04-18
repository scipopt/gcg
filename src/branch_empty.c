/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
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
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_empty.h"
#include "relax_gcg.h"
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
#include "pub_gcgvar.h"

#define BRANCHRULE_NAME          "empty"
#define BRANCHRULE_DESC          "empty branching in generic column generation"
#define BRANCHRULE_PRIORITY      1000000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0


/*
 * Callback methods for enforcing branching constraints
 */

/* copy default SCIP branching rules to allow solving restrictions of the original problem as a subSCIP without
 * Dantzig-Wolfe decomposition
 */
static
SCIP_RETCODE GCGincludeOriginalCopyPlugins(
   SCIP* scip
   )
{
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
/** copy method for master branching rule */
static
SCIP_DECL_BRANCHCOPY(branchCopyEmpty)
{
   assert(branchrule != NULL);
   assert(scip != NULL);
   SCIPdebugMessage("pricer copy called.\n");
   SCIP_CALL( GCGincludeOriginalCopyPlugins(scip) );
   return SCIP_OKAY;
}

SCIP_RETCODE GCGcreateConsOrigbranchNode(
      SCIP* scip,
      SCIP_CONS* masterbranchchildcons
   )
{
   SCIP_NODE* child;
   SCIP_CONS*  origbranch;
   SCIP_CONS** origbranchcons;
   int norigbranchcons;
   int i;

   i = 0;

   assert(masterbranchchildcons != NULL);

   SCIP_CALL( SCIPcreateChild(scip, &child, 0.0, SCIPgetLocalTransEstimate(scip)) );

   SCIPdebugMessage("Name is %s\n", GCGconsMasterbranchGetOrigbranchConsName(masterbranchchildcons) );

   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranch, GCGconsMasterbranchGetOrigbranchConsName(masterbranchchildcons), child,
            GCGconsOrigbranchGetActiveCons(scip), GCGconsMasterbranchGetOrigbranchrule(masterbranchchildcons), GCGconsMasterbranchGetOrigbranchdata(masterbranchchildcons)) );

   if( GCGconsMasterbranchGetOrigbranchdata(masterbranchchildcons) == NULL )
      SCIPdebugMessage("origbranch with no branchdata created\n");

   SCIP_CALL( SCIPaddConsNode(scip, child, origbranch, NULL) );


   norigbranchcons = GCGconsMasterbranchGetNOrigbranchCons(masterbranchchildcons);
   origbranchcons = GCGconsMasterbranchGetOrigbranchCons(masterbranchchildcons);

   for( i=0; i<norigbranchcons; ++i)
   {
      SCIP_CALL( SCIPaddConsNode(scip, child, origbranchcons[i], NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &(origbranchcons[i])) );
   }

   if( GCGconsMasterbranchGetOrigbranchConsChgVarUbNode(masterbranchchildcons) )
   {
      SCIP_CALL( SCIPchgVarUbNode(scip, child,
         GCGconsMasterbranchGetOrigbranchConsChgVarNodeVar(masterbranchchildcons),
         GCGconsMasterbranchGetOrigbranchConsChgVarNodeBound(masterbranchchildcons)) );
   }
   if( GCGconsMasterbranchGetOrigbranchConsChgVarLbNode(masterbranchchildcons) )
   {
      SCIP_CALL( SCIPchgVarUbNode(scip, child,
         GCGconsMasterbranchGetOrigbranchConsChgVarNodeVar(masterbranchchildcons),
         GCGconsMasterbranchGetOrigbranchConsChgVarNodeBound(masterbranchchildcons)) );
   }
   if( GCGconsMasterbranchGetOrigbranchConsAddPropBoundChg(masterbranchchildcons) )
   {
      SCIP_CALL( GCGconsOrigbranchAddPropBoundChg(scip, origbranch,
            GCGconsMasterbranchGetOrigbranchConsChgVarNodeVar(masterbranchchildcons),
            GCGconsMasterbranchGetOrigbranchConsAddPropBoundChgBoundtype(masterbranchchildcons),
            GCGconsMasterbranchGetOrigbranchConsAddPropBoundChgBound(masterbranchchildcons)) );
   }
   SCIP_CALL( SCIPreleaseCons(scip, &origbranch) );

   if( norigbranchcons > 0 )
      SCIPfreeMemoryArray(GCGrelaxGetMasterprob(scip), &origbranchcons);

/*   SCIP_CALL( GCGconsMasterbranchSetOrigConsData(GCGrelaxGetMasterprob(scip), masterbranchchildcons, NULL, NULL,
      GCGconsMasterbranchGetOrigbranchdata(masterbranchchildcons), NULL, 0, FALSE, FALSE, FALSE, NULL, 0, NULL, 0) );*/


   if(SCIPnodeGetNumber(GCGconsOrigbranchGetNode(GCGconsOrigbranchGetActiveCons(scip))) != SCIPnodeGetNumber(GCGconsMasterbranchGetNode(GCGconsMasterbranchGetActiveCons(GCGrelaxGetMasterprob(scip)))))
   {
      SCIPdebugMessage("norignodes = %d; nmasternodes = %d\n",SCIPnodeGetNumber(GCGconsOrigbranchGetNode(GCGconsOrigbranchGetActiveCons(scip))), SCIPnodeGetNumber(GCGconsMasterbranchGetNode(GCGconsMasterbranchGetActiveCons(GCGrelaxGetMasterprob(scip)))));
   }

   assert(SCIPgetNNodes(scip) == SCIPgetNNodes(GCGrelaxGetMasterprob(scip)));
   /*assert(SCIPnodeGetNumber(GCGconsOrigbranchGetNode(GCGconsOrigbranchGetActiveCons(scip))) == SCIPnodeGetNumber(GCGconsMasterbranchGetNode(GCGconsMasterbranchGetActiveCons(GCGrelaxGetMasterprob(scip)))));*/

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

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}
/** branching execution method relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextEmpty)
{  /*lint --e{715}*/
   SCIP* masterscip;
   SCIP_Bool feasible;
   SCIP_CONS* masterbranchcons;
   SCIP_CONS* masterbranchchildcons;
   int nchildnodes;
   int i;

   nchildnodes = 0;
   i = 0;
   feasible = TRUE;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* check whether the current original solution is integral */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), TRUE, TRUE, TRUE, TRUE, &feasible) );
#else
   SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), FALSE, TRUE, TRUE, TRUE, &feasible) );
#endif

   if( feasible )
   {
      SCIPdebugMessage("node cut off, since origsol was feasible, solval = %f\n",
         SCIPgetSolOrigObj(scip, GCGrelaxGetCurrentOrigSol(scip)));

      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Execext method of empty branching\n");

   masterscip = GCGrelaxGetMasterprob(scip);
   assert(masterscip != NULL);

   masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);
   assert(masterbranchcons != NULL);

   nchildnodes = GCGconsMasterbranchGetNChildcons(masterbranchcons);
   if( nchildnodes <= 0 )
   {
      SCIPdebugMessage("node cut off, since there is no successor node\n");

      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   for( i=0; i<nchildnodes; ++i )
   {
      masterbranchchildcons = GCGconsMasterbranchGetChildcons(masterbranchcons, i);
      assert(masterbranchchildcons != NULL);

      SCIP_CALL( GCGcreateConsOrigbranchNode(scip, masterbranchchildcons));
   }

   *result = SCIP_BRANCHED;
   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsEmpty)
{  /*lint --e{715}*/
   SCIP* masterscip;
   SCIP_Bool feasible;
   SCIP_CONS* masterbranchcons;
   SCIP_CONS* masterbranchchildcons;
   int nchildnodes;
   int i;

   nchildnodes = 0;
   i = 0;
   feasible = TRUE;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* check whether the current original solution is integral */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), TRUE, TRUE, TRUE, TRUE, &feasible) );
#else
   SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), FALSE, TRUE, TRUE, TRUE, &feasible) );
#endif

   if( feasible )
   {
      SCIPdebugMessage("node cut off, since origsol was feasible, solval = %f\n",
         SCIPgetSolOrigObj(scip, GCGrelaxGetCurrentOrigSol(scip)));

      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Execeps method of empty branching\n");

   masterscip = GCGrelaxGetMasterprob(scip);
      assert(masterscip != NULL);

      masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);
      assert(masterbranchcons != NULL);

      nchildnodes = GCGconsMasterbranchGetNChildcons(masterbranchcons);
      if( nchildnodes <= 0 )
      {
         SCIPdebugMessage("node cut off, since there is no successor node\n");

         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      for( i=0; i<nchildnodes; ++i )
      {
         masterbranchchildcons = GCGconsMasterbranchGetChildcons(masterbranchcons, i);
         assert(masterbranchchildcons != NULL);

         SCIP_CALL( GCGcreateConsOrigbranchNode(scip, masterbranchchildcons));
      }

      *result = SCIP_BRANCHED;
      return SCIP_OKAY;
}

/*
 * branching specific interface methods
 */
/** creates the master branching rule and includes it in SCIP */
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
