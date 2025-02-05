/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file    branch_staticvar.c
 *
 * @brief   static master variable branching rule
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <scip/cons_linear.h>
#include <scip/def.h>
#include <scip/scip.h>
#include <scip/type_retcode.h>
#include <scip/type_var.h>
#include <string.h>

#include "branch_staticvar.h"
#include "cons_masterbranch.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "type_branchgcg.h"


#define BRANCHRULE_NAME        "staticvar"                       /**< name of branching rule */
#define BRANCHRULE_DESC        "static mastervariable branching" /**< short description of branching rule */
#define BRANCHRULE_PRIORITY        -100000                             /**< priority of this branching rule */
#define BRANCHRULE_MAXDEPTH        -1                            /**< maximal depth level of the branching rule */
#define BRANCHRULE_MAXBOUNDDIST    1.0                           /**< maximal relative distance from current node's
                                                                      dual bound to primal bound compared to best node's
                                                                      dual bound for applying branching */


/*
 * Data structures
 */
typedef enum {
   DownBranch = 0,
   UpBranch = 1
} BranchType;

struct GCG_BranchData
{
   SCIP_VAR*             mastervar;          /**< mastervariable to branch on **/
   BranchType            branchtype;         /**< type of branch **/
   SCIP_Real             bound;              /**< chosen bound of the mastervariable **/
};

/*
 * Local methods
 */

static
SCIP_RETCODE createChildNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_VAR*             mastervar,          /**< mastervariable to branch on */
   SCIP_Real             bound               /**< chosen bound of the mastervariable */
   )
{
   SCIP* masterscip;
   GCG_BRANCHDATA* downbranchdata;
   GCG_BRANCHDATA* upbranchdata;
   SCIP_NODE* downchild;
   SCIP_NODE* upchild;
   SCIP_CONS* downcons;
   SCIP_CONS* upcons;
   char downname[SCIP_MAXSTRLEN];
   char upname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(mastervar != NULL);
   assert(GCGvarGetBlock(mastervar) == -1);
   assert(SCIPvarGetType(mastervar) <= SCIP_VARTYPE_INTEGER);
   assert(!SCIPisFeasIntegral(scip, SCIPvarGetLPSol(mastervar)));

   masterscip = GCGmasterGetOrigprob(scip);
   assert(masterscip != NULL);

   SCIPdebugMessage("createChildNodes: mastervar = %s, bound = %.2f\n", SCIPvarGetName(mastervar), bound);

   /* create down branch */
   SCIP_CALL( SCIPallocBlockMemory(scip, &downbranchdata) );
   downbranchdata->mastervar = mastervar;
   downbranchdata->branchtype = DownBranch;
   downbranchdata->bound = SCIPfloor(masterscip, bound);
   SCIP_CALL( SCIPcreateChild(masterscip, &downchild, 0.0, SCIPgetLocalTransEstimate(masterscip)));
   (void) SCIPsnprintf(downname, SCIP_MAXSTRLEN, "down(%s,%.2f)", SCIPvarGetName(downbranchdata->mastervar), downbranchdata->bound);
   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &downcons, downname, downchild,
      GCGconsMasterbranchGetActiveCons(masterscip), branchrule, downbranchdata, NULL, 0, 0) );
   SCIP_CALL( SCIPaddConsNode(masterscip, downchild, downcons, NULL) );
   SCIP_CALL( SCIPcreateConsLinear(scip, &downcons, downname, 0, NULL, NULL, -SCIPinfinity(scip), downbranchdata->bound, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddCoefLinear(scip, downcons, downbranchdata->mastervar, 1.0) );
   SCIP_CALL( SCIPreleaseCons(masterscip, &downcons) );

   /* create up branch */
   SCIP_CALL( SCIPallocBlockMemory(scip, &upbranchdata) );
   upbranchdata->mastervar = mastervar;
   upbranchdata->branchtype = UpBranch;
   upbranchdata->bound = SCIPceil(masterscip, bound);
   SCIP_CALL( SCIPcreateChild(masterscip, &upchild, 0.0, SCIPgetLocalTransEstimate(masterscip)));
   (void) SCIPsnprintf(upname, SCIP_MAXSTRLEN, "up(%s,%.2f)", SCIPvarGetName(upbranchdata->mastervar), upbranchdata->bound);
   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &upcons, upname, upchild,
      GCGconsMasterbranchGetActiveCons(masterscip), branchrule, upbranchdata, NULL, 0, 0) );
   SCIP_CALL( SCIPaddConsNode(masterscip, upchild, upcons, NULL) );
   SCIP_CALL( SCIPcreateConsLinear(scip, &upcons, upname, 0, NULL, NULL, upbranchdata->bound, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddCoefLinear(scip, upcons, upbranchdata->mastervar, 1.0) );
   SCIP_CALL( SCIPreleaseCons(masterscip, &upcons) );

   return SCIP_OKAY;
}


/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
#define branchCopyStaticVar NULL

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeStaticVar NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitStaticVar NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolStaticVar NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolStaticVar NULL


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpStaticVar)
{
   SCIP* origscip;
   SCIP_VAR** branchcands;
   int nbranchcands;
   int i;
   SCIP_VAR* chosenvar;
   SCIP_Real solval;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   SCIPdebugMessage("Execlp method of static mastervariable branching\n");

   *result = SCIP_DIDNOTRUN;

   if( GCGcurrentNodeIsGeneric(scip) )
   {
      SCIPdebugMessage("Not executing static mastervar branching, node was branched by generic branchrule\n");
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   if( GCGrelaxIsOrigSolFeasible(origscip) )
   {
      SCIPdebugMessage("node cut off, since origsol was feasible, solval = %f\n",
         SCIPgetSolOrigObj(origscip, GCGrelaxGetCurrentOrigSol(origscip)));
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPgetLPBranchCands(scip, &branchcands, NULL, NULL, &nbranchcands, NULL, NULL) );
   chosenvar = NULL;
   for( i=0; i < nbranchcands; i++ ) {
      SCIP_VAR* var = branchcands[i];

      // var must be static variable (blocknr == -1)
      if( GCGvarGetBlock(var) != -1 )
         continue;

      // var must be integral
      if( SCIPvarGetType(var) > SCIP_VARTYPE_INTEGER )
         continue;

      // solution of var must be fractional
      solval = SCIPvarGetLPSol(var);
      if( SCIPisFeasIntegral(scip, solval) )
         continue;

      chosenvar = var;
      break;
   }

   if( chosenvar == NULL )
   {
      SCIPdebugMessage("No fractional static variable found\n");
      return SCIP_OKAY;
   }

   // branch on chosenvar <= floor(solval), chosenvar >= ceil(solval)
   SCIP_CALL( createChildNodes(scip, branchrule, chosenvar, solval) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextStaticVar)
{
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsStaticVar)
{
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}

/*
 * GCG specific branching rule callbacks
 */

/** activation method for branchrule, called when a node in the master problem is activated,
 *  should perform changes to the current node's problem due to the branchdata
 */
#define branchActiveMasterStaticVar NULL


/** deactivation method for branchrule, called when a node in the master problem is deactivated,
 *  should undo changes to the current node's problem due to the branchdata
 */
#define branchDeactiveMasterStaticVar NULL

/** propagation method for branchrule, called when a node in the master problem is propagated,
 *  should perform propagation at the current node due to the branchdata
 */
 #define branchPropMasterStaticVar NULL

/** method for branchrule, called when the master LP is solved at one node,
 *  can store pseudocosts for the branching decisions
 */
#define branchMasterSolvedStaticVar NULL

/** frees branching data of an origbranch constraint (called when the origbranch constraint is deleted) */
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteStaticVar)
{
   assert(scip != NULL);
   assert(branchdata != NULL);

   SCIPdebugMessage("branchDataDeleteStaticVar: (%s %s %.2f)\n", SCIPvarGetName((*branchdata)->mastervar), ((*branchdata)->branchtype == DownBranch ? "<=" : ">="), (*branchdata)->bound);

   SCIPfreeBlockMemory(scip, branchdata);
   *branchdata = NULL;

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitStaticVar)
{
   SCIP* origprob;

   origprob = GCGmasterGetOrigprob(scip);
   assert(branchrule != NULL);
   assert(origprob != NULL);

   SCIP_CALL( GCGrelaxIncludeBranchrule(origprob, branchrule, branchActiveMasterStaticVar,
         branchDeactiveMasterStaticVar, branchPropMasterStaticVar, branchMasterSolvedStaticVar, branchDataDeleteStaticVar,
         NULL, NULL) );

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the xyz branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleStaticVar(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
         BRANCHRULE_MAXBOUNDDIST,
         branchCopyStaticVar, branchFreeStaticVar, branchInitStaticVar, branchExitStaticVar, branchInitsolStaticVar, branchExitsolStaticVar,
         branchExeclpStaticVar, branchExecextStaticVar, branchExecpsStaticVar,
         branchruledata) );

   return SCIP_OKAY;
}
