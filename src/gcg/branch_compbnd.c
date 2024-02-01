/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       */
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

/**@file    branch_compbnd.c
 *
 * @brief   branching rule based on vanderbeck's component bound branching
 * @author  Til Mohr
 * @author Marcel Schmickerath
 * @author Martin Bergner
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <scip/def.h>
#include <scip/scip.h>
#include <string.h>

#include "branch_compbnd.h"
#include "cons_integralorig.h"
#include "cons_masterbranch.h"
#include "gcg.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "type_branchgcg.h"

#include "scip/branch_allfullstrong.h"
#include "scip/branch_fullstrong.h"
#include "scip/branch_inference.h"
#include "scip/branch_mostinf.h"
#include "scip/branch_leastinf.h"
#include "scip/branch_pscost.h"
#include "scip/branch_random.h"
#include "scip/branch_relpscost.h"
#include "scip/cons_linear.h"
#include "scip/nodesel_bfs.h"
#include "scip/nodesel_dfs.h"
#include "scip/nodesel_estimate.h"
#include "scip/nodesel_hybridestim.h"
#include "scip/nodesel_restartdfs.h"


#define BRANCHRULE_NAME            "compbnd"                      /**< name of branching rule */
#define BRANCHRULE_DESC            "component bound branching"    /**< short description of branching rule */
#define BRANCHRULE_PRIORITY        0                              /**< priority of this branching rule */
#define BRANCHRULE_MAXDEPTH        -1                             /**< maximal depth level of the branching rule */
#define BRANCHRULE_MAXBOUNDDIST    1.0                            /**< maximal relative distance from current node's
                                                                   dual bound to primal bound compared to best node's
                                                                   dual bound for applying branching */

#define EVENTHDLR_NAME         "compbndbranchvaradd"
#define EVENTHDLR_DESC         "event handler for adding a new generated mastervar into the right branching constraints by using component bound branching"


/*
 * Data structures
 */

/** branching data */
struct GCG_BranchData
{
   GCG_BRANCH_TYPE       branchtype;         /**< type of branching */
   SCIP_CONS*            mastercons;         /**< constraint enforcing the branching restriction in the master problem */
   SCIP_Real             constant;           /**< constant value of the branching constraint in the master problem - either lhs or rhs, depending on branchtype */
   GCG_COMPSEQUENCE*     B;                  /**< component bound sequence which induce the current branching constraint */
   int                   Bsize;              /**< size of the component bound sequence B */
   int                   consblocknr;        /**< id of the pricing problem (or block) to which this branching constraint belongs */
   int                   nvars;              /**< number of master variables the last time the node has been visited - neccessary to later include newly generated master variables */
};

/** set of component bounds in separate */
struct GCG_Record
{
   GCG_COMPSEQUENCE**   record;              /**< array of component bound sequences in record */
   int                  recordsize;          /**< number of component bound sequences in record */
   int                  recordcapacity;      /**< capacity of record */
   int*                 sequencesizes;       /**< array of sizes of component bound sequences */
   int*                 capacities;          /**< array of capacities of component bound sequences */
};
typedef struct GCG_Record GCG_RECORD;

/*
 * Local methods
 */

 /** initialize branchdata at the node */
static
SCIP_RETCODE initNodeBranchdata(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_BRANCHDATA**      nodebranchdata,     /**< branching data to set */
   GCG_BRANCH_TYPE       branchtype,         /**< type of branch to generate */
   SCIP_Real             constant,           /**< constant value of the branching constraint in the master problem - either lhs or rhs, depending on branchtype */
   int                   blocknr             /**< block we are branching in */
   )
{
   SCIP_CALL( SCIPallocBlockMemory(scip, nodebranchdata) );

   (*nodebranchdata)->branchtype = branchtype;
   (*nodebranchdata)->consblocknr = blocknr;
   (*nodebranchdata)->mastercons = NULL;
   (*nodebranchdata)->constant = constant;
   (*nodebranchdata)->B = NULL;
   (*nodebranchdata)->Bsize = 0;
   (*nodebranchdata)->nvars = 0;

   return SCIP_OKAY;
}

/** computes the generator of mastervar for the entry in origvar
 * @return entry of the generator corresponding to origvar */
static
SCIP_Real getGeneratorEntry(
   SCIP_VAR*             mastervar,          /**< current mastervariable */
   SCIP_VAR*             origvar             /**< corresponding origvar */
   )
{
   int i;
   SCIP_VAR** origvars;
   SCIP_Real* origvals;
   int norigvars;

   assert(mastervar != NULL);
   assert(origvar != NULL);

   origvars = GCGmasterVarGetOrigvars(mastervar);
   norigvars = GCGmasterVarGetNOrigvars(mastervar);
   origvals = GCGmasterVarGetOrigvals(mastervar);

   for( i = 0; i < norigvars; ++i )
   {
      if( origvars[i] == origvar )
      {
         return origvals[i];
      }
   }

   return 0.0;
}

/** whether a master variable is in B or not */
static
SCIP_Bool isMasterVarInB(
   SCIP_VAR*             mastervar,          /**< master variable to check */
   GCG_COMPSEQUENCE*     B,                  /**< component bound sequence to check */
   int                   Bsize               /**< size of B */
   )
{
   int i;

   assert(mastervar != NULL);
   assert(B != NULL);
   assert(Bsize > 0);

   for( i = 0; i < Bsize; ++i )
   {
      SCIP_Real generatorentry = getGeneratorEntry(mastervar, B[i].component);
      if ( (B[i].sense == GCG_COMPSENSE_GE && SCIPisLT(NULL, generatorentry, B[i].bound)) ||
           (B[i].sense == GCG_COMPSENSE_LE && SCIPisGT(NULL, generatorentry, B[i].bound)) )
      {
         return FALSE;
      }
   }

   return TRUE;
}

 /** for given component bound sequence S, create |S|+1 Vanderbeck branching nodes */
static
SCIP_RETCODE createChildNodesCompBnd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   GCG_COMPSEQUENCE*     B,                  /**< Component Bound Sequence defining the nodes */
   int                   Bsize,              /**< size of B */
   int                   blocknr,            /**< number of the block */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching call */
   )
{
   SCIP*  masterscip;
   int i;
   int p;
   int identicalBlocks;
   SCIP_Real constantSum;
   int nmastervars;
   int nbranchcands;
   int nchildnodes;
   SCIP_VAR** mastervars;
   SCIP_VAR** branchcands;
   GCG_BRANCHDATA* downBranchData;
   GCG_BRANCHDATA* upBranchData;

   assert(scip != NULL);
   assert(Bsize > 0);
   assert(B != NULL);

   nchildnodes = 0;

   identicalBlocks = GCGgetNIdenticalBlocks(scip, blocknr);
   SCIPdebugMessage("Component bound branching rule Node creation for blocknr %d with %d identical blocks \n", blocknr, identicalBlocks);


   /*  get variable data of the master problem */
   masterscip = GCGgetMasterprob(scip);
   assert(masterscip != NULL);
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(nmastervars >= 0);

   SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL, NULL) );

   /* determine the constant value of the master constraint, i.e. the rhs for the down branch, and lhs for the up branch */
   constantSum = 0;
   for( i = 0; i < nmastervars; ++i )
   {
      if( GCGisMasterVarInBlock(mastervars[i], blocknr) && isMasterVarInB(mastervars[i], B, Bsize) )
      {
         constantSum += SCIPgetSolVal(masterscip, NULL, mastervars[i]);
      }
   }
   // sanity check: the sum must be fractional, otherwise something went wrong during seperation, i.e. B is incorrect
   assert(!SCIPisFeasIntegral(scip, constantSum));

   /* create two nodes */
   SCIPdebugMessage("Component bound branching rule: creating 2 nodes\n");
   SCIP_CALL( initNodeBranchdata(scip, &downBranchData, GCG_BRANCH_DOWN, SCIPfloor(scip, constantSum), blocknr) );
   SCIP_CALL( initNodeBranchdata(scip, &upBranchData, GCG_BRANCH_UP, SCIPceil(scip, constantSum), blocknr) );

   return SCIP_OKAY;
}

/** prepares information for using the generic branching scheme */
SCIP_RETCODE GCGbranchCompBndInitbranch(
   SCIP*                 masterscip,              /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,              /**< branching rule */
   SCIP_RESULT*          result                   /**< pointer to store the result of the branching call */
   )
{
   SCIP* origscip;
   SCIP_VAR** branchcands;
   SCIP_VAR** allorigvars;
   SCIP_VAR** mastervars;
   int nmastervars;
   SCIP_CONS* masterbranchcons;
   int nbranchcands;
   GCG_BRANCHDATA* branchdata;
   SCIP_VAR* mastervar;
   SCIP_Real mastervarValue;
   GCG_COMPSEQUENCE* B;
   int Bsize;
   int blocknr;
   int i;
   int j;
   int allnorigvars;

   blocknr = -2;
   B = NULL;

   assert(masterscip != NULL);

   SCIPdebugMessage("get information for component bound branching\n");

   origscip = GCGmasterGetOrigprob(masterscip);

   assert(origscip != NULL);
   SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL, NULL) );

   SCIP_CALL( SCIPgetVarsData(origscip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   /* in case original problem contains continuous variables, there are no branching cands */
   assert(nbranchcands > 0 || SCIPgetNContVars(origscip) > 0);
   mastervar = NULL;

   /* 1. Determine in what block we are branching. We select the first available block,
    *     i.e. the first block that contains a branching candidate, starting from the master block.
    */
   /* in case of continuous origvar look for "fractional" blocks using the representation (currentorigsol) in the original problem */
   if(SCIPgetNContVars(origscip) > 0)
   {
      int norigvars;
      SCIP_VAR** origvars;
      SCIP_VAR* origvar;

      norigvars = SCIPgetNVars(origscip);
      origvars = SCIPgetVars(origscip);

      nbranchcands = SCIPgetNVars(masterscip);
      branchcands = SCIPgetVars(masterscip);

      assert(nbranchcands > 0);

      for( i = 0; i < norigvars; ++i )
      {
         int k;
         SCIP_Bool checked;
         origvar = origvars[i];

         if( SCIPvarGetType(origvar) > SCIP_VARTYPE_INTEGER )
            continue;

         if( SCIPisIntegral(origscip, SCIPgetSolVal(origscip, GCGrelaxGetCurrentOrigSol(origscip), origvar)) )
            continue;

         blocknr = GCGgetBlockRepresentative(origscip, GCGvarGetBlock(origvar));

         SCIPdebugMessage("Variable %s belonging to block %d with representative %d is not integral!\n", SCIPvarGetName(origvar), GCGvarGetBlock(origvar), blocknr);

         if( blocknr == -1 )
         {
            assert(GCGoriginalVarGetNMastervars(origvar) == 1);
            mastervar = GCGoriginalVarGetMastervars(origvar)[0];
            break;
         }

         break;
      }
   } else
   {
      /* loop over all branching candidates */
      for( i = 0; i < nbranchcands; ++i )
      {
         int k;
         mastervar = branchcands[i];
         assert(GCGvarIsMaster(mastervar));

         /* if we have a master variable, we branch on it */
         if( GCGvarGetBlock(mastervar) == -1 )
         {
            assert(!GCGmasterVarIsArtificial(mastervar));
            blocknr = -1;
            break;
         }

         /* else, check if the candidate is in an procing block */
         for( j = 0; j < GCGgetNPricingprobs(origscip); ++j )
         {
            if( GCGisMasterVarInBlock(mastervar, j) )
            {
               blocknr = j;
               break;
            }
         }
      }
   }

   if( blocknr < -1 )
   {
      SCIPdebugMessage("Generic branching rule could not find variables to branch on!\n");
      SCIP_Bool rays;
      SCIP_CALL( GCGpricerExistRays(masterscip, &rays) );
      if( rays )
         SCIPwarningMessage(masterscip, "Generic branching is not compatible with unbounded problems!\n");
      return SCIP_ERROR;
   }

   /* a special case; branch on copy of an origvar directly */
#ifdef SCIP_DISABLED_CODE
   // TODO-TIL: Generic Branching handles this special case, but I have no idea why and for what purpose.
   if( blocknr == -1 )
   {
      assert(!GCGmasterVarIsLinking(mastervar));
      SCIPdebugMessage("branching on master variable\n");
      SCIP_CALL( branchDirectlyOnMastervar(origscip, mastervar, branchrule) );  // ownership of mastervar moved. TODO-TIL: check if mastervar is modified
      return SCIP_OKAY;
   }
#endif

   masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);
   SCIPdebugMessage("branching in block %d \n", blocknr);

   /* 2. Call to seperation algorithm to find a suitable B to branch on in the current block.*/
   // TODO-TIL: Call to unimplemented Seperation algorithm
   assert(Bsize > 0);
   assert(B != NULL);

   /* 3. Create the child nodes. */
   SCIP_CALL( createChildNodesCompBnd(masterscip, branchrule, B, Bsize, blocknr, result) );

   // TODO-TIL: Check if freeing B is necessary
   SCIPfreeBufferArrayNull(origscip, &B);
   SCIPdebugMessage("free F\n");

   return SCIP_OKAY;
}

/*
 * Callback methods of event handler
 *
 * The event handler is necessary to react on variable additions to the master problem.
 * If the new master variable satisfies all the component bound constraints, it must be added to the branching constraint.
 */

/* define not used callback as NULL*/
#define branchFreeCompBnd NULL
#define branchExitCompBnd NULL
#define branchInitsolCompBnd NULL
#define branchExitsolCompBnd NULL

/** adds a variable to a branching constraint */
static
SCIP_RETCODE addVarToMasterbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             mastervar,          /**< the variable to add */
   GCG_BRANCHDATA*       branchdata,         /**< branching data structure where the variable should be added */
   SCIP_Bool*            added               /**< whether the variable was added */
)
{
   int i;

   assert(scip != NULL);
   assert(mastervar != NULL);
   assert(branchdata != NULL);
   assert(added != NULL);

   *added = FALSE;

   // TODO-TIL: What the hell is this doing? I assume block -1 is master? what is block -3?
   if( GCGvarGetBlock(mastervar) == -1 || branchdata->consblocknr == -3 || !GCGisMasterVarInBlock(mastervar, branchdata->consblocknr) )
      return SCIP_OKAY;

   SCIPdebugMessage("consSsize = %d\n", branchdata->Bsize);

   if( isMasterVarInB(mastervar, branchdata->B, branchdata->Bsize) )
   {
      SCIPdebugMessage("mastervar is added\n");
      SCIP_CALL( SCIPaddCoefLinear(scip, branchdata->mastercons, mastervar, 1.0) );
      *added = TRUE;
   }

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolCompBndbranchvaradd)
{  /*lint --e{715}*/
    assert(scip != NULL);
    assert(eventhdlr != NULL);
    assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

    /* notify SCIP that your event handler wants to react on the event type */
    SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, NULL) );

    return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolCompBndbranchvaradd)
{  /*lint --e{715}*/
    assert(scip != NULL);
    assert(eventhdlr != NULL);
    assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

    /* notify SCIP that your event handler wants to drop the event type */
    SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, -1) );

    return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecCompBndbranchvaradd)
{  /*lint --e{715}*/
   SCIP* origscip;
   SCIP_CONS* masterbranchcons;
   SCIP_CONS* parentcons;
   SCIP_VAR* mastervar;
   SCIP_VAR** allorigvars;
   SCIP_VAR** mastervars;
   GCG_BRANCHDATA* branchdata;
   int allnorigvars;
   int nmastervars;

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(scip != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

   mastervar = SCIPeventGetVar(event);
   if( !GCGvarIsMaster(mastervar) )
      return SCIP_OKAY;

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   masterbranchcons = GCGconsMasterbranchGetActiveCons(scip);
   assert(masterbranchcons != NULL);

   /* if branch rule is not component bound, abort */
   if( !GCGisBranchruleCompBnd(GCGconsMasterbranchGetBranchrule(masterbranchcons)) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetVarsData(origscip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   parentcons = masterbranchcons;
   branchdata = GCGconsMasterbranchGetBranchdata(parentcons);


   if( GCGvarIsMaster(mastervar) && GCGconsMasterbranchGetBranchrule(parentcons) != NULL )
   {
      SCIP_Bool added = FALSE;
      SCIPdebugMessage("Mastervar <%s>\n", SCIPvarGetName(mastervar));
      while( parentcons != NULL && branchdata != NULL
            && branchdata->B != NULL && branchdata->Bsize > 0 )
      {
         if( GCGconsMasterbranchGetBranchrule(parentcons) == NULL || strcmp(SCIPbranchruleGetName(GCGconsMasterbranchGetBranchrule(parentcons)), "generic") != 0 )
            break;

         assert(branchdata != NULL);


         if( (branchdata->consblocknr != GCGvarGetBlock(mastervar) && GCGvarGetBlock(mastervar) != -1 )
               || (GCGvarGetBlock(mastervar) == -1 && !GCGmasterVarIsLinking(mastervar)) )
         {
            parentcons = GCGconsMasterbranchGetParentcons(parentcons);

            if( parentcons != NULL )
               branchdata = GCGconsMasterbranchGetBranchdata(parentcons);

            continue;
         }

         SCIP_CALL( addVarToMasterbranch(scip, mastervar, branchdata, &added) );

         parentcons = GCGconsMasterbranchGetParentcons(parentcons);
         branchdata = GCGconsMasterbranchGetBranchdata(parentcons);
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeCompBnd NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitCompBnd NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolCompBnd NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolCompBnd NULL


/** branching execution method for fractional LP solutions */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHEXECLP(branchExeclpCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExeclpCompBnd NULL
#endif


static
SCIP_DECL_BRANCHEXECEXT(branchExecextCompBnd)
{  /*lint --e{715}*/
   SCIPdebugMessage("Execext method of generic branching\n");

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/* from branch_master */
static
SCIP_RETCODE GCGincludeMasterCopyPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPincludeNodeselBfs(scip) );
   SCIP_CALL( SCIPincludeNodeselDfs(scip) );
   SCIP_CALL( SCIPincludeNodeselEstimate(scip) );
   SCIP_CALL( SCIPincludeNodeselHybridestim(scip) );
   SCIP_CALL( SCIPincludeNodeselRestartdfs(scip) );
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
/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyCompBnd)
{
   assert(branchrule != NULL);
   assert(scip != NULL);
   SCIP_CALL( GCGincludeMasterCopyPlugins(scip) );
   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsCompBnd)
{  /*lint --e{715}*/
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of component bound branching\n");

   if( SCIPisStopped(scip) )
   {
      SCIPwarningMessage(scip, "No branching could be created, solving process cannot be restarted...\n" );

      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   else
   {
      SCIPerrorMessage("This method is not implemented, aborting since we cannot recover!\n");
      SCIPdialogMessage(scip, NULL, "Due to numerical issues, the problem could not be solved.\n");
      SCIPdialogMessage(scip, NULL, "You can try to disable discretization and aggregation and resolve the problem.\n");

      *result = SCIP_DIDNOTRUN;
      return SCIP_ERROR;
   }
}

/*
 * GCG specific branching rule callbacks
 */

/** activation method for branchrule, called when a node in the master problem is activated,
 *  should perform changes to the current node's problem due to the branchdata
 */
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterCompBnd)
{
   SCIP_VAR** mastervars;
   int nmastervars;
   int i;
   int nvarsadded;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   SCIPdebugMessage("branchActiveMasterCompBnd: Block %d, Ssize %d\n", branchdata->consblocknr, branchdata->Bsize);

   if( branchdata->nvars >= SCIPgetNVars(scip) )
      return SCIP_OKAY;

   nmastervars = SCIPgetNVars(scip);
   mastervars = SCIPgetVars(scip);
   nvarsadded = 0;

   for( i = branchdata->nvars; i < nmastervars; ++i )
   {
      SCIP_Bool added = FALSE;
      assert(mastervars[i] != NULL);
      assert(GCGvarIsMaster(mastervars[i]));

      SCIP_CALL( addVarToMasterbranch(scip, mastervars[i], branchdata, &added) );
      if( added )
         ++nvarsadded;

   }
   SCIPdebugMessage("%d/%d vars added with contant=%g\n", nvarsadded, nmastervars-branchdata->nvars, branchdata->constant);
   branchdata->nvars = nmastervars;

   return SCIP_OKAY;
}


/** deactivation method for branchrule, called when a node in the master problem is deactivated,
 *  should undo changes to the current node's problem due to the branchdata
 */
static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterCompBnd)
{
   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   SCIPdebugMessage("branchDeactiveMasterCompBnd: Block %d, Ssize %d\n", branchdata->consblocknr, branchdata->Bsize);

   /* set number of variables since last call */
   branchdata->nvars = SCIPgetNVars(scip);
   return SCIP_OKAY;
}

/** propagation method for branchrule, called when a node in the master problem is propagated,
 *  should perform propagation at the current node due to the branchdata
 */
static
GCG_DECL_BRANCHPROPMASTER(branchPropMasterCompBnd)
{
   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);
   assert(branchdata->B != NULL);

   *result = SCIP_DIDNOTFIND;
   return SCIP_OKAY;
}

/** method for branchrule, called when the master LP is solved at one node,
 *  can store pseudocosts for the branching decisions
 */
#define branchMasterSolvedCompBnd NULL

/** frees branching data of an origbranch constraint (called when the origbranch constraint is deleted) */
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteCompBnd)
{
   assert(scip != NULL);
   assert(branchdata != NULL);

   if( *branchdata == NULL )
   {
      SCIPdebugMessage("branchDataDeleteCompBnd: cannot delete empty branchdata\n");

      return SCIP_OKAY;
   }

   if( (*branchdata)->mastercons != NULL )
   {
      SCIPdebugMessage("branchDataDeleteCompBnd: child blocknr %d, %s\n", (*branchdata)->consblocknr,
         SCIPconsGetName((*branchdata)->mastercons) );
   }
   else
   {
      SCIPdebugMessage("branchDataDeleteCompBnd: child blocknr %d, empty mastercons\n", (*branchdata)->consblocknr);
   }

   /* release constraint that enforces the branching decision */
   if( (*branchdata)->mastercons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(GCGgetMasterprob(scip), &(*branchdata)->mastercons) );
      (*branchdata)->mastercons = NULL;
   }

   if( (*branchdata)->B != NULL && (*branchdata)->Bsize > 0 )
   {
      // TODO-TIL: I just changed the parameter to (*branchdata)->Bsize from (*branchdata)->maxconsSsize, without knowing what I am doing
      SCIPfreeBlockMemoryArrayNull(scip, &((*branchdata)->B), (*branchdata)->Bsize);
      (*branchdata)->B = NULL;
      (*branchdata)->Bsize = 0;
   }

   SCIPfreeBlockMemoryNull(scip, branchdata);
   *branchdata = NULL;

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitCompBnd)
{
   SCIP* origscip;

   origscip = GCGmasterGetOrigprob(scip);
   assert(branchrule != NULL);
   assert(origscip != NULL);

   SCIPdebugMessage("Init method of Vanderbecks generic branching\n");

   SCIP_CALL( GCGrelaxIncludeBranchrule(origscip, branchrule, branchActiveMasterCompBnd,
         branchDeactiveMasterCompBnd, branchPropMasterCompBnd, branchMasterSolvedCompBnd, branchDataDeleteCompBnd) );

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the compbnd branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleCompBnd(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create compbnd branching rule data */
   branchruledata = NULL;

   SCIPdebugMessage("Include method of component bound branching\n");

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
         BRANCHRULE_MAXBOUNDDIST,
         branchCopyCompBnd, branchFreeCompBnd, branchInitCompBnd, branchExitCompBnd, branchInitsolCompBnd, branchExitsolCompBnd,
         branchExeclpCompBnd, branchExecextCompBnd, branchExecpsCompBnd,
         branchruledata) );

   /* include event handler for adding generated mastervars to the branching constraints */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, eventInitsolCompBndbranchvaradd, eventExitsolCompBndbranchvaradd,
         NULL, eventExecCompBndbranchvaradd,
         NULL) );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   SCIP_CALL( GCGconsIntegralorigAddBranchrule(scip, branchrule) );

   return SCIP_OKAY;
}

/** returns true when the branch rule is the generic branchrule */
SCIP_Bool GCGisBranchruleCompBnd(
   SCIP_BRANCHRULE*      branchrule          /**< branchrule to check */
)
{
   return (branchrule != NULL) && (strcmp(BRANCHRULE_NAME, SCIPbranchruleGetName(branchrule)) == 0);
}
