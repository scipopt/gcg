/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"
//#define SCIP_DEBUG
/**@file   branch_gcg.c
 * @ingroup BRANCHINGRULES
 * @brief  most infeasible LP branching rule
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_gcg.h"
#include "relax_gcg.h"
#include "cons_origbranch.h"


#define BRANCHRULE_NAME          "gcg"
#define BRANCHRULE_DESC          "branching for generic column generation, doing most infeasible branching"
#define BRANCHRULE_PRIORITY      100
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0




/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeGcg NULL


/** initialization method of branching rule (called after problem was transformed) */
#define branchInitGcg NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitGcg NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolGcg NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolGcg NULL


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpGcg)
{  
   SCIPdebugMessage("Execlp method of gcg branching\n");

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsGcg)
{
   SCIP_SOL* currentsol;
   SCIP_VAR** vars;
   int nvars;
   int nbinvars;
   int nintvars;

   SCIP_Real infeasibility;
   SCIP_Real score;
   SCIP_Real obj;
   SCIP_Real bestscore;
   SCIP_Real bestobj;
   int bestcand;
   int i;

   SCIP_NODE* childup;
   SCIP_NODE* childdown;
   SCIP_CONS* consup;
   SCIP_CONS* consdown;

   SCIP_CONS* origbranchup;
   SCIP_CONS* origbranchdown;


   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of gcg branching\n");

   /* get current sol */
   currentsol = GCGrelaxGetCurrentOrigSol(scip);

   /* get the variables of the original problem and the numbers of variable types */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* search for an integer variable with fractional value */
   for ( i = 0; i < nbinvars + nintvars; i++ )
   {
      assert(SCIPvarGetType(vars[i]) == (i < nbinvars ? SCIP_VARTYPE_BINARY : SCIP_VARTYPE_INTEGER));
      if ( !SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, currentsol, vars[i])))
      {
         SCIPdebugMessage("Var %s has fractional value in current solution: %f\n", SCIPvarGetName(vars[i]), SCIPgetSolVal(scip, currentsol, vars[i]));
         break;
      }
   }

   /* create the b&b-tree child-nodes of the current node */
   SCIP_CALL( SCIPcreateChild(scip, &childup, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &childdown, 0.0, SCIPgetLocalTransEstimate(scip)) );

   /* create corresponding constraints */
   SCIP_CALL( SCIPcreateConsLinear(scip, &consup, "branch_up", 0, NULL, NULL, 
         SCIPceil(scip, SCIPgetSolVal(scip, currentsol, vars[i])), SCIPinfinity(scip), 
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcreateConsLinear(scip, &consdown, "branch_down", 0, NULL, NULL, 
         -1.0 * SCIPinfinity(scip), SCIPfloor(scip, SCIPgetSolVal(scip, currentsol, vars[i])),  
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddCoefLinear(scip, consup, vars[i], 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, consdown, vars[i], 1.0) );

   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchup, "branchup", consup, vars[i], GCG_CONSSENSE_GE, 
         SCIPceil(scip, SCIPgetSolVal(scip, currentsol, vars[i]))) );
   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchdown, "branchdown", consdown, vars[i], GCG_CONSSENSE_LE, 
         SCIPfloor(scip, SCIPgetSolVal(scip, currentsol, vars[i]))) );

   /* add constraints to nodes */
   SCIP_CALL( SCIPaddConsNode(scip, childup, consup, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdown, consdown, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childup, origbranchup, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdown, origbranchdown, NULL) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &consup) );
   SCIP_CALL( SCIPreleaseCons(scip, &consdown) );
   SCIP_CALL( SCIPreleaseCons(scip, &origbranchup) );
   SCIP_CALL( SCIPreleaseCons(scip, &origbranchdown) );

      
   *result = SCIP_BRANCHED;


   return SCIP_OKAY;
}




/*
 * branching specific interface methods
 */

/** creates the most infeasible LP braching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleGcg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create inference branching rule data */
   branchruledata = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchFreeGcg, branchInitGcg, branchExitGcg, branchInitsolGcg, branchExitsolGcg, 
         branchExeclpGcg, branchExecpsGcg,
         branchruledata) );

   return SCIP_OKAY;
}
