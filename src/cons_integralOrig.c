/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   cons_integralOrig.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for the integrality constraint in the original problem
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_integralOrig.h"
#include "pricer_gcg.h"

#include "struct_vardata.h"
#include "cons_masterbranch.h"


#define CONSHDLR_NAME          "integralOrig"
#define CONSHDLR_DESC          "integrality constraint"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY      1000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY     1000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */




/*
 * Callback methods
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#define consFreeIntegralOrig NULL


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitIntegralOrig NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitIntegralOrig NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreIntegralOrig NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreIntegralOrig NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolIntegralOrig NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#define consExitsolIntegralOrig NULL


/** frees specific constraint data */
#define consDeleteIntegralOrig NULL


/** transforms constraint data into data belonging to the transformed problem */ 
#define consTransIntegralOrig NULL


/** LP initialization method of constraint handler */
#define consInitlpIntegralOrig NULL


/** separation method of constraint handler for LP solutions */
#define consSepalpIntegralOrig NULL


/** separation method of constraint handler for arbitrary primal solutions */
#define consSepasolIntegralOrig NULL


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpIntegralOrig)
{  
   SCIP* origprob;
   SCIP_VAR** origvars;
   SCIP_VARDATA* vardata;
   int norigvars;
   SCIP_Real solval;
   SCIP_Bool discretization;
   int v;
   int i;

   SCIP_NODE* child1;
   SCIP_NODE* child2;
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(conss == NULL);
   assert(nconss == 0);
   assert(result != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   SCIPdebugMessage("LP solution enforcing method of integralOrig constraint\n");

   *result = SCIP_FEASIBLE;

   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( discretization )
   {
      return SCIP_OKAY;
   }

   origvars = SCIPgetOrigVars(origprob);
   norigvars = SCIPgetNOrigVars(origprob);

   for( v = 0; v < norigvars; v++ )
   {
      if( SCIPvarGetType(origvars[v]) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      solval = 0;
      vardata = SCIPvarGetData(origvars[v]);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->data.origvardata.mastervars != NULL);
      assert(vardata->data.origvardata.mastervals != NULL);
      assert(vardata->data.origvardata.nmastervars >= 0);

      for( i = 0; i < vardata->data.origvardata.nmastervars; i++ )
      {
         solval += vardata->data.origvardata.mastervals[i] 
            * SCIPgetSolVal(scip, NULL, vardata->data.origvardata.mastervars[i]);
      }
      if( !SCIPisFeasIntegral(scip, solval) )
      {

         /* create the b&b-tree child-nodes of the current node */
         SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
         SCIP_CALL( SCIPcreateChild(scip, &child2, 0.0, SCIPgetLocalTransEstimate(scip)) );
         
         SCIP_CALL( GCGcreateConsMasterbranch(scip, &cons1, child1, GCGconsMasterbranchGetActiveCons(scip)) );
         SCIP_CALL( GCGcreateConsMasterbranch(scip, &cons2, child2, GCGconsMasterbranchGetActiveCons(scip)) );
         
         SCIP_CALL( SCIPaddConsNode(scip, child1, cons1, NULL) );
         SCIP_CALL( SCIPaddConsNode(scip, child2, cons2, NULL) );
         
         /* release constraints */
         SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons2) );

         *result = SCIP_BRANCHED;

         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;

}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsIntegralOrig)
{  
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(conss == NULL);
   assert(nconss == 0);
   assert(result != NULL);

   SCIPdebugMessage("Enfops method of integralOrig constraint: %d fractional variables\n", SCIPgetNLPBranchCands(scip));
   printf("Enfops method of integralOrig constraint!\n");

   /* call branching methods */
   SCIP_CALL( SCIPbranchLP(scip, result) );

   assert(*result == SCIP_BRANCHED);

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckIntegralOrig)
{  
   SCIP* origprob;
   SCIP_VAR** origvars;
   int norigvars;
   SCIP_VARDATA* vardata;
   SCIP_Real solval;
   SCIP_Bool discretization;
   int v;
   int i;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   SCIPdebugMessage("Check method of integralOrig constraint\n");

   *result = SCIP_FEASIBLE;

   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( discretization )
   {
      return SCIP_OKAY;
   }

   origvars = SCIPgetOrigVars(origprob);
   norigvars = SCIPgetNOrigVars(origprob);

   for( v = 0; v < norigvars && *result == SCIP_FEASIBLE; v++ )
   {
      if( SCIPvarGetType(origvars[v]) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      solval = 0;
      vardata = SCIPvarGetData(origvars[v]);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->data.origvardata.mastervars != NULL);
      assert(vardata->data.origvardata.mastervals != NULL);
      assert(vardata->data.origvardata.nmastervars >= 0);

      for ( i = 0; i < vardata->data.origvardata.nmastervars; i++ )
      {
         solval += vardata->data.origvardata.mastervals[i] 
            * SCIPgetSolVal(scip, sol, vardata->data.origvardata.mastervars[i]);
      }
      if( !SCIPisFeasIntegral(scip, solval) )
      {
         *result = SCIP_INFEASIBLE;
         
         if( printreason )
         {
            SCIPinfoMessage(scip, NULL, "violation: integrality condition of variable <%s> = %.15g\n", 
               SCIPvarGetName(origvars[v]), solval);
         }
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#define consPropIntegralOrig NULL

/** presolving method of constraint handler */
#define consPresolIntegralOrig NULL


/** propagation conflict resolving method of constraint handler */
#define consRespropIntegralOrig NULL


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockIntegralOrig)
{  
   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveIntegralOrig NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveIntegralOrig NULL


/** constraint enabling notification method of constraint handler */
#define consEnableIntegralOrig NULL


/** constraint disabling notification method of constraint handler */
#define consDisableIntegralOrig NULL

/** constraint display method of constraint handler */
#define consPrintIntegralOrig NULL

/** constraint copying method of constraint handler */
#define consCopyIntegralOrig NULL

/** constraint parsing method of constraint handler */
#define consParseIntegralOrig NULL


/*
 * constraint specific interface methods
 */

/** creates the handler for integrality constraint and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrIntegralOrig(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create integral constraint handler data */
   conshdlrdata = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeIntegralOrig, consInitIntegralOrig, consExitIntegralOrig, 
         consInitpreIntegralOrig, consExitpreIntegralOrig, consInitsolIntegralOrig, consExitsolIntegralOrig,
         consDeleteIntegralOrig, consTransIntegralOrig, consInitlpIntegralOrig,
         consSepalpIntegralOrig, consSepasolIntegralOrig, consEnfolpIntegralOrig, consEnfopsIntegralOrig, consCheckIntegralOrig, 
         consPropIntegralOrig, consPresolIntegralOrig, consRespropIntegralOrig, consLockIntegralOrig,
         consActiveIntegralOrig, consDeactiveIntegralOrig, 
         consEnableIntegralOrig, consDisableIntegralOrig,
         consPrintIntegralOrig, consCopyIntegralOrig, consParseIntegralOrig,
         conshdlrdata) );

   return SCIP_OKAY;
}
