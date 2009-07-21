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
//#define SCIP_DEBUG
/**@file   cons_infeas.c
 * @brief  constraint handler for storing the branching decisions at each node of the tree
 * @author Gerald Gamrath
 *
 */

#include <assert.h>
#include <string.h>

#include "scip/type_cons.h"
#include "cons_infeas.h"
#include "cons_masterbranch.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "infeas"
#define CONSHDLR_DESC          "store branching decision at nodes of the tree constraint handler"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   1000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              * propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/*
 * Callback methods
 */


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpInfeas)
{
   SCIP_NODE* child1;
   SCIP_NODE* child2;
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Enfolp method of infeas constraint handler.\n");

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


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsInfeas)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   SCIPerrorMessage("consEnfopsInfeas() called - this should not happen!\n");

   /* do nothing */
   *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckInfeas)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   SCIPdebugMessage("Check method of infeas constraint handler.\n");

   /* do nothing */
   *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockInfeas)
{
   return SCIP_OKAY;
}


/* define not used callbacks as NULL */
#define consFreeInfeas NULL
#define consInitsolInfeas NULL
#define consExitsolInfeas NULL
#define consDeleteInfeas NULL
#define consActiveInfeas NULL
#define consDeactiveInfeas NULL
#define consPropInfeas NULL
#define consPresolInfeas NULL
#define consRespropInfeas NULL
#define consInitInfeas NULL
#define consExitInfeas NULL
#define consInitpreInfeas NULL
#define consExitpreInfeas NULL
#define consTransInfeas NULL
#define consInitlpInfeas NULL
#define consSepalpInfeas NULL
#define consSepasolInfeas NULL
#define consEnableInfeas NULL
#define consDisableInfeas NULL
#define consPrintInfeas NULL
#define consCopyInfeas NULL
#define consParseInfeas NULL


/*
 * interface methods
 */


/** creates the handler for infeas constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrInfeas(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrData;

   SCIPdebugMessage("Including branch orig constraint handler.\n");

   conshdlrData = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeInfeas, consInitInfeas, consExitInfeas,
         consInitpreInfeas, consExitpreInfeas, consInitsolInfeas, consExitsolInfeas,
         consDeleteInfeas, consTransInfeas, consInitlpInfeas,
         consSepalpInfeas, consSepasolInfeas, consEnfolpInfeas, consEnfopsInfeas, consCheckInfeas,
         consPropInfeas, consPresolInfeas, consRespropInfeas, consLockInfeas,
         consActiveInfeas, consDeactiveInfeas,
         consEnableInfeas, consDisableInfeas,
         consPrintInfeas, consCopyInfeas, consParseInfeas, 
         conshdlrData) );

   return SCIP_OKAY;
}

/* ----------------------------------- external methods -------------------------- */



