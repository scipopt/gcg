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
/**@file   cons_integralOrig.c
 * @brief  constraint handler for storing the graph at each node of the tree
 * @author Gerald Gamrath
 *
 */

#include <assert.h>
#include <string.h>

#include "scip/type_cons.h"
#include "cons_integralOrig.h"
#include "probdata_gcg.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "integralOrig"
#define CONSHDLR_DESC          "storing graph at nodes of the tree constraint handler"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              * propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */


/** constraint data for storing graph constraints */
struct SCIP_ConsData
{
   
};



/** constraint handler data */
struct SCIP_ConshdlrData
{
   int nblocks;
   
};


/*
 * Local methods
 */



/*
 * Callback methods
 */


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeIntegralOrig)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   SCIPdebugMessage("freeing store graph constraint handler\n");

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolIntegralOrig)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONS*         cons;
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolIntegralOrig)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteIntegralOrig)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(*consdata != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   
   SCIPdebugMessage("Deleting integral orig constraint: <%s>.\n", SCIPconsGetName(cons));

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpIntegralOrig)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsIntegralOrig)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckIntegralOrig)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockIntegralOrig)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   SCIPdebugMessage("Locking method for integral orig constraint: <%s>.\n", SCIPconsGetName(cons));

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveIntegralOrig)
{
   return SCIP_OKAY;
}



/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveIntegralOrig)
{
   return SCIP_OKAY;
}



/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropIntegralOrig)
{
   return SCIP_OKAY;
}

/* define not used callbacks as NULL */
#define consPresolIntegralOrig NULL
#define consRespropIntegralOrig NULL
#define consInitIntegralOrig NULL
#define consExitIntegralOrig NULL
#define consInitpreIntegralOrig NULL
#define consExitpreIntegralOrig NULL
#define consTransIntegralOrig NULL
#define consInitlpIntegralOrig NULL
#define consSepalpIntegralOrig NULL
#define consSepasolIntegralOrig NULL
#define consEnableIntegralOrig NULL
#define consDisableIntegralOrig NULL
#define consPrintIntegralOrig NULL


/*
 * interface methods
 */


/** creates the handler for integralOrig constraints and includes it in SCIP */
SCIP_RETCODE COLORincludeConshdlrIntegralOrig(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrData;

   SCIPdebugMessage("Including graph storage constraint handler.\n");

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrData) );

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
         consPrintIntegralOrig,
         conshdlrData) );

   return SCIP_OKAY;
}


/** creates and captures a integralOrig constraint, uses knowledge of the B&B-father*/
SCIP_RETCODE COLORcreateConsIntegralOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name                /**< name of constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);

   /* find the integralOrig constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("integralOrig constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   SCIPdebugMessage("Creating store graph constraint: <%s(%d,%d)>. \n", name, (node1+1), (node2+1));

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
         TRUE, FALSE, TRUE, FALSE, TRUE) );
   
   return SCIP_OKAY;
}




/* ----------------------------------- external methods -------------------------- */

