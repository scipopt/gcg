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

/**@file   probdata_gcg.c
 * @brief  problem data for generic column generation
 * @author Gerald Gamrath
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "probdata_gcg.h"
//#include "scip/cons_linear.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <assert.h>


#define MAX_LINELEN       65536


struct SCIP_ProbData
{
   SCIP*            origprob;           /* the original problem */
   SCIP**           pricingprobs;       /* the array of pricing problems */
   int              npricingprobs;      /* the number of pricing problems */
   SCIP_CONS**      masterconss;        /* array of cons in the original problem that belong to the master problem */
   int              maxmasterconss;     /* length of the array mastercons */
   int              nmasterconss;       /* number of constraints saved in mastercons */
   SCIP_CONS***     pricingconss;       /* array of arrays of cons in the original problem belonging to pricing problem nr. i*/
   int*             maxpricingconss;    /* length of pricingcons at position i*/
   int*             npricingconss;      /* number of constraints belonging to pricing problem i */
};


/*
 * Local methods
 */

static
SCIP_RETCODE ensureSizeMasterConss(
   SCIP*                 scip,
   SCIP_PROBDATA*        probdata,
   int                   size
   )
{
   assert(scip != NULL);
   assert(probdata != NULL);
   assert(probdata->masterconss != NULL);
   
   if ( probdata->maxmasterconss < size )
   {
      probdata->maxmasterconss = MAX(probdata->maxmasterconss + 5, size);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(probdata->masterconss), probdata->maxmasterconss) );
   }
   assert(probdata->maxmasterconss >= size);

   return SCIP_OKAY;
}

static
SCIP_RETCODE ensureSizePricingConss(
   SCIP*                 scip,
   SCIP_PROBDATA*        probdata,
   int                   index,
   int                   size
   )
{
   assert(scip != NULL);
   assert(probdata != NULL);
   assert(0 <= index && index < probdata->npricingprobs);
   assert(probdata->pricingconss != NULL);
   assert(probdata->pricingconss[index] != NULL);
   
   if ( probdata->maxpricingconss[index] < size )
   {
      probdata->maxpricingconss[index] = MAX(probdata->maxpricingconss[index]+5, size);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(probdata->pricingconss[index]), probdata->maxpricingconss[index]) );
   }
   assert(probdata->maxpricingconss[index] >= size);

   return SCIP_OKAY;
}

/*
 * Callback methods of probdata
 */

/** transforms the problem */
static
SCIP_DECL_PROBTRANS(probtransGcg)
{
   assert(scip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   (*targetdata)->origprob = sourcedata->origprob;        

   return SCIP_OKAY;
}


/** deletes the transformed problem */
static
SCIP_DECL_PROBDELTRANS(probdeltransGcg)
{
   return SCIP_OKAY;
}


#define probinitsolGcg NULL
#define probexitsolGcg NULL


static
SCIP_DECL_PROBDELORIG(probdelorigGcg)
{
   assert(probdata != NULL);
   assert(*probdata != NULL);

   /* free original problem */
   SCIP_CALL( SCIPfreeTransform((*probdata)->origprob) );
   SCIP_CALL( SCIPfree(&((*probdata)->origprob)) );

   return SCIP_OKAY;
}




/*
 * probdata specific interface methods
 */

/** sets up the problem data */
SCIP_RETCODE SCIPcreateProbGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */
   int                   npricingprobs       /**< number of pricing problems */
   )
{
   char* probname;
   SCIP_PROBDATA* probdata;
   int i;
   
   printf("Creating problem: %s\n", name);
   
   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, &probdata) );

   /* initialize data */
   probdata->npricingprobs = npricingprobs;
   probdata->maxmasterconss = 5;


   /* array of constraints belonging to the master problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->masterconss), probdata->maxmasterconss) );
   
   /* array for saving of constraints belonging to one of the pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->pricingconss), npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->maxpricingconss), npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->npricingconss), npricingprobs) );
   for ( i = 0; i < npricingprobs; i++ )
   {
      probdata->maxpricingconss[i] = 5;
      probdata->npricingconss[i] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->pricingconss[i]), probdata->maxpricingconss[i]) );
   }
   
   /* initializing the scip data structure for the original problem */  
   SCIP_CALL( SCIPcreate(&(probdata->origprob)) );
   SCIP_CALL( SCIPincludeDefaultPlugins(probdata->origprob) );
 
   /* disable conflict analysis */
   SCIP_CALL( SCIPsetBoolParam(probdata->origprob, "conflict/useprop", FALSE) ); 
   SCIP_CALL( SCIPsetBoolParam(probdata->origprob, "conflict/useinflp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(probdata->origprob, "conflict/useboundlp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(probdata->origprob, "conflict/usesb", FALSE) ); 
   SCIP_CALL( SCIPsetBoolParam(probdata->origprob, "conflict/usepseudo", FALSE) );

   SCIP_CALL( SCIPcreateProb(probdata->origprob, "origprob", NULL, NULL, NULL, NULL, NULL, NULL) );
   
   /* initialize the pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->pricingprobs), npricingprobs) );
   for ( i = 0; i < npricingprobs; i++ )
   {
      //SCIP_CALL( SCIPallocBufferArray(scip, probname, MAX_LINELEN) );
      //(void) SCIPsnprintf(probname, MAX_LINELEN, "pricing_block_%n", i);

      /* initializing the scip data structure for the original problem */  
      SCIP_CALL( SCIPcreate(&(probdata->pricingprobs[i])) );
      SCIP_CALL( SCIPincludeDefaultPlugins(probdata->pricingprobs[i]) );
 
      /* disable conflict analysis */
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/useprop", FALSE) ); 
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/useinflp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/useboundlp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/usesb", FALSE) ); 
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/usepseudo", FALSE) );

      SCIP_CALL( SCIPcreateProb(probdata->pricingprobs[i], "origprob", NULL, NULL, NULL, NULL, NULL, NULL) );
   }

   /* create problem in SCIP */
   SCIP_CALL( SCIPcreateProb(scip, name, probdelorigGcg, probtransGcg, probdeltransGcg, 
         probinitsolGcg, probexitsolGcg, probdata) );


   return SCIP_OKAY;
}


/* ----------------------------------- external methods -------------------------- */


/** creates and captures a linear constraint */
SCIP_RETCODE GCGcreateConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< Is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   int                   pricingprobnr       /**< number of the pricing problem this constraint belongs to,
                                              *   -1 if the constraint should appear in the master problem? */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS* cons;

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);

   /* add constraint to the original program */
   SCIP_CALL( SCIPcreateConsLinear(probdata->origprob, &cons, name, nvars, vars, vals, lhs, rhs, initial, 
         separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   SCIP_CALL( SCIPaddCons(probdata->origprob, cons) );
   SCIPdebugMessage("added constraint %s to the original problem\n", name);
   //SCIP_CALL( SCIPreleaseCons(probdata->origprob, &cons) );

   /* if constraint belongs to the master problem, create it for the scip instance */
   if ( pricingprobnr == -1 )
   {
      /* add constraint to array masterconss */
      SCIP_CALL( ensureSizeMasterConss(scip, probdata, probdata->nmasterconss+1) );
      probdata->masterconss[probdata->nmasterconss] = cons;
      probdata->nmasterconss++;
      
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nvars, vars, vals, lhs, rhs, initial, 
            separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMessage("added constraint %s to the master problem\n", name);
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* constraint belongs to one of the pricing problems */
   else
   {
      /* add constraint to array pricingconss */
      SCIP_CALL( ensureSizePricingConss(scip, probdata, pricingprobnr, probdata->npricingconss[pricingprobnr]+1) );
      probdata->pricingconss[pricingprobnr][probdata->npricingconss[pricingprobnr]] = cons;
      probdata->npricingconss[pricingprobnr]++;
      
   }

   return SCIP_OKAY;

}


/** creates a variable of the original program */
SCIP_RETCODE GCGcreateVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable           /**< is var's column removable from the LP (due to aging or cleanup)? */
   )
{
   SCIP_PROBDATA* probdata;

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);
   assert(var != NULL);
   assert(lb <= ub);

   SCIP_CALL( SCIPcreateVar(probdata->origprob, var, name, lb, ub, obj, vartype, initial, removable, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}


/** adds variable to the original problem */
SCIP_RETCODE GCGaddOriginalVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to add */
   )
{
   SCIP_PROBDATA* probdata;

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPaddVar(probdata->origprob, var) );

   /* because the variable was added to the problem, it is captured by SCIP and we can safely release it right now
    * without making the returned *var invalid
    */
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   return SCIP_OKAY;
}


SCIP* GCGprobGetOrigprob(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   
   probdata = SCIPgetProbData(scip);
   
   assert(probdata != NULL);
   
   return probdata->origprob;
}
