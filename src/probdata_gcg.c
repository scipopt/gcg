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
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <assert.h>


#define MAX_LINELEN       65536
#define STARTMAXMASTERVARS 10


struct SCIP_ProbData
{
   SCIP*            origprob;           /* the original problem */
   SCIP**           pricingprobs;       /* the array of pricing problems */
   int              npricingprobs;      /* the number of pricing problems */
   SCIP_CONS**      masterconss;        /* array of cons in the master problem */
   int              maxmasterconss;     /* length of the array mastercons */
   int              nmasterconss;       /* number of constraints saved in mastercons */
   SCIP_CONS**      origmasterconss;    /* array of cons in the original problem that belong to the master problem */
   int              maxorigmasterconss; /* length of the array origmastercons */
   int              norigmasterconss;   /* number of constraints saved in origmastercons */
   SCIP_CONS***     pricingconss;       /* array of arrays of cons in the original problem belonging to pricing problem nr. i*/
   int*             maxpricingconss;    /* length of pricingcons at position i*/
   int*             npricingconss;      /* number of constraints belonging to pricing problem i */
   SCIP_CONS**      convconss;          /* array of convexity constraints, one for each block */
};



/*
 * Vardata methods
 */

static
SCIP_DECL_VARDELORIG(gcgvardelorig)
{
   if ( (*vardata)->vartype == GCG_VARTYPE_ORIGINAL )
   {
      SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.mastervars));
      SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.mastervals));
   }
   SCIPfreeBlockMemory(scip, vardata);


   return SCIP_OKAY;
}  


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
   assert(probdata->maxmasterconss == probdata->maxorigmasterconss);
   assert(probdata->nmasterconss == probdata->norigmasterconss);

   if ( probdata->maxmasterconss < size )
   {
      probdata->maxmasterconss = MAX(probdata->maxmasterconss + 5, size);
      probdata->maxorigmasterconss = probdata->maxmasterconss;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(probdata->masterconss), probdata->maxmasterconss) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(probdata->origmasterconss), probdata->maxorigmasterconss) );
   }
   assert(probdata->maxmasterconss >= size);
   assert(probdata->maxorigmasterconss >= size);

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
   int i;
   int j;
   SCIP_VAR** vars;
   int nvars;
   SCIP_Real* vals;
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, targetdata) );

   (*targetdata)->origprob = sourcedata->origprob;
   (*targetdata)->pricingprobs = sourcedata->pricingprobs;
   (*targetdata)->npricingprobs = sourcedata->npricingprobs;
   (*targetdata)->origmasterconss = sourcedata->origmasterconss;
   (*targetdata)->maxorigmasterconss = sourcedata->maxorigmasterconss;
   (*targetdata)->norigmasterconss = sourcedata->norigmasterconss;
   (*targetdata)->pricingconss = sourcedata->pricingconss;
   (*targetdata)->maxpricingconss = sourcedata->maxpricingconss;
   (*targetdata)->npricingconss = sourcedata->npricingconss;

   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->masterconss), sourcedata->maxmasterconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->convconss), sourcedata->npricingprobs) );
   (*targetdata)->maxmasterconss = sourcedata->maxmasterconss;
   (*targetdata)->nmasterconss = sourcedata->nmasterconss;
   
   SCIP_CALL( SCIPtransformConss(scip, sourcedata->nmasterconss, 
         sourcedata->masterconss, (*targetdata)->masterconss) );

   SCIP_CALL( SCIPtransformConss(scip, sourcedata->npricingprobs, 
         sourcedata->convconss, (*targetdata)->convconss) );

   vars = SCIPgetVars((*targetdata)->origprob);
   nvars = SCIPgetNVars((*targetdata)->origprob);
   for ( i = 0; i < nvars; i++ )
   {
      vardata = SCIPvarGetData(vars[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);

      SCIP_CALL( SCIPallocBlockMemoryArray((*targetdata)->origprob, &(vardata->data.origvardata.coefs), (*targetdata)->nmasterconss) );
      vardata->data.origvardata.ncoefs = (*targetdata)->nmasterconss;
      for ( j = 0; j < vardata->data.origvardata.ncoefs; j++ )
      {
         vardata->data.origvardata.coefs[j] = 0;
      }
   }

   /* save coefs in the vardata */
   for ( i = 0; i < (*targetdata)->norigmasterconss; i++ )
   {
      vars = SCIPgetVarsLinear((*targetdata)->origprob, (*targetdata)->origmasterconss[i]);
      nvars = SCIPgetNVarsLinear((*targetdata)->origprob, (*targetdata)->origmasterconss[i]);
      vals = SCIPgetValsLinear((*targetdata)->origprob, (*targetdata)->origmasterconss[i]);
      for ( j = 0; j < nvars; j++ )
      {
         vardata = SCIPvarGetData(vars[j]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         assert(vardata->data.origvardata.coefs != NULL);
         vardata->data.origvardata.coefs[i] = vals[j];
      }

   }

   SCIPpresolve((*targetdata)->origprob);
   //if ( SCIPisObjIntegral((*targetdata)->origprob) )
   //   SCIP_CALL( SCIPsetObjIntegral(scip) );
         
   return SCIP_OKAY;
}


/** deletes the transformed problem */
static
SCIP_DECL_PROBDELTRANS(probdeltransGcg)
{
   int i;
   SCIP_VAR** vars;
   int nvars;
   SCIP_VARDATA* vardata;

   assert(scip !=  NULL);
   assert(probdata != NULL);
   assert(*probdata != NULL);

   for ( i = 0; i < (*probdata)->nmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->masterconss[i]) );
   }
   for ( i = 0; i < (*probdata)->npricingprobs; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->convconss[i]) );
   }


   SCIPfreeTransform((*probdata)->origprob);
   vars = SCIPgetVars((*probdata)->origprob);
   nvars = SCIPgetNVars((*probdata)->origprob);
   for ( i = 0; i < nvars; i++ )
   {
      vardata = SCIPvarGetData(vars[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);

      SCIPfreeBlockMemoryArray((*probdata)->origprob, &(vardata->data.origvardata.coefs), (*probdata)->nmasterconss);
   }
   
   SCIPfreeMemoryArray(scip, &(*probdata)->masterconss);
   SCIPfreeMemoryArray(scip, &(*probdata)->convconss);

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


#define probinitsolGcg NULL
#define probexitsolGcg NULL


static
SCIP_DECL_PROBDELORIG(probdelorigGcg)
{
   int i;
   int j;

   assert(probdata != NULL);
   assert(*probdata != NULL);

   for ( i = 0; i < (*probdata)->npricingprobs; i++ )
   {
      for ( j = 0; j < (*probdata)->npricingconss[i]; j++ )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &((*probdata)->pricingconss[i][j])) );
      }
   }
   for ( i = 0; i < (*probdata)->norigmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->origmasterconss[i]) );
   }
   for ( i = 0; i < (*probdata)->nmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->masterconss[i]) );
   }
   for ( i = 0; i < (*probdata)->npricingprobs; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->convconss[i]) );
   }

   /* free pricing problems */
   for ( i = (*probdata)->npricingprobs - 1; i >= 0 ; i-- )
   {
      SCIP_CALL( SCIPfreeTransform((*probdata)->pricingprobs[i]) );
      SCIP_CALL( SCIPfree(&((*probdata)->pricingprobs[i])) );
   }

   /* free original problem */
   SCIP_CALL( SCIPfreeTransform((*probdata)->origprob) );
   SCIP_CALL( SCIPfree(&((*probdata)->origprob)) );

   /* free array for constraints */
   for ( i = (*probdata)->npricingprobs-1; i >= 0; i-- )
   {
      SCIPfreeMemoryArray(scip, &((*probdata)->pricingconss[i]));
   }
   SCIPfreeMemoryArray(scip, &((*probdata)->pricingconss));
   SCIPfreeMemoryArray(scip, &((*probdata)->maxpricingconss));
   SCIPfreeMemoryArray(scip, &((*probdata)->npricingconss));

   SCIPfreeMemoryArray(scip, &((*probdata)->masterconss));
   SCIPfreeMemoryArray(scip, &((*probdata)->origmasterconss));

   SCIPfreeMemoryArray(scip, &((*probdata)->pricingprobs));

   SCIPfreeMemoryArray(scip, &((*probdata)->convconss));

   SCIPfreeMemory(scip, probdata);

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
   probdata->nmasterconss = 0;
   probdata->maxorigmasterconss = 5;
   probdata->norigmasterconss = 0;


   /* arrays of constraints belonging to the master problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->masterconss), probdata->maxmasterconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->origmasterconss), probdata->maxorigmasterconss) );
   
   /* array for saving constraints belonging to one of the pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->pricingconss), npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->maxpricingconss), npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->npricingconss), npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->convconss), npricingprobs) );
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
   SCIP_CALL( SCIPsetIntParam(probdata->origprob, "presolving/probing/maxrounds", 0) );

   SCIP_CALL( SCIPcreateProb(probdata->origprob, "origprob", NULL, NULL, NULL, NULL, NULL, NULL) );


   SCIP_CALL( SCIPallocBufferArray(scip, &probname, 25) );   
   /* initialize the pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->pricingprobs), npricingprobs) );
   for ( i = 0; i < npricingprobs; i++ )
   {
      /* initializing the scip data structure for the original problem */  
      SCIP_CALL( SCIPcreate(&(probdata->pricingprobs[i])) );
      SCIP_CALL( SCIPincludeDefaultPlugins(probdata->pricingprobs[i]) );
 
      /* disable conflict analysis */
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/useprop", FALSE) ); 
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/useinflp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/useboundlp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/usesb", FALSE) ); 
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "conflict/usepseudo", FALSE) );
      /* disable output to console */
      SCIP_CALL( SCIPsetIntParam(probdata->pricingprobs[i], "display/verblevel", 0) );
      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(probdata->pricingprobs[i], "misc/catchctrlc", FALSE) );

      (void) SCIPsnprintf(probname, 25, "pricing_block_%d", i);
      SCIP_CALL( SCIPcreateProb(probdata->pricingprobs[i], probname, NULL, NULL, NULL, NULL, NULL, NULL) );
   }

   SCIPfreeBufferArray(scip, &probname);

   /* create problem in SCIP */
   SCIP_CALL( SCIPcreateProb(scip, name, probdelorigGcg, probtransGcg, probdeltransGcg, 
         probinitsolGcg, probexitsolGcg, probdata) );

   SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, "gcg")) );

   return SCIP_OKAY;
}


/* ----------------------------------- external methods -------------------------- */

/** create the convexity constraints belonging to the pricing blocks */
SCIP_RETCODE GCGprobCreateConvConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   char* consname;
   int i;

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &consname, 25) );

   /* create convexity constraint for the blocks in the master problem */
   for ( i = 0; i < probdata->npricingprobs; i++ )
   {

      (void) SCIPsnprintf(consname, 25, "conv_block_%d", i);

      SCIP_CALL( SCIPcreateConsLinear(scip, &(probdata->convconss[i]), consname, 0, NULL, NULL, 1, 1, TRUE, 
            TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, probdata->convconss[i]) );
   }

   SCIPfreeBufferArray(scip, &consname);
   
   return SCIP_OKAY;
   
}
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
   SCIP_VARDATA* vardata;
   int i;

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
      /* add constraint to array origmasterconss */
      SCIP_CALL( ensureSizeMasterConss(scip, probdata, probdata->norigmasterconss+1) );
      probdata->origmasterconss[probdata->norigmasterconss] = cons;
      probdata->norigmasterconss++;

      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 0, NULL, NULL, lhs, rhs, initial, 
            separate, enforce, check, propagate, local, TRUE, dynamic, removable, stickingatnode) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMessage("added constraint %s to the master problem\n", name);
      //SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* add constraint to array masterconss */
      //SCIP_CALL( ensureSizeMasterConss(scip, probdata, probdata->nmasterconss+1) );
      probdata->masterconss[probdata->nmasterconss] = cons;
      probdata->nmasterconss++;
   }

   /* constraint belongs to one of the pricing problems */
   else
   {
      SCIP_VAR* var;

      /* add constraint to array pricingconss */
      SCIP_CALL( ensureSizePricingConss(scip, probdata, pricingprobnr, probdata->npricingconss[pricingprobnr]+1) );
      probdata->pricingconss[pricingprobnr][probdata->npricingconss[pricingprobnr]] = cons;
      probdata->npricingconss[pricingprobnr]++;

      /* create (empty) constraint in the pricing problem */
      SCIP_CALL( SCIPcreateConsLinear(probdata->pricingprobs[pricingprobnr], &cons, name, 0, NULL, NULL, lhs, 
            rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, modifiable, FALSE, FALSE, 
            FALSE) );
      SCIP_CALL( SCIPaddCons(probdata->pricingprobs[pricingprobnr], cons) );
      
      /* add variables to the constraint */
      for ( i = 0; i < nvars; i++ )
      {
         assert(!SCIPisFeasZero(scip, vals[i]));
         vardata = SCIPvarGetData(vars[i]);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         assert( (vardata->data.origvardata.pricingvar == NULL) == (vardata->blocknr == -1) );
         if ( vardata->blocknr == -1 )
         {
            SCIP_CALL( GCGcreatePricingVar(scip, &var, vars[i], pricingprobnr) );
         }
         assert(vardata->data.origvardata.pricingvar != NULL);
         SCIP_CALL( SCIPaddCoefLinear(probdata->pricingprobs[pricingprobnr], cons, vardata->data.origvardata.pricingvar, vals[i]) );

      }
      
      //SCIP_CALL( SCIPprintCons(probdata->pricingprobs[pricingprobnr], cons, NULL) );
      
   }

   return SCIP_OKAY;

}


/** creates a variable of the original program */
SCIP_RETCODE GCGcreateOrigVar(
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
   SCIP_VARDATA* vardata;

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);
   assert(var != NULL);
   assert(lb <= ub);

   SCIP_CALL( SCIPallocBlockMemory(probdata->origprob, &vardata) );
   vardata->vartype = GCG_VARTYPE_ORIGINAL;
   vardata->blocknr = -1;
   vardata->data.origvardata.pricingvar = NULL;
   vardata->data.origvardata.ncoefs = 0;
   vardata->data.origvardata.nmastervars = 0;
   vardata->data.origvardata.maxmastervars = STARTMAXMASTERVARS;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.origvardata.mastervars), 
         vardata->data.origvardata.maxmastervars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.origvardata.mastervals), 
         vardata->data.origvardata.maxmastervars) );

   SCIP_CALL( SCIPcreateVar(probdata->origprob, var, name, lb, ub, obj, vartype, initial, removable, 
         gcgvardelorig, NULL, NULL, vardata) );

   return SCIP_OKAY;
}

SCIP_RETCODE GCGchgOrigVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_VARDATA* vardata;

   probdata = SCIPgetProbData(scip);

   vardata = SCIPvarGetData(var);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   
   if ( vardata->data.origvardata.pricingvar != NULL )
   {
      SCIP_CALL( SCIPchgVarUb(probdata->pricingprobs[vardata->blocknr], 
            vardata->data.origvardata.pricingvar, newbound) );
   }

   SCIP_CALL( SCIPchgVarUb(probdata->origprob, var, newbound) );

   return SCIP_OKAY;
}


SCIP_RETCODE GCGchgOrigVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_VARDATA* vardata;

   probdata = SCIPgetProbData(scip);

   vardata = SCIPvarGetData(var);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);

   if ( vardata->data.origvardata.pricingvar != NULL )
   {
      SCIP_CALL( SCIPchgVarLb(probdata->pricingprobs[vardata->blocknr], 
            vardata->data.origvardata.pricingvar, newbound) );
   }

   SCIP_CALL( SCIPchgVarLb(probdata->origprob, var, newbound) );

   return SCIP_OKAY;
}


SCIP_RETCODE GCGchgOrigVarType(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_VARTYPE          vartype             /**< new type of variable */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_VARDATA* vardata;

   probdata = SCIPgetProbData(scip);

   vardata = SCIPvarGetData(var);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);

   if ( vardata->data.origvardata.pricingvar != NULL )
   {
      SCIP_CALL( SCIPchgVarType(probdata->pricingprobs[vardata->blocknr], 
            vardata->data.origvardata.pricingvar, vartype) );
   }

   SCIP_CALL( SCIPchgVarType(probdata->origprob, var, vartype) );

   return SCIP_OKAY;
}



/** creates a variable of a pricing problem program */
SCIP_RETCODE GCGcreatePricingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   SCIP_VAR*             origvar,            /**< corresponding variable in the original program */
   int                   pricingprobnr       /**< number of the pricing problem to which the variable belongs */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_VARDATA* vardata;
   SCIP_VARDATA* origvardata;

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);
   assert(origvar != NULL);
   assert(pricingprobnr >= 0 && pricingprobnr < probdata->npricingprobs);
   
   origvardata = SCIPvarGetData(origvar);
   
   assert(origvardata != NULL);
   assert(origvardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(origvardata->data.origvardata.pricingvar == NULL);
   assert(origvardata->blocknr == -1);


   SCIP_CALL( SCIPallocBlockMemory(probdata->pricingprobs[pricingprobnr], &vardata) );
   vardata->vartype = GCG_VARTYPE_PRICING;
   vardata->blocknr = pricingprobnr;
   vardata->data.pricingvardata.origvar = origvar;

   //printf("var = %s, ub = %f, type = %d\n", SCIPvarGetName(origvar), SCIPvarGetUbGlobal(origvar), SCIPvarGetType(origvar));

   SCIP_CALL( SCIPcreateVar(probdata->pricingprobs[pricingprobnr], var, SCIPvarGetName(origvar), 
         SCIPvarGetLbGlobal(origvar), SCIPvarGetUbGlobal(origvar), 0, SCIPvarGetType(origvar), 
         TRUE, FALSE, gcgvardelorig, NULL, NULL, vardata) );

   origvardata->data.origvardata.pricingvar = *var;
   origvardata->blocknr = pricingprobnr;

   SCIP_CALL( SCIPaddVar(probdata->pricingprobs[pricingprobnr], *var) );

   /* because the variable was added to the problem, it is captured by SCIP and we can safely release it right now
    * without making the returned *var invalid
    */
   SCIP_CALL( SCIPreleaseVar(probdata->pricingprobs[pricingprobnr], var) );
   
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
   SCIP_CALL( SCIPreleaseVar(probdata->origprob, &var) );

   return SCIP_OKAY;
}


SCIP* GCGprobGetOrigprob(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   
   probdata = SCIPgetProbData(scip);
   
   assert(probdata != NULL);
   
   return probdata->origprob;
}

SCIP* GCGprobGetPricingprob(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   )
{
   SCIP_PROBDATA* probdata;
   
   assert(scip != NULL);
   
   probdata = SCIPgetProbData(scip);
   
   assert(probdata != NULL);

   return probdata->pricingprobs[pricingprobnr];
}

int GCGprobGetNPricingprobs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   
   probdata = SCIPgetProbData(scip);
   
   assert(probdata != NULL);
   
   return probdata->npricingprobs;
}

void GCGprobGetMasterConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          conss,
   int*                  nconss
   )
{
   SCIP_PROBDATA* probdata;
   
   probdata = SCIPgetProbData(scip);
   
   assert(probdata != NULL);
   assert(nconss != NULL);
   assert(probdata->masterconss != NULL);
   assert(probdata->nmasterconss >= 0);

   *conss = probdata->masterconss;   
   *nconss = probdata->nmasterconss;
}


int GCGprobGetNMasterConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   
   probdata = SCIPgetProbData(scip);
   
   assert(probdata != NULL);
   assert(probdata->nmasterconss >= 0);

   return probdata->nmasterconss;
}



void GCGprobGetOrigMasterConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          conss,
   int*                  nconss
   )
{
   SCIP_PROBDATA* probdata;
   
   probdata = SCIPgetProbData(scip);
   
   assert(probdata != NULL);
   assert(nconss != NULL);
   assert(probdata->origmasterconss != NULL);
   assert(probdata->norigmasterconss >= 0);

   *conss = probdata->origmasterconss;   
   *nconss = probdata->norigmasterconss;

}

SCIP_CONS* GCGprobGetConvCons(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   pricingprobnr
   )
{
   SCIP_PROBDATA* probdata;
   
   probdata = SCIPgetProbData(scip);
   
   assert(probdata != NULL);
   assert(probdata->convconss != NULL);
   assert(pricingprobnr >= 0 && pricingprobnr < probdata->npricingprobs);


   return probdata->convconss[pricingprobnr];
}

