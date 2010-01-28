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
//#define SCIP_DEBUG
//#define DEBUG_PRICING
//#define DEBUG_PRICING_ALL_OUTPUT
//#define CHECKNEWVAR
//#define CHECKVARBOUNDS
/**@file   pricer_gcg.c
 * @ingroup PRICERS
 * @brief  pricer for generic column generation, solves the pricing problem as a MIP
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "pricer_gcg.h"
#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"
#include "scip/scip.h"
#include "sepa_master.h"
#include "struct_vardata.h"
#include "struct_solver.h"
#include "type_solver.h"
#include <stdio.h>
#include <stdlib.h>


#define PRICER_NAME            "gcg"
#define PRICER_DESC            "pricer for gcg"
#define PRICER_PRIORITY        5000000
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

#define EPS 0.0001

#define DEFAULT_MAXVARSROUNDFARKAS 1
#define DEFAULT_MAXVARSROUNDREDCOST 100
#define DEFAULT_MAXSUCCESSFULMIPSREDCOST INT_MAX
#define DEFAULT_MAXROUNDSREDCOST INT_MAX
#define DEFAULT_MAXSOLSPROB INT_MAX
#define DEFAULT_USEHEURPRICING FALSE
#define DEFAULT_ONLYPOSCONV FALSE
#define DEFAULT_ABORTPRICING TRUE
#define DEFAULT_ONLYBEST FALSE
#define DEFAULT_SUCCESSFULMIPSREL 1.0
#define DEFAULT_DISPINFOS FALSE
#define DEFAULT_SORTING 0

#define MAXBEST 1000


/*
 * Data structures
 */


/** variable pricer data */
struct SCIP_PricerData
{
   int npricingprobs;              /* number of pricing problems */
   SCIP** pricingprobs;            /* pointers to the pricing problems */
   SCIP_Real* dualsolconv;         /* array of dual solutions for the convexity constraints */
   SCIP* origprob;                 /* the original program */
   SCIP_Real* solvals;             /* solution values of variables in the pricing problems */
   int* nvarsprob;                 /* number of variables created by the pricing probs */
   SCIP_Longint currnodenr;
   SCIP_HASHMAP* mapcons2idx;
   SCIP_Real* tmpconvdualsol;
   int* permu;
   int npricingprobsnotnull;

   SCIP_Real** bestsolvals;
   SCIP_VAR*** bestsolvars;
   int* nbestsolvars;
   int* prob;
   SCIP_Real* redcost;
   int nbestsols;
   int maxbestsols;
   int maxvars;

   SCIP_VAR** pricedvars;
   int npricedvars;
   int maxpricedvars;

   SCIP_Real probfactor;

   /** variables used for statistics */
   SCIP_CLOCK* redcostclock;
   SCIP_CLOCK* redcostsolveclock;
   SCIP_CLOCK* farkasclock;
   SCIP_CLOCK* farkassolveclock;
   SCIP_CLOCK* freeclock;
   SCIP_CLOCK* transformclock;
   int solvedsubmipsoptimal;
   int solvedsubmipsheur;
   int calls;
   int farkascalls;
   int redcostcalls;

   /* solver data */
   GCG_SOLVER**     solvers;
   int              nsolvers;

   /** parameter values */
   SCIP_VARTYPE vartype;           /* vartype of created master variables */
   int maxvarsroundfarkas;
   int maxvarsroundredcost;
   int maxsuccessfulmipsredcost;
   int maxroundsredcost;
   int maxsolsprob;
   int nroundsredcost;
   int sorting;
   SCIP_Bool useheurpricing;
   SCIP_Bool onlyposconv;
   SCIP_Bool abortpricing;
   SCIP_Bool onlybest;
   SCIP_Bool dispinfos;
   SCIP_Real successfulmipsrel;
};



/*
 * Vardata methods
 */

static
SCIP_DECL_VARDELTRANS(gcgvardeltrans)
{
   assert((*vardata)->vartype == GCG_VARTYPE_MASTER);
   SCIPfreeBlockMemoryArray(scip, &((*vardata)->data.mastervardata.origvals), (*vardata)->data.mastervardata.norigvars);
   SCIPfreeBlockMemoryArray(scip, &((*vardata)->data.mastervardata.origvars), (*vardata)->data.mastervardata.norigvars);
   
   SCIPfreeBlockMemory(scip, vardata);

   return SCIP_OKAY;
}

static
SCIP_DECL_PARAMCHGD(paramChgdOnlybestMaxvars)
{
   SCIP_PARAMDATA* paramdata;
   SCIP_PRICERDATA* pricerdata;
   int i;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);
   pricerdata = (SCIP_PRICERDATA*) paramdata;


   if( SCIPgetStage(scip) <= SCIP_STAGE_PRESOLVED )
      return SCIP_OKAY;

   /* free array if not needed anymore */
   if ( !pricerdata->onlybest && pricerdata->maxbestsols > 0 )
   {
      assert(pricerdata->bestsolvars != NULL);
      assert(pricerdata->bestsolvals != NULL);
      assert(pricerdata->nbestsolvars != NULL);
      assert(pricerdata->redcost != NULL );
      assert(pricerdata->prob != NULL );

      for( i = 0; i < pricerdata->maxbestsols; i++ )
      {
         assert(pricerdata->bestsolvars[i] != NULL);
         assert(pricerdata->bestsolvals[i] != NULL);
         
         SCIPfreeMemoryArray(scip, &(pricerdata->bestsolvars[i]));
         SCIPfreeMemoryArray(scip, &(pricerdata->bestsolvals[i]));
      }

      SCIPfreeMemoryArray(scip, &pricerdata->bestsolvars);
      SCIPfreeMemoryArray(scip, &pricerdata->bestsolvals);
      SCIPfreeMemoryArray(scip, &pricerdata->nbestsolvars);
      SCIPfreeMemoryArray(scip, &pricerdata->redcost);
      SCIPfreeMemoryArray(scip, &pricerdata->prob);

      pricerdata->bestsolvars = NULL;
      pricerdata->bestsolvals = NULL;
      pricerdata->nbestsolvars = NULL;
      pricerdata->maxbestsols = 0;
      pricerdata->nbestsols = 0;
   }

   /* create array */
   if ( pricerdata->onlybest && pricerdata->maxbestsols == 0 && pricerdata->maxvarsroundredcost <= MAXBEST)
   {
      assert(pricerdata->bestsolvars == NULL);
      assert(pricerdata->bestsolvals == NULL);
      assert(pricerdata->nbestsolvars == NULL);
      assert(pricerdata->redcost == NULL);
      assert(pricerdata->prob == NULL);

      pricerdata->maxbestsols = pricerdata->maxvarsroundredcost;

      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->bestsolvars, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->bestsolvals, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->nbestsolvars, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->redcost, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->prob, pricerdata->maxbestsols) );

      for( i = 0; i < pricerdata->maxbestsols; i++ )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->bestsolvars[i]), pricerdata->maxvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->bestsolvals[i]), pricerdata->maxvars) );
         pricerdata->nbestsolvars[i] = 0;
      }

      pricerdata->nbestsols = 0;
   }

   /* change size of array */
   if ( pricerdata->onlybest && pricerdata->maxbestsols != 0 )
   {
      assert(pricerdata->bestsolvars != NULL);
      assert(pricerdata->bestsolvals != NULL);
      assert(pricerdata->nbestsolvars != NULL);
      assert(pricerdata->redcost != NULL);
      assert(pricerdata->prob != NULL);

      SCIP_CALL( SCIPreallocMemoryArray(scip, &pricerdata->bestsolvars, pricerdata->maxvarsroundredcost) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &pricerdata->bestsolvals, pricerdata->maxvarsroundredcost) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &pricerdata->nbestsolvars, pricerdata->maxvarsroundredcost) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &pricerdata->redcost, pricerdata->maxvarsroundredcost) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &pricerdata->prob, pricerdata->maxvarsroundredcost) );

      for( i = pricerdata->maxbestsols; i < pricerdata->maxvarsroundredcost; i++ )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->bestsolvars[i]), pricerdata->maxvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->bestsolvals[i]), pricerdata->maxvars) );
         pricerdata->nbestsolvars[i] = 0;
      }

      pricerdata->maxbestsols = pricerdata->maxvarsroundredcost;
   }

   printf("paramchanged\n");

   return SCIP_OKAY;
}



/*
 * Local methods
 */

/* ensures size of pricedvars array */
static
SCIP_RETCODE ensureSizePricedvars(
   SCIP*                 scip,
   SCIP_PRICERDATA*      pricerdata,
   int                   size
   )
{
   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert(pricerdata->pricedvars != NULL);

   if ( pricerdata->maxpricedvars < size )
   {
      pricerdata->maxpricedvars = MAX( 2 * pricerdata->maxpricedvars, size );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(pricerdata->pricedvars), pricerdata->maxpricedvars) );
   }
   assert(pricerdata->maxpricedvars >= size);

   return SCIP_OKAY;
}


/* ensures size of solvers array */
static
SCIP_RETCODE ensureSizeSolvers(
   SCIP*                 scip,
   SCIP_PRICERDATA*      pricerdata
   )
{
   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));   

   if ( pricerdata->nsolvers == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->solvers), 1) );
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(pricerdata->solvers), pricerdata->nsolvers+1) );
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE solversFree(
   SCIP*                 scip,
   SCIP_PRICERDATA*      pricerdata
   )
{
   int i;

   assert(scip != NULL);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));   
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverfree != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverfree(scip, pricerdata->solvers[i]) );

         BMSfreeMemoryArray(&pricerdata->solvers[i]->name);
         BMSfreeMemoryArray(&pricerdata->solvers[i]->description);
   
         SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->solvers[i]->optfarkasclock)) );
         SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->solvers[i]->optredcostclock)) );
         SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->solvers[i]->heurfarkasclock)) );
         SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->solvers[i]->heurredcostclock)) );

         SCIPfreeMemory(scip, &(pricerdata->solvers[i]));
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE solversInit(
   SCIP*                 scip,
   SCIP_PRICERDATA*      pricerdata
   )
{
   int i;

   assert(scip != NULL);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));   
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverinit != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverinit(scip, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE solversExit(
   SCIP*                 scip,
   SCIP_PRICERDATA*      pricerdata
   )
{
   int i;

   assert(scip != NULL);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));   
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverexit != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverexit(scip, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE solversInitsol(
   SCIP*                 scip,
   SCIP_PRICERDATA*      pricerdata
   )
{
   int i;

   assert(scip != NULL);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));   
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverinitsol != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverinitsol(scip, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE solversExitsol(
   SCIP*                 scip,
   SCIP_PRICERDATA*      pricerdata
   )
{
   int i;

   assert(scip != NULL);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));   
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverexitsol != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverexitsol(scip, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}




static
SCIP_RETCODE solvePricingProblem(
   SCIP*                 scip,
   SCIP_PRICERDATA*      pricerdata,
   int                   prob,
   GCG_PRICETYPE         pricetype,
   SCIP_VAR****          solvars,
   SCIP_Real***          solvals,
   int**                 nsolvars,
   int*                  nsols,
   SCIP_STATUS*          status
   )
{
   int i;

   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert(pricerdata->pricingprobs[prob] != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));   
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solversolve != NULL )
      {
         if( pricetype == GCG_PRICETYPE_FARKAS )
         {
            SCIP_CALL( SCIPstartClock(scip, pricerdata->solvers[i]->optfarkasclock));
         }
         else
         {
            SCIP_CALL( SCIPstartClock(scip, pricerdata->solvers[i]->optredcostclock));
         }

         SCIP_CALL( pricerdata->solvers[i]->solversolve(scip, pricerdata->solvers[i], 
               pricerdata->pricingprobs[prob], prob, 
               solvars, solvals, nsolvars, nsols, status) );

         if( pricetype == GCG_PRICETYPE_FARKAS )
         {
            SCIP_CALL( SCIPstopClock(scip, pricerdata->solvers[i]->optfarkasclock));
            if( *status != SCIP_STATUS_UNKNOWN )
               pricerdata->solvers[i]->optfarkascalls++;
         }
         else
         {
            SCIP_CALL( SCIPstopClock(scip, pricerdata->solvers[i]->optredcostclock));
            if( *status != SCIP_STATUS_UNKNOWN )
               pricerdata->solvers[i]->optredcostcalls++;
         }

         if( *status == SCIP_STATUS_OPTIMAL )
            break;
      }

   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE solvePricingProblemHeur(
   SCIP*                 scip,
   SCIP_PRICERDATA*      pricerdata,
   int                   prob,
   GCG_PRICETYPE         pricetype,
   SCIP_VAR****          solvars,
   SCIP_Real***          solvals,
   int**                 nsolvars,
   int*                  nsols,
   SCIP_STATUS*          status
   )
{
   int i;

   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert(pricerdata->pricingprobs[prob] != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));   
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solversolve != NULL )
      {
         if( pricetype == GCG_PRICETYPE_FARKAS )
         {
            SCIP_CALL( SCIPstartClock(scip, pricerdata->solvers[i]->heurfarkasclock));
         }
         else
         {
            SCIP_CALL( SCIPstartClock(scip, pricerdata->solvers[i]->heurredcostclock));
         }

         SCIP_CALL( pricerdata->solvers[i]->solversolveheur(scip, pricerdata->solvers[i], 
               pricerdata->pricingprobs[prob], prob, 
               solvars, solvals, nsolvars, nsols, status) );


         if( pricetype == GCG_PRICETYPE_FARKAS )
         {
            SCIP_CALL( SCIPstopClock(scip, pricerdata->solvers[i]->heurfarkasclock));
            if( *status != SCIP_STATUS_UNKNOWN )
               pricerdata->solvers[i]->heurfarkascalls++;
         }
         else
         {
            SCIP_CALL( SCIPstopClock(scip, pricerdata->solvers[i]->heurredcostclock));
            if( *status != SCIP_STATUS_UNKNOWN )
               pricerdata->solvers[i]->heurredcostcalls++;
         }

         if( *status == SCIP_STATUS_OPTIMAL )
            break;
      }

      if( *status == SCIP_STATUS_OPTIMAL )
         break;
   }

   return SCIP_OKAY;
}


#ifdef CHECKNEWVAR
static
SCIP_RETCODE checkNewVar(
   SCIP*                 scip,
   SCIP_VAR*             newvar,
   SCIP_Real             redcost,
   SCIP_Real             dualsolconv
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;
   int i;

   SCIP_VARDATA* newvardata;
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(newvar != NULL);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* compare newvar with all existing variables */
   for( v = 0; v < nvars; v++ )
   {
      assert(vars[v] != NULL);

      /* vars[v] is the new variable itself */
      if( vars[v] == newvar )
         continue;

      /* vars[v] has a different objective function value, may not be equal to newvar */
      if( !SCIPisEQ(scip, SCIPvarGetObj(vars[v]), SCIPvarGetObj(newvar)) )
         continue;

      vardata = SCIPvarGetData(vars[v]);
      newvardata = SCIPvarGetData(newvar);
      assert(vardata != NULL);
      assert(newvardata != NULL);

      /* vars[v] belongs to a different block, may not be equal to newvar */
      if( vardata->blocknr != newvardata->blocknr )
         continue;

      assert(vardata->vartype == GCG_VARTYPE_MASTER);
      assert(newvardata->vartype == GCG_VARTYPE_MASTER);

      /* vars[v] belongs to a different block, may not be equal to newvar */
      if( vardata->data.mastervardata.norigvars != newvardata->data.mastervardata.norigvars )
         continue;

      /* compare the parts of the original variables contained in vars[i] and newvar */
      for( i = 0; i < vardata->data.mastervardata.norigvars; i++ )
      {
         /* original variables are not equal */
         if( vardata->data.mastervardata.origvars[i] != newvardata->data.mastervardata.origvars[i] )
            break;
         /* the values of the original variable are not equal */
         if( !SCIPisEQ(scip, vardata->data.mastervardata.origvals[i], newvardata->data.mastervardata.origvals[i]) )
            break;
      }

      if( i == vardata->data.mastervardata.norigvars )
      {
         printf("var %s is equal to var %s! solval = %f, ub = %g, lazyub = %g, redcost = %f, lpsolstat = %d, dualsolconv = %f\n", 
            SCIPvarGetName(newvar), SCIPvarGetName(vars[v]), SCIPgetSolVal(scip, NULL, vars[v]), 
            SCIPvarGetUbLocal(vars[v]), SCIPvarGetUbLazy(vars[v]), redcost, SCIPgetLPSolstat(scip), dualsolconv);
      }
   }

   return SCIP_OKAY;
}
#endif

#ifdef CHECKVARBOUNDS
static
SCIP_RETCODE checkVarBounds(
   SCIP*                 scip
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;
   SCIP_VARDATA* vardata;
   SCIP* origscip;

   assert(scip != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   SCIP_CALL( SCIPgetVarsData(origscip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* check whether the corresponding pricing MIP has the same bound for the variable */
   for( v = 0; v < nvars; v++ )
   {
      assert(vars[v] != NULL);

      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      assert(vardata->data.origvardata.pricingvar != NULL || vardata->blocknr == -1);

      if( vardata->blocknr == -1 )
         continue;

      if( !GCGrelaxIsPricingprobRelevant(origscip, vardata->blocknr) || GCGrelaxGetNIdenticalBlocks(origscip, vardata->blocknr) != 1 ) 
         continue;

      if( SCIPvarGetUbLocal(vars[v]) != SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar) )
      {
         printf("var %s: orig upper bound = %g, pricing upper bound = %g, global orig upper bound = %g!\n", 
            SCIPvarGetName(vars[v]), SCIPvarGetUbLocal(vars[v]), SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar),
            SCIPvarGetUbGlobal(vars[v]));
      }
      if( SCIPvarGetLbLocal(vars[v]) != SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar) )
      {
         printf("var %s: orig lower bound = %g, pricing lower bound = %g, global orig lower bound = %g!\n", 
            SCIPvarGetName(vars[v]), SCIPvarGetLbLocal(vars[v]), SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar),
            SCIPvarGetLbGlobal(vars[v]) );
      }

      assert(SCIPvarGetUbLocal(vars[v]) == SCIPvarGetUbLocal(vardata->data.origvardata.pricingvar));
      assert(SCIPvarGetLbLocal(vars[v]) == SCIPvarGetLbLocal(vardata->data.origvardata.pricingvar));
   }

   return SCIP_OKAY;
}
#endif

/* informs an original variable, that a variable in the master problem was created, 
 * that contains a part of the original variable.
 * Saves this information in the original variable's data */
SCIP_RETCODE GCGpricerAddMasterVarToOrigVar(
   SCIP*                 scip,
   SCIP_VAR*             origvar,
   SCIP_VAR*             var,
   SCIP_Real             val
   )
{
   SCIP_VARDATA* vardata;
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(origvar != NULL);
   assert(var != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);   

   vardata = SCIPvarGetData(origvar);
   assert(vardata != NULL);
   assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata->data.origvardata.mastervars != NULL);
   assert(vardata->data.origvardata.mastervals != NULL);
   assert(vardata->data.origvardata.nmastervars >= 0);
   assert(vardata->data.origvardata.maxmastervars >= vardata->data.origvardata.nmastervars);

   /* realloc mastervars array of the original variable, if needed */
   if( vardata->data.origvardata.maxmastervars == vardata->data.origvardata.nmastervars )
   {
      SCIP_CALL( SCIPreallocMemoryArray(pricerdata->origprob, &(vardata->data.origvardata.mastervars),
            2*vardata->data.origvardata.maxmastervars) );
      SCIP_CALL( SCIPreallocMemoryArray(pricerdata->origprob, &(vardata->data.origvardata.mastervals),
            2*vardata->data.origvardata.maxmastervars) );
      SCIPdebugMessage("mastervars array of var %s resized from %d to %d\n", SCIPvarGetName(origvar), 
         vardata->data.origvardata.maxmastervars, 2*vardata->data.origvardata.maxmastervars);
      vardata->data.origvardata.maxmastervars = 2*vardata->data.origvardata.maxmastervars;
   }
   /* add information to the original variable's vardata */
   vardata->data.origvardata.mastervars[vardata->data.origvardata.nmastervars] = var;
   vardata->data.origvardata.mastervals[vardata->data.origvardata.nmastervars] = val;
   vardata->data.origvardata.nmastervars++;

   return SCIP_OKAY;
}


static
SCIP_RETCODE setPricingObjs(
   SCIP*                 scip,
   GCG_PRICETYPE         pricetype
   )
{
   SCIP* origprob;
   SCIP_VARDATA* vardata;
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS** origconss;
   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_VAR** probvars;
   int nprobvars;

   SCIP_ROW** mastercuts;
   int nmastercuts;
   SCIP_ROW** origcuts;
   int norigcuts;
   SCIP_COL** cols;
   SCIP_Real* consvals;
   SCIP_Real dualsol;

   SCIP_VAR** consvars;
   int nconsvars;

   int i;
   int j;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   /* get the constraints of the master problem and the corresponding constraints in the original problem */
   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);
   origconss = GCGrelaxGetLinearOrigMasterConss(origprob);

   /* set objective value of all variables in the pricing problems to 0 (for farkas pricing) /
    * to the original objective of the variable (for redcost pricing) */
   for( i = 0; i < pricerdata->npricingprobs; i++)
   {
      if( pricerdata->pricingprobs[i] == NULL )
         continue;
      probvars = SCIPgetVars(pricerdata->pricingprobs[i]);
      nprobvars = SCIPgetNVars(pricerdata->pricingprobs[i]);

      for( j = 0; j < nprobvars; j++ )
      {
         if( pricetype == GCG_PRICETYPE_FARKAS )
         {
            SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], 0) );
         }
         else 
         {
            vardata = SCIPvarGetData(probvars[j]);
            assert(vardata->vartype == GCG_VARTYPE_PRICING);
            assert(vardata->blocknr == i);
            assert(vardata->data.pricingvardata.origvars != NULL);
            assert(vardata->data.pricingvardata.origvars[0] != NULL);
            SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], 
                  SCIPvarGetObj(vardata->data.pricingvardata.origvars[0])) );
         }  
      }
   }

   /* compute reduced cost and update objectives in the pricing problems */
   for( i = 0; i < nmasterconss; i++ )
   {
      /* farkas pricing */
      if( pricetype == GCG_PRICETYPE_REDCOST )
         dualsol = SCIPgetDualsolLinear(scip, masterconss[i]);
      /* redcost pricing */
      else 
      {
         assert(pricetype == GCG_PRICETYPE_FARKAS);
         dualsol = SCIPgetDualfarkasLinear(scip, masterconss[i]);
      }
      if( !SCIPisZero(scip, dualsol) )
      {
         /* for all variables in the constraint, modify the objective of the corresponding variable in a pricing problem */
         consvars = SCIPgetVarsLinear(origprob, origconss[i]);
         consvals = SCIPgetValsLinear(origprob, origconss[i]);
         nconsvars = SCIPgetNVarsLinear(origprob, origconss[i]);
         for( j = 0; j < nconsvars; j++ )
         {
            vardata = SCIPvarGetData(consvars[j]);
            assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
            if( vardata->blocknr != -1 && pricerdata->pricingprobs[vardata->blocknr] != NULL )
            {
               assert(vardata->data.origvardata.pricingvar != NULL);
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[vardata->blocknr], 
                     vardata->data.origvardata.pricingvar, -1.0 * dualsol * consvals[j]) );
            }
         }
      }
   }

   /* get the cuts of the master problem and the corresponding cuts in the original problem */
   mastercuts = GCGsepaGetMastercuts(scip);
   nmastercuts = GCGsepaGetNMastercuts(scip);
   origcuts = GCGsepaGetOrigcuts(scip);
   norigcuts = GCGsepaGetNOrigcuts(scip);
   
   assert(mastercuts != NULL);
   assert(origcuts != NULL);
   assert(norigcuts == nmastercuts);

   /* compute reduced cost and update objectives in the pricing problems */
   for( i = 0; i < nmastercuts; i++ )
   {
      /* farkas pricing */
      if( pricetype == GCG_PRICETYPE_REDCOST )
         dualsol = SCIProwGetDualsol(mastercuts[i]);
      /* redcost pricing */
      else
      {    
         assert(pricetype == GCG_PRICETYPE_FARKAS);
         dualsol = SCIProwGetDualfarkas(mastercuts[i]);
      }
      if( !SCIPisZero(scip, dualsol) )
      {
         /* get columns and vals of the cut */
         nconsvars = SCIProwGetNNonz(origcuts[i]);
         cols = SCIProwGetCols(origcuts[i]);
         consvals = SCIProwGetVals(origcuts[i]);

         /* get the variables corresponding to the columns in the cut */
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
         for( j = 0; j < nconsvars; j++ )
            consvars[j] = SCIPcolGetVar(cols[j]);

         /* for all variables in the cut, modify the objective of the corresponding variable in a pricing problem */
         for( j = 0; j < nconsvars; j++ )
         {
            vardata = SCIPvarGetData(consvars[j]);
            assert(vardata != NULL);
            assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
            if( vardata->blocknr != -1 && pricerdata->pricingprobs[vardata->blocknr] != NULL )
            {
               assert(vardata->data.origvardata.pricingvar != NULL);
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[vardata->blocknr], 
                     vardata->data.origvardata.pricingvar, -1.0 * dualsol * consvals[j]) );
            }
         }
         SCIPfreeBufferArray(scip, &consvars);
      }
   }

   /* get dual solutions / farkas values of the convexity constraints */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {

      assert( GCGrelaxIsPricingprobRelevant(origprob, i) 
         == (GCGrelaxGetConvCons(origprob, i) != NULL) );
      if( !GCGrelaxIsPricingprobRelevant(origprob, i) )
      {
         pricerdata->dualsolconv[i] = -1.0 * SCIPinfinity(scip);
         continue;
      }  
      if( pricetype == GCG_PRICETYPE_REDCOST )
         pricerdata->dualsolconv[i] = SCIPgetDualsolLinear(scip, GCGrelaxGetConvCons(origprob, i));
      else
      {
         assert(pricetype == GCG_PRICETYPE_FARKAS);
         pricerdata->dualsolconv[i] = SCIPgetDualfarkasLinear(scip, GCGrelaxGetConvCons(origprob, i));
      }
   }

   return SCIP_OKAY;
}

/** creates a new master variable corresponding to the given solution and problem */
static
SCIP_RETCODE createNewMasterVar(
   SCIP*                 scip,
   SCIP_VAR**            solvars,
   SCIP_Real*            solvals,
   int                   nsolvars,
   int                   prob,
   SCIP_Bool             checkonlybest,
   SCIP_Bool*            added
   )
{
   SCIP* origprob;
   SCIP_VARDATA* vardata;
   SCIP_VARDATA* newvardata;
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_Real* mastercoefs;
   char varname[SCIP_MAXSTRLEN];

   SCIP_ROW** mastercuts;
   int nmastercuts;
   SCIP_ROW** origcuts;
   int norigcuts;
   SCIP_COL** cols;
   SCIP_Real conscoef;
   SCIP_VAR* var;
   SCIP_Real* consvals;
   SCIP_Real objcoeff;
   SCIP_VAR* newvar;

   SCIP_Real objvalue;

   SCIP_CONS* linkcons;
   int c;
   int idx;

   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(solvars != NULL);
   assert(solvals != NULL);
   assert(nsolvars >= 0);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);
   
   objvalue = 0.0;

   /* for integer variable, round the solvals */
   for( i = 0; i < nsolvars; i++ )
   {
      if( SCIPvarGetType(solvars[i]) != SCIP_VARTYPE_CONTINUOUS )
      {
         assert(SCIPisEQ(scip, solvals[i], SCIPfloor(scip, solvals[i])));
         solvals[i] = SCIPfloor(scip, solvals[i]);
      }
      objvalue += solvals[i] * SCIPvarGetObj(solvars[i]);
   }

   if( SCIPisRelGE(scip, objvalue, pricerdata->dualsolconv[prob]) )
   {
      *added = FALSE;

      return SCIP_OKAY;
   }
   else
   {
      *added = TRUE;
   }

   SCIPdebugMessage("found var with redcost %g (objvalue = %g, dualsol =%g)\n", 
      objvalue - pricerdata->dualsolconv[prob], objvalue, pricerdata->dualsolconv[prob]);
   //printf("found var with redcost %f (objvalue = %g, dualsol =%g), checkonlybest = %d\n", 
   //   objvalue - pricerdata->dualsolconv[prob], objvalue, pricerdata->dualsolconv[prob], checkonlybest);

   if( checkonlybest && pricerdata->onlybest && pricerdata->maxbestsols > 0 )
   {
      int pos;

      for( pos = pricerdata->nbestsols - 1; pos >= 0 && pricerdata->redcost[pos] > 
              objvalue - pricerdata->dualsolconv[prob]; pos-- )
      {
         if( pos < pricerdata->maxbestsols - 1 )
         {
            pricerdata->prob[pos+1] = pricerdata->prob[pos];
            pricerdata->redcost[pos+1] = pricerdata->redcost[pos];
            pricerdata->nbestsolvars[pos+1] = pricerdata->nbestsolvars[pos];
            for( i = 0; i < pricerdata->nbestsolvars[pos]; i++ )
            {
               pricerdata->bestsolvars[pos+1][i] = pricerdata->bestsolvars[pos][i];
               pricerdata->bestsolvals[pos+1][i] = pricerdata->bestsolvals[pos][i];
            }
         }
         else
         {
            pricerdata->nbestsols--;
         }
      }
      pos++;

      if( pos != pricerdata->maxbestsols )
      {
         pricerdata->prob[pos] = prob;
         pricerdata->redcost[pos] = objvalue - pricerdata->dualsolconv[prob];
         pricerdata->nbestsolvars[pos] = nsolvars;
         
         for( i = 0; i < pricerdata->nbestsolvars[pos]; i++ )
         {
            pricerdata->bestsolvars[pos][i] = solvars[i];
            pricerdata->bestsolvals[pos][i] = solvals[i];
         }
         pricerdata->nbestsols++;
      }

      return SCIP_OKAY;
   }

            
   /* create data for the new variable in the master problem */
   SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
   newvardata->vartype = GCG_VARTYPE_MASTER;
   newvardata->blocknr = prob;

   /* compute objective coefficient of the variable */
   objcoeff = 0;
   for( i = 0; i < nsolvars; i++ )
   {
      if( !SCIPisZero(scip, solvals[i]) )
      {
         vardata = SCIPvarGetData(solvars[i]);
         assert(vardata->vartype == GCG_VARTYPE_PRICING);
         assert(vardata->data.pricingvardata.origvars != NULL);
         assert(vardata->data.pricingvardata.origvars[0] != NULL);
         /* add quota of original variable's objcoef to the master variable's coef */
         objcoeff += solvals[i] * SCIPvarGetObj(vardata->data.pricingvardata.origvars[0]);
      }
   }

   (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "p_%d_%d", prob, pricerdata->nvarsprob[prob]);
   pricerdata->nvarsprob[prob]++;
   
   /* create variable in the master problem */
   SCIP_CALL( SCIPcreateVar(scip, &newvar, varname, 0, INT_MAX /*GCGrelaxGetNIdenticalBlocks(origprob, prob)*/, 
         objcoeff, pricerdata->vartype, TRUE, TRUE, NULL, NULL, gcgvardeltrans, newvardata) );

   SCIPdebugMessage("found var %s with redcost %f!\n", SCIPvarGetName(newvar), 
      objvalue - pricerdata->dualsolconv[prob]);

   /* count number of non-zeros */
   newvardata->data.mastervardata.norigvars = 0;

   for( i = 0; i < nsolvars; i++ )
      if( !SCIPisZero(scip, solvals[i]) )
         newvardata->data.mastervardata.norigvars++;

   if( newvardata->data.mastervardata.norigvars > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvars), newvardata->data.mastervardata.norigvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(newvardata->data.mastervardata.origvals), newvardata->data.mastervardata.norigvars) );
   }
   else
   {
      newvardata->data.mastervardata.origvars = NULL;
      newvardata->data.mastervardata.origvals = NULL;
   }
               
   /* number of original variables yet saved in mastervardata */
   j = 0;

   /* update variable datas */
   for( i = 0; i < nsolvars; i++ )
   {
      if( !SCIPisZero(scip, solvals[i]) )
      {
         vardata = SCIPvarGetData(solvars[i]);
         assert(vardata->vartype == GCG_VARTYPE_PRICING);
         assert(vardata->data.pricingvardata.origvars != NULL);
         assert(vardata->data.pricingvardata.origvars[0] != NULL);
         /* save in the master problem variable's data the quota of the corresponding original variable */
         newvardata->data.mastervardata.origvars[j] = vardata->data.pricingvardata.origvars[0];
         newvardata->data.mastervardata.origvals[j] = solvals[i];
         /* save the quota in the original variable's data */
         SCIP_CALL( GCGpricerAddMasterVarToOrigVar(scip, vardata->data.pricingvardata.origvars[0], newvar, solvals[i]) );
         j++;
      }
   }
   assert(j == newvardata->data.mastervardata.norigvars);

   /* add variable */
   SCIP_CALL( SCIPaddPricedVar(scip, newvar, pricerdata->dualsolconv[prob] - objvalue) );

   SCIP_CALL( SCIPcaptureVar(scip, newvar) );
   SCIP_CALL( ensureSizePricedvars(scip, pricerdata, pricerdata->npricedvars + 1) );
   pricerdata->pricedvars[pricerdata->npricedvars] = newvar;
   pricerdata->npricedvars++;

   SCIP_CALL( SCIPallocBufferArray(scip, &mastercoefs, nmasterconss) );
   BMSclearMemoryArray(mastercoefs, nmasterconss);

   /* compute coef of the variable in the master constraints */
   for( i = 0; i < nsolvars; i++ )
   {
      if( !SCIPisZero(scip, solvals[i]) )
      {
         vardata = SCIPvarGetData(solvars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_PRICING);
         assert(vardata->data.pricingvardata.origvars != NULL);
         assert(vardata->data.pricingvardata.origvars[0] != NULL);
         vardata = SCIPvarGetData(vardata->data.pricingvardata.origvars[0]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         assert(vardata->data.origvardata.coefs != NULL || vardata->data.origvardata.ncoefs == 0);
                  
         /* for each coef, add coef * solval to the coef of the new variable for the corresponding constraint */
         for( c = 0; c < vardata->data.origvardata.ncoefs; c++ )
         {
            assert(!SCIPisZero(scip, vardata->data.origvardata.coefs[c]));
            SCIP_CALL( SCIPgetTransformedCons(scip, vardata->data.origvardata.linkconss[c], 
                  &linkcons) );

            idx = (int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, linkcons);
            assert(0 <= idx && idx < nmasterconss);
            assert(masterconss[idx] == linkcons);
            mastercoefs[idx] += vardata->data.origvardata.coefs[c] * solvals[i];
         }
                  
      }
   }
   /* add the variables to the master constraints */
   for( i = 0; i < nmasterconss; i++ )
   {
      if( !SCIPisZero(scip, mastercoefs[i]) )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, masterconss[i], newvar, mastercoefs[i]) );
      }
   }

   /* get the cuts of the master problem and the corresponding cuts in the original problem */
   mastercuts = GCGsepaGetMastercuts(scip);
   nmastercuts = GCGsepaGetNMastercuts(scip);
   origcuts = GCGsepaGetOrigcuts(scip);
   norigcuts = GCGsepaGetNOrigcuts(scip);

   assert(mastercuts != NULL);
   assert(origcuts != NULL);
   assert(norigcuts == nmastercuts);

   /* compute coef of the variable in the cuts and add it to the cuts */
   for( i = 0; i < nmastercuts; i++ )
   {
      /* get columns of the cut and their coefficients */
      cols = SCIProwGetCols(origcuts[i]);
      consvals = SCIProwGetVals(origcuts[i]);

      conscoef = 0;

      for( j = 0; j < SCIProwGetNNonz(origcuts[i]); j++ )
      {
         var = SCIPcolGetVar(cols[j]);
         vardata = SCIPvarGetData(var);

         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
         if( vardata->blocknr == prob )
         {
            assert(vardata->data.origvardata.pricingvar != NULL);

            for( k = 0; k < nsolvars; k++ )
               if( solvars[k] == vardata->data.origvardata.pricingvar )
               {
                  conscoef += ( consvals[j] * solvals[k] );
                  break;
               }
         }

      }
               
      if( !SCIPisZero(scip, conscoef) )
      {
         SCIP_CALL( SCIPaddVarToRow(scip , mastercuts[i], newvar, conscoef) );
         //printf("new variable has coef = %f in cut %s:\n", conscoef, SCIProwGetName(mastercuts[i]));
         //SCIP_CALL( SCIPprintRow(origprob, mastercuts[i], NULL) );
      }
   }

   /* add variable to convexity constraint */
   SCIP_CALL( SCIPaddCoefLinear(scip, GCGrelaxGetConvCons(origprob, prob), newvar, 1) );

   SCIPfreeBufferArray(scip, &mastercoefs);

#ifdef CHECKNEWVAR
   /* check whether the created variable already existed */
   SCIP_CALL( checkNewVar(scip, newvar, SCIPgetSolOrigObj(pricerdata->pricingprobs[prob], sol) - pricerdata->dualsolconv[prob], pricerdata->dualsolconv[prob]) );
#endif

   SCIPreleaseVar(scip, &newvar);

   return SCIP_OKAY;
}


/* performs the pricing routine, gets the type of pricing that should be done: farkas or redcost pricing */
static
SCIP_RETCODE performPricing(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_PRICER*          pricer,             /* the pricer */
   GCG_PRICETYPE         pricetype,          /* type of the pricing */
   SCIP_RESULT*          result,             /* result pointer */
   SCIP_Real*            lowerbound          /* lowerbound pointer */
   )
{
   SCIP_PRICERDATA* pricerdata;            /* the data of the pricer */
   SCIP* origprob;

   int i;
   int j;
   int prob;

   int nfoundvars;
   int nfoundvarsprob;
   int successfulmips;

   SCIP_Real timelimit;

   SCIP_Real bestsolval;
   SCIP_Real bestredcost;
   SCIP_Bool bestredcostvalid;
   SCIP_Bool added;

   SCIP_STATUS status;

   SCIP_VAR*** solvars;
   SCIP_Real** solvals;
   int* nsolvars;
   int nsols;

   SCIP_Real maxtmpconvdualsol;
   SCIP_Real sumtmpconvdualsol;

   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   assert(result != NULL || pricetype == GCG_PRICETYPE_FARKAS);
   assert(lowerbound != NULL || pricetype == GCG_PRICETYPE_FARKAS);

   if( pricerdata->dispinfos )
   {
      //if( pricetype == GCG_PRICETYPE_REDCOST || SCIPgetNVars(scip) % 50 == 0 )
         printf("nvars = %d, current lowerbound = %g, time = %f, node = %lld\n", SCIPgetNVars(scip), 
            SCIPgetLPObjval(scip), SCIPgetSolvingTime(scip), SCIPgetNNodes(scip));
   }

   /* check whether pricing can be aborted: if objective value is always integral
    * and the current node's current lowerbound rounded up equals the 
    * current lp objective value rounded up we don't need to continue pricing
    * since the best possible feasible solution must have at least this value
    */
   if( pricerdata->abortpricing && SCIPisObjIntegral(scip) && pricetype == GCG_PRICETYPE_REDCOST 
      && SCIPceil(scip, SCIPgetNodeDualbound(scip, SCIPgetCurrentNode(scip))) 
      == SCIPceil(scip, SCIPgetLPObjval(scip)) && SCIPgetNNodes(scip) > 1 )
   {
      if( pricerdata->dispinfos )
         printf("pricing aborted due to integral objective: node LB = %g, LP obj = %g\n", 
            SCIPgetNodeDualbound(scip, SCIPgetCurrentNode(scip)), SCIPgetLPObjval(scip));

      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( pricetype == GCG_PRICETYPE_REDCOST )
   {
      pricerdata->redcostcalls++;
      *result = SCIP_SUCCESS;
   }
   if( pricetype == GCG_PRICETYPE_FARKAS )
      pricerdata->farkascalls++;

   pricerdata->calls++;
   nfoundvars = 0;
   successfulmips = 0;

#ifdef CHECKVARBOUNDS
   SCIP_CALL( checkVarBounds(scip) );
#endif
   /* set objectives of the variables in the pricing sub-MIPs */
   SCIP_CALL( setPricingObjs(scip, pricetype) );

   maxtmpconvdualsol = 0;
   sumtmpconvdualsol = 0;

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if( GCGrelaxIsPricingprobRelevant(origprob, i) )
      {
         if( maxtmpconvdualsol < pricerdata->tmpconvdualsol[i] )
            maxtmpconvdualsol = pricerdata->tmpconvdualsol[i];
         if( maxtmpconvdualsol < - pricerdata->tmpconvdualsol[i] )
            maxtmpconvdualsol = - pricerdata->tmpconvdualsol[i];
         sumtmpconvdualsol += pricerdata->tmpconvdualsol[i];
      }
   }


   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      assert(GCGrelaxIsPricingprobRelevant(origprob, i)
         == (pricerdata->dualsolconv[i] != -1.0 * SCIPinfinity(scip)));

      pricerdata->permu[i] = i;
      pricerdata->tmpconvdualsol[i] = pricerdata->dualsolconv[i];

      if( pricerdata->sorting == 1 && pricetype == GCG_PRICETYPE_REDCOST )
      {
         pricerdata->tmpconvdualsol[i] = (1.5 * maxtmpconvdualsol + pricerdata->tmpconvdualsol[i]) / sqrt(pricerdata->nvarsprob[i] + 1);
      }
      if( pricerdata->sorting == 2 && pricetype == GCG_PRICETYPE_REDCOST )
      {
         pricerdata->tmpconvdualsol[i] = (0.5 * maxtmpconvdualsol + pricerdata->tmpconvdualsol[i]) / sqrt(pricerdata->nvarsprob[i] + 1);
      }
      if( pricerdata->sorting == 3 && pricetype == GCG_PRICETYPE_REDCOST )
      {
         pricerdata->tmpconvdualsol[i] = ( sumtmpconvdualsol/pricerdata->npricingprobsnotnull + pricerdata->tmpconvdualsol[i]) / sqrt(pricerdata->nvarsprob[i] + 1);
      }
   }
   SCIPsortDownRealInt(pricerdata->tmpconvdualsol, pricerdata->permu, pricerdata->npricingprobs);

   bestredcost = 0.0;
   bestredcostvalid = FALSE;

   if( pricerdata->useheurpricing )
   {
      SCIPdebugMessage("heuristcal pricing\n");
   
      /* solve the pricing MIPs heuristically and check whether solutions 
       * corresponding to variables with negative reduced costs where found 
       */
      for( i = 0; i < pricerdata->npricingprobs && (pricetype == GCG_PRICETYPE_FARKAS || ((pricerdata->onlybest ||
                  nfoundvars < pricerdata->maxvarsroundredcost) && successfulmips < pricerdata->maxsuccessfulmipsredcost
               && successfulmips < pricerdata->successfulmipsrel * pricerdata->npricingprobsnotnull))
              && (nfoundvars == 0 || pricerdata->dualsolconv[pricerdata->permu[i]] > 0 || !pricerdata->onlyposconv )
              && (pricetype == GCG_PRICETYPE_REDCOST || nfoundvars < pricerdata->maxvarsroundfarkas); i++)
      {
         prob = pricerdata->permu[i];

         if( pricerdata->pricingprobs[prob] == NULL )
            continue;

         /* set objective limit, such that only solutions with negative reduced costs are accepted */
         SCIP_CALL( SCIPsetObjlimit(pricerdata->pricingprobs[prob], pricerdata->dualsolconv[prob]) );

         /* set time limit */
         SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
         if( !SCIPisInfinity(scip, timelimit) )
         {
            if( timelimit - SCIPgetTotalTime(scip) > 0 )
            {
               SCIP_CALL( SCIPsetRealParam(pricerdata->pricingprobs[prob], "limits/time", 
                     timelimit - SCIPgetTotalTime(scip)) );
            }
            else
            {
               *result = SCIP_DIDNOTRUN;

               return SCIP_OKAY;
            }
         }

         pricerdata->solvedsubmipsheur++;

         SCIP_CALL( solvePricingProblemHeur(scip, pricerdata, prob, pricetype, &solvars, &solvals, 
               &nsolvars, &nsols, &status) );

         //printf("Pricingprob %d has found %d sols!\n", prob, nsols);

         nfoundvarsprob = 0;

         for( j = 0; j < nsols && nfoundvarsprob <= pricerdata->maxsolsprob &&
                 (pricetype == GCG_PRICETYPE_REDCOST || nfoundvars < pricerdata->maxvarsroundfarkas)
                 && (pricetype == GCG_PRICETYPE_FARKAS || nfoundvars < pricerdata->maxvarsroundredcost 
                    || pricerdata->onlybest); j++ )
         {
            /* create new variable, compute objective function value and add it to the master constraints and cuts it belongs to */
            SCIP_CALL( createNewMasterVar(scip, solvars[j], solvals[j], nsolvars[j], prob, 
                  pricetype == GCG_PRICETYPE_REDCOST, &added) );
            
            if ( added )
            {
               nfoundvars++;
               nfoundvarsprob++;

               if( nfoundvarsprob == 1 )
                  successfulmips++;
            }
         }
      }
      for( j = 0; j < pricerdata->npricingprobs; j++ )
         if( pricerdata->pricingprobs[j] != NULL 
            && SCIPgetStage(pricerdata->pricingprobs[j]) > SCIP_STAGE_PROBLEM)
         {
            SCIP_CALL( SCIPstartClock(scip, pricerdata->freeclock) );
            SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[j]) );
            SCIP_CALL( SCIPstopClock(scip, pricerdata->freeclock) );
         }
   }

   /* if no variables were found so far, solve the pricing MIPs to optimality and check whether
    * solutions corresponding to variables with negative reduced costs where found
    */
   if( nfoundvars == 0 )
   {
      SCIPdebugMessage("optimal pricing\n");

      bestredcostvalid = ( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL ? TRUE : FALSE );

      for( i = 0; i < pricerdata->npricingprobs && (pricetype == GCG_PRICETYPE_FARKAS || (( pricerdata->onlybest || 
                  nfoundvars < pricerdata->maxvarsroundredcost) && successfulmips < pricerdata->maxsuccessfulmipsredcost
               && successfulmips < pricerdata->successfulmipsrel * pricerdata->npricingprobsnotnull))
              && (nfoundvars == 0 || pricerdata->dualsolconv[pricerdata->permu[i]] > 0 || !pricerdata->onlyposconv)
              && (pricetype == GCG_PRICETYPE_REDCOST || nfoundvars < pricerdata->maxvarsroundfarkas); i++)
      {
         prob = pricerdata->permu[i];

         if( pricerdata->pricingprobs[prob] == NULL )
            continue;

         /* set objective limit, such that only solutions with negative reduced costs are accepted */
         //SCIP_CALL( SCIPsetObjlimit(pricerdata->pricingprobs[prob], pricerdata->dualsolconv[prob]) );

         /* set time limit */
         SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
         if( !SCIPisInfinity(scip, timelimit) )
         {
            if( timelimit - SCIPgetTotalTime(scip) > 0 )
            {
               SCIP_CALL( SCIPsetRealParam(pricerdata->pricingprobs[prob], "limits/time", 
                     timelimit - SCIPgetTotalTime(scip)) );
            }
            else
            {
               *result = SCIP_DIDNOTRUN;
               bestredcostvalid = FALSE;
               break;
            }
         }

         SCIP_CALL( solvePricingProblem(scip, pricerdata, prob, pricetype, &solvars, &solvals, 
               &nsolvars, &nsols, &status) );

         pricerdata->solvedsubmipsoptimal++;

         if( nsols > 0 )
         {
            /* compute the ojective value of the best solution */
            bestsolval = 0.0;
            for( j = 0; j < nsolvars[0]; j++ )
            {
               if( SCIPvarGetType(solvars[0][j]) != SCIP_VARTYPE_CONTINUOUS )
               {
                  assert(SCIPisEQ(scip, solvals[0][j], SCIPfloor(scip, solvals[0][j])));
                  solvals[0][j] = SCIPfloor(scip, solvals[0][j]);
               }
               bestsolval += solvals[0][j] * SCIPvarGetObj(solvars[0][j]);
            }

            if( SCIPisSumNegative(scip, bestsolval - pricerdata->dualsolconv[prob]) )
               bestredcost += GCGrelaxGetNIdenticalBlocks(origprob, prob) * 
                  (bestsolval - pricerdata->dualsolconv[prob]);
         }

         if( status != SCIP_STATUS_OPTIMAL )
            bestredcostvalid = FALSE;

         nfoundvarsprob = 0;

         for( j = 0; j < nsols && nfoundvarsprob <= pricerdata->maxsolsprob &&
                 (pricetype == GCG_PRICETYPE_REDCOST || nfoundvars < pricerdata->maxvarsroundfarkas)
                 && (pricetype == GCG_PRICETYPE_FARKAS || nfoundvars < pricerdata->maxvarsroundredcost
                    || pricerdata->onlybest); j++ )
         {
            /* create new variable, compute objective function value and add it to the master constraints and cuts it belongs to */
            SCIP_CALL( createNewMasterVar(scip, solvars[j], solvals[j], nsolvars[j], prob, 
                  pricetype == GCG_PRICETYPE_REDCOST, &added) );
            
            if ( added )
            {
               nfoundvars++;
               nfoundvarsprob++;
               
               if( nfoundvarsprob == 1 )
                  successfulmips++;
            }
         }

      }
   }

   if( pricerdata->onlybest && pricerdata->maxbestsols > 0 && pricerdata->nbestsols > 0 )
   {
      for( j = 0; j < pricerdata->nbestsols; j++ )
      {
         /* create new variable, compute objective function value and add it to the master constraints and cuts it belongs to */
         SCIP_CALL( createNewMasterVar(scip, pricerdata->bestsolvars[j], pricerdata->bestsolvals[j], 
               pricerdata->nbestsolvars[j], pricerdata->prob[j], FALSE, &added) );
         assert(added);
      }
      pricerdata->nbestsols = 0;
   }

   for( j = i; j < pricerdata->npricingprobs && bestredcostvalid; j++ )
      if( pricerdata->pricingprobs[pricerdata->permu[j]] != NULL )
         bestredcostvalid = FALSE;

   for( j = 0; j < pricerdata->npricingprobs; j++ )
      if( pricerdata->pricingprobs[j] != NULL 
         && SCIPgetStage(pricerdata->pricingprobs[j]) > SCIP_STAGE_PROBLEM)
         {
            SCIP_CALL( SCIPstartClock(scip, pricerdata->freeclock) );
            SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[j]) );
            SCIP_CALL( SCIPstopClock(scip, pricerdata->freeclock) );
         }

   if( pricetype == GCG_PRICETYPE_REDCOST && bestredcostvalid )
   {
      assert(lowerbound != NULL);
      if( pricerdata->dispinfos )
         printf("lower bound = %g, bestredcost = %g\n", SCIPgetLPObjval(scip) + bestredcost, bestredcost);

      *lowerbound = SCIPgetLPObjval(scip) + bestredcost;
   }

   SCIPdebugMessage("Pricing: found %d new vars\n", nfoundvars);

   return SCIP_OKAY;
}



/*
 * Callback methods of variable pricer
 */


/** destructor of variable pricer to free user data (called when SCIP is exiting) */

static
SCIP_DECL_PRICERFREE(pricerFreeGcg)
{ 
   SCIP_PRICERDATA* pricerdata;  

   assert(scip != NULL);
  
   /* get pricerdata */
   pricerdata = SCIPpricerGetData(pricer);

   SCIP_CALL( solversFree(scip, pricerdata) );

   SCIPfreeMemoryArray(scip, &pricerdata->solvers);

   /* free memory for pricerdata*/
   if( pricerdata != NULL )
   {
      SCIPfreeMemory(scip, &pricerdata);
   }
   
   SCIPpricerSetData(pricer, NULL);
   return SCIP_OKAY;
}



/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
static
SCIP_DECL_PRICERINITSOL(pricerInitsolGcg)
{  
   SCIP_PRICERDATA* pricerdata;
   int i;
   SCIP* origprob;
   SCIP_VAR** vars;
   int nvars;
   int v;
   SCIP_VARDATA* vardata;
   SCIP_Bool discretization;
   SCIP_CONS** masterconss;
   int nmasterconss;


   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   pricerdata->currnodenr = -1;

   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);

   /* init array containing all pricing problems */
   pricerdata->npricingprobs = GCGrelaxGetNPricingprobs(origprob);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->pricingprobs), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nvarsprob), pricerdata->npricingprobs) );
   pricerdata->npricingprobsnotnull = 0;

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if( GCGrelaxIsPricingprobRelevant(origprob, i) )
      {
         pricerdata->pricingprobs[i] = GCGrelaxGetPricingprob(origprob, i);
         pricerdata->npricingprobsnotnull++;
      }
      else
      {
         pricerdata->pricingprobs[i] = NULL;
      }
      pricerdata->nvarsprob[i] = 0;
   }
   
   /* alloc memory for arrays of reduced cost */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->dualsolconv), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->tmpconvdualsol), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->permu), pricerdata->npricingprobs) );

   /* alloc memory for solution values of variables in pricing problems */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->solvals), SCIPgetNOrigVars(origprob)) );

   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->redcostclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->farkasclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->freeclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->transformclock)) );

   pricerdata->solvedsubmipsoptimal = 0;
   pricerdata->solvedsubmipsheur = 0;
   pricerdata->calls = 0;
   pricerdata->redcostcalls = 0;
   pricerdata->farkascalls = 0;

   /* set variable type for master variables */
   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( discretization )
   {
      pricerdata->vartype = SCIP_VARTYPE_INTEGER;
   }
   else
   {
      pricerdata->vartype = SCIP_VARTYPE_CONTINUOUS;
   }

   /* for variables in the original problem that do not belong to any block, 
    * create the corresponding variable in the master problem */
   vars = SCIPgetVars(origprob);
   nvars = SCIPgetNVars(origprob);
   for( v = 0; v < nvars; v++ )
   {
      vardata = SCIPvarGetData(vars[v]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_ORIGINAL);
      if( vardata->blocknr == -1 )
      {
         SCIP_VAR* newvar;
         SCIP_VARDATA* newvardata;

         assert(vardata->data.origvardata.pricingvar == NULL);

         SCIPdebugMessage("var %s is in no block!\n", SCIPvarGetName(vars[v]));
         //printf("var %s is in no block!\n", SCIPvarGetName(vars[v]));

         /* create vardata */
         SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
         newvardata->vartype = GCG_VARTYPE_MASTER;
         newvardata->blocknr = -1;
         newvardata->data.mastervardata.norigvars = 2;

         /* save corresoponding origvar */
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, 
               &(newvardata->data.mastervardata.origvars), 2) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, 
               &(newvardata->data.mastervardata.origvals), 2) );
         newvardata->data.mastervardata.origvars[0] = vars[v];
         newvardata->data.mastervardata.origvals[0] = 1.0;
         newvardata->data.mastervardata.origvars[1] = vars[v];
         newvardata->data.mastervardata.origvals[1] = 0.0;

         /* create variable in the master problem */
         SCIP_CALL( SCIPcreateVar(scip, &newvar, SCIPvarGetName(vars[v]), 
               SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v]), SCIPvarGetObj(vars[v]), SCIPvarGetType(vars[v]), 
               TRUE, TRUE, NULL, NULL, gcgvardeltrans, newvardata) );
         SCIPaddVar(scip, newvar);

         //SCIPchgVarUbLazy(scip, newvar, SCIPvarGetUbGlobal(vars[v]));
         //SCIPchgVarLbLazy(scip, newvar, SCIPvarGetLbGlobal(vars[v]));

         SCIP_CALL( GCGpricerAddMasterVarToOrigVar(scip, vars[v], newvar, 1.0) );

         /* add variable in the master to the master constraints it belongs to */
         for( i = 0; i < vardata->data.origvardata.ncoefs; i++ )
         {
            SCIP_CONS* linkcons;
            assert(!SCIPisZero(scip, vardata->data.origvardata.coefs[i]));
            SCIP_CALL( SCIPgetTransformedCons(scip, vardata->data.origvardata.linkconss[i], 
                  &linkcons) );

            SCIP_CALL( SCIPaddCoefLinear(scip, linkcons, 
                  newvar, vardata->data.origvardata.coefs[i]) );
         }
         SCIPreleaseVar(scip, &newvar);

      }
   }

   SCIP_CALL( SCIPhashmapCreate(&(pricerdata->mapcons2idx), SCIPblkmem(scip), 10 * nmasterconss) );
   for( i = 0; i < nmasterconss; i++ )
   {
      //printf("add cons %s to hashmap: pointer %p\n", SCIPconsGetName(masterconss[i]), masterconss[i]);
      SCIP_CALL( SCIPhashmapInsert(pricerdata->mapcons2idx, masterconss[i], (void*)(size_t)i) );
      assert((int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, masterconss[i]) == i);
   }

   /* create onlybest array, if needed */
   pricerdata->maxvars = -1;
   for( i = 0; i < GCGrelaxGetNPricingprobs(origprob); i++ )
   {
      if( SCIPgetNVars(GCGrelaxGetPricingprob(origprob, i)) > pricerdata->maxvars )
         pricerdata->maxvars = SCIPgetNVars(GCGrelaxGetPricingprob(origprob, i));
   }

   if ( pricerdata->onlybest && pricerdata->maxvarsroundredcost <= MAXBEST)
   {
      pricerdata->maxbestsols = pricerdata->maxvarsroundredcost;

      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->bestsolvars, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->bestsolvals, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->nbestsolvars, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->redcost, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->prob, pricerdata->maxbestsols) );

      for( i = 0; i < pricerdata->maxbestsols; i++ )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->bestsolvars[i]), pricerdata->maxvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->bestsolvals[i]), pricerdata->maxvars) );
         pricerdata->nbestsolvars[i] = 0;
      }

      pricerdata->nbestsols = 0;
   }
   else
   {
      pricerdata->bestsolvars = NULL;
      pricerdata->bestsolvals = NULL;
      pricerdata->nbestsolvars = NULL;
      pricerdata->redcost = NULL;
      pricerdata->prob = NULL;
      pricerdata->maxbestsols = 0;
      pricerdata->nbestsols = 0;
   }

   pricerdata->npricedvars = 0;
   pricerdata->maxpricedvars = 50;
   SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->pricedvars, pricerdata->maxpricedvars) );

   SCIP_CALL( solversInitsol(scip, pricerdata) );

   return SCIP_OKAY;
}



/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolGcg)
{  
   SCIP_PRICERDATA* pricerdata;
   int i;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   SCIPhashmapFree(&(pricerdata->mapcons2idx));

   /* free bestvars array if not needed anymore */
   if ( pricerdata->onlybest && pricerdata->maxbestsols > 0 )
   {
      assert(pricerdata->bestsolvars != NULL);
      assert(pricerdata->bestsolvals != NULL);
      assert(pricerdata->nbestsolvars != NULL);
      assert(pricerdata->redcost != NULL);
      assert(pricerdata->prob != NULL);

      for( i = 0; i < pricerdata->maxbestsols; i++ )
      {
         assert(pricerdata->bestsolvars[i] != NULL);
         assert(pricerdata->bestsolvals[i] != NULL);
         
         SCIPfreeMemoryArray(scip, &(pricerdata->bestsolvars[i]));
         SCIPfreeMemoryArray(scip, &(pricerdata->bestsolvals[i]));
      }

      SCIPfreeMemoryArray(scip, &pricerdata->bestsolvars);
      SCIPfreeMemoryArray(scip, &pricerdata->bestsolvals);
      SCIPfreeMemoryArray(scip, &pricerdata->nbestsolvars);
      SCIPfreeMemoryArray(scip, &pricerdata->redcost);
      SCIPfreeMemoryArray(scip, &pricerdata->prob);


      pricerdata->bestsolvars = NULL;
      pricerdata->bestsolvals = NULL;
      pricerdata->nbestsolvars = NULL;
      pricerdata->redcost = NULL;
      pricerdata->prob = NULL;
      pricerdata->maxbestsols = 0;
      pricerdata->nbestsols = 0;
   }

   
   SCIPfreeMemoryArray(scip, &(pricerdata->pricingprobs));
   SCIPfreeMemoryArray(scip, &(pricerdata->dualsolconv));
   SCIPfreeMemoryArray(scip, &(pricerdata->tmpconvdualsol));
   SCIPfreeMemoryArray(scip, &(pricerdata->permu));
   SCIPfreeMemoryArray(scip, &(pricerdata->solvals));
   SCIPfreeMemoryArray(scip, &(pricerdata->nvarsprob));

   //SCIPfreeMemoryArray(scip, &(pricerdata->dualsolconv));

   for( i = 0; i < pricerdata->npricedvars; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &pricerdata->pricedvars[i]) );
   }
   SCIPfreeMemoryArray(scip, &pricerdata->pricedvars);

   printf("calls = %d\n", pricerdata->calls);
   printf("solved sub-MIPs heur = %d\n", pricerdata->solvedsubmipsheur);
   printf("solved sub-MIPs optimal = %d\n", pricerdata->solvedsubmipsoptimal);
   printf("farkas calls = %d, redcost calls = %d\n", pricerdata->farkascalls, pricerdata->redcostcalls);
   printf("time for farkas pricing (total): %f\n", SCIPgetClockTime(scip, pricerdata->farkasclock));
   printf("time for redcost pricing (total): %f\n", SCIPgetClockTime(scip, pricerdata->redcostclock));
   printf("time for transformation: %f\n", SCIPgetClockTime(scip, pricerdata->transformclock));
   printf("time for freeing sub-MIPs: %f\n", SCIPgetClockTime(scip, pricerdata->freeclock));

   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->redcostclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->farkasclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->freeclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->transformclock)) );

   SCIP_CALL( solversExitsol(scip, pricerdata) );
   
   return SCIP_OKAY;
}




/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostGcg)
{  
   SCIP_RETCODE retcode;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   assert(pricerdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetNTotalNodes(scip) == pricerdata->currnodenr )
   {
      pricerdata->nroundsredcost++;
   }
   else
   {
      pricerdata->currnodenr = SCIPgetNTotalNodes(scip);
      pricerdata->nroundsredcost = 0;
   }

   if( pricerdata->nroundsredcost >= pricerdata->maxroundsredcost && pricerdata->currnodenr != 1)
   {
      SCIPdebugMessage("pricing aborted at node %lld\n", pricerdata->currnodenr);
      return SCIP_OKAY;
   }

   *result = SCIP_SUCCESS;

   SCIP_CALL( SCIPstartClock(scip, pricerdata->redcostclock) );

   //printf("pricerredcost\n");
   retcode = performPricing(scip, pricer, GCG_PRICETYPE_REDCOST, result, lowerbound);

   SCIP_CALL( SCIPstopClock(scip, pricerdata->redcostclock) );

   return retcode;
}




static
SCIP_DECL_PRICERFARKAS(pricerFarkasGcg)
{  
   SCIP_RETCODE retcode;
   SCIP_PRICERDATA* pricerdata;
   
   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   assert(pricerdata != NULL);

   SCIP_CALL( SCIPstartClock(scip, pricerdata->farkasclock) );

   //printf("pricerfarkas\n");
   retcode = performPricing(scip, pricer, GCG_PRICETYPE_FARKAS, NULL, NULL);

   SCIP_CALL( SCIPstopClock(scip, pricerdata->farkasclock) );

   return retcode;
}

/* define not used callbacks as NULL */
#define pricerInitGcg NULL
#define pricerExitGcg NULL


/*
 * variable pricer specific interface methods
 */

/** creates the gcg variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origprob            /**< SCIP data structure of the original problem */
   )
{
   SCIP_PRICERDATA* pricerdata;

   SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );
   pricerdata->origprob = origprob;

   /* initialize solvers array */
   pricerdata->solvers = NULL;
   pricerdata->nsolvers = 0;


   /* include variable pricer */
   SCIP_CALL( SCIPincludePricer(scip, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerFreeGcg, pricerInitGcg, pricerExitGcg, 
         pricerInitsolGcg, pricerExitsolGcg, pricerRedcostGcg, pricerFarkasGcg,
         pricerdata) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxsuccessfulmipsredcost",
         "maximal number of pricing mips leading to new variables solved solved in one redcost pricing round",
         &pricerdata->maxsuccessfulmipsredcost, FALSE, DEFAULT_MAXSUCCESSFULMIPSREDCOST, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxvarsroundredcost",
         "maximal number of variables created in one redcost pricing round",
         &pricerdata->maxvarsroundredcost, FALSE, DEFAULT_MAXVARSROUNDREDCOST, 1, INT_MAX, 
         paramChgdOnlybestMaxvars, (SCIP_PARAMDATA*)pricerdata) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxvarsroundfarkas",
         "maximal number of variables created in one farkas pricing round",
         &pricerdata->maxvarsroundfarkas, FALSE, DEFAULT_MAXVARSROUNDFARKAS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxroundsredcost",
         "maximal number of pricing rounds per node after the root node",
         &pricerdata->maxroundsredcost, FALSE, DEFAULT_MAXROUNDSREDCOST, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxsolsprob",
         "maximal number of variables added for each block in a pricinground",
         &pricerdata->maxsolsprob, FALSE, DEFAULT_MAXSOLSPROB, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/useheurpricing",
         "should pricing be performed heuristically before solving the MIPs to optimality?",
         &pricerdata->useheurpricing, TRUE, DEFAULT_USEHEURPRICING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/onlyposconv",
         "should only pricing problems be solved with a positive dualsol of the convexity constraint, if possible?",
         &pricerdata->onlyposconv, TRUE, DEFAULT_ONLYPOSCONV, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/abortpricing",
         "should pricing be aborted due to integral objective function?",
         &pricerdata->abortpricing, TRUE, DEFAULT_ABORTPRICING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/onlybest",
         "should only the best variables (TRUE) be added in case of a maxvarsround limit or the first ones (FALSE)?",
         &pricerdata->onlybest, TRUE, DEFAULT_ONLYBEST, paramChgdOnlybestMaxvars, (SCIP_PARAMDATA*)pricerdata) );

   SCIP_CALL( SCIPaddRealParam(pricerdata->origprob, "pricing/masterpricer/successfulsubmipsrel",
         "part of the submips that are solved and lead to new variables before pricing round is aborted? (1.0 = solve all pricing MIPs)",
         &pricerdata->successfulmipsrel, FALSE, DEFAULT_SUCCESSFULMIPSREL, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/dispinfos",
         "should additional informations concerning the pricing process be displayed?",
         &pricerdata->dispinfos, FALSE, DEFAULT_DISPINFOS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/sorting",
         "which sorting method should be used to sort the pricing problems?",
         &pricerdata->sorting, FALSE, DEFAULT_SORTING, 0, 5, NULL, NULL) );


   return SCIP_OKAY;
}

/** returns the pointer to the scip instance representing the original problem */
SCIP* GCGpricerGetOrigprob(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->origprob;
}

/** returns the array of variables that were priced in during the solving process */
SCIP_VAR** GCGpricerGetPricedvars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->pricedvars;;
}

/** returns the number of variables that were priced in during the solving process */
int GCGpricerGetNPricedvars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->npricedvars;;
}


/** adds the given constraint and the given position to the hashmap of the pricer */
SCIP_RETCODE GCGpricerAddMasterconsToHashmap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint that should be added */
   int                   pos                 /**< the position of the constraint in the relaxator's masterconss array */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(pos >= 0);
   
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   SCIP_CALL( SCIPhashmapInsert(pricerdata->mapcons2idx, cons, (void*)(size_t)pos) );
   assert((int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, cons) == pos);

   SCIPdebugMessage("Added cons %s (%p) to hashmap with index %d\n", SCIPconsGetName(cons), cons, pos);
   
   return SCIP_OKAY;
}

/** includes a solver into the pricer data */
SCIP_RETCODE GCGpricerIncludeSolver(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,
   const char*           description,
   int                   priority,
   GCG_DECL_SOLVERSOLVE  ((*solversolve)),   /**<  solving method for solver */
   GCG_DECL_SOLVERSOLVEHEUR((*solveheur)),   /**<  heuristic solving method for solver */
   GCG_DECL_SOLVERFREE   ((*solverfree)),
   GCG_DECL_SOLVERINIT   ((*solverinit)),
   GCG_DECL_SOLVEREXIT   ((*solverexit)),
   GCG_DECL_SOLVERINITSOL((*solverinitsol)),
   GCG_DECL_SOLVEREXITSOL((*solverexitsol)),
   GCG_SOLVERDATA*       solverdata
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   int pos;
   
   assert(scip != NULL);
   
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   SCIP_CALL( ensureSizeSolvers(scip, pricerdata) );

   /* solvers array is sorted decreasingly wrt. the priority, find right position and shift solvers with smaller priority */
   pos = pricerdata->nsolvers;
   while( pos >= 1 && pricerdata->solvers[pos-1]->priority < priority )
   {
      pricerdata->solvers[pos] = pricerdata->solvers[pos-1];
      pos--;
   }
   SCIP_CALL( SCIPallocMemory(scip, &(pricerdata->solvers[pos])) );

   SCIP_ALLOC( BMSduplicateMemoryArray(&pricerdata->solvers[pos]->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&pricerdata->solvers[pos]->description, description, strlen(description)+1) );

   pricerdata->solvers[pos]->priority = priority;
   pricerdata->solvers[pos]->solversolve = solversolve;
   pricerdata->solvers[pos]->solversolveheur = solveheur;
   pricerdata->solvers[pos]->solverfree = solverfree;
   pricerdata->solvers[pos]->solverinit = solverinit;
   pricerdata->solvers[pos]->solverexit = solverexit;
   pricerdata->solvers[pos]->solverinitsol = solverinitsol;
   pricerdata->solvers[pos]->solverexitsol = solverexitsol;
   pricerdata->solvers[pos]->solverdata = solverdata;


   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->solvers[pos]->optfarkasclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->solvers[pos]->optredcostclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->solvers[pos]->heurfarkasclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->solvers[pos]->heurredcostclock)) );

   pricerdata->solvers[pos]->optfarkascalls = 0;
   pricerdata->solvers[pos]->optredcostcalls = 0;
   pricerdata->solvers[pos]->heurfarkascalls = 0;
   pricerdata->solvers[pos]->heurredcostcalls = 0;

   pricerdata->nsolvers++;

   return SCIP_OKAY;
}


GCG_SOLVERDATA* GCGpricerGetSolverdata(
   SCIP*                 scip,
   GCG_SOLVER*           solver
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   
   assert(scip != NULL);
   assert(solver != NULL);
   
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));   
   assert(pricerdata->nsolvers > 0);

   return solver->solverdata;
}

void GCGpricerSetSolverdata(
   SCIP*                 scip,
   GCG_SOLVER*          solver,
   GCG_SOLVERDATA*       solverdata
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   
   assert(scip != NULL);
   assert(solver != NULL);
   
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));   
   assert(pricerdata->nsolvers > 0);

   solver->solverdata = solverdata;
}




void GCGpricerPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   int i;

   assert(scip != NULL);
   
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /**@todo add constraint statistics: how many constraints (instead of cuts) have been added? */
   SCIPmessageFPrintInfo(file, "Pricing Solver     : #HeurFarkas  #OptFarkas  #HeurRedcost #OptRedcost Time: HeurFarkas  OptFarkas  HeurRedcost OptRedcost\n");

   for( i = 0; i < pricerdata->nsolvers; ++i )
   {
      GCG_SOLVER* solver;
      solver = pricerdata->solvers[i];
      assert(solver != NULL);

      SCIPmessageFPrintInfo(file, "  %-17.17s:", solver->name);
      SCIPmessageFPrintInfo(file, " %11d %11d   %11d %11d       %10.2f %10.2f   %10.2f %10.2f \n",
         solver->heurfarkascalls, solver->optfarkascalls, 
         solver->heurredcostcalls, solver->optredcostcalls, 
         SCIPgetClockTime(scip, solver->heurfarkasclock),
         SCIPgetClockTime(scip, solver->optfarkasclock),
         SCIPgetClockTime(scip, solver->heurredcostclock),
         SCIPgetClockTime(scip, solver->optredcostclock));
   }
}
