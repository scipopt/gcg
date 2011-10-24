/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* #define SCIP_DEBUG */

/**@file   pricer_gcg.c
 * @ingroup PRICERS
 * @brief  pricer for generic column generation
 * @author Gerald Gamrath
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"

#include "pricer_gcg.h"
#include "sepa_master.h"

#include "relax_gcg.h"
#include "struct_solver.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"

#define PRICER_NAME            "gcg"
#define PRICER_DESC            "pricer for gcg"
#define PRICER_PRIORITY        5000000
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

#define DEFAULT_MAXVARSROUNDFARKAS 10
#define DEFAULT_MAXVARSROUNDREDCOSTROOT 100
#define DEFAULT_MAXVARSROUNDREDCOST 100
#define DEFAULT_MAXSUCCESSFULMIPSREDCOST INT_MAX
#define DEFAULT_MAXROUNDSREDCOST INT_MAX
#define DEFAULT_MAXSOLSPROB INT_MAX
#define DEFAULT_USEHEURPRICING FALSE
#define DEFAULT_ONLYPOSCONV FALSE
#define DEFAULT_ABORTPRICINGINT TRUE
#define DEFAULT_ABORTPRICINGGAP 0.00
#define DEFAULT_USEINTERBOUNDS TRUE
#define DEFAULT_ONLYBEST FALSE
#define DEFAULT_SUCCESSFULMIPSREL 1.0
#define DEFAULT_MIPSRELREDCOSTROOT 1.0
#define DEFAULT_MIPSRELREDCOST 1.0
#define DEFAULT_MIPSRELFARKAS 1.0
#define DEFAULT_DISPINFOS FALSE
#define DEFAULT_SORTING 2

#define MAXBEST 1000


#define GCGpricerPrintInfo(scip,pricerdata, ...) do { \
   if( pricerdata->dispinfos ) { \
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,__VA_ARGS__);\
   } else {\
      SCIPdebugMessage(__VA_ARGS__); \
   }\
   }while(FALSE)

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
   int* npointsprob;               /* number of variables representing points created by the pricing probs */
   int* nraysprob;                 /* number of variables representing rays created by the pricing probs */
   SCIP_Longint currnodenr;
   SCIP_HASHMAP* mapcons2idx;
   SCIP_Real* score;
   int* permu;
   int npricingprobsnotnull;

   SCIP_Real** bestsolvals;
   SCIP_VAR*** bestsolvars;
   int* nbestsolvars;
   SCIP_Bool* bestsolisray;
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
   int maxvarsroundredcostroot;
   int maxsuccessfulmipsredcost;
   int maxroundsredcost;
   int maxsolsprob;
   int nroundsredcost;
   int sorting;
   SCIP_Bool useheurpricing;
   SCIP_Bool onlyposconv;
   SCIP_Bool abortpricingint;
   SCIP_Bool useinterbounds;
   SCIP_Bool onlybest;
   SCIP_Bool dispinfos;
   SCIP_Real successfulmipsrel;
   SCIP_Real mipsrelredcost;
   SCIP_Real mipsrelredcostroot;
   SCIP_Real mipsrelfarkas;
   SCIP_Real abortpricinggap;
};


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
   if( !pricerdata->onlybest && pricerdata->maxbestsols > 0 )
   {
      assert(pricerdata->bestsolvars != NULL);
      assert(pricerdata->bestsolvals != NULL);
      assert(pricerdata->nbestsolvars != NULL);
      assert(pricerdata->bestsolisray != NULL);
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
      SCIPfreeMemoryArray(scip, &pricerdata->bestsolisray);
      SCIPfreeMemoryArray(scip, &pricerdata->redcost);
      SCIPfreeMemoryArray(scip, &pricerdata->prob);

      pricerdata->bestsolvars = NULL;
      pricerdata->bestsolvals = NULL;
      pricerdata->nbestsolvars = NULL;
      pricerdata->bestsolisray = NULL;
      pricerdata->maxbestsols = 0;
      pricerdata->nbestsols = 0;
   }

   /* create array */
   if( pricerdata->onlybest && pricerdata->maxbestsols == 0 && pricerdata->maxvarsroundredcost <= MAXBEST)
   {
      assert(pricerdata->bestsolvars == NULL);
      assert(pricerdata->bestsolvals == NULL);
      assert(pricerdata->nbestsolvars == NULL);
      assert(pricerdata->bestsolisray == NULL);
      assert(pricerdata->redcost == NULL);
      assert(pricerdata->prob == NULL);

      pricerdata->maxbestsols = pricerdata->maxvarsroundredcost;

      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->bestsolvars, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->bestsolvals, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->nbestsolvars, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->bestsolisray, pricerdata->maxbestsols) );
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
   if( pricerdata->onlybest && pricerdata->maxbestsols != 0 )
   {
      assert(pricerdata->bestsolvars != NULL);
      assert(pricerdata->bestsolvals != NULL);
      assert(pricerdata->nbestsolvars != NULL);
      assert(pricerdata->bestsolisray != NULL);
      assert(pricerdata->redcost != NULL);
      assert(pricerdata->prob != NULL);

      SCIP_CALL( SCIPreallocMemoryArray(scip, &pricerdata->bestsolvars, pricerdata->maxvarsroundredcost) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &pricerdata->bestsolvals, pricerdata->maxvarsroundredcost) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &pricerdata->nbestsolvars, pricerdata->maxvarsroundredcost) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &pricerdata->bestsolisray, pricerdata->maxvarsroundredcost) );
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

    SCIPdebugMessage("paramchanged\n");

   return SCIP_OKAY;
}



/*
 * Local methods
 */

/** returns TRUE or FALSE, depending whether we are in the root node or not */
static SCIP_Bool isRootNode(
   SCIP* scip                           /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   return (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip));
}


/* ensures size of pricedvars array */
static
SCIP_RETCODE ensureSizePricedvars(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,        /**< Pricerdata data structure */
   int                   size
   )
{
   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert(pricerdata->pricedvars != NULL);

   if( pricerdata->maxpricedvars < size )
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
   SCIP*                 scip,              /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata         /**< Pricerdata data structure  */
   )
{
   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));

   if( pricerdata->nsolvers == 0 )
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
   SCIP*                 scip,              /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata         /**< Pricerdata data structure  */
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
   SCIP*                 scip,              /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata         /**< Pricerdata data structure  */
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
   SCIP*                 scip,              /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata         /**< Pricerdata data structure  */
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
   SCIP*                 scip,              /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata         /**< Pricerdata data structure  */
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
   SCIP*                 scip,              /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata         /**< Pricerdata data structure  */
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
   SCIP_Bool**           solisray,
   int*                  nsols,
   SCIP_STATUS*          status
   )
{
   int i;
   SCIP_Real timelimit;

   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert(pricerdata->pricingprobs[prob] != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   *status = SCIP_STATUS_UNKNOWN;

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {

      /*get time limit */
      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      if( !SCIPisInfinity(scip, timelimit) && timelimit - SCIPgetSolvingTime(scip) < 0 )
      {
         *nsols = 0;
         *status = SCIP_STATUS_TIMELIMIT;
      }

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
               solvars, solvals, nsolvars, solisray, nsols, status) );

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

         if( *status == SCIP_STATUS_OPTIMAL || *status == SCIP_STATUS_UNBOUNDED )
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
   SCIP_Bool**           solisray,
   int*                  nsols,
   SCIP_STATUS*          status
   )
{
   int i;
   SCIP_Real timelimit;

   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert(pricerdata->pricingprobs[prob] != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   *status = SCIP_STATUS_UNKNOWN;

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      /*get time limit */
      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      if( !SCIPisInfinity(scip, timelimit) && timelimit - SCIPgetSolvingTime(scip) < 1 )
      {
         *nsols = 0;
         *status = SCIP_STATUS_TIMELIMIT;
      }

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
               solvars, solvals, nsolvars, solisray, nsols, status) );


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

         if( *status == SCIP_STATUS_OPTIMAL || *status == SCIP_STATUS_UNBOUNDED )
            break;
      }

   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE setPricingObjs(
   SCIP*                 scip,              /**< SCIP data structure            */
   GCG_PRICETYPE         pricetype          /**< Farkas or Reduced cost pricing */
   )
{
   SCIP* origprob;
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
            SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], 0.0) );
         }
         else
         {
            SCIP_VAR* origvar;
            origvar = GCGpricingVarGetOrigvars(probvars[j])[0];

            assert(GCGvarGetBlock(probvars[j]) == i);

            if( GCGvarIsLinking(origvar) )
            {
               SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], 0.0) );
            }
            else
            {
               assert( GCGvarGetBlock(origvar) == i);
               SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], SCIPvarGetObj(origvar)) );
            }
         }
      }
   }

   /* compute reduced cost for linking variable constraints and update objectives in the pricing problems */
   for( i = 0; i < pricerdata->npricingprobs; i++)
   {
      if( pricerdata->pricingprobs[i] == NULL )
         continue;
      probvars = SCIPgetVars(pricerdata->pricingprobs[i]);
      nprobvars = SCIPgetNVars(pricerdata->pricingprobs[i]);

      for( j = 0; j < nprobvars; j++ )
      {
         SCIP_VAR* origvar;
         SCIP_VAR** pricingvars;
         SCIP_CONS** linkconss;

         origvar = GCGpricingVarGetOrigvars(probvars[j])[0];

         assert(GCGvarIsPricing(probvars[j]));
         assert(GCGvarGetBlock(probvars[j]) == i);

         if( !GCGvarIsLinking(origvar) )
            continue;

         pricingvars = GCGlinkingVarGetPricingVars(origvar);
         linkconss = GCGlinkingVarGetLinkingConss(origvar);
         assert(pricingvars[i] == probvars[j]);
         assert(linkconss[i] != NULL);

         /* redcost pricing */
         if( pricetype == GCG_PRICETYPE_REDCOST )
            dualsol = SCIPgetDualsolLinear(scip, linkconss[i]);
         /* farkas pricing */
         else
         {
            assert(pricetype == GCG_PRICETYPE_FARKAS);
            dualsol = SCIPgetDualfarkasLinear(scip, linkconss[i]);
         }

         /* add dual solution value to the pricing variable:
          * lambda variables get coef -1 in linking constraints --> add dualsol */
         SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[i], probvars[j], dualsol) );
      }
   }

   /* compute reduced cost and update objectives in the pricing problems */
   for( i = 0; i < nmasterconss; i++ )
   {
      /* redcost pricing */
      if( pricetype == GCG_PRICETYPE_REDCOST )
         dualsol = SCIPgetDualsolLinear(scip, masterconss[i]);
      /* farkas pricing */
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
            int blocknr;
            blocknr = GCGvarGetBlock(consvars[j]);
            assert(GCGvarIsOriginal(consvars[j]));
            /* nothing to be done if variable belongs to redundant block or
             * variable was directly transferred to the master
             * or variable is linking variable (which means, the directly transferred copy is part of the master cons) */
            if( blocknr >= 0 && pricerdata->pricingprobs[blocknr] != NULL )
            {
               assert(GCGoriginalVarGetPricingVar(consvars[j]) != NULL);
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[blocknr],
                     GCGoriginalVarGetPricingVar(consvars[j]), -1.0 * dualsol * consvals[j]) );
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
            int blocknr;
            blocknr = GCGvarGetBlock(consvars[j]);
            assert(GCGvarIsOriginal(consvars[j]));
            /* nothing to be done if variable belongs to redundant block or
             * variable was directly transferred to the master
             * or variable is linking variable (which means, the directly transferred copy is part of the master cut) */
            if( blocknr >= 0 && pricerdata->pricingprobs[blocknr] != NULL )
            {
               assert(GCGoriginalVarGetPricingVar(consvars[j]) != NULL);
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[blocknr],
                     GCGoriginalVarGetPricingVar(consvars[j]), -1.0 * dualsol * consvals[j]) );
            }
         }
         SCIPfreeBufferArray(scip, &consvars);
      }
   }

   /* get dual solutions / farkas values of the convexity constraints */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {

      assert( GCGrelaxIsPricingprobRelevant(origprob, i) == (GCGrelaxGetConvCons(origprob, i) != NULL) );
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

/** add master variable to all constraints */
static
SCIP_RETCODE addVariableToMasterconstraints(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,        /**< Pricerdata data structure */
   SCIP_VAR*             newvar,            /**< The new variable to add */
   int                   prob,              /**< number of the pricing problem the solution belongs to */
   SCIP_VAR**            solvars,           /**< array of variables with non-zero value in the solution of the pricing problem */
   SCIP_Real*            solvals,           /**< array of values in the solution of the pricing problem for variables in array solvars*/
   int                   nsolvars           /**< number of variables in array solvars */
   )
{
   int i;
   int c;
   int idx;

   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_Real* mastercoefs;
   SCIP_CONS* linkcons;

   nmasterconss = GCGrelaxGetNMasterConss(pricerdata->origprob);
   masterconss = GCGrelaxGetMasterConss(pricerdata->origprob);

   SCIP_CALL( SCIPallocBufferArray(scip, &mastercoefs, nmasterconss) );
   BMSclearMemoryArray(mastercoefs, nmasterconss);

   /* compute coef of the variable in the master constraints */
   for( i = 0; i < nsolvars; i++ )
   {
      if( !SCIPisZero(scip, solvals[i]) )
      {
         SCIP_CONS** linkconss;
         SCIP_VAR** origvars;
         SCIP_Real* coefs;
         int ncoefs;

         assert(GCGvarIsPricing(solvars[i]));
         origvars = GCGpricingVarGetOrigvars(solvars[i]);
         assert(GCGvarIsOriginal(origvars[0]));

         coefs = GCGoriginalVarGetCoefs(origvars[0]);
         ncoefs = GCGoriginalVarGetNCoefs(origvars[0]);
         assert(!SCIPisInfinity(scip, solvals[i]));

         /* original variable is a linking variable, just add it to the linkcons */
         if( GCGvarIsLinking(origvars[0]) )
         {
            SCIP_VAR** pricingvars;
            linkconss = GCGlinkingVarGetLinkingConss(origvars[0]);
            pricingvars = GCGlinkingVarGetPricingVars(origvars[0]);
            assert(pricingvars[prob] == solvars[i]);
            assert(linkconss[prob] != NULL);
            SCIP_CALL( SCIPaddCoefLinear(scip, linkconss[prob], newvar, -solvals[i]) );
            continue;
         }

         /* for each coef, add coef * solval to the coef of the new variable for the corresponding constraint */
         for( c = 0; c < ncoefs; c++ )
         {
            linkconss = GCGoriginalVarGetLinkingCons(origvars[0]);
            assert(!SCIPisZero(scip, coefs[c]));
            SCIP_CALL( SCIPgetTransformedCons(scip, linkconss[c], &linkcons) );

            idx = (int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, linkcons);
            assert(0 <= idx && idx < nmasterconss);
            assert(masterconss[idx] == linkcons);
            mastercoefs[idx] += coefs[c] * solvals[i];
         }

      }
   }

   /* add the variable to the master constraints */
   for( i = 0; i < nmasterconss; i++ )
   {
      if( !SCIPisZero(scip, mastercoefs[i]) )
      {
         assert(!SCIPisInfinity(scip, mastercoefs[i]));
         SCIP_CALL( SCIPaddCoefLinear(scip, masterconss[i], newvar, mastercoefs[i]) );
      }
   }

   SCIPfreeBufferArray(scip, &mastercoefs);
   return SCIP_OKAY;
}



/** add variable with computed coefficients to the master cuts */
static
SCIP_RETCODE addVariableToMastercuts(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_VAR*             newvar,            /**< The new variable to add */
   int                   prob,              /**< number of the pricing problem the solution belongs to */
   SCIP_VAR**            solvars,           /**< array of variables with non-zero value in the solution of the pricing problem */
   SCIP_Real*            solvals,           /**< array of values in the solution of the pricing problem for variables in array solvars*/
   int                   nsolvars           /**< number of variables in array solvars */
   )
{
   SCIP_ROW** mastercuts;
   int nmastercuts;
   SCIP_ROW** origcuts;
   int norigcuts;

   SCIP_COL** cols;
   SCIP_Real conscoef;
   SCIP_VAR* var;
   SCIP_Real* consvals;

   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(newvar != NULL);
   assert(solvars != NULL);
   assert(solvals != NULL);

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
         int blocknr;
         var = SCIPcolGetVar(cols[j]);
         blocknr = GCGvarGetBlock(var);
         assert(GCGvarIsOriginal(var));

         /* if the belongs to the same block and is no linking variable, update the coef */
         if( blocknr == prob )
         {
            for( k = 0; k < nsolvars; k++ )
            {
               if( solvars[k] == GCGoriginalVarGetPricingVar(var) )
               {
                  conscoef += ( consvals[j] * solvals[k] );
                  break;
               }
            }
         }

      }

      if( !SCIPisZero(scip, conscoef) )
         SCIP_CALL( SCIPaddVarToRow(scip , mastercuts[i], newvar, conscoef) );
   }
   return SCIP_OKAY;
}

/** creates a new master variable corresponding to the given solution and problem */
static
SCIP_RETCODE createNewMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
   SCIP_Real*            solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
   int                   nsolvars,           /**< number of variables in array solvars */
   SCIP_Bool             solisray,           /**< is the solution a ray? */
   int                   prob,               /**< number of the pricing problem the solution belongs to */
   SCIP_Bool             checkonlybest,      /**< should the method just store the solution if it is among the best ones?
                                              *   call this method later with checkonlybest = FALSE to add the best ones */
   SCIP_Bool             force,              /**< should the given variable be added also if it has non-negative reduced cost? */
   SCIP_Bool*            added,              /**< pointer to store whether the variable was successfully added */
   SCIP_VAR**            addedvar            /**< pointer to store the created variable */
   )
{
   SCIP* origprob;
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   char varname[SCIP_MAXSTRLEN];

   SCIP_Real objcoeff;
   SCIP_VAR* newvar;

   SCIP_Real objvalue;
   SCIP_Real redcost;

   int i;

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

   if( addedvar != NULL )
      *addedvar = NULL;

   objvalue = 0.0;
   redcost = 0.0;

   if( !force )
   {
      /* compute the objective function value of the solution */
      for( i = 0; i < nsolvars; i++ )
         objvalue += solvals[i] * SCIPvarGetObj(solvars[i]);

      /* compute reduced cost of variable (i.e. subtract dual solution of convexity constraint, if solution corresponds to a point) */
      redcost = ( solisray ? objvalue : objvalue - pricerdata->dualsolconv[prob]);

      if( !SCIPisSumNegative(scip, redcost) )
      {
         *added = FALSE;

         return SCIP_OKAY;
      }
   }

   *added = TRUE;

   SCIPdebugMessage("found var with redcost %g (objvalue = %g, dualsol =%g)\n", redcost, objvalue, pricerdata->dualsolconv[prob]);
   //printf("found var with redcost %f (objvalue = %g, dualsol =%g), checkonlybest = %d\n",
   //   redcost, objvalue, pricerdata->dualsolconv[prob], checkonlybest);

   if( checkonlybest && pricerdata->onlybest && pricerdata->maxbestsols > 0 )
   {
      int pos;

      for( pos = pricerdata->nbestsols - 1; pos >= 0 && pricerdata->redcost[pos] > redcost; pos-- )
      {
         if( pos < pricerdata->maxbestsols - 1 )
         {
            pricerdata->prob[pos+1] = pricerdata->prob[pos];
            pricerdata->redcost[pos+1] = pricerdata->redcost[pos];
            pricerdata->nbestsolvars[pos+1] = pricerdata->nbestsolvars[pos];
            pricerdata->bestsolisray[pos+1] = pricerdata->bestsolisray[pos];
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
         pricerdata->redcost[pos] = redcost;
         pricerdata->nbestsolvars[pos] = nsolvars;
         pricerdata->bestsolisray[pos] = solisray;

         for( i = 0; i < pricerdata->nbestsolvars[pos]; i++ )
         {
            pricerdata->bestsolvars[pos][i] = solvars[i];
            pricerdata->bestsolvals[pos][i] = solvals[i];
         }
         pricerdata->nbestsols++;
      }

      return SCIP_OKAY;
   }

   /* compute objective coefficient of the variable */
   objcoeff = 0;
   for( i = 0; i < nsolvars; i++ )
   {
      if( !SCIPisZero(scip, solvals[i]) )
      {
         SCIP_VAR* origvar;

         assert(GCGvarIsPricing(solvars[i]));
         origvar = GCGpricingVarGetOrigvars(solvars[i])[0];

         /* original variable is linking variable --> directly transferred master variable got the full obj,
          * priced-in variables get no objective value for this origvar */
         if( GCGvarIsLinking(origvar) )
            continue;

         /* add quota of original variable's objcoef to the master variable's coef */
         objcoeff += solvals[i] * SCIPvarGetObj(origvar);
      }
   }

   if( SCIPisInfinity(scip, objcoeff) )
   {
      SCIPwarningMessage("variable with infinite objective value found in pricing, change objective to SCIPinfinity()/2\n");
      objcoeff = SCIPinfinity(scip) / 2;
   }

   if( solisray )
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "r_%d_%d", prob, pricerdata->nraysprob[prob]);
      pricerdata->nraysprob[prob]++;
   }
   else
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "p_%d_%d", prob, pricerdata->npointsprob[prob]);
      pricerdata->npointsprob[prob]++;
   }

   SCIP_CALL( GCGcreateMasterVar(scip, pricerdata->pricingprobs[prob], &newvar, varname, objcoeff,
         pricerdata->vartype, solisray, prob, nsolvars, solvals, solvars));

   SCIPdebugMessage("found var %s with redcost %f!\n", SCIPvarGetName(newvar), redcost);

   /* add variable */
   if( !force )
   {
      SCIP_CALL( SCIPaddPricedVar(scip, newvar, pricerdata->dualsolconv[prob] - objvalue) );
   }
   else
   {
      SCIP_CALL( SCIPaddVar(scip, newvar) );
   }

   SCIP_CALL( SCIPcaptureVar(scip, newvar) );
   SCIP_CALL( ensureSizePricedvars(scip, pricerdata, pricerdata->npricedvars + 1) );
   pricerdata->pricedvars[pricerdata->npricedvars] = newvar;
   pricerdata->npricedvars++;

   SCIP_CALL(addVariableToMasterconstraints(scip, pricerdata, newvar, prob, solvars, solvals, nsolvars));

   SCIP_CALL(addVariableToMastercuts(scip, newvar, prob, solvars, solvals, nsolvars));

   /* add variable to convexity constraint */
   if( !solisray )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, GCGrelaxGetConvCons(origprob, prob), newvar, 1.0) );
   }

   if( addedvar != NULL )
      *addedvar = newvar;

   SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

   return SCIP_OKAY;
}

/** compute the ojective value of the best solution */
static
SCIP_Real computeSolObjValue(
   SCIP*       scip,          /**< SCIP data structure                  */
   int         nsolvars,      /**< number of variables in the solution  */
   double*     solvals,       /**< Array of solution values             */
   SCIP_VAR**  solvars        /**< Array of solution variables          */
)
{
   int j;
   SCIP_Real bestsolval;
   assert(scip != NULL);
   assert(nsolvars >= 0);
   assert(solvals != NULL);
   assert(solvars != NULL);

   bestsolval = 0.0;

   for ( j = 0; j < nsolvars; j++ )
   {
      /** @todo: round solution values??? */
      assert(solvars[j] != NULL);
      bestsolval += solvals[j] * SCIPvarGetObj(solvars[j]);
   }
   return bestsolval;

}

/**
 * check whether pricing can be aborted:
 * if objective value is always integral and the current node's current
 * lowerbound rounded up equals the current lp objective value rounded
 * up we don't need to continue pricing since the best possible feasible
 * solution must have at least this value
 */
static
SCIP_Bool canPricingBeAborted(
   SCIP*             scip,          /**< SCIP data structure */
   SCIP_PRICERDATA*  pricerdata     /**< */
   )
{
   SCIP_Bool canabort;
   assert(scip != NULL);
   assert(pricerdata != NULL);
   canabort = FALSE;
   if ( pricerdata->abortpricingint && SCIPisObjIntegral(scip)
      && SCIPisEQ(scip, SCIPceil(scip, SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip))), SCIPceil(scip, SCIPgetLPObjval(scip))) /* && SCIPgetNNodes(scip) > 1 ??????*/)
   {
      GCGpricerPrintInfo(scip, pricerdata,
            "pricing aborted due to integral objective: node LB = %g, LP obj = %g\n",
            SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip)), SCIPgetLPObjval(scip));

      canabort = TRUE;
   }
   if ( pricerdata->abortpricinggap > 0 )
   {
      SCIP_Real gap;
      gap = ABS((SCIPgetLPObjval(scip) - SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip)))/SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip)));

      if ( gap < pricerdata->abortpricinggap )
      {
         GCGpricerPrintInfo(scip, pricerdata,
               "pricing aborted due to small gap: node LB = %g, LP obj = %g, gap = %g\n",
               SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip)), SCIPgetLPObjval(scip), gap);
         canabort = TRUE;
      }
   }

   return canabort;
}

static
void sortPricingProblemsByScore(SCIP_PRICERDATA *pricerdata)
{
   int i;
   assert(pricerdata != NULL);
    /* TODO: sort w.r.t. other measures? Don't sort in Farkas pricing? Randomized? */
   for (i = 0; i < pricerdata->npricingprobs; i++)
   {
      pricerdata->permu[i] = i;
      if (pricerdata->sorting == 1)
         pricerdata->score[i] = pricerdata->dualsolconv[i];
      else if (pricerdata->sorting == 2)
         pricerdata->score[i] = -(0.2 * pricerdata->npointsprob[i] + pricerdata->nraysprob[i]);
   }

   if (pricerdata->sorting > 0)
      SCIPsortDownRealInt(pricerdata->score, pricerdata->permu, pricerdata->npricingprobs);
}

/** performs the pricing routine, gets the type of pricing that should be done: farkas or redcost pricing */
static
SCIP_RETCODE performPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer,             /**< the pricer */
   GCG_PRICETYPE         pricetype,          /**< type of the pricing */
   SCIP_RESULT*          result,             /**< result pointer */
   SCIP_Real*            lowerbound          /**< lowerbound pointer */
   )
{
   SCIP_PRICERDATA* pricerdata;            /* the data of the pricer */
   SCIP* origprob;

   int i;
   int j;
   int prob;
   int solvedmips;

   int nfoundvars;
   int nfoundvarsprob;
   int successfulmips;

   SCIP_Real timelimit;

   SCIP_Real bestsolval;
   SCIP_Real bestredcost;
   SCIP_Bool bestredcostvalid;
   SCIP_Bool added;
   SCIP_Bool root;
   SCIP_Bool duringheurpricing;

   SCIP_STATUS status;

   SCIP_VAR*** solvars;
   SCIP_Real** solvals;
   int* nsolvars;
   SCIP_Bool* solisray;
   int nsols;

   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   assert(result != NULL || pricetype == GCG_PRICETYPE_FARKAS);
   assert(lowerbound != NULL || pricetype == GCG_PRICETYPE_FARKAS);

   if(lowerbound != NULL)
      *lowerbound = -SCIPinfinity(scip);

   duringheurpricing = pricerdata->useheurpricing;
   root = isRootNode(scip);

   GCGpricerPrintInfo(scip, pricerdata, "nvars = %d, current LP objval = %g, time = %f, node = %lld\n",
         SCIPgetNVars(scip), SCIPgetLPObjval(scip), SCIPgetSolvingTime(scip), SCIPgetNNodes(scip));

   if(pricetype == GCG_PRICETYPE_REDCOST)
   {
      if( canPricingBeAborted(scip, pricerdata) )
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }

      pricerdata->redcostcalls++;
      *result = SCIP_SUCCESS;
   }
   else if( pricetype == GCG_PRICETYPE_FARKAS )
   {
      pricerdata->farkascalls++;
   }

   pricerdata->calls++;
   nfoundvars = 0;
   successfulmips = 0;
   solvedmips = 0;

   /* set objectives of the variables in the pricing sub-MIPs */
   SCIP_CALL( setPricingObjs(scip, pricetype) );

   sortPricingProblemsByScore(pricerdata);

   bestredcost = 0.0;
   bestredcostvalid = FALSE;

   i = 0; /* to make lint happy */
   if( pricerdata->useheurpricing )
   {
      SCIPdebugMessage("heuristical pricing\n");

      /* solve the pricing MIPs heuristically and check whether solutions
       * corresponding to variables with negative reduced costs where found
       */
      for( i = 0; i < pricerdata->npricingprobs && (pricetype == GCG_PRICETYPE_FARKAS || ((pricerdata->onlybest ||
                  nfoundvars < pricerdata->maxvarsroundredcost) && successfulmips < pricerdata->maxsuccessfulmipsredcost
               && successfulmips < pricerdata->successfulmipsrel * pricerdata->npricingprobsnotnull
               && (nfoundvars == 0 || solvedmips < pricerdata->mipsrelredcost * pricerdata->npricingprobsnotnull )))
              && (nfoundvars == 0 || pricerdata->dualsolconv[pricerdata->permu[i]] > 0 || !pricerdata->onlyposconv )
              && (pricetype == GCG_PRICETYPE_REDCOST || (nfoundvars < pricerdata->maxvarsroundfarkas
                    && (nfoundvars == 0 || solvedmips < pricerdata->mipsrelfarkas * pricerdata->npricingprobsnotnull))); i++)
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
            if( timelimit - SCIPgetSolvingTime(scip) > 0 )
            {
               SCIP_CALL( SCIPsetRealParam(pricerdata->pricingprobs[prob], "limits/time",
                     timelimit - SCIPgetSolvingTime(scip)) );
               SCIPdebugMessage("Tilim for pricing %d is %f\n", prob, timelimit- SCIPgetSolvingTime(scip) + 5);
            }
            else
            {
               SCIPdebugMessage("Tilim for pricing %d is < 0\n", prob);
               if( pricetype == GCG_PRICETYPE_REDCOST )
                  *result = SCIP_DIDNOTRUN;

               return SCIP_OKAY;
            }
         }

         pricerdata->solvedsubmipsheur++;
         solvedmips++;

         SCIP_CALL( solvePricingProblemHeur(scip, pricerdata, prob, pricetype, &solvars, &solvals,
               &nsolvars, &solisray, &nsols, &status) );

         //printf("Pricingprob %d has found %d sols!\n", prob, nsols);

         nfoundvarsprob = 0;

         for( j = 0; j < nsols && nfoundvarsprob <= pricerdata->maxsolsprob &&
                 (pricetype == GCG_PRICETYPE_REDCOST || nfoundvars < pricerdata->maxvarsroundfarkas)
                 && (pricetype == GCG_PRICETYPE_FARKAS || nfoundvars < pricerdata->maxvarsroundredcost
                    || pricerdata->onlybest); j++ )
         {
            /* create new variable, compute objective function value and add it to the master constraints and cuts it belongs to */
            SCIP_CALL( createNewMasterVar(scip, solvars[j], solvals[j], nsolvars[j], solisray[j], prob,
                  pricetype == GCG_PRICETYPE_REDCOST, FALSE, &added, NULL) );

            if( added )
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
              (((nfoundvars < pricerdata->maxvarsroundredcostroot) || !root ) && ((nfoundvars < pricerdata->maxvarsroundredcost) || root)))
               && successfulmips < pricerdata->maxsuccessfulmipsredcost
               && successfulmips < pricerdata->successfulmipsrel * pricerdata->npricingprobsnotnull
               && (nfoundvars == 0 || ( (root || solvedmips < pricerdata->mipsrelredcost * pricerdata->npricingprobsnotnull)
                   && (!root || solvedmips < pricerdata->mipsrelredcostroot * pricerdata->npricingprobsnotnull) ) )))
              && (nfoundvars == 0 || pricerdata->dualsolconv[pricerdata->permu[i]] > 0 || !pricerdata->onlyposconv)
              && (pricetype == GCG_PRICETYPE_REDCOST || (nfoundvars < pricerdata->maxvarsroundfarkas
                    && (nfoundvars == 0 || solvedmips < pricerdata->mipsrelfarkas * pricerdata->npricingprobsnotnull))); i++)
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
            if( timelimit - SCIPgetSolvingTime(scip) > 0 )
            {
               SCIP_CALL( SCIPsetRealParam(pricerdata->pricingprobs[prob], "limits/time",
                     timelimit - SCIPgetSolvingTime(scip)) );
               SCIPdebugMessage("Tilim for pricing %d is %f\n", prob, timelimit- SCIPgetSolvingTime(scip) + 5);
            }
            else
            {
               SCIPdebugMessage("Tilim for pricing %d is < 0\n", prob);
               if( pricetype == GCG_PRICETYPE_REDCOST )
                  *result = SCIP_DIDNOTRUN;

               bestredcostvalid = FALSE;
               break;
            }
         }

         SCIP_CALL( solvePricingProblem(scip, pricerdata, prob, pricetype, &solvars, &solvals,
               &nsolvars, &solisray, &nsols, &status) );

         pricerdata->solvedsubmipsoptimal++;
         solvedmips++;

         if( nsols > 0 )
         {
            /* compute the ojective value of the best solution */
            bestsolval = 0.0;

            for( j = 0; j < nsolvars[0]; j++ )
            {
               /* TODO: round solution values??? */
               bestsolval += solvals[0][j] * SCIPvarGetObj(solvars[0][j]);
            }

            /* TODO: ensure that the first solution is really the best one and that its objective value is the best reduced cost */
            if( SCIPisSumNegative(scip, bestsolval - pricerdata->dualsolconv[prob]) )
               bestredcost += GCGrelaxGetNIdenticalBlocks(origprob, prob) *
                  (bestsolval - pricerdata->dualsolconv[prob]);
         }

         if( status != SCIP_STATUS_OPTIMAL )
            bestredcostvalid = FALSE;

         nfoundvarsprob = 0;

         for( j = 0; j < nsols && nfoundvarsprob <= pricerdata->maxsolsprob &&
                 (pricetype == GCG_PRICETYPE_REDCOST || nfoundvars < pricerdata->maxvarsroundfarkas)
                 && (pricetype == GCG_PRICETYPE_FARKAS || ((nfoundvars < pricerdata->maxvarsroundredcost || root ) && (nfoundvars < pricerdata->maxvarsroundredcostroot || !root))
                       || pricerdata->onlybest); j++ )
         {
            /* create new variable, compute objective function value and add it to the master constraints and cuts it belongs to */
            SCIP_CALL( createNewMasterVar(scip, solvars[j], solvals[j], nsolvars[j], solisray[j], prob,
                  pricetype == GCG_PRICETYPE_REDCOST, FALSE, &added, NULL) );

            if( added )
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
               pricerdata->nbestsolvars[j], pricerdata->bestsolisray[j], pricerdata->prob[j], FALSE, FALSE, &added, NULL) );
         assert(added);
      }
      pricerdata->nbestsols = 0;
   }

   /** @todo: perhaps solve remaining pricing problems, if only few left? */
   /** @todo: solve all pricing problems all k iterations? */
   /* this makes sure that if a pricing problem has not been solved, the langrangian bound cannot be calculated */
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

   if( pricetype == GCG_PRICETYPE_REDCOST && bestredcostvalid && pricerdata->useinterbounds && duringheurpricing == FALSE)
   {
      assert(lowerbound != NULL);
      GCGpricerPrintInfo(scip, pricerdata, "lower bound = %g, bestredcost = %g\n", SCIPgetLPObjval(scip) + bestredcost, bestredcost);

      *lowerbound = SCIPgetLPObjval(scip) + bestredcost;
   }

   SCIPdebugMessage("%s pricing: found %d new vars\n", (pricetype == GCG_PRICETYPE_REDCOST ? "Redcost" : "Farkas"), nfoundvars);
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


/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitGcg)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   /* get pricerdata */
   pricerdata = SCIPpricerGetData(pricer);

   SCIP_CALL( solversInit(scip, pricerdata) );

   return SCIP_OKAY;
}


/** deinitialization method of variable pricer (called before transformed problem is freed) */
static
SCIP_DECL_PRICEREXIT(pricerExitGcg)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   /* get pricerdata */
   pricerdata = SCIPpricerGetData(pricer);

   SCIP_CALL( solversExit(scip, pricerdata) );

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
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->npointsprob), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nraysprob), pricerdata->npricingprobs) );
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
      pricerdata->npointsprob[i] = 0;
      pricerdata->nraysprob[i] = 0;
   }

   /* alloc memory for arrays of reduced cost */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->dualsolconv), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->score), pricerdata->npricingprobs) );
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
      SCIP_Real* coefs;
      int blocknr;
      int ncoefs;

      blocknr = GCGvarGetBlock(vars[v]);
      coefs = GCGoriginalVarGetCoefs(vars[v]);
      ncoefs = GCGoriginalVarGetNCoefs(vars[v]);

      assert(GCGvarIsOriginal(vars[v]));
      if( blocknr < 0 )
      {
         SCIP_CONS** linkconss;
         SCIP_VAR* newvar;
         linkconss = GCGoriginalVarGetLinkingCons(vars[v]);

         SCIP_CALL( GCGcreateInitialMasterVar(scip, vars[v], &newvar) );
         SCIP_CALL( SCIPaddVar(scip, newvar) );

         SCIP_CALL( GCGoriginalVarAddMasterVar(scip, vars[v], newvar, 1.0) );

         /* add variable in the master to the master constraints it belongs to */
         for( i = 0; i < ncoefs; i++ )
         {
            SCIP_CONS* linkcons;
            assert(!SCIPisZero(scip, coefs[i]));
            SCIP_CALL( SCIPgetTransformedCons(scip, linkconss[i], &linkcons) );

            SCIP_CALL( SCIPaddCoefLinear(scip, linkcons, newvar, coefs[i]) );
         }

         /* we copied a linking variable into the master, add it to the linkcons */
         if( GCGvarIsLinking(vars[v]) )
         {
            SCIP_CONS** linkingconss;
            linkingconss = GCGlinkingVarGetLinkingConss(vars[v]);

            for( i = 0; i < pricerdata->npricingprobs; i++ )
            {
               if( linkingconss[i] != NULL )
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, linkingconss[i], newvar, 1.0) );
               }
            }
         }

         SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

      }
   }

   SCIP_CALL( SCIPhashmapCreate(&(pricerdata->mapcons2idx), SCIPblkmem(scip), 10 * nmasterconss +1) );
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

   if( pricerdata->onlybest && pricerdata->maxvarsroundredcost <= MAXBEST)
   {
      pricerdata->maxbestsols = pricerdata->maxvarsroundredcost;

      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->bestsolvars, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->bestsolvals, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->nbestsolvars, pricerdata->maxbestsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->bestsolisray, pricerdata->maxbestsols) );
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
      pricerdata->bestsolisray = NULL;
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
   if( pricerdata->onlybest && pricerdata->maxbestsols > 0 )
   {
      assert(pricerdata->bestsolvars != NULL);
      assert(pricerdata->bestsolvals != NULL);
      assert(pricerdata->nbestsolvars != NULL);
      assert(pricerdata->bestsolisray != NULL);
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
      SCIPfreeMemoryArray(scip, &pricerdata->bestsolisray);
      SCIPfreeMemoryArray(scip, &pricerdata->redcost);
      SCIPfreeMemoryArray(scip, &pricerdata->prob);


      pricerdata->bestsolvars = NULL;
      pricerdata->bestsolvals = NULL;
      pricerdata->nbestsolvars = NULL;
      pricerdata->bestsolisray = NULL;
      pricerdata->redcost = NULL;
      pricerdata->prob = NULL;
      pricerdata->maxbestsols = 0;
      pricerdata->nbestsols = 0;
   }


   SCIPfreeMemoryArray(scip, &(pricerdata->pricingprobs));
   SCIPfreeMemoryArray(scip, &(pricerdata->dualsolconv));
   SCIPfreeMemoryArray(scip, &(pricerdata->score));
   SCIPfreeMemoryArray(scip, &(pricerdata->permu));
   SCIPfreeMemoryArray(scip, &(pricerdata->solvals));
   SCIPfreeMemoryArray(scip, &(pricerdata->npointsprob));
   SCIPfreeMemoryArray(scip, &(pricerdata->nraysprob));

   //SCIPfreeMemoryArray(scip, &(pricerdata->dualsolconv));

   for( i = 0; i < pricerdata->npricedvars; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &pricerdata->pricedvars[i]) );
   }
   SCIPfreeMemoryArray(scip, &pricerdata->pricedvars);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "calls = %d\n", pricerdata->calls);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "solved sub-MIPs heur = %d\n", pricerdata->solvedsubmipsheur);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "solved sub-MIPs optimal = %d\n", pricerdata->solvedsubmipsoptimal);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "farkas calls = %d, redcost calls = %d\n", pricerdata->farkascalls, pricerdata->redcostcalls);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "time for farkas pricing (total): %f\n", SCIPgetClockTime(scip, pricerdata->farkasclock));
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "time for redcost pricing (total): %f\n", SCIPgetClockTime(scip, pricerdata->redcostclock));
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "time for transformation: %f\n", SCIPgetClockTime(scip, pricerdata->transformclock));
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "time for freeing sub-MIPs: %f\n", SCIPgetClockTime(scip, pricerdata->freeclock));

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

   if( pricerdata->redcostcalls == 0 )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Starting reduced cost pricing...\n");

   /* update number of reduced cost pricing rounds at the current node */
   if( SCIPgetNNodes(scip) == pricerdata->currnodenr )
   {
      pricerdata->nroundsredcost++;
   }
   else
   {
      pricerdata->currnodenr = SCIPgetNNodes(scip);
      pricerdata->nroundsredcost = 0;
   }

   /* if the number of reduced cost pricing rounds at the current node exceeds the limit (and we are not at the root), stop pricing;
    * we always stop pricing, if the maximum number of reduced cost rounds is set to 0
    */
   if( pricerdata->maxroundsredcost == 0 || (pricerdata->nroundsredcost >= pricerdata->maxroundsredcost && pricerdata->currnodenr != 1) )
   {
      SCIPdebugMessage("pricing aborted at node %lld\n", pricerdata->currnodenr);
      return SCIP_OKAY;
   }

   *result = SCIP_SUCCESS;

   /* perform pricing */
   SCIP_CALL( SCIPstartClock(scip, pricerdata->redcostclock) );
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

#define pricerCopyGcg NULL

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
         pricerCopyGcg, pricerFreeGcg, pricerInitGcg, pricerExitGcg,
         pricerInitsolGcg, pricerExitsolGcg, pricerRedcostGcg, pricerFarkasGcg,
         pricerdata) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxsuccessfulmipsredcost",
         "maximal number of pricing mips leading to new variables solved solved in one redcost pricing round",
         &pricerdata->maxsuccessfulmipsredcost, FALSE, DEFAULT_MAXSUCCESSFULMIPSREDCOST, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxvarsroundredcost",
         "maximal number of variables created in one redcost pricing round",
         &pricerdata->maxvarsroundredcost, FALSE, DEFAULT_MAXVARSROUNDREDCOST, 0, INT_MAX,
         paramChgdOnlybestMaxvars, (SCIP_PARAMDATA*)pricerdata) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxvarsroundredcostroot",
         "maximal number of variables created in one redcost pricing round at the root node",
         &pricerdata->maxvarsroundredcostroot, FALSE, DEFAULT_MAXVARSROUNDREDCOSTROOT, 0, INT_MAX,
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

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/abortpricingint",
         "should pricing be aborted due to integral objective function?",
         &pricerdata->abortpricingint, TRUE, DEFAULT_ABORTPRICINGINT, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(pricerdata->origprob, "pricing/masterpricer/abortpricinggap",
         "should pricing be aborted due to small gap between dual bound and RMP objective?",
         &pricerdata->abortpricinggap, TRUE, DEFAULT_ABORTPRICINGGAP, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/useinterbounds",
         "should lagrangean intermediate dual bounds be computed and used?",
         &pricerdata->useinterbounds, TRUE, DEFAULT_USEINTERBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/onlybest",
         "should only the best variables (TRUE) be added in case of a maxvarsround limit or the first ones (FALSE)?",
         &pricerdata->onlybest, TRUE, DEFAULT_ONLYBEST, paramChgdOnlybestMaxvars, (SCIP_PARAMDATA*)pricerdata) );

   SCIP_CALL( SCIPaddRealParam(pricerdata->origprob, "pricing/masterpricer/successfulsubmipsrel",
         "part of the submips that are solved and lead to new variables before pricing round is aborted? (1.0 = solve all pricing MIPs)",
         &pricerdata->successfulmipsrel, FALSE, DEFAULT_SUCCESSFULMIPSREL, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(pricerdata->origprob, "pricing/masterpricer/mipsrelredcostroot",
         "part of the submips that are solved before redcost pricing round is aborted at the root node, if variables have been found yed? (1.0 = solve all pricing MIPs)",
         &pricerdata->mipsrelredcostroot, FALSE, DEFAULT_MIPSRELREDCOSTROOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(pricerdata->origprob, "pricing/masterpricer/mipsrelredcost",
         "part of the submips that are solved before redcost pricing round is aborted, if variables have been found yed? (1.0 = solve all pricing MIPs)",
         &pricerdata->mipsrelredcost, FALSE, DEFAULT_MIPSRELREDCOST, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(pricerdata->origprob, "pricing/masterpricer/mipsrelfarkas",
         "part of the submips that are solved before Farkas pricing round is aborted, if variables have been found yed? (1.0 = solve all pricing MIPs)",
         &pricerdata->mipsrelfarkas, FALSE, DEFAULT_MIPSRELFARKAS, 0.0, 1.0, NULL, NULL) );

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

   return pricerdata->pricedvars;
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

   return pricerdata->npricedvars;
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


/** transfers a primal solution of the original problem into the master variable space,
 *  i.e. creates one master variable for each block and adds the solution to the master problem  */
SCIP_RETCODE GCGpricerTransOrigSolToMasterVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             origsol             /**< the solution that should be transferred */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
#if 1
   SCIP_SOL* mastersol;
#endif
   SCIP_VAR* newvar;
   SCIP* origprob;
   SCIP_Bool added;
   int prob;
   int i;

   SCIP_VAR** origvars;
   SCIP_Real* origsolvals;
   int norigvars;

   SCIP_VAR*** pricingvars;
   SCIP_Real** pricingvals;
   int* npricingvars;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   /* now compute coefficients of the master variables in the master constraint */
   origvars = SCIPgetVars(origprob);
   norigvars = SCIPgetNVars(origprob);

   /* allocate memory for storing variables and solution values from the solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &origsolvals, norigvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pricingvars, pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pricingvals, pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &npricingvars, pricerdata->npricingprobs) );

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      npricingvars[i] = 0;
      if(pricerdata->pricingprobs[i] == NULL)
      {
         continue;
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &pricingvars[i], SCIPgetNVars(pricerdata->pricingprobs[i])) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pricingvals[i], SCIPgetNVars(pricerdata->pricingprobs[i])) );
   }

   /* get solution values */
   SCIP_CALL( SCIPgetSolVals(scip, origsol, norigvars, origvars, origsolvals) );

#if 1
   SCIP_CALL( SCIPcreateSol(scip, &mastersol, NULL) );
#endif
   /* store variables and solutions into arrays */
   for( i = 0; i < norigvars; i++ )
   {
      int blocknr;
      assert(GCGvarIsOriginal(origvars[i]));
      blocknr = GCGvarGetBlock(origvars[i]);
      assert(GCGoriginalVarGetPricingVar(origvars[i]) != NULL || blocknr < 0);

      if( blocknr >= 0 )
      {
         prob = blocknr;
         if(pricerdata->pricingprobs[prob] == NULL)
         {
            continue;
         }

         if( !SCIPisZero(scip, origsolvals[i]) )
         {
            pricingvars[prob][npricingvars[prob]] = GCGoriginalVarGetPricingVar(origvars[i]);
            pricingvals[prob][npricingvars[prob]] = origsolvals[i];
            npricingvars[prob]++;
         }
      }
#if 1
      else
      {
         assert(GCGoriginalVarGetNMastervars(origvars[i]) == 1);
         assert(GCGoriginalVarGetMastervars(origvars[i])[0] != NULL);
         SCIP_CALL( SCIPsetSolVal(scip, mastersol, GCGoriginalVarGetMastervars(origvars[i])[0], origsolvals[i]) );
      }
#endif
   }

   /* create variables in the master problem */
   for( prob = 0; prob < pricerdata->npricingprobs; prob++ )
   {
      if(pricerdata->pricingprobs[prob] == NULL)
      {
         continue;
      }
      SCIP_CALL( createNewMasterVar(scip, pricingvars[prob], pricingvals[prob], npricingvars[prob], FALSE, prob, FALSE, TRUE, &added, &newvar) );
      assert(added);
#if 1
      SCIP_CALL( SCIPsetSolVal(scip, mastersol, newvar, 1.0 * GCGrelaxGetNIdenticalBlocks(pricerdata->origprob, prob)) );
#endif
   }

#if 1
   SCIP_CALL( SCIPtrySolFree(scip, &mastersol, TRUE, TRUE, TRUE, TRUE, &added) );
#endif

   /* free memory for storing variables and solution values from the solution */

   for( i = pricerdata->npricingprobs - 1; i>= 0; i-- )
   {
      if(pricerdata->pricingprobs[i] == NULL)
      {
         continue;
      }
      SCIPfreeBufferArray(scip, &pricingvals[i]);
      SCIPfreeBufferArray(scip, &pricingvars[i]);
   }

   SCIPfreeBufferArray(scip, &npricingvars);
   SCIPfreeBufferArray(scip, &pricingvals);
   SCIPfreeBufferArray(scip, &pricingvars);
   SCIPfreeBufferArray(scip, &origsolvals);



   return SCIP_OKAY;
}
