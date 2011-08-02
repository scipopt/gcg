/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"
//#define DEBUG_PRICING_ALL_OUTPUT
//#define SCIP_DEBUG
/**@file   solver_mip.c
 * @brief  mip solver for pricing problems
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <pthread.h>
#include <scip_misc.h>
#include "solver_mip.h"
#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"
#include "type_solver.h"
#include "struct_solverinfo.h"
#include "struct_vardata.h"
#include "pricer_gcg.h"



#define SOLVER_NAME          "mip"
#define SOLVER_DESC          "mip solver for pricing problems"
#define SOLVER_PRIORITY      0

#define DEFAULT_CHECKSOLS TRUE



struct GCG_SolData
{
   SCIP_Real** solvals;
   SCIP_VAR*** solvars;
   SCIP_Real*  tmpsolvals;
   int* nsolvars;
   SCIP_Bool* solisray;
   int nsols;
   int maxvars;
};

typedef struct GCG_SolData GCG_SOLDATA;

struct GCG_SolverData
{
   SCIP* origprob;
   GCG_SOLDATA** soldata;
   GCG_SOLVERINFO *solverinfo;

   SCIP_Bool checksols;
};



/* ensures size of solution arrays */
static
SCIP_RETCODE ensureSizeSolvars(
   SCIP*                 scip,
   GCG_SOLDATA*       soldata,
   int                   nsols
   )
{
   int i;

   assert(scip != NULL);
   assert(soldata != NULL);

   if( soldata->nsols < nsols )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(soldata->nsolvars), nsols) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(soldata->solisray), nsols) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(soldata->solvars), nsols) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(soldata->solvals), nsols) );

      for( i = soldata->nsols; i < nsols; i++ )
      {
         soldata->nsolvars[i] = 0;
         SCIP_CALL( SCIPallocMemoryArray(scip, &(soldata->solvars[i]), soldata->maxvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(soldata->solvals[i]), soldata->maxvars) );
      }

      soldata->nsols = nsols;

   }

   return SCIP_OKAY;
}

/** checks whether the given solution is equal to one of the former solutions in the sols array */
static
SCIP_RETCODE checkSolNew(
   SCIP*                 scip,
   SCIP*                 pricingprob,
   SCIP_SOL**            sols,
   int                   idx,
   SCIP_Bool*            isnew
   )
{
   SCIP* origprob;
   SCIP_VAR** probvars;
   int nprobvars;
   SCIP_Real* newvals;

   int s;
   int i;

   assert(scip != NULL);
   assert(pricingprob != NULL);
   assert(sols != NULL);
   assert(sols[idx] != NULL);
   assert(isnew != NULL);

   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   probvars = SCIPgetVars(pricingprob);
   nprobvars = SCIPgetNVars(pricingprob);

   *isnew = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &newvals, nprobvars) );

   SCIP_CALL( SCIPgetSolVals(pricingprob, sols[idx], nprobvars, probvars, newvals) );

   for( s = 0; s < idx && *isnew == TRUE; s++ )
   {
      assert(sols[s] != NULL);
      /* TODO: ensure that the solutions are sorted  */
      /*assert(SCIPisLE(scip, SCIPgetSolOrigObj(pricingprob, sols[s]), SCIPgetSolOrigObj(pricingprob, sols[idx])) 
       *|| ABS(SCIPgetSolOrigObj(pricingprob, sols[s])) > 1e+15 * SCIPepsilon(pricingprob));*/
      if( !SCIPisEQ(scip, SCIPgetSolOrigObj(pricingprob, sols[s]), SCIPgetSolOrigObj(pricingprob, sols[idx])) )
         continue;

      if( SCIPsolGetOrigin(sols[s]) != SCIP_SOLORIGIN_ORIGINAL && SCIPsolGetOrigin(sols[idx]) != SCIP_SOLORIGIN_ORIGINAL )
         continue;

      for( i = 0; i < nprobvars; i++ )
         if( !SCIPisEQ(scip, SCIPgetSolVal(pricingprob, sols[s], probvars[i]), newvals[i]) )
            break;

      if( i == nprobvars )
         *isnew = FALSE;
   }

   SCIPfreeBufferArray(scip, &newvals);

   return SCIP_OKAY;
}



/*
 * Callback methods for pricing problem solver
 */

static
GCG_DECL_SOLVERFREE(solverFreeMip)
{
   GCG_SOLVERDATA* solverdata;

   assert(scip != NULL);
   assert(solver != NULL);

   solverdata = GCGpricerGetSolverdata(scip, solver);
   assert(solverdata != NULL);

   SCIPfreeMemory(scip, &solverdata);

   GCGpricerSetSolverdata(scip, solver, NULL);

   return SCIP_OKAY;
}

static
GCG_DECL_SOLVERINITSOL(solverInitsolMip)
{
   GCG_SOLVERDATA* solverdata;
   int i;
   int j;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(solverinfo != NULL);

   solverdata = GCGpricerGetSolverdata(scip, solver);
   assert(solverdata != NULL);

   solverdata->solverinfo = solverinfo;
   SCIP_CALL( SCIPallocMemoryArray(scip, &solverdata->soldata, GCGrelaxGetNPricingprobs(solverdata->origprob)));
   SCIPdebugMessage("Pricingprobs: %d\n",  GCGrelaxGetNPricingprobs(solverdata->origprob));
   assert(solverdata->soldata != NULL);
   // solverdata->maxvars = -1;
   for( i = 0; i < GCGrelaxGetNPricingprobs(solverdata->origprob); i++ )
   {
      SCIP_CALL(SCIPallocMemory(scip, &solverdata->soldata[i]));
      assert(solverdata->soldata[i] != NULL);
//      if( SCIPgetNVars(GCGrelaxGetPricingprob(solverdata->origprob, i)) > solverdata->maxvars )
      solverdata->soldata[i]->maxvars = SCIPgetNVars(GCGrelaxGetPricingprob(solverdata->origprob, i));
   }
   for( i = 0; i < GCGrelaxGetNPricingprobs(solverdata->origprob); i++ )
   {
      assert(solverdata->soldata[i] != NULL);
      solverdata->soldata[i]->nsols = 10;
      SCIP_CALL( SCIPallocMemoryArray(scip, &(solverdata->soldata[i]->nsolvars), solverdata->soldata[i]->nsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(solverdata->soldata[i]->solisray), solverdata->soldata[i]->nsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(solverdata->soldata[i]->solvars), solverdata->soldata[i]->nsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(solverdata->soldata[i]->solvals), solverdata->soldata[i]->nsols) );

      for( j = 0; j < solverdata->soldata[i]->nsols; j++ )
      {
         solverdata->soldata[i]->nsolvars[j] = 0;
         SCIP_CALL( SCIPallocMemoryArray(scip, &(solverdata->soldata[i]->solvars[j]), solverdata->soldata[i]->maxvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(solverdata->soldata[i]->solvals[j]), solverdata->soldata[i]->maxvars) );
      }

      SCIP_CALL( SCIPallocMemoryArray(scip, &(solverdata->soldata[i]->tmpsolvals), solverdata->soldata[i]->maxvars) );
   }
   return SCIP_OKAY;
}

static
GCG_DECL_SOLVEREXITSOL(solverExitsolMip)
{
   GCG_SOLVERDATA* solverdata;
   int i;
   int j;

   assert(scip != NULL);
   assert(solver != NULL);

   solverdata = GCGpricerGetSolverdata(scip, solver);
   assert(solverdata != NULL);
   for( i = 0; i < GCGrelaxGetNPricingprobs(solverdata->origprob); i++ )
   {
      for( j = 0; j < solverdata->soldata[i]->nsols; j++ )
      {
         SCIPfreeMemoryArray(scip, &(solverdata->soldata[i]->solvars[j]));
         SCIPfreeMemoryArray(scip, &(solverdata->soldata[i]->solvals[j]));
      }

      SCIPfreeMemoryArray(scip, &(solverdata->soldata[i]->tmpsolvals));

      SCIPfreeMemoryArray(scip, &(solverdata->soldata[i]->nsolvars));
      SCIPfreeMemoryArray(scip, &(solverdata->soldata[i]->solisray));
      SCIPfreeMemoryArray(scip, &(solverdata->soldata[i]->solvars));
      SCIPfreeMemoryArray(scip, &(solverdata->soldata[i]->solvals));
      SCIPfreeMemory(scip, &solverdata->soldata[i]);

   }
   SCIPfreeMemoryArray(scip, &solverdata->soldata);
   return SCIP_OKAY;
}

#define solverInitMip NULL
#define solverExitMip NULL


static
GCG_DECL_SOLVERSOLVE(solverSolveMip)
{
   GCG_SOLVERDATA* solverdata;
   SCIP_SOL** probsols;
   int nprobsols;

   SCIP_VAR** probvars;
   int nprobvars;

   SCIP_Bool newsol;

   int s;
   int i;

#ifdef DEBUG_PRICING_ALL_OUTPUT
   //if( pricetype == GCG_PRICETYPE_REDCOST )
   {
      //char probname[SCIP_MAXSTRLEN];
      //(void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "pricingmip_%d_%d_vars.lp", probnr, SCIPgetNVars(scip));
      //SCIP_CALL( SCIPwriteOrigProblem(pricingprob, probname, NULL, FALSE) );
      
      SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", SCIP_VERBLEVEL_HIGH) );
      
      //SCIP_CALL( SCIPwriteParams(pricingprob, "pricing.set", TRUE, TRUE) );
      
   }
#endif
   solverdata = GCGpricerGetSolverdata(scip, solver);
   assert(solverdata != NULL);
   assert(solverdata->soldata != NULL);
   assert(solverdata->solverinfo != NULL);

//   pthread_mutex_lock(&solverdata->solverinfo->update_count);
//   assert(solverdata->solverinfo->count >= 0);
//   ++(solverdata->solverinfo->count);
//   pthread_mutex_unlock(&solverdata->solverinfo->update_count);

   SCIP_CALL( SCIPtransformProb(pricingprob) );

   /* presolve the pricing submip */
   if( SCIPgetStage(pricingprob) < SCIP_STAGE_PRESOLVING )
   {
      SCIP_CALL( SCIPpresolve(pricingprob) );
   }

#if 0
   SCIP_CALL( SCIPwriteOrigProblem(pricingprob, "pricing.lp", NULL, FALSE) );
   SCIP_CALL( SCIPwriteOrigProblem(pricingprob, "pricing.cip", NULL, FALSE) );
   SCIP_CALL( SCIPwriteParams(pricingprob, "pricing.set", FALSE, FALSE) );
#endif
   /* solve the pricing submip */
   SCIP_CALL( SCIPsolve(pricingprob) );

   /* all SCIP statuses handled so far */
   assert( SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_GAPLIMIT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT 
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFEASIBLE
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT 
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_UNBOUNDED 
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFORUNBD );
   /** @todo handle userinterrupt: set result pointer (and lowerbound), handle other solution states */

   //printf("MIP pricing solver: status = %d\n", SCIPgetStatus(pricingprob));

   if( SCIPgetStatus(pricingprob) == SCIP_STATUS_UNBOUNDED 
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFORUNBD )
   {
      /* the pricing problem was declared to be (infeasible or) unbounded, but SCIP did not compute a primal ray;
       * this occurs when presolving detected (infeasibility or) unboundedness; since we need a primal ray to create 
       * the corresponding variable, we disable presolving and resolve the problem to get the primal ray out of the LP */
      if( !SCIPhasPrimalRay(pricingprob) )
      {
         SCIP_CALL( SCIPfreeTransform(pricingprob) );
         
         SCIP_CALL( SCIPsetIntParam(pricingprob, "presolving/maxrounds", 0) );
         SCIP_CALL( SCIPtransformProb(pricingprob) );
         /* solve the pricing submip */
         SCIP_CALL( SCIPsolve(pricingprob) );
      }
      assert(SCIPhasPrimalRay(pricingprob)
         || SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT 
         || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT );
   }

   //printf("MIP pricing solver: status = %d\n", SCIPgetStatus(pricingprob));

   if( SCIPgetStatus(pricingprob) == SCIP_STATUS_UNBOUNDED 
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFORUNBD )
   {
      assert(SCIPhasPrimalRay(pricingprob));

      probvars  = SCIPgetOrigVars(pricingprob);
      nprobvars = SCIPgetNOrigVars(pricingprob);

      solverdata->soldata[probnr]->nsolvars[0] = 0;

      /* store the primal ray values */
      for( i = 0; i < nprobvars; i++ )
      {
         if( SCIPisZero(scip, SCIPgetPrimalRayVal(pricingprob, probvars[i])) )
            continue;

         assert(!SCIPisInfinity(scip, SCIPgetPrimalRayVal(pricingprob, probvars[i])) && !SCIPisInfinity(scip, -SCIPgetPrimalRayVal(pricingprob, probvars[i])));

         solverdata->soldata[probnr]->solvars[0][solverdata->soldata[probnr]->nsolvars[0]] = probvars[i];
         solverdata->soldata[probnr]->solvals[0][solverdata->soldata[probnr]->nsolvars[0]] = SCIPgetPrimalRayVal(pricingprob, probvars[i]);
         solverdata->soldata[probnr]->nsolvars[0]++;

         SCIPdebugMessage("%s: %g\n", SCIPvarGetName(probvars[i]), SCIPgetPrimalRayVal(pricingprob, probvars[i]));
      }
      solverdata->soldata[probnr]->solisray[0] = TRUE;
      *solvars = solverdata->soldata[probnr]->solvars;
      *solvals = solverdata->soldata[probnr]->solvals;
      *nsolvars = solverdata->soldata[probnr]->nsolvars;
      *solisray = solverdata->soldata[probnr]->solisray;
      *nsols = 1;
      *result = SCIP_STATUS_UNBOUNDED;

      SCIPdebugMessage("pricingproblem has an unbounded ray!\n");
   }
   else if( SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT )
   {
      *solvars = solverdata->soldata[probnr]->solvars;
      *solvals = solverdata->soldata[probnr]->solvals;
      *nsolvars = solverdata->soldata[probnr]->nsolvars;
      *solisray = solverdata->soldata[probnr]->solisray;
      *nsols = 0;
      *result = SCIPgetStatus(pricingprob);
   } 
   else
   {
      /* get variables of the pricing problem */
      probvars = SCIPgetOrigVars(pricingprob);
      nprobvars = SCIPgetNOrigVars(pricingprob);

      nprobsols = SCIPgetNSols(pricingprob);
      probsols = SCIPgetSols(pricingprob);

      *nsols = 0;
      PTHREAD_CALL(pthread_mutex_lock(&solverdata->solverinfo->access_masterscip));
      SCIP_CALL( ensureSizeSolvars(scip, solverdata->soldata[probnr], nprobsols) ); //TODO: BUG: make sure mutex gets unlocked in every case!
      PTHREAD_CALL(pthread_mutex_unlock(&solverdata->solverinfo->access_masterscip));
      for( s = 0; s < nprobsols; s++ )
      {
         SCIP_Bool feasible;
         //printf("MIP pricing solver: status = %d\n", SCIPgetStatus(pricingprob));
         //printf("Objective value of solution: %g\n", SCIPgetSolOrigObj(pricingprob, probsols[s]));
         if( SCIPisInfinity(pricingprob, -SCIPgetSolOrigObj(pricingprob, probsols[s])) )
         {
           SCIPdebugMessage("unbounded solution\n");
         }
         //#ifndef NDEBUG
         SCIP_CALL( SCIPcheckSolOrig(pricingprob, probsols[s], &feasible, TRUE, TRUE) );
         //assert(feasible);
         //#endif

         /* check whether the solution is equal to one of the previous solutions */
         if( solverdata->checksols )
         {
            PTHREAD_CALL(pthread_mutex_lock(&solverdata->solverinfo->access_masterscip));
            SCIP_CALL( checkSolNew(scip, pricingprob, probsols, s, &newsol) );
            PTHREAD_CALL(pthread_mutex_unlock(&solverdata->solverinfo->access_masterscip));

            if( !newsol )
               continue;
         }

         solverdata->soldata[probnr]->nsolvars[*nsols] = 0;
         solverdata->soldata[probnr]->solisray[*nsols] = FALSE;

         SCIP_CALL( SCIPgetSolVals(pricingprob, probsols[s], nprobvars, probvars, solverdata->soldata[probnr]->tmpsolvals) );

#if 0         
         /* for integer variables, round the solution values */
         for( i = 0; i < nprobvars; i++ )
         {
            if( SCIPisZero(scip, solverdata->soldata[probnr]->tmpsolvals[i]) )
               continue;
            if( SCIPvarGetType(probvars[i]) != SCIP_VARTYPE_CONTINUOUS )
            {
               assert(SCIPisEQ(scip, solverdata->soldata[probnr]->tmpsolvals[i], SCIPfeasFloor(scip, solverdata->soldata[probnr]->tmpsolvals[i])));
               solverdata->soldata[probnr]->tmpsolvals[i] = SCIPfeasFloor(scip, solverdata->soldata[probnr]->tmpsolvals[i]);
            }
         }
#endif

         /* try to handle nearly unbounded solutions, that are only finite due to numerical troubles:
          * ATTENTION: some of the methods below are called where this is normally not allowed in SCIP
          */
#if 0
         {
            SCIP_SOL* tmpsol;
            tmpsol = NULL;

            /* check the solution values for infinity */
            for( i = 0; i < nprobvars; i++ )
            {
               if( SCIPisInfinity(scip, solverdata->soldata[probnr]->tmpsolvals[i]) )
               {
                  SCIP_VARDATA* vardata;

                  if( tmpsol == NULL )
                  {
                     SCIP_CALL( SCIPcreateSolCopy(pricingprob, &tmpsol, probsols[s]) );
                  }

                  SCIPwarningMessage("found solution for pricing problem with variable, that has solution value +infinity - try to reduce this value\n");
                  vardata = SCIPvarGetData(probvars[i]);
                  assert(vardata->vartype == GCG_VARTYPE_PRICING);
                  assert(vardata->data.pricingvardata.origvars != NULL);
                  assert(vardata->data.pricingvardata.origvars[0] != NULL);

                  SCIP_CALL( SCIPsetSolVal(pricingprob, tmpsol, probvars[i], 
                        SCIPinfinity(scip) / (2 * SCIPvarGetObj(vardata->data.pricingvardata.origvars[0]))) );

                  SCIP_CALL( SCIPcheckSol(pricingprob, tmpsol, TRUE, TRUE, TRUE, TRUE, &feasible) );
                  if( feasible )
                  {
                     solverdata->soldata[probnr]->tmpsolvals[i] = SCIPinfinity(scip) / (2 * SCIPvarGetObj(vardata->data.pricingvardata.origvars[0]));
                     SCIPwarningMessage("--> reduced value of variable to %g\n", solverdata->soldata[probnr]->tmpsolvals[i]);
                  }
                  else
                  {
                     SCIPwarningMessage("--> reducing value of variable to %g caused infeasibility\n", solverdata->soldata[probnr]->tmpsolvals[i]);
                  }
               }
            }
         
            if( tmpsol != NULL )
            {
               SCIP_CALL( SCIPfreeSol(pricingprob, &tmpsol) );
            }
         }
#endif
        

         /* store the solution values */
         for( i = 0; i < nprobvars; i++ )
         {
            if( SCIPisZero(scip, solverdata->soldata[probnr]->tmpsolvals[i]) )
               continue;

            solverdata->soldata[probnr]->solvars[*nsols][solverdata->soldata[probnr]->nsolvars[*nsols]] = probvars[i];
            solverdata->soldata[probnr]->solvals[*nsols][solverdata->soldata[probnr]->nsolvars[*nsols]] = solverdata->soldata[probnr]->tmpsolvals[i];
            solverdata->soldata[probnr]->nsolvars[*nsols]++;
         }
         
         *nsols = *nsols + 1;
      }
         
      *solvars = solverdata->soldata[probnr]->solvars;
      *solvals = solverdata->soldata[probnr]->solvals;
      *nsolvars = solverdata->soldata[probnr]->nsolvars;
      *solisray = solverdata->soldata[probnr]->solisray;

      *result = SCIP_STATUS_OPTIMAL;
      SCIPdebugMessage("pricingproblem found %d sols!\n", *nsols);
   }

#ifdef DEBUG_PRICING_ALL_OUTPUT
   SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", 0) );
   SCIP_CALL( SCIPprintStatistics(pricingprob, NULL) );
#endif

   PTHREAD_CALL(pthread_mutex_lock(&solverdata->solverinfo->update_count));
   --(solverdata->solverinfo->count);
   solverdata->solverinfo->queue[solverdata->solverinfo->nqueueentries] = probnr;
   ++(solverdata->solverinfo->nqueueentries);
   SCIPdebugMessage("Adding thread %d to the queue, we have now %d entries (%d sols).\n", probnr, solverdata->solverinfo->nqueueentries, solverdata->soldata[probnr]->nsols );
   assert(solverdata->solverinfo->count >= 0);
   PTHREAD_CALL(pthread_mutex_unlock(&solverdata->solverinfo->update_count));
   PTHREAD_CALL(pthread_cond_signal(&solverdata->solverinfo->update_cond));
   return SCIP_OKAY;
}


static
GCG_DECL_SOLVERSOLVEHEUR(solverSolveHeurMip)
{
   GCG_SOLVERDATA* solverdata;
   SCIP_SOL** probsols;
   int nprobsols;

   SCIP_VAR** probvars;
   int nprobvars;

   SCIP_Bool newsol;

   int s;
   int i;

#ifdef DEBUG_PRICING_ALL_OUTPUT
   SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", SCIP_VERBLEVEL_HIGH) );
#endif

   solverdata = GCGpricerGetSolverdata(scip, solver);
   assert(solverdata != NULL);

   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/stallnodes", 100) );
   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/nodes", 1000) );
   SCIP_CALL( SCIPsetRealParam(pricingprob, "limits/gap", 0.2) );
   //SCIP_CALL( SCIPsetIntParam(pricingprob, "limits/bestsol", 5) );

   SCIP_CALL( SCIPtransformProb(pricingprob) );

   /* presolve the pricing submip */
   if( SCIPgetStage(pricingprob) < SCIP_STAGE_PRESOLVING )
   {
      SCIP_CALL( SCIPpresolve(pricingprob) );
   }

   /* solve the pricing submip */
   SCIP_CALL( SCIPsolve(pricingprob) );
   
   /* so far, the pricing problem should be solved to optimality */
   assert( SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_GAPLIMIT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT 
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_INFEASIBLE
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT );
   /** @todo handle userinterrupt: set result pointer (and lowerbound), handle other solution states */
   
   if( SCIPgetStatus(pricingprob) == SCIP_STATUS_USERINTERRUPT
      || SCIPgetStatus(pricingprob) == SCIP_STATUS_TIMELIMIT )
   {
      *solvars = solverdata->soldata[probnr]->solvars;
      *solvals = solverdata->soldata[probnr]->solvals;
      *nsolvars = solverdata->soldata[probnr]->nsolvars;
      *solisray = solverdata->soldata[probnr]->solisray;
      *nsols = 0;
      *result = SCIP_STATUS_UNKNOWN;
   } 
   else
   {
      /* get variables of the pricing problem */
      probvars = SCIPgetOrigVars(pricingprob);
      nprobvars = SCIPgetNOrigVars(pricingprob);

      nprobsols = SCIPgetNSols(pricingprob);
      probsols = SCIPgetSols(pricingprob);

      *nsols = 0;

      SCIP_CALL( ensureSizeSolvars(scip, solverdata->soldata[probnr], nprobsols) );

      for( s = 0; s < nprobsols; s++ )
      {
         //#ifndef NDEBUG
         SCIP_Bool feasible;
         SCIP_CALL( SCIPcheckSolOrig(pricingprob, probsols[s], &feasible, TRUE, TRUE) );
         if( !feasible )
            SCIPdebugMessage("heur: %s\n", SCIPheurGetName(SCIPsolGetHeur(probsols[s])));
         assert(feasible);
         //#endif
         
         if( solverdata->checksols )
         {
            SCIP_CALL( checkSolNew(scip, pricingprob, probsols, s, &newsol) );
            
            if( !newsol )
               continue;
         }

         solverdata->soldata[probnr]->nsolvars[*nsols] = 0;
         solverdata->soldata[probnr]->solisray[*nsols] = FALSE;

         SCIP_CALL( SCIPgetSolVals(pricingprob, probsols[s], nprobvars, probvars, solverdata->soldata[probnr]->tmpsolvals) );
         
         /* for integer variable, round the solvals */
         for( i = 0; i < nprobvars; i++ )
         {
            if( SCIPisZero(scip, solverdata->soldata[probnr]->tmpsolvals[i]) )
               continue;

            solverdata->soldata[probnr]->solvars[*nsols][solverdata->soldata[probnr]->nsolvars[*nsols]] = probvars[i];

            if( SCIPvarGetType(probvars[i]) != SCIP_VARTYPE_CONTINUOUS )
            {
               assert(SCIPisEQ(scip, solverdata->soldata[probnr]->tmpsolvals[i], SCIPfeasFloor(scip, solverdata->soldata[probnr]->tmpsolvals[i])));
               solverdata->soldata[probnr]->solvals[*nsols][solverdata->soldata[probnr]->nsolvars[*nsols]]
                  = SCIPfeasFloor(scip, solverdata->soldata[probnr]->tmpsolvals[i]);
            }
            else
               solverdata->soldata[probnr]->solvals[*nsols][solverdata->soldata[probnr]->nsolvars[*nsols]] = solverdata->soldata[probnr]->tmpsolvals[i];

            solverdata->soldata[probnr]->nsolvars[*nsols]++;
         }
         
         *nsols = *nsols + 1;
      }
         
      *solvars = solverdata->soldata[probnr]->solvars;
      *solvals = solverdata->soldata[probnr]->solvals;
      *nsolvars = solverdata->soldata[probnr]->nsolvars;
      *solisray = solverdata->soldata[probnr]->solisray;

      *result = SCIP_STATUS_OPTIMAL;
      //printf("Pricingprob %d has found %d sols!\n", prob, nsols);
   }

#ifdef DEBUG_PRICING_ALL_OUTPUT
   //if( pricetype == GCG_PRICETYPE_REDCOST )
   {
      SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", 0) );
      SCIP_CALL( SCIPprintStatistics(pricingprob, NULL) );
   }
#endif

   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/stallnodes", -1) );
   SCIP_CALL( SCIPsetLongintParam(pricingprob, "limits/nodes", -1) ); 
   SCIP_CALL( SCIPsetRealParam(pricingprob, "limits/gap", 0.0) ); 
   SCIP_CALL( SCIPsetIntParam(pricingprob, "limits/bestsol", -1) ); 

   return SCIP_OKAY;
}

/** creates the most infeasible LP braching rule and includes it in SCIP */
SCIP_RETCODE GCGincludeSolverMip(
   SCIP*                 scip                /**< SCIP data structure */
   )
{   
   GCG_SOLVERDATA* data;

   SCIP_CALL( SCIPallocMemory( scip, &data) );
   data->solverinfo = NULL;
   data->origprob = GCGpricerGetOrigprob(scip);

   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY, 
         solverSolveMip, solverSolveHeurMip, solverFreeMip, solverInitMip, solverExitMip,
         solverInitsolMip, solverExitsolMip, data) );

   SCIP_CALL( SCIPaddBoolParam(data->origprob, "pricingsolver/mip/checksols",
         "should solutions of the pricing MIPs be checked for duplicity?",
         &data->checksols, TRUE, DEFAULT_CHECKSOLS, NULL, NULL) );


   return SCIP_OKAY;
}
