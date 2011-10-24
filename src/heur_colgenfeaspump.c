/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_colgenfeaspump.c
 * @ingroup PRIMALHEURISTICS
 * @brief  column generation based feasibility pump primal heuristic
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* toggle debug mode */
// #define SCIP_DEBUG

#include <assert.h>
#include <string.h>

#include "heur_colgenfeaspump.h"
#include "pricer_gcg.h"
#include "pub_gcgvar.h"
#include "relax_gcg.h"

#include "scip/scipdefplugins.h"


#define HEUR_NAME             "colgenfeaspump"
#define HEUR_DESC             "column generation based feasibility pump"
#define HEUR_DISPCHAR         'G'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXLPITERQUOT    0.01   /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS     1000   /**< additional number of allowed LP iterations */
#define DEFAULT_CYCLELENGTH      20
#define DEFAULT_MAXLOOPS         100    /**< maximal number of pumping rounds */
#define DEFAULT_MAXSTALLLOOPS    10     /**< maximal number of pumping rounds without fractionality improvement (-1: no limit) */
#define DEFAULT_OBJFACTOR        0.95
#define DEFAULT_SHIFTRATE        0.05   /**< percentage of variables to be shifted in case of a 1-cycle */

#define MINLPITER                5000  /**< minimal number of LP iterations allowed in each LP solving call */




/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   /* parameters */
   int               cyclelength;
   SCIP_Real         maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   int               maxlpiterofs;       /**< additional number of allowed LP iterations */
   int               maxloops;           /**< maximal number of pumping rounds */
   int               maxstallloops;      /**< maximal number of pumping rounds without fractionality improvement (-1: no limit) */
   SCIP_Real         objfactor;
   SCIP_Real         shiftrate;          /**< percentage of variables to be shifted in case of a 1-cycle */

   /* statistics */
   SCIP_Longint      nlpiterations;      /**< number of LP iterations used in this heuristic */
   int               nsuccess;           /**< number of runs that produced at least one feasible solution */
};




/*
 * Local methods
 */

#if 0
/** set objective values for master diving LP */
static
SCIP_RETCODE setMasterDivingObjectives(
   SCIP*                 scip,
   SCIP_HEURDATA*        heurdata
   )
{
   SCIP* masterprob;

   SCIP_VAR** origvars;
   int norigvars;
   SCIP_VAR** mastervars;
   SCIP_Real* masterobjs;
   int nmastervars;

   int i;

   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   assert(heurdata != NULL);

   SCIP_CALL( SCIPgetVarsData(scip, &origvars, &norigvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   assert(origvars != NULL);
   assert(norigvars > 0);
   assert(mastervars != NULL);
   assert(nmastervars > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &masterobjs, nmastervars) );

   GCGrelaxTransformOrigvalsToMastervals(scip, origvars, heurdata->origobjs, norigvars, mastervars, masterobjs, nmastervars);

   for( i = 0; i < nmastervars; i++ )
   {
      SCIP_CALL( SCIPchgVarObjDive(masterprob, mastervars[i], masterobjs[i]) );
   }

   SCIPfreeBufferArray(scip, &masterobjs);

   return SCIP_OKAY;
}
#endif

/** solve a diving LP on the master problem */
static
SCIP_RETCODE performDivingOnMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   maxlpiterations,    /**< maximum number of LP iterations allowed */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of used LP iterations */
   SCIP_Bool*            lperror            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
//   SCIP_Bool*            cutoff              /**< pointer to store whether the diving direction is infeasible */
   )
{
   SCIP* masterscip;
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_Bool feasible;
   SCIP_Longint oldnlpiters;
//   SCIP_Longint nodelimit;

   assert(scip != NULL);

   masterscip = GCGrelaxGetMasterprob(scip);
   assert(masterscip != NULL);

//   SCIP_CALL( SCIPgetLongintParam(masterscip, "limits/nodes", &nodelimit) );
//   SCIP_CALL( SCIPsetLongintParam(masterscip, "limits/nodes", nodelimit + 1) );

   //printf("before LP solving\n");

   oldnlpiters = SCIPgetNLPIterations(masterscip);
   SCIP_CALL( SCIPsolveDiveLP(masterscip, maxlpiterations, lperror) );
   lpsolstat = SCIPgetLPSolstat(masterscip);

   //printf("after LP solving\n");

//   SCIP_CALL( SCIPsetLongintParam(masterscip, "limits/nodes", nodelimit) );

   *nlpiterations = SCIPgetNLPIterations(masterscip) - oldnlpiters;

   //printf("lperror = %d\n", (*lperror));

   if( !(*lperror) )
   {
      //printf("lpsolstat = %d, isRelax = %d\n", lpsolstat, SCIPisLPRelax(masterscip));
      /* get LP solution status */
//      *cutoff = lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT || lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE;
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )//&& SCIPisLPRelax(masterscip) )
         SCIP_CALL( GCGrelaxUpdateCurrentSol(scip, &feasible) );
   }
   else
   {
      SCIPdebugMessage("something went wrong, an lp error occurred\n");
   }

   return SCIP_OKAY;
}

/** solve the subproblems with a distance objective function */
static
SCIP_RETCODE solvePricingProblems(
      SCIP*             scip,
      SCIP_Real         alpha,
      SCIP_SOL*         relaxsol,
      SCIP_SOL*         sol,
      SCIP_Real*        pricingobjs,
      int*              varcounter
      )
{
   SCIP* masterprob;

   SCIP_VAR** origvars;
   int norigvars;
   SCIP_VAR* origvar;
   SCIP_Real solval;
   SCIP_Real frac;
   SCIP_Real newobjcoeff;
   int idx;

   int npricingprobs;
   SCIP* pricingprob;
   SCIP_VAR** subvars;
   int nsubvars;
   int nbinvars;
   int nintvars;
   SCIP_SOL* subsol;

   int i;
   int j;
   int v;

   /* get master problem and number of pricing problems */
   masterprob = GCGrelaxGetMasterprob(scip);
   npricingprobs = GCGrelaxGetNPricingprobs(scip);

   for( i = 0; i < npricingprobs; i++ )
   {
      pricingprob = GCGrelaxGetPricingprob(scip, i);
      assert(pricingprob == NULL || GCGrelaxGetNIdenticalBlocks(scip, i) > 0);
      assert(pricingprob != NULL || GCGrelaxGetNIdenticalBlocks(scip, i) == 0);

      /* If the pricing problem is relevant, i. e. not represented by an identical one,
       * solve it */
      if( pricingprob != NULL )
      {
         int nidenticalblocks;
         SCIP_CALL( SCIPgetVarsData(pricingprob, &subvars, &nsubvars, &nbinvars, &nintvars, NULL, NULL) );
         nidenticalblocks = GCGrelaxGetNIdenticalBlocks(scip, i);

         /* The pricing problem may represent a number of other pricing problems
          * (in case of identical blocks); in that case, it has to be solved
          * once for each block */
         for( j = 0; j < nidenticalblocks; j++ )
         {
            /* change objective function values */
            for( v = 0; v < nbinvars + nintvars; v++ )
            {
               assert(GCGvarIsPricing(subvars[v]));
               origvars = GCGpricingVarGetOrigvars(subvars[v]);
               norigvars = GCGpricingVarGetNOrigvars(subvars[v]);
               assert(j < norigvars);

               origvar = origvars[j];
               idx = SCIPvarGetProbindex(origvar);
               solval = SCIPgetSolVal(scip, relaxsol, origvar);
               frac = SCIPfeasFrac(scip, solval);

               /* variables which are already integral, are treated separately */
               if( SCIPisFeasZero(scip, frac) )
               {
                  SCIP_Real lb;
                  SCIP_Real ub;

                  /* variables at their bounds should be kept there */
                  lb = SCIPvarGetLbLocal(origvar);
                  ub = SCIPvarGetUbLocal(origvar);
                  if( SCIPisFeasEQ(scip, solval, lb) )
                     newobjcoeff = 1.0;
                  else if( SCIPisFeasEQ(scip, solval, ub) )
                     newobjcoeff = -1.0;
                  else
                     newobjcoeff = 0.0;
               }
               else
               {
                  if( frac > 0.5 )
                     newobjcoeff = -1.0;
                  else
                     newobjcoeff = 1.0;
               }

               SCIP_CALL( SCIPchgVarObj(pricingprob, subvars[v], newobjcoeff) );
               pricingobjs[idx] = newobjcoeff;
               SCIP_CALL( SCIPsetSolVal(scip, sol, origvar, 0.0) );
            }
            for( v = nbinvars + nintvars; v < nsubvars; v++ )
            {
               assert(GCGvarIsPricing(subvars[v]));
               origvars = GCGpricingVarGetOrigvars(subvars[v]);
               norigvars = GCGpricingVarGetNOrigvars(subvars[v]);
               assert(j < norigvars);

               origvar = origvars[j];
               idx = SCIPvarGetProbindex(origvar);
               SCIP_CALL( SCIPchgVarObj(pricingprob, subvars[v], 0.0) );
               pricingobjs[idx] = newobjcoeff;
               SCIP_CALL( SCIPsetSolVal(scip, sol, origvar, 0.0) );
            }

            /* solve subproblem for current block */
//            GCGpricerSolveSinglePricingProblem(masterprob, i, &solvars, &solvals, &nsolvars, &status);
            SCIP_CALL( SCIPsolve(pricingprob) );
            subsol = SCIPgetBestSol(pricingprob);

            /* set solution values of corresponding block in current working solution */
            for( v = 0; v < nsubvars; v++ )
            {
               assert(GCGvarIsPricing(subvars[v]));
               origvars = GCGpricingVarGetOrigvars(subvars[v]);
               norigvars = GCGpricingVarGetNOrigvars(subvars[v]);
               assert(j < norigvars);

               solval = SCIPgetSolVal(pricingprob, subsol, subvars[v]);

               /* solution values which should be integral may not be integral due to numerics;
                * in that case, round them */
               if( SCIPvarGetType(subvars[v]) != SCIP_VARTYPE_CONTINUOUS )
               {
                  assert(SCIPisEQ(scip, solval, SCIPfloor(scip, solval)));
                  solval = SCIPfloor(scip, solval);
               }

               SCIP_CALL( SCIPsetSolVal(scip, sol, origvars[j], solval) );
            }

//            /* add column to the subSCIP of the master problem */
//            SCIP_CALL( addVarToMastersub(scip, mastersub, i, nsolvars, solvars, solvals, subvars, nsubvars, consmapfw, varcounter) );

            /* free pricing problem s. t. it can be solved again */
            SCIP_CALL( SCIPfreeTransform(pricingprob) );
         }
      }
   }

   return SCIP_OKAY;
}

/** check if there are cycles, i. e. if a solution has already been visited before */
static
SCIP_RETCODE checkCycles(
      SCIP*             scip,
      int               cyclelength,
      int               nloops,
      SCIP_SOL*         sol,
      SCIP_Real         alpha,
      SCIP_SOL**        lastsols,
      SCIP_Real*        lastalphas,
      int*              cycle
      )
{
   SCIP_VAR** vars;
   int nvars;

   SCIP_Real solval1;
   SCIP_Real solval2;

   int i;
   int j;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   // TODO: also take alphas into account

   *cycle = -1;
   for( i = 0; i < MIN(cyclelength, nloops - 1); i++ )
   {
      for( j = 0; j < nvars; j++ )
      {
         solval1 = SCIPgetSolVal(scip, sol, vars[j]);
         solval2 = SCIPgetSolVal(scip, lastsols[i], vars[j]);
         if( !SCIPisFeasEQ(scip, solval1, solval2) )
            break;
      }
      if( j == nvars )
      {
         *cycle = i;
         break;
      }
   }

   return SCIP_OKAY;
}

/** shift a solution in case of a 1-cycle */
// TODO: how to calculate scores; in particular, a reasonable weighting between nLocks and nVarshifts
static
SCIP_RETCODE shiftSol(
      SCIP*             scip,
      SCIP_SOL*         sol,
      SCIP_Real         shiftrate,
      SCIP_Real*        pricingobjs,
      int*              nshifts
      )
{
   SCIP_VAR** vars;
   int nbinvars;
   int nintvars;
   int nvars;

   int maxshifts;
   int* varshifts;
   int score;
   int minscore;
   SCIP_VAR* var;
   SCIP_VAR* pricingvar;
   SCIP_VAR* shiftvar;
   int idx;
   SCIP_Bool increase;

   int i;
   int j;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   SCIP_CALL( SCIPallocBufferArray(scip, &varshifts, nvars) );
   for( i = 0; i < nvars; i++ )
      varshifts[i] = 0;

   *nshifts = 0;

   maxshifts = shiftrate * nvars;
   minscore = INT_MAX;
   shiftvar = NULL;
   for( i = 0; i < maxshifts; i++ )
   {
      for( j = 0; j < nvars; j++ )
      {
         if( pricingobjs[j] == 0 || varshifts[j] == INT_MAX)
            continue;

         var = vars[j];
         assert(GCGvarIsOriginal(var));
         pricingvar = GCGoriginalVarGetPricingVar(var);

         score = pricingobjs[j] > 0 ? SCIPvarGetNLocksUp(pricingvar) + varshifts[j] : SCIPvarGetNLocksDown(pricingvar) + varshifts[j];
         if( score < minscore )
         {
            minscore = score;
            shiftvar = var;
            idx = j;
            increase = pricingobjs[j] > 0;
         }
      }

      if( shiftvar == NULL )
         return SCIP_OKAY;

      if( increase )
         SCIP_CALL( SCIPincSolVal(scip, sol, shiftvar, 1.0) );
      else
         SCIP_CALL( SCIPincSolVal(scip, sol, shiftvar, -1.0) );
      varshifts[idx]++;

//      SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, TRUE, &feasible) );
//      if( !feasible )
//      {
//         if( increase )
//            SCIP_CALL( SCIPincSolVal(scip, sol, shiftvar, -1.0) );
//         else
//            SCIP_CALL( SCIPincSolVal(scip, sol, shiftvar, 1.0) );
//         varshifts[j] = INT_MAX;
//      }
//      else
         (*nshifts)++;
   }

   SCIPfreeBufferArray(scip, &varshifts);

   return SCIP_OKAY;
}




/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#define heurCopyColgenfeaspump NULL

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeColgenfeaspump)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
SCIP_DECL_HEURINIT(heurInitColgenfeaspump)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of colgenfeaspump primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitColgenfeaspump NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitColgenfeaspump)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of colgenfeaspump primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitColgenfeaspump NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolColgenfeaspump)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of colgenfeaspump primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolColgenfeaspump NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolColgenfeaspump)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of colgenfeaspump primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolColgenfeaspump NULL
#endif


/** calculates an adjusted maximal number of LP iterations */
static
SCIP_Longint adjustedMaxNLPIterations(
   SCIP_Longint          maxnlpiterations,   /**< regular maximal number of LP iterations */
   SCIP_Longint          nsolsfound,         /**< total number of solutions found so far by SCIP */
   int                   nstallloops         /**< current number of stalling rounds */
   )
{
   if( nstallloops <= 1 )
   {
      if( nsolsfound == 0 )
         return 4*maxnlpiterations;
      else
         return 2*maxnlpiterations;
   }
   else
      return maxnlpiterations;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecColgenfeaspump)
{  /*lint --e{715}*/
   SCIP* masterprob;
   SCIP_HEURDATA* heurdata;
   SCIP_LPSOLSTAT lpsolstat;  /* status of the LP solution */
   SCIP_RETCODE retcode;
   SCIP_SOL** lastsols;
   SCIP_SOL* relaxsol;
   SCIP_SOL* sol;
   SCIP_SOL* tmpsol;
   SCIP_VAR** mastervars;
   SCIP_VAR** vars;
   SCIP_Bool lperror;
   SCIP_Bool success;
   SCIP_Real alpha;
   SCIP_Real objnorm;         /* Euclidean norm of the objective function, used for scaling */
   SCIP_Real scalingfactor;   /* factor to scale the original objective function with */
   SCIP_Real objfactor;
   SCIP_Real* lastalphas;
   SCIP_Real* pricingobjs;
   SCIP_Real* solvals;
   SCIP_Longint nlpiterations;    /* number of LP iterations done during one pumping round */
   SCIP_Longint maxnlpiterations; /* maximum number of LP iterations fpr this heuristic */
   SCIP_Longint nsolsfound;       /* number of solutions found by this heuristic */
   SCIP_Longint ncalls;           /* number of calls of this heuristic */
   int bestnfracs;
   int cycle;
   int maxloops;
   int maxstallloops;
   int nfracs;
   int nloops;
   int nmastervars;
   int npricingprobs;
   int nshifts;
   int nstallloops;
   int nvars;
   int varcounter;

   int i;
   int j;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get master problem */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   /* get number of pricing problems */
   npricingprobs = GCGrelaxGetNPricingprobs(scip);

   /* get original variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* get master variables' data */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetStage(masterprob) > SCIP_STAGE_SOLVING || SCIPgetLPSolstat(masterprob) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("Not executing CG Feaspump: master LP not solved to optimality.\n");
      return SCIP_OKAY;
   }

   assert(SCIPhasCurrentNodeLP(masterprob));

   /* only call heuristic, if the LP solution is basic (which allows fast resolve in diving) */
   if( !SCIPisLPSolBasic(masterprob) )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(masterprob) == SCIPgetNNodes(masterprob) && SCIPgetDepth(masterprob) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* only call the column generation feasibility pump once at the root */
   if( SCIPgetDepth(scip) == 0 && SCIPheurGetNCalls(heur) > 0 )
      return SCIP_OKAY;

   /* @todo for some reason, the heuristic is sometimes called with an invalid relaxation solution;
    *       in that case, don't execute it */
   if( !SCIPisRelaxSolValid(scip) )
   {
      SCIPdebugMessage("not executing colgen feaspump: invalid relaxation solution (should not happen!)\n");
      return SCIP_OKAY;
   }

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNLPIterations(scip) + SCIPgetNLPIterations(masterprob);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = 10*SCIPheurGetNBestSolsFound(heur) + heurdata->nsuccess;
   maxnlpiterations = (SCIP_Longint)((1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * nlpiterations);
   maxnlpiterations += heurdata->maxlpiterofs;

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;

   /* at the first root call, allow more iterations if there is no feasible solution yet */
   if( SCIPheurGetNCalls(heur) == 0 && SCIPgetNSolsFound(scip) == 0 && SCIPgetDepth(scip) == 0 )
      maxnlpiterations += nlpiterations;

   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + MINLPITER);

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("executing Column Generation Feasibility Pump ...\n");

   /* calculate factor by which alpha is decreased */
   if( heurdata->objfactor == 1.0 )
      objfactor = MIN(1.0 - 0.1 / (SCIP_Real)(1 + SCIPgetNBestSolsFound(scip)), 0.999);
   else
      objfactor = heurdata->objfactor;

   /* calculate maximal number of loops */
   maxloops = (heurdata->maxloops == -1 ? INT_MAX : heurdata->maxloops);
   maxstallloops = (heurdata->maxstallloops == -1 ? INT_MAX : heurdata->maxstallloops);

   /* allocate further memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &pricingobjs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );

   /* allocate memory for cycle handling */
   SCIP_CALL( SCIPallocBufferArray(scip, &lastsols, heurdata->cyclelength) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lastalphas, heurdata->cyclelength) );
   for( i = 0; i < heurdata->cyclelength; i++ )
   {
      SCIP_CALL( SCIPcreateSol(scip, &lastsols[i], heur) );
   }

   /* initialize working solutions */
   SCIP_CALL( SCIPcreateSol(scip, &relaxsol, heur) );
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   /* start diving */
   SCIP_CALL( SCIPstartDive(masterprob) );

   /* lp was solved optimal */
   lperror = FALSE;
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;

   /* change objectives in the master problem */
   /* scale distance function and original objective to the same norm */
   objnorm = SCIPgetObjNorm(scip);
   objnorm = MAX(objnorm, 1.0);
   scalingfactor = SQRT((SCIP_Real)(SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip))) / objnorm;

   nfracs = SCIPgetNExternBranchCands(scip);
   bestnfracs = nfracs;
   nloops = 0;
   nstallloops = 0;
   alpha = 1.0;
   cycle = -1;
   varcounter = 0;

   while( nfracs > 0
      && heurdata->nlpiterations < adjustedMaxNLPIterations(maxnlpiterations, nsolsfound, nstallloops)
      && nloops < maxloops && nstallloops < maxstallloops
      && !SCIPisStopped(scip) )
   {
      SCIP_Longint nlpiterationsleft;
      int iterlimit;

      nloops++;
      alpha *= objfactor;

      SCIPdebugMessage("CG Feasibility Pump loop %d: %d fractional variables (alpha: %.4f, stall: %d/%d)\n",
         nloops, nfracs, alpha, nstallloops, maxstallloops);

      /* create solution from diving LP and try to round it */
      SCIP_CALL( SCIPlinkRelaxSol(scip, relaxsol) );
      SCIP_CALL( SCIProundSol(scip, relaxsol, &success) );

      /* if the rounded solution is feasible and better, add it to SCIP */
      if( success )
      {
         SCIPdebugMessage("colgen feaspump found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, relaxsol));
         //         SCIP_CALL( SCIPtrySol(scip, relaxsol, FALSE, FALSE, FALSE, FALSE, &success) );
#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPtrySol(scip, relaxsol, TRUE, TRUE, TRUE, TRUE, &success) );
#else
         SCIP_CALL( SCIPtrySol(scip, relaxsol, FALSE, TRUE, TRUE, TRUE, &success) );
#endif

         if( success )
         {
            SCIPdebugMessage(" -> solution was feasible and good enough\n");
            *result = SCIP_FOUNDSOL;
         }
      }

      SCIP_CALL( SCIPlinkRelaxSol(scip, relaxsol) );

      /* solve all pricing problems and store the result in the current working solution */
      SCIPdebugMessage(" -> solving pricing problem\n");
      SCIP_CALL( solvePricingProblems(scip, alpha, relaxsol, sol, pricingobjs, &varcounter) );
      SCIPdebugMessage(" -> new solution, obj=%g\n", SCIPgetSolOrigObj(scip, sol));

      /* check for cycles */
      SCIP_CALL( checkCycles(scip, heurdata->cyclelength, nloops, sol, alpha, lastsols, lastalphas, &cycle) );

      if( cycle >= 0 )
      {
         SCIPdebugMessage(" -> cycle of length %d detected, shift variables\n", cycle + 1);
         SCIP_CALL( shiftSol(scip, sol, heurdata->shiftrate, pricingobjs, &nshifts) );
         if( nshifts > 0 )
         {
            SCIPdebugMessage(" -> %d shiftings performed\n", nshifts);
         }
         else
         {
            SCIPdebugMessage(" -> no shifting performed - change alpha\n");
            alpha /= objfactor * objfactor;
            alpha = MIN(alpha, 1.0);
            nstallloops++;
            continue;
         }
      }
//      else if( cycle > 0)
//      {
//         SCIPdebugMessage(" -> cycle of length %d detected, change alpha\n", cycle + 1);
//         alpha /= objfactor * objfactor;
//         alpha = MIN(alpha, 1.0);
//         nstallloops++;
//         continue;
//      }

      /* try to add obtained pricing solution to the solution pool; if it is feasible, then stop */
      SCIP_CALL( SCIPtrySol(scip, sol, FALSE, TRUE, FALSE, TRUE, &success) );
      if( success )
      {
         SCIPdebugMessage(" -> solving pricing problem yielded feasible solution.\n");
         *result = SCIP_FOUNDSOL;
         break;
      }

      /* change objective coefficients in master problem */
      for( i = 0; i < nmastervars; i++ )
      {
         SCIP_VAR* mastervar;
         SCIP_VAR** origvars;
         SCIP_Real* origvals;
         int norigvars;

         SCIP_VAR* origvar;
         SCIP_VAR* pricingvar;
         // SCIP_VARTYPE vartype;
         SCIP_Real masterobjcoeff;
         SCIP_Real newobjcoeff;
         SCIP_Real origval;
         SCIP_Real relaxval;
         SCIP_Real solval;

         mastervar = mastervars[i];
         masterobjcoeff = 0.0;
         origvars = GCGmasterVarGetOrigvars(mastervar);
         origvals = GCGmasterVarGetOrigvals(mastervar);
         norigvars = GCGmasterVarGetNOrigvars(mastervar);

         /* for each original variable contained in mastervar, compute newobjcoeff
          * and add origval * newobjcoeff to objective coefficient of mastervar */
         for( j = 0; j < norigvars; j++ )
         {
            origvar = origvars[j];
            origval = origvals[j];
            assert(GCGvarIsOriginal(origvar));

            pricingvar = GCGoriginalVarGetPricingVar(origvar);
            relaxval = SCIPgetSolVal(scip, relaxsol, origvar);
            solval = SCIPgetSolVal(scip, sol, origvar);
            //            vartype = SCIPvarGetType(origvar);

            /* compute new objective coefficient for original variable */
            if( SCIPisGT(scip, solval, relaxval) )
               newobjcoeff = -alpha + (1.0 - alpha) * (SCIPvarGetNLocksDown(pricingvar) + 1) / (SCIP_Real) (SCIPvarGetNLocksDown(pricingvar) + SCIPvarGetNLocksUp(pricingvar) + 1);
            else
               newobjcoeff = alpha + (1.0 - alpha) * (SCIPvarGetNLocksUp(pricingvar) + 1) / (SCIP_Real) (SCIPvarGetNLocksDown(pricingvar) + SCIPvarGetNLocksUp(pricingvar) + 1);

            /* transfer objective coeff to master variable */
            masterobjcoeff += origval * newobjcoeff;
         }

         SCIP_CALL( SCIPchgVarObjDive(masterprob, mastervar, masterobjcoeff) );
      }

      /* the LP with the new (distance) objective is solved */
      nlpiterations = SCIPgetNLPIterations(scip) + SCIPgetNLPIterations(masterprob);
      nlpiterationsleft = adjustedMaxNLPIterations(maxnlpiterations, nsolsfound, nstallloops) - heurdata->nlpiterations;
      iterlimit = MAX((int)nlpiterationsleft, MINLPITER);
      SCIPdebugMessage(" -> solve LP with iteration limit %d\n", iterlimit);

//      SCIP_CALL( setMasterDivingObjectives(scip, heurdata) );
      retcode = performDivingOnMaster(scip, iterlimit, &nlpiterations, &lperror);
      lpsolstat = SCIPgetLPSolstat(masterprob);

      /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * Hence in optimized mode, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
       */
      if( retcode != SCIP_OKAY )
      {
#ifndef NDEBUG
         if( lpsolstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
         {
            SCIP_CALL( retcode );
         }
#endif
         SCIPwarningMessage("Error while solving LP in Colgen Feaspump heuristic; LP solve terminated with code <%d>\n", retcode);
         SCIPwarningMessage("This does not affect the remaining solution procedure --> continue\n");
      }

      /* update iteration count */
      heurdata->nlpiterations += nlpiterations;
      SCIPdebugMessage(" -> number of iterations: %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", lperror=%u, lpsolstat=%d\n",
            heurdata->nlpiterations, adjustedMaxNLPIterations(maxnlpiterations, nsolsfound, nstallloops), lperror, lpsolstat);

      /* check whether LP was solved optimal */
      if( lperror || lpsolstat != SCIP_LPSOLSTAT_OPTIMAL )
         break;

      /* swap the last solutions, i.e. store the pricing solution into the lastsols array */
      tmpsol = lastsols[heurdata->cyclelength-1];
      for( i = heurdata->cyclelength-1; i > 0; i-- )
      {
         lastsols[i] = lastsols[i-1];
         lastalphas[i] = lastalphas[i-1];
      }
      lastsols[0] = sol;
      lastalphas[0] = alpha;
      sol = tmpsol;

      /* check for improvement in number of fractionals */
      nfracs = SCIPgetNExternBranchCands(scip);
      if( nfracs < bestnfracs )
      {
         bestnfracs = nfracs;
         nstallloops = 0;
      }
      else
         nstallloops++;

      SCIPdebugMessage(" -> loop finished: %d fractional variables (stall: %d/%d, iterations: %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT")\n",
            nfracs, nstallloops, maxstallloops, heurdata->nlpiterations, adjustedMaxNLPIterations(maxnlpiterations, nsolsfound, nstallloops));
   }

   /* try final solution, if no more fractional variables are left */
   if( nfracs == 0 && !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      success = FALSE;

      SCIP_CALL( SCIPlinkRelaxSol(scip, relaxsol) );
      SCIPdebugMessage("colgen feaspump found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, relaxsol));

//      SCIP_CALL( SCIPtrySol(scip, relaxsol, FALSE, FALSE, FALSE, FALSE, &success) );
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPtrySol(scip, relaxsol, TRUE, TRUE, TRUE, TRUE, &success) );
#else
      SCIP_CALL( SCIPtrySol(scip, relaxsol, FALSE, TRUE, TRUE, TRUE, &success) );
#endif

      if( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* end diving */
   SCIP_CALL( SCIPendDive(masterprob) );

   /* free memory */
   SCIPfreeBufferArray(scip, &pricingobjs);
   SCIPfreeBufferArray(scip, &solvals);

   /* free working solutions */
   SCIPfreeSol(scip, &relaxsol);
   SCIPfreeSol(scip, &sol);

   /* free memory for cycle handling */
   for( i = 0; i < heurdata->cyclelength; i++ )
   {
      SCIPfreeSol(scip, &lastsols[i]);
   }
   SCIPfreeBufferArray(scip, &lastsols);
   SCIPfreeBufferArray(scip, &lastalphas);

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the colgenfeaspump primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurColgenfeaspump(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create colgenfeaspump primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyColgenfeaspump,
         heurFreeColgenfeaspump, heurInitColgenfeaspump, heurExitColgenfeaspump,
         heurInitsolColgenfeaspump, heurExitsolColgenfeaspump, heurExecColgenfeaspump,
         heurdata) );

   /* add colgenfeaspump primal heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/"HEUR_NAME"/maxlpiterquot",
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, FALSE, DEFAULT_MAXLPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/"HEUR_NAME"/maxlpiterofs",
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, FALSE, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/cyclelength",
         "maximum length of cycles to be checked explicitly in each round",
         &heurdata->cyclelength, TRUE, DEFAULT_CYCLELENGTH, 1, 100, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxloops",
         "maximal number of pumping rounds (-1: no limit)",
         &heurdata->maxloops, TRUE, DEFAULT_MAXLOOPS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxstallloops",
         "maximal number of pumping rounds without fractionality improvement (-1: no limit)",
         &heurdata->maxstallloops, TRUE, DEFAULT_MAXSTALLLOOPS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/objfactor",
         "factor by which the regard of the objective is decreased in each round",
         &heurdata->objfactor, FALSE, DEFAULT_OBJFACTOR, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/shiftrate",
         "percentage of variables to be shifted in case of a 1-cycle",
         &heurdata->shiftrate, TRUE, DEFAULT_SHIFTRATE, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
