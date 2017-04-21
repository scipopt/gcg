/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   benders_gcg.c
 * @brief  GCG Benders' decomposition algorithm
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "benders_gcg.h"

#include "gcg.h"

#include "relax_gcg.h"
#include "struct_solver.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"
#include "cons_masterbranch.h"

#include "scip/cons_linear.h"
#include "scip/pub_var.h"
#include "scip/pub_benders.h"
#include "scip/bendersdefcuts.h"

#define BENDERS_NAME            "gcg"
#define BENDERS_DESC            "Benders' decomposition template"
#define BENDERS_PRIORITY        0
#define BENDERS_CUTLP        TRUE            /**< should Benders' cut be generated for LP solutions */
#define BENDERS_CUTPSEUDO    TRUE            /**< should Benders' cut be generated for pseudo solutions */
#define BENDERS_CUTRELAX     TRUE            /**< should Benders' cut be generated for relaxation solutions */


#define DEFAULT_USEHEURSOLVING           FALSE      /**< should heuristic solving be used */
#define DEFAULT_ABORTSOLVINGBOUND        FALSE      /**< should the subproblem solve be aborted when exceeding the upper bound */
#define DEFAULT_DISPINFOS                FALSE      /**< should additional information be displayed */
#define DEFAULT_DISABLECUTOFF            2          /**< should the cutoffbound be applied in master LP solving? (0: on, 1:off, 2:auto) */
#define DEFAULT_SORTING                  2          /**< default sorting method for pricing mips
                                                     *    0 :   order of pricing problems
                                                     *    1 :   according to dual solution of convexity constraint
                                                     *    2 :   according to reliability from previous round)
                                                     */
#define DEFAULT_THREADS                  0          /**< number of threads (0 is OpenMP default) */
#define DEFAULT_EAGERFREQ                10         /**< frequency at which all pricingproblems should be solved (0 to disable) */


#define SUBPROBLEM_STAT_ARRAYLEN_TIME 1024   /**< length of the array for Time histogram representation */
#define SUBPROBLEM_STAT_BUCKETSIZE_TIME 10   /**< size of the buckets for Time histogram representation */
#define SUBPROBLEM_STAT_ARRAYLEN_CUTS 1024   /**< length of the array for foundVars histogram representation */
#define SUBPROBLEM_STAT_BUCKETSIZE_CUTS 1    /**< size of the buckets for foundVars histogram representation */


/*
 * Data structures
 */

/* TODO: fill in the necessary Benders' decomposition data */

/** Benders' decomposition data */
struct SCIP_BendersData
{
   SCIP*                 origprob;           /**< the SCIP instance of the original problem */
   SCIP_Real*            subprobobjvals;     /**< the objective values of each pricing problem in the current iteration */
   SCIP_Longint          currnodenr;         /**< current node number in the masterproblem*/
   SCIP_HASHMAP*         mapcons2idx;        /**< hashmap mapping constraints to their index in the conss array */
   SCIP_Real*            score;              /**< score of the pricing problem problems */

   /** variables used for statistics */
   SCIP_CLOCK*           freeclock;          /**< time for freeing pricing problems */
   SCIP_CLOCK*           transformclock;     /**< time for transforming pricing problems */
   int                   solvedsubmipsoptimal; /**< number of optimal pricing runs */
   int                   solvedsubmipsheur;  /**< number of heuristical pricing runs*/
   int                   calls;              /**< number of total pricing calls */
   SCIP_Longint          pricingiters;       /**< sum of all pricing simplex iterations */

   /* solver data */
   GCG_SOLVER**          solvers;            /**< pricing solvers array */
   int                   nsolvers;           /**< number of pricing solvers */

   /* event handler */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler */

   /** parameter values */
   int                   sorting;            /**< how should the subproblems be sorted */
   SCIP_Bool             useheursolving;     /**< should heuristic solving be used? */
   SCIP_Bool             abortsolvebound;    /**< should the subproblem solve be aborted upon when exceeding the current upper bound? */
   SCIP_Bool             dispinfos;          /**< should subproblem solving information be displayed? */
   int                   disablecutoff;      /**< should the cutoffbound be applied in master LP solving (0: on, 1:off, 2:auto)? */
   int                   eagerfreq;          /**< frequency at which all pricingproblems should be solved */
   int                   threads;            /**< the number of threads used to solve the subproblems */

   /** statistics */
   int                   subproblemcalls;    /**< the number of calls to the subproblem */
   int*                  subproblemcallsdist;/**< Calls of each subproblem */
   double*               subproblemtimedist; /**< Time spent solving each subproblem */

   int*                  nodetimehist;       /**< Histogram of nodetime distribution */

   int                   eagerage;           /**< iterations since last eager iteration */
};




/*
 * Local methods
 */

/* this function performs the necessary operations at the first call to the constraint handler. This is to get around
 * the SCIP stages. See scip_farkas and scip_redcost in pricer_gcg.cpp for another example. */
static
SCIP_RETCODE bendersCallOperations(
   SCIP*                 masterprob,
   SCIP_BENDERS*         benders
   )
{
   assert(masterprob != NULL);
   assert(benders != NULL);

   /* this is the same trick that is used in pricer_gcg.cpp. If both DW and BD are used for the same problem, this may
    * need to be changed. */
   if( SCIPbendersGetNCalls(benders) == 0 )
   {
      SCIP_CALL( GCGconsMasterbranchAddRootCons(masterprob) );
   }

   return SCIP_OKAY;
}

/* returns the objective coefficient for the given variable */
static
SCIP_Real varGetObj(
      SCIP_VAR* var
   )
{
   SCIP_VAR* origvar;
   assert(var != NULL);

   origvar = GCGpricingVarGetOrigvars(var)[0];

   if( GCGoriginalVarIsLinking(origvar) )
      return 0.0;
   else
      return SCIPvarGetObj(origvar);
}

/* Initialises the objective function for all subproblems. */
static
SCIP_RETCODE setSubproblemObjs(
   SCIP_BENDERS*        benders            /**< the benders' decomposition constraint handler */
   )
{
   SCIP* subproblem;
   SCIP_VAR** probvars;
   int nprobvars;
   int nsubproblems;
   int i;
   int j;

   assert(benders != NULL);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* set objective value of all variables in the pricing problems to 0 (for farkas pricing) /
    * to the original objective of the variable (for redcost pricing)
    */
   for( i = 0; i < nsubproblems; i++ )
   {
      subproblem = SCIPbendersSubproblem(benders, i);
      if( subproblem == NULL )
         continue;

      probvars = SCIPgetVars(subproblem);
      nprobvars = SCIPgetNVars(subproblem);

      for( j = 0; j < nprobvars; j++ )
      {
         assert(GCGvarGetBlock(probvars[j]) == i);
         assert(GCGoriginalVarIsLinking(GCGpricingVarGetOrigvars(probvars[j])[0]) || (GCGvarGetBlock(GCGpricingVarGetOrigvars(probvars[j])[0]) == i));

         SCIP_CALL( SCIPchgVarObj(subproblem, probvars[j], varGetObj(probvars[j])));

         SCIPdebugMessage("pricingobj var <%s> %f\n", SCIPvarGetName(probvars[j]), varGetObj(probvars[j]));
      }
   }

   return SCIP_OKAY;
}

/* this function fixes the linking variables to the value from the master problem */
static
SCIP_RETCODE setupSubproblem(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   SCIP_SOL*             sol,                /**< the current incumbent solution to check. Can be NULL is generating cuts for LP solutions */
   int                   probnumber          /**< the subproblem number that is beting set up */
   )
{
   SCIP_VAR** mastervars;
   SCIP_VAR** masterfixedvars;
   SCIP_VAR** pricingvars;       /* the pricing vars that are related to the linking vars */
   SCIP_Real fixedval;
   int nmastervars;
   int nmasterfixedvars;
   int i;


   assert(masterprob != NULL);

   mastervars = SCIPgetVars(masterprob);
   masterfixedvars = SCIPgetFixedVars(masterprob);
   nmastervars = SCIPgetNVars(masterprob);
   nmasterfixedvars = SCIPgetNFixedVars(masterprob);


   for( i = 0; i < nmastervars + nmasterfixedvars; i++ )
   {
      SCIP_VAR* var;

      if( i < nmastervars )
         var = mastervars[i];
      else
         var = masterfixedvars[i - nmastervars];


      /* check whether the master variable is a linking variable
       * In the case that the variable is a linking variable, then all corresponding pricing variables will be fixed */
      if( GCGmasterVarIsLinking(var) )
      {
         /* storing the value of the linking variable from the given solution.
          * If the current solution is NULL, then the LP or Pseudo solution will be returned */
         fixedval = SCIPgetSolVal(masterprob, sol, var);

         /* collecting all pricing variables that are associated with this linking variable */
         pricingvars = GCGlinkingVarGetPricingVars(GCGmasterVarGetOrigvars(var)[0]);

         if( pricingvars[probnumber] != NULL )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;

            SCIP_CALL( SCIPfixVar(SCIPbendersSubproblem(benders, probnumber), pricingvars[probnumber], fixedval,
                  &infeasible, &fixed) );

            assert(!infeasible);
            assert(fixed);
         }
      }
   }

   return SCIP_OKAY;
}

/* solving the subproblems to generated Benders' cuts */
static
SCIP_RETCODE solveSubproblem(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   SCIP_SOL*             sol,                /**< the current solution, can be null if an lp of pseudo solution is checked */
   SCIP_Real*            objval,             /**< the sum of the objective function values over all subproblems */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible          /**< a flag to indicate whether all subproblems are feasible */
   )
{
   SCIP* subproblem;
   SCIP_SOL* bestsol;

   /* previous parameter settings */
   int prevCutoffParam;
   int prevPropMaxroundsParam;
   int prevPropMaxroundsRootParam;
   char prevInitAlgParam;
   char prevResolveAlgParam;
   SCIP_Bool prevConfParam;
   SCIP_Bool prevDualParam;

   assert(masterprob != NULL);
   assert(benders != NULL);
   assert(objval != NULL);
   assert(infeasible != NULL);

   (*objval) = 0;
   (*infeasible) = FALSE;


   /* TODO: This should be solved just as an LP, so as a MIP. There is too much overhead with the MIP.
    * Need to change status check for checking the LP. */
   subproblem = SCIPbendersSubproblem(benders, probnumber);

   /* modifying all of the parameters */

   /* Do we have to disable presolving? If yes, we have to store all presolving parameters. */
   SCIPsetPresolving(subproblem, SCIP_PARAMSETTING_OFF, TRUE);

   /* Disabling heuristics so that the problem is not trivially solved */
   SCIPsetHeuristics(subproblem, SCIP_PARAMSETTING_OFF, TRUE);

   /* store parameters that are changed for the generation of the subproblem cuts */
   SCIP_CALL( SCIPgetBoolParam(subproblem, "conflict/enable", &prevConfParam) );
   SCIPsetParam(subproblem, "conflict/enable", FALSE);

   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/disablecutoff", &prevCutoffParam) );
   SCIPsetIntParam(subproblem, "lp/disablecutoff", 1);

   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/initalgorithm", &prevInitAlgParam) );
   SCIPsetCharParam(subproblem, "lp/initalgorithm", 'd');
   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/resolvealgorithm", &prevResolveAlgParam) );
   SCIPsetCharParam(subproblem, "lp/resolvealgorithm", 'd');

   SCIP_CALL( SCIPgetBoolParam(subproblem, "misc/alwaysgetduals", &prevDualParam) );
   SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", TRUE);

   //SCIPinfoMessage(subproblem, NULL, "Pricing problem %d\n", i);
   //SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NORMAL) );
   //SCIP_CALL( SCIPsetBoolParam(subproblem, "display/lpinfo", TRUE) );

   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxrounds", &prevPropMaxroundsParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxrounds", 0) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxroundsroot", &prevPropMaxroundsRootParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxroundsroot", 0) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "constraints/linear/propfreq", -1) );

   SCIP_CALL( SCIPsolve(subproblem) );

   bestsol = SCIPgetBestSol(subproblem);

   if( SCIPgetStatus(subproblem) == SCIP_STATUS_OPTIMAL )
   {
      (*objval) += SCIPgetSolTransObj(subproblem, bestsol);
   }
   else if( SCIPgetStatus(subproblem) == SCIP_STATUS_INFEASIBLE )
   {
      (*infeasible) = TRUE;
      (*objval) = SCIPinfinity(masterprob);
   }
   else
      assert(FALSE);

   //SCIP_CALL( SCIPprintStatistics(subprob, NULL) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
   /* resetting the parameter settings to the previous state */
   SCIPsetPresolving(subproblem, SCIP_PARAMSETTING_DEFAULT, TRUE);
   SCIPsetHeuristics(subproblem, SCIP_PARAMSETTING_DEFAULT, TRUE);
   SCIPsetBoolParam(subproblem, "conflict/enable", prevConfParam);
   SCIPsetIntParam(subproblem, "lp/disablecutoff", prevCutoffParam);
   SCIPsetCharParam(subproblem, "lp/initalgorithm", prevInitAlgParam);
   SCIPsetCharParam(subproblem, "lp/resolvealgorithm", prevResolveAlgParam);
   SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", prevDualParam);
   SCIPsetIntParam(subproblem, "propagating/maxrounds", prevPropMaxroundsParam);
   SCIPsetIntParam(subproblem, "propagating/maxroundsroot", prevPropMaxroundsRootParam);

   return SCIP_OKAY;
}

/** sets the values the given variables in the original problem */
static
SCIP_RETCODE setOriginalProblemValues(
   SCIP*                 origprob,           /**< the SCIP instance of the original problem */
   SCIP_SOL*             origsol,            /**< the solution for the original problem */
   SCIP_VAR**            vars,               /**< the variables from the decomposed problem */
   SCIP_Real*            vals,               /**< the solution values of the given problem */
   int                   nvars,              /**< the number of variables */
   SCIP_Bool             master              /**< are the variables from the master problem */
   )
{
   SCIP_VAR** origvars;
   SCIP_Real* origvals;
   int norigvars;
   int i;

   assert(origprob != NULL);
   assert(vars != NULL);
   assert(vals != NULL);

   /* looping through all variables to update the values in the original solution */
   for( i = 0; i < nvars; i++ )
   {
      norigvars = master ? GCGmasterVarGetNOrigvars(vars[i]) : GCGpricingVarGetNOrigvars(vars[i]);
      if( norigvars > 0 )
      {
         origvars = master ? GCGmasterVarGetOrigvars(vars[i]) : GCGpricingVarGetOrigvars(vars[i]);

         if( master )
            origvals = GCGmasterVarGetOrigvals(vars[i]);

         /* all master variables should be associated with a single original variable. This is because no reformulation has
          * been performed. */
         assert(norigvars == 1);
         assert(!master || origvals[0] == 1.0);

         assert((master && GCGvarIsMaster(vars[i])) || (!master && GCGvarIsPricing(vars[i])));

         assert(!SCIPisInfinity(origprob, vals[i]));

         SCIP_CALL( SCIPsetSolVal(origprob, origsol, origvars[0], vals[i]) );
      }
   }

   return SCIP_OKAY;
}



/** creates an original problem solution from the master and subproblem solutions */
static
SCIP_RETCODE createOriginalProblemSolution(
   SCIP*                 masterprob,
   SCIP_BENDERS*         benders
   )
{
   /* Don't know whether we need to consider the fixed variables. Must check. */
   SCIP* origprob;
   SCIP* subproblem;
   SCIP_BENDERSDATA* bendersdata;
   SCIP_SOL* origsol;
   SCIP_SOL* bestsol;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nsubproblems;
   int nvars;
   int i;
   SCIP_Bool stored;

   assert(masterprob != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);

   assert(bendersdata != NULL);

   origprob = bendersdata->origprob;

   /* creating the original problem */
   SCIP_CALL( SCIPcreateSol(origprob, &origsol, GCGrelaxGetProbingheur(origprob)) );

   /* setting the values of the master variables in the original solution */

   /* getting the variable data for the master variables */
   SCIP_CALL( SCIPgetVarsData(masterprob, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(vars != NULL);

   /* getting the best solution from the master problem */
   bestsol = SCIPgetBestSol(masterprob);

   SCIP_CALL( SCIPallocBufferArray(masterprob, &vals, nvars) );
   SCIP_CALL( SCIPgetSolVals(masterprob, bestsol, nvars, vars, vals) );

   /* setting the values using the master problem solution */
   SCIP_CALL( setOriginalProblemValues(origprob, origsol, vars, vals, nvars, TRUE) );

   /* freeing the values buffer array for use for the pricing problems */
   SCIPfreeBufferArray(masterprob, &vals);

   /* setting the values of the subproblem variables in the original solution */
   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* looping through all subproblems */
   for( i = 0; i < nsubproblems; i++ )
   {
      subproblem = SCIPbendersSubproblem(benders, i);

      /* getting the variable data for the master variables */
      SCIP_CALL( SCIPgetVarsData(subproblem, &vars, &nvars, NULL, NULL, NULL, NULL) );
      assert(vars != NULL);

      /* getting the best solution from the master problem */
      bestsol = SCIPgetBestSol(subproblem);

      SCIP_CALL( SCIPallocBufferArray(subproblem, &vals, nvars) );
      SCIP_CALL( SCIPgetSolVals(subproblem, bestsol, nvars, vars, vals) );

      /* setting the values using the master problem solution */
      SCIP_CALL( setOriginalProblemValues(origprob, origsol, vars, vals, nvars, FALSE) );

      /* freeing the values buffer array for use for the pricing problems */
      SCIPfreeBufferArray(subproblem, &vals);
   }

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPtrySol(origprob, origsol, TRUE, TRUE, TRUE, TRUE, TRUE, &stored) );
#else
   SCIP_CALL( SCIPtrySol(origprob, origsol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );
#endif
   if( !stored )
   {

      SCIP_CALL( SCIPcheckSolOrig(origprob, origsol, &stored, TRUE, TRUE) );
   }
   /** @bug The solution doesn't have to be accepted, numerics might bite us, so the transformation might fail.
    *  A remedy could be: Round the values or propagate changes or call a heuristic to fix it.
    */
   SCIP_CALL( SCIPfreeSol(origprob, &origsol) );

   if( stored )
      SCIPdebugMessage("  updated current best primal feasible solution.\n");


   return SCIP_OKAY;
}




/*
 * Callback methods for Benders' decomposition
 */

/* TODO: Implement all necessary Benders' decomposition methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for benders plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BENDERSCOPY(bendersCopyGcg)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of gcg Benders' decompostion not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersCopyGcg NULL
#endif

/** destructor of Benders' decomposition to free user data (called when SCIP is exiting) */
#if 1
static
SCIP_DECL_BENDERSFREE(bendersFreeGcg)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);

   if( bendersdata != NULL )
   {
      SCIPfreeMemory(scip, &bendersdata);
   }

   return SCIP_OKAY;
}
#else
#define bendersFreeGcg NULL
#endif


/** initialization method of Benders' decomposition (called after problem was transformed) */
#if 0
static
SCIP_DECL_BENDERSINIT(bendersInitGcg)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of gcg Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitGcg NULL
#endif


/** deinitialization method of Benders' decomposition (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BENDERSEXIT(bendersExitGcg)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of gcg Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitGcg NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_BENDERSINITPRE(bendersInitpreGcg)
{  /*lint --e{715}*/
   int nsubproblems;
   int i;

   assert(scip != NULL);
   assert(benders != NULL);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   for( i = 0; i < nsubproblems; i++ )
      SCIP_CALL( GCGaddDataAuxiliaryVar(scip, SCIPbendersGetAuxiliaryVar(benders, i), i) );

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_BENDERSEXITPRE(bendersExitpreGcg)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of gcg Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitpreGcg NULL
#endif


/** solving process initialization method of Benders' decomposition (called when branch and bound process is about to begin) */
static
SCIP_DECL_BENDERSINITSOL(bendersInitsolGcg)
{  /*lint --e{715}*/
   SCIP* origprob;
   SCIP_BENDERSDATA* bendersdata;
   SCIP_CONS** masterconss;
   int nmasterconss;
   int nsubproblems;
   int origverblevel;
   int i;

   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);

   /* TODO: I don't think that I need to check for the scip_. This parameter is only provided in classes related to the
    * pricer. */
   /*
   assert(scip == scip_);
   */
   assert(bendersdata != NULL);

   origprob = bendersdata->origprob;

   /* at the beginning, the output of the master problem gets the same verbosity level
    * as the output of the original problem */
   SCIP_CALL( SCIPgetIntParam(origprob, "display/verblevel", &origverblevel) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", origverblevel) );

   bendersdata->currnodenr = -1;
   bendersdata->eagerage = 0;

   nmasterconss = GCGgetNMasterConss(origprob);
   masterconss = GCGgetMasterConss(origprob);

   /* init array containing all pricing problems */
   nsubproblems = SCIPbendersGetNSubproblems(benders);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(bendersdata->subprobobjvals), nsubproblems) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(bendersdata->subproblemcallsdist), nsubproblems) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(bendersdata->subproblemtimedist), nsubproblems) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(bendersdata->nodetimehist), SUBPROBLEM_STAT_ARRAYLEN_TIME) ); /*lint !e506*/

   BMSclearMemoryArray(bendersdata->nodetimehist, SUBPROBLEM_STAT_ARRAYLEN_TIME);

   for( i = 0; i < nsubproblems; i++ )
   {
      bendersdata->subprobobjvals[i] = SCIPinfinity(scip);
      bendersdata->subproblemcallsdist[i] = 0;
      bendersdata->subproblemtimedist[i] = 0;

      /* TODO: Not sure whether the relevance check is necessary for Benders'. Must check */
      if( TRUE || GCGisPricingprobRelevant(origprob, i) )
      {
         SCIPbendersAddSubproblem(benders, GCGgetPricingprob(origprob, i));
      }
#if 0
      else
      {
         bendersdata->pricingprobs[i] = NULL;
      }
#endif
   }

   /* alloc memory for arrays of reduced cost */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(bendersdata->score), nsubproblems) );

   SCIP_CALL( SCIPcreateCPUClock(scip, &(bendersdata->freeclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(bendersdata->transformclock)) );

   bendersdata->solvedsubmipsoptimal = 0;
   bendersdata->solvedsubmipsheur = 0;
   bendersdata->calls = 0;
   bendersdata->pricingiters = 0;

   /* TODO: Get a parameter for the type of cut to be added. Not sure what this means at the moment, but will work this
    * out. */

   SCIP_CALL( SCIPhashmapCreate(&(bendersdata->mapcons2idx), SCIPblkmem(scip), 10 * nmasterconss +1) );
   for( i = 0; i < nmasterconss; i++ )
   {
      SCIP_CALL( SCIPhashmapInsert(bendersdata->mapcons2idx, masterconss[i], (void*)(size_t)i) );
      assert((int)(size_t)SCIPhashmapGetImage(bendersdata->mapcons2idx, masterconss[i]) == i); /*lint !e507*/
   }


   /* setting the objective coefficients for the subproblems.
    * This is required because the variables are added to the pricing problems with a zero coefficient. In the DW
    * context, this is appropriate because the objective coefficients are constantly changing. In the BD context, the
    * objective coefficients are static, so they only need to be updated once. */
   SCIP_CALL( setSubproblemObjs(benders) );

   /* TODO: It is not clear that this function is needed. Only a slight modification is needed for the solver interface,
    * i.e. allowing the return of dual and farkas solutions, so it could work. If the solver interface will be used for
    * Benders' decomposition, the solversInitsol() will need to become a global function. */
   /*
    * SCIP_CALL( solversInitsol() );
    */


   /* TODO: It is not clear whether this function is needed. Will need to check what is done here. */
   /*
    * SCIP_CALL( SCIPactivateEventHdlrDisplay(scip_) );
    */

   return SCIP_OKAY;
}


/** solving process deinitialization method of Benders' decomposition (called before branch and bound process data is freed) */
static
SCIP_DECL_BENDERSEXITSOL(bendersExitsolGcg)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;
   int nsubproblems;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   SCIPhashmapFree(&(bendersdata->mapcons2idx));

   SCIP_CALL( SCIPfreeClock(scip, &(bendersdata->transformclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(bendersdata->freeclock)) );

   SCIPfreeBlockMemoryArray(scip, &(bendersdata->score), nsubproblems);

   SCIPfreeMemoryArray(scip, &(bendersdata->nodetimehist));

   SCIPfreeBlockMemoryArray(scip, &(bendersdata->subproblemtimedist), nsubproblems);
   SCIPfreeBlockMemoryArray(scip, &(bendersdata->subproblemcallsdist), nsubproblems);
   SCIPfreeBlockMemoryArray(scip, &(bendersdata->subprobobjvals), nsubproblems);

   return SCIP_OKAY;
}


/** mapping method between the master problem variables and the subproblem variables of Benders' decomposition */
static
SCIP_DECL_BENDERSGETMASTERVAR(bendersGetmastervarGcg)
{  /*lint --e{715}*/
   SCIP_VAR* origvar;
   SCIP_VAR* mastervar;

   assert(scip != NULL);
   assert(benders != NULL);
   assert(var != NULL);

   /* if there is no corresponding master variable for the input variable, then NULL is returned */
   mastervar = NULL;

   /* getting the original variable for the given pricing variable */
   origvar = GCGpricingVarGetOrigvars(var)[0];

   /* checking whether the original variable is a linking variable.
    * If this is the case, then the corresponding master variable is added to the generated cut.
    * If the pricing variable is not a linking variable, then the farkas dual value is added to the lhs */
   if( GCGoriginalVarIsLinking(origvar) )
   {
      mastervar = GCGoriginalVarGetMastervars(origvar)[0];
      mastervar = SCIPvarGetProbvar(mastervar);
   }


   return mastervar;
}


/** the execution method for Benders' decomposition */
static
SCIP_DECL_BENDERSEXEC(bendersExecGcg)
{  /*lint --e{715}*/

   SCIP_CALL( bendersCallOperations(scip, benders) );

   return SCIP_OKAY;
}


/** the post-execution method for Benders' decomposition */
static
SCIP_DECL_BENDERSPOSTSOLVE(bendersPostsolveGcg)
{  /*lint --e{715}*/

   /* creates a solution to the original problem */
   if( !infeasible )
      SCIP_CALL( createOriginalProblemSolution(scip, benders) );

   return SCIP_OKAY;
}


/** the subproblem solving method for Benders' decomposition */
static
SCIP_DECL_BENDERSSOLVESUB(bendersSolvesubGcg)
{  /*lint --e{715}*/
   SCIP_Real objval;

   assert(scip != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));
   assert(SCIPbendersSubproblem(benders, probnumber) != NULL);

   /* TODO: Need to get the solution. Is this provided or passed through the constraint handler? */

   /* setting up the subproblem to be solved */
   SCIP_CALL( setupSubproblem(scip, benders, sol, probnumber) );

   /* solves the subproblem */
   SCIP_CALL( solveSubproblem(scip, benders, sol, &objval, probnumber, infeasible) );

   return SCIP_OKAY;
}


/** the subproblem freeing method for Benders' decomposition */
static
SCIP_DECL_BENDERSFREESUB(bendersFreesubGcg)
{  /*lint --e{715}*/

   /* freeing the transform of the subproblems so that it can be updated in the next iteration. */
   /* TODO: Consider using reoptimisation. In this case you will need to free the reoptimization problem. */
   SCIP_CALL( SCIPfreeTransform(SCIPbendersSubproblem(benders, probnumber)) );

   return SCIP_OKAY;
}


/*
 * Benders' decomposition specific interface methods
 */

/** creates the gcg Benders' decomposition and includes it in SCIP */
SCIP_RETCODE SCIPincludeBendersGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origprob,           /**< the SCIP instance of the original problem */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   )
{
   SCIP_BENDERSDATA* bendersdata;
   SCIP_BENDERS* benders;

   /* create gcg Benders' decomposition data */
   SCIP_CALL( SCIPallocMemory(scip, &bendersdata) );
   bendersdata->origprob = origprob;
   bendersdata->solvers = NULL;
   bendersdata->nsolvers = 0;
   bendersdata->nodetimehist = NULL;

   benders = NULL;


   /* include Benders' decomposition */
#if 0
   /* use SCIPincludeBenders() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBenders(scip, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, nsubproblems,
         bendersCopyGcg, bendersFreeGcg, bendersInitGcg, bendersExitGcg, bendersInitpreGcg, bendersExitpreGcg,
         bendersInitsolGcg, bendersExitsolGcg, bendersGetmastervarGcg, bendersExecGcg,
         bendersPostsolveGcg, bendersFreesubGcg, bendersdata) );
#else
   /* use SCIPincludeBendersBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBendersBasic(scip, &benders, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, nsubproblems,
         BENDERS_CUTLP, BENDERS_CUTPSEUDO, BENDERS_CUTRELAX, bendersGetmastervarGcg, bendersExecGcg,
         bendersSolvesubGcg, bendersFreesubGcg, bendersdata) );
   assert(benders != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBendersCopy(scip, benders, bendersCopyGcg) );
   SCIP_CALL( SCIPsetBendersFree(scip, benders, bendersFreeGcg) );
   SCIP_CALL( SCIPsetBendersInit(scip, benders, bendersInitGcg) );
   SCIP_CALL( SCIPsetBendersExit(scip, benders, bendersExitGcg) );
   SCIP_CALL( SCIPsetBendersInitpre(scip, benders, bendersInitpreGcg) );
   SCIP_CALL( SCIPsetBendersExitpre(scip, benders, bendersExitpreGcg) );
   SCIP_CALL( SCIPsetBendersInitsol(scip, benders, bendersInitsolGcg) );
   SCIP_CALL( SCIPsetBendersExitsol(scip, benders, bendersExitsolGcg) );
   SCIP_CALL( SCIPsetBendersPostsolve(scip, benders, bendersPostsolveGcg) );
#endif

   /* including the default cuts for Benders' decomposition */
   SCIP_CALL( SCIPincludeBendersDefaultCuts(scip, benders) );

   /* add gcg Benders' decomposition parameters */
   /* parameters for the Benders' decomposition constraint handler */
   SCIP_CALL( SCIPaddBoolParam(origprob, "benders/gcg/useheursolving",
         "should subproblem solving be performed heuristically before solving the LPs to optimality?",
         &bendersdata->useheursolving, TRUE, DEFAULT_USEHEURSOLVING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "benders/gcg/abortsolvingbound",
         "should solving be aborted when the objective function is less than the current upper bound?",
         &bendersdata->abortsolvebound, TRUE, DEFAULT_ABORTSOLVINGBOUND, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "benders/gcg/dispinfos",
         "should additional informations concerning the subproblem solving process be displayed?",
         &bendersdata->dispinfos, FALSE, DEFAULT_DISPINFOS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "benders/gcg/sorting",
         "which sorting method should be used to sort the subproblems problems (0 = order of pricing problems, 1 = according to dual solution of convexity constraint, 2 = according to reliability from previous round)",
         &bendersdata->sorting, FALSE, DEFAULT_SORTING, 0, 5, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "benders/gcg/threads",
         "how many threads should be used to concurrently solve the subprolems (0 to guess threads by OpenMP)",
         &bendersdata->threads, FALSE, DEFAULT_THREADS, 0, 4096, NULL, NULL) );

   //SCIP_CALL( SCIPsetIntParam(scip, "lp/disablecutoff", DEFAULT_DISABLECUTOFF) );
   //SCIP_CALL( SCIPaddIntParam(origprob, "benders/gcg/disablecutoff",
         //"should the cutoffbound be applied in master LP solving (0: on, 1:off, 2:auto)?",
         //&bendersdata->disablecutoff, FALSE, DEFAULT_DISABLECUTOFF, 0, 2, paramChgdDisablecutoff, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "benders/gcg/eagerfreq",
            "frequency at which all subproblems should be solved (0 to disable)",
            &bendersdata->eagerfreq, FALSE, DEFAULT_EAGERFREQ, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
