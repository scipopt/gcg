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
//#define SCIP_DEBUG
#include <assert.h>
#include <string.h>

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

#define BENDERS_NAME                "gcg"
#define BENDERS_DESC                "Benders' decomposition for the Generic Column Generation package"
#define BENDERS_PRIORITY         1000
#define BENDERS_CUTLP            TRUE   /**< should Benders' cut be generated for LP solutions */
#define BENDERS_CUTPSEUDO        TRUE   /**< should Benders' cut be generated for pseudo solutions */
#define BENDERS_CUTRELAX         TRUE   /**< should Benders' cut be generated for relaxation solutions */
#define BENDERS_SHAREAUXVARS    FALSE   /**< should this Benders' share the highest priority Benders' aux vars */

#define SUBPROBLEM_STAT_ARRAYLEN_TIME 1024   /**< length of the array for Time histogram representation */
#define SUBPROBLEM_STAT_BUCKETSIZE_TIME 10   /**< size of the buckets for Time histogram representation */
#define SUBPROBLEM_STAT_ARRAYLEN_CUTS 1024   /**< length of the array for foundVars histogram representation */
#define SUBPROBLEM_STAT_BUCKETSIZE_CUTS 1    /**< size of the buckets for foundVars histogram representation */

#define LARGE_VALUE  10000    /**< a large value that is used to create an artificial solution */

/*
 * Data structures
 */

/* TODO: fill in the necessary Benders' decomposition data */

/** Benders' decomposition data */
struct SCIP_BendersData
{
   SCIP*                 origprob;           /**< the SCIP instance of the original problem */
   SCIP_SOL*             relaxsol;           /**< the solution to the original problem related to the relaxation */
};




/*
 * Local methods
 */

/* returns the objective coefficient for the given variable */
static
SCIP_Real varGetObj(
   SCIP_VAR*             var
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
   SCIP_BENDERS*         benders,            /**< the benders' decomposition constraint handler */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_VAR** probvars;
   int nprobvars;
   int i;

   assert(benders != NULL);

   /* changing the variable */
   subproblem = SCIPbendersSubproblem(benders, probnumber);
   assert(subproblem != NULL);

   probvars = SCIPgetVars(subproblem);
   nprobvars = SCIPgetNVars(subproblem);

   for( i = 0; i < nprobvars; i++ )
   {
      assert(GCGvarGetBlock(probvars[i]) == probnumber);
      assert(GCGoriginalVarIsLinking(GCGpricingVarGetOrigvars(probvars[i])[0]) || (GCGvarGetBlock(GCGpricingVarGetOrigvars(probvars[i])[0]) == probnumber));

      SCIP_CALL( SCIPchgVarObj(subproblem, probvars[i], varGetObj(probvars[i])));

      SCIPdebugMessage("pricingobj var <%s> %f\n", SCIPvarGetName(probvars[i]), varGetObj(probvars[i]));
   }

   return SCIP_OKAY;
}

/** sets the values the given variables in the original problem */
static
SCIP_RETCODE setOriginalProblemValues(
   SCIP*                 origprob,           /**< the SCIP instance of the original problem */
   SCIP*                 masterprob,         /**< the Benders' master problem */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SOL*             origsol,            /**< the solution for the original problem */
   SCIP_VAR**            vars,               /**< the variables from the decomposed problem */
   SCIP_Real*            vals,               /**< the solution values of the given problem, can be NULL */
   int                   nvars,              /**< the number of variables */
   SCIP_Bool             master,             /**< are the variables from the master problem */
   SCIP_Bool             artificial          /**< should an artifical (possibly infeasible) solution be created to
                                                  generate branching candidates */
   )
{
   SCIP_VAR** origvars;
   SCIP_Real val;
   int norigvars;
   int i;

#ifndef NDEBUG
   SCIP_Real* origvals;
#endif

   assert(origprob != NULL);
   assert(masterprob != NULL);
   assert(benders != NULL);
   assert(vars != NULL);
   assert(vals != NULL || artificial);

   /* looping through all variables to update the values in the original solution */
   for( i = 0; i < nvars; i++ )
   {
      norigvars = master ? GCGmasterVarGetNOrigvars(vars[i]) : GCGpricingVarGetNOrigvars(vars[i]);
      if( norigvars > 0 )
      {
         SCIP_VAR* mastervar;

         origvars = master ? GCGmasterVarGetOrigvars(vars[i]) : GCGpricingVarGetOrigvars(vars[i]);

#ifndef NDEBUG
         if( master )
            origvals = GCGmasterVarGetOrigvals(vars[i]);
#endif

         /* all master variables should be associated with a single original variable. This is because no reformulation has
          * been performed. */
         assert(norigvars == 1);
         assert(!master || origvals[0] == 1.0);

         assert((master && GCGvarIsMaster(vars[i])) || (!master && GCGvarIsPricing(vars[i])));

         /* for all variables that are from the subproblems, they are set to their bounds if the solution is being
          * created to identify branching candidates. */
         if( !master && artificial )
         {
            if( SCIPisNegative(origprob, SCIPvarGetObj(origvars[0])) )
            {
               val = SCIPvarGetLbGlobal(origvars[0]);
               if( SCIPisInfinity(origprob, -val) )
                  val = -LARGE_VALUE;
            }
            else
            {
               val = SCIPvarGetUbGlobal(origvars[0]);
               if( SCIPisInfinity(origprob, val) )
                  val = LARGE_VALUE;
            }
         }
         else
            val = vals[i];

         /* identifying whether the variable is a master problem variable. The variable is a master problem variable if
          * there is a mapping from the subproblem to the master problem */
         if( master )
            mastervar = vars[i];
         else
         {
            mastervar = NULL;
            SCIP_CALL( SCIPgetBendersMasterVar(masterprob, benders, vars[i], &mastervar) );
         }

         assert(!SCIPisInfinity(origprob, val));

         SCIPdebugMsg(masterprob, "setting the value of <%s> (dw variable <%s>) to %g in the original solution. "
            "Variable type: %d\n", SCIPvarGetName(origvars[0]), SCIPvarGetName(vars[i]), val, master);

         /* only update the solution value of master variables if an artificial solution is being created. */
         if( master || mastervar == NULL )
            SCIP_CALL( SCIPsetSolVal(origprob, origsol, origvars[0], val) );
      }
   }

   return SCIP_OKAY;
}



/** creates an original problem solution from the master and subproblem solutions */
static
SCIP_RETCODE createOriginalProblemSolution(
   SCIP*                 masterprob,         /**< the Benders' master problem */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SOL*             sol,                /**< the solution passed to the Benders' decomposition subproblems. */
   SCIP_Bool             artificial          /**< should an artifical (possibly infeasible) solution be created to
                                                  generate branching candidates */
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
   SCIP_CALL( SCIPallocBufferArray(masterprob, &vals, nvars) );
   SCIP_CALL( SCIPgetSolVals(masterprob, sol, nvars, vars, vals) );

   /* setting the values using the master problem solution */
   SCIP_CALL( setOriginalProblemValues(origprob, masterprob, benders, origsol, vars, vals, nvars, TRUE, artificial) );

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
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintSol(subproblem, bestsol, NULL, FALSE) );
#endif

      /* the branching candidates come from the master problem solution. However, we need a full solution to pass to the
       * original problem to find the branching candidate. So the subproblem variables are set to their bounds, creating
       * a possibly infeasible solution, but with the fractional master problem variables. */
      if( artificial )
      {
         /* setting the values of the subproblem variables to their bounds. */
         SCIP_CALL( setOriginalProblemValues(origprob, masterprob, benders, origsol, vars, NULL, nvars, FALSE, artificial) );
      }
      else
      {
         SCIP_CALL( SCIPallocBufferArray(subproblem, &vals, nvars) );
         SCIP_CALL( SCIPgetSolVals(subproblem, bestsol, nvars, vars, vals) );

         /* setting the values using the master problem solution */
         SCIP_CALL( setOriginalProblemValues(origprob, masterprob, benders, origsol, vars, vals, nvars, FALSE, artificial) );

         /* freeing the values buffer array for use for the pricing problems */
         SCIPfreeBufferArray(subproblem, &vals);
      }
   }

   /* if the solution is NULL, then the solution comes from the relaxation. Thus, it should be stored in the
    * bendersdata. When it is not null, then solution comes from a heuristic. So this solution should be passed to the
    * solution storage. */
   if( sol != NULL )
   {
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
   }
   else
   {
      if( bendersdata->relaxsol != NULL )
         SCIP_CALL( SCIPfreeSol(origprob, &bendersdata->relaxsol) );

      bendersdata->relaxsol = origsol;
   }

   //printf("Master Solution\n");
   //SCIP_CALL( SCIPprintSol(masterprob, sol, 0, 0) );
   //printf("Original Solution\n");
   //SCIP_CALL( SCIPprintSol(origprob, origsol, 0, 0) );

   return SCIP_OKAY;
}




/*
 * Callback methods for Benders' decomposition
 */

/* TODO: Implement all necessary Benders' decomposition methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of Benders' decomposition to free user data (called when SCIP is exiting) */
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


/** solving process deinitialization method of Benders' decomposition (called before branch and bound process data is freed) */
static
SCIP_DECL_BENDERSEXITSOL(bendersExitsolGcg)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);

   /* freeing the relaxation solution */
   if( bendersdata->relaxsol != NULL )
      SCIP_CALL( SCIPfreeSol(bendersdata->origprob, &bendersdata->relaxsol) );

   return SCIP_OKAY;
}


/** mapping method between the master problem variables and the subproblem variables of Benders' decomposition */
static
SCIP_DECL_BENDERSGETVAR(bendersGetvarGcg)
{  /*lint --e{715}*/
   SCIP_VAR* origvar;

   assert(scip != NULL);
   assert(benders != NULL);
   assert(var != NULL);
   assert(mappedvar != NULL);

   /* if there is no corresponding master variable for the input variable, then NULL is returned */
   (*mappedvar) = NULL;

   /* if the probnumber is -1, then the master variable is requested.
    * if the probnumber >= 0, then the subproblem variable is requested. */
   if( probnumber == -1 )
   {
      /* getting the original variable for the given pricing variable */
      origvar = GCGpricingVarGetOrigvars(var)[0];

      /* checking whether the original variable is a linking variable.
       * If this is the case, then the corresponding master variable is added to the generated cut.
       * If the pricing variable is not a linking variable, then the farkas dual value is added to the lhs */
      if( GCGoriginalVarIsLinking(origvar) )
      {
         (*mappedvar) = GCGoriginalVarGetMastervars(origvar)[0];
         //(*mappedvar) = SCIPvarGetProbvar((*mappedvar));
      }
   }
   else
   {
      assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

      /* getting the original variable for the given pricing variable */
      origvar = GCGmasterVarGetOrigvars(var)[0];

      /* checking whether the original variable is a linking variable.
       * If this is the case, then the corresponding master variable is added to the generated cut.
       * If the pricing variable is not a linking variable, then the farkas dual value is added to the lhs */
      if( GCGoriginalVarIsLinking(origvar) )
      {
         (*mappedvar) = GCGlinkingVarGetPricingVars(origvar)[probnumber];
         //if( (*mappedvar) != NULL )
            //(*mappedvar) = SCIPvarGetProbvar((*mappedvar));
      }
   }

   return SCIP_OKAY;
}


/** the post-execution method for Benders' decomposition */
static
SCIP_DECL_BENDERSPOSTSOLVE(bendersPostsolveGcg)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);

   /* creates a solution to the original problem */
#ifdef SCIP_DEBUG
   SCIPdebugMessage("The master problem solution.\n");
   SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
   if( type == SCIP_BENDERSENFOTYPE_LP && !infeasible )
   {
      /* if the problem was found to be infeasible, then an artificial solution is created. */
      SCIP_CALL( createOriginalProblemSolution(scip, benders, sol, infeasible) );
      SCIP_CALL( GCGrelaxUpdateCurrentSol(bendersdata->origprob) );
   }

   return SCIP_OKAY;
}



/** the method for creating the Benders' decomposition subproblem. This method is called during the initialisation stage
 *  (after the master problem was transformed)
 *
 *  This method must create the SCIP instance for the subproblem and add the required variables and constraints. In
 *  addition, the settings required for the solving the problem must be set here. However, some settings will be
 *  overridden by the standard solving method included in the Benders' decomposition framework. If a special solving
 *  method is desired, the user can implement the bendersSolvesubDefault callback.
 */
static
SCIP_DECL_BENDERSCREATESUB(bendersCreatesubGcg)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;
   SCIP* origprob;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);

   origprob = bendersdata->origprob;

   /* TODO: Not sure whether the relevance check is necessary for Benders'. Must check */
   if( TRUE || GCGisPricingprobRelevant(origprob, probnumber) )
   {
      SCIPaddBendersSubproblem(scip, benders, GCGgetPricingprob(origprob, probnumber));

      /* setting the objective coefficients for the subproblems.
       * This is required because the variables are added to the pricing problems with a zero coefficient. In the DW
       * context, this is appropriate because the objective coefficients are constantly changing. In the BD context, the
       * objective coefficients are static, so they only need to be updated once. */
      SCIP_CALL( setSubproblemObjs(benders, probnumber) );
   }


   return SCIP_OKAY;
}



/*
 * Benders' decomposition specific interface methods
 */

/** creates the gcg Benders' decomposition and includes it in SCIP */
SCIP_RETCODE SCIPincludeBendersGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origprob            /**< the SCIP instance of the original problem */
   )
{
   SCIP_BENDERSDATA* bendersdata;
   SCIP_BENDERS* benders;

   /* create gcg Benders' decomposition data */
   SCIP_CALL( SCIPallocMemory(scip, &bendersdata) );
   bendersdata->origprob = origprob;
   bendersdata->relaxsol = NULL;

   benders = NULL;


   /* include Benders' decomposition */
   SCIP_CALL( SCIPincludeBendersBasic(scip, &benders, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY,
         BENDERS_CUTLP, BENDERS_CUTPSEUDO, BENDERS_CUTRELAX, BENDERS_SHAREAUXVARS, bendersGetvarGcg,
         bendersCreatesubGcg, bendersdata) );
   assert(benders != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBendersFree(scip, benders, bendersFreeGcg) );
   SCIP_CALL( SCIPsetBendersInitpre(scip, benders, bendersInitpreGcg) );
   SCIP_CALL( SCIPsetBendersExitsol(scip, benders, bendersExitsolGcg) );
   SCIP_CALL( SCIPsetBendersPostsolve(scip, benders, bendersPostsolveGcg) );

   /* including the default cuts for Benders' decomposition */
   SCIP_CALL( SCIPincludeBendersDefaultCuts(scip, benders) );

   return SCIP_OKAY;
}



/** returns the last relaxation solution */
SCIP_SOL* SCIPbendersGetRelaxSol(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   )
{
   SCIP_BENDERSDATA* bendersdata;

   assert(benders != NULL);
   assert(strcmp(SCIPbendersGetName(benders), BENDERS_NAME) == 0);

   bendersdata = SCIPbendersGetData(benders);

   return bendersdata->relaxsol;
}


/** returns the original problem for the given master problem */
SCIP* GCGbendersGetOrigprob(
   SCIP*                 masterprob          /**< the master problem SCIP instance */
   )
{
   SCIP_BENDERS* benders;
   SCIP_BENDERSDATA* bendersdata;

   assert(masterprob != NULL);

   benders = SCIPfindBenders(masterprob, BENDERS_NAME);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);

   return bendersdata->origprob;
}
