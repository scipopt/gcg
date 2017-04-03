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

/**@file   cons_benders.c
 * @brief  constraint handler for benders decomposition
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include "scip/scip.h"
#include "gcg.h"

#include "cons_benders.h"

#include "relax_gcg.h"
#include "struct_solver.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"

#include "scip/cons_linear.h"
#include "scip/pub_var.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "benders"
#define CONSHDLR_DESC          "constraint handler to execute Benders' Decomposition"
#define CONSHDLR_ENFOPRIORITY    10000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY   10000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */


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


/* TODO: (optional) enable linear or nonlinear constraint upgrading */
#if 0
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#define LINCONSUPGD_PRIORITY          0 /**< priority of the constraint handler for upgrading of linear constraints */
#define NONLINCONSUPGD_PRIORITY       0 /**< priority of the constraint handler for upgrading of nonlinear constraints */
#endif

#define SUBPROBLEM_STAT_ARRAYLEN_TIME 1024                /**< length of the array for Time histogram representation */
#define SUBPROBLEM_STAT_BUCKETSIZE_TIME 10                /**< size of the buckets for Time histogram representation */
#define SUBPROBLEM_STAT_ARRAYLEN_CUTS 1024                /**< length of the array for foundVars histogram representation */
#define SUBPROBLEM_STAT_BUCKETSIZE_CUTS 1                 /**< size of the buckets for foundVars histogram representation */

/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for benders constraints */
//struct SCIP_ConsData
//{
//};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP*                 origprob;           /**< the SCIP instance of the original problem */
   int                   npricingprobs;      /**< number of pricing problems */
   SCIP**                pricingprobs;       /**< pointers to the pricing problems */
   SCIP_Real*            pricingobjvals;     /**< the objective values of each pricing problem in the current iteration */
   int*                  noptimalityprob;    /**< number of optimality cuts created by the subproblem */
   int*                  nfeasibilityprob;    /**< number of feasibility cuts created by the subproblem */
   SCIP_Longint          currnodenr;         /**< current node number in the masterproblem*/
   SCIP_HASHMAP*         mapcons2idx;        /**< hashmap mapping constraints to their index in the conss array */
   SCIP_Real*            score;              /**< score of the pricing problem problems */
   int*                  permu;              /**< current permutation of the pricing problems */
   int                   npricingprobsnotnull; /**< number of non-Null pricing problems*/

   SCIP_VAR**            auxiliaryvars;      /**< the auxiliary variables added to the master problem */
   SCIP_CONS***          optimalitycuts;    /**< array of all optimality cuts */
   SCIP_CONS***          feasibilitycuts;    /**< array of all feasibility cuts */
   int*                  noptimalitycuts;    /**< number of optimality cuts */
   int*                  nfeasibilitycuts;   /**< number of feasibility cuts */
   int*                  maxoptimalitycuts;  /**< maximal number of optimality cuts */
   int*                  maxfeasibilitycuts;  /**< maximal number of optimality cuts */

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
   int*                  nfeasibilitycutsdist;/**< Feasibility cuts found in each subproblem */
   int*                  noptimalitycutsdist;/**< Optimality cuts found in each subproblem */
   double*               subproblemtimedist; /**< Time spent solving each subproblem */

   int*                  nodetimehist;       /**< Histogram of nodetime distribution */
   int*                  optimalitycutshist; /**< Histogram of found optimality cuts distribution */
   int*                  feasibilitycutshist;/**< Histogram of found feasibility cuts distribution */

   int                   eagerage;           /**< iterations since last eager iteration */
};


/*
 * Local methods
 */

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
   SCIP_CONSHDLR*        conshdlr            /**< the benders' decomposition constraint handler */
   )
{
   SCIP* pricingprob;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR** probvars;
   int nprobvars;
   int i;
   int j;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* set objective value of all variables in the pricing problems to 0 (for farkas pricing) /
    * to the original objective of the variable (for redcost pricing)
    */
   for( i = 0; i < conshdlrdata->npricingprobs; i++ )
   {
      pricingprob = conshdlrdata->pricingprobs[i];
      if( pricingprob == NULL )
         continue;

      probvars = SCIPgetVars(pricingprob);
      nprobvars = SCIPgetNVars(pricingprob);

      for( j = 0; j < nprobvars; j++ )
      {
         assert(GCGvarGetBlock(probvars[j]) == i);
         assert(GCGoriginalVarIsLinking(GCGpricingVarGetOrigvars(probvars[j])[0]) || (GCGvarGetBlock(GCGpricingVarGetOrigvars(probvars[j])[0]) == i));

         SCIP_CALL( SCIPchgVarObj(pricingprob, probvars[j], varGetObj(probvars[j])));

         SCIPdebugMessage("pricingobj var <%s> %f\n", SCIPvarGetName(probvars[j]), varGetObj(probvars[j]));
      }
   }

   return SCIP_OKAY;
}

/* ensures the size of the optimality cuts array is large enough for a particular problem */
static
SCIP_RETCODE ensureSizeOptimalityCutsArray(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< the constraint handler data */
   int                   probnumber,         /**< the pricing problem number */
   int                   size                /**< needed size */
   )
{
   assert(masterprob != NULL);
   assert(conshdlrdata != NULL);
   assert(probnumber >= 0 && probnumber < conshdlrdata->npricingprobs);

   if( conshdlrdata->maxoptimalitycuts[probnumber] < size )
   {
      int oldsize = conshdlrdata->maxoptimalitycuts[probnumber];
      conshdlrdata->maxoptimalitycuts[probnumber] = SCIPcalcMemGrowSize(masterprob, size);
      SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &conshdlrdata->optimalitycuts[probnumber], oldsize,
            conshdlrdata->maxoptimalitycuts[probnumber]) );
   }
   assert(conshdlrdata->maxoptimalitycuts[probnumber] >= size);

   return SCIP_OKAY;
}

/* ensures the size of the feasibility cuts array is large enough for a particular problem */
static
SCIP_RETCODE ensureSizeFeasibilityCutsArray(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< the constraint handler data */
   int                   probnumber,         /**< the pricing problem number */
   int                   size                /**< needed size */
   )
{
   assert(masterprob != NULL);
   assert(conshdlrdata != NULL);
   assert(probnumber >= 0 && probnumber < conshdlrdata->npricingprobs);

   if( conshdlrdata->maxfeasibilitycuts[probnumber] < size )
   {
      int oldsize = conshdlrdata->maxfeasibilitycuts[probnumber];
      conshdlrdata->maxfeasibilitycuts[probnumber] = SCIPcalcMemGrowSize(masterprob, size);
      SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &conshdlrdata->feasibilitycuts[probnumber], oldsize,
            conshdlrdata->maxfeasibilitycuts[probnumber]) );
   }
   assert(conshdlrdata->maxfeasibilitycuts[probnumber] >= size);

   return SCIP_OKAY;
}

/* this function fixes the linking variables to the value from the master problem */
static
SCIP_RETCODE setupSubproblems(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP**                subproblems,        /**< the SCIP instances of the subproblems */
   SCIP_SOL*             sol,                /**< the current incumbent solution to check. Can be NULL is generating cuts for LP solutions */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< the constraint handler data */
   )
{
   SCIP_VAR** mastervars;
   SCIP_VAR** masterfixedvars;
   SCIP_VAR** pricingvars;       /* the pricing vars that are related to the linking vars */
   SCIP_Real fixedval;
   int nmastervars;
   int nmasterfixedvars;
   int npricingprobs;
   int i;
   int j;


   assert(masterprob != NULL);
   assert(subproblems != NULL);

   mastervars = SCIPgetVars(masterprob);
   masterfixedvars = SCIPgetFixedVars(masterprob);
   nmastervars = SCIPgetNVars(masterprob);
   nmasterfixedvars = SCIPgetNFixedVars(masterprob);

   npricingprobs = conshdlrdata->npricingprobs;

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

         /* fixing the corresponding pricing variables in all subproblems */
         for( j = 0; j < npricingprobs; j++ )
         {
            if( pricingvars[j] != NULL )
            {
               SCIP_Bool infeasible;
               SCIP_Bool fixed;

               SCIP_CALL( SCIPfixVar(subproblems[j], pricingvars[j], fixedval, &infeasible, &fixed) );

               assert(!infeasible);
               assert(fixed);
            }
         }
      }
   }

   return SCIP_OKAY;
}

/* computing as standard Benders' optimality cut from the dual solutions of the LP */
static
SCIP_RETCODE computeStandardOptimalityCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 pricingprob,        /**< the SCIP instance of the pricing problem */
   SCIP_CONSHDLR*        conshdlr,           /**< the benders' decomposition constraint handler */
   SCIP_CONS*            cut                 /**< the cut that is generated from the pricing problem */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** fixedvars;
   SCIP_CONS** conss;
   int nvars;
   int nfixedvars;
   int nconss;
   SCIP_Real dualsol;
   SCIP_Real addval;    /* the value that must be added to the lhs */
   SCIP_Real lhs;       /* the left hand side of the cut */
   int i;

#ifndef NDEBUG
   SCIP_Real verifyobj = 0;
   SCIP_Real checkobj = 0;
#endif

   assert(masterprob != NULL);
   assert(pricingprob != NULL);
   assert(conshdlr != NULL);
   assert(cut != NULL);

   nvars = SCIPgetNVars(pricingprob);
   vars = SCIPgetVars(pricingprob);
   nfixedvars = SCIPgetNFixedVars(pricingprob);
   fixedvars = SCIPgetFixedVars(pricingprob);

   nconss = SCIPgetNConss(pricingprob);
   conss = SCIPgetConss(pricingprob);


   /* looping over all constraints and setting the coefficients of the cut */
   for( i = 0; i < nconss; i++ )
   {
      dualsol = GCGconsGetDualsol(pricingprob, conss[i]);

      assert( !SCIPisInfinity(pricingprob, dualsol) && !SCIPisInfinity(pricingprob, -dualsol) );

      if( SCIPisZero(pricingprob, dualsol) )
         continue;

      lhs = SCIPgetLhsLinear(masterprob, cut);


      if( SCIPisPositive(pricingprob, dualsol) )
         addval = dualsol*GCGconsGetLhs(pricingprob, conss[i]);
      else if( SCIPisNegative(pricingprob, dualsol) )
         addval = dualsol*GCGconsGetRhs(pricingprob, conss[i]);

      lhs += addval;

      /* Update the lhs of the cut */
      SCIP_CALL( SCIPchgLhsLinear(masterprob, cut, lhs) );
   }

   /* looping over all variables to update the coefficients in the computed cut. */
   for( i = 0; i < nvars + nfixedvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_VAR* origvar;
      SCIP_Real redcost;


      if( i < nvars )
         var = vars[i];
      else
         var = fixedvars[i - nvars];

      origvar = GCGpricingVarGetOrigvars(var)[0];

      var = SCIPvarGetProbvar(var);

      redcost = SCIPgetVarRedcost(pricingprob, var);

#ifndef NDEBUG
      checkobj += SCIPvarGetUnchangedObj(var)*SCIPvarGetSol(var, TRUE);
#endif

      /* checking whether the original variable is a linking variable.
       * If this is the case, then the corresponding master variable is added to the generated cut.
       * If the pricing variable is not a linking variable, then the farkas dual value is added to the lhs */
      if( GCGoriginalVarIsLinking(origvar) )
      {
         SCIP_VAR* mastervar;
         SCIP_Real coef;

         /* getting the array of the master variables corresponding to an original variable. It is assumed that the
          * first master variable in this array is the copy of the linking variable in the master problem. */
         mastervar = GCGoriginalVarGetMastervars(origvar)[0];
         mastervar = SCIPvarGetProbvar(mastervar);

         coef = -1.0*(SCIPvarGetObj(var) + redcost);

#ifndef NDEBUG
         verifyobj -= SCIPvarGetObj(var)*SCIPvarGetSol(var, TRUE);
#endif

         SCIP_CALL( SCIPaddCoefLinear(masterprob, cut, mastervar, coef) );
      }
      else
      {
         if( !SCIPisZero(pricingprob, redcost) )
         {
             addval = 0;

             /* get current lhs of the subproblem cut */
             lhs = SCIPgetLhsLinear(masterprob, cut);

             if( SCIPisPositive(pricingprob, redcost) )
                addval = redcost*SCIPvarGetLbLocal(var);
             else if( SCIPisNegative(pricingprob, redcost) )
                addval = redcost*SCIPvarGetUbLocal(var);

             lhs += addval;

             /* Update lhs */
             SCIP_CALL( SCIPchgLhsLinear(masterprob, cut, lhs) );
         }
      }
   }

#ifndef NDEBUG
   lhs = SCIPgetLhsLinear(masterprob, cut);
   verifyobj += lhs;

   /* need to generate a solution to verify the cut. Not sure how to do this yet. */
   //SCIP_CALL( GCGcreateSolFromGcgCol(pricingprob, &sol, probnr, gcgcol) );
   //verifyobj -= SCIPgetActivityLinear(masterprob, cut, sol);
#endif

   //assert(SCIPisFeasEQ(origprob, SCIPgetLPObjval(origprob), verifyobj ));

   assert(cut != NULL);

   return SCIP_OKAY;
}

/* adds an optimality cut to the master problem and updates the conshdlr data */
static
SCIP_RETCODE addOptimalityCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< the constraint handler data */
   SCIP_CONS*            cut,                /**< the cut that is generated from the pricing problem */
   int                   probnumber          /**< the number of the pricing problem */
   )
{
   assert(masterprob != NULL);
   assert(conshdlrdata != NULL);
   assert(cut != NULL);
   assert(probnumber >= 0 && probnumber <= conshdlrdata->npricingprobs);

   SCIP_CALL( ensureSizeOptimalityCutsArray(masterprob, conshdlrdata, probnumber, conshdlrdata->maxoptimalitycuts[probnumber] + 1) );

   conshdlrdata->optimalitycuts[probnumber][conshdlrdata->noptimalitycuts[probnumber]] = cut;

   SCIP_CALL( SCIPcaptureCons(masterprob, cut) );

   ++(conshdlrdata->noptimalitycuts[probnumber]);

   return SCIP_OKAY;
}

/* computing as standard Benders' feasibility cut from the dual solutions of the LP */
/* NOTE: The cut must be created before being passed to this function */
static
SCIP_RETCODE computeStandardFeasibilityCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 pricingprob,        /**< the SCIP instance of the pricing problem */
   SCIP_CONSHDLR*        conshdlr,           /**< the benders' decomposition constraint handler */
   SCIP_CONS*            cut                 /**< the cut that is generated from the pricing problem */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** fixedvars;
   SCIP_CONS** conss;
   int nvars;
   int nfixedvars;
   int nconss;
   SCIP_Real dualsol;
   SCIP_Real lhs;       /* the left hand side of the cut */
   SCIP_Real addval;    /* the value that must be added to the lhs */
   int i;

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   //SCIP_SOL* pricingsol;
   SCIP_Real* farkascoefs;    // the coefficients of the farkas proof
   SCIP_Real farkasact = 0;   // the activities of the farkas proof
   SCIP_Real farkaslhs = 0;   // the lhs of the farkas proof
   int nconsvars;
   int j;

   //pricingsol = SCIPgetBestSol(pricingprob)

   assert(masterprob != NULL);
   assert(pricingprob != NULL);
   assert(conshdlr != NULL);
   assert(cut != NULL);

   nvars = SCIPgetNVars(pricingprob);
   vars = SCIPgetVars(pricingprob);
   nfixedvars = SCIPgetNFixedVars(pricingprob);
   fixedvars = SCIPgetFixedVars(pricingprob);

   nconss = SCIPgetNConss(pricingprob);
   conss = SCIPgetConss(pricingprob);

   SCIP_CALL( SCIPallocBufferArray(pricingprob, &farkascoefs, nvars + nfixedvars) );
   for( i = 0; i < nvars + nfixedvars; i++ )
      farkascoefs[i] = 0;

   /* looping over all constraints and setting the coefficients of the cut */
   for( i = 0; i < nconss; i++ )
   {
      dualsol = GCGconsGetDualfarkas(pricingprob, conss[i]);

      if( SCIPisZero(pricingprob, dualsol) )
         continue;

      lhs = SCIPgetLhsLinear(masterprob, cut);


      if( SCIPisPositive(pricingprob, dualsol) )
         addval = dualsol*GCGconsGetLhs(pricingprob, conss[i]);
      else if( SCIPisNegative(pricingprob, dualsol) )
         addval = dualsol*GCGconsGetRhs(pricingprob, conss[i]);

      lhs += addval;

      /* Update the lhs of the cut */
      SCIP_CALL( SCIPchgLhsLinear(masterprob, cut, lhs) );

      farkaslhs += addval;

      nconsvars = GCGconsGetNVars(pricingprob, conss[i]);
      SCIP_CALL( SCIPallocBufferArray(pricingprob, &consvars, nconsvars) );
      SCIP_CALL( SCIPallocBufferArray(pricingprob, &consvals, nconsvars) );
      SCIP_CALL( GCGconsGetVars(pricingprob, conss[i], consvars, nconsvars) );
      SCIP_CALL( GCGconsGetVals(pricingprob, conss[i], consvals, nconsvars) );

      /* loop over all variables with non-zero coefficient */
      for( j = 0; j < nconsvars; j++ )
      {
         SCIP_VAR* consvar;
         SCIP_Real consval;

         consvar = consvars[j];
         consval = consvals[j];

         /* TODO: Do we need the problem variable? */
         consvar = SCIPvarGetProbvar(consvars[j]);

         //assert(!GCGoriginalVarIsLinking(consvar));

         /* update the coefficient in the farkas activity */
         farkascoefs[SCIPvarGetProbindex(consvar)] += dualsol * consval;
      }

      SCIPfreeBufferArray(pricingprob, &consvars);
      SCIPfreeBufferArray(pricingprob, &consvals);
   }

   /* looping over all variables to update the coefficients in the computed cut. */
   for( i = 0; i < nvars + nfixedvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_VAR* origvar;

      if( i < nvars )
         var = vars[i];
      else
         var = fixedvars[i - nvars];

      origvar = GCGpricingVarGetOrigvars(var)[0];

      var = SCIPvarGetProbvar(var);

      dualsol = farkascoefs[SCIPvarGetProbindex(var)];

      if( SCIPisZero(masterprob, dualsol) )
         continue;

      /* checking whether the original variable is a linking variable.
       * If this is the case, then the corresponding master variable is added to the generated cut.
       * If the pricing variable is not a linking variable, then the farkas dual value is added to the lhs */
      if( GCGoriginalVarIsLinking(origvar) )
      {
         SCIP_VAR* mastervar;

         /* getting the array of the master variables corresponding to an original variable. It is assumed that the
          * first master variable in this array is the copy of the linking variable in the master problem. */
         mastervar = GCGoriginalVarGetMastervars(origvar)[0];

         mastervar = SCIPvarGetProbvar(mastervar);
         var = SCIPvarGetProbvar(var);

         SCIP_CALL( SCIPaddCoefLinear(masterprob, cut, mastervar, dualsol) );
      }
      else
      {
          addval = 0;

          /* get current lhs of the subproblem cut */
          lhs = SCIPgetLhsLinear(masterprob, cut);

          if( SCIPisPositive(pricingprob, dualsol) )
             addval = dualsol*SCIPvarGetUbLocal(var);
          else if( SCIPisNegative(pricingprob, dualsol) )
             addval = dualsol*SCIPvarGetLbLocal(var);

          lhs -= addval;

          /* Update lhs */
          SCIP_CALL( SCIPchgLhsLinear(masterprob, cut, lhs) );
      }
   }

#ifndef NDEBUG
   for( i = 0; i < nvars; ++i )
   {
       SCIP_VAR* var;

       var = vars[i];
       var = SCIPvarGetProbvar(var);

       assert(SCIPvarIsTransformed(var));

       if( farkascoefs[SCIPvarGetProbindex(var)] > 0.0 )
          farkasact += farkascoefs[SCIPvarGetProbindex(var)]*SCIPvarGetUbLocal(var);
       else
          farkasact += farkascoefs[SCIPvarGetProbindex(var)]*SCIPvarGetLbLocal(var);
   }
#endif


   assert(cut != NULL);


#ifndef NDEBUG
   /* TODO: Not sure about how to generate the solution for the first assert. Need to check */
   //assert(SCIPgetActivityLinear(masterprob, cut, pricingsol) < SCIPgetLhsLinear(masterprob, cut));
   assert(farkasact < farkaslhs);
   SCIPfreeBufferArray(pricingprob, &farkascoefs);
#endif


   return SCIP_OKAY;
}

/* adds an feasibility cut to the master problem and updates the conshdlr data */
static
SCIP_RETCODE addFeasibilityCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< the constraint handler data */
   SCIP_CONS*            cut,                /**< the cut that is generated from the pricing problem */
   int                   probnumber          /**< the number of the pricing problem */
   )
{
   assert(masterprob != NULL);
   assert(conshdlrdata != NULL);
   assert(cut != NULL);
   assert(probnumber >= 0 && probnumber <= conshdlrdata->npricingprobs);

   SCIP_CALL( ensureSizeFeasibilityCutsArray(masterprob, conshdlrdata, probnumber, conshdlrdata->maxfeasibilitycuts[probnumber] + 1) );

   conshdlrdata->feasibilitycuts[probnumber][conshdlrdata->nfeasibilitycuts[probnumber]] = cut;

   SCIP_CALL( SCIPcaptureCons(masterprob, cut) );

   ++(conshdlrdata->nfeasibilitycuts[probnumber]);

   return SCIP_OKAY;
}

/* adds the auxiliary variable to the generated cut. If this is the first optimality cut for the subproblem, then the
 * auxiliary variable is first created and added to the master problem. */
static
SCIP_RETCODE addAuxiliaryVariableToCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_CONSHDLR*        conshdlr,           /**< the benders' decomposition constraint handler */
   SCIP_CONS*            cut,                /**< the cut that is generated from the pricing problem */
   int                   probnumber,         /**< the number of the pricing problem */
   SCIP_Bool*            optimal             /**< flag to indicate whether the current subproblem is optimal for the master */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR* auxiliaryvar;
   SCIP_SOL* bestsol;               /* the current solution to the master problem */
   char varname[SCIP_MAXSTRLEN];    /* the name of the auxiliary variable */

   assert(masterprob != NULL);
   assert(conshdlr != NULL);
   assert(cut != NULL);

   (*optimal) = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* retrieving the best solution to check the value of the auxiliary variable against the subproblem solution */
   bestsol = SCIPgetBestSol(masterprob);

   if( conshdlrdata->noptimalitycuts[probnumber] == 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      /* getting the generic branching var add event to drop the var add event while the auxiliary variable is added. */
      eventhdlr = SCIPfindEventhdlr(masterprob, "genericbranchvaradd");

      /* dropping the varadded event of the generic branching event handler */
      SCIP_CALL( SCIPdropEvent(masterprob, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, -1) );

      /* if no optimality cuts have been added for this subproblem, then the auxiliary variable will be created and
       * added */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "auxiliaryvar_%d", probnumber );
      SCIP_CALL( GCGcreateMasterVar(masterprob, conshdlrdata->origprob, conshdlrdata->pricingprobs[probnumber], &auxiliaryvar,
            varname, 1.0, SCIP_VARTYPE_CONTINUOUS, FALSE, probnumber, 0, NULL, NULL));

      SCIP_CALL( SCIPaddVar(masterprob, auxiliaryvar) );

      conshdlrdata->auxiliaryvars[probnumber] = auxiliaryvar;

      SCIP_CALL( SCIPreleaseVar(masterprob, &auxiliaryvar) );

      /* informing the generic branching event handler to catch var added events */
      SCIP_CALL( SCIPcatchEvent(masterprob, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, NULL) );
   }
   else
   {
      SCIP_Real auxiliaryvarval;

      auxiliaryvarval = SCIPgetSolVal(masterprob, bestsol, conshdlrdata->auxiliaryvars[probnumber]);

      /* if the value of the auxiliary variable in the master problem is greater or equal to the subproblem objective,
       * then a cut is not added by the subproblem.
       * TODO: Need to use a epsilon tolerance for this check. */
      if( SCIPisGE(masterprob, auxiliaryvarval, conshdlrdata->pricingobjvals[probnumber]) )
      {
         (*optimal) = TRUE;
         return SCIP_OKAY;
      }
   }

   /* adding the auxiliary variable to the generated cut */
   SCIP_CALL( SCIPaddCoefLinear(masterprob, cut, conshdlrdata->auxiliaryvars[probnumber], 1.0) );

   return SCIP_OKAY;
}

/* generates and applies Benders' cuts */
static
SCIP_RETCODE generateAndApplyBendersCuts(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 pricingprob,        /**< the SCIP instance of the pricing problem */
   SCIP_CONSHDLR*        conshdlr,           /**< the benders' decomposition constraint handler */
   int                   probnumber,         /**< the number of the pricing problem */
   SCIP_RESULT*          result              /**< the result from solving the subproblems */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cut;                  /* the cut that will be generated from the solution to the pricing problem */
   char cutname[SCIP_MAXSTRLEN];    /* the name of the generated cut */
   SCIP_Bool optimal;               /* flag to indicate whether the current subproblem is optimal for the master */

   assert(masterprob != NULL);
   assert(pricingprob != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* setting the name of the generated cut */
   if( SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL )
      (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "optimalitycut_%d_%d", probnumber, conshdlrdata->noptimalitycuts[probnumber] );
   else if( SCIPgetStatus(pricingprob) == SCIP_STATUS_INFEASIBLE )
      (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "feasibilitycut_%d_%d", probnumber, conshdlrdata->nfeasibilitycuts[probnumber] );
   else
      assert(FALSE);

   /* creating the constraint for the cut */
   SCIP_CALL( SCIPcreateConsBasicLinear(masterprob, &cut, cutname, 0, NULL, NULL, 0.0, SCIPinfinity(masterprob)) );

   /* checking the status of the problem to determine whether a feasibility or optimality cut should be generated. */
   if( SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL )
   {
      SCIP_CALL( computeStandardOptimalityCut(masterprob, pricingprob, conshdlr, cut) );
      SCIP_CALL( addAuxiliaryVariableToCut(masterprob, conshdlr, cut, probnumber, &optimal) );

      /* if the current subproblem is optimal for the master, then we do not add a constraint. */
      if( optimal )
      {
         SCIPinfoMessage(masterprob, NULL, "No cut added for subproblem %d\n", probnumber);
         SCIP_CALL( SCIPreleaseCons(masterprob, &cut) );
         return SCIP_OKAY;
      }

      SCIP_CALL( addOptimalityCut(masterprob, conshdlrdata, cut, probnumber) );
   }
   else if( SCIPgetStatus(pricingprob) == SCIP_STATUS_INFEASIBLE )
   {
      if( SCIPgetNLPIterations(pricingprob) == 0 )
         SCIPinfoMessage(masterprob, NULL, "No iterations in pricing problem %d\n", probnumber);

      SCIP_CALL( computeStandardFeasibilityCut(masterprob, pricingprob, conshdlr, cut) );
      SCIP_CALL( addFeasibilityCut(masterprob, conshdlrdata, cut, probnumber) );
   }
   else
      assert(FALSE);

   SCIP_CALL( SCIPprintCons(masterprob, cut, NULL) );
   SCIPinfoMessage(masterprob, NULL, "\n");

   /* adding the constraint to the master problem */
   SCIP_CALL( SCIPaddCons(masterprob, cut) );

   SCIP_CALL( SCIPreleaseCons(masterprob, &cut) );

   (*result) = SCIP_CONSADDED;


   return SCIP_OKAY;
}

/* solving the subproblems to generated Benders' cuts */
static
SCIP_RETCODE solveSubproblems(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_CONSHDLR*        conshdlr,           /**< the benders' decomposition constraint handler */
   SCIP_SOL*             sol,                /**< the current solution, can be null if an lp of pseudo solution is checked */
   SCIP_Real*            objval,             /**< the sum of the objective function values over all subproblems */
   SCIP_Bool*            infeasible,         /**< a flag to indicate whether all subproblems are feasible */
   SCIP_Bool             conscheck           /**< is this function called from conscheck */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int npricingprobs;
   int i;

   /* previous parameter settings */
   int prevCutoffParam;
   int prevPropMaxroundsParam;
   int prevPropMaxroundsRootParam;
   char prevInitAlgParam;
   char prevResolveAlgParam;
   SCIP_Bool prevConfParam;
   SCIP_Bool prevDualParam;

   assert(masterprob != NULL);
   assert(conshdlr != NULL);
   assert(objval != NULL);
   assert(infeasible != NULL);

   (*objval) = 0;
   (*infeasible) = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   npricingprobs = conshdlrdata->npricingprobs;

   SCIP_CALL( SCIPprintSol(masterprob, sol, FALSE, FALSE) );

   /* setting up the subproblems by fixing variables based upon the master problem solution */
   SCIP_CALL( setupSubproblems(masterprob, conshdlrdata->pricingprobs, sol, conshdlrdata) );

   /* TODO: This should be solved just as an LP, so as a MIP. There is too much overhead with the MIP.
    * Need to change status check for checking the LP. */
   for( i = 0; i < npricingprobs; i++ )
   {
      SCIP* pricingprob;
      SCIP_SOL* bestsol;

      pricingprob = conshdlrdata->pricingprobs[i];

      /* modifying all of the parameters */

      /* Do we have to disable presolving? If yes, we have to store all presolving parameters. */
      SCIPsetPresolving(pricingprob, SCIP_PARAMSETTING_OFF, TRUE);

      /* Disabling heuristics so that the problem is not trivially solved */
      SCIPsetHeuristics(pricingprob, SCIP_PARAMSETTING_OFF, TRUE);

      /* store parameters that are changed for the generation of the subproblem cuts */
      SCIP_CALL( SCIPgetBoolParam(pricingprob, "conflict/enable", &prevConfParam) );
      SCIPsetParam(pricingprob, "conflict/enable", FALSE);

      SCIP_CALL( SCIPgetIntParam(pricingprob, "lp/disablecutoff", &prevCutoffParam) );
      SCIPsetIntParam(pricingprob, "lp/disablecutoff", 1);

      SCIP_CALL( SCIPgetCharParam(pricingprob, "lp/initalgorithm", &prevInitAlgParam) );
      SCIPsetCharParam(pricingprob, "lp/initalgorithm", 'd');
      SCIP_CALL( SCIPgetCharParam(pricingprob, "lp/resolvealgorithm", &prevResolveAlgParam) );
      SCIPsetCharParam(pricingprob, "lp/resolvealgorithm", 'd');

      SCIP_CALL( SCIPgetBoolParam(pricingprob, "misc/alwaysgetduals", &prevDualParam) );
      SCIPsetBoolParam(pricingprob, "misc/alwaysgetduals", TRUE);

      //SCIPinfoMessage(pricingprob, NULL, "Pricing problem %d\n", i);
      //SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", (int)SCIP_VERBLEVEL_NORMAL) );
      //SCIP_CALL( SCIPsetBoolParam(pricingprob, "display/lpinfo", TRUE) );

      SCIP_CALL( SCIPgetIntParam(pricingprob, "propagating/maxrounds", &prevPropMaxroundsParam) );
      SCIP_CALL( SCIPsetIntParam(pricingprob, "propagating/maxrounds", 0) );
      SCIP_CALL( SCIPgetIntParam(pricingprob, "propagating/maxroundsroot", &prevPropMaxroundsRootParam) );
      SCIP_CALL( SCIPsetIntParam(pricingprob, "propagating/maxroundsroot", 0) );

      SCIP_CALL( SCIPsetIntParam(pricingprob, "constraints/linear/propfreq", -1) );

      SCIP_CALL( SCIPsolve(pricingprob) );

      bestsol = SCIPgetBestSol(pricingprob);

      if( SCIPgetStatus(pricingprob) == SCIP_STATUS_OPTIMAL )
      {
         (*objval) += SCIPgetSolTransObj(pricingprob, bestsol);
         conshdlrdata->pricingobjvals[i] = SCIPgetSolTransObj(pricingprob, bestsol);
      }
      else if( SCIPgetStatus(pricingprob) == SCIP_STATUS_INFEASIBLE )
      {
         (*infeasible) = TRUE;
         conshdlrdata->pricingobjvals[i] = SCIPinfinity(masterprob);
         /* since the subproblem solves are called from conscheck, then we abort on an infeasible subproblem */
         if( conscheck )
         {
            /* if one subproblem is infeasible, then the current solution is infeasible */
            (*objval) = SCIPinfinity(masterprob);
            break;
         }
      }
      else
         assert(FALSE);

      //SCIP_CALL( SCIPprintStatistics(pricingprob, NULL) );

      SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
      /* resetting the parameter settings to the previous state */
      SCIPsetPresolving(pricingprob, SCIP_PARAMSETTING_DEFAULT, TRUE);
      SCIPsetHeuristics(pricingprob, SCIP_PARAMSETTING_DEFAULT, TRUE);
      SCIPsetBoolParam(pricingprob, "conflict/enable", prevConfParam);
      SCIPsetIntParam(pricingprob, "lp/disablecutoff", prevCutoffParam);
      SCIPsetCharParam(pricingprob, "lp/initalgorithm", prevInitAlgParam);
      SCIPsetCharParam(pricingprob, "lp/resolvealgorithm", prevResolveAlgParam);
      SCIPsetBoolParam(pricingprob, "misc/alwaysgetduals", prevDualParam);
      SCIPsetIntParam(pricingprob, "propagating/maxrounds", prevPropMaxroundsParam);
      SCIPsetIntParam(pricingprob, "propagating/maxroundsroot", prevPropMaxroundsRootParam);
   }

   return SCIP_OKAY;
}

/** free pricing problems */
static
SCIP_RETCODE freePricingProblems(
   SCIP_CONSHDLRDATA*    conshdlrdata
   )
{
   int j;
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->pricingprobs != NULL);

   for( j = 0; j < conshdlrdata->npricingprobs; j++ )
   {
      if( conshdlrdata->pricingprobs[j] != NULL
         && SCIPgetStage(conshdlrdata->pricingprobs[j]) > SCIP_STAGE_PROBLEM)
      {
         SCIP_CALL( SCIPfreeTransform(conshdlrdata->pricingprobs[j]) );
      }
   }

   return SCIP_OKAY;
}
#if 0
/* applies the generated cut to the master problem*/
static
SCIP_RETCODE applyCut(
   //SCIP_BENDERS_CUTTYPE  cuttype
   )
{
   return SCIP_OKAY;
}

/* apply Magnanti-Wong strengthening of the dual solutions from the optimal LP */
static
SCIP_RETCODE applyMagnantiWongDualStrengthening(
   )
{
   return SCIP_OKAY;
}
#endif


/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyBenders NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* free memory for conshdlrdata*/
   if( conshdlrdata != NULL )
   {
      SCIPfreeMemory(scip, &conshdlrdata);
   }

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   /* A possible implementation from gcg_pricer.cpp */
   assert(scip == scip_);
   assert(reducedcostpricing != NULL);
   assert(farkaspricing != NULL);

   SCIP_CALL( solversInit() );

   SCIP_CALL( reducedcostpricing->resetCalls() );
   SCIP_CALL( farkaspricing->resetCalls() );

   return SCIP_OKAY;
}
#else
#define consInitBenders NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitBenders NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreBenders NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreBenders NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolBenders)
{  /*lint --e{715}*/
   SCIP* origprob;
   int i;
   SCIP_CONS** masterconss;
   int nmasterconss;
   int origverblevel;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* TODO: I don't think that I need to check for the scip_. This parameter is only provided in classes related to the
    * pricer. */
   /*
   assert(scip == scip_);
   */
   assert(conshdlr != NULL);
   assert(conshdlrdata != NULL);

   origprob = conshdlrdata->origprob;

   /* at the beginning, the output of the master problem gets the same verbosity level
    * as the output of the original problem */
   SCIP_CALL( SCIPgetIntParam(origprob, "display/verblevel", &origverblevel) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", origverblevel) );

   conshdlrdata->currnodenr = -1;
   conshdlrdata->eagerage = 0;

   nmasterconss = GCGgetNMasterConss(origprob);
   masterconss = GCGgetMasterConss(origprob);

   /* init array containing all pricing problems */
   conshdlrdata->npricingprobs = GCGgetNPricingprobs(origprob);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->pricingprobs), conshdlrdata->npricingprobs) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->pricingobjvals), conshdlrdata->npricingprobs) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->noptimalityprob), conshdlrdata->npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->nfeasibilityprob), conshdlrdata->npricingprobs) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->subproblemcallsdist), conshdlrdata->npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->nfeasibilitycutsdist), conshdlrdata->npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->noptimalitycutsdist), conshdlrdata->npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->subproblemtimedist), conshdlrdata->npricingprobs) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrdata->nodetimehist), SUBPROBLEM_STAT_ARRAYLEN_TIME) ); /*lint !e506*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrdata->optimalitycutshist), SUBPROBLEM_STAT_ARRAYLEN_CUTS) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrdata->feasibilitycutshist), SUBPROBLEM_STAT_ARRAYLEN_CUTS) );

   BMSclearMemoryArray(conshdlrdata->nodetimehist, SUBPROBLEM_STAT_ARRAYLEN_TIME);
   BMSclearMemoryArray(conshdlrdata->optimalitycutshist, SUBPROBLEM_STAT_ARRAYLEN_CUTS);
   BMSclearMemoryArray(conshdlrdata->feasibilitycutshist, SUBPROBLEM_STAT_ARRAYLEN_CUTS);

   conshdlrdata->npricingprobsnotnull = 0;

   for( i = 0; i < conshdlrdata->npricingprobs; i++ )
   {
      conshdlrdata->pricingobjvals[i] = SCIPinfinity(scip);
      conshdlrdata->subproblemcallsdist[i] = 0;
      conshdlrdata->nfeasibilitycutsdist[i] = 0;
      conshdlrdata->noptimalitycutsdist[i] = 0;
      conshdlrdata->subproblemtimedist[i] = 0;

      if( GCGisPricingprobRelevant(origprob, i) )
      {
         conshdlrdata->pricingprobs[i] = GCGgetPricingprob(origprob, i);
         conshdlrdata->npricingprobsnotnull++;
      }
      else
      {
         conshdlrdata->pricingprobs[i] = NULL;
      }
      conshdlrdata->noptimalityprob[i] = 0;
      conshdlrdata->nfeasibilityprob[i] = 0;
   }

   /* alloc memory for arrays of reduced cost */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->score), conshdlrdata->npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->permu), conshdlrdata->npricingprobs) );

   SCIP_CALL( SCIPcreateCPUClock(scip, &(conshdlrdata->freeclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(conshdlrdata->transformclock)) );

   conshdlrdata->solvedsubmipsoptimal = 0;
   conshdlrdata->solvedsubmipsheur = 0;
   conshdlrdata->calls = 0;
   conshdlrdata->pricingiters = 0;

   /* TODO: Get a parameter for the type of cut to be added. Not sure what this means at the moment, but will work this
    * out. */

   SCIP_CALL( SCIPhashmapCreate(&(conshdlrdata->mapcons2idx), SCIPblkmem(scip), 10 * nmasterconss +1) );
   for( i = 0; i < nmasterconss; i++ )
   {
      SCIP_CALL( SCIPhashmapInsert(conshdlrdata->mapcons2idx, masterconss[i], (void*)(size_t)i) );
      assert((int)(size_t)SCIPhashmapGetImage(conshdlrdata->mapcons2idx, masterconss[i]) == i); /*lint !e507*/
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->auxiliaryvars, conshdlrdata->npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->optimalitycuts, conshdlrdata->npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->feasibilitycuts, conshdlrdata->npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->noptimalitycuts, conshdlrdata->npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->nfeasibilitycuts, conshdlrdata->npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->maxoptimalitycuts, conshdlrdata->npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->maxfeasibilitycuts, conshdlrdata->npricingprobs) );

   for( i = 0; i < conshdlrdata->npricingprobs; i++ )
   {
      conshdlrdata->noptimalitycuts[i] = 0;
      conshdlrdata->maxoptimalitycuts[i] = 50;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->optimalitycuts[i], conshdlrdata->maxoptimalitycuts[i]) );

      conshdlrdata->nfeasibilitycuts[i] = 0;
      conshdlrdata->maxfeasibilitycuts[i] = 50;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->feasibilitycuts[i], conshdlrdata->maxfeasibilitycuts[i]) );
   }


   /* setting the objective coefficients for the subproblems.
    * This is required because the variables are added to the pricing problems with a zero coefficient. In the DW
    * context, this is appropriate because the objective coefficients are constantly changing. In the BD context, the
    * objective coefficients are static, so they only need to be updated once. */
   SCIP_CALL( setSubproblemObjs(conshdlr) );

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


#define SCIP_DECL_CONSEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, SCIP_Bool restart)
/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   for( i = conshdlrdata->npricingprobs - 1; i >= 0; i++ )
   {
      SCIPfreeBlockMemoryArray(scip, &conshdlrdata->feasibilitycuts[i], conshdlrdata->maxfeasibilitycuts[i]);
      SCIPfreeBlockMemoryArray(scip, &conshdlrdata->optimalitycuts[i], conshdlrdata->maxoptimalitycuts[i]);
   }

   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->maxfeasibilitycuts, conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->maxoptimalitycuts, conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->nfeasibilitycuts, conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->noptimalitycuts, conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->feasibilitycuts, conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->optimalitycuts, conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->auxiliaryvars, conshdlrdata->npricingprobs);

   SCIPhashmapFree(&(conshdlrdata->mapcons2idx));

   SCIP_CALL( SCIPfreeClock(scip, &(conshdlrdata->transformclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(conshdlrdata->freeclock)) );

   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->permu), conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->score), conshdlrdata->npricingprobs);

   SCIPfreeMemoryArray(scip, &(conshdlrdata->feasibilitycutshist));
   SCIPfreeMemoryArray(scip, &(conshdlrdata->optimalitycutshist));
   SCIPfreeMemoryArray(scip, &(conshdlrdata->nodetimehist));

   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->subproblemtimedist), conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->noptimalitycutsdist), conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->nfeasibilitycutsdist), conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->subproblemcallsdist), conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->nfeasibilityprob), conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->noptimalityprob), conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->pricingobjvals), conshdlrdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->pricingprobs), conshdlrdata->npricingprobs);

   return SCIP_OKAY;
}


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteBenders NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransBenders NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpBenders NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 0
static
SCIP_DECL_CONSSEPALP(consSepalpBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpBenders NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolBenders NULL
#endif


//#define SCIP_DECL_CONSENFOLP(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, int nusefulconss, \
      //SCIP_Bool solinfeasible, SCIP_RESULT* result)
/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   SCIP_Real objval;
   SCIP_Bool infeasible;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   (*result) = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   SCIP_CALL( solveSubproblems(scip, conshdlr, NULL, &objval, &infeasible, FALSE) );

   for( i = 0; i < conshdlrdata->npricingprobs; i++ )
      SCIP_CALL( generateAndApplyBendersCuts(scip, conshdlrdata->pricingprobs[i], conshdlr, i, result) );

   SCIP_CALL( freePricingProblems(conshdlrdata) );

   return SCIP_OKAY;
}


//#define SCIP_DECL_CONSENFORELAX(x) SCIP_RETCODE x (SCIP* scip, SCIP_SOL* sol, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, int nusefulconss, \
      //SCIP_Bool solinfeasible, SCIP_RESULT* result)
/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   SCIP_Real objval;
   SCIP_Bool infeasible;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(sol != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   SCIP_CALL( solveSubproblems(scip, conshdlr, sol, &objval, &infeasible, FALSE) );

   for( i = 0; i < conshdlrdata->npricingprobs; i++ )
      SCIP_CALL( generateAndApplyBendersCuts(scip, conshdlrdata->pricingprobs[i], conshdlr, i, result) );

   SCIP_CALL( freePricingProblems(conshdlrdata) );

   return SCIP_OKAY;
}


//#define SCIP_DECL_CONSENFOPS(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, int nusefulconss, \
      //SCIP_Bool solinfeasible, SCIP_Bool objinfeasible, SCIP_RESULT* result)
/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   SCIP_Real objval;
   SCIP_Bool infeasible;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   SCIP_CALL( solveSubproblems(scip, conshdlr, NULL, &objval, &infeasible, FALSE) );

   for( i = 0; i < conshdlrdata->npricingprobs; i++ )
      SCIP_CALL( generateAndApplyBendersCuts(scip, conshdlrdata->pricingprobs[i], conshdlr, i, result) );

   SCIP_CALL( freePricingProblems(conshdlrdata) );

   return SCIP_OKAY;
}


/* The define is kept as a comment so I know what is being passed to this function */
#if 0
#define SCIP_DECL_CONSCHECK(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, SCIP_SOL* sol, \
      SCIP_Bool checkintegrality, SCIP_Bool checklprows, SCIP_Bool printreason, SCIP_RESULT* result)
#endif
/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckBenders)
{  /*lint --e{715}*/

   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real objval;
   SCIP_Bool infeasible;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(sol != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   SCIP_CALL( solveSubproblems(scip, conshdlr, sol, &objval, &infeasible, TRUE) );

   if( infeasible )
      (*result) = SCIP_INFEASIBLE;
   else  /* in the else case, we need to update the objective function value with objval */
      (*result) = SCIP_FEASIBLE;

   SCIP_CALL( freePricingProblems(conshdlrdata) );

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
SCIP_DECL_CONSPROP(consPropBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropBenders NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolBenders NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropBenders NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockBenders)
{  /*lint --e{715}*/
   //SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   //SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveBenders NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveBenders NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableBenders NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableBenders NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsBenders NULL
#endif


/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintBenders NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyBenders NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseBenders NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsBenders NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsBenders NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsBenders NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for benders constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origprob            /**< SCIP data structure of the original problem */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr;

   /* create benders constraint handler data */
   conshdlrdata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->origprob = origprob;
   conshdlrdata->solvers = NULL;
   conshdlrdata->nsolvers = 0;
   conshdlrdata->nodetimehist = NULL;
   conshdlrdata->optimalitycutshist = NULL;
   conshdlrdata->feasibilitycutshist = NULL;

   conshdlr = NULL;

   /* include constraint handler */
#if 0
   /* use SCIPincludeConshdlr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING,
         conshdlrCopyBenders,
         consFreeBenders, consInitBenders, consExitBenders,
         consInitpreBenders, consExitpreBenders, consInitsolBenders, consExitsolBenders,
         consDeleteBenders, consTransBenders, consInitlpBenders,
         consSepalpBenders, consSepasolBenders, consEnfolpBenders, consEnforelaxBenders, consEnfopsBenders, consCheckBenders,
         consPropBenders, consPresolBenders, consRespropBenders, consLockBenders,
         consActiveBenders, consDeactiveBenders,
         consEnableBenders, consDisableBenders, consDelvarsBenders,
         consPrintBenders, consCopyBenders, consParseBenders,
         consGetVarsBenders, consGetNVarsBenders, consGetDiveBdChgsBenders, conshdlrdata) );
#else
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpBenders, consEnfopsBenders, consCheckBenders, consLockBenders,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveBenders) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyBenders, consCopyBenders) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveBenders) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteBenders) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsBenders) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableBenders) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableBenders) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitBenders) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreBenders) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolBenders) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeBenders) );
   SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsBenders) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsBenders) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsBenders) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitBenders) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreBenders) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolBenders) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpBenders) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseBenders) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolBenders, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintBenders) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropBenders, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropBenders) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpBenders, consSepasolBenders, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransBenders) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxBenders) );

#endif


   /* parameters for the Benders' decomposition constraint handler */
   SCIP_CALL( SCIPaddBoolParam(origprob, "benders/subproblem/useheursolving",
         "should subproblem solving be performed heuristically before solving the LPs to optimality?",
         &conshdlrdata->useheursolving, TRUE, DEFAULT_USEHEURSOLVING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "benders/subproblem/abortsolvingbound",
         "should solving be aborted when the objective function is less than the current upper bound?",
         &conshdlrdata->abortsolvebound, TRUE, DEFAULT_ABORTSOLVINGBOUND, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "benders/subproblem/dispinfos",
         "should additional informations concerning the subproblem solving process be displayed?",
         &conshdlrdata->dispinfos, FALSE, DEFAULT_DISPINFOS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "benders/subproblem/sorting",
         "which sorting method should be used to sort the subproblems problems (0 = order of pricing problems, 1 = according to dual solution of convexity constraint, 2 = according to reliability from previous round)",
         &conshdlrdata->sorting, FALSE, DEFAULT_SORTING, 0, 5, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "benders/subproblem/threads",
         "how many threads should be used to concurrently solve the subprolems (0 to guess threads by OpenMP)",
         &conshdlrdata->threads, FALSE, DEFAULT_THREADS, 0, 4096, NULL, NULL) );

   //SCIP_CALL( SCIPsetIntParam(scip, "lp/disablecutoff", DEFAULT_DISABLECUTOFF) );
   //SCIP_CALL( SCIPaddIntParam(origprob, "benders/subproblem/disablecutoff",
         //"should the cutoffbound be applied in master LP solving (0: on, 1:off, 2:auto)?",
         //&conshdlrdata->disablecutoff, FALSE, DEFAULT_DISABLECUTOFF, 0, 2, paramChgdDisablecutoff, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "benders/subproblem/eagerfreq",
            "frequency at which all subproblems should be solved (0 to disable)",
            &conshdlrdata->eagerfreq, FALSE, DEFAULT_EAGERFREQ, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a benders constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
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
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsBenders() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the benders constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("benders constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   /* TODO: create and store constraint specific data here */

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a benders constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsBenders(scip, cons, name, nvars, vars, coefs, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
