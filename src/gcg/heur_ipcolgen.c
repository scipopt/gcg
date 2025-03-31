/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_ipcolgen.c
 * @brief  The integer programming column generation heuristic
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "gcg/heur_ipcolgen.h"
#include "gcg/gcg.h"
#include "gcg/pricer_gcg.h"
#include "gcg/scip_misc.h"

#include "scip/scipdefplugins.h"
#include "scip/scip_event.h"

#define HEUR_NAME             "ipcolgen"
#define HEUR_DESC             "A destroy and repair heuristic for the master problem that uses a modified pricing problem"
#define HEUR_DISPCHAR         'I'
#define HEUR_PRIORITY         -1110000
#define HEUR_FREQ             -1        /**< experiments shows that the best frequency is 5 */
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_DURINGPRICINGLOOP | SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE

#define DEFAULT_MAXNODES      5000LL    /**< maximum number of nodes to regard in the repair problem                 */
#define DEFAULT_MINFIXINGRATE 0.2       /**< minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_MINIMPROVE    0.01      /**< factor by which restricted master should at least improve the incumbent */
#define DEFAULT_MINNODES      50LL      /**< minimum number of nodes to regard in the repair problem                 */
#define DEFAULT_NODESOFS      500LL     /**< number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.5       /**< repair problem nodes in relation to nodes of the original problem       */
#define DEFAULT_LPLIMFAC      2.0       /* factor by which the limit on the number of LP depends on the node limit  */
#define DEFAULT_USELPROWS     FALSE     /**< should repair problem be created out of the rows in the LP rows,
                                          *  otherwise, the copy constructor of the constraints handlers are used*/
#define DEFAULT_COPYCUTS      TRUE      /**< if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                          *  of the original scip be copied to constraints of the subscip        */
#define DEFAULT_SOLVEAUXPROB  TRUE      /**< should an auxiliary problem be solved to find improving solutions */
#define DEFAULT_DUALWEIGHT    0.25      /**< the default value for the dual weight in the pricing objective */
#define DEFAULT_INITDYNAMICPEN 0.1      /**< the default value for the initial dynamic penalty for the pricing objective */
#define DEFAULT_BIGM          1234.56   /**< the big-M penalty for the pricing objective */
#define DEFAULT_WAITNEWSOL    TRUE      /**< should the heuristic wait until a new solution is found before executing */
#define DEFAULT_MININITIALGAP 0.5       /**< the minimum optimality gap necessary prior to the first call of the heuristic.
                                          *  The results from the paper show that this could be reduced to 0.25 */
#define DEFAULT_CALLSPERNODE  4         /**< the number of times that the heuristic is called per node */
#define DEFAULT_MAXITER       10        /**< the maximum number of weighted pricing iterations */
#define DEFAULT_NOIMPROVEITER 3         /**< the number of weighted pricing iterations without primal improvement.
                                          *  Experiments showed significant performance degradation when increasing this value */
#define DEFAULT_RINSFIXING    FALSE     /**< should a RINS style fixing be used for the repair master problem */

#define DEFAULT_RANDSEED      31        /**< initial random seed */

#define DEFAULT_ARRAYSIZE     100
#define PREVLPOBJSIZE         5

/* event handler properties */
#define EVENTHDLR_NAME         "IPColGen"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"

/* pricing callback properties */
#define PRICINGCB_NAME         "IPColGen"
#define PRICINGCB_DESC         "pre- and post-pricing methods for adding relaxation and dual weights to the pricing problem objective"
#define PRICINGCB_PRIORITY     1000000
#define PRICINGCB_EXCLUSIVE    TRUE

/*
 * Data structures
 */

enum SCIP_HeurConsType
{
   CONSTYPE_SETPACK     = 0,
   CONSTYPE_SETCOVER    = 1,
   CONSTYPE_SETPART     = 2,
   CONSTYPE_INVALID     = 3,
   CONSTYPE_NCONSTYPES  = 4
};
typedef enum SCIP_HeurConsType SCIP_HEURCONSTYPE;

/** primal heuristic data */
struct SCIP_HeurData
{
   GCG*                  gcg;                /**< GCG data structure */
   SCIP_SOL*             lastsol;            /**< the last incumbent solution used as reference point */
   SCIP_SOL*             partialsol;         /**< the partial solution that is used for the weighted pricing */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_NODE*            prevnode;           /**< the node where this heuristic was previously called */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the repair problem */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the repair problem */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   SCIP_Longint          usednodes;          /**< nodes already used by repair problem in earlier calls */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed */
   SCIP_Real             minimprove;         /**< factor by which ipcolgen should at least improve the incumbent */
   SCIP_Real             nodesquot;          /**< repair problem nodes in relation to nodes of the original problem */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP for the repair problem,
                                              *   for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   SCIP_Real             dualweight;         /**< the default value for the dual weight in the pricing objective */
   SCIP_Real             initdynamicpen;     /**< the initial value for the dynamic penalty for the pricing objective */
   SCIP_Bool             uselprows;          /**< should the repair problem be created out of the rows in the LP rows? */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in repair problem? */
   SCIP_Bool             waitnewsol;         /**< should the heuristic wait until a new solution is found before executing */
   SCIP_Bool             solveauxproblem;    /**< should the repair problem be solved to find improving solutions */
   SCIP_Real             mininitialgap;      /**< the minimum initial gap required before the first call of the heuristic */
   int                   callspernode;       /**< the number of times that the heuristic is called per node */
   int                   maxiter;            /**< the maximum number of weighted pricing iterations */
   int                   noimproveiter;      /**< the number of weighted pricing iterations without primal improvement */
   SCIP_Bool             rinsfixing;         /**< should a RINS style fixing be used for the repair master problem */

   /* parameters used for controlling the execution of the ipcolgen heuristic */
   SCIP_Bool             inheur;             /**< flag to avoid recursive calls to the heuristic */
   int                   numexec;            /**< the number of times the heuristic is executed */
   SCIP_Bool             firstexec;          /**< flag to indicate whether the first execution of the heuristic has been performed */
   SCIP_Longint          prevpricingiter;    /**< the number of pricing iterations performed before this node */
   SCIP_Real*            prevlpobjs;         /**< the LP objectives during the previous call to this heuristic */
   int                   nprevlpobjs;        /**< the number of previous LP objectives stored */
   int                   prevlpobjssize;     /**< the size of previous LP objectives array */
   int                   prevnsolsfound;     /**< the previous number of solutions found */
   SCIP_Longint          prevnnodes;         /**< the number of nodes processed in the previous call to the heuristic */
   SCIP_Longint          nwaitnodes;         /**< the number of nodes to wait until the heuristic is called next */

   /* statistics */
   SCIP_Real             firstcallgap;       /**< the optimality gap when the heuristic is first called */
   SCIP_Real             firstcallabsgap;    /**< the absolute gap when the heuristic is first called */

   /* the data required for the pricing callback functions */
   SCIP_SOL*             bestsol;            /**< the current best solution during the pricing loop */
   SCIP_Real*            penalties;          /**< the penalties applied to each constraint in the pricing objective */
   int*                  penaltiesids;       /**< the ids corresponding to the constraints to apply the penalties */
   IPC_PENALTYTYPE*      penaltytypes;       /**< the type of penalty, this is used for the penalty adjustment */
   int                   npenalties;         /**< the number of penalties */
   int                   penaltiessize;      /**< the length of the penalties array */
   int                   nmastervars;        /**< the number of variables in the master problem */
   int                   npricingprobs;      /**< the number of pricing problems */
   int                   nfixedvars;         /**< the number of variables fixed to form the partial solution */
   int                   infeascount;        /**< the number of iterations that the partial solution is infeasible */
   SCIP_Bool             abort;              /**< should the weighted pricing be aborted */
};

/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process each time the LP is solved. This is used to interrupt the solving process of the
 * repair problem
 */
static
SCIP_DECL_EVENTEXEC(eventExecIpcolgen)
{
   SCIP_HEURDATA* heurdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_LPSOLVED);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata != NULL);

   /* interrupt solution process of sub-SCIP */
   if( SCIPgetNLPs(scip) > heurdata->lplimfac * heurdata->nodelimit )
   {
      SCIPdebugMsg(scip, "interrupt after  %" SCIP_LONGINT_FORMAT " LPs\n",SCIPgetNLPs(scip));
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}


/** creates a starting solution for the heuristic if no solution has been previously found
 *
 *  The start solution is created by greedily adding the most recently generated columns. Once a column is added for
 *  each pricing problem, then the start solution has been created.
 */
static
SCIP_RETCODE createStartSolution(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP_SOL*             startsol,           /**< the starting solution that will be populated with variables */
   SCIP_Bool*            success             /**< pointer to store whether creating the starting solution was successful */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;
   int nvars;
   int npricingprobs;

   SCIP_Bool* blockvaradded;
   int nvarsadded;
   int i;

   assert(gcg != NULL);
   assert(startsol != NULL);

   scip = GCGgetMasterprob(gcg);

   /* getting the master problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   npricingprobs = GCGgetNPricingprobs(gcg);
   SCIP_CALL( SCIPallocClearBufferArray(scip, &blockvaradded, npricingprobs) );

   /* adding a variable for each block to the solution */
   nvarsadded = 0;
   for( i = nvars - 1; i >= 0 && nvarsadded < npricingprobs; i-- )
   {
      SCIP_Real varval;
      int blocknum;

      /* getting the block number for the variable */
      blocknum = GCGvarGetBlock(vars[i]);
      assert(blocknum >= 0 && blocknum < npricingprobs);

      /* if the block number is negative, then the variable is a master-only variable. Thus, it can be ignored. */
      if( blocknum < 0 )
         continue;

      /* if a variable has been added for the block, then the variable is ignored */
      if( blockvaradded[blocknum] )
         continue;

      /* adding the variable to the solution by setting the value to the variables upper bound. If the upper bound is
       * infinity, then the variable is set to 1.0
       */
      varval = SCIPvarGetUbGlobal(vars[i]);
      if( SCIPisFeasGE(scip, varval, SCIPinfinity(scip)) )
         varval = 1.0;

      SCIP_CALL( SCIPsetSolVal(scip, startsol, vars[i], varval) );

      blockvaradded[blocknum] = TRUE;
      nvarsadded++;
   }

   /* at least half of the blocks must have a variable added to the solution */
   (*success) = (nvarsadded >= npricingprobs/2.0);

   /* freeing the buffer memory */
   SCIPfreeBufferArray(scip, &blockvaradded);

   return SCIP_OKAY;
}


/** creates a partial solution by destroying a source complete solution */
static
SCIP_RETCODE createPartialSolution(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_SOL*             partialsol,         /**< the new partial solution */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   SCIP_Real             minfixingrate,      /**< percentage of integer variables that have to be fixed */
   int                   nblocks,            /**< the number of blocks in the problem */
   int*                  nfixedvars,         /**< returns the number of variables fixed to zero in the partial solution */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully  */
   )
{
   SCIP_VAR** mastervars;
   SCIP_Real mastersolval;
   int nmastervars;

   SCIP_Real fixingrate;
   int fixingcounter;
   int i;

   assert(scip != NULL);
   assert(partialsol != NULL);

   /* get variable data of the master problem */
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   i = 0;
   fixingrate = 0.0;
   fixingcounter = 0;

   /* setting variables to zero variables in the partial solution */
   while( i < nmastervars && fixingrate < minfixingrate)
   {
      int varidx;

      varidx = SCIPrandomGetInt(randnumgen, 0, nmastervars - 1);

      mastersolval = SCIPgetSolVal(scip, partialsol, mastervars[varidx]);

      /* if variable takes a non-zero value in the master solution, then it is fixed to zero.
       * the master variable must be a priced variable
       */
      if( GCGvarGetBlock(mastervars[varidx]) >= 0 && !SCIPisFeasZero(scip, mastersolval) )
      {
         SCIPdebugMsg(scip, "Fixing <%s> to zero (%g)\n", SCIPvarGetName(mastervars[varidx]), mastersolval);
         SCIP_CALL( SCIPsetSolVal(scip, partialsol, mastervars[varidx], 0.0) );
         fixingcounter++;
      }

      i++;
      fixingrate = fixingcounter / (SCIP_Real)(MAX(nblocks, 1));
   }

   (*nfixedvars) = fixingcounter;

   /* abort, if no variables are fixed (which should not happen) */
   if( fixingcounter == 0 )
   {
      SCIPdebugMessage(" -> no master variables fixed, not solving problem.\n");
      *success = FALSE;
      return SCIP_OKAY;
   }

   SCIPdebugMessage(" -> %d out of %d (%.2f percent) blocks fixed.\n", fixingcounter, nblocks, fixingrate * 100.0);

   (*success) = TRUE;

   return SCIP_OKAY;
}

/** adds a penalty to the dynamic penalties */
static
SCIP_RETCODE addPenalty(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_Real             addpenalty,         /**< the penalty to be added */
   int                   addpenaltyid,       /**< the constraint id for the added penalty */
   IPC_PENALTYTYPE       addpenaltytype,     /**< the type of the added penalty */
   SCIP_Real**           penalties,          /**< the relaxation penaltys */
   int**                 penaltiesids,       /**< the relaxtion penalty constraint ids */
   IPC_PENALTYTYPE**     penaltytypes,       /**< the penalty types */
   int*                  npenalties,         /**< the number of dynamic penalties */
   int*                  penaltiessize       /**< the number of elements in the penalties array */
   )
{
   assert(scip != NULL);
   assert(penalties != NULL);
   assert(penaltiesids != NULL);
   assert(npenalties != NULL);
   assert(penaltiessize != NULL);

   if( *npenalties >= *penaltiessize )
   {
      int oldsize = *penaltiessize;
      *penaltiessize = SCIPcalcMemGrowSize(scip, *penaltiessize + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, penalties, oldsize, *penaltiessize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, penaltiesids, oldsize, *penaltiessize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, penaltytypes, oldsize, *penaltiessize) );
   }
   assert(*npenalties < *penaltiessize);

   (*penalties)[*npenalties] = addpenalty;
   (*penaltiesids)[*npenalties] = addpenaltyid;
   (*penaltytypes)[*npenalties] = addpenaltytype;
   (*npenalties)++;

   return SCIP_OKAY;
}

/** checks whether the constraint is of a valid type and the sign of the penalty is returned */
static
SCIP_HEURCONSTYPE getConstraintType(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_CONS*            cons                /**< the constraint to check for validity */
   )
{
   consType ctype;

   /* getting the constraint type */
   ctype = GCGconsGetType(scip, cons);

   if( ctype == setpacking )
      return CONSTYPE_SETPACK;
   if( ctype == setcovering )
      return CONSTYPE_SETCOVER;
   if( ctype == setpartitioning )
      return CONSTYPE_SETPART;

   if( ctype == linear )
   {
      SCIP_Real* consvals;
      int nconsvars;
      int i;

      consvals = SCIPgetValsLinear(scip, cons);
      nconsvars = SCIPgetNVarsLinear(scip, cons);

      for( i = 0; i < nconsvars; i++ )
      {
         if( !(SCIPisZero(scip, consvals[i]) || SCIPisEQ(scip, consvals[i], 1.0)) )
            return CONSTYPE_INVALID;
      }

      if( SCIPisInfinity(scip, -SCIPgetLhsLinear(scip, cons)) && SCIPisEQ(scip, SCIPgetRhsLinear(scip, cons), 1.0) )
         return CONSTYPE_SETPACK;

      if( SCIPisInfinity(scip, SCIPgetRhsLinear(scip, cons)) && SCIPisEQ(scip, SCIPgetLhsLinear(scip, cons), 1.0) )
         return CONSTYPE_SETCOVER;

      if( SCIPisEQ(scip, SCIPgetLhsLinear(scip, cons), 1.0) && SCIPisEQ(scip, SCIPgetRhsLinear(scip, cons), 1.0) )
         return CONSTYPE_SETPART;
   }

   return CONSTYPE_INVALID;
}

/** using the current solution, compute the dynamic penalties to apply to the column generation pricing problem */
static
SCIP_RETCODE computeDynamicPenalties(
   GCG*                  gcg,               /**< the SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution used to set up the repair problem */
   SCIP_Real             initdynamicpen,     /**< the value of the initial dynamic penalty */
   SCIP_Real**           penalties,          /**< the dynamic penalties, or NULL if the penalties are not computed */
   int**                 penaltiesids,       /**< the dynamic penalties constraint ids, or NULL if the penalties are not computed */
   IPC_PENALTYTYPE**     penaltytypes,       /**< the penalty types */
   int*                  npenalties,         /**< the number of dynamic penalties, or NULL if the penalties are not computed */
   int*                  penaltiessize,      /**< the number of elements in the dynamic penalties array, or NULL if the penalties are not computed */
   SCIP_Bool             repairsol           /**< should variables be fixed to zero to repair the solution */
   )
{
   SCIP* scip;
   SCIP_CONS** masterconss;
   int nmasterconss;
   int i;
   int j;

   assert(gcg != NULL);
   assert(sol != NULL);

   scip = GCGgetMasterprob(gcg);

   nmasterconss = GCGgetNMasterConss(gcg);
   masterconss = GCGgetMasterConss(gcg);

   /* setting the dynamic penalties per constraint. The penalties can only be applied to set covering or set packing
    * constraints. This is because all constraint coefficients are positive and the LHS != RHS.
    * TODO: extend this to handling more general constraints
    */
   for( i = 0; i < nmasterconss; i++ )
   {
      SCIP_HEURCONSTYPE constype;

      constype = getConstraintType(scip, masterconss[i]);

      /* if the constraint is of type SET COVERING or SET PACKING, then we add penalties to the pricing objective.
       * Otherwise, no penalties will be added
       */
      if( constype == CONSTYPE_SETPACK || constype == CONSTYPE_SETCOVER || constype == CONSTYPE_SETPART )
      {
         SCIP_VAR** consvars;
         int nconsvars;
         SCIP_Real sum;
         SCIP_Real conspenalty;

         nconsvars = GCGconsGetNVars(scip, masterconss[i]);
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );

         if( nconsvars > 0 )
         {
            SCIP_CALL( GCGconsGetVars(scip, masterconss[i], consvars, nconsvars) );
         }

         /* computing the activity of the constraint */
         sum = 0.0;
         for( j = 0; j < nconsvars; j++ )
            sum += SCIPgetSolVal(scip, sol, consvars[j]);

         SCIPdebugMsg(scip, "constraint <%s>  activity: %g\n", SCIPconsGetName(masterconss[i]), sum);

         /* if the solution should be repaired, then any set packing or set partitioning constraint with an activity
          * currently exceeding 1.0 will have all variables fixed to zero
          */
         if( repairsol && (constype == CONSTYPE_SETPACK || constype == CONSTYPE_SETPART) && SCIPisSumGT(scip, sum, 1.0) )
         {
            SCIP_Real solval;

            j = 0;
            while( SCIPisSumGT(scip, sum, 1.0) && j < nconsvars )
            {
               solval = SCIPgetSolVal(scip, sol, consvars[j]);
               if( SCIPisGE(scip, solval, 1.0) )
               {
                  SCIP_CALL( SCIPsetSolVal(scip, sol, consvars[j], 0.0) );
                  sum -= solval;
               }
               j++;
            }
         }

         /* for set packing, if the constraint is satisfied, then we must penalise any non-zeros in that constraint.
          * for set covering, if the constraint is violated, then we must reward any non-zeros in that constraint
          */
         if( penalties != NULL )
         {
            assert(penaltiesids != NULL);
            assert(penaltytypes != NULL);
            assert(npenalties != NULL);
            assert(penaltiessize != NULL);

            conspenalty = initdynamicpen*MAX(sum, 1.0);
            if( SCIPisSumGE(scip, sum, 1.0) )
            {
               if( constype == CONSTYPE_SETPART || constype == CONSTYPE_SETPACK )
               {
                  SCIP_CALL( addPenalty(scip, DEFAULT_BIGM, i, PENALTYTYPE_BIGM, penalties, penaltiesids, penaltytypes,
                        npenalties, penaltiessize) );
               }
            }
            else
            {
               if( constype == CONSTYPE_SETPACK )
               {
                  SCIP_CALL( addPenalty(scip, conspenalty, i, PENALTYTYPE_SETPACK, penalties, penaltiesids, penaltytypes,
                        npenalties, penaltiessize) );
               }
               else if( constype == CONSTYPE_SETCOVER )
               {
                  SCIP_CALL( addPenalty(scip, -1.0*conspenalty, i, PENALTYTYPE_SETCOVER, penalties, penaltiesids, penaltytypes,
                        npenalties, penaltiessize) );
               }
               else if( constype == CONSTYPE_SETPART )
               {
                  SCIP_CALL( addPenalty(scip, conspenalty, i, PENALTYTYPE_SETPACK, penalties, penaltiesids, penaltytypes,
                        npenalties, penaltiessize) );
                  SCIP_CALL( addPenalty(scip, -1.0*conspenalty, i, PENALTYTYPE_SETCOVER, penalties, penaltiesids, penaltytypes,
                        npenalties, penaltiessize) );
               }
            }
         }

         SCIPfreeBufferArray(scip, &consvars);
      }
   }

   return SCIP_OKAY;
}

/** adjusts the dynamic penalties w.r.t to the added columns */
static
SCIP_RETCODE adjustDynamicPenalties(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP_HEURDATA*        heurdata            /**< the heuristic data */
   )
{
   SCIP* scip;
   SCIP_SOL* augsol;
   SCIP_VAR** mastervars;
   SCIP_CONS** masterconss;
   int nmastervars;
   int i;
   int j;

   SCIP_Bool penaltychanged;
   SCIP_Bool abort;
   SCIP_Bool infinitepenalty;

   assert(gcg != NULL);
   assert(heurdata->penalties != NULL);
   assert(heurdata->penaltiesids != NULL);
   assert(heurdata->penaltytypes != NULL);

   scip = GCGgetMasterprob(gcg);

   /* getting the master variables. It is assumed that the variables at the end of the array correspond to the newly
    * added variables
    */
   mastervars = SCIPgetVars(scip);
   nmastervars = SCIPgetNVars(scip);

   masterconss = GCGgetMasterConss(gcg);

   /* if no columns have been added to the master problem, then no penalty update will be performed */
   if( heurdata->nmastervars >= nmastervars )
      return SCIP_OKAY;

   /* creating the new augmentation solution */
   SCIP_CALL( SCIPcreateSol(scip, &augsol, NULL) );

   /* setting the values in the augmentation and partial solution */
   for( i = heurdata->nmastervars; i < nmastervars; i++ )
   {
      SCIP_CALL( SCIPsetSolVal(scip, augsol, mastervars[i], 1.0) );
   }

   /* scanning all of the constraints that have been previously identified and computing their activity with respect to
    * the new priced columns
    */
   abort = TRUE;
   penaltychanged = FALSE;
   infinitepenalty = FALSE;
   for( i = 0; i < heurdata->npenalties; i++ )
   {
      IPC_PENALTYTYPE penaltytype;
      SCIP_CONS* cons;
      SCIP_VAR** consvars;
      int nconsvars;
      SCIP_Real sum;
      SCIP_Real factor;

      penaltytype = heurdata->penaltytypes[i];

      /* if the penalty is a BIGM penalty, then the penalty is not adjusted */
      if( penaltytype == PENALTYTYPE_BIGM )
         continue;

      assert(heurdata->penaltiesids[i] <= GCGgetNMasterConss(gcg));
      cons = masterconss[heurdata->penaltiesids[i]];

      nconsvars = GCGconsGetNVars(scip, cons);
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );

      if( nconsvars > 0 )
      {
         SCIP_CALL( GCGconsGetVars(scip, cons, consvars, nconsvars) );
      }

      /* computing the activity of the constraint */
      sum = 0.0;
      for( j = 0; j < nconsvars; j++ )
         sum += SCIPgetSolVal(scip, augsol, consvars[j]);

      SCIPdebugMsg(scip, "adjusting penalties -- constraint <%s>  activity: %g\n", SCIPconsGetName(cons), sum);

      /* if the activity is positive, this means that the new variable is covering the constraints. As such, the penalties
       * must be adjusted
       * for set packing, the penalty is increased
       * for set covering, the penalty is decreased
       */
      factor = 3*MAX(sum, 1.0);
      if( SCIPisSumGE(scip, sum, 1.0) )
      {
         if( penaltytype == PENALTYTYPE_SETPACK )
         {
            heurdata->penalties[i] *= factor;
            penaltychanged = TRUE;
         }
         else if( penaltytype == PENALTYTYPE_SETCOVER )
         {
            heurdata->penalties[i] /= factor;
            penaltychanged = TRUE;
         }
      }
      else
      {
         if( penaltytype == PENALTYTYPE_SETPACK )
         {
            heurdata->penalties[i] /= factor;
            penaltychanged = TRUE;
         }
         else if( penaltytype == PENALTYTYPE_SETCOVER )
         {
            heurdata->penalties[i] *= factor;
            penaltychanged = TRUE;
         }
      }

      SCIPdebugMsg(scip, "new penalty for constraint <%s>: %g\n", SCIPconsGetName(cons), heurdata->penalties[i]);

      if ( !SCIPisFeasZero(scip, heurdata->penalties[i]) )
         abort = FALSE;
      else
         heurdata->penalties[i] = 0.0;

      if( SCIPisGT(scip, heurdata->penalties[i], 1e+5) || SCIPisGT(scip, -heurdata->penalties[i], 1e+5) )
         infinitepenalty = TRUE;

      SCIPfreeBufferArray(scip, &consvars);
   }

   /* freeing the solution */
   SCIP_CALL( SCIPfreeSol(scip, &augsol) );

   /* if no penalty has been changed, then the weighted pricing is aborted */
   heurdata->abort = abort || !penaltychanged || infinitepenalty;

   return SCIP_OKAY;
}


/* ---------------- Callback methods of the pricing callback plugin ---------------- */

/** the pre-pricing method of the pricing callback technique
 *
 *  This method is called immediately before the pricing is performed in the GCG pricer. At this point, it is possible
 *  to modify the solving data is used within the pricing for new variables. Any data that is modified should be
 *  reverted in the post-pricing method.
 */
static
GCG_DECL_PRICINGCBPREPRICING(pricingcbPrepricingIpcolgen)
{
#ifdef SCIP_DEBUG
   int i;
#endif
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;

   assert(gcg != NULL);
   assert(pricer != NULL);

   scip = GCGgetMasterprob(gcg);

   (*result) = SCIP_DIDNOTRUN;

   /* the callback is only executed during redcost pricing */
   if( type != GCG_PRICETYPE_REDCOST )
      return SCIP_OKAY;

   /* getting the IP colgen heuristic */
   heur = SCIPfindHeur(scip, HEUR_NAME);
   assert(heur != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* this callback can only be executed during the heuristic */
   if( !heurdata->inheur )
      return SCIP_OKAY;

   /* aborting the weighted pricing */
   if( heurdata->abort )
   {
      (*abort) = TRUE;
      return SCIP_OKAY;
   }

   /* setting the weight and dynamic penalties for the dual values and the master constraints, respectively */
   GCGsetPricingObjDualWeight(gcg, heurdata->dualweight);
   GCGsetPricingObjRelaxWeight(gcg, heurdata->penalties, heurdata->penaltiesids, heurdata->npenalties);

#ifdef SCIP_DEBUG
   for( i = 0; i < heurdata->npenalties; i++ )
   {
      if( heurdata->penalties[i] != 0.0 )
         SCIPdebugMsg(scip, "penalty[%d](%d): %g\n", i, heurdata->penaltiesids[i], heurdata->penalties[i]);
   }
#endif

   (*result) = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** the post-pricing method of the pricing callback technique
 *
 *  This method is called immediately after the pricing is performed in the GCG pricer. This method should be used to
 *  revert any changes made to in the pre-pricing method.
 */
static
GCG_DECL_PRICINGCBPOSTPRICING(pricingcbPostpricingIpcolgen)
{
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;

   assert(gcg != NULL);
   assert(pricer != NULL);

   scip = GCGgetMasterprob(gcg);
   (*result) = SCIP_DIDNOTRUN;

   /* the callback is only executed during redcost pricing */
   if( type != GCG_PRICETYPE_REDCOST )
      return SCIP_OKAY;

   /* getting the IP colgen heuristic */
   heur = SCIPfindHeur(scip, HEUR_NAME);
   assert(heur != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* this callback can only be executed during the heuristic */
   if( !heurdata->inheur )
      return SCIP_OKAY;

   SCIP_CALL( adjustDynamicPenalties(gcg, heurdata) );

   /* resetting the dual weight and the dynamic penalties */
   GCGsetPricingObjDualWeight(gcg, 1.0);
   GCGsetPricingObjRelaxWeight(gcg, NULL, NULL, 0);

   heurdata->nmastervars = SCIPgetNVars(scip);

   (*result) = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * Local methods
 */

/** pricing in new variables with a weighted objective function for the original master scip instance. The new variables
 *  will not form part of any stored solution, so they will not be fixed when setting up the repair problem
 */
static
SCIP_RETCODE performWeightedPricing(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP_HEURDATA*        heurdata,           /**< the heuristic data */
   int                   nblocks             /**< the number of blocks in the problem */
   )
{
   SCIP* scip;
   SCIP* origprob;
   int npricingiter;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;
   SCIP_Bool lpsolved;
   SCIP_Real lpobjval;
#ifdef SCIP_DEBUG
   int i;
#endif

   /* pricing parameter settings */
   int maxcolsroundredcostroot;
   int maxcolsroundredcost;
   int maxroundsredcost;
   int heurpricingiters;
   char sorting;

   assert(gcg != NULL);

   scip = GCGgetMasterprob(gcg);
   origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   SCIPdebugMsg(scip, "Performing weighted pricing\n");

   /* the number of pricing iterations is dependent on the number of pricing problems */
   npricingiter = nblocks;

   /* setting a flag to indicate that the next solve is from within the heuristic */
   heurdata->inheur = TRUE;

   heurdata->abort = FALSE;

   heurdata->infeascount = 0;

   /* storing the previously set pricing parameter settings */
   SCIP_CALL( SCIPgetIntParam(origprob, "pricing/masterpricer/maxcolsroundredcostroot", &maxcolsroundredcostroot) );
   SCIP_CALL( SCIPgetIntParam(origprob, "pricing/masterpricer/maxcolsroundredcost", &maxcolsroundredcost) );
   SCIP_CALL( SCIPgetIntParam(origprob, "pricing/masterpricer/maxroundsredcost", &maxroundsredcost) );
   SCIP_CALL( SCIPgetIntParam(origprob, "pricing/masterpricer/heurpricingiters", &heurpricingiters) );
   SCIP_CALL( SCIPgetCharParam(origprob, "pricing/masterpricer/sorting", &sorting) );

   /* setting the pricing parameters for the repair phase of the heuristic */
   SCIP_CALL( SCIPsetIntParam(origprob, "pricing/masterpricer/maxcolsroundredcostroot", 1) );
   SCIP_CALL( SCIPsetIntParam(origprob, "pricing/masterpricer/maxcolsroundredcost", 1) );
   SCIP_CALL( SCIPsetIntParam(origprob, "pricing/masterpricer/maxroundsredcost", 1) );
   SCIP_CALL( SCIPsetIntParam(origprob, "pricing/masterpricer/heurpricingiters", INT_MAX) );
   SCIP_CALL( SCIPsetCharParam(origprob, "pricing/masterpricer/sorting", 'd') );

#ifdef SCIP_MOREDEBUG
   SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", TRUE) );
#endif

   /* enabling the pricing callback plugin and marking at as exclusive */
   GCGpricingcbSetEnabled(GCGpricerFindPricingcb(gcg, PRICINGCB_NAME), TRUE);
   GCGpricingcbSetExclusive(GCGpricerFindPricingcb(gcg, PRICINGCB_NAME), TRUE);

   /* solving the restricted master problem with the alternative pricing objective function */
   SCIP_CALL( GCGrelaxPerformProbingWithPricing(gcg, npricingiter, NULL, NULL, &lpobjval, &lpsolved,  &lperror, &cutoff) );
   assert(!lperror);

   /* disabling the pricing callback plugin */
   GCGpricingcbSetEnabled(GCGpricerFindPricingcb(gcg, PRICINGCB_NAME), FALSE);
   GCGpricingcbSetExclusive(GCGpricerFindPricingcb(gcg, PRICINGCB_NAME), FALSE);

   /* resetting the pricing parameter settings */
   SCIP_CALL( SCIPsetIntParam(origprob, "pricing/masterpricer/maxcolsroundredcostroot", maxcolsroundredcostroot) );
   SCIP_CALL( SCIPsetIntParam(origprob, "pricing/masterpricer/maxcolsroundredcost", maxcolsroundredcost) );
   SCIP_CALL( SCIPsetIntParam(origprob, "pricing/masterpricer/maxroundsredcost", maxroundsredcost) );
   SCIP_CALL( SCIPsetIntParam(origprob, "pricing/masterpricer/heurpricingiters", heurpricingiters) );
   SCIP_CALL( SCIPsetCharParam(origprob, "pricing/masterpricer/sorting", sorting) );

#ifdef SCIP_MOREDEBUG
   SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", FALSE) );
#endif

   /* after the solve, the recursion flag can be unset */
   heurdata->inheur = FALSE;

   return SCIP_OKAY;
}

/** destroys and repairs the current best solution. This is the main loop of the ipcolgen heuristic. */
static
SCIP_RETCODE destroyAndRepairSolution(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP_HEURDATA*        heurdata,           /**< the heuristic data */
   SCIP_SOL*             partialsol,         /**< the partial solution that will be updated */
   int                   npricingprobs,      /**< the number of pricing problems */
   SCIP_Bool*            terminate           /**< indicates whether the heuristic must be terminated */
   )
{
   SCIP* scip;
   SCIP_Bool success;
   SCIP_Bool updatesol;
   int i;

   assert(gcg != NULL);
   assert(heurdata != NULL);

   scip = GCGgetMasterprob(gcg);

   SCIPdebugMsg(scip, "Destroy and repair\n");

   heurdata->bestsol = SCIPgetBestSol(scip);

   /* setting the number of relaxation penalties to zero to overwrite the previous penalties */
   heurdata->npenalties = 0;
   heurdata->nmastervars = SCIPgetNVars(scip);
   heurdata->npricingprobs = npricingprobs;

   updatesol = TRUE;

   for( i = 0; i < 3; i++ )
   {
      if( updatesol )
      {
         SCIP_Real initdynamicpen = heurdata->initdynamicpen;

         SCIP_CALL( createPartialSolution(scip, partialsol, heurdata->randnumgen, heurdata->minfixingrate,
               npricingprobs, &heurdata->nfixedvars, &success) );

         /* setting the number of penalties to zero to avoid the reallocation of the memory arrays */
         heurdata->npenalties = 0;

         /* computing the relaxation penalties */
         SCIP_CALL( computeDynamicPenalties(gcg, partialsol, initdynamicpen, &heurdata->penalties,
               &heurdata->penaltiesids, &heurdata->penaltytypes, &heurdata->npenalties, &heurdata->penaltiessize, FALSE) );

         updatesol = FALSE;
      }

      /* if a sufficient number of fixings were not performed, then the heuristic will exit */
      if( !success || !heurdata->npenalties )
      {
         (*terminate) = TRUE;
         break;
      }

      /* pricing in new columns that will be used to find new solutions when solving the repair problem */
      SCIP_CALL( performWeightedPricing(gcg, heurdata, npricingprobs) );

      /* if the best solution is updated, then the partial solution needs to be updated */
      if( heurdata->bestsol != SCIPgetBestSol(scip) )
      {
         SCIP_CALL( SCIPfreeSol(scip, &partialsol) );
         heurdata->bestsol = SCIPgetBestSol(scip);
         SCIP_CALL( SCIPcreateSolCopy(scip, &partialsol, heurdata->bestsol) );
         SCIP_CALL( SCIPunlinkSol(scip, partialsol) );

         updatesol = TRUE;
      }
   }

   /* resetting the number of relaxation penalties to zero */
   heurdata->npenalties = 0;

   return SCIP_OKAY;
}

/** creates the sub-SCIP for the repair problem that is solved to find improving solutions */
static
SCIP_RETCODE createRepairProblem(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP**                repairprob,         /**< the SCIP data structure for the repair problem */
   SCIP_HEURDATA*        heurdata,           /**< the heuristic data */
   SCIP_VAR**            repairprobvars      /**< an array to store the variables of the repair problem */
   )
{
   SCIP_HASHMAP* varmapfw;                   /* mapping of master variables to repair problem variables */
   SCIP_VAR** mastervars;
   int nmastervars;
   int i;

   assert(scip != NULL);

   /* get variable data of the master problem */
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   /* initializing the repair problem */
   SCIP_CALL( SCIPcreate(repairprob) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem((*repairprob)), nmastervars) );

   if( heurdata->uselprows )
   {
      char probname[SCIP_MAXSTRLEN];

      /* copy all plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins((*repairprob)) );

      /* get name of the original problem and add the string "_repair" */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_repair", SCIPgetProbName(scip));

      /* create the repair problem */
      SCIP_CALL( SCIPcreateProb((*repairprob), probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* copy all variables */
      SCIP_CALL( SCIPcopyVars(scip, (*repairprob), varmapfw, NULL, NULL, NULL, 0, TRUE) );
   }
   else
   {
      SCIP_Bool valid;

      valid = FALSE;

      SCIP_CALL( SCIPcopy(scip, (*repairprob), varmapfw, NULL, "repairprob", TRUE, FALSE, FALSE, TRUE, &valid) ); /** @todo check for thread safeness */

      if( heurdata->copycuts )
      {
         /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, (*repairprob), varmapfw, NULL, TRUE, NULL) );
      }

      SCIPdebugMessage("Copying the SCIP instance was %scomplete.\n", valid ? "" : "not ");
   }

   for( i = 0; i < nmastervars; i++ )
      repairprobvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, mastervars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   return SCIP_OKAY;
}

/** set up the repair problem by fixing variables based on the provide solution values */
static
SCIP_RETCODE setupRepairProblem(
   SCIP*                 scip,               /**< SCIP data structure for master problem */
   SCIP*                 repairprob,         /**< SCIP data structure for repair problem */
   SCIP_VAR**            repairprobvars,     /**< the variables of the repair problem */
   SCIP_SOL*             sol,                /**< the source solution */
   SCIP_Bool             uselprows,          /**< should repair problem be created out of the rows in the LP rows? */
   SCIP_Bool             rinsfixing          /**< should a RINS-style fixing be used for the master problem */
   )
{
   SCIP_VAR** mastervars;
   int nmastervars;
   int i;

   /* get variable data of the master problem */
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   /* fix variables based upon sol or the current LP solution in the repair problem */
   for( i = 0; i < nmastervars; i++ )
   {
      SCIP_Real mastersolval;

      mastersolval = SCIPgetSolVal(scip, sol, mastervars[i]);

      if( rinsfixing )
      {
         SCIP_Real lpsolval;

         lpsolval = SCIPgetSolVal(scip, NULL, mastervars[i]);

         /* if the variable takes the same value in the primal solution and the relaxation, then it is fixed */
         if( GCGvarGetBlock(mastervars[i]) >= 0 && !SCIPisFeasEQ(scip, lpsolval, mastersolval) )
         {
            SCIP_CALL( SCIPchgVarUbGlobal(repairprob, repairprobvars[i], mastersolval) );
            SCIP_CALL( SCIPchgVarLbGlobal(repairprob, repairprobvars[i], mastersolval) );
         }
      }
      else
      {
         /* if variable takes a non-zero value in the master solution, then it is fixed to that value */
         if( GCGvarGetBlock(mastervars[i]) >= 0 && !SCIPisZero(scip, mastersolval) )
         {
            SCIP_CALL( SCIPchgVarUbGlobal(repairprob, repairprobvars[i], mastersolval) );
            SCIP_CALL( SCIPchgVarLbGlobal(repairprob, repairprobvars[i], mastersolval) );
         }
      }
   }

   if( uselprows )
   {
      SCIP_ROW** rows;                          /* original scip rows                         */
      int nrows;

      /* get the rows and their number */
      SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

      /* copy all rows to linear constraints */
      for( i = 0; i < nrows; i++ )
      {
         SCIP_CONS* cons;
         SCIP_VAR** consvars;
         SCIP_COL** cols;
         SCIP_Real constant;
         SCIP_Real lhs;
         SCIP_Real rhs;
         SCIP_Real* vals;
         int nnonz;
         int j;

         /* ignore rows that are only locally valid */
         if( SCIProwIsLocal(rows[i]) )
            continue;

         /* get the row's data */
         constant = SCIProwGetConstant(rows[i]);
         lhs = SCIProwGetLhs(rows[i]) - constant;
         rhs = SCIProwGetRhs(rows[i]) - constant;
         vals = SCIProwGetVals(rows[i]);
         nnonz = SCIProwGetNNonz(rows[i]);
         cols = SCIProwGetCols(rows[i]);

         assert(lhs <= rhs);

         /* allocate memory array to be filled with the corresponding repair problem variables */
         SCIP_CALL( SCIPallocBufferArray(repairprob, &consvars, nnonz) );
         for( j = 0; j < nnonz; j++ )
            consvars[j] = repairprobvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

         /* create a new linear constraint and add it to the repair problem */
         SCIP_CALL( SCIPcreateConsLinear(repairprob, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(repairprob, cons) );
         SCIP_CALL( SCIPreleaseCons(repairprob, &cons) );

         /* free temporary memory */
         SCIPfreeBufferArray(repairprob, &consvars);
      }
   }

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by translating the solution of the repair problem */
static
SCIP_RETCODE createNewSol(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP*                 repairprob,         /**< SCIP structure of repair problem */
   SCIP_VAR**            repairprobvars,     /**< the variables of the restricted master problem */
   SCIP_HEUR*            heur,               /**< Restricted Master heuristic structure */
   SCIP_SOL*             repairprobsol,      /**< solution of the restricted master problem */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP* origprob;
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_VAR** mastervars;
   int        nvars;
   int        nmastervars;
   SCIP_Real* repairprobvals;
   SCIP_SOL*  newmastersol;

   assert(gcg != NULL);
   assert(repairprob != NULL);
   assert(repairprobvars != NULL);
   assert(repairprobsol != NULL);

   scip = GCGgetMasterprob(gcg);
   origprob = GCGgetOrigprob(gcg);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(origprob, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(nmastervars == SCIPgetNOrigVars(repairprob));

   SCIP_CALL( SCIPallocBufferArray(scip, &repairprobvals, nmastervars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(repairprob, repairprobsol, nmastervars, repairprobvars, repairprobvals) );

   /* create new solution for the master problem */
   SCIP_CALL( SCIPcreateSol(scip, &newmastersol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newmastersol, nmastervars, mastervars, repairprobvals) );

   /* add solution to the master problem, GCG will translate it and add it to the original problem*/
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintSol(scip, newmastersol, NULL, FALSE) );
   SCIP_CALL( SCIPtrySolFree(scip, &newmastersol, TRUE, TRUE, TRUE, TRUE, TRUE, success) );
#else
   SCIP_CALL( SCIPtrySolFree(scip, &newmastersol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );
#endif
   if( *success == FALSE )
   {
      SCIPdebugMessage("WARNING: original solution feasible, but no solution has been added to master problem.\n");
   }

   SCIPfreeBufferArray(scip, &repairprobvals);

   return SCIP_OKAY;
}


/** solve repair problem */
static
SCIP_RETCODE solveRepairProblem(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP*                 repairprob,         /**< SCIP data structure of the repair problem */
   SCIP_HEUR*            heur,               /**< the ipcolgen heuristic structure */
   SCIP_HEURDATA*        heurdata,           /**< the heuristic data */
   SCIP_VAR**            repairprobvars,     /**< the variables of the repair problem */
   SCIP_Longint          nnodes,             /**< the number of nodes provided to the heuristic */
   SCIP_RESULT*          result              /**< the result from performing the heuristic */
   )
{
   SCIP* origprob;
   SCIP* scip;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_Real cutoff;
   SCIP_Real upperbound;
   SCIP_SOL** repairprobsols;
   int nrepairprobsols;
   SCIP_Bool success;

   int i;

#ifdef NDEBUG
   SCIP_RETCODE retstat;
#endif

   assert(repairprob != NULL);
   assert(heurdata != NULL);

   scip = GCGgetMasterprob(gcg);
   origprob = GCGgetOrigprob(gcg);

   eventhdlr = NULL;
   /* create event handler for LP events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(repairprob, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecIpcolgen, NULL) );
   if( eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* do not abort repair problem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(repairprob, "misc/catchctrlc", FALSE) );

   /* disable output to console */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(repairprob, "display/verblevel", SCIP_VERBLEVEL_FULL) );
#else
   SCIP_CALL( SCIPsetIntParam(repairprob, "display/verblevel", SCIP_VERBLEVEL_NONE) );
   SCIP_CALL( SCIPsetBoolParam(repairprob, "timing/statistictiming", FALSE) );
#endif

   /* set limits for the repair problem */
   SCIP_CALL( SCIPcopyLimits(scip, repairprob) );
   SCIP_CALL( SCIPsetLongintParam(repairprob, "limits/stallnodes", MAX(10, nnodes/10)) );
   SCIP_CALL( SCIPsetLongintParam(repairprob, "limits/nodes", nnodes) );

   /* forbid recursive call of heuristics solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(repairprob, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(repairprob, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(repairprob, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(scip, "estimate") != NULL )
   {
      SCIP_CALL( SCIPsetIntParam(repairprob, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(scip, "inference") != NULL )
   {
      SCIP_CALL( SCIPsetIntParam(repairprob, "branching/inference/priority", INT_MAX/4) );
   }

   /* disable conflict analysis */
   if( !SCIPisParamFixed(repairprob, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(repairprob, "conflict/enable", FALSE) );
   }

   /* setting an objective cutoff with respect to the best solution */
   assert(!SCIPisInfinity(origprob,SCIPgetUpperbound(origprob)));

   upperbound = SCIPgetUpperbound(origprob) - SCIPsumepsilon(origprob);

   if( !SCIPisInfinity(origprob,-1.0*SCIPgetLowerbound(origprob)) )
   {
      cutoff = (1-heurdata->minimprove)*SCIPgetUpperbound(origprob) + heurdata->minimprove*SCIPgetLowerbound(origprob);
   }
   else
   {
      if( SCIPgetUpperbound ( origprob ) >= 0 )
         cutoff = ( 1 - heurdata->minimprove ) * SCIPgetUpperbound ( origprob );
      else
         cutoff = ( 1 + heurdata->minimprove ) * SCIPgetUpperbound ( origprob );
   }
   cutoff = MIN(upperbound, cutoff);
   SCIP_CALL( SCIPsetObjlimit(repairprob, cutoff) );

   /* solve the repair problem */

   /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
    * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
#ifdef NDEBUG
   retstat = SCIPpresolve(repairprob);
   if( retstat != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while presolving subMIP in IPColGen; IPColGen terminated with code <%d>\n",retstat);
   }
#else
   SCIP_CALL( SCIPpresolve(repairprob) );
#endif

   SCIPdebugMessage("presolved the IPColGen repair problem: %d vars, %d cons\n", SCIPgetNVars(repairprob), SCIPgetNConss(repairprob));

   SCIPdebugMessage("solving the IPColGen repair problem: maxnodes=%"SCIP_LONGINT_FORMAT"\n", nnodes);

   /* catching the LP solved events */
   SCIP_CALL( SCIPcatchEvent(repairprob, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

#ifdef NDEBUG
   retstat = SCIPsolve(repairprob);
   if( retstat != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while solving subMIP in IPColGen; IPColGen terminated with code <%d>\n",retstat);
   }
#else
   SCIP_CALL( SCIPsolve(repairprob) );
#endif

   /* drop LP events of sub-SCIP */
   SCIP_CALL( SCIPdropEvent(repairprob, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

   /* print solving statistics of repair problem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(repairprob, NULL) ) );

   SCIPdebugMessage(" -> %d feasible solution(s) found.\n", SCIPgetNSols(repairprob));

   heurdata->usednodes += SCIPgetNNodes(repairprob);

   /* check, whether a solution was found;
    * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
    */
   nrepairprobsols = SCIPgetNSols(repairprob);
   repairprobsols = SCIPgetSols(repairprob);
   success = FALSE;
   for( i = 0; i < nrepairprobsols && !success; ++i )
   {
      SCIP_CALL( createNewSol(gcg, repairprob, repairprobvars, heur, repairprobsols[i], &success) );
   }

   if( success )
      *result = SCIP_FOUNDSOL;

   return SCIP_OKAY;
}


/** checks conditions to determine whether the ipcolgen heuristic should be executed.
 *  Prior to the first call of the ipcolgen heuristic, the average LP objective value is computed to identify the point
 *  when the tailing-off effect starts. This is when the heuristic is first called.
 */
static
SCIP_RETCODE executeHeuristic(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP_HEURDATA*        heurdata,           /**< the heuristic data */
   SCIP_Bool*            execute             /**< returns whether the heuristic should be executed */
   )
{
   SCIP* scip;
   int maxdepth;
   int nblocks;
   int i;

   assert(gcg != NULL);
   assert(heurdata != NULL);
   assert(execute != NULL);

   (*execute) = FALSE;

   scip = GCGgetMasterprob(gcg);

   nblocks = GCGgetNPricingprobs(gcg);

   /* if the previously called node is different to the current node, then the number of pricing iterations must be
    * reset
    */
   if( SCIPgetCurrentNode(scip) != heurdata->prevnode )
   {
      heurdata->prevpricingiter = SCIPpricerGetNCalls(SCIPfindPricer(scip, "gcg"));
      heurdata->numexec = 0;
      heurdata->prevnode = SCIPgetCurrentNode(scip);
   }

   /* storing the number of solutions found before the heuristic is called */
   heurdata->prevnsolsfound = SCIPgetNSols(scip);

   /* a solution must exist before the heuristic can be executed */
   if( SCIPgetBestSol(scip) == NULL )
      return SCIP_OKAY;

   /* imposing the maximum depth setting */
   maxdepth = SCIPheurGetMaxdepth(SCIPfindHeur(scip, HEUR_NAME));
   if( maxdepth >= 0 && SCIPgetDepth(scip) > maxdepth )
      return SCIP_OKAY;

   /* the heuristic will only be called if the current gap is large enough. */
   if( !heurdata->firstexec && SCIPgetGap(scip) < heurdata->mininitialgap )
      return SCIP_OKAY;

   /* the heuristic is only called if a new incumbent solution is found */
   if( heurdata->waitnewsol && SCIPgetBestSol(scip) != NULL && SCIPgetBestSol(scip) == heurdata->lastsol
      && SCIPgetDepth(scip) == 0 )
      return SCIP_OKAY;

   /* only run the heuristic if enough nodes have been processed since the last call */
   if( SCIPgetDepth(scip) > 0 && (SCIPgetNNodes(scip) - heurdata->prevnnodes) < heurdata->nwaitnodes )
      return SCIP_OKAY;

   /* we require at least n blocks of pricing iterations to be performed before the first execution of the heuristic */
   if( SCIPpricerGetNCalls(SCIPfindPricer(scip, "gcg")) < nblocks )
      return SCIP_OKAY;

   /* restricting the number of executions */
   if( heurdata->numexec >= heurdata->callspernode )
      return SCIP_OKAY;

   /* if the heuristic has been executed once, then we require nblocks*0.1 pricing iterations */
   if( heurdata->firstexec )
   {
      if( SCIPpricerGetNCalls(SCIPfindPricer(scip, "gcg")) - heurdata->prevpricingiter > nblocks*0.1 )
         (*execute) = TRUE;

      return SCIP_OKAY;
   }

   /* collecting the previous LP objective values after we have reached the pricer calls threshold */
   if( heurdata->nprevlpobjs == heurdata->prevlpobjssize )
   {
      for( i = 0; i < heurdata->prevlpobjssize - 1; i++ )
      {
         heurdata->prevlpobjs[i] = heurdata->prevlpobjs[i +  1];
      }
   }
   else
      heurdata->nprevlpobjs++;

   heurdata->prevlpobjs[heurdata->nprevlpobjs - 1] = SCIPgetLPObjval(scip);

   /* if enough previous LP objectives have been collected, then average difference is computed and used to assess the
    * execution of the heuristic
    */
   if( heurdata->nprevlpobjs == heurdata->prevlpobjssize )
   {
      SCIP_Real avgobj = 0.0;
      int objcount = 0;

      for( i = 0; i < heurdata->nprevlpobjs; i++ )
      {
         if( !SCIPisZero(scip, heurdata->prevlpobjs[i]) && !SCIPisInfinity(scip, heurdata->prevlpobjs[i]) )
         {
            avgobj += heurdata->prevlpobjs[i];
            objcount++;
         }
      }

      if( objcount < heurdata->nprevlpobjs )
         return SCIP_OKAY;

      avgobj = avgobj/(SCIP_Real) objcount;
      SCIPdebugMsg(scip, "Average objective: %g Current Objective: %g Relative Difference: %g\n", avgobj,
         SCIPgetLPObjval(scip), (avgobj - SCIPgetLPObjval(scip))/REALABS(SCIPgetLPObjval(scip)));

      /* if the average objective and the current object are equal, this could indicate that the column generation has
       * stalled. So we ignore this situation
       */
      if( SCIPisFeasEQ(scip, avgobj, SCIPgetLPObjval(scip)) )
         return SCIP_OKAY;

      /* the threshold for starting the heuristic is if the average difference in the lp objectives is 0.01% of the LP
       * objective value
       */
      if( !SCIPisZero(scip, avgobj)
         && SCIPisLT(scip, (avgobj - SCIPgetLPObjval(scip))/REALABS(SCIPgetLPObjval(scip)), 0.0001) )
      {
         (*execute) = TRUE;
      }
   }

   if( (*execute) && !heurdata->firstexec )
   {
      heurdata->firstexec = TRUE;

      /* storing first executing statistics */
      heurdata->firstcallgap = SCIPgetGap(scip);
      heurdata->firstcallabsgap = SCIPgetPrimalbound(scip) - SCIPgetDualbound(scip);
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#define heurCopyIPcolgen NULL

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeIPcolgen)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   assert(scip == GCGgetDwMasterprob(heurdata->gcg));

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitIPcolgen)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->lastsol = NULL;
   heurdata->inheur = FALSE;
   heurdata->numexec = 0;

   heurdata->prevnode = NULL;

   heurdata->firstexec = FALSE;
   heurdata->prevpricingiter = 0;

   heurdata->prevnsolsfound = 0;
   heurdata->prevnnodes = 0;
   heurdata->nwaitnodes = 0;

   heurdata->prevlpobjssize = PREVLPOBJSIZE;
   heurdata->nprevlpobjs = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->prevlpobjs, heurdata->prevlpobjssize) );

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen,
         DEFAULT_RANDSEED, TRUE) );

   /* initialising the dynamic penalties arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->penalties, DEFAULT_ARRAYSIZE) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->penaltiesids, DEFAULT_ARRAYSIZE) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->penaltytypes, DEFAULT_ARRAYSIZE) );
   heurdata->npenalties = 0;
   heurdata->penaltiessize = DEFAULT_ARRAYSIZE;

   /* initialising the statistics */
   heurdata->firstcallgap = SCIPinfinity(scip);
   heurdata->firstcallabsgap = SCIPinfinity(scip);

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitIPcolgen)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* freeing the relaxation penalties arrays */
   SCIPfreeBlockMemoryArray(scip, &heurdata->penaltytypes, heurdata->penaltiessize);
   SCIPfreeBlockMemoryArray(scip, &heurdata->penaltiesids, heurdata->penaltiessize);
   SCIPfreeBlockMemoryArray(scip, &heurdata->penalties, heurdata->penaltiessize);

   SCIPfreeBlockMemoryArray(scip, &heurdata->prevlpobjs, heurdata->prevlpobjssize);

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolIPcolgen NULL


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolIPcolgen NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecIPcolgen)
{  /*lint --e{715}*/
   GCG* gcg;
   SCIP* origprob;
   SCIP_HEURDATA* heurdata;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Bool discretization;
   SCIP_Bool success;
   SCIP_Bool execute;
   SCIP_Longint nnodes;

   SCIP* repairprob;
   SCIP_VAR** repairprobvars;

   SCIP_SOL* bestsol;
   SCIP_RESULT locresult;
   int iter;
   int solcount;
   int noimprove;
   int maxiter;
   int maxsol;
   int maxnoimprove;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   gcg = heurdata->gcg;
   assert(scip == GCGgetDwMasterprob(gcg));

   /* get original problem */
   origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   maxiter = heurdata->maxiter;
   maxsol = 200;
   maxnoimprove = heurdata->noimproveiter;

   *result = SCIP_DIDNOTRUN;

   /* checking whether the current call is a recursive call to the heuristic. If this heuristic is called from the
    * repair problem, then it is not executed.
    */
   if( heurdata->inheur )
      return SCIP_OKAY;

   /* this heuristic works only for the discretization approach */
   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( !discretization )
      return SCIP_OKAY;

   *result = SCIP_DELAYED;

   /* checking whether the heuristic should be executed */
   SCIP_CALL( executeHeuristic(gcg, heurdata, &execute) );

   if( !execute )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward the heuristic if it succeeded often */
   nnodes = (SCIP_Longint)(nnodes * (SCIPheurGetNSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nnodes -= (SCIP_Longint)(100.0 * SCIPheurGetNCalls(heur));  /* count the setup costs for the sub-MIP as 100 nodes */
   nnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nnodes -= heurdata->usednodes;
   nnodes = MIN(nnodes, heurdata->maxnodes);
   heurdata->nodelimit = nnodes;

   /* check whether we have enough nodes left to call IPColGen */
   if( nnodes < heurdata->minnodes )
      return SCIP_OKAY;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(origprob, "limits/time", &timelimit) );
   if( !SCIPisInfinity(origprob, timelimit) )
      timelimit -= SCIPgetSolvingTime(origprob);
   SCIP_CALL( SCIPgetRealParam(origprob, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(origprob, memorylimit) )
      memorylimit -= SCIPgetMemUsed(origprob)/1048576.0;
   if( timelimit < 10.0 || memorylimit <= 0.0 )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Executing IPColGen ...\n");

   *result = SCIP_DIDNOTFIND;

   /* creating the partial solution for setting the dynamic penalties and fixings in the repair problem */
   bestsol = SCIPgetBestSol(scip);

   /* if the best solution doesn't exist, then a starting solution is created */
   if( bestsol == NULL )
   {
      SCIP_SOL* startsol;

      /* creating the starting solution */
      SCIP_CALL( SCIPcreateSol(scip, &startsol, NULL) );

      success = FALSE;
      SCIP_CALL( createStartSolution(gcg, startsol, &success) );

      if( !success )
      {
         SCIP_CALL( SCIPfreeSol(scip, &startsol) );
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPcreateSolCopy(scip, &heurdata->partialsol, startsol) );
      SCIP_CALL( SCIPunlinkSol(scip, heurdata->partialsol) );

      SCIP_CALL( SCIPfreeSol(scip, &startsol) );
   }
   else
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &heurdata->partialsol, bestsol) );
      SCIP_CALL( SCIPunlinkSol(scip, heurdata->partialsol) );
   }

   iter = 0;
   solcount = 0;
   noimprove = 0;
   while( iter < maxiter && solcount < maxsol && noimprove < maxnoimprove )
   {
      SCIP_Bool terminate;

      /* entering probing mode for the master problem */
      SCIP_CALL( GCGrelaxStartProbing(gcg, NULL) );

      locresult = SCIP_DIDNOTFIND;

      /* performing the destroy and repair of the best solution */
      SCIP_CALL( destroyAndRepairSolution(gcg, heurdata, heurdata->partialsol, GCGgetNPricingprobs(gcg), &terminate) );

      terminate = FALSE;
      if( terminate )
      {
         SCIP_CALL( GCGrelaxEndProbing(gcg) );

         /* the partial solution is no longer needed */
         SCIP_CALL( SCIPfreeSol(scip, &heurdata->partialsol) );

         goto TERMINATE;
      }

      if( heurdata->solveauxproblem )
      {
         success = FALSE;

         /* allocating memory for the repair problem variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &repairprobvars, SCIPgetNVars(scip)) );

         /* creates the repair problem to find improving solutions */
         SCIP_CALL( createRepairProblem(scip, &repairprob, heurdata, repairprobvars) );

         /* set up repair problem by fixing variables based on the partial solution */
         SCIP_CALL( setupRepairProblem(scip, repairprob, repairprobvars, heurdata->partialsol, heurdata->uselprows,
               heurdata->rinsfixing) );
         SCIPdebugMessage("IPColGen repair problem: %d vars, %d cons, success=%u\n", SCIPgetNVars(repairprob),
            SCIPgetNConss(repairprob), success);

         /* solve the repair problem to find improving integer solutions */
         SCIP_CALL( solveRepairProblem(gcg, repairprob, heur, heurdata, repairprobvars, nnodes,
               &locresult) );

         /* free the repair problem */
         SCIPfreeBufferArray(scip, &repairprobvars);
         SCIP_CALL( SCIPfree(&repairprob) );
      }

      /* incrementing the iteration counter */
      iter++;

      /* incrementing the solution counter */
      if( bestsol != SCIPgetBestSol(scip) )
      {
         SCIP_CALL( SCIPfreeSol(scip, &heurdata->partialsol) );
         bestsol = SCIPgetBestSol(scip);
         SCIP_CALL( SCIPcreateSolCopy(scip, &heurdata->partialsol, bestsol) );
         SCIP_CALL( SCIPunlinkSol(scip, heurdata->partialsol) );

         solcount++;
         noimprove = 0;
      }
      else
         noimprove++;

      (*result) = MAX(locresult, (*result));

      /* ending probing mode for the master problem */
      SCIP_CALL( GCGrelaxEndProbing(gcg) );
   }

   /* the partial solution is freed between every iteration of the algorithm.
    * NOTE: this is probably not very efficient and it may be possible to free the solution only when the incumbent is
    * updated
    */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->partialsol) );

   SCIPdebugMessage("Finished IPColGen ...\n");

TERMINATE:
   /* if a new solution is found, then the number of waiting nodes is reset. Otherwise, it is increased */
   if( heurdata->prevnsolsfound < SCIPgetNSols(scip) )
      heurdata->nwaitnodes = 0;
   else
      heurdata->nwaitnodes = (heurdata->nwaitnodes + 1)*10;

   /* updating the last stored solution */
   heurdata->lastsol = SCIPgetBestSol(scip);

   /* incrementing the number of executions */
   heurdata->numexec++;

   /* updating the number of pricing iterations */
   heurdata->prevpricingiter = SCIPpricerGetNCalls(SCIPfindPricer(scip, "gcg"));

   /* updating the number of nodes processed */
   heurdata->prevnnodes = SCIPgetNNodes(scip);

   return SCIP_OKAY;
}




/*
 * primal heuristic specific interface methods
 */

/** creates the IPColGen heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurIPcolgen(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_HEURDATA* heurdata;

   scip = GCGgetDwMasterprob(gcg);

   /* create IPColGen data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heurdata->gcg = gcg;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyIPcolgen, heurFreeIPcolgen, heurInitIPcolgen, heurExitIPcolgen,
         heurInitsolIPcolgen, heurExitsolIPcolgen, heurExecIPcolgen,
         heurdata) );

   /* include the pricing callback plugin */
   SCIP_CALL( GCGpricerIncludePricingcb(gcg, PRICINGCB_NAME, PRICINGCB_DESC, PRICINGCB_PRIORITY,
         NULL, NULL, NULL, NULL, NULL, pricingcbPrepricingIpcolgen, pricingcbPostpricingIpcolgen, NULL) );

   /* add IPColGen parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minfixingrate",
         "minimum percentage of integer variables that have to be fixed",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the repair problem",
         &heurdata->maxnodes,  TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start IPColGen",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of repair problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which IPColGen should at least improve the incumbent  ",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/uselprows",
         "should the repair problem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in repair problem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/dualweight",
         "the weight for the dual values in the pricing problem objective",
         &heurdata->dualweight, FALSE, DEFAULT_DUALWEIGHT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/initdynamicpen",
         "the initial dynamic penalty for the master constraints in the pricing problem objective",
         &heurdata->initdynamicpen, FALSE, DEFAULT_INITDYNAMICPEN, 0.0, 1e+10, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/waitnewsol",
         "should the heuristic wait until a new solution is found before executing",
         &heurdata->waitnewsol, TRUE, DEFAULT_WAITNEWSOL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/solveauxproblem",
         "should an auxiliary problem be solved to find improving solutions",
         &heurdata->solveauxproblem, TRUE, DEFAULT_SOLVEAUXPROB, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/mininitialgap",
         "the minimum initial gap that is necessary before the first call of the heuristic",
         &heurdata->mininitialgap, FALSE, DEFAULT_MININITIALGAP, 0.0, 1e+10, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/callspernode",
         "the maximum number of times that the heuristic is called in each node",
         &heurdata->callspernode, TRUE, DEFAULT_CALLSPERNODE, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxiter",
         "the maximum number of weighted pricing iterations",
         &heurdata->maxiter, TRUE, DEFAULT_MAXITER, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/noimproveiter",
         "the maximum number of weighted pricing iterations without primal improvement",
         &heurdata->noimproveiter, TRUE, DEFAULT_NOIMPROVEITER, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/rinsfixing",
         "should a RINS-style fixing be used for the repair master problem",
         &heurdata->rinsfixing, TRUE, DEFAULT_RINSFIXING, NULL, NULL) );

   return SCIP_OKAY;
}
