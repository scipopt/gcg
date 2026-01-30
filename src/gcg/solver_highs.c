/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file solver_highs.c
 * @brief  highs solver for pricing problems
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define DEBUG_PRICING_ALL_OUTPUT */

#include <math.h>
#include <string.h>
#include <stdint.h>

#include "scip/scip.h"
#include "gcg/gcg.h"
#include "gcg/solver_highs.h"
#include "gcg/pub_solver.h"
#include "gcg/pub_gcgcol.h"

#include "gcg/pricer_gcg.h"
#include "gcg/scip_misc.h"

#include "interfaces/highs_c_api.h"

#define CHECK_ZERO(x) { int _restat_;                                   \
      if( (_restat_ = (x)) != 0 )                                       \
      {                                                                 \
         SCIPerrorMessage("Error in pricing solver: HIGHS returned %d\n", _restat_); \
         retval = SCIP_INVALIDRESULT;                                   \
         goto TERMINATE;                                                \
      }                                                                 \
   }

#define CHECK_SOLVER_RUN(x) { int _restat_;                             \
      if( (_restat_ = (x)) != 0 && (_restat_ = (x) != 1) )              \
      {                                                                 \
         SCIPerrorMessage("Error in pricing solver: HIGHS returned %d\n", _restat_); \
         retval = SCIP_INVALIDRESULT;                                   \
         goto TERMINATE;                                                \
      }                                                                 \
   }

#define SOLVER_NAME                  "highs"
#define SOLVER_DESC                  "highs solver for pricing problems"
#define SOLVER_PRIORITY              10
#define SOLVER_HEURENABLED           FALSE   /**< indicates whether the heuristic solving method of the solver should be enabled */
#define SOLVER_EXACTENABLED          FALSE   /**< indicates whether the exact solving method of the solver should be enabled */

#define DEFAULT_CHECKSOLS            TRUE    /**< should solutions of the pricing MIPs be checked for duplicity? */
#define DEFAULT_THREADS              1       /**< number of threads the HIGHS pricing solver is allowed to use (0: automatic) */
#define DEFAULT_STARTNODELIMIT       1000LL  /**< start node limit for heuristic pricing */
#define DEFAULT_STARTGAPLIMIT        0.2     /**< start gap limit for heuristic pricing */
#define DEFAULT_STARTSOLLIMIT        10LL    /**< start solution limit for heuristic pricing */
#define DEFAULT_NODELIMITFAC         1.25    /**< factor by which to increase node limit for heuristic pricing (1.0: add start limit) */
#define DEFAULT_STALLNODELIMITFAC    1.25    /**< factor by which to increase stalling node limit for heuristic pricing */
#define DEFAULT_GAPLIMITFAC          0.8     /**< factor by which to decrease gap limit for heuristic pricing (1.0: subtract start limit) */
#define DEFAULT_SOLLIMITFAC          1.5     /**< factor by which to increase solution limit for heuristic pricing (1.0: add start limit) */


/** pricing solver data */
struct GCG_SolverData
{
   GCG*                  gcg;                /**< GCG instance */
   SCIP**                pricingprobs;       /**< array storing the SCIP instances for all pricing problems */
   int                   npricingprobs;      /**< number of pricing problems */
   void**                highsptr;            /**< array of pointers to Highs instances */
   int*                  nupdates;           /**< array storing the number of updates for all of the pricing problems */
   SCIP_Longint*         curnodelimit;       /**< current node limit per pricing problem */
   SCIP_Real*            curgaplimit;        /**< current gap limit per pricing problem */
   SCIP_Longint*         cursollimit;        /**< current solution limit per pricing problem */
   /**
    *  information about the basic pricing problem (without potential branching constraints)
    */
   SCIP_VAR***           pricingvars;        /**< array storing the variables of the pricing problems */
   int**                 pricingvartypes;    /**< the variable types of the variables in the highs instances */
   SCIP_CONS***          pricingconss;       /**< array storing the constraints of the pricing problems */
   int*                  npricingvars;       /**< array storing the number of variables of the pricing problems */
   int*                  nbasicpricingconss; /**< array storing the basic number of constraints of the pricing problems */
   /**
    *  parameters
    */
   SCIP_Bool             checksols;          /**< should solutions of the pricing MIPs be checked for duplicity? */
   int                   threads;            /**< number of threads the HIGHS pricing solver is allowed to use (0: automatic) */
   SCIP_Longint          startnodelimit;     /**< start node limit for heuristic pricing */
   SCIP_Real             startgaplimit;      /**< start gap limit for heuristic pricing */
   SCIP_Longint          startsollimit;      /**< start solution limit for heuristic pricing */
   SCIP_Real             nodelimitfac;       /**< factor by which to increase node limit for heuristic pricing (1.0: add start limit) */
   SCIP_Real             gaplimitfac;        /**< factor by which to decrease gap limit for heuristic pricing (1.0: subtract start limit) */
   SCIP_Real             sollimitfac;        /**< factor by which to increase solution limit for heuristic pricing (1.0: add start limit) */

   SCIP_Bool*            ismip;              /**< array storing whether the subproblem is a MIP. */
};

/*
 * local methods
 */

/** creates a HIGHS environment and builds the pricing problem */
static
SCIP_RETCODE buildProblem(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVERDATA*       solverdata,         /**< solver data structure */
   SCIP*                 pricingprob,        /**< pricing problem */
   int                   probnr              /**< problem number */
   )
{
   SCIP* scip;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_VAR* var;

   /* the data to supply to Highs */
   int sense;
   double* varobj;
   double* varlb;
   double* varub;
   int* vartype;
   double* rowupper;
   double* rowlower;
   int* astart;
   int* aindex;
   double* avalue;

   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_VARTYPE type;
   SCIP_RETCODE retval;
   int nconss = 0;
   int nvars = 0;
   int varidx;
   int nconsvars;
   int nnonzeros;
   int idx;
   int c;
   int i;
   int v;

   scip = GCGgetMasterprob(gcg);

   /* open HIGHS environment and create problem */
   solverdata->highsptr[probnr] = Highs_create();
   solverdata->pricingprobs[probnr] = pricingprob;

   retval = SCIP_OKAY;

   /* set parameters */
   CHECK_ZERO( Highs_setDoubleOptionValue(solverdata->highsptr[probnr],
           "mip_rel_gap", 0.0) );
   CHECK_ZERO( Highs_setDoubleOptionValue(solverdata->highsptr[probnr],
           "mip_abs_gap", 0.0) );
   CHECK_ZERO( Highs_setDoubleOptionValue(solverdata->highsptr[probnr],
           "mip_feasibility_tolerance", SCIPfeastol(pricingprob)) );
   CHECK_ZERO( Highs_setIntOptionValue(solverdata->highsptr[probnr],
           "threads", solverdata->threads) );
   CHECK_ZERO( Highs_setBoolOptionValue(solverdata->highsptr[probnr],
         "output_flag", FALSE) );
#ifdef DEBUG_PRICING_ALL_OUTPUT
   CHECK_ZERO( Highs_setBoolOptionValue(solverdata->highsptr[probnr],
           "output_flag", TRUE) );
#endif

   /* set objective sense */
   assert(SCIPgetObjsense(pricingprob) == SCIP_OBJSENSE_MINIMIZE);
   sense = 1; // minimisation

   conss = SCIPgetOrigConss(pricingprob);
   nconss = SCIPgetNOrigConss(pricingprob);
   vars = SCIPgetOrigVars(pricingprob);
   nvars = SCIPgetNOrigVars(pricingprob);

   /* arrays for storing the basic constraints and variables */
   solverdata->npricingvars[probnr] = nvars;
   solverdata->nbasicpricingconss[probnr] = nconss;

   SCIP_CALL( SCIPallocMemoryArray(scip, &solverdata->pricingvars[probnr], nvars) ); /*lint !e866*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &solverdata->pricingvartypes[probnr], nvars) ); /*lint !e866*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &solverdata->pricingconss[probnr], nconss) ); /*lint !e866*/

   /* temporary memory for storing all data about the variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &varobj, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vartype, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varlb, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varub, nvars) );

   /* temporary memory for storing data about the constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &rowlower, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowupper, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &astart, nconss + 1) );

   /* collect information about variables: bounds, objective function, name, type */
   for( i = 0; i < nvars; i++ )
   {
      var = vars[i];
      varidx = SCIPvarGetIndex(var);
      assert(0 <= varidx);
      assert(varidx < nvars);
      solverdata->pricingvars[probnr][varidx] = var;
      SCIP_CALL( SCIPcaptureVar(pricingprob, var) );

      varlb[varidx] = SCIPvarGetLbLocal(var);
      varub[varidx] = SCIPvarGetUbLocal(var);
      varobj[varidx] = SCIPvarGetObj(var);

      type = SCIPvarGetType(var);

      switch( type )
      {
      case SCIP_VARTYPE_CONTINUOUS:
         vartype[varidx] = 0;
         break;
      case SCIP_VARTYPE_BINARY:
      case SCIP_VARTYPE_INTEGER:
         vartype[varidx] = 1;
         break;
      default:
         SCIPerrorMessage("invalid variable type\n");
         return SCIP_INVALIDDATA;
      }

      /* storing whether the problem is a MIP */
      if( vartype[varidx] > 0 )
         solverdata->ismip[probnr] = TRUE;

      /* storing the pricing variable type */
      solverdata->pricingvartypes[probnr][varidx] = vartype[varidx];
   }

   /* collect right hand sides and ranges of the constraints, count total number of nonzeros */
   nnonzeros = 0;
   for( c = 0; c < nconss; ++c )
   {
      solverdata->pricingconss[probnr][c] = conss[c];
      SCIP_CALL( SCIPcaptureCons(pricingprob, conss[c]) );

      nnonzeros += GCGconsGetNVars(scip, conss[c]);
      lhs = GCGconsGetLhs(pricingprob, conss[c]);
      rhs = GCGconsGetRhs(pricingprob, conss[c]);

      if( SCIPisInfinity(scip, -lhs) )
         rowlower[c] = -Highs_getInfinity(solverdata->highsptr[probnr]);
      else
         rowlower[c] = lhs;

      if( SCIPisInfinity(scip, rhs) )
         rowupper[c] = Highs_getInfinity(solverdata->highsptr[probnr]);
      else
         rowupper[c] = rhs;
   }

   /* temporary memory for storing data about coefficients in the constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars) );

   /* temporary memory for storing the CSR matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &aindex, nnonzeros) );
   SCIP_CALL( SCIPallocBufferArray(scip, &avalue, nnonzeros) );

   /* collect nonzeros */
   for( c = 0, idx = 0; c < nconss; ++c )
   {
      /* storing the start index for the constraint */
      astart[c] = idx;

      nconsvars = GCGconsGetNVars(scip, conss[c]);

      /* while this should not happen, sometimes a subproblem is detected that
       * has empty constraints.
       */
      if (nconsvars == 0)
         continue;

      SCIP_CALL( GCGconsGetVals(pricingprob, conss[c], consvals, nvars) );
      SCIP_CALL( GCGconsGetVars(pricingprob, conss[c], consvars, nvars) );

      /* get coefficients */
      for( v = 0; v < nconsvars; ++v )
      {
         aindex[idx] = SCIPvarGetIndex(consvars[v]);
         avalue[idx] = (double)consvals[v];
         idx++;
      }
   }
   assert(idx == nnonzeros);
   astart[nconss] = nnonzeros;

   /* passing the data to Highs to build the MIP model
    * NOTE: the format is row-wise, where HiGHS typically expects col-wise. The value for row-wise is 2,
    * col-wise is 1*/
   CHECK_ZERO( Highs_passMip(solverdata->highsptr[probnr], nvars, nconss, nnonzeros, 2, sense, 0, varobj, varlb,
         varub, rowlower, rowupper, astart, aindex, avalue, vartype) );

#ifdef WRITEPROBLEMS
   {
      char filename[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "highs-%s.lp", SCIPgetProbName(pricingprob));
      SCIPinfoMessage(pricingprob, NULL, "print pricing problem to %s\n", filename);
      CHECK_ZERO( Highs_writeModel(solverdata->highsptr[probnr], filename) );
   }
#endif

 TERMINATE:
   /* free temporary memory */
   if( avalue != NULL )
   {
      SCIPfreeBufferArray(scip, &avalue);
      SCIPfreeBufferArray(scip, &aindex);

      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);

      SCIPfreeBufferArray(scip, &astart);
      SCIPfreeBufferArray(scip, &rowupper);
      SCIPfreeBufferArray(scip, &rowlower);

      SCIPfreeBufferArray(scip, &varub);
      SCIPfreeBufferArray(scip, &varlb);
      SCIPfreeBufferArray(scip, &vartype);
      SCIPfreeBufferArray(scip, &varobj);
   }

   return retval;
}


/** updates bounds and objective coefficients of variables in the given pricing problem */
static
SCIP_RETCODE updateVars(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_SOLVERDATA*       solverdata,         /**< solver data structure */
   SCIP*                 pricingprob,        /**< pricing problem */
   int                   probnr,             /**< problem number */
   SCIP_Bool             varobjschanged,     /**< have the objective coefficients changed? */
   SCIP_Bool             varbndschanged      /**< have the lower and upper bounds changed? */
   )
{
   SCIP_VAR** vars;
   double* varobj;
   double* collower;
   double* colupper;
   int* objidx;
   int* updatevaridx;
   int nvars;
   int npricingvars;

   SCIP_RETCODE retval;
   int i;

   vars = SCIPgetOrigVars(pricingprob);
   nvars = SCIPgetNOrigVars(pricingprob);
   npricingvars = solverdata->npricingvars[probnr];

   assert(npricingvars == nvars);
   assert(npricingvars == Highs_getNumCol(solverdata->highsptr[probnr]));

   retval = SCIP_OKAY;

   if( varobjschanged )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &objidx, npricingvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varobj, npricingvars) );
   }

   if( varbndschanged )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &updatevaridx, npricingvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &collower, npricingvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &colupper, npricingvars) );
   }

   /* get new bounds and objective coefficients of variables */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_VAR* origvar;
      SCIP_VAR* var;
      int varidx;

      origvar = vars[i];
      varidx = SCIPvarGetIndex(origvar);
      assert(0 <= varidx);
      assert(varidx < npricingvars);

      if( SCIPgetStage(pricingprob) >= SCIP_STAGE_TRANSFORMED )
         var = SCIPvarGetTransVar(origvar);
      else
         var = origvar;

      if( varbndschanged )
      {
         updatevaridx[(size_t)varidx] = varidx;
         collower[(size_t)varidx] = (double) SCIPvarGetLbGlobal(var);
         colupper[(size_t)varidx] = (double) SCIPvarGetUbGlobal(var);
      }

      if( varobjschanged )
      {
         objidx[varidx] = varidx;
         varobj[varidx] = SCIPvarGetObj(origvar);
      }
   }

   /* update bounds and objective coefficient of basic variables */
   if( varbndschanged )
   {
      CHECK_ZERO( Highs_changeColsBoundsBySet(solverdata->highsptr[probnr], nvars, updatevaridx, collower, colupper) );
   }
   if( varobjschanged )
   {
      CHECK_ZERO( Highs_changeColsCostBySet(solverdata->highsptr[probnr], nvars, objidx, varobj) );
   }

TERMINATE:
   if( varbndschanged )
   {
      SCIPfreeBufferArray(scip, &colupper);
      SCIPfreeBufferArray(scip, &collower);
      SCIPfreeBufferArray(scip, &updatevaridx);
   }
   if( varobjschanged)
   {
      SCIPfreeBufferArray(scip, &varobj);
      SCIPfreeBufferArray(scip, &objidx);
   }

   return retval;
}

/** updates branching constraints in the given pricing problem */
static
SCIP_RETCODE updateBranchingConss(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_SOLVERDATA*       solverdata,         /**< solver data structure */
   SCIP*                 pricingprob,        /**< pricing problem */
   int                   probnr              /**< problem number */
   )
{
   SCIP_CONS** conss;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;

   /* data for the rows */
   double* newrowlower;
   double* newrowupper;
   int* newstart;
   int* newindex;
   double* newvalue;

   int nconss;
   int nbasicpricingconss;
   int nhighsrows;
   int nnewconss = 0;
   int nnonzeros;
   int nvars;

   SCIP_RETCODE retval;
   int idx;
   int c;

   conss = SCIPgetOrigConss(pricingprob);
   nconss = SCIPgetNOrigConss(pricingprob);
   nbasicpricingconss = solverdata->nbasicpricingconss[probnr];

   nvars = SCIPgetNOrigVars(pricingprob);

   retval = SCIP_OKAY;

   nhighsrows = Highs_getNumRow(solverdata->highsptr[probnr]);

   if( nbasicpricingconss < nhighsrows )
   {
      CHECK_ZERO( Highs_deleteRowsByRange(solverdata->highsptr[probnr], nbasicpricingconss, nhighsrows - 1) );
   }

   nnewconss = nconss - nbasicpricingconss;

   if( nnewconss == 0 )
      return retval;

   /* temporary arrays for storing data about new constraints */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &newrowlower, nnewconss) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &newrowupper, nnewconss) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &newstart, nnewconss + 1) );

   /* get information about new constraints */
   nnonzeros = 0;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_Real lhs;
      SCIP_Real rhs;
      int nconsvars;
      int considx;

      /* we assume that nothing changed about the basic constraints */
      if( c < nbasicpricingconss )
      {
         assert(conss[c] == solverdata->pricingconss[probnr][c]);
         continue;
      }

      considx = c - nbasicpricingconss;
      assert(considx >= 0);

      nconsvars = GCGconsGetNVars(scip, conss[c]);
      lhs = GCGconsGetLhs(pricingprob, conss[c]);
      rhs = GCGconsGetRhs(pricingprob, conss[c]);

      if( SCIPisInfinity(scip, -lhs) )
         newrowlower[considx] = -Highs_getInfinity(solverdata->highsptr[probnr]);
      else
         newrowlower[considx] = lhs;

      if( SCIPisInfinity(scip, rhs) )
         newrowupper[considx] = Highs_getInfinity(solverdata->highsptr[probnr]);
      else
         newrowupper[considx] = rhs;

      nnonzeros += nconsvars;
   }

   /* temporary arrays for getting variables and coefficients in new constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars) );

   /* temporary arrays for storing data about new constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &newindex, nnonzeros) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newvalue, nnonzeros) );

   /* collect coefficients in new constriants */
   for( c = 0, idx = 0; c < nconss; ++c )
   {
      int considx;
      int nconsvars;
      int v;

      /* we assume that nothing changed about the basic constraints */
      if( c < nbasicpricingconss )
      {
         assert(conss[c] == solverdata->pricingconss[probnr][c]);
         continue;
      }

      considx = c - nbasicpricingconss;
      assert(considx >= 0);

      newstart[considx] = idx;
      nconsvars = GCGconsGetNVars(scip, conss[c]);
      SCIP_CALL( GCGconsGetVars(pricingprob, conss[c], consvars, nvars) );
      SCIP_CALL( GCGconsGetVals(pricingprob, conss[c], consvals, nvars) );

      /* get coefficients */
      for( v = 0; v < nconsvars; ++v )
      {
         newindex[idx] = SCIPvarGetIndex(consvars[v]);
         newvalue[idx] = (double)consvals[v];
         idx++;
      }
   }
   assert(idx == nnonzeros);
   newstart[nconss - nbasicpricingconss] = idx;

   /* add new constraints */
   CHECK_ZERO( Highs_addRows(solverdata->highsptr[probnr], nnewconss, newrowlower, newrowupper, nnonzeros, newstart,
         newindex, newvalue) );

TERMINATE:
   if( nnewconss > 0 )
   {
      SCIPfreeBufferArray(scip, &newvalue);
      SCIPfreeBufferArray(scip, &newindex);

      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);

      SCIPfreeBufferArray(scip, &newstart);
      SCIPfreeBufferArray(scip, &newrowupper);
      SCIPfreeBufferArray(scip, &newrowlower);
   }

   return retval;
}


/** solves the pricing problem with HIGHS */
static
SCIP_RETCODE solveHighs(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVERDATA*       solverdata,         /**< solver data structure */
   SCIP*                 pricingprob,        /**< pricing problem */
   int                   probnr,             /**< problem number */
   SCIP_Real             dualsolconv,        /**< dual solution value of the corresponding convexity constraint */
   SCIP_Real*            lowerbound,         /**< pointer to store lower bound */
   int*                  ncols,              /**< pointer to store number of columns */
   GCG_PRICINGSTATUS*    status              /**< pointer to store the pricing status */
   )
{ /*lint -e715*/
   SCIP* scip;
   GCG_COL* col;
   SCIP_RETCODE retval;
   SCIP_Bool predisabled = FALSE;
   double* highssolvals;
   double upperbound;
   double* primsol = NULL;
   int numcols;
   int modelstatus;
   int runretval;

   SCIP_SOL* sol;
   SCIP_Bool feasible;

   scip = GCGgetMasterprob(gcg);
   *ncols = 0;
   *status = GCG_PRICINGSTATUS_UNKNOWN;
   upperbound = SCIPinfinity(pricingprob);

   retval = SCIP_OKAY;

   numcols = Highs_getNumCol(solverdata->highsptr[probnr]);
   assert(numcols == SCIPgetNOrigVars(pricingprob));

   SCIP_CALL( SCIPallocBufferArray(scip, &highssolvals, numcols) );
 SOLVEAGAIN:
   /* the optimization call */
   runretval = Highs_run(solverdata->highsptr[probnr]);
   CHECK_SOLVER_RUN(runretval);

   /* get model status from Highs */
   modelstatus = Highs_getModelStatus(solverdata->highsptr[probnr]);

   /* handle HIGHS solution status */
   switch( modelstatus )
   {
      /* pricing problem was solved to optimality.
       * NOTE: the optimal status is also returned when a gap limit is reached. As such, the heuristic pricing with
       * gap limits needs to be handled here.
       */
      case 7: /* HighsModelStatus::kOptimal */
      {
         double mipgap;
         assert(runretval == 0);

         /* getting the MIP gap for the solution */
         CHECK_ZERO( Highs_getDoubleInfoValue(solverdata->highsptr[probnr], "mip_gap", &mipgap) );

         /* if the gap is 0, then we have an optimal solution */
         if( !solverdata->ismip[probnr] || SCIPisFeasZero(scip, mipgap) )
         {
            *status = GCG_PRICINGSTATUS_OPTIMAL;
         }
         else
         {
            if( solverdata->gaplimitfac < 1.0 )
               solverdata->curgaplimit[probnr] *= solverdata->gaplimitfac;
            else
               solverdata->curgaplimit[probnr] = MAX(solverdata->curgaplimit[probnr] - solverdata->startgaplimit, 0.0);
            SCIPdebugMessage("   -> gap limit reached, decreasing to %g\n", solverdata->curgaplimit[probnr]);
            *status = GCG_PRICINGSTATUS_SOLVERLIMIT;
         }

         upperbound = Highs_getObjectiveValue(solverdata->highsptr[probnr]);
         break;
      }

      /* pricing problem was proven to be infeasible */
      case 8: /* HighsModelStatus::kInfeasible */
         assert(runretval == 0);
         *status = GCG_PRICINGSTATUS_INFEASIBLE;
         break;

      /* pricing problem is possibly unbounded */
      case 9:  /* HighsModelStatus::kUnboundedOrInfeasible */
      case 10: /* HighsModelStatus::kUnbounded */
      {
         int highsretval;
         SCIP_Bool hasprimalray;

         assert(runretval == 0);

         SCIP_CALL( SCIPallocBufferArray(scip, &primsol, numcols) );

         CHECK_ZERO( Highs_getSolution(solverdata->highsptr[probnr], primsol, NULL, NULL, NULL) );

         highsretval = Highs_getPrimalRay(solverdata->highsptr[probnr], &hasprimalray, highssolvals);

         if( highsretval != 0 )
         {
            assert(!predisabled);

            SCIPdebugMessage("   -> disable presolving in HIGHS and only solve the LP relaxation to get primal ray\n");

            CHECK_ZERO( Highs_setStringOptionValue(solverdata->highsptr[probnr], "presolve", "off") );
            CHECK_ZERO( Highs_setBoolOptionValue(solverdata->highsptr[probnr], "solve_relaxation", TRUE) );

            predisabled = TRUE;

            SCIPfreeBufferArray(scip, &primsol);

            goto SOLVEAGAIN;
         }
         else
         {
            CHECK_ZERO( highsretval );
         }

         /* since the primal ray is found by solving an LP, then it is necessary to round the integer variables in
          * the direction of the objective function.
          */
         for (int i = 0; i < numcols; i++)
         {
            /* checking for integer or binary variables. These have the variable type 1 */
            if( solverdata->pricingvartypes[probnr][i] == 1 && !SCIPisIntegral(scip, highssolvals[i]) )
            {
               if( SCIPgetObjsense(pricingprob) == SCIP_OBJSENSE_MINIMIZE )
                  highssolvals[i] = SCIPisPositive(pricingprob, SCIPvarGetObj(solverdata->pricingvars[probnr][i])) ?
                     SCIPfloor(pricingprob, highssolvals[i]) : SCIPceil(pricingprob, highssolvals[i]);
               else
                  highssolvals[i] = SCIPisNegative(pricingprob, SCIPvarGetObj(solverdata->pricingvars[probnr][i])) ?
                     SCIPfloor(pricingprob, highssolvals[i]) : SCIPceil(pricingprob, highssolvals[i]);
            }
            assert(SCIPisIntegral(scip, highssolvals[i]));
         }


         SCIP_CALL( GCGcreateGcgCol(gcg, pricingprob, &col, probnr, solverdata->pricingvars[probnr], highssolvals, numcols, TRUE, SCIPinfinity(pricingprob)) );
         SCIP_CALL( GCGpricerAddCol(gcg, col) );
         ++(*ncols);

         *status = GCG_PRICINGSTATUS_UNBOUNDED;

         SCIPfreeBufferArray(scip, &primsol);

         goto TERMINATE;
      }

      /* a heuristic pricing limit was reached and may be increased in the next round */
      /* The iteration limit model status is used for indicating a node or a solution is reached. It must be checked
       * to identify what limit actually triggered the termination. Since HiGHS doesn't record the number of
       * solutions found, we can only determine the limit from the number of processed nodes.
       */
      case 14: /* HighsModelStatus::kIterationLimit */
      {
         int64_t nodecount;

         assert(runretval == 1);

         /* getting the node and solution count */
         CHECK_ZERO( Highs_getInt64InfoValue(solverdata->highsptr[probnr], "mip_node_count", &nodecount) );

         if( nodecount > solverdata->curnodelimit[probnr] )
         {
            int solstatus;
            SCIP_Longint nodelimitfac;

            /* checking whether a solution exists. If not, then we don't know the current solution status */
            CHECK_ZERO( Highs_getIntInfoValue(solverdata->highsptr[probnr], "primal_solution_status", &solstatus) );

            if( solstatus < 2 ) /* either infeasible or no solution */
            {
               *status = GCG_PRICINGSTATUS_UNKNOWN;
               break;
            }

            nodelimitfac = solverdata->nodelimitfac;

            if( nodelimitfac > 1.0 )
               solverdata->curnodelimit[probnr] = (SCIP_Longint) (solverdata->curnodelimit[probnr] * solverdata->nodelimitfac);
            else
               solverdata->curnodelimit[probnr] += solverdata->startnodelimit;
            SCIPdebugMessage("   -> node limit reached, increasing to %"SCIP_LONGINT_FORMAT"\n", solverdata->curnodelimit[probnr]);
         }
         else
         {
            if( solverdata->sollimitfac > 1.0 )
               solverdata->cursollimit[probnr] = (SCIP_Longint) (solverdata->cursollimit[probnr] * solverdata->sollimitfac);
            else
               solverdata->cursollimit[probnr] += solverdata->startsollimit;
            SCIPdebugMessage("   -> solution limit reached, increasing to %"SCIP_LONGINT_FORMAT"\n", solverdata->cursollimit[probnr]);
         }
         *status = GCG_PRICINGSTATUS_SOLVERLIMIT;
         upperbound = Highs_getObjectiveValue(solverdata->highsptr[probnr]);
         break;
      }

      /* a limit is reached, but not handled by GCG. We check the solution status to determine whether we can use the
       * upper bound
      */
      case 13: /* HighsModelStatus::kTimeLimit */
      {
         int solstatus;

         assert(runretval == 1);

         /* checking whether a solution exists. If not, then we don't know the current solution status */
         CHECK_ZERO( Highs_getIntInfoValue(solverdata->highsptr[probnr], "primal_solution_status", &solstatus) );

         if( solstatus == 2 ) /* the solution is feasible */
         {
            upperbound = Highs_getObjectiveValue(solverdata->highsptr[probnr]);
         }
         *status = GCG_PRICINGSTATUS_UNKNOWN;
         break;
      }

      case 15: /* HighsModelStatus::kUnknown */
      {
         int solstatus;

         assert(runretval == 1);

         /* checking whether a solution exists. If not, then we don't know the current solution status */
         CHECK_ZERO( Highs_getIntInfoValue(solverdata->highsptr[probnr], "primal_solution_status", &solstatus) );

         if( solstatus == 2 ) /* the solution is feasible */
         {
            upperbound = Highs_getObjectiveValue(solverdata->highsptr[probnr]);
         }
         *status = GCG_PRICINGSTATUS_UNKNOWN;
         break;
      }

      default:
         *status = GCG_PRICINGSTATUS_UNKNOWN;
         goto TERMINATE;
   }

   /* if the pricing problem is a MIP, then we can collect the dual bound. Otherwise, if it is an LP, we expect that
    * the problem is solved to optimality, so the upperbound == lowerbound
   */
   if( solverdata->ismip[probnr] )
   {
      CHECK_ZERO( Highs_getDoubleInfoValue(solverdata->highsptr[probnr], "mip_dual_bound", lowerbound) );
   }
   else
      *lowerbound = upperbound;

   /* checking whether the lower bound is returned from Highs as inf */
   if (isinf(*lowerbound))
      *lowerbound = -SCIPinfinity(scip);

   assert(SCIPisFeasLE(scip, *lowerbound, upperbound));

   SCIPdebugMessage("   -> pricing problem %d solved: status=%d, modelstatus=%d, lowerbound=%g, upperbound=%g\n",
       probnr, *status, modelstatus, *lowerbound, upperbound);

   assert(SCIPisFeasEQ(scip, *lowerbound, upperbound) || *status != GCG_PRICINGSTATUS_OPTIMAL);

   /* extracting the best solution and checked if it has a negative reduced cost */
   CHECK_ZERO( Highs_getSolution(solverdata->highsptr[probnr], highssolvals, NULL, NULL, NULL) );

   /* creating a solution from the column generated from solving Highs */
   SCIP_CALL( SCIPcreateOrigSol(pricingprob, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVals(pricingprob, sol, numcols, solverdata->pricingvars[probnr], highssolvals) );

   feasible = FALSE;

   /* check whether the solution is feasible */
   if( !solverdata->checksols )
   {
      SCIP_CALL( SCIPcheckSolOrig(pricingprob, sol, &feasible, FALSE, FALSE) );

      /* if the optimal solution is not feasible, we return GCG_PRICINGSTATUS_UNKNOWN as status */
      if( !feasible )
         *status = GCG_PRICINGSTATUS_UNKNOWN;
   }
   else
      feasible = TRUE;

   if( feasible )
   {
      SCIP_CALL( GCGcreateGcgColFromSol(gcg, pricingprob, NULL, NULL, &col, probnr, sol, FALSE, SCIPinfinity(pricingprob)) );
      SCIP_CALL( GCGpricerAddCol(gcg, col) );
      ++(*ncols);
   }

   SCIP_CALL( SCIPfreeSol(pricingprob, &sol) );

   assert(*status != GCG_PRICINGSTATUS_OPTIMAL || *ncols > 0);
 TERMINATE:
   if( predisabled )
   {
      CHECK_ZERO( Highs_setStringOptionValue(solverdata->highsptr[probnr], "presolve", "on") );
      CHECK_ZERO( Highs_setBoolOptionValue(solverdata->highsptr[probnr], "solve_relaxation", FALSE) );
   }

   if( primsol != NULL )
      SCIPfreeBufferArray(scip, &primsol);

   SCIPfreeBufferArray(scip, &highssolvals);

   return retval;
}


/** destructor of pricing solver to free user data (called when SCIP is exiting) */
static
GCG_DECL_SOLVERFREE(solverFreeHighs)
{
   GCG_SOLVERDATA* solverdata;

   assert(gcg != NULL);
   assert(solver != NULL);

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   SCIPfreeMemory(GCGgetDwMasterprob(gcg), &solverdata);
   GCGsolverSetData(solver, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of pricing solver (called when branch and bound process is about to begin) */
static
GCG_DECL_SOLVERINITSOL(solverInitsolHighs)
{
   SCIP* scip;
   GCG_SOLVERDATA* solverdata;
   int npricingprobs;

   int i;

   assert(gcg != NULL);
   assert(solver != NULL);

   scip = GCGgetMasterprob(gcg);
   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   solverdata->npricingprobs = GCGgetNPricingprobs(gcg);
   npricingprobs = solverdata->npricingprobs;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->highsptr), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->nupdates), npricingprobs) );
   BMSclearMemoryArray(solverdata->nupdates, npricingprobs);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->pricingprobs), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->pricingvars), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->pricingvartypes), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->pricingconss), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->npricingvars), npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solverdata->nbasicpricingconss), npricingprobs) );
   BMSclearMemoryArray(solverdata->npricingvars, npricingprobs);
   BMSclearMemoryArray(solverdata->nbasicpricingconss, npricingprobs);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solverdata->curnodelimit, npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solverdata->curgaplimit, npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solverdata->cursollimit, npricingprobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solverdata->ismip, npricingprobs) );

   BMSclearMemoryArray(solverdata->ismip, npricingprobs);

   for( i = 0; i < npricingprobs; ++i )
   {
      if( GCGisPricingprobRelevant(gcg, i) )
      {
         SCIP_CALL( buildProblem(gcg, solverdata, GCGgetPricingprob(gcg, i), i) );
      }

      solverdata->curnodelimit[i] = solverdata->startnodelimit;
      solverdata->curgaplimit[i] = solverdata->startgaplimit;
      solverdata->cursollimit[i] = solverdata->startsollimit;
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of pricing solver (called before branch and bound process data is freed) */
static
GCG_DECL_SOLVEREXITSOL(solverExitsolHighs)
{
   SCIP* scip;
   GCG_SOLVERDATA* solverdata;
   SCIP_RETCODE retval;
   int npricingprobs;
   int i;
   int j;

   assert(gcg != NULL);
   assert(solver != NULL);

   scip = GCGgetMasterprob(gcg);
   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   retval = SCIP_OKAY;

   npricingprobs = GCGgetNPricingprobs(gcg);

   /* free pricing problems */
   for( i = 0; i < npricingprobs; ++i )
   {
      if( GCGisPricingprobRelevant(gcg, i) )
      {
         /* free the HiGHS instance */
         Highs_destroy(solverdata->highsptr[i]);

         if( solverdata->nbasicpricingconss[i] > 0 )
         {
            /* release stored constraints */
            for( j = 0; j < solverdata->nbasicpricingconss[i]; ++j )
            {
               SCIP_CALL( SCIPreleaseCons(solverdata->pricingprobs[i], &solverdata->pricingconss[i][j]) );
            }
            SCIPfreeMemoryArray(scip, &(solverdata->pricingconss[i]));
         }

         if( solverdata->npricingvars[i] > 0 )
         {
            /* release stored constraints */
            for( j = 0; j < solverdata->npricingvars[i]; ++j )
            {
               SCIP_CALL( SCIPreleaseVar(solverdata->pricingprobs[i], &solverdata->pricingvars[i][j]) );
            }
            SCIPfreeMemoryArray(scip, &(solverdata->pricingvartypes[i]));
            SCIPfreeMemoryArray(scip, &(solverdata->pricingvars[i]));
         }
      }
   }

   SCIPfreeBlockMemoryArray(scip, &solverdata->ismip, solverdata->npricingprobs);

   SCIPfreeBlockMemoryArray(scip, &solverdata->cursollimit, solverdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &solverdata->curgaplimit, solverdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &solverdata->curnodelimit, solverdata->npricingprobs);

   SCIPfreeBlockMemoryArray(scip, &(solverdata->nbasicpricingconss), solverdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(solverdata->npricingvars), solverdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(solverdata->pricingconss), solverdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(solverdata->pricingvartypes), solverdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(solverdata->pricingvars), solverdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(solverdata->pricingprobs), solverdata->npricingprobs);

   SCIPfreeBlockMemoryArray(scip, &(solverdata->nupdates), solverdata->npricingprobs);
   SCIPfreeBlockMemoryArray(scip, &(solverdata->highsptr), solverdata->npricingprobs);

   return retval;
}

#define solverInitHighs NULL
#define solverExitHighs NULL

/** update method for pricing solver, used to update solver specific pricing problem data */
static
GCG_DECL_SOLVERUPDATE(solverUpdateHighs)
{
   GCG_SOLVERDATA* solverdata;

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   SCIPdebugMessage("HIGHS solver -- update data for problem %d: varobjschanged = %u, varbndschanged = %u, consschanged = %u\n",
      probnr, varobjschanged, varbndschanged, consschanged);

   /* update pricing problem information */
   SCIP_CALL( updateVars(GCGgetMasterprob(gcg), solverdata, pricingprob, probnr, varobjschanged, varbndschanged) );
   if( consschanged )
   {
      SCIP_CALL( updateBranchingConss(GCGgetMasterprob(gcg), solverdata, pricingprob, probnr) );
   }

   solverdata->curnodelimit[probnr] = solverdata->startnodelimit;
   solverdata->curgaplimit[probnr] = solverdata->startgaplimit;
   solverdata->cursollimit[probnr] = solverdata->startsollimit;

#ifdef WRITEPROBLEMS
   /* Print the pricing problem after updating:
    *  * after checking variable bounds, because they change in particular when a new generic branching subproblem is considered
    *  * but not after adding new branching constraints, since objectives will be set afterwards before solving
    */
   if( varbndschanged && !consschanged )
   {
      char filename[SCIP_MAXSTRLEN];

      ++(solverdata->nupdates[probnr]);

      (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "highs-%s-%d-%d.lp", SCIPgetProbName(pricingprob), SCIPgetNNodes(scip), solverdata->nupdates[probnr]);
      SCIPinfoMessage(pricingprob, NULL, "print pricing problem to %s\n", filename);
      CHECK_ZERO( Highs_writeModel(solverdata->highsptr[probnr], filename) );
   }
#endif

   return SCIP_OKAY;
}

/** heuristic solving method of HIGHS solver */
static
GCG_DECL_SOLVERSOLVEHEUR(solverSolveHeurHighs)
{
   GCG_SOLVERDATA* solverdata;
   int ncols;
   SCIP_RETCODE retval;

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   SCIPdebugMessage("calling heuristic pricing with HIGHS for pricing problem %d\n", probnr);

   retval = SCIP_OKAY;

   /* set heuristic limits */
   CHECK_ZERO( Highs_setIntOptionValue(solverdata->highsptr[probnr], "mip_max_nodes",
         (int) solverdata->curnodelimit[probnr]) );
   CHECK_ZERO( Highs_setDoubleOptionValue(solverdata->highsptr[probnr], "mip_rel_gap",
         (double) solverdata->curgaplimit[probnr]) );
   CHECK_ZERO( Highs_setIntOptionValue(solverdata->highsptr[probnr], "mip_max_improving_sols",
         (int) solverdata->cursollimit[probnr]) );

   /* solve the pricing problem and evaluate solution */
   SCIP_CALL( solveHighs(gcg, solverdata, pricingprob, probnr, dualsolconv, lowerbound, &ncols, status) );
   assert(*status != GCG_PRICINGSTATUS_OPTIMAL || ncols > 0);

 TERMINATE:
   return retval;
}

/** solving method for pricing solver which solves the pricing problem to optimality */
static
GCG_DECL_SOLVERSOLVE(solverSolveHighs)
{
   GCG_SOLVERDATA* solverdata;
   int ncols;
   SCIP_RETCODE retval;

   assert(solver != NULL);

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   SCIPdebugMessage("calling exact pricing with HIGHS for pricing problem %d\n", probnr);

   retval = SCIP_OKAY;

   /* set limits to (infinite/zero) default values */
   CHECK_ZERO( Highs_setIntOptionValue(solverdata->highsptr[probnr], "mip_max_nodes", INT_MAX) );
   CHECK_ZERO( Highs_setDoubleOptionValue(solverdata->highsptr[probnr], "mip_rel_gap", 0.0) );
   CHECK_ZERO( Highs_setIntOptionValue(solverdata->highsptr[probnr], "mip_max_improving_sols", INT_MAX) );

   /* solve the pricing problem and evaluate solution */
   SCIP_CALL( solveHighs(gcg, solverdata, pricingprob, probnr, dualsolconv, lowerbound, &ncols, status) );
   assert(*status != GCG_PRICINGSTATUS_OPTIMAL || ncols > 0);

 TERMINATE:
   return retval;
}

/** creates the HIGHS pricing solver and includes it in GCG */
SCIP_RETCODE GCGincludeSolverHighs(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* origprob;
   GCG_SOLVERDATA* solverdata;
   char name[SCIP_MAXSTRLEN];

   SCIP_CALL( SCIPallocMemory(GCGgetDwMasterprob(gcg), &solverdata) );
   origprob = GCGgetOrigprob(gcg);

   SCIP_CALL( GCGpricerIncludeSolver(gcg, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY,
         SOLVER_HEURENABLED, SOLVER_EXACTENABLED,
         solverUpdateHighs, solverSolveHighs, solverSolveHeurHighs, solverFreeHighs, solverInitHighs,
         solverExitHighs, solverInitsolHighs, solverExitsolHighs, solverdata));

   SCIP_CALL( SCIPaddBoolParam(origprob, "pricingsolver/highs/checksols",
         "should solutions of the pricing MIPs be checked for duplicity?",
         &solverdata->checksols, TRUE, DEFAULT_CHECKSOLS, NULL, NULL));

   SCIP_CALL( SCIPaddIntParam(origprob, "pricingsolver/highs/threads",
         "number of threads the HiGHS pricing solver is allowed to use (0: automatic)",
         &solverdata->threads, TRUE, DEFAULT_THREADS, 0, INT_MAX, NULL, NULL));

   SCIP_CALL( SCIPaddLongintParam(origprob, "pricingsolver/highs/startnodelimit",
         "start node limit for heuristic pricing",
         &solverdata->startnodelimit, TRUE, DEFAULT_STARTNODELIMIT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/highs/startgaplimit",
         "start gap limit for heuristic pricing",
         &solverdata->startgaplimit, TRUE, DEFAULT_STARTGAPLIMIT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(origprob, "pricingsolver/highs/startsollimit",
         "start solution limit for heuristic pricing",
         &solverdata->startsollimit, TRUE, DEFAULT_STARTSOLLIMIT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/highs/nodelimitfac",
         "factor by which to increase node limit for heuristic pricing (1.0: add start limit)",
         &solverdata->nodelimitfac, TRUE, DEFAULT_NODELIMITFAC, 1.0, SCIPinfinity(origprob), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/highs/gaplimitfac",
         "factor by which to decrease gap limit for heuristic pricing (1.0: subtract start limit)",
         &solverdata->gaplimitfac, TRUE, DEFAULT_GAPLIMITFAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/highs/sollimitfac",
         "factor by which to increase solution limit for heuristic pricing (1.0: add start limit)",
         &solverdata->sollimitfac, TRUE, DEFAULT_SOLLIMITFAC, 1.0, SCIPinfinity(origprob), NULL, NULL) );
   
   SCIPsnprintf(name, SCIP_MAXSTRLEN, "HiGHS %s", Highs_version());
   SCIP_CALL( SCIPincludeExternalCodeInformation(origprob, name, "High performance serial and parallel solver for large-scale sparse LP, MIP, and QP models (https://highs.dev/)") );

   return SCIP_OKAY;
}
