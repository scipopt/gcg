/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//#define SCIP_DEBUG
/**@file   pricer_gcg.cpp
 * @brief  pricer for generic column generation
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @author Alexander Gross
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <cstring>

#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"

#include "pricer_gcg.h"
#include "objpricer_gcg.h"
#include "sepa_master.h"

#include "relax_gcg.h"
#include "struct_solver.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"
#include "cons_masterbranch.h"
#include "objscip/objscip.h"
#include "objpricer_gcg.h"
#include "class_pricingtype.h"
#include "scip/scip.h"


#ifdef _OPENMP
#include <omp.h>
#endif

using namespace scip;

#define PRICER_NAME            "gcg"
#define PRICER_DESC            "pricer for gcg"
#define PRICER_PRIORITY        5000000
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

#define DEFAULT_MAXSOLSPROB              INT_MAX    /**< maximal number of solution per pricing problem*/
#define DEFAULT_USEHEURPRICING           FALSE      /**< should heuristic pricing be used */
#define DEFAULT_ABORTPRICINGINT          TRUE       /**< should the pricing be aborted when integral */
#define DEFAULT_ABORTPRICINGGAP          0.00       /**< gap at which the pricing is aborted */
#define DEFAULT_SUCCESSFULMIPSREL        1.0        /**< factor of successful mips to be solved */
#define DEFAULT_DISPINFOS                FALSE      /**< should additional information be displayed */
#define DEFAULT_SORTING                  2          /**< default sorting method for pricing mips
                                                     *    0 :   order of pricing problems
                                                     *    1 :   according to dual solution of convexity constraint
                                                     *    2 :   according to reliability from previous round)
                                                     */
#define DEFAULT_THREADS                  0          /**< number of threads (0 is OpenMP default) */

#define EVENTHDLR_NAME         "probdatavardeleted"
#define EVENTHDLR_DESC         "event handler for variable deleted event"

/** small macro to simplify printing pricer information */
#define GCGpricerPrintInfo(scip,pricerdata, ...) do { \
   if( pricerdata->dispinfos ) { \
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,__VA_ARGS__);\
   } else {\
      SCIPdebugMessage(__VA_ARGS__); \
   }\
   }while( FALSE )

#define PRICER_STAT_ARRAYLEN_TIME 1024                /**< length of the array for Time histogram representation */
#define PRICER_STAT_BUCKETSIZE_TIME 10                /**< size of the buckets for Time histogram representation */
#define PRICER_STAT_ARRAYLEN_VARS 1024                /**< length of the array for foundVars histogram representation */
#define PRICER_STAT_BUCKETSIZE_VARS 1                 /**< size of the buckets for foundVars histogram representation */

/*
 * Data structures
 */


/** variable pricer data */
struct SCIP_PricerData
{
   int                   npricingprobs;      /**< number of pricing problems */
   SCIP**                pricingprobs;       /**< pointers to the pricing problems */
   SCIP_Real*            dualsolconv;        /**< array of dual solutions for the convexity constraints */
   SCIP_Real*            solvals;            /**< solution values of variables in the pricing problems */
   int*                  npointsprob;        /**< number of variables representing points created by the pricing probs */
   int*                  nraysprob;          /**< number of variables representing rays created by the pricing probs */
   SCIP_Longint          currnodenr;         /**< current node number in the masterproblem*/
   SCIP_HASHMAP*         mapcons2idx;        /**< hashmap mapping constraints to their index in the conss array */
   SCIP_Real*            score;              /**< score of the pricing problem problems */
   int*                  permu;              /**< current permutation of the pricing problems */
   int                   npricingprobsnotnull; /**< number of non-Null pricing problems*/

   SCIP_VAR**            pricedvars;         /**< array of all priced variables */
   int                   npricedvars;        /**< number of priced variables */
   int                   maxpricedvars;      /**< maximal number of priced variables */

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
   SCIP_VARTYPE          vartype;            /**< vartype of created master variables */
   int                   maxsolsprob;        /**< maximal number of solutions per pricing problem */
   int                   nroundsredcost;     /**< number of reduced cost rounds */
   int                   sorting;            /**< how should pricing problems be sorted */
   SCIP_Bool             useheurpricing;     /**< should heuristic pricing be used */
   SCIP_Bool             abortpricingint;    /**< should the pricing be aborted on integral solutions */
   SCIP_Bool             dispinfos;          /**< should pricing information be displayed*/
   SCIP_Real             successfulmipsrel;  /**< Factor of successful MIPs solved until pricing be aborted */
   SCIP_Real             abortpricinggap;    /**< Gap at which pricing should be aborted */


   /** statistics */
   int                   oldvars;            /**< Vars of last pricing iteration */
   int*                  farkascallsdist;    /**< Calls of each farkas pricing problem */
   int*                  farkasfoundvars;    /**< Found vars of each farkas pricing problem */
   double*               farkasnodetimedist; /**< Time spend in each farkas pricing problem */

   int*                  redcostcallsdist;   /**< Calls of each redcost pricing problem */
   int*                  redcostfoundvars;   /**< Found vars of each redcost pricing problem */
   double*               redcostnodetimedist; /**< Time spend in each redcost pricing problem */

   int*                  nodetimehist;       /**< Histogram of nodetime distribution */
   int*                  foundvarshist;      /**< Histogram of foundvars distribution */

   double                rootnodedegeneracy; /**< degeneracy of the root node */
   double                avgrootnodedegeneracy; /**< average degeneray of all nodes */
   int                   ndegeneracycalcs;   /**< number of observations */
};


/*
 * Callback methods of event handler
 */

/** destructor of event handler to free user data (called when SCIP is exiting) */
#define eventFreeVardeleted NULL

/** initialization method of event handler (called after problem was transformed) */
#define eventInitVardeleted NULL

/** deinitialization method of event handler (called before transformed problem is freed) */
#define eventExitVardeleted NULL

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
#define eventInitsolVardeleted NULL

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#define eventExitsolVardeleted NULL

/** frees specific event data */
#define eventDeleteVardeleted NULL

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecVardeleted)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_VAR** origvars;
   int i;

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
   assert(pricerdata != NULL);

   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARDELETED);
   var = SCIPeventGetVar(event);
   assert(var != NULL);

   SCIPdebugMessage("remove master variable %s from pricerdata and corresponding original variables\n", SCIPvarGetName(var));

   assert(GCGvarIsMaster(var));
   origvars = GCGmasterVarGetOrigvars(var);
   assert(origvars != NULL);

   /* remove master variable from corresponding pricing original variables */
   for( i = 0; i < GCGmasterVarGetNOrigvars(var); ++i )
   {
      SCIP_CALL( GCGoriginalVarRemoveMasterVar(scip, origvars[i], var) );
   }

   /* remove variable from array of stored priced variables */
   for( i = 0; i < pricerdata->npricedvars; ++i )
   {
      if( pricerdata->pricedvars[i] == var )
      {
         /* drop vardeleted event on variable */
         SCIP_CALL( SCIPdropVarEvent(scip, pricerdata->pricedvars[i], SCIP_EVENTTYPE_VARDELETED,
               pricerdata->eventhdlr, NULL, -1) );

         SCIP_CALL( SCIPreleaseVar(scip, &(pricerdata->pricedvars[i])) );
         (pricerdata->npricedvars)--;
         pricerdata->pricedvars[i] = pricerdata->pricedvars[pricerdata->npricedvars];

         break;
      }
   }
   assert(i <= pricerdata->npricedvars);
#ifndef NDEBUG
   for( ; i < pricerdata->npricedvars; ++i )
   {
      assert(pricerdata->pricedvars[i] != var);
   }
#endif

   return SCIP_OKAY;
}


/*
 * Local methods
 */


/** return TRUE or FALSE whether the master LP is solved to optimality */
SCIP_Bool ObjPricerGcg::isMasterLPOptimal()
{
   assert(GCGisMaster(scip_));

   return SCIPgetLPSolstat(scip_) == SCIP_LPSOLSTAT_OPTIMAL;
}


/** return TRUE or FALSE whether pricing problem has been solved to optimality */
SCIP_Bool  ObjPricerGcg::isPricingOptimal(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(!GCGisMaster(scip) && !GCGisOriginal(scip));

   return SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL;
}

/** ensures size of pricedvars array */
SCIP_RETCODE ObjPricerGcg::ensureSizePricedvars(
   int                   size                /**< needed size */
   )
{
   assert(pricerdata->pricedvars != NULL);

   if( pricerdata->maxpricedvars < size )
   {
      int oldsize;

      oldsize = pricerdata->maxpricedvars;
      pricerdata->maxpricedvars = SCIPcalcMemGrowSize(scip_, size);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip_, &(pricerdata->pricedvars), oldsize, pricerdata->maxpricedvars) );
   }
   assert(pricerdata->maxpricedvars >= size);

   return SCIP_OKAY;
}


/** ensures size of solvers array */
SCIP_RETCODE ObjPricerGcg::ensureSizeSolvers()
{
   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));

   if( pricerdata->nsolvers == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip_, &(pricerdata->solvers), 1) ); /*lint !e506*/
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip_, &(pricerdata->solvers), pricerdata->nsolvers+1) );
   }

   return SCIP_OKAY;
}

#ifdef ENABLESTATISTICS
/** gets the NodeTimeDistribution in the form of a histogram */
static
void GCGpricerGetNodeTimeHistogram(
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   SCIP_Real             time                /**< time the pricingproblem needed */
   )
{
   int i;

   /* 1000* because mapping milliseconds on the index i */
   i = 1000*time/PRICER_STAT_BUCKETSIZE_TIME; /*lint !e524 */

   if( i >= PRICER_STAT_ARRAYLEN_TIME )
   {
      i = PRICER_STAT_ARRAYLEN_TIME-1;
   }
   pricerdata->nodetimehist[i]++;

}


/** gets the FoundVarsDistribution in form of a histogram */
static
void GCGpricerGetFoundVarsHistogram(
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   int                   foundvars           /**< foundVars in pricingproblem */
   )
{
   int i;
   i = foundvars/PRICER_STAT_BUCKETSIZE_VARS;

   if( i >= PRICER_STAT_ARRAYLEN_VARS )
   {
      i = PRICER_STAT_ARRAYLEN_VARS-1;
   }
   pricerdata->foundvarshist[i]++;

}


/** gets the statistics of the pricingprobs like calls, foundvars and time */
static
void GCGpricerCollectStatistic(
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   GCG_PRICETYPE         type,               /**< type of pricing: optimal or heuristic */
   int                   probindex,          /**< index of the pricingproblem */
   SCIP_Real             time                /**< time the pricingproblem needed */
   )
{
   int foundvars;

   foundvars = pricerdata->npricedvars - pricerdata->oldvars;

   if( type == GCG_PRICETYPE_FARKAS )
   {

      pricerdata->farkascallsdist[probindex]++; /*Calls*/
      pricerdata->farkasfoundvars[probindex] += foundvars;
      pricerdata->farkasnodetimedist[probindex] += time;   /*Time*/

   }
   else if( type == GCG_PRICETYPE_REDCOST )
   {

      pricerdata->redcostcallsdist[probindex]++;
      pricerdata->redcostfoundvars[probindex] += foundvars;
      pricerdata->redcostnodetimedist[probindex] += time;

   }
   pricerdata->oldvars = pricerdata->npricedvars;

   GCGpricerGetNodeTimeHistogram(pricerdata, time);
   GCGpricerGetFoundVarsHistogram(pricerdata, foundvars);

}
#endif

/** frees all solvers */
SCIP_RETCODE ObjPricerGcg::solversFree()
{
   int i;

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverfree != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverfree(scip_, pricerdata->solvers[i]) );

         BMSfreeMemoryArray(&pricerdata->solvers[i]->name);
         BMSfreeMemoryArray(&pricerdata->solvers[i]->description);

         SCIP_CALL( SCIPfreeClock(scip_, &(pricerdata->solvers[i]->optfarkasclock)) );
         SCIP_CALL( SCIPfreeClock(scip_, &(pricerdata->solvers[i]->optredcostclock)) );
         SCIP_CALL( SCIPfreeClock(scip_, &(pricerdata->solvers[i]->heurfarkasclock)) );
         SCIP_CALL( SCIPfreeClock(scip_, &(pricerdata->solvers[i]->heurredcostclock)) );

         SCIPfreeMemory(scip, &(pricerdata->solvers[i]));
      }
   }

   return SCIP_OKAY;
}

/** calls the init method on all solvers */
SCIP_RETCODE ObjPricerGcg::solversInit()
{
   int i;

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverinit != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverinit(scip_, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}

/** calls the exit method on all solvers */
SCIP_RETCODE ObjPricerGcg::solversExit()
{
   int i;

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverexit != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverexit(scip_, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}

/** calls the initsol method on all solvers */
SCIP_RETCODE ObjPricerGcg::solversInitsol()
{
   int i;

   if( pricerdata->npricingprobs == 0 )
      return SCIP_OKAY;

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverinitsol != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverinitsol(scip_, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}

/** calls the exitsol method of all solvers */
SCIP_RETCODE ObjPricerGcg::solversExitsol()
{
   int i;

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverexitsol != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverexitsol(scip_, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}

/** returns the gegeneracy of the masterproblem */
SCIP_RETCODE ObjPricerGcg::computeCurrentDegeneracy(
   double*               degeneracy          /**< pointer to store degeneracy */
   )
{
   int ncols;
   int nrows;
   int i;
   int count;
   int countz;
   int colindex;
   double currentVal;
   int* indizes;
   SCIP_COL** cols;
   SCIP_VAR* var;

   assert(degeneracy != NULL);

   *degeneracy = 0.0;
   ncols = SCIPgetNLPCols(scip_);
   nrows = SCIPgetNLPRows(scip_);
   cols = SCIPgetLPCols(scip_);

   SCIP_CALL( SCIPallocMemoryArray(scip_, &indizes, ncols+nrows) );

   for( i = 0; i < ncols+nrows; i++ )
   {
      indizes[i] = 0;
   }

   /* gives indices of Columns in Basis and indices of vars in Basis */
   SCIP_CALL( SCIPgetLPBasisInd(scip_, indizes) );

   countz = 0;
   count = 0;

   for( i = 0; i < nrows; i++ )
   {
      colindex = indizes[i];
      /* is column if >0 it is column in basis, <0 is for row */
      if( colindex > 0 )
      {
         var = SCIPcolGetVar(cols[colindex]);

         currentVal = SCIPgetSolVal(scip_, NULL, var);

         if( SCIPisEQ(scip_, currentVal, 0.0) )
            countz++;

         count++;
      }
   }

   /* Degeneracy in % */
   if( count > 0 )
      *degeneracy = ((double)countz / count);

   assert(*degeneracy <= 1.0 && *degeneracy >= 0);

   SCIPfreeMemoryArray(scip_, &indizes);

   return SCIP_OKAY;
}

/** initializes the pointers to the appropriate structures */
SCIP_RETCODE ObjPricerGcg::getSolverPointers(
   GCG_SOLVER*           solver,             /**< pricing solver */
   PricingType*          pricetype,          /**< type of pricing: optimal or heuristic */
   SCIP_Bool             optimal,            /**< should the pricing problem be solved optimal or heuristically */
   SCIP_CLOCK**          clock,              /**< clock belonging to this setting */
   int**                 calls,              /**< calls belonging to this setting */
   GCG_DECL_SOLVERSOLVE((**solversolve))     /**< solving function belonging to this setting */
   )
{
   assert(solver != NULL);
   assert(clock != NULL);
   assert(calls != NULL);
   switch( optimal )
   {
   case TRUE:
      if( pricetype->getType() == GCG_PRICETYPE_FARKAS )
      {
         *clock = solver->optfarkasclock;
         *calls = &(solver->optfarkascalls);
      }
      else
      {
         *clock = solver->optredcostclock;
         *calls = &(solver->optredcostcalls);
      }
      *solversolve = solver->solversolve;
      break;
   case FALSE:
      if( pricetype->getType() == GCG_PRICETYPE_FARKAS )
      {
         *clock = solver->heurfarkasclock;
         *calls = &(solver->heurfarkascalls);
      }
      else
      {
         *clock = solver->heurredcostclock;
         *calls = &(solver->heurredcostcalls);
      }
      *solversolve = solver->solversolveheur;
      break;
   default:
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/** set subproblem timelimit */
SCIP_RETCODE ObjPricerGcg::setPricingProblemTimelimit(
   SCIP*                 pricingscip         /**< SCIP of the pricingproblem */
   )
{
   SCIP_Real timelimit;

   assert(pricingscip != NULL);

   /* set time limit */
   SCIP_CALL( SCIPgetRealParam(scip_, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip_, timelimit) )
   {
      if( timelimit - SCIPgetSolvingTime(scip_) > 0 )
      {
         SCIP_CALL( SCIPsetRealParam(pricingscip, "limits/time", timelimit - SCIPgetSolvingTime(scip_)) );
         SCIPdebugMessage("Time limit for prob <%s> is %f\n", SCIPgetProbName(pricingscip), timelimit- SCIPgetSolvingTime(scip_));
      }
      else
      {
         SCIPdebugMessage("Time limit for prob <%s> is < 0\n", SCIPgetProbName(pricingscip));
      }
   }
   return SCIP_OKAY;
}

/** set subproblem memory limit */
SCIP_RETCODE ObjPricerGcg::setPricingProblemMemorylimit(
   SCIP*                 pricingscip         /**< SCIP of the pricingproblem */
   )
{
   SCIP_Real memlimit;

   assert(pricingscip != NULL);

   assert(GCGisOriginal(origprob));

   SCIP_CALL( SCIPgetRealParam(origprob, "limits/memory", &memlimit) );

   if( !SCIPisInfinity(origprob, memlimit) )
   {
      memlimit -= SCIPgetMemUsed(origprob)/1048576.0 + GCGgetPricingprobsMemUsed(origprob) - SCIPgetMemUsed(pricingscip)/1048576.0;
      if( memlimit < 0 )
         memlimit = 0.0;
      SCIP_CALL( SCIPsetRealParam(pricingscip, "limits/memory", memlimit) );
   }

   return SCIP_OKAY;
}


/** set all pricing problem limits */
SCIP_RETCODE ObjPricerGcg::setPricingProblemLimits(
   int                   prob,               /**< index of the pricing problem */
   SCIP_Bool             optimal            /**< heuristic or optimal pricing */
   )
{
   assert(pricerdata != NULL);
   assert(prob >= 0 && prob < pricerdata->npricingprobs);

   /** @todo set objective limit, such that only solutions with negative reduced costs are accepted? */
   if( !optimal )
   {
      SCIP_CALL( SCIPsetObjlimit(pricerdata->pricingprobs[prob], pricerdata->dualsolconv[prob]) );
   }

   SCIP_CALL( setPricingProblemTimelimit(pricerdata->pricingprobs[prob]) );
   SCIP_CALL( setPricingProblemMemorylimit(pricerdata->pricingprobs[prob]) );

   return SCIP_OKAY;
}

/** solves a specific pricing problem
 * @todo simplify
 */
SCIP_RETCODE ObjPricerGcg::solvePricingProblem(
   int                   prob,               /**< index of pricing problem */
   PricingType*          pricetype,          /**< type of pricing: optimal or heuristic */
   SCIP_Bool             optimal,            /**< should the pricing problem be solved optimal or heuristically */
   SCIP_Real*            lowerbound,         /**< dual bound returned by pricing problem */
   SCIP_SOL**            sols,               /**< pointer to store solutions */
   SCIP_Bool*            solisray,           /**< array to indicate whether solution is a ray */
   int                   maxsols,            /**< size of the sols array to indicate maximum solutions */
   int*                  nsols,              /**< number of solutions */
   SCIP_STATUS*          status              /**< solution status of the pricing problem */
   )
{
   SCIP_RETCODE retcode;
   int i;

   assert(pricerdata->pricingprobs[prob] != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   *status = SCIP_STATUS_UNKNOWN;

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      SCIP_CLOCK* clock;
      int* calls;
      GCG_SOLVER* solver;
      GCG_DECL_SOLVERSOLVE((*solversolve));

      #pragma omp critical (limits)
      {
         retcode = setPricingProblemLimits(prob, optimal);
      }
      SCIP_CALL( retcode );

      solver = pricerdata->solvers[i];
      assert(solver != NULL);

      SCIP_CALL( getSolverPointers(solver, pricetype, optimal, &clock, &calls, &solversolve) );
      assert(solversolve == solver->solversolve || solversolve == solver->solversolveheur);

      /* continue if the appropriate solver is not available */
      if( solversolve == NULL )
      {
         continue;
      }
      #pragma omp critical (clock)
      {
         SCIP_CALL_ABORT( SCIPstartClock(scip_, clock) );
      }

      SCIP_CALL( solversolve(pricerdata->pricingprobs[prob], solver, prob, lowerbound, sols, solisray, maxsols, nsols, status) );

      if(optimal)
      {
         #pragma omp atomic
         pricerdata->solvedsubmipsoptimal++;
      }
      else
      {
         #pragma omp atomic
         pricerdata->solvedsubmipsheur++;
      }

      #pragma omp critical (clock)
      {
         SCIP_CALL_ABORT( SCIPstopClock(scip_, clock) );
      }

      if( *status != SCIP_STATUS_UNKNOWN )
      {
         #pragma omp atomic
         (*calls)++;
      }

      if( *status == SCIP_STATUS_OPTIMAL || *status == SCIP_STATUS_UNBOUNDED )
      {
         if( optimal )
         {

#ifdef ENABLESTATISTICS
            #pragma omp critical (collectstats)
            GCGpricerCollectStatistic(pricerdata, pricetype->getType(), prob,
               SCIPgetSolvingTime(pricerdata->pricingprobs[prob]));
#endif
            if( SCIPgetStage(pricerdata->pricingprobs[prob]) > SCIP_STAGE_SOLVING )
            {
               #pragma omp atomic
               pricerdata->pricingiters += SCIPgetNLPIterations(pricerdata->pricingprobs[prob]);
            }
         }
         break;
      }
   }

   return SCIP_OKAY;
}

/** computes the pricing problem objectives
 *  @todo this method could use more parameters as it is private
 */
SCIP_RETCODE ObjPricerGcg::setPricingObjs(
   PricingType*          pricetype           /**< Farkas or Reduced cost pricing */
   )
{
   SCIP_CONS** origconss;
   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_VAR** probvars;
   int nprobvars;

   SCIP_ROW** mastercuts;
   int nmastercuts;
   SCIP_ROW** origcuts;
   SCIP_COL** cols;
   SCIP_Real* consvals;
   SCIP_Real dualsol;

   SCIP_VAR** consvars;
   int nconsvars;

   int i;
   int j;

   /* get the constraints of the master problem and the corresponding constraints in the original problem */
   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);
   origconss = GCGrelaxGetLinearOrigMasterConss(origprob);

   /* set objective value of all variables in the pricing problems to 0 (for farkas pricing) /
    * to the original objective of the variable (for redcost pricing)
    */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if( pricerdata->pricingprobs[i] == NULL )
         continue;
      probvars = SCIPgetVars(pricerdata->pricingprobs[i]);
      nprobvars = SCIPgetNVars(pricerdata->pricingprobs[i]);

      for( j = 0; j < nprobvars; j++ )
      {
         assert(GCGvarGetBlock(probvars[j]) == i);
         assert( GCGvarIsLinking(GCGpricingVarGetOrigvars(probvars[j])[0]) || (GCGvarGetBlock(GCGpricingVarGetOrigvars(probvars[j])[0]) == i));
         SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], pricetype->varGetObj(probvars[j])));
      }
   }

   /* compute reduced cost for linking variable constraints and update objectives in the pricing problems */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if( pricerdata->pricingprobs[i] == NULL )
         continue;
      probvars = SCIPgetVars(pricerdata->pricingprobs[i]);
      nprobvars = SCIPgetNVars(pricerdata->pricingprobs[i]);

      for( j = 0; j < nprobvars; j++ )
      {
         SCIP_VAR* origvar;
         SCIP_CONS** linkconss;
#ifndef NDEBUG
         SCIP_VAR** pricingvars;
#endif
         assert(GCGvarIsPricing(probvars[j]));
         assert(GCGvarGetBlock(probvars[j]) == i);

         origvar = GCGpricingVarGetOrigvars(probvars[j])[0];
         if( !GCGvarIsLinking(origvar) )
            continue;

#ifndef NDEBUG
         pricingvars = GCGlinkingVarGetPricingVars(origvar);
#endif
         linkconss = GCGlinkingVarGetLinkingConss(origvar);
         assert(pricingvars[i] == probvars[j]);
         assert(linkconss[i] != NULL);

         dualsol = pricetype->consGetDual(scip_, linkconss[i]);

         /* add dual solution value to the pricing variable:
          * lambda variables get coef -1 in linking constraints --> add dualsol
          */
         SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[i], probvars[j], dualsol) );
      }
   }

   /* compute reduced cost and update objectives in the pricing problems */
   for( i = 0; i < nmasterconss; i++ )
   {

      dualsol = pricetype->consGetDual(scip_, masterconss[i]);

      if( !SCIPisZero(scip_, dualsol) )
      {
#ifdef PRINTDUALSOLS
         SCIPdebugMessage("mastercons <%s> dualsol: %g\n", SCIPconsGetName(masterconss[i]), dualsol);
#endif

         /* for all variables in the constraint, modify the objective of the corresponding variable in a pricing problem */
         consvars = SCIPgetVarsLinear(origprob, origconss[i]);
         consvals = SCIPgetValsLinear(origprob, origconss[i]);
         nconsvars = SCIPgetNVarsLinear(origprob, origconss[i]);
         for( j = 0; j < nconsvars; j++ )
         {
            int blocknr;
            blocknr = GCGvarGetBlock(consvars[j]);
            assert(GCGvarIsOriginal(consvars[j]));
            /* nothing to be done if variable belongs to redundant block or variable was directly transferred to the master
             * or variable is linking variable (which means, the directly transferred copy is part of the master cons)
             */
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
   mastercuts = GCGsepaGetMastercuts(scip_);
   nmastercuts = GCGsepaGetNMastercuts(scip_);
   origcuts = GCGsepaGetOrigcuts(scip_);

   assert(mastercuts != NULL);
   assert(origcuts != NULL);
   assert(GCGsepaGetNOrigcuts(scip_) == nmastercuts);

   /* compute reduced cost and update objectives in the pricing problems */
   for( i = 0; i < nmastercuts; i++ )
   {
      dualsol = pricetype->rowGetDual(mastercuts[i]);

      if( !SCIPisZero(scip_, dualsol) )
      {
         /* get columns and vals of the cut */
         nconsvars = SCIProwGetNNonz(origcuts[i]);
         cols = SCIProwGetCols(origcuts[i]);
         consvals = SCIProwGetVals(origcuts[i]);

         /* get the variables corresponding to the columns in the cut */
         SCIP_CALL( SCIPallocMemoryArray(scip, &consvars, nconsvars) );
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
         SCIPfreeMemoryArray(scip, &consvars);
      }
   }

   /* get dual solutions / farkas values of the convexity constraints */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {

      assert( GCGrelaxIsPricingprobRelevant(origprob, i) == (GCGrelaxGetConvCons(origprob, i) != NULL) );
      if( !GCGrelaxIsPricingprobRelevant(origprob, i) )
      {
         pricerdata->dualsolconv[i] = -1.0 * SCIPinfinity(scip_);
         continue;
      }

      pricerdata->dualsolconv[i] = pricetype->consGetDual(scip_, GCGrelaxGetConvCons(origprob, i));
#ifdef PRINTDUALSOLS
      if( GCGrelaxIsPricingprobRelevant(origprob, i) )
      {
         SCIPdebugMessage("convcons <%s> dualsol: %g\n", SCIPconsGetName(GCGrelaxGetConvCons(origprob, i)), pricerdata->dualsolconv[i]);
      }
#endif
   }

   return SCIP_OKAY;
}

/** add master variable to all constraints */
SCIP_RETCODE ObjPricerGcg::addVariableToMasterconstraints(
   SCIP_VAR*             newvar,             /**< The new variable to add */
   int                   prob,               /**< number of the pricing problem the solution belongs to */
   SCIP_VAR**            solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
   SCIP_Real*            solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
   int                   nsolvars            /**< number of variables in array solvars */
   )
{
   int i;
   int c;
   int idx;

   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_Real* mastercoefs;
   SCIP_CONS* linkcons;

   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);

   SCIP_CALL( SCIPallocBufferArray(scip_, &mastercoefs, nmasterconss) );
   BMSclearMemoryArray(mastercoefs, nmasterconss);

   /* compute coef of the variable in the master constraints */
   for( i = 0; i < nsolvars; i++ )
   {
      if( !SCIPisZero(scip_, solvals[i]) )
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
         assert(!SCIPisInfinity(scip_, solvals[i]));

         /* original variable is a linking variable, just add it to the linkcons */
         if( GCGvarIsLinking(origvars[0]) )
         {
#ifndef NDEBUG
            SCIP_VAR** pricingvars;
            pricingvars = GCGlinkingVarGetPricingVars(origvars[0]);
#endif
            linkconss = GCGlinkingVarGetLinkingConss(origvars[0]);

            assert(pricingvars[prob] == solvars[i]);
            assert(linkconss[prob] != NULL);
            SCIP_CALL( SCIPaddCoefLinear(scip_, linkconss[prob], newvar, -solvals[i]) );
            continue;
         }

         /* for each coef, add coef * solval to the coef of the new variable for the corresponding constraint */
         for( c = 0; c < ncoefs; c++ )
         {
            linkconss = GCGoriginalVarGetMasterconss(origvars[0]);
            assert(!SCIPisZero(scip_, coefs[c]));
            SCIP_CALL( SCIPgetTransformedCons(scip_, linkconss[c], &linkcons) );

            idx = (int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, linkcons); /*lint !e507*/
            assert(0 <= idx && idx < nmasterconss);
            assert(masterconss[idx] == linkcons);
            mastercoefs[idx] += coefs[c] * solvals[i];
         }

      }
   }

   /* add the variable to the master constraints */
   for( i = 0; i < nmasterconss; i++ )
   {
      if( !SCIPisZero(scip_, mastercoefs[i]) )
      {
         assert(!SCIPisInfinity(scip_, mastercoefs[i]) && !SCIPisInfinity(scip_, -mastercoefs[i]));
         SCIP_CALL( SCIPaddCoefLinear(scip_, masterconss[i], newvar, mastercoefs[i]) );
      }
   }

   SCIPfreeBufferArray(scip_, &mastercoefs);
   return SCIP_OKAY;
}



/** add variable with computed coefficients to the master cuts */
static
SCIP_RETCODE addVariableToMastercuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             newvar,             /**< The new variable to add */
   int                   prob,               /**< number of the pricing problem the solution belongs to */
   SCIP_VAR**            solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
   SCIP_Real*            solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
   int                   nsolvars            /**< number of variables in array solvars */
   )
{
   SCIP_ROW** mastercuts;
   int nmastercuts;
   SCIP_ROW** origcuts;

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

   assert(mastercuts != NULL);
   assert(origcuts != NULL);
   assert(GCGsepaGetNOrigcuts(scip) == nmastercuts);

   /* compute coef of the variable in the cuts and add it to the cuts */
   for( i = 0; i < nmastercuts; i++ )
   {
      if( !SCIProwIsInLP(mastercuts[i]) )
         continue;

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
            for( k = 0; k < nsolvars; k++ )
               if( solvars[k] == GCGoriginalVarGetPricingVar(var) )
               {
                  conscoef += ( consvals[j] * solvals[k] );
                  break;
               }
      }

      if( !SCIPisZero(scip, conscoef) )
         SCIP_CALL( SCIPaddVarToRow(scip , mastercuts[i], newvar, conscoef) );
   }

   return SCIP_OKAY;
}

/** adds new variable to the end of the priced variables array */
SCIP_RETCODE ObjPricerGcg::addVariableToPricedvars(
   SCIP_VAR*             newvar              /**< variable to add */
   )
{
   SCIP_CALL( ensureSizePricedvars(pricerdata->npricedvars + 1) );
   pricerdata->pricedvars[pricerdata->npricedvars] = newvar;
   pricerdata->npricedvars++;

   return SCIP_OKAY;
}

SCIP_Real ObjPricerGcg::computeRedCost(
   SCIP_VAR**            solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
   SCIP_Real*            solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
   int                   nsolvars,           /**< number of variables in array solvars */
   SCIP_Bool             solisray,           /**< is the solution a ray? */
   int                   prob               /**< number of the pricing problem the solution belongs to */
   )
{
   SCIP_Real objvalue;
   objvalue = 0.0;

   /* compute the objective function value of the solution */
   for( int i = 0; i < nsolvars; i++ )
      objvalue += solvals[i] * SCIPvarGetObj(solvars[i]);

   /* compute reduced cost of variable (i.e. subtract dual solution of convexity constraint, if solution corresponds to a point) */
   return (solisray ? objvalue : objvalue - pricerdata->dualsolconv[prob]);
}

/** creates a new master variable corresponding to the given solution and problem */
SCIP_RETCODE ObjPricerGcg::createNewMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
   SCIP_Real*            solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
   int                   nsolvars,           /**< number of variables in array solvars */
   SCIP_Bool             solisray,           /**< is the solution a ray? */
   int                   prob,               /**< number of the pricing problem the solution belongs to */
   SCIP_Bool             force,              /**< should the given variable be added also if it has non-negative reduced cost? */
   SCIP_Bool*            added,              /**< pointer to store whether the variable was successfully added */
   SCIP_VAR**            addedvar            /**< pointer to store the created variable */
   )
{
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
   assert(pricerdata != NULL);

   if( addedvar != NULL )
      *addedvar = NULL;

   objvalue = 0.0;
   redcost = 0.0;

   if( !force )
   {
      /* compute the objective function value of the solution */
      redcost = computeRedCost(solvars, solvals, nsolvars, solisray == TRUE, prob);

      if( !SCIPisSumNegative(scip, redcost) )
      {
         SCIPdebugMessage("var with redcost %g (objvalue = %g, dualsol =%g) was not added\n", redcost, objvalue, pricerdata->dualsolconv[prob]);
         *added = FALSE;

         return SCIP_OKAY;
      }
      SCIPdebugMessage("found var with redcost %g (objvalue = %g, dualsol =%g)\n", redcost, objvalue, pricerdata->dualsolconv[prob]);
   }
   else
   {
      SCIPdebugMessage("force var (objvalue = %g, dualsol =%g)\n",  objvalue, pricerdata->dualsolconv[prob]);
   }

   *added = TRUE;

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
      SCIPwarningMessage(scip, "variable with infinite objective value found in pricing, change objective to SCIPinfinity()/2\n");
      objcoeff = SCIPinfinity(scip) / 2;
   }

   if( solisray == 1 )
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
         pricerdata->vartype, solisray == 1, prob, nsolvars, solvals, solvars));

   SCIPvarMarkDeletable(newvar);

   SCIP_CALL( SCIPcatchVarEvent(scip, newvar, SCIP_EVENTTYPE_VARDELETED,
         pricerdata->eventhdlr, NULL, NULL) );


   /* add variable */
   if( !force )
   {
      SCIP_CALL( SCIPaddPricedVar(scip, newvar, pricerdata->dualsolconv[prob] - objvalue) );
   }
   else
   {
      SCIP_CALL( SCIPaddVar(scip, newvar) );
   }

   SCIP_CALL( addVariableToPricedvars(newvar) );
   SCIP_CALL( addVariableToMasterconstraints(newvar, prob, solvars, solvals, nsolvars) );
   SCIP_CALL( addVariableToMastercuts(scip, newvar, prob, solvars, solvals, nsolvars) );

   /* add variable to convexity constraint */
   if( !(solisray == 1))
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, GCGrelaxGetConvCons(origprob, prob), newvar, 1.0) );
   }

   if( addedvar != NULL )
   {

      *addedvar = newvar;
   }

   GCGupdateVarStatistics(scip, origprob, newvar, redcost);

   return SCIP_OKAY;
}

/**
 * check whether pricing can be aborted:
 * if objective value is always integral and the current node's current
 * lowerbound rounded up equals the current lp objective value rounded
 * up we don't need to continue pricing since the best possible feasible
 * solution must have at least this value
 */
SCIP_Bool  ObjPricerGcg::canPricingBeAborted()
{
   SCIP_Bool canabort;

   canabort = FALSE;

   if( pricerdata->abortpricingint && SCIPisObjIntegral(scip_)
      && SCIPisEQ(scip_, SCIPceil(scip_, SCIPgetNodeLowerbound(scip_, SCIPgetCurrentNode(scip_))), SCIPceil(scip_, SCIPgetLPObjval(scip_))) /* && SCIPgetNNodes(scip) > 1 ??????*/)
   {
      GCGpricerPrintInfo(scip_, pricerdata, "pricing aborted due to integral objective: node LB = %g, LP obj = %g\n",
            SCIPgetNodeLowerbound(scip_, SCIPgetCurrentNode(scip_)), SCIPgetLPObjval(scip_));

      canabort = TRUE;
   }

   if( !canabort && pricerdata->abortpricinggap > 0 )
   {
      SCIP_Real gap;
      gap = (SCIPgetLPObjval(scip_) - SCIPgetNodeLowerbound(scip_, SCIPgetCurrentNode(scip_)))/SCIPgetNodeLowerbound(scip_, SCIPgetCurrentNode(scip_));
      gap = ABS(gap);

      if( gap < pricerdata->abortpricinggap )
      {
         GCGpricerPrintInfo(scip_, pricerdata, "pricing aborted due to small gap: node LB = %g, LP obj = %g, gap = %g\n",
               SCIPgetNodeLowerbound(scip_, SCIPgetCurrentNode(scip_)), SCIPgetLPObjval(scip_), gap);

         canabort = TRUE;
      }
   }

   return canabort;
}

/** sorts pricing problems according to their score */
void ObjPricerGcg::sortPricingProblemsByScore()
{
   int i;
   assert(pricerdata != NULL);
    /** @todo sort w.r.t. other measures? Don't sort in Farkas pricing? Randomized? */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      pricerdata->permu[i] = i;
      switch( pricerdata->sorting )
      {
      case 1:
         pricerdata->score[i] = pricerdata->dualsolconv[i];
         break;
      case 2:
         pricerdata->score[i] = -(0.2 * pricerdata->npointsprob[i] + pricerdata->nraysprob[i]);
         break;
      default:
         pricerdata->score[i] = 0.0;
         break;
      }
   }

   if( pricerdata->sorting > 0 )
      SCIPsortDownRealInt(pricerdata->score, pricerdata->permu, pricerdata->npricingprobs);
}


/** returns whether pricing can be aborted */
SCIP_Bool ObjPricerGcg::abortPricing(
   PricingType*          pricetype,          /**< type of pricing*/
   int                   nfoundvars,         /**< number of variables found so far */
   int                   solvedmips,         /**< number of MIPS solved so far */
   int                   successfulmips,     /**< number of sucessful mips solved so far */
   SCIP_Bool             optimal             /**< optimal or heuristic pricing */
)
{
   if( optimal )
      return pricetype->canOptimalPricingBeAborted(nfoundvars, solvedmips, successfulmips, pricerdata->successfulmipsrel, pricerdata->npricingprobsnotnull);
   else
      return pricetype->canHeuristicPricingBeAborted(nfoundvars, solvedmips, successfulmips, pricerdata->successfulmipsrel, pricerdata->npricingprobsnotnull);
}


/** free pricing problems */
SCIP_RETCODE ObjPricerGcg::freePricingProblems()
{
   int j;
   assert(pricerdata != NULL);
   assert(pricerdata->pricingprobs != NULL);

   for( j = 0; j < pricerdata->npricingprobs; j++ )
      if( pricerdata->pricingprobs[j] != NULL
         && SCIPgetStage(pricerdata->pricingprobs[j]) > SCIP_STAGE_PROBLEM)
         {
            SCIP_CALL( SCIPstartClock(scip_, pricerdata->freeclock) );
            SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[j]) );
            SCIP_CALL( SCIPstopClock(scip_, pricerdata->freeclock) );
         }

   return SCIP_OKAY;
}

int ObjPricerGcg::countPricedVariables(
   int& prob,
   SCIP_SOL** sols,
   int nsols,
   SCIP_Bool* solisray
   )
{
   SCIP_VAR** solvars;
   SCIP_Real* solvals;
   int nsolvars;
   int nfoundvars;

   nfoundvars = 0;
   solvars = SCIPgetOrigVars(pricerdata->pricingprobs[prob]);
   nsolvars = SCIPgetNOrigVars(pricerdata->pricingprobs[prob]);
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &solvals, nsolvars));
   for( int j = 0; j < nsols; ++j )
   {
      SCIPdebugMessage("Solution %d of prob %d (%p)\n", j, prob, sols[j]);
      SCIP_CALL_ABORT( SCIPgetSolVals(pricerdata->pricingprobs[prob], sols[j], nsolvars, solvars, solvals) );
      if ( SCIPisNegative(scip_, computeRedCost(solvars, solvals, nsolvars, solisray[j], prob)) )
      {
         nfoundvars += 1;
      }
   }
   SCIPfreeMemoryArray(scip_, &solvals);
   return nfoundvars;
}

/** performs optimal or farkas pricing */
SCIP_RETCODE ObjPricerGcg::performPricing(
   PricingType*          pricetype,          /**< type of pricing */
   SCIP_Bool             optimal,            /**< heuristic or optimal pricing */
   SCIP_RESULT*          result,             /**< result pointer */
   int*                  nfoundvars,         /**< pointer to store number of found variables */
   SCIP_Real*            bestredcost,        /**< pointer to store reduced cost */
   SCIP_Bool*            bestredcostvalid    /**< pointer to store whether the reduced cost returned is valid */
   )
{
   int i;
   int j;

   int solvedmips;
   int successfulmips;
   int nfoundvarsprob;
   SCIP_Bool added;

   SCIP_Real pricinglowerbound;

   SCIP_SOL*** sols;
   int* nsols;
   SCIP_Bool** solisray;
   int maxsols;
   SCIP_RETCODE retcode;

   assert(pricetype->getType() == GCG_PRICETYPE_FARKAS || result != NULL);

   assert(nfoundvars != NULL);
   assert(bestredcost != NULL);
   assert(bestredcostvalid != NULL);

   SCIPdebugMessage("%s pricing\n", optimal ? "optimal" : "heuristic");

   solvedmips = 0;
   successfulmips = 0;
   retcode = SCIP_OKAY;
   *nfoundvars = 0;

   *bestredcost = 0.0;
   *bestredcostvalid = ( SCIPgetLPSolstat(scip_) == SCIP_LPSOLSTAT_OPTIMAL && optimal ? TRUE : FALSE );
   pricinglowerbound = SCIPinfinity(scip_);

   maxsols = MAX(MAX(farkaspricing->getMaxvarsround(),reducedcostpricing->getMaxvarsround()),reducedcostpricing->getMaxvarsroundroot());

   SCIP_CALL( SCIPallocMemoryArray(scip_, &nsols, pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip_, &sols, pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip_, &solisray, pricerdata->npricingprobs) );

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip_, &(solisray[i]), maxsols) );
      SCIP_CALL( SCIPallocMemoryArray(scip_, &(sols[i]), maxsols) );
   }

#ifdef _OPENMP
      if( threads > 0 )
         omp_set_num_threads(threads);
#endif

   #pragma omp parallel for private(j)
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {

      for( j = 0; j < maxsols; ++j )
      {
         solisray[i][j] = FALSE;
         sols[i][j] = NULL;
      }
   }

   #pragma omp parallel for ordered firstprivate(pricinglowerbound) shared(retcode, optimal, solisray, sols, nsols, maxsols,pricetype,bestredcost,bestredcostvalid,nfoundvars,successfulmips) reduction(+:solvedmips)
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      int prob;
      SCIP_STATUS status;
      SCIP_RETCODE private_retcode;

      int nvarsfound = *nfoundvars;
      status = SCIP_STATUS_UNKNOWN;
      prob = pricerdata->permu[i];

      if( pricerdata->pricingprobs[prob] == NULL || retcode != SCIP_OKAY)
         goto done;

      if( abortPricing(pricetype, *nfoundvars, solvedmips, successfulmips, optimal) )
         goto done;

      #pragma omp ordered
      {
         private_retcode = solvePricingProblem(prob, pricetype, optimal, &pricinglowerbound, sols[prob], solisray[prob], maxsols, &nsols[prob], &status);

         #pragma omp critical (retcode)
         retcode = private_retcode;

         #pragma omp atomic
         *nfoundvars += countPricedVariables(prob, sols[prob], nsols[prob], solisray[prob] );

         if(nvarsfound < *nfoundvars)
         {
            #pragma omp atomic
            successfulmips += 1;
         }
      }

      if( optimal )
      {
        if( !SCIPisInfinity(scip_, pricinglowerbound) && isPricingOptimal(pricerdata->pricingprobs[prob]) && isMasterLPOptimal() )
            assert(!SCIPisSumPositive(scip_, pricinglowerbound - pricerdata->dualsolconv[prob]));

         #pragma omp atomic
         *bestredcost += GCGrelaxGetNIdenticalBlocks(origprob, prob) * (pricinglowerbound - pricerdata->dualsolconv[prob]);
      }

      #pragma omp atomic
      solvedmips++;

   done:
      ;
   }

   SCIP_CALL( retcode );

#ifdef _OPENMP
   SCIPdebugMessage("We are here with currently %d threads.\n", omp_get_num_threads());
#endif

   /** @todo perhaps solve remaining pricing problems, if only few left? */
   /** @todo solve all pricing problems all k iterations? */
   *nfoundvars = 0;

   for( i = 0; i < pricerdata->npricingprobs; ++i )
   {
      int prob;
      prob = pricerdata->permu[i];
      if( pricerdata->pricingprobs[prob] == NULL )
         continue;

      if( !isPricingOptimal(pricerdata->pricingprobs[prob]) )
      {
         SCIPdebugMessage("Pricing prob %d was not solved to optimality, reduced cost invalid\n", prob);
         *bestredcostvalid = FALSE;
      }

      if( SCIPgetStage(pricerdata->pricingprobs[prob]) <= SCIP_STAGE_PROBLEM )
         continue;

      nfoundvarsprob = 0;

      for( j = 0; j < nsols[prob]; ++j )
      {
         SCIP_VAR** solvars;
         SCIP_Real* solvals;
         int nsolvars;

         /** add variable only if we cannot abort */
         if( (nfoundvarsprob <= pricerdata->maxsolsprob &&
             (pricetype->getType() == GCG_PRICETYPE_REDCOST || *nfoundvars < pricetype->getMaxvarsround()) &&
             (pricetype->getType() == GCG_PRICETYPE_FARKAS || ((*nfoundvars < pricetype->getMaxvarsround() || isRootNode(scip_) ) &&
             (*nfoundvars <reducedcostpricing->getMaxvarsroundroot() || !isRootNode(scip_))))) )
         {
            solvars = SCIPgetOrigVars(pricerdata->pricingprobs[prob]);
            nsolvars = SCIPgetNOrigVars(pricerdata->pricingprobs[prob]);
            SCIP_CALL( SCIPallocMemoryArray(scip, &solvals, nsolvars) );
            SCIPdebugMessage("Solution %d of prob %d (%p)\n", j, prob, sols[prob][j]);
            SCIP_CALL( SCIPgetSolVals(pricerdata->pricingprobs[prob], sols[prob][j], nsolvars, solvars, solvals) );

            /* create new variable, compute objective function value and add it to the master constraints and cuts it belongs to */
            SCIP_CALL( createNewMasterVar(scip_, solvars, solvals, nsolvars, solisray[prob][j], prob,
                  FALSE, &added, NULL) );

            if( added )
            {
               ++(*nfoundvars);
               nfoundvarsprob++;
            }
            SCIPfreeMemoryArray(scip, &solvals);
         }

         if( solisray[prob][j] )
         {
            SCIP_CALL( SCIPfreeSol(pricerdata->pricingprobs[prob], &sols[prob][j]) );
            SCIPdebugMessage("Freeing solution %d of prob %d.\n", j, prob);
         }
      }

      if( nfoundvarsprob > 0 )
         successfulmips++;

   }

   for( i = 0; i < pricerdata->npricingprobs; ++i )
   {
      SCIPdebugMessage("Freeing information for pricing %d.\n", i);
      SCIPfreeMemoryArray(scip, &(sols[i]));
      SCIPfreeMemoryArray(scip, &(solisray[i]));
   }

   SCIPfreeMemoryArray(scip, &solisray);
   SCIPfreeMemoryArray(scip, &sols);
   SCIPfreeMemoryArray(scip, &nsols);

   /* free the pricingproblems if they exist and need to be freed */
   SCIP_CALL( freePricingProblems() );

   return SCIP_OKAY;
}

/** performs the pricing routine, gets the type of pricing that should be done: farkas or redcost pricing */
SCIP_RETCODE ObjPricerGcg::priceNewVariables(
   PricingType*          pricetype,          /**< type of the pricing */
   SCIP_RESULT*          result,             /**< result pointer */
   SCIP_Real*            lowerbound          /**< lowerbound pointer */
   )
{
   int nfoundvars;
   double degeneracy;

   SCIP_Real bestredcost;
   SCIP_Bool bestredcostvalid;

   assert(result != NULL || pricetype->getType() == GCG_PRICETYPE_FARKAS);
   assert(lowerbound != NULL || pricetype->getType() == GCG_PRICETYPE_FARKAS);

   if( lowerbound != NULL )
      *lowerbound = -SCIPinfinity(scip_);

   GCGpricerPrintInfo(scip_, pricerdata, "nvars = %d, current LP objval = %g, time = %f, node = %lld\n",
         SCIPgetNVars(scip_), SCIPgetLPObjval(scip_), SCIPgetSolvingTime(scip_), SCIPgetNNodes(scip_));

   if( pricetype->getType() == GCG_PRICETYPE_REDCOST )
   {
      assert(result != NULL);

      /* terminate early, if applicable */
      if( canPricingBeAborted() )
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }

      *result = SCIP_SUCCESS;
   }

   pricetype->incCalls();

   pricerdata->calls++;
   nfoundvars = 0;

   /* set objectives of the variables in the pricing sub-MIPs */
   SCIP_CALL( setPricingObjs(pricetype) );

   sortPricingProblemsByScore();

   bestredcost = 0.0;
   bestredcostvalid = TRUE;

   if( pricerdata->useheurpricing )
   {
      SCIP_CALL( performPricing(pricetype, FALSE, result, &nfoundvars, &bestredcost, &bestredcostvalid) );
   }

   /* if no variables were found so far, solve the pricing MIPs to optimality and check whether
    * solutions corresponding to variables with negative reduced costs where found
    */
   if( nfoundvars == 0 )
   {
      SCIP_CALL( performPricing(pricetype, TRUE, result, &nfoundvars, &bestredcost, &bestredcostvalid) );
   }

   if( nfoundvars == 0 && isRootNode(scip_) )
   {
      SCIP_CALL( SCIPsetIntParam(scip_, "display/verblevel", 0) );
   }

   if( pricetype->getType() == GCG_PRICETYPE_REDCOST && bestredcostvalid )
   {
      assert(lowerbound != NULL);
      GCGpricerPrintInfo(scip_, pricerdata, "lower bound = %g, bestredcost = %g\n", SCIPgetLPObjval(scip_) + bestredcost, bestredcost);

      *lowerbound = SCIPgetLPObjval(scip_) + bestredcost;
   }

   SCIPdebugMessage("%s pricing: found %d new vars\n", (pricetype->getType() == GCG_PRICETYPE_REDCOST ? "Redcost" : "Farkas"), nfoundvars);

   if( isRootNode(scip_) && pricetype->getType() == GCG_PRICETYPE_REDCOST && pricetype->getCalls() > 0 )
   {
      SCIP_CALL( computeCurrentDegeneracy(&degeneracy) );

      pricerdata->rootnodedegeneracy = degeneracy;

      /* Complicated calculation for numerical stability:
       *     E[\sum_{i=1}^n x_i] = (E[\sum_{i=1}^{n-1} x_i]*(n-1) + x_n)/n
       *     E[\sum_{i=1}^n x_i] = E[\sum_{i=1}^{n-1} x_i]*(n-1)/n + x_n/n
       * <=> E[\sum_{i=1}^n x_i] = E[\sum_{i=1}^{n-1} x_i]-E[\sum_{i=1}^{n-1} x_i]/n + x_n/n
       * <=> E_n = E_{n-1} - E_{n-1}/n + x_n/n
       * <=> E -= E/n - x_n/n
       */
      ++pricerdata->ndegeneracycalcs;
      pricerdata->avgrootnodedegeneracy -= pricerdata->avgrootnodedegeneracy/(pricerdata->ndegeneracycalcs) - degeneracy/(pricerdata->ndegeneracycalcs);
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of variable pricer
 */

 ObjPricerGcg::ObjPricerGcg(
    SCIP*              scip,               /**< SCIP data structure */
    SCIP*              origscip,           /**< SCIP data structure of original problem */
    const char*        name,               /**< name of variable pricer */
    const char*        desc,               /**< description of variable pricer */
    int                priority,           /**< priority of the variable pricer */
    SCIP_Bool          delay,
    SCIP_PRICERDATA*   p_pricerdata
    ) : ObjPricer(scip, name, desc, priority, delay)
 {

    assert(origscip!= NULL);
    pricerdata = p_pricerdata;
    origprob = origscip;
 }

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
SCIP_DECL_PRICERFREE(ObjPricerGcg::scip_free)
{
   assert(scip == scip_);
   SCIP_CALL( solversFree() );

   SCIPfreeMemoryArray(scip, &pricerdata->solvers);

   /* free memory for pricerdata*/
   if( pricerdata != NULL )
   {
      SCIPfreeMemory(scip, &pricerdata);
   }

   delete reducedcostpricing;
   delete farkaspricing;

   SCIPpricerSetData(pricer, NULL);
   return SCIP_OKAY;
}


/** initialization method of variable pricer (called after problem was transformed) */
SCIP_DECL_PRICERINIT(ObjPricerGcg::scip_init)
{
   assert(scip == scip_);
   SCIP_CALL( solversInit() );

   SCIP_CALL( reducedcostpricing->resetCalls() );
   SCIP_CALL( farkaspricing->resetCalls() );

   return SCIP_OKAY;
}


/** deinitialization method of variable pricer (called before transformed problem is freed) */
SCIP_DECL_PRICEREXIT(ObjPricerGcg::scip_exit)
{
   assert(scip == scip_);
   SCIP_CALL( solversExit() );

   return SCIP_OKAY;
}


/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
SCIP_DECL_PRICERINITSOL(ObjPricerGcg::scip_initsol)
{
   int i;
   int norigvars;
   SCIP_Bool discretization;
   SCIP_CONS** masterconss;
   int nmasterconss;
   int origverblevel;

   assert(scip == scip_);
   assert(pricer != NULL);

   /* at the beginning, the output of the master problem gets the same verbosity level
    * as the output of the original problem */
   SCIP_CALL( SCIPgetIntParam(origprob, "display/verblevel", &origverblevel) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", origverblevel) );

   pricerdata->currnodenr = -1;

   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);

   /* init array containing all pricing problems */
   pricerdata->npricingprobs = GCGrelaxGetNPricingprobs(origprob);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->pricingprobs), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->npointsprob), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nraysprob), pricerdata->npricingprobs) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->farkascallsdist), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->farkasfoundvars), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->farkasnodetimedist), pricerdata->npricingprobs) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcostcallsdist), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcostfoundvars), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcostnodetimedist), pricerdata->npricingprobs) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nodetimehist), PRICER_STAT_ARRAYLEN_TIME) ); /*lint !e506*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->foundvarshist), PRICER_STAT_ARRAYLEN_VARS) ); /*lint !e506*/

   BMSclearMemoryArray(pricerdata->nodetimehist, PRICER_STAT_ARRAYLEN_TIME);
   BMSclearMemoryArray(pricerdata->foundvarshist, PRICER_STAT_ARRAYLEN_VARS);

   pricerdata->oldvars=0;

   pricerdata->npricingprobsnotnull = 0;

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {

      pricerdata->farkascallsdist[i] = 0;
      pricerdata->farkasfoundvars[i] = 0;
      pricerdata->farkasnodetimedist[i] = 0;
      pricerdata->redcostcallsdist[i] = 0;
      pricerdata->redcostfoundvars[i] = 0;
      pricerdata->redcostnodetimedist[i]= 0;


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
   norigvars = SCIPgetNOrigVars(origprob);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->solvals), norigvars) );

   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->freeclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->transformclock)) );

   pricerdata->solvedsubmipsoptimal = 0;
   pricerdata->solvedsubmipsheur = 0;
   pricerdata->calls = 0;
   pricerdata->pricingiters = 0;

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

   SCIP_CALL( SCIPhashmapCreate(&(pricerdata->mapcons2idx), SCIPblkmem(scip), 10 * nmasterconss +1) );
   for( i = 0; i < nmasterconss; i++ )
   {
      SCIP_CALL( SCIPhashmapInsert(pricerdata->mapcons2idx, masterconss[i], (void*)(size_t)i) );
      assert((int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, masterconss[i]) == i); /*lint !e507*/
   }

   pricerdata->npricedvars = 0;
   pricerdata->maxpricedvars = 50;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &pricerdata->pricedvars, pricerdata->maxpricedvars) );

   pricerdata->rootnodedegeneracy = 0.0;
   pricerdata->avgrootnodedegeneracy = 0.0;
   pricerdata->ndegeneracycalcs = 0;

   SCIP_CALL( solversInitsol() );

   return SCIP_OKAY;
}


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
SCIP_DECL_PRICEREXITSOL(ObjPricerGcg::scip_exitsol)
{
   int i;

   assert(scip == scip_);
   assert(pricer != NULL);
   assert(pricerdata != NULL);

   SCIPhashmapFree(&(pricerdata->mapcons2idx));

   SCIPfreeMemoryArray(scip, &(pricerdata->pricingprobs));
   SCIPfreeMemoryArray(scip, &(pricerdata->dualsolconv));
   SCIPfreeMemoryArray(scip, &(pricerdata->score));
   SCIPfreeMemoryArray(scip, &(pricerdata->permu));
   SCIPfreeMemoryArray(scip, &(pricerdata->solvals));
   SCIPfreeMemoryArray(scip, &(pricerdata->npointsprob));
   SCIPfreeMemoryArray(scip, &(pricerdata->nraysprob));

   SCIPfreeMemoryArray(scip, &(pricerdata->farkascallsdist));
   SCIPfreeMemoryArray(scip, &(pricerdata->farkasfoundvars));
   SCIPfreeMemoryArray(scip, &(pricerdata->farkasnodetimedist));

   SCIPfreeMemoryArray(scip, &(pricerdata->redcostcallsdist));
   SCIPfreeMemoryArray(scip, &(pricerdata->redcostfoundvars));
   SCIPfreeMemoryArray(scip, &(pricerdata->redcostnodetimedist));

   SCIPfreeMemoryArray(scip, &(pricerdata->nodetimehist));
   SCIPfreeMemoryArray(scip, &(pricerdata->foundvarshist));

   pricerdata->nodetimehist = NULL;
   pricerdata->foundvarshist = NULL;

   for( i = 0; i < pricerdata->npricedvars; i++ )
   {
      SCIP_CALL( SCIPdropVarEvent(scip, pricerdata->pricedvars[i], SCIP_EVENTTYPE_VARDELETED,
            pricerdata->eventhdlr, NULL, -1) );

      SCIP_CALL( SCIPreleaseVar(scip, &pricerdata->pricedvars[i]) );
   }
   SCIPfreeBlockMemoryArray(scip, &pricerdata->pricedvars, pricerdata->maxpricedvars);

   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->freeclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->transformclock)) );

   SCIP_CALL( solversExitsol() );

   return SCIP_OKAY;
}


/** reduced cost pricing method of variable pricer for feasible LPs */
SCIP_DECL_PRICERREDCOST(ObjPricerGcg::scip_redcost)
{
   SCIP_RETCODE retcode;

   assert(scip == scip_);
   assert(pricer != NULL);
   assert(pricerdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( reducedcostpricing->getCalls() == 0 )
   {
      /** @todo This is just a workaround around SCIP stages! */
      if( farkaspricing->getCalls() == 0 )
      {
         SCIP_CALL( SCIPconsMasterbranchAddRootCons(scip) );
      }
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Starting reduced cost pricing...\n");
   }
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
   if( reducedcostpricing->getMaxrounds() == 0 || (pricerdata->nroundsredcost >= reducedcostpricing->getMaxrounds() && pricerdata->currnodenr != 1) )
   {
      SCIPdebugMessage("pricing aborted at node %lld\n", pricerdata->currnodenr);
      return SCIP_OKAY;
   }

   *result = SCIP_SUCCESS;

   /* perform pricing */
   SCIP_CALL( reducedcostpricing->startClock() );
   retcode = priceNewVariables(reducedcostpricing, result, lowerbound);
   SCIP_CALL( reducedcostpricing->startClock() );

   return retcode;
}

/** farcas pricing method of variable pricer for infeasible LPs */
SCIP_DECL_PRICERFARKAS(ObjPricerGcg::scip_farkas)
{
   SCIP_RETCODE retcode;
   SCIP_SOL** origsols;
   int norigsols;
   int i;

   assert(scip == scip_);
   assert(pricer != NULL);
   assert(pricerdata != NULL);

   /** @todo This is just a workaround around SCIP stages! */
   if(reducedcostpricing->getCalls() == 0  && farkaspricing->getCalls() == 0 )
   {
      SCIP_CALL( SCIPconsMasterbranchAddRootCons(scip) );
   }

   /* get solutions from the original problem */
   origsols = SCIPgetSols(origprob);
   norigsols = SCIPgetNSols(origprob);
   assert(norigsols >= 0);

   /* Add already known solutions for the original problem to the master variable space */
   /** @todo This is just a workaround around SCIP stages! */
   if( farkaspricing->getCalls() == 0 )
   {
      for( i = 0; i < norigsols; ++i )
      {
         assert(origsols[i] != NULL);
         SCIP_CALL( GCGpricerTransOrigSolToMasterVars(scip, origsols[i]) );
      }
      /* return if we transferred solutions as the master should be feasible */
      if( norigsols > 0 )
      {
         farkaspricing->incCalls();
         return SCIP_OKAY;
      }
   }

   SCIP_CALL( farkaspricing->startClock() );
   retcode = priceNewVariables(farkaspricing, NULL, NULL);
   SCIP_CALL( farkaspricing->stopClock() );

   return retcode;
}

void ObjPricerGcg::createPricingTypes()
{
   farkaspricing = new FarkasPricing(scip_);
   reducedcostpricing = new ReducedCostPricing(scip_);
}

/*
 * variable pricer specific interface methods
 */

/** creates the GCG variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origprob            /**< SCIP data structure of the original problem */
   )
{
   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;

   SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );

   /* initialize solvers array */
   pricerdata->solvers = NULL;
   pricerdata->nsolvers = 0;
   pricerdata->nodetimehist = NULL;
   pricerdata->foundvarshist = NULL;

   pricer = new ObjPricerGcg(scip, origprob, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY, pricerdata);
   /* include variable pricer */
   SCIP_CALL( SCIPincludeObjPricer(scip, pricer, TRUE) );

   pricer->createPricingTypes();

   /* include event handler into master SCIP */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, eventFreeVardeleted, eventInitVardeleted, eventExitVardeleted,
         eventInitsolVardeleted, eventExitsolVardeleted, eventDeleteVardeleted, eventExecVardeleted,
         NULL) );

   pricerdata->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/maxsolsprob",
         "maximal number of variables added for each block in a pricinground",
         &pricerdata->maxsolsprob, FALSE, DEFAULT_MAXSOLSPROB, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "pricing/masterpricer/useheurpricing",
         "should pricing be performed heuristically before solving the MIPs to optimality?",
         &pricerdata->useheurpricing, TRUE, DEFAULT_USEHEURPRICING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "pricing/masterpricer/abortpricingint",
         "should pricing be aborted due to integral objective function?",
         &pricerdata->abortpricingint, TRUE, DEFAULT_ABORTPRICINGINT, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricing/masterpricer/abortpricinggap",
         "should pricing be aborted due to small gap between dual bound and RMP objective?",
         &pricerdata->abortpricinggap, TRUE, DEFAULT_ABORTPRICINGGAP, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricing/masterpricer/successfulsubmipsrel",
         "part of the submips that are solved and lead to new variables before pricing round is aborted? (1.0 = solve all pricing MIPs)",
         &pricerdata->successfulmipsrel, FALSE, DEFAULT_SUCCESSFULMIPSREL, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "pricing/masterpricer/dispinfos",
         "should additional informations concerning the pricing process be displayed?",
         &pricerdata->dispinfos, FALSE, DEFAULT_DISPINFOS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/sorting",
         "which sorting method should be used to sort the pricing problems (0 = order of pricing problems, 1 = according to dual solution of convexity constraint, 2 = according to reliability from previous round)",
         &pricerdata->sorting, FALSE, DEFAULT_SORTING, 0, 5, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/threads",
         "how many threads should be used to concurrently solve the pricing problem (0 to guess threads by OpenMP)",
         &ObjPricerGcg::threads, FALSE, DEFAULT_THREADS, 0, 4096, NULL, NULL) );

   return SCIP_OKAY;
}

/** returns the pointer to the scip instance representing the original problem */
SCIP* GCGpricerGetOrigprob(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   ObjPricerGcg* pricer;

   assert(scip != NULL);

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   return pricer->getOrigprob();
}

/** returns the array of variables that were priced in during the solving process */
SCIP_VAR** GCGpricerGetPricedvars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
   assert(pricerdata != NULL);

   return pricerdata->pricedvars;
}

/** returns the number of variables that were priced in during the solving process */
int GCGpricerGetNPricedvars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
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
   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(pos >= 0);

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
   assert(pricerdata != NULL);

   SCIP_CALL( SCIPhashmapInsert(pricerdata->mapcons2idx, cons, (void*)(size_t)pos) );
   assert((int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, cons) == pos); /*lint !e507*/

   SCIPdebugMessage("Added cons %s (%p) to hashmap with index %d\n", SCIPconsGetName(cons), cons, pos);

   return SCIP_OKAY;
}

/** includes a solver into the pricer data */
SCIP_RETCODE GCGpricerIncludeSolver(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of solver */
   const char*           description,        /**< description of solver */
   int                   priority,           /**< priority of solver */
   GCG_DECL_SOLVERSOLVE  ((*solversolve)),   /**< solving method for solver */
   GCG_DECL_SOLVERSOLVEHEUR((*solveheur)),   /**< heuristic solving method for solver */
   GCG_DECL_SOLVERFREE   ((*solverfree)),    /**< free method of solver */
   GCG_DECL_SOLVERINIT   ((*solverinit)),    /**< init method of solver */
   GCG_DECL_SOLVEREXIT   ((*solverexit)),    /**< exit method of solver */
   GCG_DECL_SOLVERINITSOL((*solverinitsol)), /**< initsol method of solver */
   GCG_DECL_SOLVEREXITSOL((*solverexitsol)), /**< exitsol method of solver */
   GCG_SOLVERDATA*       solverdata          /**< solverdata data structure */
   )
{
   int pos;

   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
   assert(pricerdata != NULL);


   SCIP_CALL( pricer->ensureSizeSolvers() );

   /* solvers array is sorted decreasingly wrt. the priority, find right position and shift solvers with smaller priority */
   pos = pricerdata->nsolvers;
   while( pos >= 1 && pricerdata->solvers[pos-1]->priority < priority )
   {
      pricerdata->solvers[pos] = pricerdata->solvers[pos-1];
      pos--;
   }
   SCIP_CALL( SCIPallocMemory(scip, &(pricerdata->solvers[pos])) ); /*lint !e866*/

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

/** returns the solverdata of a solver */
GCG_SOLVERDATA* GCGsolverGetSolverdata(
   GCG_SOLVER*           solver              /**< pointer so solver */
   )
{
   assert(solver != NULL);

   return solver->solverdata;
}

/** sets solver data of specific solver */
void GCGsolverSetSolverdata(
   GCG_SOLVER*           solver,             /**< pointer to solver  */
   GCG_SOLVERDATA*       solverdata          /**< solverdata data structure */
   )
{
   assert(solver != NULL);

   solver->solverdata = solverdata;
}

/** writes out a list of all pricing problem solvers */
void GCGpricerPrintListOfSolvers(
   SCIP*                 scip                /**< SCIP data structure */
   )
{

   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;
   int nsolvers;
   int i;

   assert(scip != NULL);

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));

   nsolvers = pricerdata->nsolvers;

   SCIPdialogMessage(scip, NULL, " solver               priority description\n --------------       -------- -----------\n");

   for( i = 0; i < nsolvers; ++i )
   {
      SCIPdialogMessage(scip, NULL,  " %-20s", pricerdata->solvers[i]->name);
      SCIPdialogMessage(scip, NULL,  " %8d", pricerdata->solvers[i]->priority);
      SCIPdialogMessage(scip, NULL,  " %s\n", pricerdata->solvers[i]->description);
   }
}

/** prints pricer statistics */
void GCGpricerPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;
   int i;
   double start;
   double end;

   assert(scip != NULL);

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
   assert(pricerdata != NULL);

   /**@todo add constraint statistics: how many constraints (instead of cuts) have been added? */
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Pricing Solver     : #HeurFarkas  #OptFarkas  #HeurRedcost #OptRedcost Time: HeurFarkas  OptFarkas  HeurRedcost OptRedcost\n");

   for( i = 0; i < pricerdata->nsolvers; ++i )
   {
      GCG_SOLVER* solver;
      solver = pricerdata->solvers[i];
      assert(solver != NULL);

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  %-17.17s:", solver->name);
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " %11d %11d   %11d %11d       %10.2f %10.2f   %10.2f %10.2f \n",
         solver->heurfarkascalls, solver->optfarkascalls,
         solver->heurredcostcalls, solver->optredcostcalls,
         SCIPgetClockTime(scip, solver->heurfarkasclock),
         SCIPgetClockTime(scip, solver->optfarkasclock),
         SCIPgetClockTime(scip, solver->heurredcostclock),
         SCIPgetClockTime(scip, solver->optredcostclock));
   }

   /* print of Pricing Statistics */

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Farkas pricing Statistic:\nno.\t#Calls\t\t#Vars\t\ttime(s)\n");

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d  \t %d \t\t %d \t\t %.2f \n", i, pricerdata->farkascallsdist[i],
         pricerdata->farkasfoundvars[i], pricerdata->farkasnodetimedist[i]);

   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Reduced Cost pricing Statistic:\nno.\t#Calls\t\t#Vars\t\ttime(s)\n");

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d  \t %d \t\t %d \t\t %.2f \n", i, pricerdata->redcostcallsdist[i],
         pricerdata->redcostfoundvars[i], pricerdata->redcostnodetimedist[i]);

   }

   /* print of Histogram Buckets !=0      */

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Histogram Time\n");
   for( i = 0; i < PRICER_STAT_ARRAYLEN_TIME; i++ )
   {
      start = (1.0 * i * PRICER_STAT_BUCKETSIZE_TIME)/1000.0;
      end = start + PRICER_STAT_BUCKETSIZE_TIME/1000.0;

      if( pricerdata->nodetimehist[i] != 0 )
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "From\t%.4f\t-\t%.4f\ts:\t\t%d \n", start, end, pricerdata->nodetimehist[i]);
   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Histogram Found Vars\n");

   for( i = 0; i < PRICER_STAT_ARRAYLEN_VARS; i++ )
   {
      start = i * PRICER_STAT_BUCKETSIZE_VARS;
      end = start + PRICER_STAT_BUCKETSIZE_VARS;

      if( pricerdata->foundvarshist[i] != 0 )
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "From\t%.0f\t-\t%.0f\tvars:\t\t%d \n", start, end, pricerdata->foundvarshist[i]);
   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Pricing Summary:\n");
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Calls                            : %d\n", pricerdata->calls);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Farkas Pricing Calls             : %d\n", pricer->getFarkasPricing()->getCalls());
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Farkas Pricing Time              : %f\n", pricer->getFarkasPricing()->getClockTime());
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Reduced Cost Pricing Calls       : %d\n", pricer->getReducedCostPricing()->getCalls());
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Reduced Cost Pricing Time        : %f\n", pricer->getReducedCostPricing()->getClockTime());
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Solved subMIPs Heuristic Pricing : %d\n", pricerdata->solvedsubmipsheur);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Solved subMIPs Optimal Pricing   : %d\n", pricerdata->solvedsubmipsoptimal);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Time for transformation          : %f\n", SCIPgetClockTime(scip, pricerdata->transformclock));
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Time for freeing subMIPs         : %f\n", SCIPgetClockTime(scip, pricerdata->freeclock));

}


/** transfers a primal solution of the original problem into the master variable space,
 *  i.e. creates one master variable for each block and adds the solution to the master problem  */
SCIP_RETCODE GCGpricerTransOrigSolToMasterVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             origsol             /**< the solution that should be transferred */
   )
{
   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_SOL* mastersol;
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

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
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
      if( pricerdata->pricingprobs[i] == NULL )
         continue;

      SCIP_CALL( SCIPallocBufferArray(scip, &pricingvars[i], SCIPgetNVars(pricerdata->pricingprobs[i])) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &pricingvals[i], SCIPgetNVars(pricerdata->pricingprobs[i])) ); /*lint !e866*/
   }

   /* get solution values */
   SCIP_CALL( SCIPgetSolVals(scip, origsol, norigvars, origvars, origsolvals) );
   SCIP_CALL( SCIPcreateSol(scip, &mastersol, SCIPgetSolHeur(origprob, origsol)) );

   /* store variables and solutions into arrays */
   for( i = 0; i < norigvars; i++ )
   {
      int blocknr;
      assert(GCGvarIsOriginal(origvars[i]));
      blocknr = GCGvarGetBlock(origvars[i]);
      assert(blocknr < 0 || GCGoriginalVarGetPricingVar(origvars[i]) != NULL);

      if( blocknr >= 0 )
      {
         prob = blocknr;
         if( pricerdata->pricingprobs[prob] == NULL )
            continue;

         if( !SCIPisZero(scip, origsolvals[i]) )
         {
            pricingvars[prob][npricingvars[prob]] = GCGoriginalVarGetPricingVar(origvars[i]);
            pricingvals[prob][npricingvars[prob]] = origsolvals[i];
            npricingvars[prob]++;
         }
      }
      else
      {
         assert((GCGoriginalVarGetNMastervars(origvars[i]) == 1) || (GCGvarIsLinking(origvars[i])));
         assert(GCGoriginalVarGetMastervars(origvars[i])[0] != NULL);
         SCIP_CALL( SCIPsetSolVal(scip, mastersol, GCGoriginalVarGetMastervars(origvars[i])[0], origsolvals[i]) );
      }
   }

   /* create variables in the master problem */
   for( prob = 0; prob < pricerdata->npricingprobs; prob++ )
   {
      if( pricerdata->pricingprobs[prob] == NULL )
      {
         continue;
      }
      SCIP_CALL( pricer->createNewMasterVar(scip, pricingvars[prob], pricingvals[prob], npricingvars[prob], FALSE, prob, TRUE, &added, &newvar) );
      assert(added);

      SCIP_CALL( SCIPsetSolVal(scip, mastersol, newvar, 1.0 * GCGrelaxGetNIdenticalBlocks(origprob, prob)) );
   }

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPtrySolFree(scip, &mastersol, TRUE, TRUE, TRUE, TRUE, &added) );
#else
   SCIP_CALL( SCIPtrySolFree(scip, &mastersol, FALSE, TRUE, TRUE, TRUE, &added) );
#endif

   /* free memory for storing variables and solution values from the solution */

   for( i = pricerdata->npricingprobs - 1; i>= 0; i-- )
   {
      if( pricerdata->pricingprobs[i] == NULL )
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


/** create initial master variables */
SCIP_RETCODE GCGpricerCreateInitialMastervars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;
   int i;
   SCIP* origprob;
   SCIP_VAR** vars;
   int nvars;
   int npricingprobs;
   int v;

   assert(scip != NULL);

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
   assert(pricerdata != NULL);

   origprob = pricer->getOrigprob();
   assert(origprob != NULL);

   npricingprobs = GCGrelaxGetNPricingprobs(origprob);
   assert(npricingprobs >= 0);

   /* for variables in the original problem that do not belong to any block,
    * create the corresponding variable in the master problem
    */
   vars = SCIPgetVars(origprob);
   nvars = SCIPgetNVars(origprob);
   for( v = 0; v < nvars; v++ )
   {
      SCIP_Real* coefs;
      int blocknr;
      int ncoefs;
      SCIP_VAR* var;

      /* var = SCIPvarGetProbvar(vars[v]); */
      var = vars[v];
      blocknr = GCGvarGetBlock(var);
      coefs = GCGoriginalVarGetCoefs(var);
      ncoefs = GCGoriginalVarGetNCoefs(var);

      assert(GCGvarIsOriginal(var));
      if( blocknr < 0 )
      {
         SCIP_CONS** linkconss;
         SCIP_VAR* newvar;

         SCIP_CALL( GCGcreateInitialMasterVar(scip, var, &newvar) );
         SCIP_CALL( SCIPaddVar(scip, newvar) );

         SCIP_CALL( GCGoriginalVarAddMasterVar(scip, var, newvar, 1.0) );

         linkconss = GCGoriginalVarGetMasterconss(var);

         /* add variable in the master to the master constraints it belongs to */
         for( i = 0; i < ncoefs; i++ )
         {
            assert(!SCIPisZero(scip, coefs[i]));
            /*            SCIP_CALL( SCIPgetTransformedCons(scip, linkconss[i], &linkcons) );*/

            SCIP_CALL( SCIPaddCoefLinear(scip, linkconss[i], newvar, coefs[i]) );
         }

         /* we copied a linking variable into the master, add it to the linkcons */
         if( GCGvarIsLinking(var) )
         {
            SCIP_CONS** linkingconss;
            linkingconss = GCGlinkingVarGetLinkingConss(var);

            for( i = 0; i < npricingprobs; i++ )
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
   return SCIP_OKAY;
}

/** get root node degeneracy */
SCIP_Real GCGpricerGetDegeneracy(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
   assert(pricerdata != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) <= SCIP_STAGE_SOLVING && isRootNode(scip) )
   {
      return pricerdata->avgrootnodedegeneracy;
   }
   else
      return SCIPinfinity(scip);
}

/* get number of iterations in pricing problems */
SCIP_Longint GCGpricerGetPricingSimplexIters(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
   assert(pricerdata != NULL);

   return pricerdata->pricingiters;
}

/** print simplex iteration statistics */
extern
SCIP_RETCODE GCGpricerPrintSimplexIters(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   ObjPricerGcg* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = static_cast<ObjPricerGcg*>(SCIPfindObjPricer(scip, PRICER_NAME));
   assert(pricer != NULL);

   pricerdata = pricer->getPricerdata();
   assert(pricerdata != NULL);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Simplex iterations :       iter\n");
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  Master LP        : %10lld\n", SCIPgetNLPIterations(scip));
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  Pricing LP       : %10lld\n", pricerdata->pricingiters);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  Original LP      : %10lld\n", SCIPgetNLPIterations(pricer->getOrigprob()));

   return SCIP_OKAY;
}