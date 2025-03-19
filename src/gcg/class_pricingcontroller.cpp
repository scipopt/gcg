/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   class_pricingcontroller.cpp
 * @brief  pricing controller managing the pricing strategy
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/class_pricingcontroller.h"
#include "gcg/pricer_gcg.h"
#include "gcg/class_pricingtype.h"
#include "gcg/gcg.h"
#include "gcg/scip_misc.h"
#include "gcg/branch_generic.h"
#include "gcg/cons_masterbranch.h"
#include "gcg/pub_gcgpqueue.h"
#include "gcg/pub_pricingjob.h"
#include "gcg/pub_pricingprob.h"
#include "gcg/pub_solver.h"
#include "gcg/pricingjob.h"
#include "gcg/pricingprob.h"

#include "scip/scip.h"
#include "objscip/objscip.h"

#ifdef _OPENMP
#include "gcg/struct_locks.h"
#endif

#include <exception>

#define DEFAULT_HEURPRICINGITERS         1          /**< maximum number of heuristic pricing iterations per pricing call and problem */
#define DEFAULT_MAXHEURDEPTH             -1         /**< maximum depth at which heuristic pricing should be performed (-1 for infinity) */
#define DEFAULT_SORTING                  'r'        /**< order by which the pricing problems should be sorted:
                                                     *    'i'ndices
                                                     *    'd'ual solutions of convexity constraints
                                                     *    'r'eliability from all previous rounds
                                                     *    reliability from the 'l'ast nroundscol rounds
                                                     */
#define DEFAULT_NROUNDSCOL               15
#define DEFAULT_CHUNKSIZE                INT_MAX    /**< maximal number of pricing problems to be solved during one pricing loop */
#define DEFAULT_EAGERFREQ                10         /**< frequency at which all pricingproblems should be solved (0 to disable) */

#define SCIP_CALL_EXC(x)   do                                                                                 \
                       {                                                                                      \
                          SCIP_RETCODE _retcode_;                                                             \
                          if( (_retcode_ = (x)) !=  SCIP_OKAY )                                               \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _retcode_);                    \
                             throw std::exception();                                                          \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )


namespace gcg {

Pricingcontroller::Pricingcontroller(
   GCG*                   gcgstruct
   )
{
   gcg = gcgstruct;
   scip_ = GCGgetMasterprob(gcg);
   origprob = GCGgetOrigprob(gcg);
   pricingprobs = NULL;
   npricingprobs = 0;
   pricingjobs = NULL;
   npricingjobs = 0;

   sorting = DEFAULT_SORTING;
   nroundscol = DEFAULT_NROUNDSCOL;
   chunksize = DEFAULT_CHUNKSIZE;
   eagerfreq = DEFAULT_EAGERFREQ;

   pqueue = NULL;
   maxniters = INT_MAX;
   nchunks = 1;
   curchunk = 0;

   pricingtype_ = NULL;

   eagerage = 0;
   nsolvedprobs = 0;

#ifdef _OPENMP
   locks = GCGgetLocks(gcg);
#endif
}

Pricingcontroller::~Pricingcontroller()
{
}

SCIP_RETCODE Pricingcontroller::addParameters()
{
   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/heurpricingiters",
         "maximum number of heuristic pricing iterations per pricing call and problem",
         &heurpricingiters, FALSE, DEFAULT_HEURPRICINGITERS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/maxheurdepth",
         "maximum depth at which heuristic pricing should be performed (-1 for infinity)",
         &maxheurdepth, FALSE, DEFAULT_MAXHEURDEPTH, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(origprob, "pricing/masterpricer/sorting",
         "order by which the pricing problems should be sorted ('i'ndices, 'd'ual solutions of convexity constraints, 'r'eliability from previous rounds, reliability from the 'l'ast nroundscol rounds)",
         &sorting, FALSE, DEFAULT_SORTING, "dilr", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/nroundscol",
         "number of previous pricing rounds for which the number of improving columns should be counted",
         &nroundscol, TRUE, DEFAULT_NROUNDSCOL, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/chunksize",
         "maximal number of pricing problems to be solved during one pricing loop",
         &chunksize, TRUE, DEFAULT_CHUNKSIZE, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/eagerfreq",
         "frequency at which all pricingproblems should be solved (0 to disable)",
         &eagerfreq, FALSE, DEFAULT_EAGERFREQ, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** comparison operator for pricing jobs w.r.t. their solution priority */
SCIP_DECL_SORTPTRCOMP(Pricingcontroller::comparePricingjobs)
{
   GCG_PRICINGJOB* pricingjob1;
   GCG_PRICINGJOB* pricingjob2;
   GCG_PRICINGPROB* pricingprob1;
   GCG_PRICINGPROB* pricingprob2;
   GCG_SOLVER* solver1;
   GCG_SOLVER* solver2;

   pricingjob1 = (GCG_PRICINGJOB*) elem1;
   pricingjob2 = (GCG_PRICINGJOB*) elem2;

   pricingprob1 = GCGpricingjobGetPricingprob(pricingjob1);
   pricingprob2 = GCGpricingjobGetPricingprob(pricingjob2);

   solver1 = GCGpricingjobGetSolver(pricingjob1);
   solver2 = GCGpricingjobGetSolver(pricingjob2);

   /* preliminary order of sorting:
    *  * priority of pricing solvers
    *  * heuristic before exact
    *  * score
    */

   if( GCGsolverGetPriority(solver1) > GCGsolverGetPriority(solver2) )
      return -1;
   else if( GCGsolverGetPriority(solver1) < GCGsolverGetPriority(solver2) )
      return 1;

   if( GCGpricingjobIsHeuristic(pricingjob1) && GCGpricingjobIsHeuristic(pricingjob2) )
   {
      if( GCGpricingjobGetNHeurIters(pricingjob1) < GCGpricingjobGetNHeurIters(pricingjob2) )
         return -1;
      else if( GCGpricingjobGetNHeurIters(pricingjob1) > GCGpricingjobGetNHeurIters(pricingjob2) )
         return 1;
   }

   if( GCGpricingjobIsHeuristic(pricingjob1) != GCGpricingjobIsHeuristic(pricingjob2) )
   {
      if( GCGpricingjobIsHeuristic(pricingjob1) )
         return -1;
      else if( GCGpricingjobIsHeuristic(pricingjob2) )
         return 1;
   }

   if( GCGpricingjobGetScore(pricingjob1) > GCGpricingjobGetScore(pricingjob2) )
      return -1;
   else if( GCGpricingjobGetScore(pricingjob1) < GCGpricingjobGetScore(pricingjob2) )
      return 1;

   // @todo: preliminary tie breaking by pricing problem index
   if( GCGpricingprobGetProbnr(pricingprob1) < GCGpricingprobGetProbnr(pricingprob2) )
      return -1;
   else
      return 1;

   return 0;
}

/** for each pricing problem, get its corresponding generic branching constraints */
SCIP_RETCODE Pricingcontroller::getGenericBranchconss()
{
   /* get current branching rule */
   SCIP_CONS* branchcons = GCGconsMasterbranchGetActiveCons(gcg);
   SCIP_BRANCHRULE* branchrule = GCGconsMasterbranchGetBranchrule(branchcons);

   assert(branchcons != NULL);
   assert(SCIPnodeGetDepth(GCGconsMasterbranchGetNode(branchcons)) == 0 || branchrule != NULL || SCIPinProbing(scip_));

   while( GCGisBranchruleGeneric(branchrule) )
   {
      GCG_BRANCHDATA* branchdata;
      SCIP_CONS* mastercons;
      int consblocknr;

      int i;

      branchdata = GCGconsMasterbranchGetBranchdata(branchcons);
      assert(branchdata != NULL);

      mastercons = GCGbranchGenericBranchdataGetMastercons(branchdata);
      consblocknr = GCGbranchGenericBranchdataGetConsblocknr(branchdata);
      assert(mastercons != NULL);
      assert(consblocknr >= 0 || consblocknr == -3);

      if( consblocknr >= 0 )
      {
         for( i = 0; i < npricingprobs; ++i )
         {
            /* search for the pricing problem to which the generic branching decision belongs */
            if( consblocknr == GCGpricingprobGetProbnr(pricingprobs[i]) )
            {
               SCIP_CALL( GCGpricingprobAddGenericBranchData(gcg, pricingprobs[i], branchcons,
                  pricingtype_->consGetDual(mastercons)) );
               break;
            }
         }
         assert(i < npricingprobs);
      }

      branchcons = GCGconsMasterbranchGetParentcons(branchcons);
      branchrule = GCGconsMasterbranchGetBranchrule(branchcons);
   }

   return SCIP_OKAY;
}

/** check if a pricing problem is done */
SCIP_Bool Pricingcontroller::pricingprobIsDone(
   GCG_PRICINGPROB*      pricingprob        /**< pricing problem structure */
   ) const
{
   return GCGpricingprobGetNImpCols(pricingprob) > 0
      || (GCGpricingprobGetStatus(pricingprob) == GCG_PRICINGSTATUS_OPTIMAL && GCGpricingprobGetBranchconsIdx(pricingprob) == 0)
      || GCGpricingprobGetStatus(pricingprob) == GCG_PRICINGSTATUS_INFEASIBLE
      || GCGpricingprobGetStatus(pricingprob) == GCG_PRICINGSTATUS_UNBOUNDED;
}

/** check whether the next generic branching constraint of a pricing problem must be considered */
SCIP_Bool Pricingcontroller::pricingprobNeedsNextBranchingcons(
   GCG_PRICINGPROB*      pricingprob        /**< pricing problem structure */
   ) const
{
   return GCGpricingprobGetNImpCols(pricingprob) == 0
      && GCGpricingprobGetStatus(pricingprob) == GCG_PRICINGSTATUS_OPTIMAL
      && GCGpricingprobGetBranchconsIdx(pricingprob) > 0;
}

SCIP_RETCODE Pricingcontroller::initSol()
{
   int nblocks = GCGgetNPricingprobs(gcg);
   GCG_SOLVER** solvers = GCGpricerGetSolvers(gcg);
   int nsolvers = GCGpricerGetNSolvers(gcg);
   int actchunksize = MIN(chunksize, GCGgetNRelPricingprobs(gcg));

   npricingprobs = 0;
   npricingjobs = 0;
   nchunks = (int) SCIPceil(scip_, (SCIP_Real) GCGgetNRelPricingprobs(gcg) / actchunksize);
   curchunk = nchunks - 1;
   eagerage = 0;

   /* create pricing problem and pricing job data structures */
   maxpricingprobs = SCIPcalcMemGrowSize(scip_, GCGgetNRelPricingprobs(gcg));
   maxpricingjobs = SCIPcalcMemGrowSize(scip_, GCGgetNRelPricingprobs(gcg) * nsolvers);
   SCIP_CALL_EXC( SCIPallocBlockMemoryArray(scip_, &pricingprobs, maxpricingprobs) );
   SCIP_CALL_EXC( SCIPallocBlockMemoryArray(scip_, &pricingjobs, maxpricingjobs) );
   for( int i = 0; i < nblocks; ++i )
   {
      if( GCGisPricingprobRelevant(gcg, i) )
      {
         SCIP_CALL_EXC( GCGpricingprobCreate(gcg, &pricingprobs[npricingprobs], GCGgetPricingprob(gcg, i), i, nroundscol) );

         for( int j = 0; j < nsolvers; ++j )
         {
            if( GCGsolverIsHeurEnabled(solvers[j]) || GCGsolverIsExactEnabled(solvers[j]) )
            {
               SCIP_CALL_EXC( GCGpricingjobCreate(gcg, &pricingjobs[npricingjobs], pricingprobs[npricingprobs], solvers[j], npricingprobs / actchunksize) );
               ++npricingjobs;
               break;
            }
         }
         ++npricingprobs;
      }
   }

#ifdef _OPENMP
   SCIP_CALL_EXC( GCGpqueueCreate(scip_, &pqueue, npricingjobs, comparePricingjobs, &locks->memorylock) );
#else
   SCIP_CALL_EXC( GCGpqueueCreate(scip_, &pqueue, npricingjobs, comparePricingjobs, NULL) );
#endif

   return SCIP_OKAY;
}

SCIP_RETCODE Pricingcontroller::exitSol()
{
   GCGpqueueFree(&pqueue);

   for( int i = 0; i < npricingprobs; ++i )
   {
      GCGpricingprobFree(gcg, &pricingprobs[i]);
   }
   for( int i = 0; i < npricingjobs; ++i )
   {
      GCGpricingjobFree(gcg, &pricingjobs[i]);
   }
   SCIPfreeBlockMemoryArray(scip_, &pricingprobs, maxpricingprobs);
   SCIPfreeBlockMemoryArray(scip_, &pricingjobs, maxpricingjobs);

   return SCIP_OKAY;
}

/** pricing initialization, called right at the beginning of pricing */
SCIP_RETCODE Pricingcontroller::initPricing(
   PricingType*          pricingtype         /**< type of pricing */
   )
{
   SCIP_Longint tmpmaxniters;

   pricingtype_ = pricingtype;

   /* move chunk index forward */
   curchunk = (curchunk + 1) % nchunks;
   startchunk = curchunk;

   for( int i = 0; i < npricingprobs; ++i )
      GCGpricingprobInitPricing(pricingprobs[i]);

   SCIP_CALL( getGenericBranchconss() );

   /* calculate maximal possible number of pricing iterations per mis-pricing iteration */
   tmpmaxniters = 0;
   for( int i = 0; i < npricingprobs; ++i )
      tmpmaxniters += GCGpricerGetNSolvers(gcg) * ((SCIP_Longint) heurpricingiters + 1) * (GCGpricingprobGetNGenericBranchconss(pricingprobs[i]) + 1);
   maxniters = (int) MIN(tmpmaxniters, 16383);

   SCIPdebugMessage("initialize pricing, chunk = %d/%d\n", curchunk+1, nchunks);

   return SCIP_OKAY;
}

/** pricing deinitialization, called when pricing is finished */
void Pricingcontroller::exitPricing()
{
   for( int i = npricingprobs-1; i >= 0; --i )
      GCGpricingprobExitPricing(pricingprobs[i], nroundscol);

   pricingtype_ = NULL;
}

/** setup the priority queue (done once per stabilization round): add all pricing jobs to be performed */
SCIP_RETCODE Pricingcontroller::setupPriorityQueue(
   SCIP_Real*            dualsolconv         /**< dual solution values / Farkas coefficients of convexity constraints */
   )
{
   SCIPdebugMessage("Setup pricing queue, chunk = %d/%d\n", curchunk+1, nchunks);

   GCGpqueueClear(pqueue);

   /* reset pricing problems */
   for( int i = 0; i < npricingprobs; ++i )
      GCGpricingprobReset(gcg, pricingprobs[i]);

   for( int i = 0; i < npricingjobs; ++i )
   {
      int probnr = GCGpricingprobGetProbnr(GCGpricingjobGetPricingprob(pricingjobs[i]));

      SCIP_CALL_EXC( GCGpricingjobSetup(gcg, pricingjobs[i],
         (heurpricingiters > 0
            && (maxheurdepth == -1 || SCIPnodeGetDepth(SCIPgetCurrentNode(scip_)) <= maxheurdepth)),
         sorting, nroundscol, dualsolconv[probnr], GCGpricerGetNPointsProb(gcg, probnr), GCGpricerGetNRaysProb(gcg, probnr)) );

      if( GCGpricingjobGetChunk(pricingjobs[i]) == curchunk )
      {
         SCIP_CALL_EXC( GCGpqueueInsert(pqueue, (void*) pricingjobs[i]) );
      }
   }

   nsolvedprobs = 0;

   return SCIP_OKAY;
}

/** get the next pricing job to be performed */
GCG_PRICINGJOB* Pricingcontroller::getNextPricingjob()
{
   GCG_PRICINGJOB* pricingjob = NULL;

   if( SCIPisStopped(scip_) )
      return NULL;

   do
   {
      pricingjob = (GCG_PRICINGJOB*) GCGpqueueRemove(pqueue);
   }
   while( pricingjob != NULL && pricingprobIsDone(GCGpricingjobGetPricingprob(pricingjob)) );

   return pricingjob;
}

/** add the information that the next branching constraint must be added,
 * and for the pricing job, reset heuristic pricing counter and flag
 */
SCIP_RETCODE Pricingcontroller::pricingprobNextBranchcons(
   GCG_PRICINGPROB*      pricingprob         /**< pricing problem structure */
   )
{
   int i;
   SCIP_RETCODE retcode;

   GCGpricingprobNextBranchcons(pricingprob);

   /* reset heuristic pricing counter and flag for every corresponding pricing job */
   if( heurpricingiters > 0 )
      for( i = 0; i < npricingjobs; ++i )
         if( GCGpricingjobGetPricingprob(pricingjobs[i]) == pricingprob )
            GCGpricingjobResetHeuristic(pricingjobs[i]);

   /* re-sort the priority queue */
   GCG_SET_LOCK(&locks->pricinglock);
   retcode = GCGpqueueResort(pqueue);
   GCG_UNSET_LOCK(&locks->pricinglock);
   
   SCIP_CALL( retcode );

   return SCIP_OKAY;
}

/** set an individual time limit for a pricing job */
SCIP_RETCODE Pricingcontroller::setPricingjobTimelimit(
   GCG_PRICINGJOB*       pricingjob          /**< pricing job */
   )
{
   SCIP* pricingscip = GCGpricingprobGetPricingscip(GCGpricingjobGetPricingprob(pricingjob));
   SCIP_Real mastertimelimit;
   SCIP_Real timelimit;

   SCIP_CALL( SCIPgetRealParam(scip_, "limits/time", &mastertimelimit) );

   /* do not give pricing job more time than is left for solving the master problem */
   timelimit = MAX(0, mastertimelimit - SCIPgetSolvingTime(scip_));

   SCIP_CALL( SCIPsetRealParam(pricingscip, "limits/time", timelimit) );

   return SCIP_OKAY;
}

/** update solution information of a pricing problem */
void Pricingcontroller::updatePricingprob(
   GCG_PRICINGPROB*      pricingprob,        /**< pricing problem structure */
   GCG_PRICINGSTATUS     status,             /**< new pricing status */
   SCIP_Real             lowerbound,         /**< new lower bound */
   int                   nimpcols            /**< number of new improving columns */
   )
{
   GCGpricingprobUpdate(gcg, pricingprob, status, lowerbound, nimpcols);
}

/** decide whether a pricing job must be treated again */
void Pricingcontroller::evaluatePricingjob(
   GCG_PRICINGJOB*       pricingjob,        /**< pricing job */
   GCG_PRICINGSTATUS     status             /**< status of pricing job */
   )
{
   GCG_PRICINGPROB* pricingprob = GCGpricingjobGetPricingprob(pricingjob);
   SCIP_Bool heuristic = GCGpricingjobIsHeuristic(pricingjob);

   /* Go to the next heuristic pricing iteration */
   if( heuristic )
      GCGpricingjobIncreaseNHeurIters(pricingjob);

   /* If the pricing job has not yielded any improving column, possibly solve it again;
    * increase at least one of its limits, or solve it exactly if it was solved heuristically before
    */
   // @todo: update score of pricing job
   if( !pricingprobIsDone(pricingprob) )
   {
      SCIPdebugMessage("Solving problem %d with <%s> has not yielded improving columns.\n",
         GCGpricingprobGetProbnr(pricingprob), GCGsolverGetName(GCGpricingjobGetSolver(pricingjob)));

      if( heuristic && status != GCG_PRICINGSTATUS_OPTIMAL && status != GCG_PRICINGSTATUS_NOTAPPLICABLE )
      {
         assert(status == GCG_PRICINGSTATUS_UNKNOWN || status == GCG_PRICINGSTATUS_SOLVERLIMIT);

         if( status != GCG_PRICINGSTATUS_SOLVERLIMIT || GCGpricingjobGetNHeurIters(pricingjob) >= heurpricingiters )
         {
            GCGpricingjobSetExact(pricingjob);
            SCIPdebugMessage("  -> set exact\n");
         }
         else
         {
            SCIPdebugMessage("  -> increase a limit\n");
         }
         
         GCG_SET_LOCK(&locks->pricinglock);
         SCIP_CALL_EXC( GCGpqueueInsert(pqueue, (void*) pricingjob) );
         GCG_UNSET_LOCK(&locks->pricinglock);

         return;
      }

      if( pricingprobNeedsNextBranchingcons(pricingprob) )
      {

         SCIPdebugMessage("  -> consider next generic branching constraint.\n");

         SCIP_CALL_EXC( pricingprobNextBranchcons(pricingprob) );

         GCG_SET_LOCK(&locks->pricinglock);
         SCIP_CALL_EXC( GCGpqueueInsert(pqueue, (void*) pricingjob) );
         GCG_UNSET_LOCK(&locks->pricinglock);

         return;
      }

      GCGpricingjobNextSolver(gcg, pricingjob);
      if( heurpricingiters > 0 )
         GCGpricingjobResetHeuristic(pricingjob);
      if( GCGpricingjobGetSolver(pricingjob) != NULL )
      {
         SCIPdebugMessage("  -> use another solver\n");
         
         GCG_SET_LOCK(&locks->pricinglock);
         SCIP_CALL_EXC( GCGpqueueInsert(pqueue, (void*) pricingjob) );
         GCG_UNSET_LOCK(&locks->pricinglock);
      }
   }
   else
   {
      #pragma omp atomic update
      ++nsolvedprobs;
   }
}

/** collect solution results from all pricing problems */
void Pricingcontroller::collectResults(
   GCG_COL**             bestcols,           /**< best found columns per pricing problem */
   SCIP_Bool*            infeasible,         /**< pointer to store whether pricing is infeasible */
   SCIP_Bool*            optimal,            /**< pointer to store whether all pricing problems were solved to optimality */
   SCIP_Real*            bestobjvals,        /**< array to store best lower bounds */
   SCIP_Real*            beststabobj,        /**< pointer to store total lower bound */
   SCIP_Real*            bestredcost,        /**< pointer to store best total reduced cost */
   SCIP_Bool*            bestredcostvalid    /**< pointer to store whether best reduced cost is valid */
   )
{
   int nblocks = GCGgetNPricingprobs(gcg);
   SCIP_Bool foundcols = FALSE;

   /* initializations */
   *infeasible = FALSE;
   *optimal = TRUE;
   *beststabobj = 0.0;
   *bestredcost = 0.0;
   for( int i = 0; i < nblocks; ++i )
      bestobjvals[i] = -SCIPinfinity(scip_);

   for( int i = 0; i < npricingprobs; ++i )
   {
      int probnr = GCGpricingprobGetProbnr(pricingprobs[i]);
      int nidentblocks = GCGgetNIdenticalBlocks(gcg, probnr);
      SCIP_Real lowerbound = GCGpricingprobGetLowerbound(pricingprobs[i]);

      /* check infeasibility */
      *infeasible |= GCGpricingprobGetStatus(pricingprobs[i]) == GCG_PRICINGSTATUS_INFEASIBLE;

      /* check optimality */
      *optimal &= GCGpricingprobGetStatus(pricingprobs[i]) == GCG_PRICINGSTATUS_OPTIMAL;

      if( GCGpricingprobGetNImpCols(pricingprobs[i]) > 0 )
         foundcols = TRUE;

      /* update lower bound information */
      bestobjvals[probnr] = SCIPisInfinity(scip_, ABS(lowerbound)) ? lowerbound : nidentblocks * lowerbound;
      if( SCIPisInfinity(scip_, -lowerbound) )
         *beststabobj = -SCIPinfinity(scip_);
      else if( !SCIPisInfinity(scip_, -(*beststabobj)) )
         *beststabobj += bestobjvals[probnr];

      if( bestcols[probnr] != NULL )
         *bestredcost += GCGcolGetRedcost(bestcols[probnr]) * nidentblocks;
   }

   *infeasible |= (pricingtype_->getType() == GCG_PRICETYPE_FARKAS && *optimal && !foundcols);
   *bestredcostvalid &= foundcols || *optimal;
}

/** check if the next chunk of pricing problems is to be used */
SCIP_Bool Pricingcontroller::checkNextChunk()
{
   int nextchunk = (curchunk + 1) % nchunks;

   if( nextchunk == startchunk )
   {
      SCIPdebugMessage("not considering next chunk.\n");
      return FALSE;
   }
   else
   {
      SCIPdebugMessage("need considering next chunk = %d/%d\n", nextchunk+1, nchunks);
      curchunk = nextchunk;
      return TRUE;
   }
}

/** decide whether the pricing loop can be aborted */
SCIP_Bool Pricingcontroller::canPricingloopBeAborted(
   PricingType*          pricingtype,        /**< type of pricing (reduced cost or Farkas) */
   int                   nfoundcols,         /**< number of negative reduced cost columns found so far */
   int                   nsuccessfulprobs    /**< number of pricing problems solved successfully so far */
   ) const
{
   int nrelpricingprobs = GCGgetNRelPricingprobs(gcg);
   int _nsolvedprobs;

   if( eagerage == eagerfreq )
      return FALSE;

   #pragma omp atomic read
   _nsolvedprobs = nsolvedprobs;

   return !((nfoundcols < pricingtype->getMaxcolsround())
         && nsuccessfulprobs < pricingtype->getMaxsuccessfulprobs()
         && nsuccessfulprobs < pricingtype->getRelmaxsuccessfulprobs() * nrelpricingprobs
         && (nfoundcols == 0 || _nsolvedprobs < pricingtype->getRelmaxprobs() * nrelpricingprobs));
}

void Pricingcontroller::resetEagerage()
{
   eagerage = 0;
}

void Pricingcontroller::increaseEagerage()
{
   if( eagerfreq > 0 )
   {
      eagerage++;
   }
}

/** for a given problem index, get the corresponding pricing problem (or NULL, if it does not exist) */
GCG_PRICINGPROB* Pricingcontroller::getPricingprob(
   int                   probnr              /**< index of the pricing problem */
   )
{
   for( int i = 0; i < npricingprobs; ++i )
      if( GCGpricingprobGetProbnr(pricingprobs[i]) == probnr )
         return pricingprobs[i];

   return NULL;
}

/** get maximal possible number of pricing iterations */
int Pricingcontroller::getMaxNIters() const
{
   return maxniters;
}

} /* namespace gcg */
