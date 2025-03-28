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

/**@file   solver.c
 * @brief  methods for GCG pricing solvers
 * @author Henri Lotze
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/pub_solver.h"
#include "gcg/solver.h"
#include "gcg/struct_solver.h"

#include "gcg/gcg.h"
#include "gcg/pricer_gcg.h"

#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif


/** compares two solvers w. r. t. their priorities */
SCIP_DECL_SORTPTRCOMP(GCGsolverComp)
{  /*lint --e{715}*/
   GCG_SOLVER* solver1 = (GCG_SOLVER*) elem1;
   GCG_SOLVER* solver2 = (GCG_SOLVER*) elem2;

   assert(solver1 != NULL);
   assert(solver2 != NULL);

   return solver2->priority - solver1->priority; /* prefer higher priorities */
}

/** creates a GCG pricing solver */
SCIP_RETCODE GCGsolverCreate(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER**          solver,             /**< pointer to pricing solver data structure */
   const char*           name,               /**< name of solver */
   const char*           desc,               /**< description of solver */
   int                   priority,           /**< priority of solver */
   SCIP_Bool             heurenabled,        /**< flag to indicate whether heuristic solving method of the solver is enabled */
   SCIP_Bool             exactenabled,        /**< flag to indicate whether exact solving method of the solver is enabled */
   GCG_DECL_SOLVERUPDATE((*solverupdate)),   /**< update method for solver */
   GCG_DECL_SOLVERSOLVE  ((*solversolve)),   /**< solving method for solver */
   GCG_DECL_SOLVERSOLVEHEUR((*solveheur)),   /**< heuristic solving method for solver */
   GCG_DECL_SOLVERFREE   ((*solverfree)),    /**< free method of solver */
   GCG_DECL_SOLVERINIT   ((*solverinit)),    /**< init method of solver */
   GCG_DECL_SOLVEREXIT   ((*solverexit)),    /**< exit method of solver */
   GCG_DECL_SOLVERINITSOL((*solverinitsol)), /**< initsol method of solver */
   GCG_DECL_SOLVEREXITSOL((*solverexitsol)), /**< exitsol method of solver */
   GCG_SOLVERDATA*       solverdata          /**< pricing solver data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];
   SCIP* masterprob;

   assert(gcg != NULL);
   assert(solver != NULL);

   masterprob = GCGgetMasterprob(gcg);

   if( solveheur == NULL && solversolve == NULL )
   {
      SCIPwarningMessage(masterprob, "Solver <%s> has neither heuristic nor exact solving method and will not be included.\n", name);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocMemory(masterprob, solver) ); /*lint !e866*/

   SCIP_ALLOC( BMSduplicateMemoryArray(&(*solver)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*solver)->desc, desc, strlen(desc)+1) );

   (*solver)->solverupdate = solverupdate;
   (*solver)->solversolve = solversolve;
   (*solver)->solversolveheur = solveheur;
   (*solver)->solverfree = solverfree;
   (*solver)->solverinit = solverinit;
   (*solver)->solverexit = solverexit;
   (*solver)->solverinitsol = solverinitsol;
   (*solver)->solverexitsol = solverexitsol;
   (*solver)->solverdata = solverdata;

   SCIP_CALL( SCIPcreateCPUClock(masterprob, &((*solver)->optfarkasclock)) );
   SCIP_CALL( SCIPcreateCPUClock(masterprob, &((*solver)->optredcostclock)) );
   SCIP_CALL( SCIPcreateCPUClock(masterprob, &((*solver)->heurfarkasclock)) );
   SCIP_CALL( SCIPcreateCPUClock(masterprob, &((*solver)->heurredcostclock)) );

   (*solver)->optfarkascalls = 0;
   (*solver)->optredcostcalls = 0;
   (*solver)->heurfarkascalls = 0;
   (*solver)->heurredcostcalls = 0;

   if( solveheur != NULL )
   {
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "pricingsolver/%s/heurenabled", name);
      (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "flag to indicate whether heuristic solving method of solver <%s> is enabled", name);
      SCIP_CALL( SCIPaddBoolParam(GCGgetOrigprob(gcg), paramname, paramdesc,
                     &((*solver)->heurenabled), FALSE, heurenabled, NULL, NULL));
   }
   else
      (*solver)->heurenabled = FALSE;

   if( solversolve != NULL )
   {
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "pricingsolver/%s/exactenabled", name);
      (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "flag to indicate whether exact solving method of solver <%s> is enabled", name);
      SCIP_CALL( SCIPaddBoolParam(GCGgetOrigprob(gcg), paramname, paramdesc,
                     &((*solver)->exactenabled), FALSE, exactenabled, NULL, NULL));
   }
   else
      (*solver)->exactenabled = FALSE;

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "pricingsolver/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of solver <%s>", name);
   SCIP_CALL( SCIPaddIntParam(GCGgetOrigprob(gcg), paramname, paramdesc,
                  &((*solver)->priority), FALSE, priority, INT_MIN/4, INT_MAX/4, NULL, NULL));

   return SCIP_OKAY;
}

/** calls destructor and frees memory of GCG pricing solver */
SCIP_RETCODE GCGsolverFree(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER**          solver              /**< pointer to pricing solver data structure */
   )
{
   SCIP* masterprob;
   assert(gcg != NULL);
   assert(solver != NULL);
   assert(*solver != NULL);

   masterprob = GCGgetMasterprob(gcg);

   if( (*solver)->solverfree != NULL )
   {
      SCIP_CALL( (*solver)->solverfree(gcg, *solver) );
   }

   BMSfreeMemoryArray(&(*solver)->name);
   BMSfreeMemoryArray(&(*solver)->desc);

   SCIP_CALL( SCIPfreeClock(masterprob, &((*solver)->optfarkasclock)) );
   SCIP_CALL( SCIPfreeClock(masterprob, &((*solver)->optredcostclock)) );
   SCIP_CALL( SCIPfreeClock(masterprob, &((*solver)->heurfarkasclock)) );
   SCIP_CALL( SCIPfreeClock(masterprob, &((*solver)->heurredcostclock)) );

   SCIPfreeMemory(masterprob, solver);

   return SCIP_OKAY;
}

/** initializes GCG pricing solver */
SCIP_RETCODE GCGsolverInit(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   SCIP_Bool resetstat;
   SCIP* masterprob = GCGgetMasterprob(gcg);

   assert(gcg != NULL);
   assert(solver != NULL);

   SCIP_CALL( SCIPgetBoolParam(masterprob, "misc/resetstat", &resetstat) );

   if( resetstat )
   {
      SCIP_CALL( SCIPresetClock(masterprob, solver->optfarkasclock) );
      SCIP_CALL( SCIPresetClock(masterprob, solver->optredcostclock) );
      SCIP_CALL( SCIPresetClock(masterprob, solver->heurfarkasclock) );
      SCIP_CALL( SCIPresetClock(masterprob, solver->heurredcostclock) );

      solver->optfarkascalls = 0;
      solver->optredcostcalls = 0;
      solver->heurfarkascalls = 0;
      solver->heurredcostcalls = 0;
   }

   if( solver->solverinit != NULL )
   {
      SCIP_CALL( solver->solverinit(gcg, solver) );
   }

   return SCIP_OKAY;
}

/** calls exit method of GCG pricing solver */
SCIP_RETCODE GCGsolverExit(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(gcg != NULL);
   assert(solver != NULL);

   if( solver->solverexit != NULL )
   {
      SCIP_CALL( solver->solverexit(gcg, solver) );
   }

   return SCIP_OKAY;
}

/** calls solving process initialization method of GCG pricing solver */
SCIP_RETCODE GCGsolverInitsol(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(gcg != NULL);
   assert(solver != NULL);

   if( solver->solverinitsol != NULL )
   {
      SCIP_CALL( solver->solverinitsol(gcg, solver) );
   }

   return SCIP_OKAY;
}

/** calls solving process deinitialization method of GCG pricing solver */
SCIP_RETCODE GCGsolverExitsol(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(gcg != NULL);
   assert(solver != NULL);

   if( solver->solverexitsol != NULL )
   {
      SCIP_CALL( solver->solverexitsol(gcg, solver) );
   }

   return SCIP_OKAY;
}

/** calls update method of GCG pricing solver */
SCIP_RETCODE GCGsolverUpdate(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP*                 pricingprob,        /**< the pricing problem that should be solved */
   GCG_SOLVER*           solver,             /**< pricing solver */
   int                   probnr,             /**< number of the pricing problem */
   SCIP_Bool             varobjschanged,     /**< have the objective coefficients changed? */
   SCIP_Bool             varbndschanged,     /**< have the lower and upper bounds changed? */
   SCIP_Bool             consschanged        /**< have the constraints changed? */
   )
{
   assert(pricingprob != NULL);
   assert(solver != NULL);

   if( solver->solverupdate != NULL && (varobjschanged || varbndschanged || consschanged) )
   {
      SCIP_CALL( solver->solverupdate(gcg, pricingprob, solver, probnr, varobjschanged, varbndschanged, consschanged) );
   }

   return SCIP_OKAY;
}

/** calls heuristic or exact solving method of GCG pricing solver
 * @note This method has to be threadsafe!
 */
SCIP_RETCODE GCGsolverSolve(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP*                 pricingprob,        /**< the pricing problem that should be solved */
   GCG_SOLVER*           solver,             /**< pricing solver */
   SCIP_Bool             redcost,            /**< is reduced cost (TRUE) or Farkas (FALSE) pricing performed? */
   SCIP_Bool             heuristic,          /**< shall the pricing problem be solved heuristically? */
   int                   probnr,             /**< number of the pricing problem */
   SCIP_Real             dualsolconv,        /**< dual solution of the corresponding convexity constraint */
   SCIP_Real*            lowerbound,         /**< pointer to store lower bound of pricing problem */
   GCG_PRICINGSTATUS*    status,             /**< pointer to store the returned pricing status */
   SCIP_Bool*            solved              /**< pointer to store whether the solution method was called */
   )
{
   SCIP_CLOCK* clock;
   SCIP* masterprob;

   assert(gcg != NULL);
   assert(pricingprob != NULL);
   assert(solver != NULL);
   assert(lowerbound != NULL);
   assert(status != NULL);

   masterprob = GCGgetMasterprob(gcg);

   *solved = FALSE;

   if( heuristic )
   {
      if( solver->heurenabled )
      {
         assert(solver->solversolveheur != NULL);

         if( redcost )
            clock = solver->heurredcostclock;
         else
            clock = solver->heurfarkasclock;

#ifdef _OPENMP
         if( omp_get_num_threads() == 1 )
            SCIP_CALL_ABORT( SCIPstartClock(masterprob, clock) );
#else
         SCIP_CALL_ABORT( SCIPstartClock(masterprob, clock) );
#endif

         SCIP_CALL( solver->solversolveheur(gcg, pricingprob, solver, probnr, dualsolconv, lowerbound, status) );
         *solved = TRUE;

#ifdef _OPENMP
         if( omp_get_num_threads() == 1 )
            SCIP_CALL_ABORT( SCIPstopClock(masterprob, clock) );
#else
         SCIP_CALL_ABORT( SCIPstopClock(masterprob, clock) );
#endif
      }
   }
   else
   {
      if( solver->exactenabled )
      {
         assert(solver->solversolve != NULL);

         if( redcost )
            clock = solver->optredcostclock;
         else
            clock = solver->optfarkasclock;

#ifdef _OPENMP
         if( omp_get_num_threads() == 1 )
            SCIP_CALL_ABORT( SCIPstartClock(masterprob, clock) );
#else
         SCIP_CALL_ABORT( SCIPstartClock(masterprob, clock) );
#endif

         SCIP_CALL( solver->solversolve(gcg, pricingprob, solver, probnr, dualsolconv, lowerbound, status) );
         *solved = TRUE;

#ifdef _OPENMP
         if( omp_get_num_threads() == 1 )
            SCIP_CALL_ABORT( SCIPstopClock(masterprob, clock) );
#else
          SCIP_CALL_ABORT( SCIPstopClock(masterprob, clock) );
#endif
      }
   }

   if( *status != GCG_PRICINGSTATUS_NOTAPPLICABLE && *solved )
   {
      if( redcost )
         if( heuristic )
            #pragma omp atomic update
            ++solver->heurredcostcalls;
         else
            #pragma omp atomic update
            ++solver->optredcostcalls;
      else
         if( heuristic )
            #pragma omp atomic update
            ++solver->heurfarkascalls;
         else
            #pragma omp atomic update
            ++solver->optfarkascalls;
   }

   return SCIP_OKAY;
}

/** gets user data of GCG pricing solver */
GCG_SOLVERDATA* GCGsolverGetData(
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(solver != NULL);

   return solver->solverdata;
}

/** sets user data of GCG pricing solver */
void GCGsolverSetData(
   GCG_SOLVER*           solver,             /**< pricing solver */
   GCG_SOLVERDATA*       solverdata          /**< pricing solver data */
   )
{
   assert(solver != NULL);

   solver->solverdata = solverdata;
}

/** gets name of GCG pricing solver */
const char* GCGsolverGetName(
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(solver != NULL);

   return solver->name;
}

/** gets description of GCG pricing solver */
const char* GCGsolverGetDesc(
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(solver != NULL);
   
   return solver->desc;
}

/** gets priority of GCG pricing solver */
int GCGsolverGetPriority(
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(solver != NULL);
   
   return solver->priority;
}

/** gets whether heuristic solving method of GCG pricing solver is enabled */
SCIP_Bool GCGsolverIsHeurEnabled(
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(solver != NULL);

   return solver->heurenabled;
}

/** gets whether exact solving method of GCG pricing solver is enabled */
SCIP_Bool GCGsolverIsExactEnabled(
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(solver != NULL);
   
   return solver->exactenabled;
}

/** gets number of exact Farkas pricing calls of pricing solver */
int GCGsolverGetOptFarkasCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(solver != NULL);

   return solver->optfarkascalls;
}

/** gets number of exact reduced cost pricing calls of pricing solver */
int GCGsolverGetOptRedcostCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(solver != NULL);

   return solver->optredcostcalls;
}

/** gets number of heuristic Farkas pricing calls of pricing solver */
int GCGsolverGetHeurFarkasCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(solver != NULL);

   return solver->heurfarkascalls;
}

/** gets number of heuristic reduced cost pricing calls of pricing solver */
int GCGsolverGetHeurRedcostCalls(
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(solver != NULL);

   return solver->heurredcostcalls;
}

/** gets exact Farkas pricing time of pricing solver */
SCIP_Real GCGsolverGetOptFarkasTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(gcg != NULL);
   assert(solver != NULL);

   return SCIPgetClockTime(GCGgetMasterprob(gcg), solver->optfarkasclock);
}

/** gets exact reduced cost pricing time of pricing solver */
SCIP_Real GCGsolverGetOptRedcostTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(gcg != NULL);
   assert(solver != NULL);

   return SCIPgetClockTime(GCGgetMasterprob(gcg), solver->optredcostclock);
}

/** gets heuristic Farkas pricing time of pricing solver */
SCIP_Real GCGsolverGetHeurFarkasTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(gcg != NULL);
   assert(solver != NULL);

   return SCIPgetClockTime(GCGgetMasterprob(gcg), solver->heurfarkasclock);
}

/** gets heuristic reduced cost pricing time of pricing solver */
SCIP_Real GCGsolverGetHeurRedcostTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_SOLVER*           solver              /**< pricing solver */
   )
{
   assert(gcg != NULL);
   assert(solver != NULL);

   return SCIPgetClockTime(GCGgetMasterprob(gcg), solver->heurredcostclock);
}
