/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2012 Operations Research, RWTH Aachen University       */
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

/**@file   objpricer_gcg.h
 * @ingroup PUBLICMETHODS
 * @brief  GCG variable pricer
 * @author Martin Bergner
 * @ingroup PRICERS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "objscip/objscip.h"
#include "class_pricingtype.h"

#ifndef __SCIP_OBJPRICER_GCG__
#define __SCIP_OBJPRICER_GCG__


class ObjPricerGcg : public scip::ObjPricer
{
public:
   /*lint --e{1540}*/

   SCIP*                 origprob;           /**< the original program */
   SCIP_PRICERDATA *pricerdata;
   static int threads;

   /** default constructor */
   ObjPricerGcg(
         SCIP* scip, /**< SCIP data structure */
         SCIP*              origscip,           /**< SCIP data structure of original problem */
         const char* name, /**< name of variable pricer */
         const char* desc, /**< description of variable pricer */
         int priority, /**< priority of the variable pricer */
         unsigned int delay, /**< should the pricer be delayed until no other pricers or already existing*/
         SCIP_PRICERDATA *pricerdata
   );
   /** destructor */
   virtual ~ObjPricerGcg()
   {};

   /** destructor of variable pricer to free user data (called when SCIP is exiting) */
   virtual SCIP_DECL_PRICERFREE(scip_free);
   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_DECL_PRICERINIT(scip_init);
   /** deinitialization method of variable pricer (called before transformed problem is freed) */
   virtual SCIP_DECL_PRICEREXIT(scip_exit);
   /** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
   virtual SCIP_DECL_PRICERINITSOL(scip_initsol);
   /** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
   virtual SCIP_DECL_PRICEREXITSOL(scip_exitsol);
   /** reduced cost pricing method of variable pricer for feasible LPs */
   virtual SCIP_DECL_PRICERREDCOST(scip_redcost);
   /** farkas pricing method of variable pricer for infeasible LPs */
   virtual SCIP_DECL_PRICERFARKAS(scip_farkas);
   inline SCIP_PRICERDATA *getPricerdata()
   {
      return pricerdata;
   }

   ;
   /** computes the pricing problem objectives */
   SCIP_RETCODE setPricingObjs(
         PricingType *pricetype);

   /** performs the pricing routine, gets the type of pricing that should be done: farkas or redcost pricing */
   SCIP_RETCODE priceNewVariables(
         PricingType*   pricetype,          /**< type of the pricing */
         SCIP_RESULT*   result,             /**< result pointer */
         double*        lowerbound          /**< lower bound of pricingproblems */
      );

   /** creates a new master variable corresponding to the given solution and problem */
   SCIP_RETCODE createNewMasterVar(
         SCIP*          scip,               /**< SCIP data structure */
         SCIP_VAR**     solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
         double*        solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
         int            nsolvars,           /**< number of variables in array solvars */
         unsigned int   solisray,           /**< is the solution a ray? */
         int            prob,               /**< number of the pricing problem the solution belongs to */
         unsigned int   force,              /**< should the given variable be added also if it has non-negative reduced cost? */
         unsigned int*  added,              /**< pointer to store whether the variable was successfully added */
         SCIP_VAR**     addedvar
      );

   /** performs optimal or farkas pricing */
   SCIP_RETCODE performPricing(
         PricingType *pricetype, /**< type of pricing */
         unsigned int optimal, /**< heuristic or optimal pricing */
         SCIP_RESULT *result, /**< result pointer */
         int * nfoundvars, /**< pointer to store number of found variables */
         double* bestredcost, /**< pointer to store reduced cost */
         unsigned int* bestredcostvalid);

    FarkasPricing *getFarkasPricing() const
    {
        return farkaspricing;
    }

    ReducedCostPricing *getReducedCostPricing() const
    {
        return reducedcostpricing;
    }

    /** ensures size of solvers array */
    SCIP_RETCODE ensureSizeSolvers();

   SCIP* getOrigprob() const
   {
      return origprob;
   }

private:
   ReducedCostPricing *reducedcostpricing;
   FarkasPricing *farkaspricing;
   //PricingMode *pricingmode;

   /** free pricing problems */
   SCIP_RETCODE freePricingProblems();

   /** returns whether pricing can be aborted */
   SCIP_Bool abortPricing(
      PricingType*          pricetype,          /**< type of pricing*/
      int                   nfoundvars,         /**< number of variables found so far */
      int                   solvedmips,         /**< number of MIPS solved so far */
      int                   successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Bool             optimal             /**< optimal or heuristic pricing */
   );

   SCIP_Real  computeRedCost(
      SCIP_VAR**            solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
      SCIP_Real*            solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
      int                   nsolvars,           /**< number of variables in array solvars */
      SCIP_Bool             solisray,           /**< is the solution a ray? */
      int                   prob               /**< number of the pricing problem the solution belongs to */
      );

   int countPricedVariables(
      int& prob,
      SCIP_SOL** sols,
      int nsols,
      SCIP_Bool* solisray
      );

   /** return TRUE or FALSE whether the master LP is solved to optimality */
   SCIP_Bool isMasterLPOptimal();


   /** return TRUE or FALSE whether pricing problem has been solved to optimality */
   SCIP_Bool  isPricingOptimal(
      SCIP*                 scip                /**< SCIP data structure */
      );

   /** ensures size of pricedvars array */
   SCIP_RETCODE ensureSizePricedvars(
      int                   size                /**< needed size */
      );

   /** adds new variable to the end of the priced variables array */
   SCIP_RETCODE addVariableToPricedvars(
      SCIP_VAR*             newvar              /**< variable to add */
      );

   /** add master variable to all constraints */
   SCIP_RETCODE addVariableToMasterconstraints(
      SCIP_VAR*             newvar,             /**< The new variable to add */
      int                   prob,               /**< number of the pricing problem the solution belongs to */
      SCIP_VAR**            solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
      SCIP_Real*            solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
      int                   nsolvars            /**< number of variables in array solvars */
      );

   /**
    * check whether pricing can be aborted:
    * if objective value is always integral and the current node's current
    * lowerbound rounded up equals the current lp objective value rounded
    * up we don't need to continue pricing since the best possible feasible
    * solution must have at least this value
    */
   SCIP_Bool canPricingBeAborted();

   /** sorts pricing problems according to their score */
   void sortPricingProblemsByScore();

   /** returns the gegeneracy of the masterproblem */
   SCIP_RETCODE computeCurrentDegeneracy(
      double*               degeneracy          /**< pointer to store degeneracy */
      );

   /** initializes the pointers to the appropriate structures */
   SCIP_RETCODE getSolverPointers(
      GCG_SOLVER*           solver,             /**< pricing solver */
      PricingType*          pricetype,          /**< type of pricing: optimal or heuristic */
      SCIP_Bool             optimal,            /**< should the pricing problem be solved optimal or heuristically */
      SCIP_CLOCK**          clock,              /**< clock belonging to this setting */
      int**                 calls,              /**< calls belonging to this setting */
      GCG_DECL_SOLVERSOLVE((**solversolve))     /**< solving function belonging to this setting */
      );

   /** set subproblem timelimit */
   SCIP_RETCODE setPricingProblemTimelimit(
      SCIP*                 pricingscip         /**< SCIP of the pricingproblem */
      );

   /** set subproblem memory limit */
   SCIP_RETCODE setPricingProblemMemorylimit(
      SCIP*                 pricingscip         /**< SCIP of the pricingproblem */
      );

   /** set all pricing problem limits */
   SCIP_RETCODE setPricingProblemLimits(
      int                   prob,               /**< index of the pricing problem */
      SCIP_Bool             optimal            /**< heuristic or optimal pricing */
      );

   /** solves a specific pricing problem
    * @todo simplify
    */
   SCIP_RETCODE solvePricingProblem(
      int                   prob,               /**< index of pricing problem */
      PricingType*          pricetype,          /**< type of pricing: optimal or heuristic */
      SCIP_Bool             optimal,            /**< should the pricing problem be solved optimal or heuristically */
      SCIP_Real*            lowerbound,         /**< dual bound returned by pricing problem */
      SCIP_SOL**            sols,               /**< pointer to store solutions */
      SCIP_Bool*            solisray,           /**< array to indicate whether solution is a ray */
      int                   maxsols,            /**< size of the sols array to indicate maximum solutions */
      int*                  nsols,              /**< number of solutions */
      SCIP_STATUS*          status              /**< solution status of the pricing problem */
      );

   /** frees all solvers */
   SCIP_RETCODE solversFree();

   /** calls the init method on all solvers */
   SCIP_RETCODE solversInit();

   /** calls the exit method on all solvers */
   SCIP_RETCODE solversExit();

   /** calls the initsol method on all solvers */
   SCIP_RETCODE solversInitsol();

   /** calls the exitsol method of all solvers */
   SCIP_RETCODE solversExitsol();
};

int ObjPricerGcg::threads;
#endif
