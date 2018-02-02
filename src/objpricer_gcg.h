/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

#ifndef GCG_OBJPRICER_GCG_H_
#define GCG_OBJPRICER_GCG_H_

#include "objscip/objscip.h"
#include "class_pricingtype.h"
#include "class_stabilization.h"
#include "class_colpool.h"
#include "pub_gcgcol.h"

using gcg::Stabilization;
using gcg::Colpool;

class ObjPricerGcg : public scip::ObjPricer
{
public:
   /*lint --e{1540}*/

   SCIP*              origprob;           /**< the original program */
   SCIP_PRICERDATA *pricerdata;           /**< pricerdata data structure */
   Colpool *colpool;                      /**< column pool */
   static int threads;

   /** default constructor */
   ObjPricerGcg(
         SCIP* scip, /**< SCIP data structure */
         SCIP*              origscip,           /**< SCIP data structure of original problem */
         const char* name, /**< name of variable pricer */
         const char* desc, /**< description of variable pricer */
         int priority, /**< priority of the variable pricer */
         unsigned int delay, /**< should the pricer be delayed until no other pricers or already existing*/
         SCIP_PRICERDATA *pricerdata /**< pricerdata data structure */
   );
   /** destructor */
   virtual ~ObjPricerGcg()
   {}

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

   /** computes the pricing problem objectives */
   SCIP_RETCODE setPricingObjs(
      PricingType*          pricetype,          /**< Farkas or Reduced cost pricing */
      SCIP_Bool             stabilize           /**< do we use stabilization ? */
   );

   /** update reduced cost of columns in column pool */
   SCIP_RETCODE updateRedcostColumnPool(
      PricingType*          pricetype           /**< type of pricing: reduced cost or Farkas */
      );
   /** method to price new columns from Column Pool */
   SCIP_RETCODE priceColumnPool(
      PricingType*          pricetype,          /**< type of pricing: reduced cost or Farkas */
      int*                  pnfoundvars         /**< pointer to store number of priced variables */
      );

   /** performs the pricing routine, gets the type of pricing that should be done: farkas or redcost pricing */
   SCIP_RETCODE priceNewVariables(
      PricingType*       pricetype,          /**< type of the pricing */
      SCIP_RESULT*       result,             /**< result pointer */
      double*            lowerbound          /**< lower bound of pricingproblems */
   );

   /** creates a new master variable corresponding to the given solution and problem */
   SCIP_RETCODE createNewMasterVar(
      SCIP*              scip,               /**< SCIP data structure */
      PricingType*       pricetype,          /**< type of the pricing */
      SCIP_SOL*          sol,                /**< solution to compute reduced cost for */
      SCIP_VAR**         solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
      double*            solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
      int                nsolvars,           /**< number of variables in array solvars */
      unsigned int       solisray,           /**< is the solution a ray? */
      int                prob,               /**< number of the pricing problem the solution belongs to */
      unsigned int       force,              /**< should the given variable be added also if it has non-negative reduced cost? */
      unsigned int*      added,              /**< pointer to store whether the variable was successfully added */
      SCIP_VAR**         addedvar            /**< pointer to store the created variable */
   );

   /** creates a new master variable corresponding to the given gcg column */
   SCIP_RETCODE createNewMasterVarFromGcgCol(
      SCIP*                 scip,               /**< SCIP data structure */
      PricingType*          pricetype,          /**< type of pricing */
      GCG_COL*              gcgcol,             /**< GCG column data structure */
      SCIP_Bool             force,              /**< should the given variable be added also if it has non-negative reduced cost? */
      SCIP_Bool*            added,              /**< pointer to store whether the variable was successfully added */
      SCIP_VAR**            addedvar            /**< pointer to store the created variable */
   );

   /** performs optimal or farkas pricing */
   SCIP_RETCODE performPricing(
      PricingType*   pricetype,          /**< type of pricing */
      unsigned int   optimal,            /**< heuristic or optimal pricing */
      SCIP_RESULT*   result,             /**< result pointer */
      int*           nfoundvars,         /**< pointer to store number of found variables */
      SCIP_Real*     lowerbound,         /**< pointer to store lowerbound obtained due to lagrange bound */
      SCIP_Bool*     bestredcostvalid    /**< pointer to store if bestredcost are valid (pp solvedoptimal) */
   );

   const FarkasPricing *getFarkasPricing() const
   {
      return farkaspricing;
   }

   const ReducedCostPricing *getReducedCostPricing() const
   {
      return reducedcostpricing;
   }

   ReducedCostPricing *getReducedCostPricingNonConst()
   {
      return reducedcostpricing;
   }

   /** ensures size of solvers array */
   SCIP_RETCODE ensureSizeSolvers();

   SCIP* getOrigprob()
   {
      return origprob;
   }

   /** create the pointers for the pricing types */
   SCIP_RETCODE createPricingTypes();

   /** create the pointers for the stabilization */
   void createStabilization();

   /** create the pointers for the colpool */
   void createColpool();

   /* computes the objective value of the current (stabilized) dual variables) in the dual program */
   SCIP_RETCODE getStabilizedDualObjectiveValue(
      SCIP_Real*         stabdualval         /**< pointer to store stabilized dual objective value */
   );
private:
   ReducedCostPricing *reducedcostpricing;
   FarkasPricing *farkaspricing;
   Stabilization *stabilization;
   /** free pricing problems */
   SCIP_RETCODE freePricingProblems();

   /** returns whether pricing can be aborted */
   SCIP_Bool abortPricing(
      PricingType*          pricetype,          /**< type of pricing*/
      int                   nfoundvars,         /**< number of variables found so far */
      int                   solvedmips,         /**< number of MIPS solved so far */
      int                   successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Bool             optimal             /**< optimal or heuristic pricing */
   ) const;

   SCIP_Real  computeRedCost(
      PricingType*          pricetype,          /**< type of pricing */
      SCIP_SOL*             sol,                /**< solution to compute reduced cost for */
      SCIP_Bool             solisray,           /**< is the solution a ray? */
      int                   prob,               /**< number of the pricing problem the solution belongs to */
      SCIP_Real*            objvalptr           /**< pointer to store the computed objective value */
   ) const;

   SCIP_Real computeRedCostGcgCol(
      PricingType*          pricetype,          /**< type of pricing */
      GCG_Col*              gcgcol,             /**< gcg column to compute reduced cost for */
      SCIP_Real*            objvalptr           /**< pointer to store the computed objective value */
      ) const;

   /** counts the number of variables with negative reduced cost */
   int countPricedVariables(
      PricingType*          pricetype,          /**< pricing type, farkas or redcost */
      int&                  prob,               /**< number of the pricing problem */
      GCG_COL**             cols,               /**< columns corresponding to solutions which should be investigated */
      int                   ncols               /**< number of columns */
   ) const;

   /** return TRUE or FALSE whether the master LP is solved to optimality */
   SCIP_Bool isMasterLPOptimal() const;


   /** return TRUE or FALSE whether pricing problem has been solved to optimality */
   SCIP_Bool  isPricingOptimal(
      SCIP*                 scip,               /**< SCIP data structure */
      SCIP_STATUS           status              /**< status of problem */
   ) const;

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
   SCIP_Bool canPricingBeAborted() const;

   /** sorts pricing problems according to their score */
   void sortPricingProblemsByScore() const;

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
   ) const;

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
      PricingType*          pricetype,          /**< type of pricing: reduced cost or Farkas */
      SCIP_Bool             optimal             /**< heuristic or optimal pricing */
   );

   /** generic method to generate feasible columns from the pricing problem
    * @note This message has to be threadsafe!
    */
   SCIP_RETCODE generateColumnsFromPricingProblem(
      int                   prob,               /**< index of pricing problem */
      PricingType*          pricetype,          /**< type of pricing: reduced cost or Farkas */
      SCIP_Bool             optimal,            /**< should the pricing problem be solved optimal or heuristically */
      SCIP_Real*            lowerbound,         /**< dual bound returned by pricing problem */
      GCG_COL**             cols,               /**< pointer to store columns corresponding to solutions */
      int                   maxcols,            /**< size of the cols array to indicate maximum columns */
      int*                  ncols,              /**< number of columns */
      SCIP_STATUS*          status              /**< solution status of the pricing problem */
   );

   /** solves a specific pricing problem
    * @todo simplify
    * @note This message has to be threadsafe!
    */
   SCIP_RETCODE solvePricingProblem(
      int                   prob,               /**< index of pricing problem */
      PricingType*          pricetype,          /**< type of pricing: reduced cost or Farkas */
      SCIP_Bool             optimal,            /**< should the pricing problem be solved optimal or heuristically */
      SCIP_Real*            lowerbound,         /**< dual bound returned by pricing problem */
      GCG_COL**             cols,               /**< pointer to store columns corresponding to solutions */
      int                   maxcols,            /**< size of the cols array to indicate maximum columns */
      int*                  ncols,              /**< number of columns */
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

   /** computes the stack of masterbranch constraints up to the last generic branching node
    * @note This message has to be threadsafe!
    */
   SCIP_RETCODE computeGenericBranchingconssStack(
      PricingType*          pricetype,          /**< type of pricing: reduced cost or Farkas */
      int                   prob,               /**< index of pricing problem */
      SCIP_CONS***          consstack,          /**< stack of branching constraints */
      int*                  nconsstack,         /**< size of the stack */
      SCIP_Real**           consduals           /**< dual values of the masterbranch solutions */
   ) const;

   /** add bounds change from constraint from the pricing problem at this node
    * @note This message has to be threadsafe!
    */
   SCIP_RETCODE addBranchingBoundChangesToPricing(
      int                   prob,               /**< index of pricing problem */
      SCIP_CONS*            branchcons          /**< branching constraints from which bound should applied */
   ) const;

   SCIP_RETCODE checkBranchingBoundChanges(
      int                   prob,               /**< index of pricing problem */
      SCIP_SOL*             sol,                /**< solution to check */
      SCIP_CONS*            branchcons,         /**< branching constraints from which bound should applied */
      SCIP_Bool*            feasible            /**< check whether the solution is feasible */
   ) const;

   /** check bounds change from constraint from the pricing problem at this node
    * @note This message has to be threadsafe!
    */
   SCIP_RETCODE checkBranchingBoundChangesGcgCol(
      GCG_COL*              gcgcol,             /**< gcg column to check */
      SCIP_CONS*            branchcons,         /**< branching constraints from which bound should applied */
      SCIP_Bool*            feasible            /**< check whether the solution is feasible */
   ) const;

};

#endif
