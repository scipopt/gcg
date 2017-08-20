/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       */
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

/**@file   solver_knapsack.c
 * @brief  knapsack solver for pricing problems
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>

#include "solver_independentset.h"
#include "scip/cons_linear.h"
#include "type_solver.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "pub_gcgcol.h"
#include "cliquer.h"
#include "graph.h" //Graph structure from cliquer library


#define SOLVER_NAME          "independentSet"
#define SOLVER_DESC          "independent set solver for pricing problems"
#define SOLVER_PRIORITY      500

#define SOLVER_ENABLED       TRUE  /**< indicates whether the solver should be enabled */

/** knapsack pricing solver needs no solverdata */
/* struct GCG_SolverData {}; */


/*
 * Local methods
 */

/** solve the pricing problem as a knapsack problem, either exactly or approximately */
static
SCIP_RETCODE solveIndependentSet(
   SCIP_Bool             exactly,            /**< should the pricing problem be solved to optimality or heuristically? */
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   GCG_SOLVER*           solver,             /**< solver data structure */
   int                   probnr,             /**< problem number */
   SCIP_Real*            lowerbound,         /**< pointer to store lower bound */
   GCG_COL**             cols,               /**< array of columns corresponding to solutions */
   int                   maxcols,            /**< size of preallocated array */
   int*                  ncols,              /**< pointer to store number of columns */
   SCIP_STATUS*          result              /**< pointer to store pricing problem status */
   )
{ /*lint -e715 */
   SCIP_CONS** constraints;
   SCIP_VAR** consVars;
   SCIP_VAR** indSetVars;
   graph_t *g;
   set_t clique;
   SCIP_Real* solvals;
   SCIP_Real* consvals;
   int nsolvars;
   SCIP_VAR** pricingprobvars;
   SCIP_Real signHelper;
   SCIP_Real biggestObj;
   int npricingprobvars;
   int nconss;
   int indSetConstraintCount;
   int indexCount;
   int unique0;
   int unique1;
   int nodeIndex0;
   int nodeIndex1;
   int coefFlag;
//   int indSetBound;
   int scalingFactor;
   int nvars;
   int debugCounter;
   SCIP_Bool retcode;
//   FILE *outputfile;

   /* check preconditions */
   assert(pricingprob != NULL);
   assert(solver != NULL);
   assert(lowerbound != NULL);
   assert(cols != NULL);
   assert(ncols != NULL);
   assert(result != NULL);

   pricingprobvars = SCIPgetVars(pricingprob);
   npricingprobvars = SCIPgetNVars(pricingprob);

   /* All variables of the problem are expected to be binary */
   if( SCIPgetNBinVars(pricingprob) < npricingprobvars ){
      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   /* Cliquer library explicitly asks for the node weights to be positive integers.
    * Additionally, the sum of node weights to be smaller than INT_MAX.
    * We restrict our scaling factor to always honor this constraint
    */
   scalingFactor = ((int) INT_MAX / npricingprobvars) - npricingprobvars;

   /* All objective values have to be negative or 0 - library restriction */
   /* During the test we also check for the biggest objective value */
   biggestObj = 0.0;
   for ( int i = 0; i < npricingprobvars; ++i )
   {
      signHelper = SCIPvarGetObj(pricingprobvars[i]);
      if( signHelper > 0 )
      {
         *result = SCIP_STATUS_UNKNOWN;
         return SCIP_OKAY;
      }
      if( signHelper < biggestObj )
      {
         biggestObj = signHelper;
      }
   }
   if(biggestObj < (-1.0))
   {
      scalingFactor = abs((int) scalingFactor / biggestObj);
   }

   /* Build complementary graph by first creating a complete graph and then deleting edges of IS constraints */
   /* Size is first chosen to be maximal and then later cropped down to the actual number of nodes */
   g = graph_new(npricingprobvars);
   for ( int i = 0; i < npricingprobvars; ++i )
   {
      for ( int j = 0; j < npricingprobvars; ++j )
      {
         if ( i != j )
         {
            GRAPH_ADD_EDGE(g,i,j);
         }
         
      }
   }
   ASSERT(graph_test(g,stderr));

   constraints = SCIPgetConss(pricingprob);
   nconss = SCIPgetNConss(pricingprob);

   SCIP_CALL( SCIPallocBufferArray(pricingprob, &indSetVars,npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &solvals, npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &consVars, npricingprobvars) );
   nsolvars = npricingprobvars;
   indexCount=0;

   for ( int i = 0; i < npricingprobvars; ++i )
   {
      indSetVars[i] = NULL;
   }

   /* Used to keep track of whether a variable of an IS constraint was seen before or not */
   unique0=1;
   unique1=1;

   /* Used to keep track of the index of IS variables in the indSetVars array to later delete edges */
   nodeIndex0=0;
   nodeIndex1=0;

   /* Used to determine whether we have a proper coupling constraint */ 
   coefFlag = 0;

   /* Used to determine the maximum size of the IS */
   //indSetBound = npricingprobvars;

   indSetConstraintCount = 0;
   for ( int i = 0; i < nconss; ++i )
   {
      assert(constraints[i] != NULL);
      consVars = SCIPgetVarsLinear(pricingprob, constraints[i]);

      /* Check if we have an IS constraint */
      if(   SCIPgetNVarsLinear(pricingprob, constraints[i]) == 2 
         && SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),(SCIP_Real)1) ) 
      {
         indSetConstraintCount++;
         /* Check if the Vars are already in our array of unique vars */
         for ( int j = 0; j < npricingprobvars; ++j )
         {
            if ( consVars[0] == indSetVars[j] )
            {
               unique0 = 0;
               nodeIndex0 = j;
            }
            if ( consVars[1] == indSetVars[j] ) 
            {
               unique1 = 0;
               nodeIndex1 = j;
            }
         }
         /* Add Vars to our array of unique vars if not yet included */
         if(unique0)
         {
            indSetVars[indexCount] = consVars[0];
            nodeIndex0 = indexCount;
            g->weights[nodeIndex0]=1+abs((int) (scalingFactor * SCIPvarGetObj(indSetVars[indexCount])));
            ++indexCount;
         }
         if(unique1)
         {
            indSetVars[indexCount] = consVars[1];
            nodeIndex1 = indexCount;
            g->weights[nodeIndex1]=1+abs((int) (scalingFactor * SCIPvarGetObj(indSetVars[indexCount])));
            ++indexCount;
         }
         unique0 = 1;
         unique1 = 1;
         
         if(GRAPH_IS_EDGE(g,nodeIndex0,nodeIndex1))
         {
            GRAPH_DEL_EDGE(g,nodeIndex0,nodeIndex1);
         }
      }
      else
      {
         SCIPgetConsNVars(pricingprob,constraints[i],&nvars,&retcode);
         consvals = SCIPgetValsLinear(pricingprob, constraints[i]);

         /* Check the coefficients of the variables in the constraint */
         for(int k = 0; k < nvars; ++k)
         {
            if(consvals[k] != 1 && !coefFlag)
            {
               //indSetBound = abs((int) consvals[k]);
               coefFlag = 1;
            }
            else if(consvals[k] != 1 && coefFlag)
            {
               /* More than one variable has a coefficient unequal to 1 */
               *result = SCIP_STATUS_UNKNOWN;
               SCIPfreeBufferArray(pricingprob,&indSetVars);
               SCIPfreeBufferArray(pricingprob,&solvals);
               SCIPfreeBufferArray(pricingprob,&consVars);
               graph_free(g);
               return SCIP_OKAY;
            }
         }
         /* Check if we have a clique constraint (rhs 1 and coefficients 1) */
         if(!coefFlag && SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),(SCIP_Real)1) )
         {
            /* TODO: Modify the graph accordingly */
         }
         /* Check if we have a coupling constraint (rhs 0) */
         else if( SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),(SCIP_Real)0) )
         {
            /* Since the indSetBound was set, we are good here.*/
         }
         else{
            *result = SCIP_STATUS_UNKNOWN;
            SCIPfreeBufferArray(pricingprob,&indSetVars);
            SCIPfreeBufferArray(pricingprob,&solvals);
            SCIPfreeBufferArray(pricingprob,&consVars);
            graph_free(g);
            return SCIP_OKAY;
         }
         /*
         SCIPdebugMessage("Number of vars in non-IS-constraint: %d \n", nvars);
         outputfile = fopen("output.txt", "a+");
         SCIPprintCons(pricingprob, constraints[i], outputfile);
         fputs("\n\n",outputfile);
         fclose(outputfile);
         */
         
      }
   }

   ASSERT(graph_test(g,stderr));
   debugCounter = 0;
   
   /* indexCount now holds the actual number of unique IS variables, thus we truncate */
   graph_resize(g,indexCount);
   if(biggestObj == 0)
   {
      clique = clique_unweighted_find_single(g,0,0,FALSE,NULL);
   }
   else
   {
      clique = clique_find_single(g,0,0,FALSE,NULL);
   }
   for ( int i = 0; i < set_size(clique); ++i )
   {
      if(SET_CONTAINS(clique,i))
      {
         solvals[i] = 1.0;
         //SCIPdebugMessage("Var index %d is in the max IS\n", i);
         ++debugCounter;         
      }
      else
      {
         solvals[i] = 0.0;
      }
   }
   SCIP_CALL( GCGcreateGcgCol(pricingprob, &cols[0], probnr, pricingprobvars, solvals, nsolvars, FALSE, SCIPinfinity(pricingprob)) );
   //*result = SCIP_STATUS_OPTIMAL;
   
   SCIPdebugMessage("Biggest objective value: %g \n", biggestObj);
   SCIPdebugMessage("Total number of vars in the max IS: %d \n", debugCounter);
   SCIPdebugMessage("Total number of unique variables: %d \n", npricingprobvars);
   SCIPdebugMessage("Total number of unique IS variables: %d \n", indexCount);
   SCIPdebugMessage("%g%% of the variables are IS variables \n", (((double)indexCount)/npricingprobvars)*100);
   SCIPdebugMessage("Total number of constraints: %d \n", nconss);
   SCIPdebugMessage("Total number of IS constraints: %d \n", indSetConstraintCount);
   SCIPdebugMessage("%g%% of the contraints are IS contraints \n\n", (((double)indSetConstraintCount)/nconss)*100);

   SCIPfreeBufferArray(pricingprob,&indSetVars);
   SCIPfreeBufferArray(pricingprob,&solvals);
   SCIPfreeBufferArray(pricingprob,&consVars);
   graph_free(g);
   return SCIP_OKAY;

}

/*
 * Callback methods for pricing problem solver
 */


#define solverFreeIndependentSet NULL
#define solverInitsolIndependentSet NULL
#define solverExitsolIndependentSet NULL
#define solverInitIndependentSet NULL
#define solverExitIndependentSet NULL


/** exact solving method for independent set solver */
static
GCG_DECL_SOLVERSOLVE(solverSolveIndependentSet)
{  /*lint --e{715}*/

   /* solve the independent set problem exactly */
   SCIP_CALL( solveIndependentSet(TRUE, pricingprob, solver, probnr, lowerbound, cols, maxcols, ncols, result) );

   return SCIP_OKAY;
}


/** heuristic solving method of independent set solver */
static
GCG_DECL_SOLVERSOLVEHEUR(solverSolveHeurIndependentSet)
{  /*lint --e{715}*/

   /* solve the independent set problem approximately */
   SCIP_CALL( solveIndependentSet(FALSE, pricingprob, solver, probnr, lowerbound, cols, maxcols, ncols, result) );

   return SCIP_OKAY;
}


/** creates the independent set solver for pricing problems and includes it in GCG */
SCIP_RETCODE GCGincludeSolverIndependentSet(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY, SOLVER_ENABLED, solverSolveIndependentSet,
         solverSolveHeurIndependentSet, solverFreeIndependentSet, solverInitIndependentSet, solverExitIndependentSet,
         solverInitsolIndependentSet, solverExitsolIndependentSet, NULL) );

   return SCIP_OKAY;
}
