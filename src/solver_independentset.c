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

/**@file   solver_independentset.c
 * @brief  independent set solver for pricing problems
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

#define SOLVER_NAME          "independentset"
#define SOLVER_DESC          "independent set solver for pricing problems"
#define SOLVER_PRIORITY      500

#define SOLVER_ENABLED       TRUE  /**< indicates whether the solver should be enabled */

#define DEFAULT_DENSITY      0.90

struct GCG_SolverData 
{
   SCIP_Real            density;             /**< graph density threshold above which to use solver */
}; 


/*
 * Local methods
 */


/* Cliquer function to preemptively stop execution, currently unused */
/*
static boolean custom_time_function(int level, int i, int n, int max, double cputime, double realtime, clique_options *opts)
{
   // cputime is time in seconds as a double 
   if( cputime > 10.0 )
   {
      return FALSE;
   }
   else
   {
      SCIPdebugMessage("Current runtime/weight: %g/%d\n",cputime,max);
      return TRUE;
   }
}
*/
/** solve the pricing problem as an independent set problem, in an approximate way */
static
SCIP_RETCODE solveIndependentSet(
   SCIP_Bool             exactly,            /**< should the pricing problem be solved to optimality or heuristically? */
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   GCG_SOLVERDATA*       solver,             /**< solver data structure */
   int                   probnr,             /**< problem number */
   SCIP_Real*            lowerbound,         /**< pointer to store lower bound */
   GCG_COL**             cols,               /**< array of columns corresponding to solutions */
   int                   maxcols,            /**< size of preallocated array */
   int*                  ncols,              /**< pointer to store number of columns */
   SCIP_STATUS*          result              /**< pointer to store pricing problem status */
   )
{ /*lint -e715 */
   SCIP_CONS** constraints;
   SCIP_VAR** consvars;
   SCIP_VAR** indsetvars;
   SCIP_VAR** pricingprobvars;
   graph_t* g;
   FILE* outputfile;
   SCIP_Real* solvals;
   SCIP_Real* consvals;
   SCIP_Real signhelper;
   SCIP_Real biggestobj;
   SCIP_Bool retcode;
   set_t clique;
   clique_options cl_opts;
   int nsolvars;
   int npricingprobvars;
   int nconss;
   int indexcount;
   int unique0;
   int unique1;
   int nodeindex0;
   int nodeindex1;
   int coefindex;
   int scalingfactor;
   int nvars;
   int nedges;
   int debugcounter;
   int indsetconstraintcount;

   assert(pricingprob != NULL);
   assert(solver != NULL);
   assert(lowerbound != NULL);
   assert(cols != NULL);
   assert(ncols != NULL);
   assert(result != NULL);

   pricingprobvars = SCIPgetVars(pricingprob);
   npricingprobvars = SCIPgetNVars(pricingprob);

   constraints = SCIPgetConss(pricingprob);
   nconss = SCIPgetNConss(pricingprob);

   SCIP_CALL( SCIPallocBufferArray(pricingprob, &indsetvars,npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &solvals, npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &consvars, npricingprobvars) );

   /* Initialize arrays to ensure data consistency */
   for( int i = 0; i < npricingprobvars; ++i )
   {
      indsetvars[i] = NULL;
      solvals[i] = -1.0; /* To later determine whether a variable was constrained */
      consvars[i] = NULL;
   }

   /* For pricer */
   nsolvars = npricingprobvars;

   /* Used to keep track of node indizes for bijection while building the graph */
   indexcount = 0;

   /* Used to keep track of whether a variable of an IS constraint was seen before or not */
   unique0 = 1;
   unique1 = 1;

   /* Used to keep track of the index of IS variables in the indsetvars array to later delete edges */
   nodeindex0 = 0;
   nodeindex1 = 0;

   /* Used to determine the kind of non-IS constraints */ 
   coefindex = -1;

   /* Counters for debugging purposes */
   indsetconstraintcount = 0;
   debugcounter = 0;

   /* All variables of the problem are expected to be binary */
   if( SCIPgetNBinVars(pricingprob) < npricingprobvars ){
      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   /* Cliquer library explicitly asks for the node weights to be positive integers.
    * Additionally, the sum of node weights to be smaller than INT_MAX.
    * We restrict our scaling factor to always honor this constraint
    */
   scalingfactor = ((int) INT_MAX / npricingprobvars) - npricingprobvars;

   /* All objective values have to be negative or 0 - library restriction */
   /* During the test we also check for the biggest objective value */
   biggestobj = 0.0;
   for( int i = 0; i < npricingprobvars; ++i )
   {
      signhelper = SCIPvarGetObj(pricingprobvars[i]);
      if( SCIPisLT(pricingprob,0,signhelper) )
      {
         *result = SCIP_STATUS_UNKNOWN;
         return SCIP_OKAY;
      }
      if( SCIPisLT(pricingprob,signhelper,biggestobj) )
      {
         biggestobj = signhelper;
         //SCIPprintVar(pricingprob, pricingprobvars[i], NULL);
      }
      //SCIPdebugMessage("Objective value %g \n", signhelper);
   }
   if( SCIPisLT(pricingprob,biggestobj,-1.0) )
   {
      /* Ensure that INT_MAX is never reached by the sum of all scaled weights */
      scalingfactor = abs((int) scalingfactor / biggestobj);
   }

   /* Build complementary graph by first creating a complete graph and then deleting edges of IS constraints */
   /* Size is first chosen to be maximal and then later cropped down to the actual number of nodes */
   g = graph_new(npricingprobvars);
   for( int i = 0; i < npricingprobvars; ++i )
   {
      for( int j = 0; j < npricingprobvars; ++j )
      {
         if( i != j )
         {
            GRAPH_ADD_EDGE(g,i,j);
         }
         
      }
   }

   /* Main loop to check the nature of each constraint */
   for( int i = 0; i < nconss; ++i )
   {
      assert(constraints[i] != NULL);
      consvars = SCIPgetVarsLinear(pricingprob, constraints[i]);

      /* Check if we have an IS constraint */
      if( SCIPgetNVarsLinear(pricingprob, constraints[i]) == 2 &&
          SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),(SCIP_Real)1) ) 
      {
         indsetconstraintcount++;
         /* Preprocessing: Constraint is only relevant for pricing if one of the variables has an objective value != 0 */
         if( SCIPvarGetObj(consvars[0]) != 0 || SCIPvarGetObj(consvars[1]) != 0 )
         {
            /* Check if the Vars are already in our array of unique vars */
            for( int j = 0; j < npricingprobvars; ++j )
            {
               if( consvars[0] == indsetvars[j] )
               {
                  unique0 = 0;
                  nodeindex0 = j;
               }
               if( consvars[1] == indsetvars[j] ) 
               {
                  unique1 = 0;
                  nodeindex1 = j;
               }
            }
            /* Add Vars to our array of unique vars if not yet included */
            if( unique0 )
            {
               indsetvars[indexcount] = consvars[0];
               nodeindex0 = indexcount;
               /*
               SCIPdebugMessage("Node Index: %d\n",indexcount);
               SCIPdebugMessage("Index: %d \n",SCIPvarGetProbindex(indsetvars[indexcount]));
               SCIPdebugMessage("Weight: %d \n", 1+abs((int) (scalingfactor * SCIPvarGetObj(indsetvars[indexcount]))));
               */
               g->weights[nodeindex0] = 1 + abs((int) (scalingfactor * SCIPvarGetObj(indsetvars[indexcount])));
               ++indexcount;
            }
            if( unique1 )
            {
               indsetvars[indexcount] = consvars[1];
               nodeindex1 = indexcount;
               /*
               SCIPdebugMessage("Node Index: %d\n",indexcount);
               SCIPdebugMessage("Index: %d\n",SCIPvarGetProbindex(indsetvars[indexcount]));
               SCIPdebugMessage("Weight: %d \n", 1+abs((int) (scalingfactor * SCIPvarGetObj(indsetvars[indexcount]))));
               */
               g->weights[nodeindex1] = 1 + abs((int) (scalingfactor * SCIPvarGetObj(indsetvars[indexcount])));
               ++indexcount;
            }
            unique0 = 1;
            unique1 = 1;
            
            if( GRAPH_IS_EDGE(g,nodeindex0,nodeindex1) )
            {
               GRAPH_DEL_EDGE(g,nodeindex0,nodeindex1);
            }
         }
      }
      else
      {
         /* The current constraint is no IS constraint */
         SCIPgetConsNVars(pricingprob,constraints[i],&nvars,&retcode);
         consvals = SCIPgetValsLinear(pricingprob, constraints[i]);
         consvars = SCIPgetVarsLinear(pricingprob, constraints[i]);

         /* Check the coefficients of the variables in the constraint */
         for( int k = 0; k < nvars; ++k )
         {
            //SCIPdebugMessage("Index: %d \n", SCIPvarGetIndex(consvars[k]));
            if( consvals[k] != 1 && (coefindex == -1) )
            {
               //indSetBound = abs((int) consvals[k]);
               //SCIPdebugMessage("\n\nCoefficient: %d \n\n", abs(consvals[k]));
               coefindex = k;
            }
            else if( consvals[k] != 1 && coefindex != -1 )
            {
               /* More than one variable has a coefficient unequal to 1 */

               
               SCIPdebugMessage("Unknown constraint in problem.");
               outputfile = fopen("output.txt", "a+");
               SCIPprintCons(pricingprob, constraints[i], outputfile);
               fputs("\n\n",outputfile);
               fclose(outputfile);
               

               SCIPfreeBufferArray(pricingprob,&indsetvars);
               SCIPfreeBufferArray(pricingprob,&solvals);
               SCIPfreeBufferArray(pricingprob,&consvars);
               graph_free(g);
               *result = SCIP_STATUS_UNKNOWN;
               return SCIP_OKAY;
            }
         }
         /* Check if we have a clique constraint (rhs 1 and coefficients 1) */
         if( !(coefindex == -1) && SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),(SCIP_Real)1) )
         {
            /* Delete the edges between all the variables of the constraint. 
               This way, at most one can be part of the maximum clique */
            for( int j = 0; j < nvars; ++j )
            {
               /* We are only interested in vars potentially relevant for pricing (!= 0) */
               if( SCIPvarGetObj(consvars[j]) != 0 )
               {
                  /* Determine nodeindex0 */
                  for( int k = 0; k < npricingprobvars; ++k )
                  {
                     if( consvars[j] == indsetvars[k] )
                     {
                        nodeindex0 = k;
                        unique0 = 0;
                     }
                  }
                  if( unique0 )
                  {
                     indsetvars[indexcount] = consvars[j];
                     nodeindex0 = indexcount;
                     g->weights[nodeindex0] = 1 + abs((int) (scalingfactor * SCIPvarGetObj(indsetvars[indexcount])));
                     ++indexcount;
                  }
                  /* Determine nodeindex1 */
                  for( int l = 0; l < nvars; ++l )
                  {
                     for( int k = 0; k < npricingprobvars; ++k )
                     {
                        if( consvars[l] == indsetvars[k] )
                        {
                           nodeindex1 = k;
                           unique1 = 0;
                        }
                     }
                     if( unique1 )
                     {
                        indsetvars[indexcount] = consvars[j];
                        nodeindex1 = indexcount;
                        g->weights[nodeindex1] = 1 + abs((int) (scalingfactor * SCIPvarGetObj(indsetvars[indexcount])));
                        ++indexcount;
                     }
                     /* Delete the edge between nodeindex0 and nodeindex1 */
                     if( nodeindex0 != nodeindex1 )
                     {
                        if( GRAPH_IS_EDGE(g,nodeindex0,nodeindex1) )
                        {
                           GRAPH_DEL_EDGE(g,nodeindex0,nodeindex1);
                        }
                     }
                     unique1 = 1;
                  }
               }
               unique0 = 1;
            }
         }
         /* Check if we have a coupling constraint (rhs 0) */
         else if( !(coefindex == -1) && SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),(SCIP_Real)0) )
         {
            /* Special case: The coupling constraint is purely decorative (coefficient >= #vars)*/
            if( abs(consvals[coefindex]) + 1 >= nvars )
            {
               solvals[SCIPvarGetProbindex(consvars[coefindex])] = 1.0;
               //SCIPdebugMessage("\n\nCoupling constraint redundant.\n\n");
            }
            /* Special case: The coefficient is -1, we set the coupling variable to 1 and treat the case like 
               a clique constraint. */
            else if( abs(consvals[coefindex]) == 1 )
            {
               solvals[SCIPvarGetProbindex(consvars[coefindex])] = 1.0;
               /* Delete the edges between all the variables of the constraint that are not the coupling variable.
               This way, at most one can be part of the maximum clique*/
               for(int j = 0; j < nvars; ++j)
               {
                  /* We are only interested in vars potentially relevant for pricing (!= 0) */
                  if( SCIPvarGetObj(consvars[j]) != 0 )
                  {
                     /* Determine nodeindex0 */
                     for( int k = 0; k < npricingprobvars; ++k )
                     {
                        if( consvars[j] == indsetvars[k] || coefindex == k )
                        {
                           nodeindex0 = k;
                           unique0 = 0;
                        }
                     }
                     if( unique0 )
                     {
                        indsetvars[indexcount] = consvars[j];
                        nodeindex0 = indexcount;
                        g->weights[nodeindex0] = 1 + abs((int) (scalingfactor * SCIPvarGetObj(indsetvars[indexcount])));
                        ++indexcount;
                     }
                     /* Determine nodeindex1 */
                     for( int l = 0; l < nvars; ++l )
                     {
                        for( int k = 0; k < npricingprobvars; ++k )
                        {
                           if( consvars[l] == indsetvars[k] || coefindex == k )
                           {
                              nodeindex1 = k;
                              unique1 = 0;
                           }
                        }
                        if( unique1 )
                        {
                           indsetvars[indexcount] = consvars[j];
                           nodeindex1 = indexcount;
                           g->weights[nodeindex1] = 1 + abs((int) (scalingfactor * SCIPvarGetObj(indsetvars[indexcount])));
                           ++indexcount;
                        }
                        if( nodeindex0 != nodeindex1 && nodeindex0 != coefindex && nodeindex1 != coefindex )
                        {
                           if( GRAPH_IS_EDGE(g,nodeindex0,nodeindex1) )
                           {
                              GRAPH_DEL_EDGE(g,nodeindex0,nodeindex1);
                           }
                        }
                        unique1 = 1;
                     }
                     unique0 = 1;
                  }
               }
            }
            else
            {
               /* Coupling coefficient is between 1 and npricingprobvars. */

               
               SCIPdebugMessage("Unknown constraint in problem.");
               outputfile = fopen("output.txt", "a+");
               SCIPprintCons(pricingprob, constraints[i], outputfile);
               fputs("\n\n",outputfile);
               fclose(outputfile);
               
               
               SCIPfreeBufferArray(pricingprob,&indsetvars);
               SCIPfreeBufferArray(pricingprob,&solvals);
               SCIPfreeBufferArray(pricingprob,&consvars);
               graph_free(g);
               set_free(clique);

               *result = SCIP_STATUS_UNKNOWN;
               return SCIP_OKAY;
            }
         }
         else{
            /* Constraint is neither a coupling nor a clique constraint */

            SCIPdebugMessage("Unknown constraint in problem.");
            outputfile = fopen("output.txt", "a+");
            SCIPprintCons(pricingprob, constraints[i], outputfile);
            fputs("\n\n",outputfile);
            fclose(outputfile);

               
            SCIPfreeBufferArray(pricingprob,&indsetvars);
            SCIPfreeBufferArray(pricingprob,&solvals);
            SCIPfreeBufferArray(pricingprob,&consvars);
            graph_free(g);
            set_free(clique);

            *result = SCIP_STATUS_UNKNOWN;
            return SCIP_OKAY;
         }
         coefindex = -1;
      }
   }

   /* Assert that the graph was built in a proper way */ 
   ASSERT(graph_test(g,NULL));

   /* Determine number of edges for graph density calculation */
   nedges = 0;
   for ( int i = 0; i < g->n; i++ )
   {
      for ( int j = 0; j < g->n; j++ )
      {
         if( SET_CONTAINS_FAST(g->edges[i],j) )
         {
            nedges++;
         }
      }
   }

   nedges /= 2;

   /* Test if the density criteria is met */
   SCIPdebugMessage("Graph density in current round: %g%% \n", (float)nedges/((float)(g->n - 1)*(g->n)/2));
   if( SCIPisLT(pricingprob, (float)nedges/((float)(g->n - 1) * (g->n) / 2), solver->density) )
   {
      SCIPdebugMessage("Density below %g%% (%g%%) \n", solver->density, (float)nedges/((float)(g->n - 1)*(g->n)/2));
      SCIPfreeBufferArray(pricingprob,&indsetvars);
      SCIPfreeBufferArray(pricingprob,&solvals);
      SCIPfreeBufferArray(pricingprob,&consvars);
      graph_free(g);

      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;      
   }
   
   /* indexcount now holds the actual number of unique IS variables, thus we truncate */
   assert(indexcount == npricingprobvars-1);

   //SCIPdebugMessage("graphweighted before resize = %u\n", graph_weighted(g));
   if( indexcount > 0 ){
      graph_resize(g,indexcount);
   }
   //SCIPdebugMessage("resized graph, indexcount = %d\n", indexcount);
   //SCIPdebugMessage("graphweighted after resize = %u\n", graph_weighted(g));

   /* Set cliquer options */
   cl_opts.reorder_function=reorder_by_default; //default: reorder_by_default
   cl_opts.reorder_map=NULL;
   cl_opts.time_function=NULL; //default: clique_print_time
   cl_opts.output=NULL;
   cl_opts.user_function=NULL;
   cl_opts.user_data=NULL;
   cl_opts.clique_list=NULL;
   cl_opts.clique_list_length=0;

   //graph_print(g);
   if( biggestobj == 0 )
   {
      /* TODO: This call is of no real value atm since all objectives are 0 */
      clique = clique_unweighted_find_single(g,0,0,FALSE,&cl_opts);
   }
   else
   {
      clique = clique_find_single(g,0,0,FALSE,&cl_opts);
   }
   //SCIPdebugMessage("found clique of size %d\n", set_size(clique));
   for( int i = 0; i < indexcount; ++i )
   {
      if( SET_CONTAINS(clique,i) && SCIPvarGetObj(indsetvars[i]) != 0 )
      {
         //SCIPdebugMessage("Objective value: %g \n",SCIPvarGetObj(indsetvars[i]));
         //SCIPdebugMessage("Var with index %d is in the max IS\n", SCIPvarGetProbindex(indsetvars[i]));
         //SCIPdebugMessage("Node Index: %d\n",i);
         solvals[SCIPvarGetProbindex(indsetvars[i])] = 1.0;
         ++debugcounter;         
      }
      else
      {
         solvals[SCIPvarGetProbindex(indsetvars[i])] = 0.0;
      }
   }

   /* There may be variables left which are unconstrained. We set these to 1 manually if they have an objective value != 0*/
   for( int i = 0; i < npricingprobvars; ++i )
   {
      if( solvals[i] < 0 )
      {
         if( SCIPvarGetObj(pricingprobvars[i]) != 0 )
         {
            solvals[i] = 1.0;
         }
         else
         {
            solvals[i] = 0.0;
         }
         
      }
   }

   //SCIPdebugMessage("Clique Weight: %d \n",clique_max_weight(g,&cl_opts));
   SCIP_CALL( GCGcreateGcgCol(pricingprob, &cols[0], probnr, pricingprobvars, solvals, nsolvars, FALSE, SCIPinfinity(pricingprob)) );
   *ncols = 1;
   *result = SCIP_STATUS_OPTIMAL;

   /*
   SCIPdebugMessage("Biggest objective value: %g \n", biggestobj);
   SCIPdebugMessage("Total number of vars in the max IS: %d \n", debugcounter);
   SCIPdebugMessage("Total number of unique variables: %d \n", npricingprobvars);
   SCIPdebugMessage("Total number of unique IS variables: %d \n", indexcount);
   SCIPdebugMessage("%g%% of the variables are IS variables \n", (((double)indexcount)/npricingprobvars)*100);
   SCIPdebugMessage("Total number of constraints: %d \n", nconss);
   SCIPdebugMessage("Total number of IS constraints: %d \n", indsetconstraintcount);
   SCIPdebugMessage("%g%% of the contraints are IS contraints \n\n", (((double)indsetconstraintcount)/nconss)*100);
   */
   SCIPfreeBufferArray(pricingprob,&indsetvars);
   SCIPfreeBufferArray(pricingprob,&solvals);
   SCIPfreeBufferArray(pricingprob,&consvars);
   graph_free(g);
   set_free(clique);
   return SCIP_OKAY;

}

/*
 * Callback methods for pricing problem solver
 */

/** destructor of pricing solver to free user data (called when SCIP is exiting) */
static
GCG_DECL_SOLVERFREE(solverFreeIndependentSet)
{
   GCG_SOLVERDATA* solverdata;

   assert(scip != NULL);
   assert(solver != NULL);

   solverdata = GCGsolverGetSolverdata(solver);
   assert(solverdata != NULL);

   SCIPfreeMemory(scip, &solverdata);

   GCGsolverSetSolverdata(solver, NULL);

   return SCIP_OKAY;
}

#define solverInitsolIndependentSet NULL
#define solverExitsolIndependentSet NULL
#define solverInitIndependentSet NULL
#define solverExitIndependentSet NULL
//#define solverSolveIndependentSet NULL

/** exact solving method for independent set solver */
static
GCG_DECL_SOLVERSOLVE(solverSolveIndependentSet)
{  /*lint --e{715}*/

   GCG_SOLVERDATA* solverdata;

   solverdata = GCGsolverGetSolverdata(solver);
   assert(solverdata != NULL);
   /* solve the independent set problem exactly */
   SCIP_CALL( solveIndependentSet(TRUE, pricingprob, solverdata, probnr, lowerbound, cols, maxcols, ncols, result) );

   return SCIP_OKAY;
}


/** heuristic solving method of independent set solver */
static
GCG_DECL_SOLVERSOLVEHEUR(solverSolveHeurIndependentSet)
{  /*lint --e{715}*/
   GCG_SOLVERDATA* solverdata;

   solverdata = GCGsolverGetSolverdata(solver);
   assert(solverdata != NULL);
   /* solve the independent set problem approximately */
   SCIP_CALL( solveIndependentSet(FALSE, pricingprob, solverdata, probnr, lowerbound, cols, maxcols, ncols, result) );

   return SCIP_OKAY;
}


/** creates the independent set solver for pricing problems and includes it in GCG */
SCIP_RETCODE GCGincludeSolverIndependentSet(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP* origprob;
   GCG_SOLVERDATA* solverdata;

   origprob = GCGmasterGetOrigprob(scip);
   SCIP_CALL( SCIPallocMemory(scip, &solverdata) );

   SCIP_CALL( GCGpricerIncludeSolver(scip, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY, SOLVER_ENABLED, solverSolveIndependentSet,
         solverSolveHeurIndependentSet, solverFreeIndependentSet, solverInitIndependentSet, solverExitIndependentSet,
         solverInitsolIndependentSet, solverExitsolIndependentSet, solverdata) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/independentset/density",
         "graph density threshold above which to use solver",
         &solverdata->density, TRUE, DEFAULT_DENSITY, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
