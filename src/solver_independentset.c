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
 * @author Henri Lotze
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <assert.h>

#include "solver_independentset.h"
#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"
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

#define DEFAULT_DENSITY      0.85

struct GCG_SolverData 
{
   SCIP_Real             density;             /**< graph density threshold above which to use solver */
}; 


/*
 * Local methods
 */
/* Add a variable to the bijection graph g and indsetvars array. Returns the index of the corresponding node in the graph. */
static int indset_addNodeToGraph(int indexcount, int npricingprobvars, int scalingfactor, SCIP_VAR** indsetvars, graph_t* g, SCIP_VAR* consvar)
{
   int unique;
   int nodeindex;
   unique = 1;

   for( int j = 0; j < indexcount; ++j )
   {
      if( consvar == indsetvars[j] )
      {
         /* Var already part of the graph */
         unique = 0;
         nodeindex = j;
      }
   }
   if( unique )
   {
      /* Var not yet part of graph, add it with its corresponding weight */
      indsetvars[indexcount] = consvar;
      if( !(SCIPvarGetObj(indsetvars[indexcount]) > 0) )
      {
         g->weights[indexcount] = 1 + abs((int) (scalingfactor * SCIPvarGetObj(indsetvars[indexcount])));
      }
      else
      {
         g->weights[indexcount] = 1;
      }
      nodeindex = indexcount;
   }
   return nodeindex;
}

/* Get the node index of a given variable in the bijection if mapped, else return -1 */
static int indset_getNodeIndex(SCIP_VAR* var, SCIP_VAR** indsetvars, int indexcount)
{
   for( int j = 0; j < indexcount; ++j )
   {
      if( var == indsetvars[j] )
      {
         return j;
      }
   }
   return -1;
}

/* Basic idea of the heuristic solver: The biggest independent set in a graph corresponds to the biggest clique
 * of the complement graph, for which we use the cliquer library to find it. We therefore transform the variables 
 * into graph nodes and delete the edge between two nodes if there is an independent set constraint involving both. 
 * By doing this, they cannot both be part of the maximum clique and thus not be both part of the independent set.
 * The correspondence between variables and graph nodes is done by a bijection using the indsetvars array:
 * The variable indsetvars[i] is the i-th node of the graph, indexcount keeps track of the next unmapped graph node.
 * Since we want to add a column with the best reduced cost, we take the objective coefficient of variables into
 * account by giving their graph nodes corresponding weights and searching for a weight-maximal clique.
 *
 * This solver is heuristic since the scaling by weight is limited by the cliquer library. In most realistic
 * scenarios, the result of this solver should be optimal.
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
   SCIP_CONS**    constraints;
   SCIP_CONS**    markedconstraints;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR**     lconsvars;
   SCIP_VAR**     vconsvars;
   SCIP_VAR**     mconsvars;
   SCIP_VAR**     indsetvars;
   SCIP_VAR**     pricingprobvars;
   SCIP_SOL*      conssol;
   graph_t*       g;
   SCIP_Real*     solvals;
   SCIP_Real*     consvals;
   SCIP_Real      candidate;
   SCIP_Real      biggestobj;
   SCIP_Bool      retcode;
   set_t          clique;
   clique_options cl_opts;
   int            markedcount;
   int*           markedvars;
   int            npricingprobvars;
   int            nconss;
   int            indexcount;
   int            nodeindex0;
   int            nodeindex1;
   int            nodeindex2;
   int            coefindex;
   int            scalingfactor;
   int            nvars;
   int            nedges;
   int            flag;
   SCIP_RESULT    consresult;

   FILE *outputcons;


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

   /* All variables of the problem are expected to be binary */
   if( SCIPgetNBinVars(pricingprob) < npricingprobvars )
   {
      SCIPdebugMessage("Exit: Nonbinary variables.\n");
      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   /* Cliquer library explicitly asks for the node weights to be positive integers.
    * Additionally, the sum of node weights needs to be smaller than INT_MAX.
    * We restrict our scaling factor to always honor this constraint
    */
   scalingfactor = ((int) INT_MAX / npricingprobvars) - npricingprobvars;

   /* Check for the biggest objective value to safely adjust the scalingfactor */
   biggestobj = 0.0;
   for( int i = 0; i < npricingprobvars; ++i )
   {
      candidate = SCIPvarGetObj(pricingprobvars[i]);
      if( SCIPisLT(pricingprob,candidate,biggestobj) )
      {
         biggestobj = candidate;
      }
   }
   if( SCIPisLT(pricingprob,biggestobj,-1.0) )
   {
      /* Ensure that INT_MAX is never reached by the sum of all scaled weights */
      scalingfactor = abs((int) scalingfactor / biggestobj);
   }

   SCIP_CALL( SCIPallocBufferArray(pricingprob,&markedconstraints,nconss) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob,&indsetvars,npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob,&solvals,npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob,&mconsvars,npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob,&markedvars,npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob,&lconsvars,npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob,&vconsvars,npricingprobvars) );

   /* Initialize arrays to ensure data consistency */
   for( int i = 0; i < npricingprobvars; ++i )
   {
      indsetvars[i] = NULL;
      solvals[i] = -1.0; /* To later determine whether a variable was constrained */
      markedvars[i] = 0.0; /* To later determine if a var has to be set to 0 */
   }
   /*

   /* Used to keep track of node indizes for bijection while building the graph */
   indexcount = 0;

   /* Used to keep track of the index of IS variables in the indsetvars array to later delete edges */
   nodeindex0 = 0;
   nodeindex1 = 0;

   /* Used to determine the kind of non-IS constraints */ 
   coefindex = -1;

   /* Use to handle a rare combination of IS and varbound constraints */
   markedcount = 0;

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
      conshdlr = SCIPconsGetHdlr(constraints[i]);
      assert(conshdlr != NULL);

      /* The constraint may not be of type 'linear' */
      if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
      {
         lconsvars = SCIPgetVarsLinear(pricingprob, constraints[i]);
         /* Check if we have an IS constraint */
         if( SCIPgetNVarsLinear(pricingprob, constraints[i]) == 2 &&
             SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),(SCIP_Real)1) ) 
         {
            /* Preprocessing: Constraint is only relevant for pricing if one of the variables has an objective value != 0 */
            if( SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[0]),0) || SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[1]),0) )
            {
               if( SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[0]),0) )
               {
                  nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, lconsvars[0]);
                  if( nodeindex0 == indexcount )
                  {
                     ++indexcount;
                  }
               }

               if( SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[1]),0) )
               {
                  nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, lconsvars[1]);
                  if( nodeindex1 == indexcount )
                  {
                     ++indexcount;
                  }
               }
               if( SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[0]),0) && SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[1]),0) )
               {
                  if( GRAPH_IS_EDGE(g,nodeindex0,nodeindex1) )
                  {
                     GRAPH_DEL_EDGE(g,nodeindex0,nodeindex1);
                  }
               }
            }
         }
         /* Handle other constraints that behave like IS constraints, i.e. cx+dy<=rhs with c+d>rhs, c>0, d>0 */
         else if( SCIPgetNVarsLinear(pricingprob,constraints[i]) == 2 )
         {
            consvals = SCIPgetValsLinear(pricingprob, constraints[i]);
            if( consvals[0] > 0 && consvals[1] > 0 && SCIPisLT(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),consvals[0] + consvals[1]) 
               && !SCIPisLT(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),consvals[0])
               && !SCIPisLT(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),consvals[1]))
            {
               /* As before, the constraint is only regarded if it is relevant for pricing */
               if( SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[0]),0) )
               {
                  nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, lconsvars[0]);
                  if( nodeindex0 == indexcount )
                  {
                     ++indexcount;
                  }
               }

               if( SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[1]),0) )
               {
                  nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, lconsvars[1]);
                  if( nodeindex1 == indexcount )
                  {
                     ++indexcount;
                  }
               }
               if( SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[0]),0) && SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[1]),0) )
               {
                  if( GRAPH_IS_EDGE(g,nodeindex0,nodeindex1) )
                  {
                     GRAPH_DEL_EDGE(g,nodeindex0,nodeindex1);
                  }
               }
            }
         }
         else
         {
            /* The current constraint is no linear IS constraint */
            SCIPgetConsNVars(pricingprob,constraints[i],&nvars,&retcode);
            consvals = SCIPgetValsLinear(pricingprob,constraints[i]);
            lconsvars = SCIPgetVarsLinear(pricingprob,constraints[i]);

            /* Check the coefficients of the variables in the constraint */
            for( int k = 0; k < nvars; ++k )
            {
               if( consvals[k] != 1 && (coefindex == -1) )
               {
                  coefindex = k;
               }
               else if( consvals[k] != 1 && coefindex != -1 )
               {
                  /* More than one variable has a coefficient unequal to 1 */
                  SCIPdebugMessage("Exit: More than one coefficient unequal 1, Iteration: %d.\n",i);
                  SCIPfreeBufferArray(pricingprob,&markedconstraints);
                  SCIPfreeBufferArray(pricingprob,&indsetvars);
                  SCIPfreeBufferArray(pricingprob,&solvals);
                  SCIPfreeBufferArray(pricingprob,&mconsvars);
                  SCIPfreeBufferArray(pricingprob,&lconsvars);
                  SCIPfreeBufferArray(pricingprob,&vconsvars);
                  SCIPfreeBufferArray(pricingprob,&markedvars);
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
                  if( SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[j]),0) )
                  {
                     /* Determine nodeindex0 */
                     nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, lconsvars[j]);
                     if( nodeindex0 == indexcount )
                     {
                        ++indexcount;
                     }
                     /* Determine nodeindex1 */
                     for( int l = j + 1; l < nvars; ++l )
                     {
                        nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, lconsvars[l]);
                        if( nodeindex1 == indexcount )
                        {
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
                     }
                  }
               }
            }
            /* Check if we have a coupling constraint (rhs 0) */
            else if( !(coefindex == -1) && SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),(SCIP_Real)0) )
            {
               /* Special case: The coupling constraint is purely decorative (coefficient + 1 coupling var >= #vars)*/
               if( abs(consvals[coefindex]) + 1 >= nvars )
               {
                  /* We cannot guarantee that there is no constraint of the form x+CouplingVar <= 1 */
                  /* If the node is part of the maximum clique, it is safe to set it to one, so we simply add it to the graph */
                  
                  nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, lconsvars[coefindex]);
                  if( nodeindex0 == indexcount )
                  {
                     ++indexcount;
                  }
                  /* We additionally have to mark the variable to later set it to one */ 
                  solvals[SCIPvarGetProbindex(lconsvars[coefindex])] = -2.0;
               }
               /* Special case: The coefficient is -1, we treat the case like a clique constraint. */
               else if( abs(consvals[coefindex]) == 1 )
               {
                  /* We cannot guarantee that there is no constraint of the form x+CouplingVar <= 1 */
                  /* If the node is part of the maximum clique, it is safe to set it to one, so we simply add it to the graph */
                  /* We additionally have to mark the variable to later set it to one */ 
                  nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, lconsvars[coefindex]);
                  if( nodeindex0 == indexcount )
                  {
                     ++indexcount;
                  }
                  /* We additionally have to mark the variable to later set it to one */ 
                  solvals[SCIPvarGetProbindex(lconsvars[coefindex])] = -2.0;

                  /* Delete the edges between all the variables of the constraint that are not the coupling variable.
                     This way, at most one can be part of the maximum clique */
                  for( int j = 0; j < nvars; ++j )
                  {
                     /* We are only interested in vars potentially relevant for pricing (!= 0) */
                     if( j != coefindex && SCIPisLT(pricingprob,SCIPvarGetObj(lconsvars[j]),0) )
                     {
                        /* Determine nodeindex0 */
                        nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, lconsvars[j]);
                        if( nodeindex0 == indexcount )
                        {
                           ++indexcount;
                        }
                        /* Determine nodeindex1 */
                        for( int l = j + 1; l < nvars; ++l )
                        {
                           if( l != coefindex )
                           {
                              nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, lconsvars[l]);
                              if( nodeindex1 == indexcount )
                              {
                                 ++indexcount;
                              }
                              if( nodeindex0 != nodeindex1 )
                              {
                                 if( GRAPH_IS_EDGE(g,nodeindex0,nodeindex1) )
                                 {
                                    GRAPH_DEL_EDGE(g,nodeindex0,nodeindex1);
                                 }
                              }
                           }
                        }
                     }
                  }
               }
               else
               {
                  /* Coupling coefficient is between 1 and npricingprobvars. */
                  SCIPdebugMessage("Exit: Coupling coefficient wrong, Iteration: %d.\n",i);
                  SCIPfreeBufferArray(pricingprob,&markedconstraints);
                  SCIPfreeBufferArray(pricingprob,&indsetvars);
                  SCIPfreeBufferArray(pricingprob,&solvals);
                  SCIPfreeBufferArray(pricingprob,&mconsvars);
                  SCIPfreeBufferArray(pricingprob,&lconsvars);
                  SCIPfreeBufferArray(pricingprob,&vconsvars);
                  SCIPfreeBufferArray(pricingprob,&markedvars);
                  graph_free(g);

                  *result = SCIP_STATUS_UNKNOWN;
                  return SCIP_OKAY;
               }
            }
            else{
               /* Constraint is neither a coupling nor a clique constraint */ 
               SCIPdebugMessage("Exit: Unknown constraint, Iteration: %d.\n",i);
               SCIPfreeBufferArray(pricingprob,&markedconstraints);
               SCIPfreeBufferArray(pricingprob,&indsetvars);
               SCIPfreeBufferArray(pricingprob,&solvals);
               SCIPfreeBufferArray(pricingprob,&mconsvars);
               SCIPfreeBufferArray(pricingprob,&lconsvars);
               SCIPfreeBufferArray(pricingprob,&vconsvars);
               SCIPfreeBufferArray(pricingprob,&markedvars);
               graph_free(g);

               *result = SCIP_STATUS_UNKNOWN;
               return SCIP_OKAY;
            }
            coefindex = -1;
         }
      }
      /* Constraint may be of type varbound: lhs <= x + c*y <= rhs */
      else if( strcmp(SCIPconshdlrGetName(conshdlr), "varbound") == 0 )
      {
         vconsvars[0] = SCIPgetVarVarbound(pricingprob,constraints[i]);
         vconsvars[1] = SCIPgetVbdvarVarbound(pricingprob,constraints[i]);
         /* Check value of rhs to be 0 and c to be <= -1 */ 
         if ( SCIPisInfinity(pricingprob, -SCIPgetLhsVarbound(pricingprob,constraints[i])) )
         {
            if( SCIPisEQ(pricingprob,SCIPgetRhsVarbound(pricingprob,constraints[i]),0) )
            {
               if( SCIPisLT(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,constraints[i]),-1) || SCIPisEQ(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,constraints[i]),-1) ) 
               {
                  /* if x may be relevant, add both x and y to graph */
                  if( SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[0]),0) )
                  {
                     if( SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[0]),0) )
                     {
                        nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, vconsvars[0]);
                        if( nodeindex0 == indexcount )
                        {
                           ++indexcount;
                        }
                     }

                     if( SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[1]),0) )
                     {
                        nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, vconsvars[1]);
                        if( nodeindex1 == indexcount )
                        {
                           ++indexcount;
                        }
                     }
                     /* It may be the case, that both the constraints x - y <= 0 and x + y <= 1 are part of the problem */
                     /* Although rare, we later ensure that we do not set x to 1 while y is set to 0 */
                     markedconstraints[markedcount] = constraints[i];
                     ++markedcount;
                     SCIPwarningMessage(pricingprob,"Added one vbd of special type \n");
                  }
                  /* If only y may be relevant, add only y to the graph */
                  else if( SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[1]),0) )
                  {
                     if( SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[1]),0) )
                     {
                        nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, vconsvars[1]);
                        if( nodeindex1 == indexcount )
                        {
                           ++indexcount;
                        }
                     }
                  }
                  /* If none of the nodes are relevant, ignore both since they will be set to 0 */
               }
               else
               {
                  /* Coefficient c of varbound is > -1 and we do not have an IS constraint*/ 
                  SCIPdebugMessage("Exit: Coefficient of Varbound wrong, Iteration: %d, Rhs:%g,Coeff:%g.\n",i,SCIPgetRhsVarbound(pricingprob,constraints[i]),SCIPgetVbdcoefVarbound(pricingprob,constraints[i]));
                  SCIPfreeBufferArray(pricingprob,&markedconstraints);
                  SCIPfreeBufferArray(pricingprob,&indsetvars);
                  SCIPfreeBufferArray(pricingprob,&solvals);
                  SCIPfreeBufferArray(pricingprob,&mconsvars);
                  SCIPfreeBufferArray(pricingprob,&lconsvars);
                  SCIPfreeBufferArray(pricingprob,&vconsvars);
                  SCIPfreeBufferArray(pricingprob,&markedvars);
                  graph_free(g);

                  *result = SCIP_STATUS_UNKNOWN;
                  return SCIP_OKAY;
               }
            }
            /* Rhs of varbound unequal to 0
             * It may still be the case that we have an IS constraint with a non-linear handler
             * The constraint may also be of the form c + 1 > rhs and c < rhs, i.e. a non-standard IS-constraint
             * We treat these cases like a regular IS constraint 
             */
            else if( (SCIPisEQ(pricingprob,SCIPgetRhsVarbound(pricingprob,constraints[i]),1) && SCIPisEQ(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,constraints[i]),1))
                  || (SCIPisLT(pricingprob,SCIPgetRhsVarbound(pricingprob,constraints[i]),SCIPgetVbdcoefVarbound(pricingprob,constraints[i]) + 1) 
                  && SCIPisLT(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,constraints[i]), SCIPgetRhsVarbound(pricingprob,constraints[i]))) )
            {
               /* Preprocessing: Constraint is only relevant for pricing if one of the variables has an objective value != 0 */
               if( SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[0]),0) || SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[1]),0) )
               {
                  if( SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[0]),0) )
                  {
                     nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, vconsvars[0]);
                     if( nodeindex0 == indexcount )
                     {
                        ++indexcount;
                     }
                  }

                  if( SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[1]),0) )
                  {
                     nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, vconsvars[1]);
                     if( nodeindex1 == indexcount )
                     {
                        ++indexcount;
                     }
                  }
                  if( SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[0]),0) && SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[1]),0) )
                  {
                     if( GRAPH_IS_EDGE(g,nodeindex0,nodeindex1) )
                     {
                        GRAPH_DEL_EDGE(g,nodeindex0,nodeindex1);
                     }
                  }
               }            
            }
            else
            {
               /* Rhs of varbound unequal to 0 and no IS constraint*/ 
               SCIPdebugMessage("Exit: Rhs of Varbound unhandled, Iteration: %d, Rhs: %g, Coeff:%g.\n",i,SCIPgetRhsVarbound(pricingprob,constraints[i]),SCIPgetVbdcoefVarbound(pricingprob,constraints[i]));
               SCIPfreeBufferArray(pricingprob,&markedconstraints);
               SCIPfreeBufferArray(pricingprob,&indsetvars);
               SCIPfreeBufferArray(pricingprob,&solvals);
               SCIPfreeBufferArray(pricingprob,&mconsvars);
               SCIPfreeBufferArray(pricingprob,&lconsvars);
               SCIPfreeBufferArray(pricingprob,&vconsvars);
               SCIPfreeBufferArray(pricingprob,&markedvars);
               graph_free(g);

               *result = SCIP_STATUS_UNKNOWN;
               return SCIP_OKAY;            
            }
         }
         /* We have a varbound constraint of type x + cy == rhs */
         else if( SCIPisEQ(pricingprob, SCIPgetLhsVarbound(pricingprob,constraints[i]), SCIPgetRhsVarbound(pricingprob,constraints[i])) )
         {
            /* If the rhs is 0 and c == -1, both variables have to be set to 0 or to 1 */ 
            if( (SCIPgetRhsVarbound(pricingprob,constraints[i]) == 0) && (SCIPgetVbdcoefVarbound(pricingprob,constraints[i]) == -1) )
            {
               /* Check if adding both variables to the solution would be worth it objective-wise:
                  This is heuristical, as a positive objective won't be weighted in the clique search */
               if( SCIPisLT(pricingprob,SCIPvarGetObj(vconsvars[0]) + SCIPvarGetObj(vconsvars[1]),0) )
               {
                  
                  nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, vconsvars[0]);
                  if( nodeindex0 == indexcount )
                  {
                     ++indexcount;
                  }

                  nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, vconsvars[1]);
                  if( nodeindex1 == indexcount )
                  {
                     ++indexcount;
                  }
                  /* Technically the edge between both can still be deleted, if both cons x-y==0 and x+y<=1 are present. */
                  /* If the edge is deleted, we later force both to be zero */
                  markedconstraints[markedcount] = constraints[i];
                  ++markedcount;
                  SCIPdebugMessage("nodeindex0: %d, nodeindex1: %d,names: %s and %s\n",nodeindex0,nodeindex1,SCIPvarGetName(vconsvars[0]),SCIPvarGetName(vconsvars[1]));
                  
               }
               else
               {  
                  /* Both variables have to be set to 0 to satisfy the constraint. */
                  SCIPdebugMessage("Marked1: %d:%d, Marked2: %d:%d,names: %s and %s\n",nodeindex0,SCIPvarGetProbindex(vconsvars[0]),nodeindex1,SCIPvarGetProbindex(vconsvars[1]),SCIPvarGetName(vconsvars[0]),SCIPvarGetName(vconsvars[1]));
                  markedconstraints[markedcount] = constraints[i];
                  ++markedcount;
                  solvals[SCIPvarGetProbindex(vconsvars[0])] = 0.0;
                  solvals[SCIPvarGetProbindex(vconsvars[1])] = 0.0;
               }
            }
            /*
            // Our right hand side is 1 but setting y to 1 would not satisfy the constraint
            else if( SCIPgetRhsVarbound(pricingprob,constraints[i]) == 1 && SCIPgetVbdcoefVarbound(pricingprob,constraints[i]) != 1 )
            {
               SCIPwarningMessage(pricingprob, "Special1.\n");
               // x has to be set to one to satisfy the constraint
               solvals[SCIPvarGetProbindex(consvars[0])] = 1.0;
               solvals[SCIPvarGetProbindex(consvars[1])] = 0.0;
               markedvars[SCIPvarGetProbindex(consvars[0])] = 1;
               // TODO: x may not be part of max clique due to a seperate IS-constraint
            }
            // Our right hand side is 1 and setting y to 1 would already satisfy the constraint 
            else if( SCIPisEQ(pricingprob,SCIPgetRhsVarbound(pricingprob,constraints[i]),1) && SCIPisEQ(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,constraints[i]),1) )
            {
               // One of the variables has to be set to one to satisfy the constraint but not both 
               // This is done by adding both to the graph and simply deleting the edge between them 
               // TODO: it may still be that neither is part of the max ind set 
               nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[0]);
               if( nodeindex0 == indexcount )
               {
                  ++indexcount;
               }

               nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[1]);
               if( nodeindex1 == indexcount )
               {
                  ++indexcount;
               }

               if( GRAPH_IS_EDGE(g,nodeindex0,nodeindex1) )
               {
                  GRAPH_DEL_EDGE(g,nodeindex0,nodeindex1);
               }
               SCIPwarningMessage(pricingprob, "Special2.\n");
            }
            */
            else
            {
               /* RHS is unequal 0 and unequal 1 */
               SCIPdebugMessage("Exit: Unhandled equality constraint, c: %g, rhs: %g.\n", SCIPgetVbdcoefVarbound(pricingprob,constraints[i]), SCIPgetRhsVarbound(pricingprob,constraints[i]));
               SCIPfreeBufferArray(pricingprob,&markedconstraints);
               SCIPfreeBufferArray(pricingprob,&indsetvars);
               SCIPfreeBufferArray(pricingprob,&solvals);
               SCIPfreeBufferArray(pricingprob,&mconsvars);
               SCIPfreeBufferArray(pricingprob,&lconsvars);
               SCIPfreeBufferArray(pricingprob,&vconsvars);
               SCIPfreeBufferArray(pricingprob,&markedvars);
               graph_free(g);

               *result = SCIP_STATUS_UNKNOWN;
               return SCIP_OKAY;   
            }
         }
         else
         {
            /* We have a varbound of type lhs <= x + c*y */
            SCIPdebugMessage("Exit: Varbound of type lhs <= x+c*y, c: %g, rhs: %g.\n", SCIPgetVbdcoefVarbound(pricingprob,constraints[i]), SCIPgetRhsVarbound(pricingprob,constraints[i]));
            SCIPdebugMessage("Constraint handler: %s\n", SCIPconshdlrGetName(conshdlr));
            SCIPfreeBufferArray(pricingprob,&markedconstraints);
            SCIPfreeBufferArray(pricingprob,&indsetvars);
            SCIPfreeBufferArray(pricingprob,&solvals);
            SCIPfreeBufferArray(pricingprob,&mconsvars);
            SCIPfreeBufferArray(pricingprob,&lconsvars);
            SCIPfreeBufferArray(pricingprob,&vconsvars);
            SCIPfreeBufferArray(pricingprob,&markedvars);
            graph_free(g);

            *result = SCIP_STATUS_UNKNOWN;
            return SCIP_OKAY;   
         }
      }
      else
      {
         /* Constraint handler neither linear nor varbound */
         SCIPdebugMessage("Exit: Unhandled constraint handler, Iteration: %d.\n",i);
         SCIPdebugMessage("Constraint handler: %s\n", SCIPconshdlrGetName(conshdlr));
         SCIPfreeBufferArray(pricingprob,&markedconstraints);
         SCIPfreeBufferArray(pricingprob,&indsetvars);
         SCIPfreeBufferArray(pricingprob,&solvals);
         SCIPfreeBufferArray(pricingprob,&mconsvars);
         SCIPfreeBufferArray(pricingprob,&lconsvars);
         SCIPfreeBufferArray(pricingprob,&vconsvars);
         SCIPfreeBufferArray(pricingprob,&markedvars);
         graph_free(g);

         *result = SCIP_STATUS_UNKNOWN;
         return SCIP_OKAY;         
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
   if( SCIPisLT(pricingprob, (float)nedges/((float)(g->n - 1) * (g->n) / 2), solver->density) )
   {
      SCIPdebugMessage("Exit: Density criteria not met,density: %g.\n",(float)nedges/((float)(g->n - 1) * (g->n) / 2));
      SCIPfreeBufferArray(pricingprob,&markedconstraints);
      SCIPfreeBufferArray(pricingprob,&indsetvars);
      SCIPfreeBufferArray(pricingprob,&solvals);
      SCIPfreeBufferArray(pricingprob,&mconsvars);
      SCIPfreeBufferArray(pricingprob,&lconsvars);
      SCIPfreeBufferArray(pricingprob,&vconsvars);
      SCIPfreeBufferArray(pricingprob,&markedvars);
      graph_free(g);

      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Graph size: %d.\n", indexcount);
   ASSERT( indexcount <= npricingprobvars );

   /* indexcount now holds the actual number of unique IS variables, thus we truncate */
   if( indexcount > 0 )
   {
      graph_resize(g,indexcount);
   }

   /* We reuse the markedvars array to ensure consistency between equality constraints */
   for( int i = 0; i < npricingprobvars; ++i )
   {
      markedvars[i] = 0.0; /* To later determine the solvals of equality linked variables */
   }

   /* Handle the case of varbound-IS combination */
   // Handle equality-Vbd-constraints
   for( int i = 0; i < markedcount; ++i )
   {
      mconsvars[0] = SCIPgetVarVarbound(pricingprob,markedconstraints[i]);
      mconsvars[1] = SCIPgetVbdvarVarbound(pricingprob,markedconstraints[i]);
      nodeindex0 = indset_getNodeIndex(mconsvars[0],indsetvars,indexcount);
      nodeindex1 = indset_getNodeIndex(mconsvars[1],indsetvars,indexcount);

      // Marked vars can be both equality and inequality constraints, here we only handle the equality ones.
      if( SCIPisEQ(pricingprob, SCIPgetLhsVarbound(pricingprob,markedconstraints[i]), SCIPgetRhsVarbound(pricingprob,markedconstraints[i])) )
      {
         // First handle general constraint checks, i.e. those that have to hold independent of which node is part of the graph
         for( int j = 0; j < nconss; ++j )
         {
            conshdlr = SCIPconsGetHdlr(constraints[j]);
            if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
            {
               if( SCIPgetNVarsLinear(pricingprob, constraints[j]) > 2 )
               {
                  // If both variables are part of a clique constraint, not both can become 1, so they both have to become 0
                  // The current constraint is no linear IS constraint
                  SCIPgetConsNVars(pricingprob,constraints[j],&nvars,&retcode);
                  consvals = SCIPgetValsLinear(pricingprob, constraints[j]);
                  lconsvars = SCIPgetVarsLinear(pricingprob, constraints[j]);

                  for( int k = 0; k < nvars; ++k )
                  {
                     if( consvals[k] != 1 && (coefindex == -1) )
                     {
                        coefindex = k;
                     }
                  }
                  /* Check if we have a clique constraint (rhs 1 and coefficients 1) */
                  if( !(coefindex == -1) && SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[j]),(SCIP_Real)1) )
                  {
                     for( int k = 0; k < nvars; ++k )
                     {
                        if( mconsvars[0] == lconsvars[k] )
                        {
                           for( int l = 0; l < nvars; ++l )
                           {
                              if( mconsvars[1] == lconsvars[l] )
                              {
                                 solvals[SCIPvarGetProbindex(lconsvars[0])] = 0.0;
                                 solvals[SCIPvarGetProbindex(lconsvars[1])] = 0.0;                       
                              }
                           }
                        }
                     }
                  } 
                  //Case Coupling constraint (rhs 0)
                  else if( !(coefindex == -1) && SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[j]),(SCIP_Real)0) )
                  {
                     // We only have to handle the case that one of the marked vars is the coupling variable.
                     if( mconsvars[0] == lconsvars[coefindex] || mconsvars[1] == lconsvars[coefindex] )
                     {
                        //@TODO for a later version: Handle this case.
                        SCIPdebugMessage("Exit: A variable of an equality constraint is the coupling variable.");
                        SCIPfreeBufferArray(pricingprob,&markedconstraints);
                        SCIPfreeBufferArray(pricingprob,&indsetvars);
                        SCIPfreeBufferArray(pricingprob,&solvals);
                        SCIPfreeBufferArray(pricingprob,&mconsvars);
                        SCIPfreeBufferArray(pricingprob,&lconsvars);
                        SCIPfreeBufferArray(pricingprob,&vconsvars);
                        SCIPfreeBufferArray(pricingprob,&markedvars);
                        graph_free(g);

                        *result = SCIP_STATUS_UNKNOWN;
                        return SCIP_OKAY;
                     }
                  }
               }
            }
            else if( strcmp(SCIPconshdlrGetName(conshdlr), "varbound") == 0 )
            {
               // Constraint has a varbound handler
               vconsvars[0] = SCIPgetVarVarbound(pricingprob,constraints[j]);
               vconsvars[1] = SCIPgetVbdvarVarbound(pricingprob,constraints[j]);
               nodeindex2 = indset_getNodeIndex(vconsvars[0],indsetvars,indexcount);
               if( SCIPisEQ(pricingprob, SCIPgetLhsVarbound(pricingprob,constraints[j]), SCIPgetRhsVarbound(pricingprob,constraints[j])) )
               {
                  /* If the rhs is 0 and c == -1, both variables have to be set to 0 or to 1 */ 
                  if( SCIPgetRhsVarbound(pricingprob,constraints[j]) == 0 && SCIPgetVbdcoefVarbound(pricingprob,constraints[j]) == -1 )
                  {
                     // We should assure that we are not handling the marked constraint, since we can't do anything with it.
                     if( markedconstraints[i] != constraints[j] )
                     {

                        if( mconsvars[0] == vconsvars[0] || mconsvars[0] == vconsvars[1] || mconsvars[1] == vconsvars[0] || mconsvars[1] == vconsvars[1] )
                        {
                           markedvars[SCIPvarGetProbindex(vconsvars[0])] = 1;
                           markedvars[SCIPvarGetProbindex(vconsvars[1])] = 1;
                           markedvars[SCIPvarGetProbindex(mconsvars[0])] = 1;
                           markedvars[SCIPvarGetProbindex(mconsvars[1])] = 1;
                        }
                     }
                  }
               }
               else if( SCIPgetRhsVarbound(pricingprob,constraints[j]) == 1 && SCIPgetVbdcoefVarbound(pricingprob,constraints[j]) != 1 )
               {
                  if( mconsvars[0] == vconsvars[1] || mconsvars[1] == vconsvars[1] )
                  {
                     // The constraint would be violated if the marked variables are assigned to 1, so we force them to 0
                     solvals[SCIPvarGetProbindex(mconsvars[0])] = 0.0;
                     solvals[SCIPvarGetProbindex(mconsvars[1])] = 0.0;
                  }
               }
            }
            // One of the variables is not part of the graph, i.e. either its objective is >= 0 or it is unconstrained.
            // Case unconstrained: Just set its solval to the val of the other variable of the constraint. (heuristic, may be positive objective and set to 1)
            // Case objective >= 0 and constrained: More complex.
            nvars = 0;
            if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
            {
               // Check if our constraint is an IS constraint
               if( SCIPgetNVarsLinear(pricingprob, constraints[j]) == 2 
               && SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[j]),(SCIP_Real)1) )
               {
                  lconsvars = SCIPgetVarsLinear(pricingprob, constraints[j]);
                  nvars = 2;
               }
            }
            else if( strcmp(SCIPconshdlrGetName(conshdlr), "varbound") == 0 )
            {
               if( (SCIPisEQ(pricingprob,SCIPgetRhsVarbound(pricingprob,constraints[j]),1) && SCIPisEQ(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,constraints[j]),1))
                  || (SCIPisLT(pricingprob,SCIPgetRhsVarbound(pricingprob,constraints[j]),SCIPgetVbdcoefVarbound(pricingprob,constraints[j]) + 1) 
                     && SCIPisLT(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,constraints[j]), SCIPgetRhsVarbound(pricingprob,constraints[j]))) )
               {
                  vconsvars[0] = SCIPgetVarVarbound(pricingprob,constraints[j]);
                  vconsvars[1] = SCIPgetVbdvarVarbound(pricingprob,constraints[j]);
                  nvars = 2;
               }
            }
            if( nvars == 2 )
            {
               if( nodeindex0 != -1 && nodeindex1 != -1 )
               {
                  if( !GRAPH_IS_EDGE(g,nodeindex0,nodeindex1) )
                  {
                     // Not both variables may be 1, so the only option left for satisfiability is to set both to 0.
                     solvals[SCIPvarGetProbindex(mconsvars[0])] = 0.0;
                     solvals[SCIPvarGetProbindex(mconsvars[1])] = 0.0;
                  }
               }
               if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
               {
                  //If the marked var that is not in the graph is set to 1, do we violate a constraint?
                  if( mconsvars[0] == lconsvars[0] || mconsvars[1] == lconsvars[0] )
                  {
                     nodeindex2 = indset_getNodeIndex(lconsvars[1],indsetvars,indexcount);
                     if( nodeindex2 != -1 )
                     {
                           // Not both variables may be 1, so the only option left for satisfiability is to set both to 0.
                           solvals[SCIPvarGetProbindex(mconsvars[0])] = 0.0;
                           solvals[SCIPvarGetProbindex(mconsvars[1])] = 0.0;
                     }
                  }
                  else if( mconsvars[0] == lconsvars[1] || mconsvars[1] == lconsvars[1] )
                  {
                     nodeindex2 = indset_getNodeIndex(lconsvars[0],indsetvars,indexcount);
                     if( nodeindex2 != -1 )
                     {
                        // Not both variables may be 1, so the only option left for satisfiability is to set both to 0.
                        solvals[SCIPvarGetProbindex(mconsvars[0])] = 0.0;
                        solvals[SCIPvarGetProbindex(mconsvars[1])] = 0.0;
                     }
                  }
               }
               else if( strcmp(SCIPconshdlrGetName(conshdlr), "varbound") == 0 )
               {
                  //If the marked var that is not in the graph is set to 1, do we violate a constraint?
                  if( mconsvars[0] == vconsvars[0] || mconsvars[1] == vconsvars[0] )
                  {
                     nodeindex2 = indset_getNodeIndex(vconsvars[1],indsetvars,indexcount);
                     if( nodeindex2 != -1 )
                     {
                        // Not both variables may be 1, so the only option left for satisfiability is to set both to 0.
                        solvals[SCIPvarGetProbindex(mconsvars[0])] = 0.0;
                        solvals[SCIPvarGetProbindex(mconsvars[1])] = 0.0;
                     }
                  }
                  else if( mconsvars[0] == vconsvars[1] || mconsvars[1] == vconsvars[1] )
                  {
                     nodeindex2 = indset_getNodeIndex(vconsvars[0],indsetvars,indexcount);
                     if( nodeindex2 != -1 )
                     {
                        // Not both variables may be 1, so the only option left for satisfiability is to set both to 0.
                        solvals[SCIPvarGetProbindex(mconsvars[0])] = 0.0;
                        solvals[SCIPvarGetProbindex(mconsvars[1])] = 0.0;
                     }
                  }
               }
            }            
         }
      }
   }

   /* Handle linked equality constraints. They all have to agree on a solval of 0 or 1 to be valid. */
   flag = -1;
   for( int j = 0; j < npricingprobvars; ++j )
   {
      if( markedvars[j] == 1 )
      {
         if( solvals[SCIPvarGetProbindex(pricingprobvars[j])] == 0 )
         {
            if( flag == 1)
            {
               SCIPwarningMessage(pricingprob, "Contradiction in forced solvals for equality constraint (found 1 when expecting 0)\n");
            }
            flag = 0;

         }
         else if( solvals[SCIPvarGetProbindex(pricingprobvars[j])] == 1 )
         {
            if( flag == 0)
            {
               SCIPwarningMessage(pricingprob, "Contradiction in forced solvals for equality constraint (found 0 when expecting 1)\n");
            }
            flag = 1;
         }
      }
   }
   // If no single variable has to be set to be set to 0, we may set the variables to 1
   if( flag == -1 )
   {
      flag = 1;
   }
   if( flag == 0 )
   {
      for( int j = 0; j < npricingprobvars; ++j )
      {
         if( markedvars[j] == 1 )
         {
            solvals[SCIPvarGetProbindex(pricingprobvars[j])] = 0.0;
         }
      }
   }
   else if( flag == 1 )
   {
      for( int j = 0; j < npricingprobvars; ++j )
      {
         if( markedvars[j] == 1 )
         {
            solvals[SCIPvarGetProbindex(pricingprobvars[j])] = 1.0;
         }
      }
   }

   /* Clean up the graph. If a variable's solval has been set to 0, it should not be part of the max clique */
   /* We enforce this by isolating the node and setting its weight to 1 */
   for( int i = 0; i < npricingprobvars; ++i )
   {
      if( solvals[SCIPvarGetProbindex(pricingprobvars[i])] == 0 )
      {
         nodeindex0 = indset_getNodeIndex(pricingprobvars[i],indsetvars,indexcount);
         /* The var is part of the graph if its index is unequal to -1 */
         if( nodeindex0 != -1 )
         {
            for( int j = 0; j < indexcount; ++j )
            {
               if( GRAPH_IS_EDGE(g,nodeindex0,j) )
               {
                  GRAPH_DEL_EDGE(g,nodeindex0,j);
               }
            }
            g->weights[nodeindex0] = 1;       
         }
      }
   }

   /* Set cliquer options */
   cl_opts.reorder_function = reorder_by_default; //default: reorder_by_default
   cl_opts.reorder_map = NULL;
   cl_opts.time_function = NULL; //default: clique_print_time
   cl_opts.output = NULL;
   cl_opts.user_function = NULL;
   cl_opts.user_data = NULL;
   cl_opts.clique_list = NULL;
   cl_opts.clique_list_length = 0;

   /* Find maximum weight cliques using the cliquer library*/
   if( biggestobj == 0 )
   {
      clique = clique_unweighted_find_single(g,0,0,FALSE,&cl_opts);
   }
   else
   {
      clique = clique_find_single(g,0,0,FALSE,&cl_opts);
   }

   /* Set all members of the maximum clique with objective coefficient <0 to 1 */
   for( int i = 0; i < indexcount; ++i )
   {
      /* Coupling variables were pre-set to -2.0, if they are part of the maximum clique, we enable them. */
      if( SET_CONTAINS(clique,i) && (SCIPisLT(pricingprob,SCIPvarGetObj(indsetvars[i]),0) || solvals[SCIPvarGetProbindex(indsetvars[i])] == -2.0) 
         && solvals[SCIPvarGetProbindex(indsetvars[i])] != 0.0 )
      {
         solvals[SCIPvarGetProbindex(indsetvars[i])] = 1.0;
         //SCIPdebugMessage("Node %d is Var %s.\n",i,SCIPvarGetName(indsetvars[i]));
      }
      else
      {
         //SCIPdebugMessage("Node %d is Var %s.\n",i,SCIPvarGetName(indsetvars[i]));
         /* We may have set some variables manually already, e.g. coupling variables */
         if( solvals[SCIPvarGetProbindex(indsetvars[i])] != 1.0)
         {
            solvals[SCIPvarGetProbindex(indsetvars[i])] = 0.0;
         }
         
      }
   }

   /* Handle the case of marked inequality constraints of type x - y <= 0 in combination with x + y <= 1 -Constraints */
   for( int i = 0; i < markedcount; ++i )
   {
      if( SCIPisEQ(pricingprob,SCIPgetRhsVarbound(pricingprob,markedconstraints[i]),0) 
         && (SCIPisLT(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,markedconstraints[i]),-1) 
            || SCIPisEQ(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,markedconstraints[i]),-1)) )
      {
         vconsvars[0] = SCIPgetVarVarbound(pricingprob,markedconstraints[i]);
         vconsvars[1] = SCIPgetVbdvarVarbound(pricingprob,markedconstraints[i]);
         if( (solvals[SCIPvarGetProbindex(vconsvars[0])] == 1) && (SCIPvarGetProbindex(vconsvars[1]) == 0) )
         {
            solvals[SCIPvarGetProbindex(vconsvars[0])] = 0;
            solvals[SCIPvarGetProbindex(vconsvars[1])] = 0;
         }
      }
   }

   /* Handle the case that there are still solvals of equality constraints that do not agree.
    * This may occur if one is unset (solval:-1) and the other one is already set (solval 0 or 1)
    */
   for( int i = 0; i < markedcount; ++i )
   {
      vconsvars[0] = SCIPgetVarVarbound(pricingprob,markedconstraints[i]);
      vconsvars[1] = SCIPgetVbdvarVarbound(pricingprob,markedconstraints[i]);

      if( solvals[SCIPvarGetProbindex(vconsvars[0])] != solvals[SCIPvarGetProbindex(vconsvars[1])] 
         && SCIPisEQ(pricingprob, SCIPgetLhsVarbound(pricingprob,markedconstraints[i]), SCIPgetRhsVarbound(pricingprob,markedconstraints[i])) 
         )
      {
         if( solvals[SCIPvarGetProbindex(vconsvars[0])] == 0 || solvals[SCIPvarGetProbindex(vconsvars[1])] == 0 )
         {
            solvals[SCIPvarGetProbindex(vconsvars[0])] = 0.0;
            solvals[SCIPvarGetProbindex(vconsvars[1])] = 0.0;
         }
         else
         {
            // One or both of the vars are unset and the other one, if not -1, is forced to be 1, thus we can set both to 1
            solvals[SCIPvarGetProbindex(vconsvars[0])] = 1.0;
            solvals[SCIPvarGetProbindex(vconsvars[1])] = 1.0;
         }
      }
   }

   /* There may be variables left which are unconstrained. We set these to 1 manually if they have an objective value != 0*/
   for( int i = 0; i < npricingprobvars; ++i )
   {
      if( solvals[i] == -1.0 )
      {
         if( SCIPisLT(pricingprob,SCIPvarGetObj(pricingprobvars[i]),0) )
         {
            solvals[i] = 1.0;
         }
         else
         {
            solvals[i] = 0.0;
         }
      }
      SCIPdebugMessage("Var %d has Name %s and value %g .\n",i,SCIPvarGetName(pricingprobvars[i]),solvals[i]);
   }
   /*
   //BEGIN Debug
   
   SCIP_CALL( SCIPcreateSol(pricingprob,&conssol,NULL) );
   SCIP_CALL( SCIPsetSolVals(pricingprob,conssol,npricingprobvars,pricingprobvars,solvals) );

   //outputgraph = fopen("writegraph.dimacs","w");
   outputcons = fopen("constraints.out","a");
   fprintf(outputcons, "Markedcount: %d\n", markedcount);
   for( int i = 0; i < nconss; ++i )
   {
      SCIPgetConsNVars(pricingprob,constraints[i],&nvars,&retcode);
      conshdlr = SCIPconsGetHdlr(constraints[i]);
      if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
      {
         lconsvars = SCIPgetVarsLinear(pricingprob, constraints[i]);
      }
      else
      {
         vconsvars[0] = SCIPgetVarVarbound(pricingprob,constraints[i]); 
         vconsvars[1] = SCIPgetVbdvarVarbound(pricingprob,constraints[i]);
      }
      
      SCIPprintCons(pricingprob,constraints[i],outputcons);
      fputs(" |",outputcons);
      for( int j = 0; j < nvars; ++j )
      {
         if(strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0)
            fprintf(outputcons," %g", solvals[SCIPvarGetProbindex(lconsvars[j])]);
         else 
            fprintf(outputcons," %g", solvals[SCIPvarGetProbindex(vconsvars[j])]);

         if( solvals[SCIPvarGetProbindex(consvars[j])] == 1 )
         {
            if( strcmp(SCIPvarGetName(consvars[j]), "pr0_t_x#36#1") == 0 || strcmp(SCIPvarGetName(consvars[j]), "pr0_t_x#37#1") == 0)
            {
               if( SCIPisFeasZero(pricingprob, SCIPvarGetUbLocal(consvars[j])) )
               {
                  printf("Var: %s, has ub %g.\n", SCIPvarGetName(consvars[j]),  SCIPvarGetUbLocal(consvars[j]));
               }

               //SCIPprintVar(pricingprob,origvar,NULL);
            }
         }

      }

      fputs(" ||",outputcons);
      for( int j = 0; j < nvars; ++j )
      {
         if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
            fprintf(outputcons," %d", indset_getNodeIndex(lconsvars[j],indsetvars,indexcount));
         else
            fprintf(outputcons," %d", indset_getNodeIndex(vconsvars[j],indsetvars,indexcount));
      }
      if( nvars == 2 )
      {
         if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
         {
            nodeindex0 = indset_getNodeIndex(lconsvars[0],indsetvars,indexcount);
            nodeindex1 = indset_getNodeIndex(lconsvars[1],indsetvars,indexcount);
         }
         else{
            nodeindex0 = indset_getNodeIndex(vconsvars[0],indsetvars,indexcount);
            nodeindex1 = indset_getNodeIndex(vconsvars[1],indsetvars,indexcount);            
         }

         if( nodeindex0 != -1 && nodeindex0 != -1 )
         {
            if( GRAPH_IS_EDGE(g,nodeindex0,nodeindex1))
            {
               fputs(" || edge",outputcons);
            }
         }

      }

      fputs(" ||",outputcons);
      for( int j = 0; j < nvars; ++j )
      {
         if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
            fprintf(outputcons," %g  ", SCIPvarGetObj(lconsvars[j]));
         else
            fprintf(outputcons," %g  ", SCIPvarGetObj(vconsvars[j]));
      }
      SCIP_CALL( SCIPcheckCons(pricingprob,constraints[i],conssol,FALSE,FALSE,TRUE,&consresult) );
      //fprintf(outputcons, "Resultcode: %d", consresult);
      if( consresult != SCIP_FEASIBLE )
      {
         fputs("||| violated",outputcons);
      }
      fputs("\n",outputcons);

   }
   //printf("\n\n");
   fputs("\n\n\n",outputcons);
   //graph_write_dimacs_ascii(g,NULL,outputgraph);
   //fclose(outputgraph);
   fclose(outputcons);

   SCIPdebugMessage("\n\n");
   
   //END DEBUG
   */
   
   /* Create a column corresponding to our clique result */
   SCIP_CALL( GCGcreateGcgCol(pricingprob, &cols[0], probnr, pricingprobvars, solvals, npricingprobvars, FALSE, SCIPinfinity(pricingprob)) );
   *ncols = 1;
   *result = SCIP_STATUS_OPTIMAL;

   SCIPfreeBufferArray(pricingprob,&markedconstraints);
   SCIPfreeBufferArray(pricingprob,&indsetvars);
   SCIPfreeBufferArray(pricingprob,&solvals);
   SCIPfreeBufferArray(pricingprob,&markedvars);
   SCIPfreeBufferArray(pricingprob,&mconsvars);
   SCIPfreeBufferArray(pricingprob,&lconsvars);
   SCIPfreeBufferArray(pricingprob,&vconsvars);
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
#define solverSolveIndependentSet NULL

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