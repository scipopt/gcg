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

#define DEFAULT_DENSITY      0.90

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

   for( int j = 0; j < npricingprobvars; ++j )
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
      g->weights[indexcount] = 1 + abs((int) (scalingfactor * SCIPvarGetObj(indsetvars[indexcount])));
      nodeindex = indexcount;
   }
   return nodeindex;
}

/* Get the node index of a given variable in the bijection if mapped, else return -1 */
static int indset_getNodeIndex(SCIP_VAR* var, SCIP_VAR** indsetvars, int npricingprobvars)
{
   for( int j = 0; j < npricingprobvars; ++j )
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
   SCIP_VAR**     consvars;
   SCIP_VAR**     indsetvars;
   SCIP_VAR**     pricingprobvars;
   graph_t*       g;
   SCIP_Real*     solvals;
   SCIP_Real*     consvals;
   SCIP_Real      signhelper;
   SCIP_Real      biggestobj;
   SCIP_Bool      retcode;
   set_t          clique;
   clique_options cl_opts;
   int            markedcount;
   int            nsolvars;
   int            npricingprobvars;
   int            nconss;
   int            indexcount;
   int            nodeindex0;
   int            nodeindex1;
   int            coefindex;
   int            scalingfactor;
   int            nvars;
   int            nedges;



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

   /* All objective values have to be negative or 0 - library restriction */
   /* During the test we also check for the biggest objective value */
   biggestobj = 0.0;
   for( int i = 0; i < npricingprobvars; ++i )
   {
      signhelper = SCIPvarGetObj(pricingprobvars[i]);
      if( SCIPisLT(pricingprob,0,signhelper) )
      {
         SCIPdebugMessage("Exit: Wrong coefficient sign.\n");
         *result = SCIP_STATUS_UNKNOWN;
         return SCIP_OKAY;
      }
      if( SCIPisLT(pricingprob,signhelper,biggestobj) )
      {
         biggestobj = signhelper;
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
   SCIP_CALL( SCIPallocBufferArray(pricingprob,&consvars,npricingprobvars) );

   /* Initialize arrays to ensure data consistency */
   for( int i = 0; i < npricingprobvars; ++i )
   {
      indsetvars[i] = NULL;
      solvals[i] = -1.0; /* To later determine whether a variable was constrained */
      consvars[i] = NULL;
   }
   for( int i = 0; i < nconss; ++i )
   {
      markedconstraints[i] = NULL;
   }

   /* For pricer */
   nsolvars = npricingprobvars;

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
         consvars = SCIPgetVarsLinear(pricingprob, constraints[i]);
         /* Check if we have an IS constraint */
         if( SCIPgetNVarsLinear(pricingprob, constraints[i]) == 2 &&
             SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),(SCIP_Real)1) ) 
         {
            /* Preprocessing: Constraint is only relevant for pricing if one of the variables has an objective value != 0 */
            if( SCIPvarGetObj(consvars[0]) != 0 || SCIPvarGetObj(consvars[1]) != 0 )
            {
               if( SCIPvarGetObj(consvars[0]) != 0 )
               {
                  nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[0]);
                  if( nodeindex0 == indexcount )
                  {
                     ++indexcount;
                  }
               }

               if( SCIPvarGetObj(consvars[1]) != 0 )
               {
                  nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[1]);
                  if( nodeindex1 == indexcount )
                  {
                     ++indexcount;
                  }
               }
               if( SCIPvarGetObj(consvars[0]) != 0 && SCIPvarGetObj(consvars[1]) != 0 )
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
            if( consvals[0] > 0 && consvals[1] > 0 && SCIPisLT(pricingprob, SCIPgetRhsLinear(pricingprob,constraints[i]),consvals[0] + consvals[1]) )
            {
               /* As before, the constraint is only regarded if it is relevant for pricing */
               if( SCIPvarGetObj(consvars[0]) != 0 )
               {
                  nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[0]);
                  if( nodeindex0 == indexcount )
                  {
                     ++indexcount;
                  }
               }

               if( SCIPvarGetObj(consvars[1]) != 0 )
               {
                  nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[1]);
                  if( nodeindex1 == indexcount )
                  {
                     ++indexcount;
                  }
               }
               if( SCIPvarGetObj(consvars[0]) != 0 && SCIPvarGetObj(consvars[1]) != 0 )
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
            consvals = SCIPgetValsLinear(pricingprob, constraints[i]);
            consvars = SCIPgetVarsLinear(pricingprob, constraints[i]);

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
                     nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[j]);
                     if( nodeindex0 == indexcount )
                     {
                        ++indexcount;
                     }
                     /* Determine nodeindex1 */
                     for( int l = j + 1; l < nvars; ++l )
                     {
                        nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[l]);
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
               /* Special case: The coupling constraint is purely decorative (coefficient >= #vars)*/
               if( abs(consvals[coefindex]) + 1 >= nvars )
               {
                  solvals[SCIPvarGetProbindex(consvars[coefindex])] = 1.0;
               }
               /* Special case: The coefficient is -1, we treat the case like a clique constraint. */
               else if( abs(consvals[coefindex]) == 1 )
               {
                  /* See if the objective coefficient of the coupling variable is != 0 */
                  if( SCIPvarGetObj(consvars[coefindex]) != 0 )
                  {
                     /* The coupling variable can always be set to 1 */
                     solvals[SCIPvarGetProbindex(consvars[coefindex])] = 1.0;
                  }
                  /* Delete the edges between all the variables of the constraint that are not the coupling variable.
                     This way, at most one can be part of the maximum clique */
                  for( int j = 0; j < nvars; ++j )
                  {
                     /* We are only interested in vars potentially relevant for pricing (!= 0) */
                     if( j != coefindex && SCIPvarGetObj(consvars[j]) != 0 )
                     {
                        /* Determine nodeindex0 */
                        nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[j]);
                        if( nodeindex0 == indexcount )
                        {
                           ++indexcount;
                        }
                        /* Determine nodeindex1 */
                        for( int l = j + 1; l < nvars; ++l )
                        {
                           if( l != coefindex )
                           {
                              nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[l]);
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
                  SCIPfreeBufferArray(pricingprob,&consvars);
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
               SCIPfreeBufferArray(pricingprob,&consvars);
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
         consvars[0] = SCIPgetVarVarbound(pricingprob,constraints[i]);
         consvars[1] = SCIPgetVbdvarVarbound(pricingprob,constraints[i]);
         /* Check value of rhs to be 0 and c to be <= -1 */ 
         if( SCIPisEQ(pricingprob,SCIPgetRhsVarbound(pricingprob,constraints[i]),0) )
         {
            if( SCIPisLT(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,constraints[i]),-1) || SCIPisEQ(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,constraints[i]),-1) ) 
            {
               /* if x may be relevant, add both x and y to graph */
               if( SCIPvarGetObj(consvars[0]) != 0 )
               {
                  if( SCIPvarGetObj(consvars[0]) != 0 )
                  {
                     nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[0]);
                     if( nodeindex0 == indexcount )
                     {
                        ++indexcount;
                     }
                  }

                  if( SCIPvarGetObj(consvars[1]) != 0 )
                  {
                     nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[1]);
                     if( nodeindex1 == indexcount )
                     {
                        ++indexcount;
                     }
                  }
                  /* It may be the case, that both the constraints x - y <= 0 and x + y <= 1 are part of the problem */
                  /* Although rare, we have to later ensure that we do not set x to 1 while y is set to 0 */
                  markedconstraints[markedcount] = constraints[i];
                  ++markedcount;
               }
               /* if only y may be relevant, add only y to the graph */
               else if( SCIPvarGetObj(consvars[1]) != 0 )
               {
                  if( SCIPvarGetObj(consvars[1]) != 0 )
                  {
                     nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[1]);
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
               SCIPfreeBufferArray(pricingprob,&consvars);
               graph_free(g);

               *result = SCIP_STATUS_UNKNOWN;
               return SCIP_OKAY;
            }
         }
         /* Rhs of varbound unequal to 0
          * It may still be the case that we have an IS constraint with a non-linear handler
          * We treat this case like a regular IS constraint 
          */
         else if( SCIPisEQ(pricingprob,SCIPgetRhsVarbound(pricingprob,constraints[i]),1) && SCIPisEQ(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,constraints[i]),1) )
         {
            /* Preprocessing: Constraint is only relevant for pricing if one of the variables has an objective value != 0 */
            if( SCIPvarGetObj(consvars[0]) != 0 || SCIPvarGetObj(consvars[1]) != 0 )
            {
               if( SCIPvarGetObj(consvars[0]) != 0 )
               {
                  nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[0]);
                  if( nodeindex0 == indexcount )
                  {
                     ++indexcount;
                  }
               }

               if( SCIPvarGetObj(consvars[1]) != 0 )
               {
                  nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[1]);
                  if( nodeindex1 == indexcount )
                  {
                     ++indexcount;
                  }
               }
               if( SCIPvarGetObj(consvars[0]) != 0 && SCIPvarGetObj(consvars[1]) != 0 )
               {
                  if( GRAPH_IS_EDGE(g,nodeindex0,nodeindex1) )
                  {
                     GRAPH_DEL_EDGE(g,nodeindex0,nodeindex1);
                  }
               }
            }            
         }
         /* Lastly, the constraint may be of the form c + 1 > rhs and c < rhs, i.e. a non-standard IS-constraint */
         else if( SCIPisLT(pricingprob,SCIPgetRhsVarbound(pricingprob,constraints[i]),SCIPgetVbdcoefVarbound(pricingprob,constraints[i]) + 1) 
               && SCIPisLT(pricingprob,SCIPgetVbdcoefVarbound(pricingprob,constraints[i]), SCIPgetRhsVarbound(pricingprob,constraints[i])) )
         {
            /* Preprocessing: Constraint is only relevant for pricing if one of the variables has an objective value != 0 */
            if( SCIPvarGetObj(consvars[0]) != 0 || SCIPvarGetObj(consvars[1]) != 0 )
            {
               if( SCIPvarGetObj(consvars[0]) != 0 )
               {
                  nodeindex0 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[0]);
                  if( nodeindex0 == indexcount )
                  {
                     ++indexcount;
                  }
               }

               if( SCIPvarGetObj(consvars[1]) != 0 )
               {
                  nodeindex1 = indset_addNodeToGraph(indexcount, npricingprobvars, scalingfactor, indsetvars, g, consvars[1]);
                  if( nodeindex1 == indexcount )
                  {
                     ++indexcount;
                  }
               }
               if( SCIPvarGetObj(consvars[0]) != 0 && SCIPvarGetObj(consvars[1]) != 0 )
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
            SCIPdebugMessage("Exit: Rhs of Varbound wrong, Iteration: %d, Rhs: %g, Coeff:%g.\n",i,SCIPgetRhsVarbound(pricingprob,constraints[i]),SCIPgetVbdcoefVarbound(pricingprob,constraints[i]));
            SCIPfreeBufferArray(pricingprob,&markedconstraints);
            SCIPfreeBufferArray(pricingprob,&indsetvars);
            SCIPfreeBufferArray(pricingprob,&solvals);
            SCIPfreeBufferArray(pricingprob,&consvars);
            graph_free(g);

            *result = SCIP_STATUS_UNKNOWN;
            return SCIP_OKAY;            
         }
      }
      else
      {
         /* Constraint handler neither linear nor varbound */
         SCIPdebugMessage("Exit: Nonlinear constraint handler, Iteration: %d.\n",i);
         SCIPdebugMessage("Constraint handler: %s\n", SCIPconshdlrGetName(conshdlr));
         SCIPfreeBufferArray(pricingprob,&markedconstraints);
         SCIPfreeBufferArray(pricingprob,&indsetvars);
         SCIPfreeBufferArray(pricingprob,&solvals);
         SCIPfreeBufferArray(pricingprob,&consvars);
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
      SCIPfreeBufferArray(pricingprob,&markedconstraints);
      SCIPfreeBufferArray(pricingprob,&indsetvars);
      SCIPfreeBufferArray(pricingprob,&solvals);
      SCIPfreeBufferArray(pricingprob,&consvars);
      graph_free(g);

      *result = SCIP_STATUS_UNKNOWN;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Graph size: %d.\n", indexcount);
   ASSERT( indexcount < npricingprobvars );

   /* indexcount now holds the actual number of unique IS variables, thus we truncate */
   if( indexcount > 0 )
   {
      graph_resize(g,indexcount);
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

   /* Set all members of the maximum clique with objective coefficient != 0 to 1 */
   for( int i = 0; i < indexcount; ++i )
   {
      if( SET_CONTAINS(clique,i) && SCIPvarGetObj(indsetvars[i]) != 0 )
      {
         solvals[SCIPvarGetProbindex(indsetvars[i])] = 1.0;
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

   /* Handle the case of varbound-IS combination */
   for( int i = 0; i < markedcount; ++i )
   {
      nodeindex0 = indset_getNodeIndex(SCIPgetVarVarbound(pricingprob,markedconstraints[i]),indsetvars,npricingprobvars);
      nodeindex1 = indset_getNodeIndex(SCIPgetVbdvarVarbound(pricingprob,markedconstraints[i]),indsetvars,npricingprobvars);   
      if( !GRAPH_IS_EDGE(g,nodeindex0,nodeindex1) )   
      {
         /* It is sufficient to ensure that only the first variable is set to 0, coefficient of varbound is negative */
         solvals[nodeindex0] = 0.0;
      }   
   }

   /* Create a column corresponding to our clique result */
   SCIP_CALL( GCGcreateGcgCol(pricingprob, &cols[0], probnr, pricingprobvars, solvals, nsolvars, FALSE, SCIPinfinity(pricingprob)) );
   *ncols = 1;
   *result = SCIP_STATUS_OPTIMAL;

   SCIPfreeBufferArray(pricingprob,&markedconstraints);
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
