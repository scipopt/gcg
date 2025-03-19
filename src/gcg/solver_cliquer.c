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

/**@file   solver_cliquer.c
 * @brief  heuristic solver for pricing problems that solves independent set problems with cliquer
 * @author Henri Lotze
 * @author Christian Puchert
 * @author Johannes Ehls
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */

#include <assert.h>

#include "gcg/solver_cliquer.h"
#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"
#include "gcg/pub_solver.h"
#include "gcg/pricer_gcg.h"
#include "gcg/relax_gcg.h"
#include "gcg/pub_gcgcol.h"
#include "gcg/pub_gcgvar.h"
#include "gcg/gcg.h"
#include "gcg/sepa_master.h"

#include "cliquer/cliquer.h"


#define SOLVER_NAME          "cliquer"
#define SOLVER_DESC          "heuristic solver for pricing problems that solves independent set problems with cliquer"
#define SOLVER_PRIORITY      150

#define SOLVER_HEURENABLED   TRUE            /**< indicates whether the solver should be enabled */
#define SOLVER_EXACTENABLED  FALSE           /**< indicates whether the solver should be enabled */

#define DEFAULT_NODELIMIT    200
#define DEFAULT_DENSITY      1.0
#define DEFAULT_DENSITYSTART 75
#define DEFAULT_USELINCUTOFF TRUE
#define DEFAULT_SLOPE        -1980
#define DEFAULT_INTERCEPT    2000
#define DEFAULT_OBJCOEFDISTR 0
#define DEFAULT_USEMULTIPL   FALSE
#define DEFAULT_CLIQUECONSTHRESH 0.5

/*
 * Data structures.
 */

/* Constraint type (combination of handler type and constraint form) to use in this solver. */
enum cliquerConsType
{
   LINEAR_IS,
   LINEAR_IS_LIKE,
   LINEAR_CLIQUE,
   LINEAR_COUPLING_DECORATIVE,
   LINEAR_COUPLING_CLIQUE,
   VARBND_SAME,
   VARBND_STD,
   VARBND_IS
};
typedef enum cliquerConsType CLIQUER_CONSTYPE;

struct GCG_SolverData
{
   SCIP_Real             density;            /**< graph density threshold above which to use solver */
   int                   densitystart;       /**< graph node threshold above which to apply density and linear cutoff */
   int                   nodelimit;          /**< graph node threshold below which to use solver */
   SCIP_Real             cliqueconsthresh;   /**< clique constraint percentage threshold below which to use solver */
   SCIP_Bool*            isnotapplicable;    /**< array tracking if solver is not applicable in root node (& no cuts) */
   int                   objcoefdistr;       /**< param deciding strategy for distributing obj coefs of coupling vars */
   SCIP_Bool             usemultiplicity;    /**< boolean for activating usage of var multiplicities for weighting */
   SCIP_Bool             uselincutoff;       /**< boolean for activating linear cutoff */
   SCIP_Real             lincutoffslope;     /**< slope for linear cutoff */
   SCIP_Real             lincutoffintercept; /**< intercept for linear cutoff */
};


/*
 * Local methods
 */

/** Returns whether the given var is linked in some way with other variables */
static INLINE
SCIP_Bool isVarLinked(
   SCIP_VAR** linkedvars,                         /**< Array of variables that are linked by eq-constraints */
   int        nlinkedvars,                        /**< Index of linkedvars array */
   SCIP_VAR*  var                                 /**< Variable whose membership in the linkedvars array is to be checked */
   )
{
   SCIP_Bool islinked;
   int       i;

   islinked = FALSE;
   for( i = 0; i < nlinkedvars; ++i )
   {
      if( linkedvars[i] == var )
         islinked = TRUE;
   }
   return islinked;
}

/** Returns whether 2 variables are linked, either simply or in a transitive way in respect to a given linkmatrix matrix.
  * Use of the wrapper function areVarsLinked(..) is recommended
  */
static
SCIP_Bool areVarsLinkedRec(
   int**          linkmatrix,                     /**< Matrix indicating which variables are linked by a node */
   int            vindex1,                        /**< Problem index of the first variable in the pair that is to be checked */
   int            vindex2,                        /**< Problem index of the second variable in the pair that is to be checked */
   int*           vartrace,                       /**< Array to keep track of which nodes have already been visited during recursion */
   int*           traceindex,                     /**< Pointer to index to keep track of the number of visited nodes during recursion */
   SCIP_VAR**     linkedvars,                     /**< Array of variables that are linked by eq-constraints */
   int            nlinkedvars                     /**< Index of linkedvars array */
   )
{
   SCIP_Bool varintrace;
   int       nextvarindex;
   int       i;
   int       j;

   /* Simple, direct link? (Matrix is symmetric) */
   if( linkmatrix[vindex1][vindex2] )
   {
      return TRUE;
   }

   /* More complex link by transitivity? */
   /* Mark current node visited by adding it to the trace. */
   vartrace[*traceindex] = vindex1;
   (*traceindex)++;
   for( i = 0; i < nlinkedvars; ++i )
   {
      nextvarindex = SCIPvarGetProbindex(linkedvars[i]);
      if( linkmatrix[vindex1][nextvarindex] )
      {
         /* To ensure termination, we have to keep track of the visited vars */
         varintrace = FALSE;
         for( j = 0; j < *traceindex; ++j )
         {
            if( vartrace[j] == nextvarindex )
            {
               varintrace = TRUE;
               break;
            }
         }
         if( !varintrace && areVarsLinkedRec(linkmatrix, nextvarindex, vindex2, vartrace, traceindex, linkedvars, nlinkedvars) )
         {
            return TRUE;
         }
      }
   }
   return FALSE;
}

/** Wrapper function for areVarsLinkedRec, mallocs and cleans up the necessary memory and passes through the result */
static
SCIP_Bool areVarsLinked(
   SCIP*          scip,                           /**< The problem instance */
   int**          linkmatrix,                     /**< Matrix indicating which variables are linked by a node */
   SCIP_VAR*      var1,                           /**< The first variable in the pair that is to be checked */
   SCIP_VAR*      var2,                           /**< The second variable in the pair that is to be checked */
   SCIP_VAR**     linkedvars,                     /**< Array of variables that are linked by eq-constraints */
   int            nlinkedvars                     /**< Index of linkedvars array */
   )
{
   int*      vartrace;
   int       traceindex;
   int       i;
   int       vindex1;
   int       vindex2;
   SCIP_Bool varslinked;

   vindex1 = SCIPvarGetProbindex(var1);
   vindex2 = SCIPvarGetProbindex(var2);

   /* We can save effort if a direct link is present */
   if( linkmatrix[vindex1][vindex2] )
   {
      return TRUE;
   }

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vartrace, nlinkedvars) );

   traceindex = 0;
   for( i = 0; i < nlinkedvars; ++i )
   {
      vartrace[i] = -1;
   }

   varslinked = areVarsLinkedRec(linkmatrix, vindex1, vindex2, vartrace, &traceindex, linkedvars, nlinkedvars);

   SCIPfreeBufferArray(scip, &vartrace);

   return varslinked;
}

/** Update transitivity in the linkmatrix matrix between 2 variables that are to be linked and all linked variables */
static
void updateVarLinks(
   SCIP*      scip,                                /**< The Problem instance */
   int**      linkmatrix,                          /**< Matrix indicating which variables are linked by a node */
   SCIP_VAR*  var1,                                /**< The first variable in the pair that is to be checked */
   SCIP_VAR*  var2,                                /**< The second variable in the pair that is to be checked */
   SCIP_VAR** linkedvars,                          /**< Array of variables that are linked by eq-constraints */
   int*       nlinkedvars                          /**< Index of linkedvars array */
   )
{
   int        varindex1, varindex2;
   int        i;
   SCIP_Bool  newvar1;
   SCIP_Bool  newvar2;

   newvar1 = TRUE;
   newvar2 = TRUE;

   /* Check if the variables are part of a link already, add them elsewise to the linkedvars array */
   for( i = 0; i < *nlinkedvars; ++i )
   {
      if( linkedvars[i] == var1 )
      {
         newvar1 = FALSE;
      }
      else if( linkedvars[i] == var2 )
      {
         newvar2 = FALSE;
      }
   }
   if( newvar1 )
   {
      linkedvars[*nlinkedvars] = var1;
      ++(*nlinkedvars);
   }
   if( newvar2 )
   {
      linkedvars[*nlinkedvars] = var2;
      ++(*nlinkedvars);
   }

   varindex1 = SCIPvarGetProbindex(var1);
   varindex2 = SCIPvarGetProbindex(var2);

   /* Variables may have not been simply linked before */
   linkmatrix[varindex1][varindex2] = 1;
   linkmatrix[varindex2][varindex1] = 1;

   /* The following loop is not necessary, as the equality graph itself is enough.
    * One might want to check performance implications of establishing cliques from all connected components.
    * This might increase speed of method areVarsLinked().
    * Or the converse: Delete the following loop to save time in this function but potentially increase time
    * consumption of areVarsLinked() which would need to traverse more nodes (increased transitivity).
    * [Test this!]
    */
   for( i = 0; i < *nlinkedvars; ++i )
   {
      /* It is sufficient to check the links between var1 and all other vars, since var1 and var2 are linked */
      if( varindex1 != SCIPvarGetProbindex(linkedvars[i]) )
      {
         if( areVarsLinked(scip, linkmatrix, var1, linkedvars[i], linkedvars, *nlinkedvars) )
         {
            /* Add links to both var1 and var2 */
            linkmatrix[varindex1][SCIPvarGetProbindex(linkedvars[i])] = 1;
            linkmatrix[SCIPvarGetProbindex(linkedvars[i])][varindex1] = 1;
            linkmatrix[varindex2][SCIPvarGetProbindex(linkedvars[i])] = 1;
            linkmatrix[SCIPvarGetProbindex(linkedvars[i])][varindex2] = 1;
         }
      }
   }
}

/** Get the node index of a given variable in a given array (for graph), else return -1 */
static INLINE
int getNodeIndex(
   SCIP_VAR*      var,                            /**< Variable for which the node index is to be determined */
   SCIP_VAR**     vararray,                     /**< Array of variables that are mapped to a node of the graph */
   int            indexcount                      /**< Number of variables that are mapped in the graph */
   )
{
   int        i;

   for( i = 0; i < indexcount; ++i )
   {
      if( var == vararray[i] )
         return i;
   }
   return -1;
}

/** Returns the node index of a given variable in the bijection or that of a linked variable, if any. */
static 
int getLinkedNodeIndex(
   SCIP*          scip,                           /**< The problem instance */
   SCIP_VAR*      var,                            /**< Variable for which the node index is to be determined */
   SCIP_VAR**     indsetvars,                     /**< Array of variables that are mapped to a node of the graph */
   int            indexcount,                     /**< Number of variables that are mapped in the graph */
   int**          linkmatrix,                     /**< Matrix indicating which variables are linked by a node */
   SCIP_VAR**     linkedvars,                     /**< Array of variables that are linked by eq-constraints */
   int            nlinkedvars                     /**< Index of linkedvars array */
   )
{
   int        nodeindex;
   int        i;

   nodeindex = getNodeIndex(var, indsetvars, indexcount);
   if( nodeindex == -1 && isVarLinked(linkedvars, nlinkedvars, var) )
   {
      for( i = 0; i < nlinkedvars; ++i )
      {
         if( linkedvars[i] != var )
         {
            if( areVarsLinked(scip, linkmatrix, var, linkedvars[i], linkedvars, nlinkedvars) )
            {
               nodeindex = getNodeIndex(linkedvars[i], indsetvars, indexcount);
               if( nodeindex != -1 )
               {
                  return nodeindex;
               }
            }
         }
      }
   }
   else
   {
      return nodeindex;
   }
   return -1;
}

/** Computes the number of reachable nodes from a given variable */
static INLINE
int countReachableVars(
   SCIP*          scip,         /**< The problem instance */
   int**          linkmatrix,   /**< Matrix indicating which variables are linked */
   SCIP_VAR*      var,          /**< The starting variable */
   SCIP_VAR**     linkedvars,   /**< Array of variables that are linked by eq-constraints */
   int            nlinkedvars   /**< Number of linked variables */
)
{
   int*       stack;           /* Stack for DFS */
   SCIP_Bool* varvisited;      /* Array to track visited nodes */
   int        stacksize;        /* Current size of stack */
   int        count;            /* Count of reachable nodes */
   int        actind;
   int        popedvarind;
   int        i;

   count = 0;
   stacksize = 0;

   actind = SCIPvarGetProbindex(var);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &stack, nlinkedvars) );
   SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &varvisited, nlinkedvars) );

   /* Start DFS from the given variable */
   stack[stacksize++] = actind;
   varvisited[getNodeIndex(var, linkedvars, nlinkedvars)] = TRUE;

   while( stacksize > 0 )
   {
      popedvarind = stack[--stacksize];         /* pop from stack */
      count++;

      /* Explore all linked variables */
      for( i = 0; i < nlinkedvars; i++ )
      {
         actind = SCIPvarGetProbindex(linkedvars[i]);
         if( linkmatrix[popedvarind][actind] && !varvisited[i] )
         {
            stack[stacksize++] = actind;
            varvisited[i] = TRUE;
         }
      }
   }

   SCIPfreeBufferArray(scip, &varvisited);
   SCIPfreeBufferArray(scip, &stack);

   return count;
}

/** Returns a representative variable for a given variable in the linked variable-bijection, if any, else NULL.*/
static INLINE
SCIP_VAR* getLinkedNodeVar(
   SCIP*          scip,                           /**< The problem instance */
   SCIP_VAR*      var,                            /**< Variable for which the node index is to be determined */
   int**          linkmatrix,                     /**< Matrix indicating which variables are linked by a node */
   SCIP_VAR**     linkedvars,                     /**< Array of variables that are linked by eq-constraints */
   int            nlinkedvars                     /**< Index of linkedvars array */
)
{
   int        i;

   for( i = 0; i < nlinkedvars; ++i )
   {
      if( linkedvars[i] == var || areVarsLinked(scip, linkmatrix, var, linkedvars[i], linkedvars, nlinkedvars) )
         return linkedvars[i];
   }
   return NULL;
}

/** Returns the index of the representative variable for a given variable in the linked variable-bijection, if any.*/
static INLINE
int getNodeIndexCouplDigraph(
   SCIP*          scip,                           /**< The problem instance */
   SCIP_VAR*      var,                            /**< Variable for which the node index is to be determined */
   int**          linkmatrix,                     /**< Matrix indicating which variables are linked by a node */
   SCIP_VAR**     linkedvars,                     /**< Array of variables that are linked by eq-constraints */
   int            nlinkedvars,                    /**< Index of linkedvars array */
   SCIP_VAR**     varsincouplings,                /**< Array of variables that are involved in couplings */
   int            nvarsincouplings                /**< Index of varsincouplings array */
   )
{
   if( isVarLinked(linkedvars, nlinkedvars, var) )
      var = getLinkedNodeVar(scip, var, linkmatrix, linkedvars, nlinkedvars);

   return (var == NULL) ? -1 : getNodeIndex(var, varsincouplings, nvarsincouplings);
}

/** Aggregates objective coefficients for linked variables in the aggrobjcoef array. */
static
void aggregateObjCoef(
   SCIP*          scip,                           /**< The problem instance */
   int**          linkmatrix,                     /**< Matrix indicating which variables are linked by a node */
   SCIP_VAR**     linkedvars,                     /**< Array of variables that are linked by eq-constraints */
   int            nlinkedvars,                    /**< Index of linkedvars array */
   SCIP_Real*     aggrobjcoef                     /**< Array of aggregated objective coefficients. */
   )
{
   int        i;
   int        j;
   SCIP_Real  aggr;
   int*       valisset;            /* 0 : unset; 1 : set; 2 : to be set. */

   SCIP_CALL_ABORT( SCIPallocClearBufferArray(scip, &valisset, nlinkedvars) );

   for( i = 0; i < nlinkedvars; i++)
   {
      if( valisset[i] == 1 )
         continue;

      aggr = SCIPvarGetObj(linkedvars[i]);
      valisset[i] = 2;

      for( j = i + 1; j < nlinkedvars; j++ )
      {
         if( !valisset[j] && areVarsLinked(scip, linkmatrix, linkedvars[i], linkedvars[j], linkedvars, nlinkedvars) )
         {
            aggr += SCIPvarGetObj(linkedvars[j]);
            valisset[j] = 2;
         }
      }

      for( j = i; j < nlinkedvars; j++ )
      {
         if( valisset[j] == 2 )
         {
            aggrobjcoef[SCIPvarGetProbindex(linkedvars[j])] = aggr;
            valisset[j] = 1;
         }
      }
   }

   SCIPfreeBufferArray(scip, &valisset);
}

/** Returns the aggregated objective coefficient. */
static INLINE
SCIP_Real getAggrObjCoef(
   SCIP_VAR*      var,                            /**< The variable for which the aggr. obj. coef. should be returned */
   int            nlinkedvars,                    /**< Index of linkedvars array */
   int            nvarsincouplings,               /**< Index of varsincouplings array */
   SCIP_Real*     aggrobjcoef                     /**< Array of aggregated objective coefficients. */
   )
{
   if( nlinkedvars > 0 || nvarsincouplings > 0 )
      return aggrobjcoef[SCIPvarGetProbindex(var)];
   return SCIPvarGetObj(var);
}

/** Set (/update) aggregated objective coefficients savely. */
static INLINE
void setAggrObjCoef(
   SCIP*          scip,                           /**< The problem instance */
   SCIP_VAR*      var,                            /**< The variable for which the aggr. obj. coef. should be updated */
   SCIP_Real      value,                          /**< Value to set for given var (and linked vars) */
   SCIP_VAR**     linkedvars,                     /**< Array of variables that are linked by eq-constraints */
   int            nlinkedvars,                    /**< Index of linkedvars array */
   int**          linkmatrix,                     /**< Matrix indicating which variables are linked by a node */
   SCIP_Real*     aggrobjcoef                     /**< Array of aggregated objective coefficients. */
)
{
   int            i;

   /* If variable is linked, we need to update all linked variables' objective coefficients too */
   if( isVarLinked(linkedvars, nlinkedvars, var) )
   {
      for( i = 0; i < nlinkedvars; i++ )
      {
         if( areVarsLinked(scip, linkmatrix, var, linkedvars[i], linkedvars, nlinkedvars) )
            aggrobjcoef[SCIPvarGetProbindex(linkedvars[i])] = value;
      }
   }
   else
   {
      aggrobjcoef[SCIPvarGetProbindex(var)] = value;
   }
}

/** Add a variable to the bijection graph g and indsetvars array. Returns the index of the corresponding node in the graph. */
static
int addVarToGraph(
   SCIP*          scip,                           /**< The problem instance */
   graph_t*       g,                              /**< Graph into which to insert the new variable as a node */
   SCIP_VAR*      consvar,                        /**< The variable that is to be assigned a node in the graph */
   int*           indexcount,                     /**< Pointer to Index of the next unassigned node in the graph */
   SCIP_Real      scalingfactor,                  /**< Factor for scaling the weight of newly mapped nodes */
   SCIP_VAR**     indsetvars,                     /**< Array that keeps track of variables that are part of the graph */
   int**          linkmatrix,                     /**< Matrix indicating which variables are linked by a node */
   SCIP_VAR**     linkedvars,                     /**< Array of variables that are linked by eq-constraints */
   int            nlinkedvars,                    /**< Index of linkedvars array */
   int            nvarsincouplings,               /**< Index of varsincouplings array */
   SCIP_Real*     aggrobjcoef                     /**< Array of aggregated objective coefficients. */
   )
{
   int nodeindex;
   SCIP_Real aggrObj;
   if( isVarLinked(linkedvars, nlinkedvars, consvar) )
   {
      nodeindex = getLinkedNodeIndex(scip, consvar, indsetvars, *indexcount, linkmatrix, linkedvars, nlinkedvars);
   }
   else
   {
      nodeindex = getNodeIndex(consvar, indsetvars, *indexcount);
   }
   if( nodeindex == -1 )
   {
      /* Var not yet part of graph, add it with its corresponding weight */
      indsetvars[*indexcount] = consvar;
      aggrObj = getAggrObjCoef(indsetvars[*indexcount], nlinkedvars, nvarsincouplings, aggrobjcoef);
      if( SCIPisLT(scip, aggrObj, 0.0) )
      {
         g->weights[*indexcount] = 1 + abs((int) (scalingfactor * aggrObj));
      }
      else
      {
         g->weights[*indexcount] = 1;
      }
      nodeindex = *indexcount;
      ++(*indexcount);
   }
   return nodeindex;
}

/** Set the solvals of a variable and of all its linked variables, if any */
static
void setLinkedSolvals(
   SCIP*       scip,                              /**< The problem instance */
   SCIP_Real*  solvals,                           /**< Array holding the current solution values of all variables of the problem */
   int**       linkmatrix,                        /**< Matrix indicating which variables are linked by a node */
   SCIP_VAR**  linkedvars,                        /**< Array of variables that are linked by eq-constraints */
   int         nlinkedvars,                       /**< Index of linkedvars array */
   SCIP_VAR*   var,                               /**< Var that may be linked that itself should be set to the value val */
   SCIP_Real   val                                /**< Value to which all linked vars are supposed to be set to */
   )
{
   int i;

   solvals[SCIPvarGetProbindex(var)] = val;

   for( i = 0; i < nlinkedvars; ++i )
   {
      if( var != linkedvars[i] )
      {
         if( areVarsLinked(scip, linkmatrix, var, linkedvars[i], linkedvars, nlinkedvars) )
         {
            solvals[SCIPvarGetProbindex(linkedvars[i])] = val;
            assert(SCIPisGE(scip, val, SCIPvarGetLbLocal(linkedvars[i])) &&
                   SCIPisLE(scip, val, SCIPvarGetUbLocal(linkedvars[i])));
         }
      }
   }
}

/** Check if the objective coefficients of the variables are already Integral */
static 
SCIP_Bool areObjectivesIntegral(
   SCIP*          scip,                           /**< The problem instance */
   SCIP_VAR**     linkedvars,                     /**< Array of variables that are linked by eq-constraints */
   int            nlinkedvars,                    /**< Index of linkedvars array */
   SCIP_VAR**     varsincouplings,                /**< Array of variables that are involved in coupling */
   int            nvarsincouplings,               /**< Index of varsincouplings array */
   SCIP_Real*     aggrobjcoef                     /**< Array of aggregated objective coefficients. */
   )
{
   SCIP_Real  objval;
   SCIP_Real  nvars;
   SCIP_VAR** vars;
   int        i;

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   for( i = 0; i < nvars; ++i )
   {
      objval = getAggrObjCoef(vars[i], nlinkedvars, nvarsincouplings, aggrobjcoef);
      if( !SCIPisZero(scip, objval-((int)objval)) )
      {
         return FALSE;
      }
   }
   return TRUE;
}

/** Scale the objective coefficients of the variables maximally s.t. they become integral and the sum of values does not exceed INT_MAX */
static 
SCIP_Real scaleRelativeToMax(
   SCIP*          scip,                           /**< The problem instance */
   SCIP_VAR**     linkedvars,                     /**< Array of variables that are linked by eq-constraints */
   int            nlinkedvars,                    /**< Index of linkedvars array */
   SCIP_VAR**     varsincouplings,                /**< Array of variables that are involved in coupling */
   int            nvarsincouplings,               /**< Index of varsincouplings array */
   SCIP_Real*     aggrobjcoef                     /**< Array of aggregated objective coefficients. */
   )
{
   SCIP_Real  scalingfactor;
   SCIP_Real  varval;
   SCIP_Real  biggestobj;
   SCIP_Real  nvars;
   SCIP_VAR** vars;
   int        i;

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   scalingfactor = (INT_MAX / nvars) - nvars;

   /* Check for the biggest objective value to safely adjust the scalingfactor */
   biggestobj = 0.0;
   for( i = 0; i < nvars; ++i )
   {
      varval = getAggrObjCoef(vars[i], nlinkedvars, nvarsincouplings, aggrobjcoef);
      if( SCIPisLT(scip, varval, biggestobj) )
      {
         biggestobj = varval;
      }
   }
   if( SCIPisLT(scip, biggestobj, -1.0) )
   {
      /* Ensure that INT_MAX is never reached by the sum of all scaled weights */
      scalingfactor = fabs(scalingfactor / biggestobj);
   }
   return scalingfactor;
}

/** Set isnotapplicable to true for given problem number if solver is applied at root node and no cuts are applied. */
static
void setProblemNotApplicable(
   SCIP*                 scip,               /**< master problem SCIP data structure */
   int                   probnr,             /**< pricing problem number */
   SCIP_Bool*            isnotapplicable     /**< array storing if solver is applicable to problems */
   )
{
   if( SCIPgetFocusDepth(scip) == 0 && GCGsepaGetNCuts(scip) == 0 )
      isnotapplicable[probnr] = TRUE;
}

/** Returns index in adjacency matrix of coupling digraph (after inserting it if not already contained). */
static INLINE
int assureInCouplingGraph(
   SCIP*                 scip,               /**< The problem instance */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in coupling */
   int*                  nvarsincouplings,   /**< Index of varsincouplings array */
   SCIP_VAR*             var,                /**< Variable whose membership in the varsincouplings array is to be checked */
   SCIP_VAR**            linkedvars,         /**< Array of variables that are linked by eq-constraints */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   int**                 linkmatrix          /**< Matrix indicating which variables are linked by eq-constraints */
   )
{
   int       i;

   /* If var is linked, we map all linked vars to one representative digraph node and return its index. */
   if( isVarLinked(linkedvars, nlinkedvars, var) )
      var = getLinkedNodeVar(scip, var, linkmatrix, linkedvars, nlinkedvars);

   for( i = 0; i < *nvarsincouplings; ++i )
   {
      if( varsincouplings[i] == var )
         return i;
   }

   varsincouplings[*nvarsincouplings] = var;
   return (*nvarsincouplings)++;
}

/** Update coupling digraph for a given coupling(like) constraint, i.e., the coupling and constraint vars.
 *  This is done by assuring all variables have a corresponding index in the graph.
 *  Then, directed edges from coupling var to each other var involved in the constraint are inserted. */
static INLINE
void updateCouplingDiGraph(
   SCIP*                 scip,               /**< The problem instance */
   SCIP_VAR**            consvars,           /**< Array of variables in constraint for which to update the graph */
   int                   nconsvars,          /**< Index of consvars array */
   SCIP_VAR*             var,                /**< Coupling / vbd variable for which to update graph */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int*                  nvarsincouplings,   /**< Index of varsincouplings array */
   SCIP_VAR**            linkedvars,         /**< Array of variables that are linked by eq-constraints */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   int**                 linkmatrix          /**< Matrix indicating which variables are linked by eq-constraints */
   )
{
   int            i;
   int            couplingvarind;
   int            consvarind;

   couplingvarind = assureInCouplingGraph(scip, varsincouplings, nvarsincouplings, var, linkedvars, nlinkedvars,
                                          linkmatrix);

   /* Insert (increment) edges for all non-coupling variables in the constraint. */
   for( i = 0; i < nconsvars; i++ )
   {
      if( consvars[i] == var )
         continue;

      consvarind  = assureInCouplingGraph(scip, varsincouplings, nvarsincouplings, consvars[i], linkedvars, nlinkedvars,
                                          linkmatrix);
      if( !couplingmatrix[couplingvarind][consvarind] )
         couplingmatrix[couplingvarind][consvarind] = 1;
   }
}

/** Recursive function - should be called through isCouplingRelevant() or other non-recursive methods.
 *  Checks if a variable for a given index in the coupling digraph is relevant for distributing the coefficients. */
static
SCIP_Bool isCouplingRelevantRec(
   SCIP*                 scip,               /**< The problem instance */
   int                   varind,             /**< Variable index for which to check if it has a successor */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int                   nvarsincouplings,   /**< Index of varsincouplings array */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   SCIP_Real*            aggrobjcoef,        /**< Array of aggregated objective coefficients. */
   int                   maxdepth            /**< Integer limiting recursion to finite depth in case of a cycle. */
)
{
   int            i;

   if( SCIPisLT(scip, getAggrObjCoef(varsincouplings[varind], nlinkedvars, nvarsincouplings, aggrobjcoef), 0) )
      return TRUE;

   if( maxdepth <= 0 )
      return FALSE;

   for( i = 0; i < nvarsincouplings; i++ )
   {
      if( couplingmatrix[varind][i] > 0 &&
          isCouplingRelevantRec(scip, i, couplingmatrix, varsincouplings, nvarsincouplings, nlinkedvars, aggrobjcoef,
                                maxdepth - 1))
      {
         return TRUE;
      }
   }
   return FALSE;
}

/** Checks if a variable for a given index in the coupling digraph is relevant for distributing the coefficients.
 *  A variable is said to be relevant iff
 *      - it has a negative objective coefficient or
 *      - in the digraph a node corresponding to a variable with negative obj. coef. is reachable in the digraph. */
static INLINE
SCIP_Bool isCouplingRelevant(
   SCIP*                 scip,               /**< The problem instance */
   int                   varind,             /**< Variable index for which to check if it has a successor */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int                   nvarsincouplings,   /**< Index of varsincouplings array */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   SCIP_Real*            aggrobjcoef         /**< Array of aggregated objective coefficients. */
   )
{
   return isCouplingRelevantRec(scip, varind, couplingmatrix, varsincouplings, nvarsincouplings, nlinkedvars,
                                aggrobjcoef, nvarsincouplings);
}

/** Same as isCouplingRelevant(), but takes a SCIP_VAR pointer instead of the digraph index. */
static INLINE
SCIP_Bool isCouplingRelevantVar(
   SCIP*                 scip,               /**< The problem instance */
   SCIP_VAR*             var,                /**< Variable for which to check if it has a successor */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int                   nvarsincouplings,   /**< Index of varsincouplings array */
   int**                 linkmatrix,         /**< Matrix indicating which variables are linked by eq-constraints */
   SCIP_VAR**            linkedvars,         /**< Array of variables that are linked by eq-constraints */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   SCIP_Real*            aggrobjcoef         /**< Array of aggregated objective coefficients. */
)
{
   int            varind;

   varind = getNodeIndexCouplDigraph(scip, var, linkmatrix, linkedvars, nlinkedvars, varsincouplings, nvarsincouplings);
   assert(varind > -1);

   return isCouplingRelevantRec(scip, varind, couplingmatrix, varsincouplings, nvarsincouplings, nlinkedvars,
                                aggrobjcoef, nvarsincouplings);
}

/** Checks in the coupling digraph if a var with a given index has a coupling relevant successor. */
static INLINE
SCIP_Bool hasSuccessorRel(
   SCIP*                 scip,               /**< The problem instance */
   int                   varind,             /**< Variable index for which to check if it has a successor */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int                   nvarsincouplings,   /**< Index of varsincouplings array */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   SCIP_Real*            aggrobjcoef         /**< Array of aggregated objective coefficients. */
   )
{
   int            i;

   for( i = 0; i < nvarsincouplings; i++ )
   {
      if( couplingmatrix[varind][i] > 0 &&
          isCouplingRelevantRec(scip, i, couplingmatrix, varsincouplings, nvarsincouplings, nlinkedvars, aggrobjcoef,
                                nvarsincouplings) )
      {
         return TRUE;
      }
   }
   return FALSE;
}

/** Get the number of coupling-relevant successors in the Coupling DiGraph. */
static INLINE
int getNSuccessorsRelevant(
   SCIP*                 scip,               /**< The problem instance */
   int                   varind,             /**< Variable index for which to count the successors */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int                   nvarsincouplings,   /**< Index of varsincouplings array */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   SCIP_Real*            aggrobjcoef,        /**< Array of aggregated objective coefficients. */
   int*                  varmultiplicities,  /**< Array holding number of vars represented by digraph node */
   SCIP_Bool             usemultiplicities   /**< Wether to use the given multiplicities */
)
{
   int            i;
   int            nsuccessors;

   nsuccessors = 0;
   for( i = 0; i < nvarsincouplings; i++ )
   {
      if( couplingmatrix[varind][i] > 0 &&
          isCouplingRelevantRec(scip, i, couplingmatrix, varsincouplings, nvarsincouplings, nlinkedvars, aggrobjcoef,
                                nvarsincouplings) )
      {
         nsuccessors += usemultiplicities ? varmultiplicities[i] : 1;
      }
   }
   return nsuccessors;
}

/** Compute entries for varmultiplicities array holding counts of represented variables per digraph node */
static INLINE
void initVarMultiplicities(
   SCIP*                 scip,               /**< The problem instance */
   int**                 linkmatrix,         /**< Matrix indicating which variables are linked */
   SCIP_VAR**            linkedvars,         /**< Array of variables that are linked by eq-constraints */
   int                   nlinkedvars,        /**< Number of linked variables */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int                   nvarsincouplings,   /**< Index of varsincouplings array */
   int*                  varmultiplicities   /**< Array holding number of vars represented by digraph node */
   )
{
   int            i;

   for( i = 0; i < nvarsincouplings; i++ )
   {
      if( isVarLinked(linkedvars, nlinkedvars, varsincouplings[i]) )
         varmultiplicities[i] = countReachableVars(scip, linkmatrix, varsincouplings[i], linkedvars, nlinkedvars);
      else
         varmultiplicities[i] = 1;
   }
}

/*
 * Idea - Distribution Strategy 1: Natural Coefficient Share Distribution.
 *
 * Distributes a variable's objective coefficient among its coupled successors using a structured
 * closed-form approximation, ensuring no feasible assignment overestimates the original objective function.
 *
 * **Step 1: Distribute the Objective Coefficient Among Constraints**
 * Given m constraints of the form x_{i,1} + ... + x_{i,n_i} ≤ c_i * y, each constraint i receives:
 *
 *      w_i = w * (c_i / sum(c_j for all j in constraints))
 *
 * **Step 2: Distribute w_i Among Relevant Variables**
 * Each relevant x_{i,j} variable in constraint i gets:
 *
 *      w_ij = w_i * (1 / n_rel_cons_vars)
 *
 * **Implementation:**
 * - Iterate over constraints, identify relevant variables, and accumulate their coefficient shares (`coefshares[j]`).
 * - Normalize using `denominator`:
 *
 *      aggrobjcoef[j] += actvarcoef * (coefshares[j] / denominator)
 *
 * **Special Cases:**
 * - **Varbound Standard & Clique Constraints**: Each relevant variable gets **1**.
 * - **Decorative Coupling Constraints** (x_1 + ... + x_n ≤ c * y, c ≥ n): Each variable gets **1 / n_rel_cons_vars**.
 *
 * This method ensures the total assigned weight remains ≤ w and that feasible selections of x_ij never exceed w.
 */

/** Coupling(-like) Variable Objective Distribution Strategy 1: Natural Coefficient Share Distribution. */
static INLINE
void objCoefDistrHeurNatural(
   SCIP*                 scip,               /**< The problem instance */
   int                   actvarind,          /**< Index of actual variable */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int                   nvarsincouplings,   /**< Index of varsincouplings array */
   int*                  varmultiplicities,  /**< Array holding number of vars represented by digraph node */
   SCIP_Bool             usemultipl,         /**< Wether varmultiplicities should be used for weighting */
   SCIP_VAR**            linkedvars,         /**< Array of variables that are linked by eq-constraints */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   int**                 linkmatrix,         /**< Matrix indicating which variables are linked by eq-constraints */
   SCIP_Real*            aggrobjcoef,        /**< Array of aggregated objective coefficients. */
   SCIP_CONS**           constraints,        /**< Array containing pointers to SCIP constraints of pricing problem */
   int                   nconss,             /**< Index of constraints array */
   int*                  couplingcoefindices,/**< Array for coupling coefficient of each constraint (if coupling) */
   CLIQUER_CONSTYPE*     cliquerconstypes    /**< Array holding constraint types (specific to this solver) */
   )
{
   int            i;
   int            j;
   int            denominator;
   int            nrelconsvars;
   int            nsuccessors;
   int            nlconsvars;
   int            couplingindex;
   int            mappedindex;
   SCIP_Real      frac;
   SCIP_Real      coeftodistr;
   SCIP_Real      actvarcoef;
   SCIP_Real*     coefshares;                /* array storing the fraction of the actvar coefficient they receive */
   SCIP_HASHMAP*  vartocoefsharemap;         /* hash storing relevant var index; image is index for coefshares array */
   SCIP_VAR**     lconsvars;
   SCIP_Bool      retcode;

   denominator = 0;

   nsuccessors = getNSuccessorsRelevant(scip, actvarind, couplingmatrix, varsincouplings, nvarsincouplings,
                                        nlinkedvars, aggrobjcoef, varmultiplicities, FALSE);
   SCIP_CALL_ABORT( SCIPallocClearBufferArray(scip, &coefshares, nsuccessors) );
   SCIP_CALL_ABORT( SCIPhashmapCreate(&vartocoefsharemap, SCIPblkmem(scip), nvarsincouplings) );

   /* Setup mapping to coefficient share array */
   j = 0;
   for( i = 0; i < nvarsincouplings; i++ )
   {
      if( couplingmatrix[actvarind][i] > 0 &&
          isCouplingRelevant(scip, i, couplingmatrix, varsincouplings, nvarsincouplings, nlinkedvars, aggrobjcoef) )
      {
         SCIPhashmapSetImageInt(vartocoefsharemap, (void*)(size_t)i, j);
         j++;
      }
   }

   assert(j == nsuccessors);

   /* Calculate coefficient shares */
   for( i = 0; i < nconss; i++ )
   {
      nrelconsvars = 0;
      switch( cliquerconstypes[i] )
      {
         case VARBND_STD:
            if( SCIPgetVbdvarVarbound(scip, constraints[i]) == varsincouplings[actvarind] &&
                isCouplingRelevantVar(scip, SCIPgetVarVarbound(scip, constraints[i]), couplingmatrix, varsincouplings,
                                      nvarsincouplings, linkmatrix, linkedvars, nlinkedvars, aggrobjcoef) )
            {
               couplingindex = getNodeIndexCouplDigraph(scip, SCIPgetVarVarbound(scip, constraints[i]), linkmatrix,
                                                        linkedvars, nlinkedvars, varsincouplings, nvarsincouplings);
               mappedindex = SCIPhashmapGetImageInt(vartocoefsharemap, (void*)(size_t)couplingindex);
               nrelconsvars += usemultipl ? varmultiplicities[couplingindex] : 1;
               coefshares[mappedindex] += usemultipl ? varmultiplicities[couplingindex] : 1;
            }
            break;
         case LINEAR_COUPLING_CLIQUE:
         case LINEAR_COUPLING_DECORATIVE:
            lconsvars = SCIPgetVarsLinear(scip, constraints[i]);
            if( lconsvars[couplingcoefindices[i]] == varsincouplings[actvarind] )
            {
               /* Get # of relevant variables in constraint */
               SCIPgetConsNVars(scip, constraints[i], &nlconsvars, &retcode);
               for( j = 0; j < nlconsvars; j++ )
               {
                  if( j == couplingcoefindices[i] )
                     continue;
                  if( isCouplingRelevantVar(scip, lconsvars[j], couplingmatrix, varsincouplings, nvarsincouplings,
                                            linkmatrix, linkedvars, nlinkedvars, aggrobjcoef))
                  {
                     couplingindex = getNodeIndexCouplDigraph(scip, lconsvars[j], linkmatrix, linkedvars,
                                                              nlinkedvars, varsincouplings, nvarsincouplings);
                     nrelconsvars += usemultipl ? varmultiplicities[couplingindex] : 1;
                  }
               }
               /* Add values to coefshares of variables in constraint */
               if( nrelconsvars > 0 )
               {
                  frac = cliquerconstypes[i] == LINEAR_COUPLING_CLIQUE ? 1 : (1 / nrelconsvars);
                  for( j = 0; j < nlconsvars; j++ )
                  {
                     if( j == couplingcoefindices[i] )
                        continue;
                     if( isCouplingRelevantVar(scip, lconsvars[j], couplingmatrix, varsincouplings, nvarsincouplings,
                                               linkmatrix, linkedvars, nlinkedvars, aggrobjcoef) )
                     {
                        couplingindex = getNodeIndexCouplDigraph(scip, lconsvars[j], linkmatrix, linkedvars,
                                                                 nlinkedvars, varsincouplings, nvarsincouplings);
                        mappedindex = SCIPhashmapGetImageInt(vartocoefsharemap, (void*)(size_t)couplingindex);

                        nrelconsvars += usemultipl ? varmultiplicities[couplingindex] : 1;
                        coefshares[mappedindex] += (usemultipl ? varmultiplicities[couplingindex] : 1) * frac;
                     }
                  }
               }
            }
            break;
         default:
            break;
      }
      denominator += nrelconsvars;
   }

   /* Distribute coeff. of act. var to all successor variables. */
   coeftodistr = getAggrObjCoef(varsincouplings[actvarind], nlinkedvars, nvarsincouplings, aggrobjcoef);
   for( i = 0; i < nvarsincouplings; i++ )
   {
      if( couplingmatrix[actvarind][i] > 0 &&
          isCouplingRelevant(scip, i, couplingmatrix, varsincouplings, nvarsincouplings, nlinkedvars, aggrobjcoef) )
      {
         frac = coefshares[SCIPhashmapGetImageInt(vartocoefsharemap, (void *) (size_t) i)];
         actvarcoef = getAggrObjCoef(varsincouplings[i], nlinkedvars, nvarsincouplings, aggrobjcoef);
         setAggrObjCoef(scip, varsincouplings[i], actvarcoef + coeftodistr * (frac / denominator),
                        linkedvars, nlinkedvars, linkmatrix, aggrobjcoef);
      }
   }

   SCIPhashmapFree(&vartocoefsharemap);
   SCIPfreeBufferArray(scip, &coefshares);
}

/*
 * Idea - Distribution Strategy 2: Independent set (IS)-based share distribution.
 *
 * This heuristic distributes the objective coefficient of a coupled variable among its relevant successor variables
 * by leveraging **maximum independent sets** in the coupling graph. The goal is to ensure a balanced coefficient
 * distribution while preventing overestimation of the objective function.
 *
 * **Core Idea:**
 * - Construct a **graph representation** where nodes correspond to relevant successor variables.
 * - Create edges between all pairs of successor variables initially.
 * - Remove edges based on constraints of type `LINEAR_COUPLING_CLIQUE`, ensuring that only truly independent
 *   variables remain connected.
 * - Compute a **maximum independent set (MIS)** in this reduced graph using the **cliquer** library.
 * - Distribute the objective coefficient **equally** among the variables in this MIS.
 *
 * This ensures no overestimation of the redistributed objective coefficient.
 *
 * Attention:  If the coupling graph is too large, the cliquer library might not solve the problem in acceptable time.
 *             Therefore, a hard limit of 200 nodes is implemented. Otherwise, no distribution is done.
 */

/** Coupling(-like) Variable Objective Distribution Strategy 2: Independent Set-based share distribution. */
static INLINE
void objCoefDistrHeurIS(
   SCIP*                 scip,               /**< The problem instance */
   int                   actvarind,          /**< Index of actual variable */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int                   nvarsincouplings,   /**< Index of varsincouplings array */
   int*                  varmultiplicities,  /**< Array holding number of vars represented by digraph node */
   SCIP_Bool             usemultipl,         /**< Wether varmultiplicities should be used for weighting */
   SCIP_VAR**            linkedvars,         /**< Array of variables that are linked by eq-constraints */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   int**                 linkmatrix,         /**< Matrix indicating which variables are linked by eq-constraints */
   SCIP_Real*            aggrobjcoef,        /**< Array of aggregated objective coefficients. */
   SCIP_CONS**           constraints,        /**< Array containing pointers to SCIP constraints of pricing problem */
   int                   nconss,             /**< Index of constraints array */
   int*                  couplingcoefindices,/**< Array for coupling coefficient of each constraint (if coupling) */
   CLIQUER_CONSTYPE*     cliquerconstypes    /**< Array holding constraint types (specific to this solver) */
)
{
   int            i;
   int            j;
   int            k;
   int            couplind1;
   int            couplind2;
   int            denominator;
   int            nsuccessors;
   int            nlconsvars;
   SCIP_Real      actvarcoef;
   SCIP_Real      coeftodistr;
   SCIP_HASHMAP*  vartosuccmap;       /* hash storing relevant var index; image is index for corresponding graph node */
   SCIP_VAR**     lconsvars;
   SCIP_Bool      retcode;
   set_t          clique;
   graph_t*       g;
   clique_options cl_opts;

   nsuccessors = getNSuccessorsRelevant(scip, actvarind, couplingmatrix, varsincouplings, nvarsincouplings,
                                        nlinkedvars, aggrobjcoef, varmultiplicities, FALSE);

   if( nsuccessors > 200 || nsuccessors == 0)
      return;

   g = graph_new(nsuccessors);
   for( i = 0; i < nsuccessors; ++i )
   {
      for( j = i + 1; j < nsuccessors; ++j )
      {
         GRAPH_ADD_EDGE(g, i, j);
      }
   }

   /* Setup mapping to coefficient share array */
   j = 0;
   SCIP_CALL_ABORT( SCIPhashmapCreate(&vartosuccmap, SCIPblkmem(scip), nvarsincouplings) );
   for( i = 0; i < nvarsincouplings; i++ )
   {
      if( couplingmatrix[actvarind][i] > 0 &&
          isCouplingRelevant(scip, i, couplingmatrix, varsincouplings, nvarsincouplings, nlinkedvars, aggrobjcoef) )
      {
         SCIPhashmapSetImageInt(vartosuccmap, (void*)(size_t)i, j);
         if( usemultipl )
            g->weights[j] = varmultiplicities[i];
         j++;
      }
   }
   assert(j == nsuccessors);

   /* Create IS problem to determine coefficient share */
   for( i = 0; i < nconss; i++ )
   {
      if( LINEAR_COUPLING_CLIQUE == cliquerconstypes[i] &&
          SCIPgetVarsLinear(scip, constraints[i])[couplingcoefindices[i]] == varsincouplings[actvarind] )
      {
         /* Delete edges between nodes of relevant variables in constraint */
         lconsvars = SCIPgetVarsLinear(scip, constraints[i]);
         SCIPgetConsNVars(scip, constraints[i], &nlconsvars, &retcode);
         for( j = 0; j < nlconsvars; j++ )
         {
            if( j == couplingcoefindices[i] )
               continue;
            if( isCouplingRelevantVar(scip, lconsvars[j], couplingmatrix, varsincouplings, nvarsincouplings,
                                      linkmatrix, linkedvars, nlinkedvars, aggrobjcoef) )
            {
               for( k = j + 1; k < nlconsvars; k++ )
               {
                  if( k == couplingcoefindices[i] )
                     continue;
                  if( isCouplingRelevantVar(scip, lconsvars[k], couplingmatrix, varsincouplings, nvarsincouplings,
                                            linkmatrix, linkedvars, nlinkedvars, aggrobjcoef) )
                  {
                     couplind1 = getNodeIndexCouplDigraph(scip, lconsvars[j], linkmatrix, linkedvars, nlinkedvars,
                                                          varsincouplings, nvarsincouplings);
                     couplind2 = getNodeIndexCouplDigraph(scip, lconsvars[k], linkmatrix, linkedvars, nlinkedvars,
                                                          varsincouplings, nvarsincouplings);
                     if( couplind1 != couplind2)
                        GRAPH_DEL_EDGE(g, SCIPhashmapGetImageInt(vartosuccmap, (void*)(size_t)couplind1),
                                       SCIPhashmapGetImageInt(vartosuccmap, (void*)(size_t)couplind2));
                  }
               }
            }
         }
      }
   }

   /* Calculate max variables set to 1 at once as the max independent set cardinality */

   /* Set cliquer options */
   cl_opts.reorder_function = reorder_by_default; /* default: reorder_by_default */
   cl_opts.reorder_map = NULL;
   cl_opts.time_function = NULL; /* default: clique_print_time */
   cl_opts.output = NULL;
   cl_opts.user_function = NULL;
   cl_opts.user_data = NULL;
   cl_opts.clique_list = NULL;
   cl_opts.clique_list_length = 0;

   /* Find maximum weight clique using the cliquer library */
   clique = clique_find_single(g, 0, 0, FALSE, &cl_opts);

   if( !usemultipl )
      denominator = set_size(clique);
   else if( !graph_weighted(g) )
      denominator = set_size(clique) * g->weights[0];
   else
   {
      i = -1;
      denominator = 0;
      while( (i = set_return_next(clique, i)) >= 0 )
         denominator += g->weights[i];
   }

   /* Distribute coeff. of act. var to all successor variables. */
   coeftodistr = getAggrObjCoef(varsincouplings[actvarind], nlinkedvars, nvarsincouplings, aggrobjcoef);
   for( i = 0; i < nvarsincouplings; i++ )
   {
      if( couplingmatrix[actvarind][i] > 0 &&
          isCouplingRelevant(scip, i, couplingmatrix, varsincouplings, nvarsincouplings, nlinkedvars, aggrobjcoef) )
      {
         actvarcoef = getAggrObjCoef(varsincouplings[i], nlinkedvars, nvarsincouplings, aggrobjcoef);
         setAggrObjCoef(scip, varsincouplings[i],
                        actvarcoef + coeftodistr * ((SCIP_Real)(usemultipl ? varmultiplicities[i] : 1.0) / denominator),
                        linkedvars, nlinkedvars, linkmatrix, aggrobjcoef);
      }
   }

   /* Free memory */
   set_free(clique);
   graph_free(g);
   SCIPhashmapFree(&vartosuccmap);
}

/*
 * Idea - Distribution Strategy 3: Uniform Coefficient Share Distribution. (Fastest heuristic implemented)
 *
 * Uniformly distributes the objective coefficient of a variable among all its relevant coupled successor variables.
 * Each successor receives an equal share, ensuring no overestimation of the redistributed objective coefficient.
 */

/** Coupling(-like) Variable Objective Distribution Strategy 3: Uniform Coefficient Share Distribution. */
static INLINE
void objCoefDistrHeurUniform(
   SCIP*                 scip,               /**< The problem instance */
   int                   actvarind,          /**< Index of actual variable */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int                   nvarsincouplings,   /**< Index of varsincouplings array */
   int*                  varmultiplicities,  /**< Array holding number of vars represented by digraph node */
   SCIP_Bool             usemultipl,         /**< Wether varmultiplicities should be used for weighting */
   SCIP_VAR**            linkedvars,         /**< Array of variables that are linked by eq-constraints */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   int**                 linkmatrix,         /**< Matrix indicating which variables are linked by eq-constraints */
   SCIP_Real*            aggrobjcoef         /**< Array of aggregated objective coefficients. */
)
{
   int            i;
   int            nsuccessors;
   SCIP_Real      actvarcoef;
   SCIP_Real      coeftodistr;

   nsuccessors = getNSuccessorsRelevant(scip, actvarind, couplingmatrix, varsincouplings, nvarsincouplings,
                                        nlinkedvars, aggrobjcoef, varmultiplicities, usemultipl);

   /* Distribute coeff. of act. var to all successor variables. */
   coeftodistr = getAggrObjCoef(varsincouplings[actvarind], nlinkedvars, nvarsincouplings, aggrobjcoef);
   for( i = 0; i < nvarsincouplings; i++ )
   {
      if( couplingmatrix[actvarind][i] > 0 &&
          isCouplingRelevant(scip, i, couplingmatrix, varsincouplings, nvarsincouplings, nlinkedvars, aggrobjcoef) )
      {
         actvarcoef = getAggrObjCoef(varsincouplings[i], nlinkedvars, nvarsincouplings, aggrobjcoef);
         setAggrObjCoef(scip, varsincouplings[i],
                        actvarcoef + coeftodistr * ((SCIP_Real)(usemultipl ? varmultiplicities[i] : 1.0) / nsuccessors),
                        linkedvars, nlinkedvars, linkmatrix, aggrobjcoef);
      }
   }
}


/** Distributes the objective coefficient of coupling(-like) vars to all other vars occurring in those constraints.
 *  -> Recursive function: Should be called through non-recursive wrapper function 'distributeObjCoef()'. */
static
void distributeObjCoefRec(
   SCIP*                 scip,               /**< The problem instance */
   int                   actvarind,          /**< Index of actual variable */
   SCIP_Bool*            isdistributed,      /**< Array storing if a coupled variable was visited */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int                   nvarsincouplings,   /**< Index of varsincouplings array */
   int*                  varmultiplicities,  /**< Array holding number of vars represented by digraph node */
   SCIP_Bool             usemultipl,         /**< Wether varmultiplicities should be used for weighting */
   SCIP_VAR**            linkedvars,         /**< Array of variables that are linked by eq-constraints */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   int**                 linkmatrix,         /**< Matrix indicating which variables are linked by eq-constraints */
   SCIP_Real*            aggrobjcoef,        /**< Array of aggregated objective coefficients. */
   SCIP_CONS**           constraints,        /**< Array containing pointers to SCIP constraints of pricing problem */
   int                   nconss,             /**< Index of constraints array */
   int*                  couplingcoefindices,/**< Array for coupling coefficient of each constraint (if coupling) */
   CLIQUER_CONSTYPE*     cliquerconstypes,   /**< Array holding constraint types (specific to this solver) */
   int                   selecteddistrheur   /**< Number of selected obj coef distribution heuristic. */
   )
{
   int            i;

   if( isdistributed[actvarind] )
   {
      /* Cycle in digraph detected. Could improve handling, but because very unlikely we just end recursion here.
       * Cyle of form: x <= y, y <= z, z <= x
       * (-> i.e. it follows x = y = z; aggegated coef. could be distributed equally among other coupled vars.) */
      return;
   }

   isdistributed[actvarind] = TRUE;               /* mark actual variable visited. */

   /* Recursive case : if var has predecessor(s) - visit unvisited predecessor(s) first. */
   for( i = 0; i < nvarsincouplings; i++ )
   {
      if( actvarind == i )
         continue;
      if( couplingmatrix[i][actvarind] > 0 && !isdistributed[i] )
         distributeObjCoefRec(scip, i, isdistributed, couplingmatrix, varsincouplings, nvarsincouplings,
                              varmultiplicities, usemultipl, linkedvars, nlinkedvars, linkmatrix, aggrobjcoef,
                              constraints, nconss, couplingcoefindices, cliquerconstypes, selecteddistrheur);
   }

   /* Base case: All predecessors are distributed. */
   /* If now has positive (aggr.) obj. coeff.: Distribute coeff. of act. var to all successor variables. */
   if( SCIPisGT(scip, getAggrObjCoef(varsincouplings[actvarind], nlinkedvars, nvarsincouplings, aggrobjcoef), 0))
   {
      switch( selecteddistrheur )
      {
         case 1:
            objCoefDistrHeurNatural(scip, actvarind, couplingmatrix, varsincouplings, nvarsincouplings,
                                    varmultiplicities, usemultipl, linkedvars, nlinkedvars, linkmatrix, aggrobjcoef,
                                    constraints, nconss, couplingcoefindices, cliquerconstypes);
            break;
         case 2:
            objCoefDistrHeurIS(scip, actvarind, couplingmatrix, varsincouplings, nvarsincouplings, varmultiplicities,
                               usemultipl, linkedvars, nlinkedvars, linkmatrix, aggrobjcoef, constraints, nconss,
                               couplingcoefindices, cliquerconstypes);
            break;
         case 3:
            objCoefDistrHeurUniform(scip, actvarind, couplingmatrix, varsincouplings, nvarsincouplings,
                                    varmultiplicities, usemultipl, linkedvars, nlinkedvars, linkmatrix, aggrobjcoef);
            break;
         default:
            break;
      }
   }
}

/*
 * Idea - Objective Coefficient Distribution of Coupling Variables Coefficient:
 *
 * ** Problem: **
 * As the cliquer solver can only handle non-negative integer weights, the objective coefficients are (besides
 * scaled) inverted and all negative coefficients are just set to 1 (heuristically!).
 * Thus, after the independent set problem is transformed (heuristically) to a clique problem, the objective
 * coefficients of coupling variables worsening the solutions objective value are not properly reflected in the
 * weights of the max weighted clique problem.
 *
 * ** Correction Attempt: **
 * The objective coefficient w of a coupling variable y (constraints of form: x_1 + ... + x_n <= c*y) can be distributed
 * among all other variables x_1, ..., x_n in the constraint(s) to get an objective that is closer to the actual one.
 *
 * ** Implementation: **
 * We implement this by building and utilizing a digraph to process the distribution hierarchically. I.e., if there
 * exist constraints of the form z_1 + z_2 <= z_3 and z_3 + z_4 <= z_5, we first want to distribute the coefficient of
 * z_5 to the ones of z_3 and z_4 and only then distribute the (aggregated) objective coefficient of z_3 to the ones of
 * z_1 and z_2.
 *
 * Furthermore, we only distribute coefficients that worsen the solution (negative ones after inversion). Also, we only
 * distribute to relevant variables, as the others should already be never chosen because there objective coefficient
 * suggests so.
 */

/** Distributes the objective coefficient of coupling(-like) vars to all other vars occurring in those constraints. */
static
void distributeObjCoef(
   SCIP*                 scip,               /**< The problem instance */
   SCIP_CONS**           constraints,        /**< Array containing pointers to SCIP constraints of pricing problem */
   int                   nconss,             /**< Index of constraints array */
   int*                  couplingcoefindices,/**< Array for coupling coefficient of each constraint (if coupling) */
   CLIQUER_CONSTYPE*     cliquerconstypes,   /**< Array holding constraint types (specific to this solver) */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are linked by a node */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   int                   nvarsincouplings,   /**< Index of varsincouplings array */
   SCIP_VAR**            linkedvars,         /**< Array of variables that are linked by eq-constraints */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   int**                 linkmatrix,         /**< Matrix indicating which variables are linked by eq-constraints */
   SCIP_Real*            aggrobjcoef,        /**< Array of aggregated objective coefficients. */
   int                   selecteddistrheur,  /**< Number of selected obj coef distribution heuristic. */
   SCIP_Bool             usemultipl          /**< Wether varmultiplicities should be used for weighting */
   )
{
   int            i;
   SCIP_Bool*     isdistributed;                   /* array storing if a coupled variable was visited */
   int*           varmultiplicities;               /* array holding the number of linked (represented) vars per node */

   /* Local memory allocation */
   varmultiplicities = NULL;
   SCIP_CALL_ABORT( SCIPallocClearBufferArray(scip, &isdistributed, nvarsincouplings) );
   if( usemultipl )
   {
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &varmultiplicities, nvarsincouplings) );
      initVarMultiplicities(scip, linkmatrix, linkedvars, nlinkedvars, couplingmatrix, varsincouplings,
                            nvarsincouplings, varmultiplicities);
   }

   for( i = 0; i < nvarsincouplings; i++ )
   {
      if( !isdistributed[i] &&
          hasSuccessorRel(scip, i, couplingmatrix, varsincouplings, nvarsincouplings, nlinkedvars, aggrobjcoef) )
      {
         distributeObjCoefRec(scip, i, isdistributed, couplingmatrix, varsincouplings, nvarsincouplings,
                              varmultiplicities, usemultipl, linkedvars, nlinkedvars, linkmatrix, aggrobjcoef,
                              constraints, nconss, couplingcoefindices, cliquerconstypes, selecteddistrheur);
         continue;
      }
      isdistributed[i] = TRUE;               /* mark actual variable visited. */
   }

   /* Free local memory */
   if( usemultipl )
      SCIPfreeBufferArray(scip, &varmultiplicities);
   SCIPfreeBufferArray(scip, &isdistributed);
}

/**
 * Determine cliquer constraint type and save it in the cliquerconstype array.
 * Also, build same-constraint equality graph and coupling digraph.
 * @return SCIP bool that is false in case a constraint is encountered that cannot be handled,
 *          true if propagation was successful.
 */
static
SCIP_Bool determineCliquerConsTypes(
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   SCIP_CONS**           constraints,        /**< Array containing pointers to SCIP constraints of pricing problem */
   SCIP_VAR**            linkedvars,         /**< Array of variables that are linked by eq-constraints */
   SCIP_VAR**            vconsvars,          /**< Array to store variables of a varbound constraint. */
   int*                  markedconsindices,  /**< Array containing pointers to SCIP constraints to mark them */
   int**                 linkmatrix,         /**< Matrix indicating which variables are linked by a node */
   int*                  couplingcoefindices,/**< Array for coupling coefficient of each constraint (if coupling) */
   int*                  nlinkedvars,        /**< Index of linkedvars array */
   int                   nconss,             /**< Index of constraints array */
   int*                  markedcount,        /**< Index of markedconstraints array */
   CLIQUER_CONSTYPE*     cliquerconstypes    /**< Array holding constraint types (specific to this solver) */
   )
{
   /* Local variables. */
   SCIP_CONSHDLR*    conshdlr;
   SCIP_Real*        consvals;
   SCIP_Bool         retcode;
   int               nvars;

   int               i;
   int               j;

   /* Loop for checking and saving the constraint types. This is done to ease the handling of the cases later on. */
   /* Also the case of the occurrence of constraints that can not be handled by the solver is covered. */
   /* Also, the equality graph is built through updating the linkmatrix every time a "same"-constraint is encountered. */
   for( i = 0; i < nconss; ++i )
   {
      assert(constraints[i] != NULL);
      conshdlr = SCIPconsGetHdlr(constraints[i]);
      assert(conshdlr != NULL);

      /* The constraint may be of type 'linear' */
      if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
      {
         consvals = SCIPgetValsLinear(pricingprob, constraints[i]);
         if( !SCIPisEQ(pricingprob, SCIPgetLhsLinear(pricingprob, constraints[i]), SCIPgetRhsLinear(pricingprob, constraints[i])) )
         {
            /* Check if we have an IS constraint */
            if( SCIPgetNVarsLinear(pricingprob, constraints[i]) == 2 && SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob, constraints[i]), 1) )
            {
               cliquerconstypes[i] = LINEAR_IS;
            }
            /* Handle other constraints that behave like IS constraints, i.e. cx+dy<=rhs with c+d>rhs, c>0, d>0 */
            else if( SCIPgetNVarsLinear(pricingprob, constraints[i]) == 2 && consvals[0] > 0 && consvals[1] > 0 &&
                     SCIPisLT(pricingprob, SCIPgetRhsLinear(pricingprob, constraints[i]), consvals[0] + consvals[1]) &&
                     !SCIPisLT(pricingprob, SCIPgetRhsLinear(pricingprob, constraints[i]), consvals[0]) &&
                     !SCIPisLT(pricingprob, SCIPgetRhsLinear(pricingprob, constraints[i]), consvals[1]) )
            {
               cliquerconstypes[i] = LINEAR_IS_LIKE;
            }
            else
            {
               /* The current constraint is no linear IS constraint */
               SCIPgetConsNVars(pricingprob, constraints[i], &nvars, &retcode);

               /* Check the coefficients of the variables in the constraint */
               for( j = 0; j < nvars; ++j )
               {
                  if( consvals[j] != 1 && (couplingcoefindices[i] == -1) )
                  {
                     couplingcoefindices[i] = j;
                  }
                  else if( consvals[j] != 1 && couplingcoefindices[i] != -1 )
                  {
                     /* More than one variable has a coefficient unequal to 1 */
                     SCIPdebugMessage("Exit: More than one coefficient unequal 1 in linear non-IS constraint.\n");
                     return FALSE;

                     /*
                      * Could handle other types of constraints similar to coupling constraints.
                      * -> E.g.: One var. coeff. < 0 and this var is fixed to 0: Others must also be fixed to 0.
                      *          Otherwise, cannot handle!
                      * -> To handle those, they must be identified and marked somehow to check if the coeff. is fixed
                      *    to 0 after propagation.
                      *    If not, the constraint cannot be handled. -> TERMINATE with GCG_PRICINGSTATUS_NOTAPPLICABLE
                      */
                  }
               }
               /* Check if we have a clique constraint (rhs 1 and coefficients 1) */
               if( (couplingcoefindices[i] == -1) && SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob, constraints[i]), 1) )
               {
                  cliquerconstypes[i] = LINEAR_CLIQUE;
               }
               /* Check if we have a coupling constraint (rhs 0) */
               else if( (couplingcoefindices[i] != -1) && SCIPisEQ(pricingprob, SCIPgetRhsLinear(pricingprob, constraints[i]), 0.0) )
               {
                  /* Special case: The coupling constraint is purely decorative (coefficient + 1 of coupling var >= #vars)*/
                  if( abs(consvals[couplingcoefindices[i]]) + 1 >= nvars )
                  {
                     cliquerconstypes[i] = LINEAR_COUPLING_DECORATIVE;
                  }
                  /* Special case: The coefficient is -1, we treat the case like a clique constraint. */
                  else if( abs(consvals[couplingcoefindices[i]]) == 1 )
                  {
                     cliquerconstypes[i] = LINEAR_COUPLING_CLIQUE;
                  }
                  else
                  {
                     /* Coupling coefficient is between 1 and npricingprobvars. */
                     SCIPdebugMessage("Exit: Coupling coefficient unhandled, coef: %g.\n", consvals[couplingcoefindices[i]]);
                     return FALSE;
                  }
               }
               else
               {
                  /* Constraint is neither a coupling nor a clique constraint */
                  SCIPdebugMessage("Exit: Unhandled linear constraint.\n");
                  return FALSE;
               }
            }
         }
         else
         {
            /* Constraint is a linear equality constraint */
            SCIPdebugMessage("Exit: Unhandled linear constraint: Equality constraint.\n");
            return FALSE;
         }
      }
         /* Constraint may be of type varbound: lhs <= x + c*y <= rhs */
      else if( strcmp(SCIPconshdlrGetName(conshdlr), "varbound") == 0 )
      {
         vconsvars[0] = SCIPgetVarVarbound(pricingprob, constraints[i]);
         vconsvars[1] = SCIPgetVbdvarVarbound(pricingprob, constraints[i]);

         /* Check for "same"-constraints present in Ryan-Foster-Branching and save the links between the variables. */
         /* These are constraints of type: constraint of type x = y (lhs = rhs = 0 and c = -1)*/
         if( SCIPisEQ(pricingprob, SCIPgetLhsVarbound(pricingprob, constraints[i]), SCIPgetRhsVarbound(pricingprob, constraints[i])) )
         {
            /* c == -1, thus variables have to become both 0 or both 1 */
            if( (SCIPgetRhsVarbound(pricingprob, constraints[i]) == 0) && (SCIPgetVbdcoefVarbound(pricingprob, constraints[i]) == -1) )
            {
               cliquerconstypes[i] = VARBND_SAME;

               /* Build the equality graph through updating the linkmatrix. */
               updateVarLinks(pricingprob, linkmatrix, vconsvars[0], vconsvars[1], linkedvars, nlinkedvars);
               /* Since the vars may not be part of the graph, we have to be able to set their solval later, thus we save the constraint index */
               markedconsindices[*markedcount] = i;
               ++(*markedcount);
            }
            else
            {
               /* RHS is unequal 0 and unequal 1 */
               SCIPdebugMessage("Exit: Unhandled equality constraint, c: %g, rhs: %g.\n",
                                SCIPgetVbdcoefVarbound(pricingprob, constraints[i]), SCIPgetRhsVarbound(pricingprob, constraints[i]));
               return FALSE;
            }
         }

         /* Check value of rhs to be 0 and of c to be <= -1 */
         if( SCIPisInfinity(pricingprob, -SCIPgetLhsVarbound(pricingprob, constraints[i])) )
         {
            if( SCIPisEQ(pricingprob, SCIPgetRhsVarbound(pricingprob, constraints[i]), 0) )
            {
               if( SCIPisLT(pricingprob, SCIPgetVbdcoefVarbound(pricingprob, constraints[i]), -1) ||
                   SCIPisEQ(pricingprob, SCIPgetVbdcoefVarbound(pricingprob, constraints[i]), -1) )
               {
                  cliquerconstypes[i] = VARBND_STD;
               }
               else
               {
                  /* Coefficient c of varbound is > -1 and we do not have an IS constraint*/
                  SCIPdebugMessage("Exit: Coefficient of Varbound unhandled Rhs: %g, Coeff: %g.\n",
                                   SCIPgetRhsVarbound(pricingprob, constraints[i]), SCIPgetVbdcoefVarbound(pricingprob, constraints[i]));
                  return FALSE;
               }
            }
            /*
             * Rhs of varbound unequal to 0.
             * It may still be the case that we have an IS constraint with a non-linear handler.
             * The constraint may also be of the form c + 1 > rhs and c < rhs, i.e. a non-standard IS-constraint.
             * We treat these cases like a regular IS constraint.
             */
            else if( (SCIPisEQ(pricingprob, SCIPgetRhsVarbound(pricingprob, constraints[i]), 1) &&
                      SCIPisEQ(pricingprob, SCIPgetVbdcoefVarbound(pricingprob, constraints[i]), 1)) ||
                     (SCIPisLT(pricingprob, SCIPgetRhsVarbound(pricingprob, constraints[i]), SCIPgetVbdcoefVarbound(pricingprob, constraints[i]) + 1) &&
                      SCIPisLT(pricingprob, SCIPgetVbdcoefVarbound(pricingprob, constraints[i]), SCIPgetRhsVarbound(pricingprob, constraints[i]))) )
            {
               cliquerconstypes[i] = VARBND_IS;
            }
            else
            {
               /* Rhs of varbound unequal to 0 and no IS constraint*/
               SCIPdebugMessage("Exit: Rhs of Varbound unhandled, Rhs: %g, Coeff:%g.\n",
                                SCIPgetRhsVarbound(pricingprob, constraints[i]), SCIPgetVbdcoefVarbound(pricingprob, constraints[i]));
               return FALSE;
            }
         }
            /* We may have a varbound constraint of type x + cy == rhs */
         else if( SCIPisEQ(pricingprob, SCIPgetLhsVarbound(pricingprob, constraints[i]), SCIPgetRhsVarbound(pricingprob, constraints[i])) )
         {
            /* If the rhs is 0 and c == -1, both variables have to be set to 0 or to 1 */
            if( !((SCIPgetRhsVarbound(pricingprob, constraints[i]) == 0) && (SCIPgetVbdcoefVarbound(pricingprob, constraints[i]) == -1)) )
            {
               /* RHS is unequal 0 and unequal 1 */
               SCIPdebugMessage("Exit: Unhandled equality constraint, c: %g, rhs: %g.\n",
                                SCIPgetVbdcoefVarbound(pricingprob, constraints[i]), SCIPgetRhsVarbound(pricingprob, constraints[i]));
               return FALSE;
            }
         }
         else
         {
            /* We have a varbound of type lhs <= x + c*y */
            SCIPdebugMessage("Exit: Varbound of type lhs <= x+c*y, c: %g, rhs: %g.\n",
                             SCIPgetVbdcoefVarbound(pricingprob, constraints[i]), SCIPgetRhsVarbound(pricingprob, constraints[i]));
            SCIPdebugMessage("Constraint handler: %s\n", SCIPconshdlrGetName(conshdlr));
            return FALSE;
         }
      }
      else
      {
         /* Constraint handler neither linear nor varbound */
         SCIPdebugMessage("Exit: Unhandled constraint handler: %s \n", SCIPconshdlrGetName(conshdlr));
         return FALSE;
      }
   }

#ifdef SCIP_DEBUG
   SCIPdebugMessage("Overview over instances constraint types:\n");
   const char* typeNames[] = {
      "LINEAR_IS", "LINEAR_IS_LIKE", "LINEAR_CLIQUE", "LINEAR_COUPLING_DECORATIVE",
      "LINEAR_COUPLING_CLIQUE", "VARBND_SAME", "VARBND_STD", "VARBND_IS"
   };
   int typeCount[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
   for( i = 0; i < nconss; i++ )
      typeCount[cliquerconstypes[i]]++;
   for( i = LINEAR_IS; i <= VARBND_IS; i++ )
   {
      SCIPdebugMessage("   Type '%s' : %i \n", typeNames[i], typeCount[i]);
   }
#endif

   /* Has no invalid constraint. */
   return TRUE;
}

/**
 * Propagate fixings of variables through constraints until the set of fixed variables is stable.
 * @return SCIP bool that is false in case the problem is infeasible, true if propagation was successful.
 */
static
SCIP_Bool propagateVariablefixings(
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   SCIP_CONS**           constraints,        /**< Array containing pointers to SCIP constraints of pricing problem */
   SCIP_VAR**            linkedvars,         /**< Array of variables that are linked by eq-constraints */
   SCIP_VAR**            vconsvars,          /**< Array to store variables of a varbound constraint. */
   SCIP_VAR**            varsincouplings,    /**< Array of variables that are involved in couplings */
   SCIP_Real*            solvals,            /**< Array holding the current solution values of all problem variables */
   int**                 couplingmatrix,     /**< Matrix indicating which variables are involved in couplings */
   int**                 linkmatrix,         /**< Matrix indicating which variables are linked by a node */
   int*                  nvarsincouplings,   /**< Index of varsincouplings array */
   int*                  couplingcoefindices,/**< Array for coupling coefficient of each constraint (if coupling) */
   int*                  consvarsfixedcount, /**< Array for counting how many variables are fixed per constraint */
   int                   nlinkedvars,        /**< Index of linkedvars array */
   int                   nconss,             /**< Index of constraints array */
   int*                  nfixedvars,         /**< Integer counting the number of currently fixed problem variables */
   CLIQUER_CONSTYPE*     cliquerconstypes    /**< Array holding constraint types (specific to this solver) */
   )
{
   /* Local variables. */
   SCIP_CONSHDLR*    conshdlr;
   SCIP_VAR**        lconsvars;
   SCIP_Bool         retcode;
   int               nvars;
   int               nvarsfixedtoone;
   int               vartoset;

   int               i;
   int               j;
   int               prevfixed;

   /* Compute implied variable fixings. */
   /* This is done by propagating the fixings already found over the constraints. */
   /* It is stopped once the set of fixed variables becomes stable across one iteration. */
   prevfixed = -1;        /* Need at least one iteration (because it is checked if linked variables appear in IS-constraint, i.e., x = y and x + y <= 1). */
   while( prevfixed < *nfixedvars ) {

      /* We still have a fixed variable to be processed. Iterate through constraints. */
      prevfixed = *nfixedvars;
      for( i = 0; i < nconss; i++ )
      {
         assert(constraints[i] != NULL);
         conshdlr = SCIPconsGetHdlr(constraints[i]);
         assert(conshdlr != NULL);

         /* Variables do not know in which constraints they appear. */
         /* Therefore, we count how many variables are fixed per constraint to skip constraints which have only fixed variables. */
         /* The constraint is checked if it is consistent with the fixings. Afterwards, the counter is updated. */
         /* This ensures every constraint is checked for consistency once before we skip it. */

         /* Check nature of the constraint */

         /* The constraint may be of type 'linear' */
         if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
         {
            lconsvars = SCIPgetVarsLinear(pricingprob, constraints[i]);

            /* Add coupling cons vars to coupling digraph. */
            if( cliquerconstypes[i] == LINEAR_COUPLING_CLIQUE || cliquerconstypes[i] == LINEAR_COUPLING_DECORATIVE )
            {
               SCIPgetConsNVars(pricingprob, constraints[i], &nvars, &retcode);
               updateCouplingDiGraph(pricingprob, lconsvars, nvars, lconsvars[couplingcoefindices[i]], couplingmatrix,
                                     varsincouplings, nvarsincouplings, linkedvars, nlinkedvars, linkmatrix);
            }

            /* If all variables are fixed, constraint can be skipped */
            if( consvarsfixedcount[i] == SCIPgetNVarsLinear(pricingprob, constraints[i]) )
               continue;

            if( cliquerconstypes[i] == LINEAR_IS || cliquerconstypes[i] == LINEAR_IS_LIKE )
            {
               /* Propagate variable fixings through IS-constraint. */
               if( solvals[SCIPvarGetProbindex(lconsvars[0])] == 1 && solvals[SCIPvarGetProbindex(lconsvars[1])] == 1 )
               {
                  /* Both variables are fixed to 1 which contradicts the IS constraint. -> Infeasible. */
                  SCIPdebugMessage("Exit: Both variables in IS-constraint fixed to 1.\n");
                  return FALSE;
               }
               else if( solvals[SCIPvarGetProbindex(lconsvars[0])] == 1 && solvals[SCIPvarGetProbindex(lconsvars[1])] == -1 )
               {
                  /* One variable0 is fixed to 1 -> fix variable1 to 0. */
                  solvals[SCIPvarGetProbindex(lconsvars[1])] = 0;
                  (*nfixedvars)++;
                  consvarsfixedcount[i] = 2;
               }
               else if( solvals[SCIPvarGetProbindex(lconsvars[0])] == -1 && solvals[SCIPvarGetProbindex(lconsvars[1])] == 1 )
               {
                  /* One variable1 is fixed to 1 -> fix variable0 to 0. */
                  solvals[SCIPvarGetProbindex(lconsvars[0])] = 0;
                  (*nfixedvars)++;
                  consvarsfixedcount[i] = 2;
               }
               else if( solvals[SCIPvarGetProbindex(lconsvars[0])] == -1 &&
                        solvals[SCIPvarGetProbindex(lconsvars[1])] == -1 &&
                        areVarsLinked(pricingprob, linkmatrix, lconsvars[0], lconsvars[1], linkedvars, nlinkedvars) )
               {
                  /* The two variables are linked and appear in an IS-constraint, i.e., x = y and x + y <= 1.
                   * -> Both variables must be fixed to 0. Thus calling the setter for one is sufficient */
                  setLinkedSolvals(pricingprob, solvals, linkmatrix, linkedvars, nlinkedvars, lconsvars[0], 0.0);
                  *nfixedvars += 2;
                  consvarsfixedcount[i] = 2;
               }
            }
            else
            {
               /* The current constraint is no linear IS constraint */
               SCIPgetConsNVars(pricingprob, constraints[i], &nvars, &retcode);
               nvarsfixedtoone = 0;

               /* Count the number of variables with a fixed value of 1. */
               for( j = 0; j < nvars; ++j )
               {
                  if( solvals[SCIPvarGetProbindex(lconsvars[j])] == 1 )
                     nvarsfixedtoone++;
               }

               if( cliquerconstypes[i] == LINEAR_CLIQUE && nvarsfixedtoone > 1 )
               {
                  /* More than one variable has a value fixed to 1 */
                  SCIPdebugMessage("Exit: More than one variable with value fixed to 1 in clique constraint.\n");
                  return FALSE;
               }
               else if( cliquerconstypes[i] == LINEAR_COUPLING_CLIQUE && nvarsfixedtoone > 2 )
               {
                  SCIPdebugMessage("Exit: To many variable values fixed to 1 in coupling constraint with coupling variable value fixed to 1.\n");
                  return FALSE;
               }
               else if( (cliquerconstypes[i] == LINEAR_COUPLING_CLIQUE || cliquerconstypes[i] == LINEAR_COUPLING_DECORATIVE) &&
                        solvals[SCIPvarGetProbindex(lconsvars[couplingcoefindices[i]])] == 0 &&
                        nvarsfixedtoone >= 1 )
               {
                  SCIPdebugMessage("Exit: To many variable values fixed to 1 in coupling constraint with coupling variable value fixed to 0.\n");
                  return FALSE;
               }
               else if( (cliquerconstypes[i] == LINEAR_CLIQUE && nvarsfixedtoone == 1) ||
                        (cliquerconstypes[i] == LINEAR_COUPLING_CLIQUE &&
                         nvarsfixedtoone == 2 &&
                         solvals[SCIPvarGetProbindex(lconsvars[couplingcoefindices[i]])] == 1) ||
                        ((cliquerconstypes[i] == LINEAR_COUPLING_CLIQUE || cliquerconstypes[i] == LINEAR_COUPLING_DECORATIVE) &&
                         solvals[SCIPvarGetProbindex(lconsvars[couplingcoefindices[i]])] == 0) )
               {
                  /* We have a clique constraint with exactly one variables value fixed to 1. */
                  /* Or a coupling constraint that can be handled like a clique constraint with exactly one variable value fixed to 1. */
                  /* Or a coupling constraint (clique or decorative) that has the coupling variable fixed to 0. */

                  /* In all these cases: All other involved variables need to be fixed to 0. */
                  for( j = 0; j < nvars; j++ )
                  {
                     /* The solvals of the other variables are either 0 or -1 */
                     /* Only fix to 0 and add to fixed variable array if value is -1 */
                     if( solvals[SCIPvarGetProbindex(lconsvars[j])] == -1 )
                     {
                        solvals[SCIPvarGetProbindex(lconsvars[j])] = 0;
                        (*nfixedvars)++;
                     }
                  }
                  consvarsfixedcount[i] = nvars; /* All variables of this constraint are fixed now. */
               }
               else if( (cliquerconstypes[i] == LINEAR_COUPLING_CLIQUE || cliquerconstypes[i] == LINEAR_COUPLING_DECORATIVE) &&
                        nvarsfixedtoone == 1 &&
                        solvals[SCIPvarGetProbindex(lconsvars[couplingcoefindices[i]])] == -1 )
               {
                  /* We have a coupling constraint with one variable (different from the coupling variable!) fixed to 1.
                   * And the coupling variable is unfixed. Then the coupling variable needs to be fixed to 1 too. */
                  solvals[SCIPvarGetProbindex(lconsvars[couplingcoefindices[i]])] = 1;
                  (*nfixedvars)++;

                  /* In case of a clique constraint, we can fix all other variables than the (now 2) fixed ones to 0. */
                  if( cliquerconstypes[i] == LINEAR_COUPLING_CLIQUE )
                  {
                     for( j = 0; j < nvars; j++ )
                     {
                        /* The solvals of the other variables are either 0 or -1 */
                        /* Only fix to 0 and add to fixed variable array if value is -1 */
                        if( solvals[SCIPvarGetProbindex(lconsvars[j])] == -1 )
                        {
                           solvals[SCIPvarGetProbindex(lconsvars[j])] = 0;
                           (*nfixedvars)++;
                        }
                     }
                     consvarsfixedcount[i] = nvars; /* All variables of this constraint are fixed now. */
                  }
               }
            }
         }
         /* Constraint may be of type varbound: lhs <= x + c*y <= rhs */
         else if( strcmp(SCIPconshdlrGetName(conshdlr), "varbound") == 0 )
         {
            vconsvars[0] = SCIPgetVarVarbound(pricingprob, constraints[i]);
            vconsvars[1] = SCIPgetVbdvarVarbound(pricingprob, constraints[i]);

            /* Add coupling-like varbound cons vars to coupling digraph. */
            if( cliquerconstypes[i] == VARBND_STD )
               updateCouplingDiGraph(pricingprob, vconsvars, 2, vconsvars[1], couplingmatrix, varsincouplings,
                                     nvarsincouplings, linkedvars, nlinkedvars, linkmatrix);

            if( consvarsfixedcount[i] == 2 )
               continue; /* If all variables are fixed, constraint can be skipped */

            if( cliquerconstypes[i] == VARBND_SAME )
            {
               /* Propagate variable fixings through same-constraint. */
               if( solvals[SCIPvarGetProbindex(vconsvars[0])] >= 0 && solvals[SCIPvarGetProbindex(vconsvars[1])] >= 0 &&
                   solvals[SCIPvarGetProbindex(vconsvars[0])] != solvals[SCIPvarGetProbindex(vconsvars[1])] )
               {
                  /* One variable is fixed to 1, the other to 0. -> Infeasible. */
                  SCIPdebugMessage("Exit: Variables in same-constraint are fixed to different values.\n");
                  return FALSE;
               }
               else if( solvals[SCIPvarGetProbindex(vconsvars[0])] >= 0 && solvals[SCIPvarGetProbindex(vconsvars[1])] == -1 )
               {
                  /* Fix (the unfixed) variable1 to the value of variable0 */
                  solvals[SCIPvarGetProbindex(vconsvars[1])] = solvals[SCIPvarGetProbindex(vconsvars[0])];
                  (*nfixedvars)++;
                  consvarsfixedcount[i] = 2; /* All variables of this constraint are fixed now. */
               }
               else if( solvals[SCIPvarGetProbindex(vconsvars[0])] == -1 && solvals[SCIPvarGetProbindex(vconsvars[1])] >= 0 )
               {
                  /* Fix (the unfixed) variable0 to the value of variable1 */
                  solvals[SCIPvarGetProbindex(vconsvars[0])] = solvals[SCIPvarGetProbindex(vconsvars[1])];
                  (*nfixedvars)++;
                  consvarsfixedcount[i] = 2; /* All variables of this constraint are fixed now. */
               }
            }
            /* From here on we may have a varbound constraint with x + c*y <= b. */
            else
            {
               if( solvals[SCIPvarGetProbindex(vconsvars[0])] == 1 &&
                   ((cliquerconstypes[i] == VARBND_STD && solvals[SCIPvarGetProbindex(vconsvars[1])] == 0) ||
                    (cliquerconstypes[i] == VARBND_IS && solvals[SCIPvarGetProbindex(vconsvars[1])] == 1)) )
               {
                  if( solvals[SCIPvarGetProbindex(vconsvars[1])] == 0 )
                     SCIPdebugMessage("Exit: x fixed to 1, y fixed to 0 in varbound constraint.\n");
                  if( solvals[SCIPvarGetProbindex(vconsvars[1])] == 1 )
                     SCIPdebugMessage("Exit: Both variables fixed to 1 in non-linear handler IS-constraint.\n");
                  return FALSE;
               }
               else if( cliquerconstypes[i] == VARBND_STD &&
                        ((solvals[SCIPvarGetProbindex(vconsvars[0])] == 1 && solvals[SCIPvarGetProbindex(vconsvars[1])] == -1) ||
                         (solvals[SCIPvarGetProbindex(vconsvars[0])] == -1 && solvals[SCIPvarGetProbindex(vconsvars[1])] == 0)) )
               {
                  /* Constraint behaving like x <= c*y, c >= 1 - and one variable is already fixed to 1. */
                  /* Variable to set and value to set the variable to. */
                  if( solvals[SCIPvarGetProbindex(vconsvars[0])] == 1 && solvals[SCIPvarGetProbindex(vconsvars[1])] == -1 )
                     vartoset = 1;        /* x is fixed to 1 and y is unset -> set y to 1. */
                  else
                     vartoset = 0;        /* y is fixed to 0 and x is unset -> set x to 0. */

                  solvals[SCIPvarGetProbindex(vconsvars[vartoset])] = vartoset;
                  (*nfixedvars)++;
                  consvarsfixedcount[i] = 2; /* All variables of this constraint are fixed now. */
               }
               else if( cliquerconstypes[i] == VARBND_IS &&
                        ((solvals[SCIPvarGetProbindex(vconsvars[0])] == 1 && solvals[SCIPvarGetProbindex(vconsvars[1])] == -1) ||
                         (solvals[SCIPvarGetProbindex(vconsvars[0])] == -1 && solvals[SCIPvarGetProbindex(vconsvars[1])] == 1)) )
               {
                  /* Constraint behaving like x + y <= 1 - and one variable is already fixed to 1. */
                  /* Variable to set and value to set the variable to. */
                  if( solvals[SCIPvarGetProbindex(vconsvars[0])] == 1 && solvals[SCIPvarGetProbindex(vconsvars[1])] == -1 )
                     vartoset = 1;        /* x is fixed to 1 and y is unset -> set y to 0. */
                  else
                     vartoset = 0;        /* y is fixed to 1 and x is unset -> set x to 0. */

                  solvals[SCIPvarGetProbindex(vconsvars[vartoset])] = 0;
                  (*nfixedvars)++;
                  consvarsfixedcount[i] = 2; /* All variables of this constraint are fixed now. */
               }
               else if( cliquerconstypes[i] == VARBND_IS &&
                        solvals[SCIPvarGetProbindex(vconsvars[0])] == -1 &&
                        solvals[SCIPvarGetProbindex(vconsvars[1])] == -1 &&
                        areVarsLinked(pricingprob, linkmatrix, vconsvars[0], vconsvars[1], linkedvars, nlinkedvars) )
               {
                  /* The two variables are linked and appear in an IS-constraint, i.e., x = y and x + y <= 1.
                   * -> Both variables must be fixed to 0. Thus calling the setter for one is sufficient */
                  setLinkedSolvals(pricingprob, solvals, linkmatrix, linkedvars, nlinkedvars, vconsvars[0], 0.0);
                  *nfixedvars += 2;
                  consvarsfixedcount[i] = 2; /* All variables of this constraint are fixed now. */
               }
            }
         }
      }
   }
   /* No conflicting variable fixings encountered. */
   return TRUE;
}

/* Basic idea of the heuristic solver: The biggest independent set in a graph corresponds to the biggest clique
 * of the complement graph, for which we use the cliquer library to find it. We therefore transform the variables 
 * into graph nodes and delete the edge between two nodes if there is an independent set constraint involving both. 
 * By doing this, they cannot both be part of the maximum clique and thus not be both part of the independent set.
 * The correspondence between variables and graph nodes is done by a bijection using the indsetvars array:
 * The variable indsetvars[i] is the i-th node of the graph, indexcount keeps track of the next unmapped graph node.
 * There is also the possibility that two variables x and y are linked with an equality constraint x-y = 0 due to
 * Ryan-Foster-Branching. In this case, all linked variables are mapped to the same node. There are functions
 * to get the corresponding node index.
 *
 * Since we want to add a column with the best reduced cost, we take the objective coefficient of variables into
 * account by giving their graph nodes corresponding weights and searching for a weight-maximal clique.
 *
 * If you would like to add the handling of more types of constraints, please note that the current code 
 * assumes that at no point edges are added to the graph, except during initialisation.
 *
 * This solver is currently able to handle the following type of constraints:
 * IS-Constraints, i.e. c*x + d*y <= 1*e
 * Coupling-Constraints, i.e. v + w + x -c*y <= 0
 * Clique-Constraints, i.e. v + w + x + y <= 1
 * Same-Constraints, i.e. x - y = 0 for varbound-constraints.
 * Vbd-constraints of type x - c*y <= 0 for c <= -1
 */

/** Solve the pricing problem as an independent set problem, in an approximate way. */
static
SCIP_RETCODE solveCliquer(
   SCIP_Bool             exactly,            /**< should the pricing problem be solved to optimality or heuristically? */
   GCG*                  gcg,                /**< GCG data structure */
   SCIP*                 pricingprob,        /**< pricing problem SCIP data structure */
   GCG_SOLVERDATA*       solver,             /**< solver data structure */
   int                   probnr,             /**< problem number */
   SCIP_Real*            lowerbound,         /**< pointer to store lower bound */
   GCG_PRICINGSTATUS*    status              /**< pointer to store pricing problem status */
   )
{ /*lint -e715 */
   SCIP_CONS**       constraints;
   SCIP_CONSHDLR*    conshdlr;
   SCIP_VAR**        lconsvars;
   SCIP_VAR**        vconsvars;
   SCIP_VAR**        indsetvars;
   SCIP_VAR**        pricingprobvars;
   SCIP_VAR**        linkedvars;
   SCIP_VAR**        couplvars;
   SCIP_Real*        solvals;
   SCIP_Real*        aggrobjcoef;
   SCIP_Real         density;
   SCIP_Real         scalingfactor;
   SCIP_Bool         retcode;
   SCIP_Bool         ismemgraphallocated;
   set_t             clique;
   graph_t*          g;
   clique_options    cl_opts;
   int**             linkmatrix;
   int**             couplingmatrix;
   int*              couplingcoefindices;
   int*              consvarsfixedcount;
   int*              consvarsfixedtozerocount;
   int*              markedconsindices;
   int               nlinkedvars;
   int               ncouplvars;
   int               npricingprobvars;
   int               nvars;
   int               nconss;
   int               nedges;
   int               markedcount;
   int               indexcount;
   int               nfixedvars;
   int               nodeindex0;
   int               nodeindex1;
   int               cliqueconscount;
   CLIQUER_CONSTYPE* cliquerconstypes;
   GCG_COL*          col;

   int               i;
   int               j;
   int               k;

   assert(gcg != NULL);
   assert(pricingprob != NULL);
   assert(solver != NULL);
   assert(lowerbound != NULL);
   assert(status != NULL);

   /* Check if solver already found itself to be not applicable to actual problem. */
   if( solver->isnotapplicable[probnr] )
   {
      SCIPdebugMessage("Exit: Solver already found to be not applicable to pricing problem %i.\n", probnr);
      *status = GCG_PRICINGSTATUS_NOTAPPLICABLE;
      return SCIP_OKAY;
   }

   pricingprobvars = SCIPgetVars(pricingprob);
   npricingprobvars = SCIPgetNVars(pricingprob);

   constraints = SCIPgetConss(pricingprob);
   nconss = SCIPgetNConss(pricingprob);

   /* All variables of the problem are expected to be binary */
   if( SCIPgetNBinVars(pricingprob) < npricingprobvars )
   {
      SCIPdebugMessage("Exit: Nonbinary variables.\n");
      *status = GCG_PRICINGSTATUS_NOTAPPLICABLE;
      setProblemNotApplicable(scip, probnr, solver->isnotapplicable);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(pricingprob, &markedconsindices, nconss) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &solvals, npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &linkedvars, npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &linkmatrix, npricingprobvars) );
   for( i = 0; i < npricingprobvars; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(pricingprob, &linkmatrix[i], npricingprobvars) );
   }
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &cliquerconstypes, nconss) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &consvarsfixedcount, nconss) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &couplingcoefindices, nconss) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &vconsvars, 2) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &couplvars, npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &couplingmatrix, npricingprobvars) );
   for( i = 0; i < npricingprobvars; ++i )
   {
      SCIP_CALL( SCIPallocClearBufferArray(pricingprob, &couplingmatrix[i], npricingprobvars) );
   }


   /* Boolean variables to keep track of what was allocated. */
   ismemgraphallocated = FALSE;

   /* Used to keep track of node indizes for bijection while building the graph */
   indexcount = 0;

   /* Use to handle a rare combination of IS and varbound constraints */
   markedcount = 0;

   /* Used to keep track of the index of the linkedvars array */
   nlinkedvars = 0;

   /* Used to keep track of the number of variables that have a fixed value */
   nfixedvars = 0;

   /* Used to keep track of the number of variables that are involved in coupling constraints (or std varbnd) */
   ncouplvars = 0;

   /* Build complementary graph by first creating a complete graph and then deleting edges of IS constraints. */
   /* Size is first chosen to be maximal and then later cropped down to the actual number of nodes. */
   /* Initialize the linkmatrix array. */
   /* Initialize the solvals array. */
   g = graph_new(npricingprobvars);
   for( i = 0; i < npricingprobvars; ++i )
   {
      for( j = 0; j < npricingprobvars; ++j )
      {
         if( i < j )
         {
            GRAPH_ADD_EDGE(g, i, j);
         }
         linkmatrix[i][j] = 0;
      }
      /* If bounds fix variables to some value, initialize solvals with this value. */
      if( SCIPisLT(pricingprob, SCIPvarGetUbLocal(pricingprobvars[i]), 1.0) )
      {
         solvals[i] = 0.0;
         nfixedvars++;
      }
      else if( SCIPisGT(pricingprob, SCIPvarGetLbLocal(pricingprobvars[i]), 0.0) )
      {
         solvals[i] = 1.0;
         nfixedvars++;
      }
      else
         solvals[i] = -1.0; /* To later determine whether a variable was constrained */
   }

   SCIPdebugMessage( "Number of variables fixed by bound (before propagation): %d (of %d).\n", nfixedvars, npricingprobvars );

   for( i = 0; i < nconss; i++ )
   {
      consvarsfixedcount[i] = 0;       /* Initialize array to count the number of fixed vars per constraint. */
      couplingcoefindices[i] = -1;     /* Initialize array to save coupling coefficient if constraint is a coupling constraint. */
   }


   /* Determine constraint types for easier handling later on.
    * Also, it is checked for constraints that cannot be handled by this solver. */
   if( !determineCliquerConsTypes(pricingprob, constraints, linkedvars, vconsvars, markedconsindices, linkmatrix,
                                  couplingcoefindices, &nlinkedvars, nconss, &markedcount, cliquerconstypes) )
   {
      /* Encountered constraint that can not be handled. */
      *status = GCG_PRICINGSTATUS_NOTAPPLICABLE;
      setProblemNotApplicable(scip, probnr, solver->isnotapplicable);
      goto TERMINATE;
   }

   /* Cliquer may perform worse than other solvers (e.g. SCIP) on problems containing many clique inequalities:
    * ->    Thus, we do not apply the solver if the percentage of clique constraints exceeds a threshold parameter. */
   if( SCIPisLT(pricingprob, solver->cliqueconsthresh, 1.0) )
   {
      cliqueconscount = 0;
      for( i = 0; i < nconss; i++ )
      {
         if( cliquerconstypes[i] == LINEAR_CLIQUE )
            cliqueconscount++;
      }

      if( SCIPisGT(pricingprob, ((SCIP_Real)cliqueconscount / nconss), solver->cliqueconsthresh) )
      {
         SCIPdebugMessage("Exit: Clique-constraint percentage threshold exceeded, clique-cons perc.: %3.f\n",
                          ((SCIP_Real)cliqueconscount / nconss));
         *status = GCG_PRICINGSTATUS_NOTAPPLICABLE;
         goto TERMINATE;
      }
   }

   /* Propagate the already fixed variables to (potentially) get more fixed variables.
    * Also, builds coupling digraph to distribute objective coefficients of coupling variables. */
   if( (nlinkedvars > 0 || nfixedvars > 0 || solver->objcoefdistr > 0) &&
       !propagateVariablefixings(pricingprob, constraints, linkedvars, vconsvars, couplvars, solvals, couplingmatrix,
                                 linkmatrix, &ncouplvars, couplingcoefindices, consvarsfixedcount, nlinkedvars,
                                 nconss, &nfixedvars, cliquerconstypes) )
   {
      /* Variables are fixed in a conflicting way. -> problem is infeasible. */
      *status = GCG_PRICINGSTATUS_INFEASIBLE;
      goto TERMINATE;
   }

   SCIPdebugMessage( "Number of variables fixed before building the graph (after propagation): %d (of %d).\n",
                     nfixedvars, npricingprobvars );

   /* Check if all variables of the pricing problem are fixed. In this case, it is the only feasible solution. */
   /* No graph needs to be built, we just can build the corresponding column. */
   if( nfixedvars == npricingprobvars )
      goto CREATECOLUMN;

   /* Allocate memory needed for building the graph and creating a column. */
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &indsetvars, npricingprobvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &consvarsfixedtozerocount, nconss) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &aggrobjcoef, npricingprobvars) );
   ismemgraphallocated = TRUE;

   SCIPdebugMessage("nlinkedvars = %i , ncouplvars = %i , coefdistrheur = %i\n", nlinkedvars, ncouplvars, solver->objcoefdistr);

   /* If any variable is linked or coupled, aggrobjcoef array needs to be initialized. */
   if( nlinkedvars > 0 || ncouplvars > 0 )
   {
      for( i = 0; i < npricingprobvars; ++i )
         aggrobjcoef[SCIPvarGetProbindex(pricingprobvars[i])] = SCIPvarGetObj(pricingprobvars[i]);
   }

   /* Before adding nodes to the graph, aggregating the objective coefficients may be necessary if "same"-constraints exist. */
   if( nlinkedvars > 0 )
      aggregateObjCoef(pricingprob, linkmatrix, linkedvars, nlinkedvars, aggrobjcoef);

   /* If there are coupling or standard varbound constraints, it may be necessary to distribute objective coefficients. */
   if( solver->objcoefdistr > 0 && ncouplvars > 0 )
      distributeObjCoef(pricingprob, constraints, nconss, couplingcoefindices, cliquerconstypes, couplingmatrix,
                        couplvars, ncouplvars, linkedvars, nlinkedvars, linkmatrix, aggrobjcoef, solver->objcoefdistr,
                        solver->usemultiplicity);

   /* Now calculate scaling factor based on maximum aggregated objective coefficient value. */

   /* Cliquer library explicitly demands the node weights to be positive integers.
    * Additionally, the sum of node weights needs to be smaller than INT_MAX.
    * We restrict our scaling factor to always honor this constraint.
    */
   if( !areObjectivesIntegral(pricingprob, linkedvars, nlinkedvars, couplvars, ncouplvars, aggrobjcoef) )
      scalingfactor = scaleRelativeToMax(pricingprob, linkedvars, nlinkedvars, couplvars, ncouplvars, aggrobjcoef);
   else
      scalingfactor = 1.0;


   /* Count number of fixed variables and fixed-to-0 variables per constraint. */
   for( i = 0; i < nconss; i++ )
   {
      consvarsfixedtozerocount[i] = 0;

      /* Skip if there are no fixed variables. */
      if( nfixedvars <= 0)
         continue;

      /* Get variables of the constraint in dependence of the constraint handler */
      switch( cliquerconstypes[i] )
      {
         case LINEAR_IS:
         case LINEAR_IS_LIKE:
         case LINEAR_CLIQUE:
         case LINEAR_COUPLING_DECORATIVE:
         case LINEAR_COUPLING_CLIQUE:
            lconsvars = SCIPgetVarsLinear(pricingprob, constraints[i]);
            SCIPgetConsNVars(pricingprob, constraints[i], &nvars, &retcode);
            break;
         case VARBND_SAME:
         case VARBND_STD:
         case VARBND_IS:
            vconsvars[0] = SCIPgetVarVarbound(pricingprob, constraints[i]);
            vconsvars[1] = SCIPgetVbdvarVarbound(pricingprob, constraints[i]);
            lconsvars = vconsvars;
            nvars = 2;
            break;
      }

      /* Count variables fixed to 0. */
      for( j = 0; j < nvars; j++ )
      {
         if( solvals[SCIPvarGetProbindex(lconsvars[j])] == 0 )
            consvarsfixedtozerocount[i]++;
      }

      if( consvarsfixedtozerocount[i] == 0 )
      {
         /* Count of fixed variables is still correct. */
         continue;
      }
      if( consvarsfixedcount[i] < nvars && consvarsfixedtozerocount[i] == nvars )
      {
         /* All variables fixed to 0. */
         consvarsfixedcount[i] = consvarsfixedtozerocount[i];
      }
      else if( consvarsfixedcount[i] < nvars )
      {
         /* Need to recount the overall number of fixed variables. */
         consvarsfixedcount[i] = 0;
         for( j = 0; j < nvars; j++ )
         {
            if( solvals[SCIPvarGetProbindex(lconsvars[j])] >= 0 )
               consvarsfixedcount[i]++;
         }
      }
   }

   /* All links have to be established first before we can add nodes to the graph, else pairs (a,b) and (c,d) would be mapped to different nodes */
   /* if link (b,c) is present but later in the list. We have to run through the constraints again as the linked variables need to be assigned to nodes */
   /* in order for the rest of the logic to work out (node indices are fetched during runtime) */
   for( i = 0; i < markedcount; ++i )
   {
      /* Since we know that all marked constraints at this point are same constraints, we can just add them to the graph. */
      vconsvars[0] = SCIPgetVarVarbound(pricingprob, constraints[markedconsindices[i]]);
      if( SCIPisLT(pricingprob, getAggrObjCoef(vconsvars[0], nlinkedvars, ncouplvars, aggrobjcoef), 0) )
      {
         nodeindex0 = addVarToGraph(pricingprob, g, vconsvars[0], &indexcount, scalingfactor, indsetvars, linkmatrix,
                                    linkedvars, nlinkedvars, ncouplvars, aggrobjcoef);
      }
   }

   /* Main loop to check the nature of each constraint and manipulate the graph accordingly (add nodes, remove edges). */
   for( i = 0; i < nconss; ++i )
   {
      assert(constraints[i] != NULL);
      conshdlr = SCIPconsGetHdlr(constraints[i]);
      assert(conshdlr != NULL);

      SCIPgetConsNVars(pricingprob, constraints[i], &nvars, &retcode);

      /* The constraint may not be of type 'linear' */
      if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
      {
         /* If all variables are fixed, constraint can be skipped. */
         if( consvarsfixedcount[i] == nvars )
            continue;

         lconsvars = SCIPgetVarsLinear(pricingprob, constraints[i]);

         if((cliquerconstypes[i] == LINEAR_IS &&
              (SCIPisLT(pricingprob, getAggrObjCoef(lconsvars[0], nlinkedvars, ncouplvars, aggrobjcoef), 0) ||
               SCIPisLT(pricingprob, getAggrObjCoef(lconsvars[1], nlinkedvars, ncouplvars, aggrobjcoef), 0))) ||
             cliquerconstypes[i] == LINEAR_IS_LIKE )
         {
            /* One variable fixed to 0 (the other is not fixed): constraint is relaxed -> continue. */
            if( consvarsfixedcount[i] == 1 && consvarsfixedtozerocount[i] == 1 )
               continue;

            /* Add variables nodes to graph if they have a negative (aggregated) objective coefficient. */
            nodeindex0 = -1;
            if( SCIPisLT(pricingprob, getAggrObjCoef(lconsvars[0], nlinkedvars, ncouplvars, aggrobjcoef), 0) )
               nodeindex0 = addVarToGraph(pricingprob, g, lconsvars[0], &indexcount, scalingfactor, indsetvars, linkmatrix,
                                          linkedvars, nlinkedvars, ncouplvars, aggrobjcoef);

            nodeindex1 = -1;
            if( SCIPisLT(pricingprob, getAggrObjCoef(lconsvars[1], nlinkedvars, ncouplvars, aggrobjcoef), 0) )
               nodeindex1 = addVarToGraph(pricingprob, g, lconsvars[1], &indexcount, scalingfactor, indsetvars, linkmatrix,
                                          linkedvars, nlinkedvars, ncouplvars, aggrobjcoef);

            /* If both vairables nodes are added and an edge exists between the two in the graph: delete this edge. */
            if( nodeindex0 >= 0 && nodeindex1 >= 0 && GRAPH_IS_EDGE(g, nodeindex0, nodeindex1) )
            {
               GRAPH_DEL_EDGE(g, nodeindex0, nodeindex1);
            }
         }
         else
         {
            /* Cases in which constraint is relaxed through fixings. -> continue. */
            if( (cliquerconstypes[i] == LINEAR_CLIQUE && consvarsfixedcount[i] == nvars - 1 && consvarsfixedtozerocount[i] == nvars - 1) ||
                (cliquerconstypes[i] == LINEAR_COUPLING_CLIQUE &&
                 consvarsfixedcount[i] == nvars - 1 &&
                 consvarsfixedtozerocount[i] == nvars - 2 &&
                 solvals[SCIPvarGetProbindex(lconsvars[couplingcoefindices[i]])] == 1.0) )
               continue;

            /* If coupling constraint, add coupling var to graph and mark it in the solvals array. */
            if( (cliquerconstypes[i] == LINEAR_COUPLING_DECORATIVE && solvals[SCIPvarGetProbindex(lconsvars[couplingcoefindices[i]])] != 1.0) ||
                cliquerconstypes[i] == LINEAR_COUPLING_CLIQUE )
            {
               /* We cannot guarantee that there is no constraint of the form x+CouplingVar <= 1 */
               /* If the node is part of the maximum clique, it is safe to set it to one, so we simply add it to the graph */
               nodeindex0 = addVarToGraph(pricingprob, g, lconsvars[couplingcoefindices[i]], &indexcount, scalingfactor, indsetvars,
                                          linkmatrix, linkedvars, nlinkedvars, ncouplvars, aggrobjcoef);

               /* We additionally have to mark the variable to later set it to one */
               if( solvals[SCIPvarGetProbindex(lconsvars[couplingcoefindices[i]])] < 0.0 )
                   solvals[SCIPvarGetProbindex(lconsvars[couplingcoefindices[i]])] = -2.0;
            }

            /*
             * If (coupling-)clique constraint, try to add all (non-coupling) variables corresponding nodes to the graph
             * if the objective coefficient is < 0 and remove edges between these nodes (if added).
             */
            if( cliquerconstypes[i] == LINEAR_CLIQUE || cliquerconstypes[i] == LINEAR_COUPLING_CLIQUE )
            {
               /* Delete the edges between all the variables of the constraint (that are not the coupling variable).
                        This way, at most one can be part of the maximum clique */
               for( j = 0; j < nvars; ++j )
               {
                  /* We are only interested in vars potentially relevant for pricing (obj < 0) */
                  if((cliquerconstypes[i] != LINEAR_COUPLING_CLIQUE || j != couplingcoefindices[i]) &&
                     SCIPisLT(pricingprob, getAggrObjCoef(lconsvars[j], nlinkedvars, ncouplvars, aggrobjcoef), 0) &&
                      solvals[SCIPvarGetProbindex(lconsvars[j])] != 0.0 )
                  {
                     /* Determine nodeindex0 */
                     nodeindex0 = addVarToGraph(pricingprob, g, lconsvars[j], &indexcount, scalingfactor, indsetvars, linkmatrix,
                                                linkedvars, nlinkedvars, ncouplvars, aggrobjcoef);

                     /* Determine nodeindex1 */
                     for( k = j + 1; k < nvars; ++k )
                     {
                        if((cliquerconstypes[i] != LINEAR_COUPLING_CLIQUE || k != couplingcoefindices[i]) &&
                           SCIPisLT(pricingprob, getAggrObjCoef(lconsvars[k], nlinkedvars, ncouplvars, aggrobjcoef), 0) &&
                            solvals[SCIPvarGetProbindex(lconsvars[k])] != 0.0 )
                        {
                           nodeindex1 = addVarToGraph(pricingprob, g, lconsvars[k], &indexcount, scalingfactor, indsetvars, linkmatrix,
                                                      linkedvars, nlinkedvars, ncouplvars, aggrobjcoef);

                           if( (nodeindex0 != nodeindex1) && GRAPH_IS_EDGE(g, nodeindex0, nodeindex1) )
                           {
                              GRAPH_DEL_EDGE(g, nodeindex0, nodeindex1);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      /* Constraint may be of type varbound: lhs <= x + c*y <= rhs */
      else if( strcmp(SCIPconshdlrGetName(conshdlr), "varbound") == 0 )
      {
         /* If all variables are fixed, constraint can be skipped */
         if( consvarsfixedcount[i] == 2 )
            continue;

         vconsvars[0] = SCIPgetVarVarbound(pricingprob, constraints[i]);
         vconsvars[1] = SCIPgetVbdvarVarbound(pricingprob, constraints[i]);

         /* Form: x <= d*y with d >= 1 */
         if( cliquerconstypes[i] == VARBND_STD )
         {
            /* If x fixed to 0 or y fixed to 1 (and other variable not fixed): constraint relaxed -> continue. */
            if( consvarsfixedcount[i] == 1 &&
                (solvals[SCIPvarGetProbindex(vconsvars[0])] == 0.0 || solvals[SCIPvarGetProbindex(vconsvars[1])] == 1.0) )
               continue;

            /* if x may be relevant, add both x and y to graph */
            if( SCIPisLT(pricingprob, getAggrObjCoef(vconsvars[0], nlinkedvars, ncouplvars, aggrobjcoef), 0) )
            {
               nodeindex0 = addVarToGraph(pricingprob, g, vconsvars[0], &indexcount, scalingfactor, indsetvars, linkmatrix,
                                          linkedvars, nlinkedvars, ncouplvars, aggrobjcoef);
               nodeindex1 = addVarToGraph(pricingprob, g, vconsvars[1], &indexcount, scalingfactor, indsetvars, linkmatrix,
                                          linkedvars, nlinkedvars, ncouplvars, aggrobjcoef);
               /* It may be the case, that both the constraints x - y <= 0 and x + y <= 1 are part of the problem */
               /* Although rare, we later ensure that we do not set x to 1 while y is set to 0 */
               markedconsindices[markedcount] = i;
               ++markedcount;
            }
            /* If only y may be relevant, add only y to the graph */
            else if( SCIPisLT(pricingprob, getAggrObjCoef(vconsvars[1], nlinkedvars, ncouplvars, aggrobjcoef), 0) )
            {
               nodeindex1 = addVarToGraph(pricingprob, g, vconsvars[1], &indexcount, scalingfactor, indsetvars, linkmatrix,
                                          linkedvars, nlinkedvars, ncouplvars, aggrobjcoef);
            }
            /* If none of the nodes are relevant, force x to be zero, since the constraint would be violated if x = 1 and y = 0 */
            else
            {
               /* This logic might not always be correct. These vairables might be set to 1 in an optimal solution if, e.g.,
                * they appear as coupling variables in other constraints. Even if they are both not "relevant".
                */
               setLinkedSolvals(pricingprob, solvals, linkmatrix, linkedvars, nlinkedvars, vconsvars[0], 0.0);
            }
         }

         /* Form: x + y <= 1 */
         if( cliquerconstypes[i] == VARBND_IS )
         {
            /* If x fixed to 0 or y fixed to 0 (and other variable unfixed): constraint relaxed -> continue. */
            if( consvarsfixedcount[i] == 1 && consvarsfixedtozerocount[i] == 1 )
               continue;

            /* Preprocessing: Constraint is only relevant for pricing if one of the variables has an objective value < 0 */
            if( SCIPisLT(pricingprob, getAggrObjCoef(vconsvars[0], nlinkedvars, ncouplvars, aggrobjcoef), 0) ||
                SCIPisLT(pricingprob, getAggrObjCoef(vconsvars[1], nlinkedvars, ncouplvars, aggrobjcoef), 0) )
            {
               nodeindex0 = -1;
               if( SCIPisLT(pricingprob, getAggrObjCoef(vconsvars[0], nlinkedvars, ncouplvars, aggrobjcoef), 0) )
                  nodeindex0 = addVarToGraph(pricingprob, g, vconsvars[0], &indexcount, scalingfactor, indsetvars, linkmatrix,
                                             linkedvars, nlinkedvars, ncouplvars, aggrobjcoef);

               nodeindex1 = -1;
               if( SCIPisLT(pricingprob, getAggrObjCoef(vconsvars[1], nlinkedvars, ncouplvars, aggrobjcoef), 0) )
                  nodeindex1 = addVarToGraph(pricingprob, g, vconsvars[1], &indexcount, scalingfactor, indsetvars, linkmatrix,
                                             linkedvars, nlinkedvars, ncouplvars, aggrobjcoef);

               if( nodeindex0 >= 0 && nodeindex1 >= 0 && GRAPH_IS_EDGE(g, nodeindex0, nodeindex1) )
               {
                  GRAPH_DEL_EDGE(g, nodeindex0, nodeindex1);
               }
            }
         }
      }
   }


   /* Assert that the graph was built in a proper way */ 
   assert(graph_test(g, NULL));

   /* Determine number of edges for graph density calculation */
   nedges = 0;
   for( i = 0; i < g->n; i++ )
   {
      for( j = 0; j < g->n; j++ )
      {
         if( SET_CONTAINS_FAST(g->edges[i], j) )
         {
            nedges++;
         }
      }
   }
   nedges /= 2;

   density = (SCIP_Real)nedges / ((SCIP_Real)(g->n - 1) * (g->n) / 2);

   SCIPdebugMessage("Problem number: %i ; Tree depth: %i ; Graph size: %d ; Graph density: %g\n",
                    probnr, SCIPgetFocusDepth(scip), indexcount, density);

   /* Test if the node threshold is respected */
   if( SCIPisGT(pricingprob, indexcount, solver->nodelimit) )
   {
      SCIPdebugMessage("Exit: Node threshold exceeded, number of nodes: %i.\n", indexcount);
      *status = GCG_PRICINGSTATUS_NOTAPPLICABLE;
      goto TERMINATE;
   }

   /* Only apply density / linear cutoff if density start threshold is exceeded. */
   if( SCIPisGT(pricingprob, indexcount, solver->densitystart) )
   {
      /* Test if the density criteria is met */
      if( SCIPisGT(pricingprob, density, solver->density) )
      {
         SCIPdebugMessage("Exit: Density criteria not met, density: %g.\n", density);
         *status = GCG_PRICINGSTATUS_NOTAPPLICABLE;
         goto TERMINATE;
      }

      /* Next, check if linear cutoff is activated, If yes, check if linear cutoff equation is met.
       * If (n)odes and (d)ensity have values n > m*d + b (with slope m and intercept b), solver is not applied.
       * Default values currently are: m = -1980, b = 1900
       */
      if( solver->uselincutoff &&
          SCIPisGT(pricingprob, indexcount, solver->lincutoffslope * density + solver->lincutoffintercept) )
      {
         SCIPdebugMessage("Exit: Linear threshold n <= m*d + b exceeded (i.e.: %i > %.1f * %.2f + %.1f).\n",
                          indexcount, solver->lincutoffslope, density, solver->lincutoffintercept);
         *status = GCG_PRICINGSTATUS_NOTAPPLICABLE;
         goto TERMINATE;
      }
   }

   ASSERT( indexcount <= npricingprobvars );

   /* indexcount now holds the actual number of unique IS variables, thus we truncate the graph */
   if( indexcount > 0 )
   {
      graph_resize(g, indexcount);
   }

   /* Clean up the graph. If a variable's solval has been set to 0, it should not be part of the max clique */
   /* We enforce this by isolating the node and setting its weight to 1 as nodes cannot be deleted */
   for( i = 0; i < npricingprobvars; ++i )
   {
      if( solvals[SCIPvarGetProbindex(pricingprobvars[i])] == 0 )
      {
         nodeindex0 = getLinkedNodeIndex(pricingprob, pricingprobvars[i], indsetvars, indexcount, linkmatrix, linkedvars, nlinkedvars);
         /* The var is part of the graph if its index is unequal to -1 */
         if( nodeindex0 != -1 )
         {
            for( j = 0; j < indexcount; ++j )
            {
               if( GRAPH_IS_EDGE(g, nodeindex0, j) )
               {
                  GRAPH_DEL_EDGE(g, nodeindex0, j);
               }
            }
            g->weights[nodeindex0] = 1;
         }
      }
   }

   /* Set cliquer options */
   cl_opts.reorder_function = reorder_by_default; /* default: reorder_by_default */
   cl_opts.reorder_map = NULL;
   cl_opts.time_function = NULL; /* default: clique_print_time */
   cl_opts.output = NULL;
   cl_opts.user_function = NULL;
   cl_opts.user_data = NULL;
   cl_opts.clique_list = NULL;
   cl_opts.clique_list_length = 0;

   /* Find maximum weight clique using the cliquer library */
   clique = clique_find_single(g, 0, 0, FALSE, &cl_opts);

   /* Set all members of the maximum clique with objective coefficient < 0 to 1 */
   for( i = 0; i < indexcount; ++i )
   {
      /* Coupling variables were pre-set to -2.0, if they are part of the maximum clique, we enable them. 
       * If we have already set a variable to 0, this was intended and should not be reverted.
       *
       * NOTE: As long as coupling variables may have positive cost but have cost of 1 in graph, the solver is of heuristic nature.
       *       The '-2.0'-marked coupling variables are set to 1, even if they render the solution be of non-negative reduced cost.
       *       The max-weighted-clique solver would need to support also negative costs to heal this.
       */
      if( SET_CONTAINS(clique, i) &&
          (SCIPisLT(pricingprob, getAggrObjCoef(indsetvars[i], nlinkedvars, ncouplvars, aggrobjcoef), 0) ||
           solvals[SCIPvarGetProbindex(indsetvars[i])] == -2.0) &&
          solvals[SCIPvarGetProbindex(indsetvars[i])] != 0.0 )
      {
         /* Set all linked variables, if any */
         setLinkedSolvals(pricingprob, solvals, linkmatrix, linkedvars, nlinkedvars, indsetvars[i], 1.0);
      }
      else
      {
         /* We may have set some variables manually already, e.g. coupling variables */
         if( solvals[SCIPvarGetProbindex(indsetvars[i])] != 1.0)
         {
            setLinkedSolvals(pricingprob, solvals, linkmatrix, linkedvars, nlinkedvars, indsetvars[i], 0.0);
         }
      }
   }

   for( i = 0; i < markedcount; ++i )
   {
      vconsvars[0] = SCIPgetVarVarbound(pricingprob, constraints[markedconsindices[i]]);
      vconsvars[1] = SCIPgetVbdvarVarbound(pricingprob, constraints[markedconsindices[i]]);

      /* Handle the case of marked inequality constraints of type x - y <= 0 in combination with x + y <= 1 -Constraints */
      if( cliquerconstypes[markedconsindices[i]] == VARBND_STD )
      {
         /* Check if a violating assignment was made and correct it */
         if( (solvals[SCIPvarGetProbindex(vconsvars[0])] == 1) && (solvals[SCIPvarGetProbindex(vconsvars[1])] == 0) )
         {
            setLinkedSolvals(pricingprob, solvals, linkmatrix, linkedvars, nlinkedvars, vconsvars[0], 0.0);
         }
      }

      /* Handle the case that there are still solvals of equality constraints that do not agree.
       * This may occur if one is unset (solval:-1) and the other one is already set (solval 0 or 1)
       */
      if( solvals[SCIPvarGetProbindex(vconsvars[0])] != solvals[SCIPvarGetProbindex(vconsvars[1])] &&
          cliquerconstypes[markedconsindices[i]] == VARBND_SAME )
      {
         if( solvals[SCIPvarGetProbindex(vconsvars[0])] == 0 || solvals[SCIPvarGetProbindex(vconsvars[1])] == 0 )
         {
            setLinkedSolvals(pricingprob, solvals, linkmatrix, linkedvars, nlinkedvars, vconsvars[0], 0.0);
         }
         else
         {
            /* One or both of the vars are unset and the other one, if not -1, is forced to be 1, thus we can set both to 1 */
            setLinkedSolvals(pricingprob, solvals, linkmatrix, linkedvars, nlinkedvars, vconsvars[0], 1.0);
         }
      }
   }

   /* There may be variables left which are unconstrained. We set these to 1 manually if they have an objective value < 0*/
   for( i = 0; i < npricingprobvars; ++i )
   {
      if( solvals[i] == -1.0 )
      {
         if( SCIPisLT(pricingprob, getAggrObjCoef(pricingprobvars[i], nlinkedvars, ncouplvars, aggrobjcoef), 0) )
         {
            solvals[i] = 1.0;
         }
         else
         {
            solvals[i] = 0.0;
         }
      }
   }

   /* Create a column corresponding to our clique result */
 CREATECOLUMN:
   SCIP_CALL( GCGcreateGcgCol(pricingprob, &col, probnr, pricingprobvars, solvals, npricingprobvars, FALSE, SCIPinfinity(pricingprob)) );
   SCIP_CALL( GCGpricerAddCol(gcg, col) );
   *status = GCG_PRICINGSTATUS_UNKNOWN;
   if( indexcount > 0 )
      set_free(clique); /* clique can only be freed if non-empty */

 TERMINATE:
   if( ismemgraphallocated )
   {
      SCIPfreeBufferArray(pricingprob, &indsetvars);
      SCIPfreeBufferArray(pricingprob, &consvarsfixedtozerocount);
      SCIPfreeBufferArray(pricingprob, &aggrobjcoef);
   }
   SCIPfreeBufferArray(pricingprob, &couplvars);
   for( i = 0; i < npricingprobvars; ++i )
   {
      SCIPfreeBufferArray(pricingprob, &couplingmatrix[i]);
   }
   SCIPfreeBufferArray(pricingprob, &couplingmatrix);
   SCIPfreeBufferArray(pricingprob, &vconsvars);
   SCIPfreeBufferArray(pricingprob, &consvarsfixedcount);
   SCIPfreeBufferArray(pricingprob, &couplingcoefindices);
   for( i = 0; i < npricingprobvars; ++i )
   {
      SCIPfreeBufferArray(pricingprob, &linkmatrix[i]);
   }
   SCIPfreeBufferArray(pricingprob, &linkmatrix);
   SCIPfreeBufferArray(pricingprob, &linkedvars);
   SCIPfreeBufferArray(pricingprob, &solvals);
   SCIPfreeBufferArray(pricingprob, &markedconsindices);
   SCIPfreeBufferArray(pricingprob, &cliquerconstypes);
   graph_free(g);

   return SCIP_OKAY;
}

/*
 * Callback methods for pricing problem solver
 */

/** destructor of pricing solver to free user data (called when SCIP is exiting) */
static
GCG_DECL_SOLVERFREE(solverFreeCliquer)
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
GCG_DECL_SOLVERINITSOL(solverInitsolCliquer)
{
   GCG_SOLVERDATA* solverdata;
   int npricingprobs;

   assert(scip != NULL);
   assert(solver != NULL);

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   /* allocate and initialize isnotapplicable array */
   npricingprobs = GCGgetNPricingprobs(GCGmasterGetOrigprob(scip));
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &solverdata->isnotapplicable, npricingprobs) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of pricing solver (called before branch and bound process data is freed) */
static
GCG_DECL_SOLVEREXITSOL(solverExitsolCliquer)
{
   GCG_SOLVERDATA* solverdata;
   int npricingprobs;

   assert(scip != NULL);
   assert(solver != NULL);

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   /* free isnotapplicable array */
   npricingprobs = GCGgetNPricingprobs(GCGmasterGetOrigprob(scip));
   SCIPfreeBlockMemoryArray(scip, &solverdata->isnotapplicable, npricingprobs);

   return SCIP_OKAY;
}

#define solverInitCliquer NULL
#define solverExitCliquer NULL
#define solverUpdateCliquer NULL
#define solverSolveCliquer NULL

/** heuristic solving method of independent set solver */
static
GCG_DECL_SOLVERSOLVEHEUR(solverSolveHeurCliquer)
{  /*lint --e{715}*/
   GCG_SOLVERDATA* solverdata;

   solverdata = GCGsolverGetData(solver);
   assert(solverdata != NULL);

   /* solve the independent set problem approximately */
   SCIP_CALL( solveCliquer(FALSE, gcg, pricingprob, solverdata, probnr, lowerbound, status) );

   return SCIP_OKAY;
}

/** creates the cliquer solver for pricing problems and includes it in GCG */
SCIP_RETCODE GCGincludeSolverCliquer(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_SOLVERDATA* solverdata;
   SCIP* origprob;

   origprob = GCGgetOrigprob(gcg);

   SCIP_CALL( SCIPallocMemory(GCGgetDwMasterprob(gcg), &solverdata) );

   SCIP_CALL( GCGpricerIncludeSolver(gcg, SOLVER_NAME, SOLVER_DESC, SOLVER_PRIORITY,
         SOLVER_HEURENABLED, SOLVER_EXACTENABLED,
         solverUpdateCliquer, solverSolveCliquer, solverSolveHeurCliquer,
         solverFreeCliquer, solverInitCliquer, solverExitCliquer,
         solverInitsolCliquer, solverExitsolCliquer, solverdata) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/cliquer/density",
          "graph density threshold below which to use solver",
          &solverdata->density, TRUE, DEFAULT_DENSITY, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricingsolver/cliquer/densitystart",
          "graph node threshold above which to apply density threshold / linear cutoff (below not applied)",
          &solverdata->densitystart, TRUE, DEFAULT_DENSITYSTART, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/cliquer/maxcliqueconsperc",
          "threshold for share of clique constraints in pricing problem below which to use solver (disabled = 1.0)",
          &solverdata->cliqueconsthresh, TRUE, DEFAULT_CLIQUECONSTHRESH, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricingsolver/cliquer/nodelimit",
          "graph node threshold below which to use solver",
          &solverdata->nodelimit, TRUE, DEFAULT_NODELIMIT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricingsolver/cliquer/objcoefdistr",
          "distribution of objective coefficients of coupling variables (disabled = 0, natural share = 1, "
          "MIS-based = 2, uniform = 3)",
          &solverdata->objcoefdistr, TRUE, DEFAULT_OBJCOEFDISTR, 0, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "pricingsolver/cliquer/usemultiplicity",
          "should the usage of multiplicity of linked variables be used to weight distributed coefficients be enabled? "
          "(only useful with objcoefdistr != 0)",
          &solverdata->usemultiplicity, TRUE, DEFAULT_USEMULTIPL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origprob, "pricingsolver/cliquer/lincutoff/enable",
          "should linear cutoff (n > m*d + b) for usage of solver <cliquer>, based on graph (d)ensity and (n)odes, "
          "be enabled?",
          &solverdata->uselincutoff, FALSE, DEFAULT_USELINCUTOFF, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/cliquer/lincutoff/slope",
          "slope m in the linear cutoff formula (n > m*d + b), with (d)ensity and (n)odes",
          &solverdata->lincutoffslope, TRUE, DEFAULT_SLOPE, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricingsolver/cliquer/lincutoff/intercept",
          "intercept b in the linear cutoff formula (n > m*d + b), with (d)ensity and (n)odes",
          &solverdata->lincutoffintercept, TRUE, DEFAULT_INTERCEPT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
