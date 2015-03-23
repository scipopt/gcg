/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2014 Operations Research, RWTH Aachen University       */
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

/**@file   heur_setcover.c
 * @brief  setcover primal heuristic
 * @author Tobias Oelschlaegel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIP_DEBUG

#include <assert.h>
#include <string.h>
#include "gcg.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "scip/clock.h"
#include "scip/cons_linear.h"
#include "heur_setcover.h"


#define HEUR_NAME             "setcover"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         '?'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define SCP_CORE_TENT_SIZE    5
#define SCP_LAMBDA_ADJUSTMENTS
#define SCP_LAMBDA_P          50
#define SCP_STOP_CRIT_ITER    300
#define SCP_STOP_CRIT_DIFF    1.0
#define SCP_STOP_CRIT_PER     0.99
#define SCP_PI_MIN            0.3
#define SCP_PI_ALPHA          1.1
#define SCP_BETA              1.0

/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Bool firststart;
   SCIP_Bool issetcover;
};

typedef struct
{
   SCIP_Real *costs;
   /* SCIP_Real costs_fixed; */
   SCIP_Bool *fixedvars;
   SCIP_Bool *covered;
   int numfixed;
   int maxrows;
   int maxcolumns;
} SCP_Lagrange;

typedef struct
{
   SCIP_HASHTABLE *corevariables;
   SCIP_HASHMAP *mapvariables; /* maps variables-indices to [0, nvariables) in array 'variables' */
   SCIP_VAR **variables;
   int *nvarconstraints; /* counts how many times a variable occurs in active constraints */
   int nvariables; /* total number of variables that can be found in active constraints */
   SCIP_Bool columnsavailable;
   int **columns;
   int nconstraints; /* total number of constraints (including inactive ones) */
   int nactiveconstraints; /* total number of active constraints */
   int maxconstraintvariables; /* greatest number of variables some constraint contains */
   SCIP_CONS **constraints;
} SCP_Core;

typedef struct
{
   SCIP_Real *u;
   SCIP_Real *x_lagrange;
   SCIP_Real *x_lagrange_global;
   SCIP_Real *x_greedy;
   SCIP_Real *subgradient;
   SCIP_Real *lagrangian_costs;
   SCIP_Real ub_greedy;
   SCIP_Real lb_lagrange;
   SCIP_Real lb_lagrange_global;
   SCIP_Real *delta;
   int *delta_index;
} SCP_Lagrange_Sol;
/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static SCIP_DECL_HASHGETKEY(hashGetKeyVar)
{
   return elem;
}


static SCIP_DECL_HASHKEYEQ(hashKeyEqVar)
{
   SCIP_VAR *var1 = (SCIP_VAR *) key1;
   SCIP_VAR *var2 = (SCIP_VAR *) key2;

   if(SCIPvarGetIndex(var1) == SCIPvarGetIndex(var2))
      return TRUE;
   else
      return FALSE;
}

static SCIP_DECL_HASHKEYVAL(hashKeyValVar)
{
   SCIP_VAR *var = (SCIP_VAR *) key;
   return (unsigned int) SCIPvarGetIndex(var);
}

static SCIP_RETCODE allocateMemoryForSolution(SCIP *scip, SCP_Lagrange *lag, SCP_Lagrange_Sol *mult)
{
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->u, lag->maxrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->subgradient, lag->maxrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->x_lagrange, lag->maxcolumns) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->x_lagrange_global, lag->maxcolumns) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->lagrangian_costs, lag->maxcolumns) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->x_greedy, lag->maxcolumns) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->delta, lag->maxcolumns) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->delta_index, lag->maxcolumns) );
   return SCIP_OKAY;
}

static SCIP_RETCODE freeMemoryForSolution(SCIP *scip, SCP_Lagrange_Sol *mult)
{
   SCIPfreeBufferArray(scip, &mult->u);
   SCIPfreeBufferArray(scip, &mult->subgradient);
   SCIPfreeBufferArray(scip, &mult->x_lagrange);
   SCIPfreeBufferArray(scip, &mult->x_lagrange_global);
   SCIPfreeBufferArray(scip, &mult->lagrangian_costs);
   SCIPfreeBufferArray(scip, &mult->x_greedy);
   SCIPfreeBufferArray(scip, &mult->delta);
   SCIPfreeBufferArray(scip, &mult->delta_index);
   return SCIP_OKAY;
}

static SCIP_RETCODE initTentativeCore(SCIP *scip, SCP_Core *core)
{
   SCIP_VAR **vars;
   SCIP_Bool success;
   int nvars;
   int i, j;

   assert(scip != NULL);
   assert(core != NULL);

   SCIP_CALL( SCIPhashtableCreate(&core->corevariables, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );
   SCIP_CALL( SCIPhashmapCreate(&core->mapvariables, SCIPblkmem(scip), SCIPcalcHashtableSize(SCIPgetNVars(scip))) );

   core->columnsavailable = FALSE;
   core->columns = NULL;
   core->nvariables = SCIPgetNVars(scip);
   core->variables = SCIPgetVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &core->nvarconstraints, core->nvariables) );
   for(i = 0; i< core->nvariables; i++)
      core->nvarconstraints[i] = 0;

   /* construct mapping of variable-indices to array 'variables' */
   for(i = 0; i < core->nvariables; i++)
   {
      int varidx = SCIPvarGetIndex(core->variables[i]);
      SCIP_CALL( SCIPhashmapInsert(core->mapvariables, (void *) (size_t) varidx, (void *) (size_t) i) );
   }

   core->nactiveconstraints = 0;
   core->maxconstraintvariables = 0;
   core->nconstraints = SCIPgetNConss(scip);
   core->constraints = SCIPgetConss(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->nvariables) );

   for(i = 0; i < core->nconstraints; i++)
   {
      if(SCIPconsIsActive(core->constraints[i]) == FALSE)
         continue;

      /* get all variables that are part of this constraint */
      SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[i], &nvars, &success) );
      if(success == FALSE)
         continue;

      SCIP_CALL( SCIPgetConsVars(scip, core->constraints[i], vars, core->nvariables, &success) );
      if(success == FALSE)
         continue;

      if(nvars > core->maxconstraintvariables)
         core->maxconstraintvariables = nvars;

      for(j = 0; j < nvars; j++)
      {
         int varidx = SCIPvarGetIndex(vars[j]);
         int varpos = (int) (size_t) SCIPhashmapGetImage(core->mapvariables, (void *) (size_t) varidx);

         if(i < SCP_CORE_TENT_SIZE)
         {
            /* add this variable to the core if it's not already in there */
            if(SCIPhashtableExists(core->corevariables, core->variables[varpos]) == FALSE)
               SCIPhashtableSafeInsert(core->corevariables, core->variables[varpos]);
         }

         /* increase the number of constraints this variable is part of */
         core->nvarconstraints[varpos]++;
      }

      core->nactiveconstraints++;
   }

   SCIPfreeBufferArray(scip, &vars);
   SCIPdebugMessage("%lli variables in the tentative core\n", SCIPhashtableGetNElements(core->corevariables));

   return SCIP_OKAY;
}

static SCIP_RETCODE computeCoreColumns(SCIP *scip, SCP_Core *core)
{
   SCIP_Bool success;
   SCIP_VAR **vars;
   int i, j, k;
   int nvars;

   assert(SCIP != NULL);
   assert(core != NULL);

   /* don't compute again if already computed */
   if(core->columnsavailable)
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &core->columns, core->nvariables) );
   for(i = 0; i < core->nvariables; i++)
   {
      core->columns[i] = NULL;

      /* check if variable is a core variable */
      if(SCIPhashtableExists(core->corevariables, core->variables[i]) == FALSE)
         continue;

      /* only allocate memory of it is part of any constraints at all (this should always be the case!) */
      if(core->nvarconstraints[i] == 0)
         continue;

      SCIP_CALL( SCIPallocBufferArray(scip, &core->columns[i], core->nvarconstraints[i]) );

      for(j = 0; j < core->nvarconstraints[i]; j++)
         core->columns[i][j] = -1;

      /* SCIPdebugMessage("allocated %i slots for variable %i\n", core->nvarconstraints[i], SCIPvarGetIndex(core->variables[i])); */
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );

   for(i = 0; i < core->nconstraints; i++)
   {
      if(SCIPconsIsActive(core->constraints[i]) == FALSE)
         continue;

      /* get all variables that are part of this constraint */
      SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[i], &nvars, &success) );
      if(success == FALSE)
         continue;

      SCIP_CALL( SCIPgetConsVars(scip, core->constraints[i], vars, core->maxconstraintvariables, &success) );
      if(success == FALSE)
         continue;

      for(j = 0; j < nvars; j++)
      {
         int varidx = SCIPvarGetIndex(vars[j]);
         int varpos = (int) (size_t) SCIPhashmapGetImage(core->mapvariables, (void *) (size_t) varidx);

         if(SCIPhashtableExists(core->corevariables, core->variables[varpos]) == FALSE)
            continue;

         /* add this constraint to the column of the variable */
         for(k = 0; k < core->nvarconstraints[varpos]; k++)
         {
            if(core->columns[varpos][k] == -1)
            {
               core->columns[varpos][k] = i;
               break;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &vars);
   core->columnsavailable = TRUE;
   return SCIP_OKAY;
}

static SCIP_RETCODE freeCore(SCIP *scip, SCP_Core *core)
{
   assert(core != NULL);
   assert(core->corevariables != NULL);
   assert(core->mapvariables != NULL);
   assert(core->nvarconstraints != NULL);

   SCIPhashtableFree(&core->corevariables);
   SCIPhashmapFree(&core->mapvariables);
   SCIPfreeBufferArray(scip, &core->nvarconstraints);

   if(core->columnsavailable == TRUE)
   {
      int i;
      for(i = 0; i < core->nvariables; i++)
      {
         if(core->columns[i] != NULL)
            SCIPfreeBufferArray(scip, &core->columns[i]);
      }
      SCIPfreeBufferArray(scip, &core->columns);
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE copySolution(SCP_Lagrange *lag, SCP_Lagrange_Sol *dest, SCP_Lagrange_Sol *source)
{
   int i;
   for(i = 0; i < lag->maxcolumns; i++)
   {
      dest->x_greedy[i] = source->x_greedy[i];
      dest->x_lagrange[i] = source->x_lagrange[i];
      dest->x_lagrange_global[i] = source->x_lagrange_global[i];
      dest->lagrangian_costs[i] = source->lagrangian_costs[i];
      dest->delta[i] = source->delta[i];
   }

   for(i = 0; i < lag->maxrows; i++)
   {
      dest->u[i] = source->u[i];
      dest->subgradient[i] = source->subgradient[i];
   }

   dest->ub_greedy = source->ub_greedy;
   dest->lb_lagrange = source->lb_lagrange;

   return SCIP_OKAY;
}

static SCIP_RETCODE markRowsCoveredByFixedVariables(SCIP *scip, SCP_Lagrange *lag)
{
   SCIP_COL **columns;
   int ncolumns;
   int i, j;

   SCIP_CALL( SCIPgetLPColsData(scip, &columns, &ncolumns) );
   for(i = 0; i < ncolumns; i++)
   {
      if(lag->fixedvars[SCIPcolGetIndex(columns[i])] == TRUE)
      {
         /* mark all rows of this column as covered */
         SCIP_ROW **rows = SCIPcolGetRows(columns[i]);
         int nrows = SCIPcolGetNNonz(columns[i]);

         for(j = 0; j < nrows; j++)
            lag->covered[SCIProwGetIndex(rows[j])] = TRUE;
      }
   }
   return SCIP_OKAY;
}

static int countUncoveredRows(SCP_Lagrange *lag, SCIP_COL *col)
{
   int count = 0;
   int i;
   int nrows = SCIPcolGetNNonz(col);
   SCIP_ROW **rows = SCIPcolGetRows(col);

   for(i = 0; i < nrows; i++)
   {
      SCIP_ROW *row = rows[i];
      if(lag->covered[SCIProwGetIndex(row)] == FALSE)
         count++;
   }

   return count;
}

static SCIP_RETCODE checkSetCover(SCIP *scip, SCP_Lagrange *lag, SCP_Lagrange_Sol *mult, SCIP_Bool *issetcover)
{
   SCIP_Longint *num_covered;
   SCIP_ROW **rows;
   int nrows;
   int i, j;

   /* TODO: use constraints instead of the LP, just need to iterate over all constraints and check that some of their variables are in the set cover */
   SCIP_CALL( SCIPallocBufferArray(scip, &num_covered, lag->maxrows) );
   for(i = 0; i < lag->maxrows; i++)
      num_covered[i] = 0;

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   for(i = 0; i < nrows; i++)
   {
      int rowidx = SCIProwGetIndex(rows[i]);
      if(lag->covered[rowidx] == FALSE)
      {
         SCIP_COL **covered_by = SCIProwGetCols(rows[i]);
         int numcols = SCIProwGetNNonz(rows[i]);

         for(j = 0; j < numcols; j++)
         {
            int colidx = SCIPcolGetIndex(covered_by[j]);
            if(mult->x_greedy[colidx] == TRUE)
               num_covered[rowidx]++;
         }
      }
   }

   *issetcover = TRUE;
   for(i = 0; (i < lag->maxrows) && *issetcover; i++)
   {
      if((lag->covered[i] == FALSE) && (num_covered[i] == 0))
         *issetcover = FALSE;
   }

   SCIPfreeBufferArray(scip, &num_covered);
   return SCIP_OKAY;
}

static SCIP_RETCODE greedySetCover(SCIP *scip, SCP_Lagrange *lag, SCP_Lagrange_Sol *mult)
{
   SCIP_COL **columns;
   int ncolumns;
   SCIP_Bool *uncovered;
   int *mu;
   SCIP_Real *gamma;
   SCIP_Longint *cover;
   int cover_size = 0;
   int i, j, k;
   int numuncovered = 0;
   SCIP_Bool issetcover = FALSE;

   mult->ub_greedy = 0.0;

   SCIP_CALL( SCIPallocBufferArray(scip, &uncovered, lag->maxrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mu, lag->maxcolumns) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cover, lag->maxcolumns) );
   SCIP_CALL( SCIPallocBufferArray(scip, &gamma, lag->maxcolumns) );

   for(i = 0; i < lag->maxrows; i++)
   {
      if(lag->covered[i] == TRUE)
      {
         uncovered[i] = FALSE;
      }
      else
      {
         uncovered[i] = TRUE;
         numuncovered++;
      }
   }

   for(i = 0; i < lag->maxcolumns; i++)
   {
      mult->x_greedy[i] = FALSE;
      mu[i] = 0.0;
      gamma[i] = lag->costs[i];
   }

   SCIP_CALL( SCIPgetLPColsData(scip, &columns, &ncolumns) );
   for(i = 0; i < ncolumns; i++)
   {
      int colidx = SCIPcolGetIndex(columns[i]);
      if(lag->fixedvars[colidx] == FALSE)
      {
         SCIP_ROW **rows = SCIPcolGetRows(columns[i]);
         int nrows = SCIPcolGetNNonz(columns[i]);

         for(j = 0; j < nrows; j++)
         {
            int rowidx = SCIProwGetIndex(rows[j]);
            if(uncovered[rowidx] == TRUE)
            {
               mu[colidx] = mu[colidx] + 1;
               gamma[colidx] = gamma[colidx] - mult->u[rowidx];
            }
         }
      }
   }

   while(numuncovered > 0)
   {
      SCIP_Bool found = FALSE;
      int minvariable = 0;
      SCIP_Real minscore = 0.0;

      /* TODO: replace this by a priority queue */
      for(i = 0; i < lag->maxcolumns; i++)
      {
         if((lag->fixedvars[i] == FALSE) && (mult->x_greedy[i] == FALSE) && (mu[i] > 0))
         {
            SCIP_Real score = (gamma[i] > 0) ? (gamma[i] / mu[i]) : (gamma[i] * mu[i]);
            if((found == FALSE) || (score < minscore))
            {
               found = TRUE;
               minvariable = i;
               minscore = score;
            }
         }
      }

      /* take 'minvariable' into the solution */
      mult->x_greedy[minvariable] = TRUE;
      cover[cover_size++] = minvariable;

      mult->ub_greedy = mult->ub_greedy + lag->costs[minvariable];

      SCIP_CALL( SCIPgetLPColsData(scip, &columns, &ncolumns) );
      for(i = 0; i < ncolumns; i++)
      {
         int colidx = SCIPcolGetIndex(columns[i]);

         /* we have found the column of this variable; there is probably a better way to do this than iterating through all columns */
         if(colidx == minvariable)
         {
            SCIP_ROW **rows = SCIPcolGetRows(columns[i]);
            int nrows = SCIPcolGetNNonz(columns[i]);

            for(j = 0; j < nrows; j++)
            {
               int rowidx = SCIProwGetIndex(rows[j]);
               if(uncovered[rowidx] == TRUE)
               {
                  SCIP_COL **covcolumns = SCIProwGetCols(rows[j]);
                  int ncovcolumns = SCIProwGetNNonz(rows[j]);

                  uncovered[rowidx] = FALSE;
                  numuncovered--;

                  /* decrease mu and gamma for all columns covering this row */
                  for(k = 0; k < ncovcolumns; k++)
                  {
                     int covidx = SCIPcolGetIndex(covcolumns[k]);
                     if((lag->fixedvars[covidx] == FALSE) && (mult->x_greedy[covidx] == FALSE))
                     {
                        mu[covidx] = mu[covidx] - 1;
                        gamma[covidx] = gamma[covidx] + mult->u[rowidx];
                     }
                  }
               }
            }
         }
      }
   }

   /* remove redundant columns from the solution. TODO: exact solution if set cover is very small */
   for(i = 0; i < cover_size; i++)
   {
      /* remove cover[i] from the set cover and test whether it is still a set cover. Leave it out if it is */
      assert(mult->x_greedy[cover[i]] == TRUE);
      mult->x_greedy[cover[i]] = FALSE;

      /* TODO: count how many times a set is covered. then we can easily check whether this is still a set cover */
      SCIP_CALL( checkSetCover(scip, lag, mult, &issetcover) );
      if(issetcover == TRUE)
         mult->ub_greedy = mult->ub_greedy - lag->costs[cover[i]];
      else
         mult->x_greedy[cover[i]] = TRUE;
   }

   /* verify that the solution actually is a set cover! */
   SCIP_CALL( checkSetCover(scip, lag, mult, &issetcover) );
   if(issetcover == FALSE)
   {
      SCIPdebugMessage("greedy algorithm did not compute a valid set cover\n");
      SCIPABORT();
   }

   SCIPfreeBufferArray(scip, &cover);
   SCIPfreeBufferArray(scip, &uncovered);
   SCIPfreeBufferArray(scip, &mu);
   SCIPfreeBufferArray(scip, &gamma);

   return SCIP_OKAY;
}

/* computes lagrangian costs for ALL columns */
static SCIP_RETCODE computeLagrangianCosts(SCIP *scip, SCP_Lagrange *lag, SCP_Lagrange_Sol *mult)
{
   SCIP_COL **columns;
   SCIP_ROW **rows;
   int ncolumns;
   int nrows;
   int i, j;

   SCIP_CALL( SCIPgetLPColsData(scip, &columns, &ncolumns) );
   for(i = 0; i < ncolumns; i++)
   {
      int colidx = SCIPcolGetIndex(columns[i]);
      SCIP_Real lagrangian_costs = lag->costs[colidx];

      rows = SCIPcolGetRows(columns[i]);
      nrows = SCIPcolGetNNonz(columns[i]);

      /* compute c_i - sum_{j \in I_i} u_j */
      for(j = 0; j < nrows; j++)
      {
         int rowidx = SCIProwGetIndex(rows[j]);
         lagrangian_costs = lagrangian_costs - mult->u[rowidx];
      }

      mult->lagrangian_costs[colidx] = lagrangian_costs;
   }

   return SCIP_OKAY;
}
/* computes an optimal solution to the lagrangian relaxation, see formulae (4), (5) in the paper */
static SCIP_RETCODE computeOptimalSolution(SCIP *scip, SCP_Lagrange *lag, SCP_Lagrange_Sol *mult)
{
   SCIP_ROW **rows;
   int nrows;
   int i;

   mult->lb_lagrange = 0.0;
   mult->lb_lagrange_global = 0.0;

   SCIP_CALL( computeLagrangianCosts(scip, lag, mult) );

   for(i = 0; i < lag->maxcolumns; i++)
   {
      mult->x_lagrange[i] = FALSE;
      mult->x_lagrange_global[i] = FALSE;

      /* x[i] = 1 if costs < 0, 0 otherwise */
      if(mult->lagrangian_costs[i] < 0.0)
      {
         mult->x_lagrange_global[i] = TRUE;
         mult->lb_lagrange_global = mult->lb_lagrange_global + mult->lagrangian_costs[i];

         if(lag->fixedvars[i] == FALSE)
         {
            mult->x_lagrange[i] = TRUE;
            mult->lb_lagrange = mult->lb_lagrange + mult->lagrangian_costs[i];
         }
      }
   }

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   for(i = 0; i < nrows; i++)
   {
      int rowidx = SCIProwGetIndex(rows[i]);

      mult->lb_lagrange_global = mult->lb_lagrange_global + mult->u[rowidx];

      if(lag->covered[rowidx] == FALSE)
         mult->lb_lagrange = mult->lb_lagrange + mult->u[rowidx];
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE computeSubgradient(SCIP *scip, SCP_Lagrange *lag, SCP_Lagrange_Sol *mult)
{
   SCIP_ROW **rows;
   int nrows;
   int i;

   /* compute subgradient: s_i(u) = 1 - sum_{j \in J_i} x_j(u) */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   for(i = 0; i < nrows; i++)
   {
      int rowidx = SCIProwGetIndex(rows[i]);
      if(lag->covered[rowidx] == FALSE)
      {
         SCIP_COL **covered_by = SCIProwGetCols(rows[i]);
         int numcols = SCIProwGetNNonz(rows[i]);
         int j;

         mult->subgradient[rowidx] = 1.0;

         for(j = 0; j < numcols; j++)
         {
            int colidx = SCIPcolGetIndex(covered_by[j]);

            assert(lag->fixedvars[colidx] == FALSE);

            mult->subgradient[rowidx] = mult->subgradient[rowidx] - mult->x_lagrange[colidx];
         }
      }
      else
         mult->subgradient[rowidx] = 0.0;
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE subgradientOptimization(SCIP *scip, SCP_Lagrange *lag, SCP_Lagrange_Sol *best_mult_lb, SCP_Lagrange_Sol *best_mult_ub)
{
   int max_iter;
   unsigned int seed;
#ifdef SCP_LAMBDA_ADJUSTMENTS
   SCIP_Real last_lb[SCP_LAMBDA_P];
   int last_data_pos = 0;
#endif
   SCP_Lagrange_Sol last_mult, next_mult, tmp_mult;
   SCIP_Real lambda = 0.1;
   SCIP_Real norm;
   SCIP_Real stop_crit_lb = 0.0;
   int iter, i;

   SCIP_CALL( allocateMemoryForSolution(scip, lag, &next_mult) );
   SCIP_CALL( allocateMemoryForSolution(scip, lag, &last_mult) );

   /* save data from best lower bound multiplier in last_mult */
   SCIP_CALL( copySolution(lag, &last_mult, best_mult_lb) );

   /* permutate best u by multiplying each entry with a uniformly random value in the range [0.9, 1.1] */
   seed = (unsigned int) SCIPround(scip, SCIPclockGetTimeOfDay());
   for(i = 0; i < lag->maxrows; i++)
   {
      if(lag->covered[i] == FALSE)
         last_mult.u[i] = SCIPgetRandomReal(0.9, 1.1, &seed) * last_mult.u[i];
   }

   /* subgradient optimization */
   max_iter = 10 * lag->maxrows;

   for(iter = 0; iter < max_iter; iter++)
   {
      SCIP_ROW **rows;
      int nrows;

      SCIP_CALL( computeSubgradient(scip, lag, &last_mult) );

      /* compute norm of the subgradient */
      norm = 0.0;
      for(i = 0; i < lag->maxrows; i++)
      {
         if(lag->covered[i] == FALSE)
            norm = norm + last_mult.subgradient[i] * last_mult.subgradient[i];
      }

      /* Held-Karp update */
      SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
      for(i = 0; i < nrows; i++)
      {
         int rowidx = SCIProwGetIndex(rows[i]);
         if(lag->covered[rowidx] == FALSE)
         {
            SCIP_Real hk = last_mult.u[rowidx] + lambda * (best_mult_ub->ub_greedy - last_mult.lb_lagrange) * last_mult.subgradient[rowidx] / norm;
            next_mult.u[rowidx] = (hk > 0) ? hk : 0.0;
         }
         else
            next_mult.u[rowidx] = 0.0;
      }

      SCIP_CALL( computeOptimalSolution(scip, lag, &next_mult) );
      SCIP_CALL( greedySetCover(scip, lag, &next_mult) );

      if(next_mult.ub_greedy < best_mult_ub->ub_greedy)
         SCIP_CALL( copySolution(lag, best_mult_ub, &next_mult) );

      if(next_mult.lb_lagrange > best_mult_lb->lb_lagrange)
         SCIP_CALL( copySolution(lag, best_mult_lb, &next_mult) );

      /*SCIPdebugMessage("so: iter %i, lb: %f, ub: %f, best lb: %f, best ub: %f\n", iter, next_mult.lb_lagrange, next_mult.ub_greedy, best_mult_lb->lb_lagrange, best_mult_ub->ub_greedy);*/

#ifdef SCP_LAMBDA_ADJUSTMENTS
      /* save last 'p' lower and upper bounds */
      last_lb[last_data_pos++] = next_mult.lb_lagrange;

      if(last_data_pos >= SCP_LAMBDA_P)
      {
         SCIP_Real max_lb = last_lb[0];
         SCIP_Real min_lb = last_lb[0];

         for(i = 1; i < SCP_LAMBDA_P; i++)
         {
            if(last_lb[i] > max_lb)
               max_lb = last_lb[i];
            if(last_lb[i] < min_lb)
               min_lb = last_lb[i];
         }

         /* NOTE: without these adjustments, the lower and upper bound is much better on the current test instance */
         /* if min_lb and max_lb differ by more than 1%, lambda is halved */
         if(max_lb - min_lb > 0.01)
            lambda = lambda / 2;

         /* if min_lb and max_lb differ by less than 0.1%, lambda is multiplied by 1.5 */
         if(max_lb - min_lb < 0.001)
            lambda = lambda * 1.5;

         last_data_pos = 0;
      }
#endif

      /* swap next_mult and last_mult */
      tmp_mult = last_mult;
      last_mult = next_mult;
      next_mult = tmp_mult;

      if(iter % SCP_STOP_CRIT_ITER == 0)
      {
         if((iter > 0) && (best_mult_lb->lb_lagrange - stop_crit_lb <= SCP_STOP_CRIT_DIFF) && (stop_crit_lb / best_mult_lb->lb_lagrange >= SCP_STOP_CRIT_PER))
            break;

         stop_crit_lb = best_mult_lb->lb_lagrange;
      }
   }

   SCIPdebugMessage("num iterations: %i\n", iter);
   SCIP_CALL( freeMemoryForSolution(scip, &next_mult) );
   SCIP_CALL( freeMemoryForSolution(scip, &last_mult) );

   return SCIP_OKAY;
}

static SCIP_RETCODE computeInitialLagrangeMultiplier(SCIP *scip, SCP_Lagrange *lag, SCP_Lagrange_Sol *mult)
{
   SCIP_ROW **rows;
   int nrows;
   int i;

   /* find first lagrangian multiplier */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   for(i = 0; i < nrows; i++)
   {
      if(lag->covered[SCIProwGetIndex(rows[i])] == FALSE)
      {
         SCIP_COL **covered_by = SCIProwGetCols(rows[i]);
         int numcols = SCIProwGetNNonz(rows[i]);
         int j;
         SCIP_Bool found = FALSE;
         SCIP_Real mincost = 0.0;

         assert(numcols > 0);

         for(j = 0; j < numcols; j++)
         {
            int idx = SCIPcolGetIndex(covered_by[j]);
            int numuncovered = countUncoveredRows(lag, covered_by[j]);
            SCIP_Real cost = lag->costs[idx] / numuncovered;

            if((found == FALSE) || (cost < mincost))
               mincost = cost;

            found = TRUE;
         }

         mult->u[SCIProwGetIndex(rows[i])] = mincost;
      }
      else
         mult->u[SCIProwGetIndex(rows[i])] = 0.0;
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE reportSolution(SCIP *scip, SCP_Lagrange *lag, SCP_Lagrange_Sol *mult, SCIP_HEUR *heur)
{
   SCIP_VAR **solvars;
   SCIP_Real *solvals;
   int nsolvars, i;
   SCIP_SOL *newsol;
   SCIP_Bool success = FALSE;
   SCIP_Bool foundsol = FALSE;
   SCIP_CONS **cons;
   int ncons;

   SCIP_CALL( SCIPgetVarsData(scip, &solvars, &nsolvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nsolvars) );

   for(i = 0; i < nsolvars; i++)
   {
      int varidx = SCIPvarGetIndex(solvars[i]);
      if((mult->x_greedy[varidx] == TRUE))
         solvals[i] = 1.0;
      else if((lag->fixedvars[varidx] == TRUE) && (SCIPisZero(scip, SCIPvarGetObj(solvars[i])) == FALSE))
         solvals[i] = 1.0;
      else
         solvals[i] = 0.0;
   }

   SCIP_CALL( SCIPsetSolVals(scip, newsol, nsolvars, solvars, solvals) );

   /* test all constraints and check if the activity is correct, adjust free variable if necessary */
   cons = SCIPgetConss(scip);
   ncons = SCIPgetNConss(scip);
   for(i = 0; (i < ncons) && (foundsol == FALSE); i++)
   {
      SCIP_Real lhs = 0.0, activity = 0.0;
      SCIP_Real *vals = NULL;
      SCIP_Bool valuesallones = FALSE;

      if(!strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])), "linear"))
      {
         lhs = SCIPgetLhsLinear(scip, cons[i]);
         /*rhs = SCIPgetRhsLinear(scip, cons[i]);*/
         activity = SCIPgetActivityLinear(scip, cons[i], newsol);
         vals = SCIPgetValsLinear(scip, cons[i]);
      }
      else if(!strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])), "logicor"))
      {
         valuesallones = TRUE;
      }
      else if(!strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])), "masterbranch"))
      {
         /* do nothing */
      }
      else
      {
         SCIPerrorMessage("constraint is '%s', can't handle this\n", SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])));
         SCIPABORT();
      }

      if(lhs > activity)
      {
         int j, nvars;
         SCIP_Bool changed = FALSE;
         SCIP_VAR **vars;

         SCIP_CALL( SCIPgetConsNVars(scip, cons[i], &nvars, &success) );
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
         SCIP_CALL( SCIPgetConsVars(scip, cons[i], vars, nvars, &success) );

         for(j = 0; (changed == FALSE) && (j < nvars); j++)
         {
            int nmastervars = GCGmasterVarGetNOrigvars(vars[j]);
            SCIP_Bool costszero = SCIPisZero(scip, SCIPvarGetObj(vars[j]));
            SCIP_Bool changevar = FALSE;

            if(costszero == TRUE)
            {
               changevar = (nmastervars == 0);
               if(changevar == FALSE)
               {
                  SCIP_Real *mastervals = GCGmasterVarGetOrigvals(vars[j]);
                  int k;

                  changevar = TRUE;
                  for(k = 0; (k < nmastervars) && (changevar == TRUE); k++)
                  {
                     if(SCIPisZero(scip, mastervals[k]) == FALSE)
                        changevar = FALSE;
                  }

                  if(changevar == TRUE)
                     SCIPdebugMessage("master variable has only original variables with zero costs\n");
               }
               else
                  SCIPdebugMessage("master variable has no original variables\n");

               if(changevar == TRUE)
               {
                  if(valuesallones)
                     SCIP_CALL( SCIPincSolVal(scip, newsol, vars[j], lhs - activity) );
                  else if(vals[j] != 0.0)
                     SCIP_CALL( SCIPincSolVal(scip, newsol, vars[j], (lhs - activity) / vals[j]) );

                  SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, &foundsol) );
                  changed = TRUE;
               }
            }
         }

         SCIPfreeBufferArray(scip, &vars);
      }

      /* TODO: take care of the case rhs < activity. Question: can this occur? */
   }

   SCIP_CALL( SCIPtrySolFree(scip, &newsol, TRUE, TRUE, TRUE, TRUE, &success) );
   SCIPfreeBufferArray(scip, &solvals);

   if(success == TRUE)
      SCIPdebugMessage("new solution found by set covering heuristic\n");

   return SCIP_OKAY;
}

static SCIP_RETCODE computeDelta(SCIP *scip, SCP_Lagrange *lag, SCP_Lagrange_Sol *mult_ub)
{
   SCIP_ROW **rows;
   SCIP_COL **columns;
   SCIP_Real *numcovered;
   int nrows;
   int ncolumns;
   int i, j;

   for(i = 0; i < lag->maxcolumns; i++)
      mult_ub->delta_index[i] = i;

   SCIP_CALL( SCIPallocBufferArray(scip, &numcovered, lag->maxrows) );

   /* 1. compute for each row by how many columns it is covered */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   for(i = 0; i < nrows; i++)
   {
      SCIP_COL **covered_by = SCIProwGetCols(rows[i]);
      int numcols = SCIProwGetNNonz(rows[i]);

      numcovered[i] = 0;
      for(j = 0; j < numcols; j++)
      {
         int idx = SCIPcolGetIndex(covered_by[j]);
         if(mult_ub->x_greedy[idx] == TRUE)
            numcovered[i] = numcovered[i] + 1;
      }
   }

   /* 2. for each column that is in the solution, compute delta */
   SCIP_CALL( SCIPgetLPColsData(scip, &columns, &ncolumns) );
   for(i = 0; i < ncolumns; i++)
   {
      int colidx = SCIPcolGetIndex(columns[i]);

      mult_ub->delta[colidx] = 0.0;

      if(mult_ub->x_greedy[colidx] == TRUE)
      {
         rows = SCIPcolGetRows(columns[i]);
         nrows = SCIPcolGetNNonz(columns[i]);

         if(mult_ub->lagrangian_costs[colidx] > 0)
            mult_ub->delta[colidx] = mult_ub->lagrangian_costs[colidx];

         for(j = 0; j < nrows; j++)
         {
            int rowidx = SCIProwGetIndex(rows[j]);
            mult_ub->delta[colidx] = mult_ub->delta[colidx] + mult_ub->u[rowidx] * (numcovered[rowidx] - 1) / numcovered[rowidx];
         }
      }
   }

   SCIPfreeBufferArray(scip, &numcovered);

   return SCIP_OKAY;
}

static SCIP_RETCODE threePhase(SCIP *scip, SCP_Lagrange *lag, SCP_Lagrange_Sol *mult_lb, SCP_Lagrange_Sol *mult_ub, SCIP_HEUR *heur)
{
   SCIP_Bool improved = TRUE;

   SCIP_CALL( computeOptimalSolution(scip, lag, mult_lb) );
   SCIP_CALL( greedySetCover(scip, lag, mult_lb) );

   while(improved == TRUE)
   {
      improved = FALSE;

      SCIP_CALL( subgradientOptimization(scip, lag, mult_lb, mult_ub) );
      SCIPdebugMessage("lb: %f, ub: %f\n", mult_lb->lb_lagrange, mult_ub->ub_greedy);

   }

   SCIP_CALL( reportSolution(scip, lag, mult_ub, heur) );

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_HEURCOPY(heurCopySetcover)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of setcover primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurCopySetcover NULL
#endif

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_HEURFREE(heurFreeSetcover)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of setcover primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurFreeSetcover NULL
#endif


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
SCIP_DECL_HEURINIT(heurInitSetcover)
{
   return SCIP_OKAY;
}
#else
#define heurInitSetcover NULL
#endif

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitSetcover)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of setcover primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitSetcover NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolSetcover)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of setcover primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolSetcover NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolSetcover)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of setcover primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolSetcover NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSetcover)
{
   SCIP *origprob;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   origprob = GCGmasterGetOrigprob(scip);
   assert(origprob != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMessage("set cover heuristic called\n");

   if(GCGisMasterSetCovering(origprob) == FALSE)
      return SCIP_OKAY;

   if(SCIPgetNVars(scip) == 0)
      return SCIP_OKAY;
   else
   {
      SCP_Core core;
      SCIP_CALL( initTentativeCore(scip, &core) );
      SCIP_CALL( computeCoreColumns(scip, &core) );
      SCIP_CALL( freeCore(scip, &core) );
   }

   if(SCIPgetNVars(scip) > 0)
   {
      SCP_Lagrange lag;
      SCP_Lagrange_Sol best_mult_lb, best_mult_ub;
      SCP_Lagrange_Sol best_total_ub;
      SCIP_VAR **vars = SCIPgetVars(scip);
      SCIP_ROW **rows;
      SCIP_Bool best_total_valid = FALSE;
      SCIP_Bool stopcft = FALSE;
      int maxval;
      int num;
      int i;
      num = 0;

      /* find maximum index of any variable */
      num = SCIPgetNVars(scip);
      maxval = 0;
      for(i = 0; i < num; i++)
      {
         SCIP_VAR *var = vars[i];
         if(SCIPvarGetIndex(var) >= maxval)
            maxval = SCIPvarGetIndex(var) + 1;
      }

      lag.maxcolumns = maxval;

      /* allocate memory to store the interval [0, maxval] */
      SCIP_CALL( SCIPallocBufferArray(scip, &lag.fixedvars, maxval) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lag.costs, maxval) );
      for(i = 0; i < maxval; i++)
      {
         int idx = SCIPvarGetIndex(vars[i]);
         lag.costs[idx] = SCIPvarGetObj(vars[i]);

         /* mark all variables with costs <= 0 as fixed */
         lag.fixedvars[i] = (lag.costs[idx] > 0) ? FALSE : TRUE;
      }

      /* mark all variables with costs <= 0 as fixed */
      for(i = 0; i < num; i++)
      {
         if(SCIPvarGetObj(vars[i]) <= 0)
         {
            lag.fixedvars[SCIPvarGetIndex(vars[i])] = TRUE;
            lag.numfixed = lag.numfixed + 1;
            /* lag.costs_fixed = lag.costs_fixed + lag.costs[SCIPvarGetIndex(vars[i])]; */
            /*SCIPdebugMessage("fixing variable %i because it has zero costs\n", SCIPvarGetIndex(vars[i]));*/
         }
      }

      /* find maximum index of any row */
      SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &num) );
      maxval = 0;
      for(i = 0; i < num; i++)
      {
         if(SCIProwGetIndex(rows[i]) >= maxval)
            maxval = SCIProwGetIndex(rows[i]) + 1;
      }
      lag.maxrows = maxval;

      SCIP_CALL( SCIPallocBufferArray(scip, &lag.covered, maxval) );

      /* mark all rows as uncovered */
      for(i = 0; i < lag.maxrows; i++)
         lag.covered[i] = FALSE;

      markRowsCoveredByFixedVariables(scip, &lag);

      SCIP_CALL( allocateMemoryForSolution(scip, &lag, &best_mult_lb) );
      SCIP_CALL( allocateMemoryForSolution(scip, &lag, &best_mult_ub) );
      SCIP_CALL( allocateMemoryForSolution(scip, &lag, &best_total_ub) );
      SCIP_CALL( computeInitialLagrangeMultiplier(scip, &lag, &best_mult_lb) );
      SCIP_CALL( greedySetCover(scip, &lag, &best_mult_lb) );
      SCIP_CALL( copySolution(&lag, &best_mult_ub, &best_mult_lb) );

      while(stopcft == FALSE)
      {
         SCIP_Bool foundbest = FALSE;
         SCIP_Real pi = 0.0;
         int numrowscovered, maxnumrowscovered;

         /* 1. define the sub-instance I_F derived from fixing x_j = 1 for all J in F */
         for(i = 0; i < lag.maxrows; i++)
            lag.covered[i] = FALSE;

         markRowsCoveredByFixedVariables(scip, &lag);

         /* 2. apply procedure 3-PHASE(x, u) to I_F */
         SCIP_CALL( threePhase(scip, &lag, &best_mult_lb, &best_mult_ub, heur) );

         /* 3. obtain the SCP solution x* associated with the partial solution x */
         for(i = 0; i < lag.maxcolumns; i++)
         {
            if(lag.fixedvars[i] == TRUE)
            {
               best_mult_ub.x_greedy[i] = TRUE;
               best_mult_ub.ub_greedy = best_mult_ub.ub_greedy + lag.costs[i];
            }
         }

         /* 5. update the best SCP solution */
         if((best_total_valid == FALSE) || (best_mult_ub.ub_greedy < best_total_ub.ub_greedy))
         {
            SCIP_CALL( copySolution(&lag, &best_total_ub, &best_mult_lb) );
            best_total_valid = TRUE;
            foundbest = TRUE;
         }

         /* 6. define the new value for the threshold pi */
         if((foundbest == TRUE) || (lag.numfixed == 0))
            pi = SCP_PI_MIN;
         else
            pi = pi * SCP_PI_ALPHA;

         /* 7. update the set of columns fixed */
         SCIP_CALL( computeDelta(scip, &lag, &best_total_ub) );

         SCIPsortRealInt(best_total_ub.delta, best_total_ub.delta_index, lag.maxcolumns);

         /* unfix all variables */
         lag.numfixed = 0;
         for(i = 0; i < lag.maxcolumns; i++)
            lag.fixedvars[i] = FALSE;

         for(i = 0; i < lag.maxrows; i++)
            lag.covered[i] = FALSE;

         numrowscovered = 0;
         maxnumrowscovered = lag.maxrows * pi; /* TODO: this needs to be nrows */

         for(i = 0; (i < lag.maxcolumns) && (numrowscovered < maxnumrowscovered); i++)
         {
            int colidx = best_total_ub.delta_index[i];
            if(best_total_ub.x_greedy[colidx] == TRUE)
            {
               /* fix this variable */

               /* find all rows that are covered by this column and check if they are already covered */
            }
         }

         /* TODO: define core */
         stopcft = TRUE;
      }

      SCIP_CALL( freeMemoryForSolution(scip, &best_mult_lb) );
      SCIP_CALL( freeMemoryForSolution(scip, &best_mult_ub) );
      SCIP_CALL( freeMemoryForSolution(scip, &best_total_ub) );
      SCIPfreeBufferArray(scip, &lag.fixedvars);
      SCIPfreeBufferArray(scip, &lag.covered);
      SCIPfreeBufferArray(scip, &lag.costs);
   }
   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the setcover primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSetcover(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata = NULL;
   SCIP_HEUR* heur;

   heur = NULL;

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopySetcover, heurFreeSetcover, heurInitSetcover, heurExitSetcover, heurInitsolSetcover, heurExitsolSetcover, heurExecSetcover,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecSetcover, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopySetcover) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeSetcover) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitSetcover) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitSetcover) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolSetcover) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolSetcover) );
#endif

   /* add setcover primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
