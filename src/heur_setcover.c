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

/* #define SCIP_DEBUG*/

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
/*#define SCP_LAMBDA_ADJUSTMENTS*/
#define SCP_LAMBDA_P          50
#define SCP_STOP_CRIT_ITER    300
#define SCP_STOP_CRIT_DIFF    1.0
#define SCP_STOP_CRIT_PER     0.99
#define SCP_PI_MIN            0.3
#define SCP_PI_ALPHA          1.1
#define SCP_BETA              1.0
#define SCP_GREEDY_SIMPLE

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
   SCIP_HASHTABLE *varsfixed;
   SCIP_HASHTABLE *rowscovered;
} SCP_Instance;

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
   int nactiveconstraints; /* total number of active constraints for which the variables can be retrieved */
   int maxconstraintvariables; /* greatest number of variables some constraint contains */
   SCIP_CONS **constraints;
} SCP_Core;

typedef struct
{
   SCIP_HASHTABLE *x_greedy_local; /* contains variables that are part of a greedy solution */
   SCIP_Real *u; /* lagrange multipliers for the rows */
   SCIP_Real *subgradient;
   SCIP_Real *lagrangian_costs_local;
   SCIP_Real *lagrangian_costs_global;
   SCIP_Real ub_greedy_local; /* bound computed by the greedy set cover algorithm for the restricted instance */
   SCIP_Real lb_lagrange_local; /* lower bound by lagrange relaxation for the restricted instance */
   SCIP_Real lb_lagrange_global; /* lower bound by lagrange relaxation for the unrestricted instance */
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


static SCIP_DECL_HASHGETKEY(hashGetKeyInt)
{
   return elem;
}


static SCIP_DECL_HASHKEYEQ(hashKeyEqInt)
{
   int var1 = (int) (size_t) key1;
   int var2 = (int) (size_t) key2;

   if(var1 == var2)
      return TRUE;
   else
      return FALSE;
}

static SCIP_DECL_HASHKEYVAL(hashKeyValInt)
{
   int var = (int) (size_t) key;
   return (unsigned int) var;
}

static SCIP_RETCODE allocateMemoryForSolution(SCIP *scip, SCP_Core *core, SCP_Lagrange_Sol *mult)
{
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->u, core->nconstraints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->subgradient, core->nconstraints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->lagrangian_costs_local, core->nvariables) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->lagrangian_costs_global, core->nvariables) );
   SCIP_CALL( SCIPhashtableCreate(&mult->x_greedy_local, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );
   return SCIP_OKAY;
}

static SCIP_RETCODE freeMemoryForSolution(SCIP *scip, SCP_Lagrange_Sol *mult)
{
   SCIPfreeBufferArray(scip, &mult->u);
   SCIPfreeBufferArray(scip, &mult->subgradient);
   SCIPfreeBufferArray(scip, &mult->lagrangian_costs_local);
   SCIPfreeBufferArray(scip, &mult->lagrangian_costs_global);
   SCIPhashtableFree(&mult->x_greedy_local);
   return SCIP_OKAY;
}

static SCIP_RETCODE copySolution(SCP_Core *core, SCP_Lagrange_Sol *dest, SCP_Lagrange_Sol *source)
{
   int i;

   SCIPhashtableRemoveAll(dest->x_greedy_local);

   for(i = 0; i < core->nvariables; i++)
   {
      dest->lagrangian_costs_local[i] = source->lagrangian_costs_local[i];
      dest->lagrangian_costs_global[i] = source->lagrangian_costs_global[i];

      if(SCIPhashtableExists(source->x_greedy_local, core->variables[i]) == TRUE)
         SCIPhashtableInsert(dest->x_greedy_local, core->variables[i]);
   }

   for(i = 0; i < core->nconstraints; i++)
   {
      dest->u[i] = source->u[i];
      dest->subgradient[i] = source->subgradient[i];
   }

   dest->lb_lagrange_global = source->lb_lagrange_global;
   dest->lb_lagrange_local = source->lb_lagrange_local;
   dest->ub_greedy_local = source->ub_greedy_local;

   return SCIP_OKAY;
}

static SCIP_RETCODE initInstance(SCIP *scip, SCP_Instance *inst)
{
   SCIP_CALL( SCIPhashtableCreate(&inst->varsfixed, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );
   SCIP_CALL( SCIPhashtableCreate(&inst->rowscovered, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyInt, hashKeyEqInt, hashKeyValInt, NULL) );
   return SCIP_OKAY;
}

static SCIP_RETCODE copyInstance(SCIP *scip, SCP_Core *core, SCP_Instance *dest, SCP_Instance *source)
{
   int i;

   for(i = 0; i < core->nvariables; i++)
   {
      if(SCIPhashtableExists(source->varsfixed, core->variables[i]) == TRUE)
         SCIPhashtableSafeInsert(dest->varsfixed, core->variables[i]);
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE freeInstance(SCIP *scip, SCP_Instance *inst)
{
   SCIPhashtableFree(&inst->varsfixed);
   SCIPhashtableFree(&inst->rowscovered);
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
      {
         SCIPdebugMessage("constraint %i (%s): cant get number of variables\n", i, SCIPconsGetName(core->constraints[i]));
         continue;
      }

      SCIP_CALL( SCIPgetConsVars(scip, core->constraints[i], vars, core->nvariables, &success) );
      if(success == FALSE)
      {
         SCIPdebugMessage("constraint %i (%s): cant get variables\n", i, SCIPconsGetName(core->constraints[i]));
         continue;
      }

      if(nvars > core->maxconstraintvariables)
         core->maxconstraintvariables = nvars;

      for(j = 0; j < nvars; j++)
      {
         int varidx = SCIPvarGetIndex(vars[j]);
         int varpos = (int) (size_t) SCIPhashmapGetImage(core->mapvariables, (void *) (size_t) varidx);

         if(j < SCP_CORE_TENT_SIZE)
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

/* adds all indices of rows to inst->rowscovered for all rows that are covered by the variables in inst->varsfixed */
static SCIP_RETCODE markRowsCoveredByFixedVariables(SCIP *scip, SCP_Core *core, SCP_Instance *inst)
{
   int i, j;

//   if(core->columnsavailable == TRUE)
//   {
//      /* iterate through all variables and check if they are fixed */
//      for(i = 0; i < core->nvariables; i++)
//      {
//         if(SCIPhashtableExists(inst->varsfixed, core->variables[i]) == FALSE)
//            continue;
//
//         if(SCIPhashtableExists(core->corevariables, core->variables[i]) == FALSE)
//         {
//            SCIPdebugMessage("fixed variable is not in the core, this should not happen!\n");
//            continue;
//         }
//
//         /* iterate through the column of this variable and mark all rows as fixed */
//         for(j = 0; j < core->nvarconstraints[i]; i++)
//         {
//            if(SCIPhashtableExists(inst->rowscovered, (void *) (size_t) core->columns[i][j]) == FALSE)
//               SCIPhashtableSafeInsert(inst->rowscovered, (void *) (size_t) core->columns[i][j]);
//         }
//      }
//   }
//   else
   {
      SCIP_VAR **vars;
      int nvars;

      SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );

      for(i = 0; i < core->nconstraints; i++)
      {
         SCIP_Bool success = FALSE;

         if(SCIPconsIsActive(core->constraints[i]) == FALSE)
            continue;

         SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[i], &nvars, &success) );
         if(success == FALSE)
            continue;

         SCIP_CALL( SCIPgetConsVars(scip, core->constraints[i], vars, core->maxconstraintvariables, &success) );
         if(success == FALSE)
            continue;

         for(j = 0; j < nvars; j++)
         {
            if(SCIPhashtableExists(inst->varsfixed, vars[j]) == TRUE)
            {
               if(SCIPhashtableExists(inst->rowscovered, (void *) (size_t) i) == FALSE)
                  SCIPhashtableInsert(inst->rowscovered, (void *) (size_t) i);
               break;
            }
         }
      }

      SCIPfreeBufferArray(scip, &vars);
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE checkSetCover(SCIP *scip, SCP_Core *core, SCP_Lagrange_Sol *mult, SCIP_Bool *issetcover)
{
   SCIP_Bool success;
   SCIP_VAR **vars;
   int nvars;
   int i, j;

   *issetcover = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );

   /* iterate through all constraints and check whether each of them constains a variable that is part of the cover */
   for(i = 0; i < core->nconstraints; i++)
   {
      SCIP_Bool rowcovered = FALSE;

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

         if(SCIPhashtableExists(core->corevariables, core->variables[varpos]) == TRUE)
         {
            rowcovered = TRUE;
            break;
         }
      }

      if(rowcovered == FALSE)
      {
         *issetcover = FALSE;
         break;
      }
   }

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

static SCIP_RETCODE greedySetCover(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *mult, int nfixvars, SCIP_Real *costsfixed)
{
   SCIP_HASHTABLE *rowscovered;
   int i, j;
   int nrowsuncovered = 0, nfixed = 0;

   mult->ub_greedy_local = 0.0;
   SCIPhashtableRemoveAll(mult->x_greedy_local);

#ifndef SCP_GREEDY_SIMPLE
   SCIPerrorMessage("sophisticated greedy algorithm is not implemented yet\n");
   SCIPABORT();
#endif

   SCIP_CALL( SCIPhashtableCreate(&rowscovered, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyInt, hashKeyEqInt, hashKeyValInt, NULL) );

   for(i = 0; i < core->nconstraints; i++)
   {
      SCIP_Bool success;
      int nvars;

      if(SCIPconsIsActive(core->constraints[i]) == FALSE)
         continue;

      /* this is actually necessary because there exist constraints where this fails, and we simply need to ignore them */
      SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[i], &nvars, &success) );
      if(success == FALSE)
         continue;

      if(SCIPhashtableExists(inst->rowscovered, (void *) (size_t) i) == TRUE)
         SCIPhashtableSafeInsert(rowscovered, (void *) (size_t) i);
      else
         nrowsuncovered++;
   }

#ifdef SCP_GREEDY_SIMPLE
   if(core->columnsavailable == FALSE)
   {
      SCIPerrorMessage("simple greedy algorithm requires core columns to be available\n");
      SCIPABORT();
   }

   while(nrowsuncovered > 0)
   {
      SCIP_Real minscore = SCIP_REAL_MAX;
      int mincolumn = -1;
      SCIP_Bool fixedvar = FALSE;

      /* compute scores for all core columns */
      for(i = 0; i < core->nvariables; i++)
      {
         int muh = 0;
         SCIP_Real gamma = SCIPvarGetObj(core->variables[i]);

         if(SCIPhashtableExists(core->corevariables, core->variables[i]) == FALSE)
            continue;
         else if(SCIPhashtableExists(inst->varsfixed, core->variables[i]) == TRUE)
            continue;

         for(j = 0; j < core->nvarconstraints[i]; j++)
         {
            if(SCIPhashtableExists(rowscovered, (void *) (size_t) core->columns[i][j]) == FALSE)
            {
               gamma -= mult->u[core->columns[i][j]];
               muh++;
            }
         }

         /* skip columns that do not cover anything */
         if(muh > 0)
         {
            SCIP_Real score = (gamma > 0) ? (gamma / ((SCIP_Real) muh)) : (gamma * muh);
            if(score < minscore)
            {
               minscore = score;
               mincolumn = i;
            }
         }
      }

      if(mincolumn == -1)
      {
         SCIPdebugMessage("there exist uncovered rows but no columns can cover them\n");
         SCIPABORT();
      }

      /* add variable 'variables[mincolumn]' to the set cover */
      mult->ub_greedy_local += SCIPvarGetObj(core->variables[mincolumn]);
      SCIPhashtableSafeInsert(mult->x_greedy_local, core->variables[mincolumn]);

      if(nfixed < nfixvars)
      {
         SCIPhashtableSafeInsert(inst->varsfixed, core->columns[mincolumn]);
         nfixed++;
         fixedvar = TRUE;
         if(costsfixed != NULL)
            *costsfixed += SCIPvarGetObj(core->variables[mincolumn]);
      }

      for(j = 0; j < core->nvarconstraints[mincolumn]; j++)
      {
         int col = core->columns[mincolumn][j];
         if(SCIPhashtableExists(rowscovered, (void *) (size_t) col) == FALSE)
         {
            SCIPhashtableSafeInsert(rowscovered, (void *) (size_t) col);
            if(fixedvar == TRUE)
               SCIPhashtableSafeInsert(inst->rowscovered, (void *) (size_t) col);
            nrowsuncovered--;
         }
      }
   }
#endif

   /* TODO: remove redundant columns */

   SCIPhashtableFree(&rowscovered);
   return SCIP_OKAY;
}

/* computes lagrangian costs for ALL columns */
static SCIP_RETCODE computeGlobalLagrangianCosts(SCIP *scip, SCP_Core *core, SCP_Lagrange_Sol *mult)
{
   SCIP_VAR **vars;
   int nvars;
   int i, j;

   /* set all lagrangian costs to objective values */
   for(i = 0; i < core->nvariables; i++)
      mult->lagrangian_costs_global[i] = SCIPvarGetObj(core->variables[i]);

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );

   for(i = 0; i < core->nconstraints; i++)
   {
      SCIP_Bool success = FALSE;

      if(SCIPconsIsActive(core->constraints[i]) == FALSE)
         continue;

      /* get all variables that are part of this constraint */
      SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[i], &nvars, &success) );
      if(success == FALSE)
         continue;

      SCIP_CALL( SCIPgetConsVars(scip, core->constraints[i], vars, core->maxconstraintvariables, &success) );
      if(success == FALSE)
         continue;

      /* for each variable: subtract u[i] from the variable's costs */
      for(j = 0; j < nvars; j++)
      {
         int varidx = SCIPvarGetIndex(vars[j]);
         int varpos = (int) (size_t) SCIPhashmapGetImage(core->mapvariables, (void *) (size_t) varidx);

         mult->lagrangian_costs_global[varpos] -= mult->u[i];
      }
   }

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/* computes lagrangian costs only for columns in the core considering only rows from a reduced instance */
static SCIP_RETCODE computeLocalLagrangianCosts(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *mult)
{
   SCIP_VAR **vars;
   int nvars;
   int i, j;

   /* set all lagrangian costs to objective values */
   for(i = 0; i < core->nvariables; i++)
      mult->lagrangian_costs_local[i] = SCIPvarGetObj(core->variables[i]);

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );

   for(i = 0; i < core->nconstraints; i++)
   {
      SCIP_Bool success = FALSE;

      /* skip rows that are not part of the reduced instance */
      if(SCIPhashtableExists(inst->rowscovered, (void *) (size_t) i) == TRUE)
         continue;

      if(SCIPconsIsActive(core->constraints[i]) == FALSE)
         continue;

      /* get all variables that are part of this constraint */
      SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[i], &nvars, &success) );
      if(success == FALSE)
         continue;

      SCIP_CALL( SCIPgetConsVars(scip, core->constraints[i], vars, core->maxconstraintvariables, &success) );
      if(success == FALSE)
         continue;

      /* for each variable: subtract u[i] from the variable's costs */
      for(j = 0; j < nvars; j++)
      {
         int varidx = SCIPvarGetIndex(vars[j]);
         int varpos = (int) (size_t) SCIPhashmapGetImage(core->mapvariables, (void *) (size_t) varidx);

         if(SCIPhashtableExists(core->corevariables, vars[j]) == TRUE)
            mult->lagrangian_costs_local[varpos] -= mult->u[i];
      }
   }

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/* computes an optimal solution to the lagrangian relaxation, see formulae (4), (5) in the paper */
static SCIP_RETCODE computeOptimalSolution(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *mult)
{
   SCIP_VAR **vars;
   SCIP_Bool success;
   int nvars;
   int i, j;

   mult->lb_lagrange_local = 0.0;
   mult->lb_lagrange_global = 0.0;

   SCIP_CALL( computeLocalLagrangianCosts(scip, core, inst, mult) );

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );

   for(i = 0; i < core->nvariables; i++)
   {
      if(SCIPhashtableExists(core->corevariables, core->variables[i]) == FALSE)
         continue;

      if(mult->lagrangian_costs_local[i] < 0.0)
      {
         if(SCIPhashtableExists(inst->varsfixed, core->variables[i]) == FALSE)
            mult->lb_lagrange_local = mult->lb_lagrange_local + mult->lagrangian_costs_local[i];
      }
   }

   for(i = 0; i < core->nconstraints; i++)
   {
      mult->subgradient[i] = 0.0;

      if(SCIPconsIsActive(core->constraints[i]) == FALSE)
         continue;

      if(SCIPhashtableExists(inst->rowscovered, (void *) (size_t) i) == TRUE)
         continue;

      SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[i], &nvars, &success) );
      if(success == FALSE)
         continue;

      SCIP_CALL( SCIPgetConsVars(scip, core->constraints[i], vars, core->maxconstraintvariables, &success) );
      if(success == FALSE)
         continue;

      mult->subgradient[i] = 1.0;

      for(j = 0; j < nvars; j++)
      {
         int varidx = SCIPvarGetIndex(vars[j]);
         int varpos = (int) (size_t) SCIPhashmapGetImage(core->mapvariables, (void *) (size_t) varidx);
         if(SCIPhashtableExists(core->corevariables, vars[j]) == FALSE)
            continue;

         if(mult->lagrangian_costs_local[varpos] < 0.0)
            mult->subgradient[i] -= 1.0;
      }

      mult->lb_lagrange_local = mult->lb_lagrange_local + mult->u[i];
   }

   SCIPfreeBufferArray(scip, &vars);
   return SCIP_OKAY;
}

static SCIP_RETCODE subgradientOptimization(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *best_mult_lb, SCP_Lagrange_Sol *best_mult_ub, SCIP_Real costsfixed)
{
   int max_iter;
   unsigned int seed;
   SCP_Lagrange_Sol last_mult, next_mult, tmp_mult;
   SCIP_Real norm, lambda = 0.1;
   int iter, i;
#ifdef SCP_LAMBDA_ADJUSTMENTS
   SCIP_Real last_lb[SCP_LAMBDA_P];
   int last_data_pos = 0;
#endif
   SCIP_Real stop_crit_lb = 0.0;

   SCIP_CALL( allocateMemoryForSolution(scip, core, &next_mult) );
   SCIP_CALL( allocateMemoryForSolution(scip, core, &last_mult) );

   /* save data from best lower bound multiplier in last_mult */
   SCIP_CALL( copySolution(core, &last_mult, best_mult_lb) );

   /* permutate best u by multiplying each entry with a uniformly random value in the range [0.9, 1.1] */
   seed = (unsigned int) SCIPround(scip, SCIPclockGetTimeOfDay());
   for(i = 0; i < core->nconstraints; i++)
   {
      if(SCIPhashtableExists(inst->rowscovered, (void *) (size_t) i) == FALSE)
         last_mult.u[i] = SCIPgetRandomReal(0.9, 1.1, &seed) * last_mult.u[i];
      else
         last_mult.u[i] = 0.0;
   }

   /* subgradient optimization */
   max_iter = 10 * core->nconstraints;

   for(iter = 0; iter < max_iter; iter++)
   {
      /* compute norm of the subgradient */
      norm = 0.0;
      for(i = 0; i < core->nconstraints; i++)
         norm = norm + last_mult.subgradient[i] * last_mult.subgradient[i]; /* we have subgradient[i] = 0.0 if row i is not to be considered */

      /* Held-Karp update */
      for(i = 0; i < core->nconstraints; i++)
      {
         SCIP_Real hk = last_mult.u[i] + lambda * (best_mult_ub->ub_greedy_local - last_mult.lb_lagrange_local) * last_mult.subgradient[i] / norm;
         next_mult.u[i] = 0.0;

         if(hk > 0.0)
            next_mult.u[i] = hk;
      }

      SCIP_CALL( computeOptimalSolution(scip, core, inst, &next_mult) );
      SCIP_CALL( greedySetCover(scip, core, inst, &next_mult, 0, NULL) );
      /*SCIPdebugMessage("local lb: %f, local ub: %f\n", next_mult.lb_lagrange_local, next_mult.ub_greedy_local);*/

      if(next_mult.ub_greedy_local + costsfixed < best_mult_ub->ub_greedy_local)
      {
         SCIP_CALL( copySolution(core, best_mult_ub, &next_mult) );
         best_mult_ub->ub_greedy_local += costsfixed;
      }

      if(next_mult.lb_lagrange_local > best_mult_lb->lb_lagrange_local)
         SCIP_CALL( copySolution(core, best_mult_lb, &next_mult) );

      /*SCIPdebugMessage("so: iter %i, lb: %f, ub: %f, best lb: %f, best ub: %f\n", iter, next_mult.lb_lagrange, next_mult.ub_greedy, best_mult_lb->lb_lagrange, best_mult_ub->ub_greedy);*/

#ifdef SCP_LAMBDA_ADJUSTMENTS
      /* save last 'p' lower and upper bounds */
      last_lb[last_data_pos++] = next_mult.lb_lagrange_local;

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
         if((iter > 0) && (best_mult_lb->lb_lagrange_local - stop_crit_lb <= SCP_STOP_CRIT_DIFF) && (stop_crit_lb / best_mult_lb->lb_lagrange_local >= SCP_STOP_CRIT_PER))
            break;

         stop_crit_lb = best_mult_lb->lb_lagrange_local;
      }
   }

   /*SCIPdebugMessage("num iterations: %i\n", iter);*/
   SCIP_CALL( freeMemoryForSolution(scip, &next_mult) );
   SCIP_CALL( freeMemoryForSolution(scip, &last_mult) );

   return SCIP_OKAY;
}

static SCIP_RETCODE computeInitialLagrangeMultiplier(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *mult)
{

   int i, j;

   if(core->columnsavailable == TRUE)
   {
      for(i = 0; i < core->nconstraints; i++)
         mult->u[i] = SCIP_REAL_MAX;

      for(i = 0; i < core->nvariables; i++)
      {
         int nuncovered = 0;

         if(SCIPhashtableExists(core->corevariables, core->variables[i]) == FALSE)
            continue;

         if(SCIPhashtableExists(inst->varsfixed, core->variables[i]) == TRUE)
            continue;

         /* count how many uncovered, active rows this column covers */
         for(j = 0; j < core->nvarconstraints[i]; j++)
         {
            int colpos = core->columns[i][j];

            /* skip inactive constraints */
            if(SCIPconsIsActive(core->constraints[colpos]) == FALSE)
               continue;

            /* skip covered rows */
            if(SCIPhashtableExists(inst->rowscovered, (void *) (size_t) colpos) == TRUE)
               continue;

            nuncovered++;
         }

         /* if this column covers any rows, update their cost if necessary */
         if(nuncovered > 0)
         {
            SCIP_Real costs = SCIPvarGetObj(core->variables[i]) / ((SCIP_Real) nuncovered);

            for(j = 0; j < core->nvarconstraints[i]; j++)
            {
               int colpos = core->columns[i][j];

               /* skip inactive constraints */
               if(SCIPconsIsActive(core->constraints[colpos]) == FALSE)
                  continue;

               /* skip covered rows */
               if(SCIPhashtableExists(inst->rowscovered, (void *) (size_t) colpos) == TRUE)
                  continue;

               if(costs < mult->u[colpos])
                  mult->u[colpos] = costs;
            }
         }
      }
   }
   else
   {
      SCIP_VAR **vars;
      int nvars;
      int *nuncoveredactive;

      SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nuncoveredactive, core->nvariables) );

      for(i = 0; i < core->nvariables; i++)
         nuncoveredactive[i] = 0;

      for(i = 0; i < core->nconstraints; i++)
      {
         SCIP_Bool success = FALSE;
         if(SCIPhashtableExists(inst->rowscovered, (void *) (size_t) i) == TRUE)
            continue;

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
            if(SCIPhashtableExists(core->corevariables, vars[j]) == FALSE)
               continue;

            if(SCIPhashtableExists(inst->varsfixed, vars[j]) == FALSE)
            {
               int varidx = SCIPvarGetIndex(vars[j]);
               int varpos = (int) (size_t) SCIPhashmapGetImage(core->mapvariables, (void *) (size_t) varidx);
               nuncoveredactive[varpos]++;
            }
         }
      }

      for(i = 0; i < core->nconstraints; i++)
      {
         SCIP_Bool found = FALSE;
         SCIP_Bool success = FALSE;

         if(SCIPhashtableExists(inst->rowscovered, (void *) (size_t) i) == TRUE)
         {
            mult->u[i] = 0.0;
            continue;
         }

         if(SCIPhashtableExists(inst->rowscovered, (void *) (size_t) i) == TRUE)
            continue;

         if(SCIPconsIsActive(core->constraints[i]) == FALSE)
            continue;

         SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[i], &nvars, &success) );
         if(success == FALSE)
            continue;

         SCIP_CALL( SCIPgetConsVars(scip, core->constraints[i], vars, core->maxconstraintvariables, &success) );
         if(success == FALSE)
            continue;

         for(j = 0; j < nvars; j++)
         {
            if(SCIPhashtableExists(core->corevariables, vars[j]) == FALSE)
               continue;
            else
            {
               int varidx = SCIPvarGetIndex(vars[j]);
               int varpos = (int) (size_t) SCIPhashmapGetImage(core->mapvariables, (void *) (size_t) varidx);
               SCIP_Real costs = SCIP_REAL_MAX;

               if(nuncoveredactive[varpos] > 0)
                  costs = SCIPvarGetObj(core->variables[varpos]) / ((SCIP_Real) nuncoveredactive[varpos]);

               if((found == FALSE) || (costs < mult->u[i]))
               {
                  found = TRUE;
                  mult->u[i] = costs;
               }
            }
         }
      }

      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &nuncoveredactive);
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE reportSolution(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *mult, SCIP_HEUR *heur)
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
      if(SCIPhashtableExists(mult->x_greedy_local, solvars[i]) == TRUE)
         solvals[i] = 1.0;
      else if((SCIPhashtableExists(inst->varsfixed, solvars[i]) == TRUE) && (SCIPisZero(scip, SCIPvarGetObj(solvars[i])) == FALSE))
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
         SCIPdebugMessage("constraint %i is a masterbranch\n", SCIPconsGetPos(cons[i]));
         continue;
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

         SCIPdebugMessage("constraint %i: left hand side is violated by %f\n", i, lhs - activity);
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
                  else
                     SCIPdebugMessage("could not adjust activity\n");

                  SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, &foundsol) );
                  changed = TRUE;
               }
            }
         }

         SCIPfreeBufferArray(scip, &vars);
      }

      if(!strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])), "linear"))
      {
         if(SCIPgetLhsLinear(scip, cons[i]) > SCIPgetActivityLinear(scip, cons[i], newsol))
            SCIPdebugMessage("activity is still smaller than lhs\n");
      }
      /* TODO: take care of the case rhs < activity. Question: can this occur? */
   }

   SCIP_CALL( SCIPtrySolFree(scip, &newsol, TRUE, TRUE, TRUE, TRUE, &success) );
   SCIPfreeBufferArray(scip, &solvals);

   if(success == TRUE)
      SCIPdebugMessage("new solution found by set covering heuristic\n");

   return SCIP_OKAY;
}

static SCIP_RETCODE threePhase(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *mult_lb, SCP_Lagrange_Sol *mult_ub, SCIP_HEUR *heur)
{
   SCP_Lagrange_Sol mult_lb_local;
   SCP_Instance subinst;
   int nfixvars, i;
   SCIP_Real costsfixed = 0.0;

   nfixvars = core->nconstraints / 200;
   if(nfixvars < 1)
      nfixvars = 1;

   SCIP_CALL( allocateMemoryForSolution(scip, core, &mult_lb_local) );

   /* we first create our own copy of the instance, as we need to mark variables as fixed until all variables are fixed */
   SCIP_CALL( initInstance(scip, &subinst) );
   SCIP_CALL( copyInstance(scip, core, &subinst, inst) );
   SCIP_CALL( markRowsCoveredByFixedVariables(scip, core, &subinst) );

   /* next, compute initial lagrange multipliers and find a first lower bound */
   SCIP_CALL( computeInitialLagrangeMultiplier(scip, core, &subinst, &mult_lb_local) );

   /* computeOptimalSolution also computes the subgradient */
   SCIP_CALL( computeOptimalSolution(scip, core, &subinst, &mult_lb_local) );
   SCIP_CALL( greedySetCover(scip, core, &subinst, &mult_lb_local, 0, NULL) );

   /* we now have a lower and upper bound in mult_lb_local for the instance 'inst'. We take these as our starting values */
   SCIP_CALL( copySolution(core, mult_ub, &mult_lb_local) );
   SCIP_CALL( copySolution(core, mult_lb, &mult_lb_local) );

   /* stop if all rows are covered by fixed variables */
   while(core->nactiveconstraints > SCIPhashtableGetNElements(subinst.rowscovered))
   {
      SCIP_CALL( subgradientOptimization(scip, core, &subinst, &mult_lb_local, mult_ub, costsfixed) );
      SCIP_CALL( greedySetCover(scip, core, &subinst, &mult_lb_local, nfixvars, &costsfixed) );

      SCIPdebugMessage("best upper bound: %f, best local lb: %f, fixed costs: %f\n", mult_ub->ub_greedy_local, mult_lb_local.lb_lagrange_local, costsfixed);
      /* mult_ub->greedy_local is an upper bound for the instance 'inst'. mult_lb_local.lb_lagrange_local is only a lower bound for the instance 'subinst' */
      if(costsfixed + mult_lb_local.lb_lagrange_local > mult_ub->ub_greedy_local)
         break;
   }

   /* in order to have a feasible set cover in 'mult_ub', we need to add all variables that were fixed */
   for(i = 0; i < core->nvariables; i++)
   {
      if(SCIPhashtableExists(inst->varsfixed, core->variables[i]) == FALSE)
         continue;

      if(SCIPhashtableExists(subinst.varsfixed, core->variables[i]) == TRUE)
         SCIPhashtableSafeInsert(mult_ub->x_greedy_local, core->variables[i]);
   }

   SCIP_CALL( freeMemoryForSolution(scip, &mult_lb_local) );
   SCIP_CALL( freeInstance(scip, &subinst) );

   return SCIP_OKAY;
}

static SCIP_RETCODE setCoveringHeuristic(SCIP *scip, SCIP_HEUR *heur)
{
   SCP_Core core;
   SCP_Instance inst;
   SCP_Lagrange_Sol best_mult_lb, best_mult_ub;
   SCP_Lagrange_Sol best_total_ub;
   SCIP_Bool stopcft = FALSE, foundsol = FALSE, success;
   SCIP_Bool computecorecolumns = TRUE;
   int i;

   /* basic setup, for each row find first five columns covering it */
   SCIP_CALL( initTentativeCore(scip, &core) );

   if(computecorecolumns == TRUE)
      SCIP_CALL( computeCoreColumns(scip, &core) );

   SCIP_CALL( initInstance(scip, &inst) );

   SCIP_CALL( allocateMemoryForSolution(scip, &core, &best_mult_lb) );
   SCIP_CALL( allocateMemoryForSolution(scip, &core, &best_mult_ub) );
   SCIP_CALL( allocateMemoryForSolution(scip, &core, &best_total_ub) );

   while(stopcft == FALSE)
   {
      /* 1. derive the reduced subinstance by marking rows covered by fixed variables */
      SCIPhashtableRemoveAll(inst.rowscovered);
      markRowsCoveredByFixedVariables(scip, &core, &inst);

      /* 2. apply procedure three-phase to find an optimal lagrange multiplier */
      SCIP_CALL( threePhase(scip, &core, &inst, &best_mult_lb, &best_mult_ub, heur) );

      /* 3. extend solution to a global solution by adding fixed variables */
      for(i = 0; i < core.nvariables; i++)
      {
         if(SCIPhashtableExists(inst.varsfixed, core.variables[i]) == TRUE)
         {
            SCIPhashtableSafeInsert(best_mult_ub.x_greedy_local, core.variables[i]);
            best_mult_ub.ub_greedy_local += SCIPvarGetObj(core.variables[i]);
         }
      }

      if((foundsol == FALSE) || (best_mult_ub.ub_greedy_local < best_total_ub.ub_greedy_local))
      {
         SCIP_CALL( copySolution(&core, &best_total_ub, &best_mult_ub) );
         foundsol = TRUE;
      }

      stopcft = TRUE;
   }

   SCIP_CALL( checkSetCover(scip, &core, &best_total_ub, &success) );
   if(success == FALSE)
      SCIPdebugMessage("final solution is not a valid set cover\n");

   SCIP_CALL( reportSolution(scip, &core, &inst, &best_total_ub, heur) );

   SCIP_CALL( freeMemoryForSolution(scip, &best_mult_lb) );
   SCIP_CALL( freeMemoryForSolution(scip, &best_mult_ub) );
   SCIP_CALL( freeMemoryForSolution(scip, &best_total_ub) );

   /*
    * SCIP_CALL( reportSolution(scip, core, inst, mult_ub, heur) );
    */

   SCIP_CALL( freeInstance(scip, &inst) );
   SCIP_CALL( freeCore(scip, &core) );

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

   SCIP_CALL( setCoveringHeuristic(scip, heur) );

   *result = SCIP_FOUNDSOL;
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
