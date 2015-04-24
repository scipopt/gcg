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

#define SCP_CORE_TENT_SIZE    10          /**< number of columns covering each row that are added to the tentative core at the beginning           */
#define SCP_LAMBDA_ADJUSTMENTS            /**< adjust step size during the subgradient phase                                                       */
#define SCP_LAMBDA_P          50          /**< number of iterations after which lambda is adjusted                                                 */
#define SCP_LAMBDA            0.1         /**< initial step size for the subgradient phase                                                         */
#define SCP_STOP_CRIT_ITER    300         /**< number of iterations of the subgradient phase after which the stopping criterion is tested again    */
#define SCP_STOP_CRIT_DIFF    1.0         /**< stop if absolute difference between best lower and upper bound is less than SCP_STOP_CRIT_DIFF, and */
#define SCP_STOP_CRIT_PER     0.99        /**<     the relative gap between best lower and upper bound is less than (1 - SCP_STOP_CRIT_PER         */
#define SCP_PI_MIN            0.3         /**< percentage of rows to be removed after fixing columns                                               */
#define SCP_PI_ALPHA          1.1         /**< increase of pi when no improvement was made, i.e. more columns will be fixed                        */
#define SCP_BETA              1.005       /**< allowed a gap between lowerbound and upper bound during the subgradient phase                       */
#define SCP_MAX_ITER          300         /**< maximum number of iterations of three-phase                                                         */
#define SCP_MAX_ITER_NO_IMP   10          /**< stop of no improvements during the last X iterations of three-phase                                 */
#define SCP_HEUR_MAX_ITER     250         /**< number of multipliers that are used for computing greedy solutions during each iteration            */

/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */

typedef struct
{
   int size;                     /**< actual number of elements stored in the queue    */
   int reserved;                 /**< number of elements for which memory is allocated */
   SCIP_Real *keys;              /**< array of keys                                    */
   int *data;                    /**< array of values                                  */
   int **positions;              /**< array of pointers to save positions of elements  */
} PQueue;

typedef struct
{
   SCIP_HASHTABLE *varsfixed;    /**< set that contains fixed variables                                  */
   SCIP_HASHTABLE *rowscovered;  /**< set that contains indices of rows covered by the fixed variables   */
   SCIP_Real costsfixed;         /**< total costs of variables that are fixed                            */
} SCP_Instance;

typedef struct
{
   SCIP_HASHTABLE *corevariables;   /**< set of core variables */
   int *listcorevariables;          /**< array of indices of core variables */
   int ncorevariables;              /**< number of core variables */
   SCIP_HASHMAP *mapvariables;      /**< maps variables-indices to [0, nvariables) in array 'variables' */
   SCIP_VAR **variables;            /**< all variables of the problem */
   int *nvarconstraints;            /**< array that contains for each variable the number of constraints it is part of */
   int nvariables;                  /**< total number of variables  */
   SCIP_Bool columnsavailable;      /**< if set then 'columns' contains the columns for all core variables */
   int **columns;                   /**< columns of core variables, NULL if not a core variable */
   SCIP_Bool rowsavailable;         /**< if set then 'rows' contains all rows reduced to core variables */
   int **rows;                      /**< rows that only contain core variables */
   int *nrowvars;                   /**< number of core variables a row contains */
   int nconstraints;                /**< total number of constraints (including inactive ones) */
   int nactiveconstraints;          /**< total number of active constraints for which the variables can be retrieved */
   int maxconstraintvariables;      /**< greatest number of variables some constraint contains */
   SCIP_CONS **constraints;         /**< all constraints of the problem */
   SCIP_Real *delta;                /**< delta values of variables */
   int *delta_pos;
   int *solgreedy;
   int nsolgreedy;
} SCP_Core;

typedef struct
{
   SCIP_HASHTABLE *x_greedy_local;     /**< contains variables that are part of a greedy solution, this is not necessarily a global solution */
   SCIP_Real *u;                       /**< lagrange multipliers for the rows                                                                */
   SCIP_Real *subgradient;
   SCIP_Real *lagrangian_costs_local;  /**< lagrangian costs (for a certain instance) when only uncovered rows are considered                */
   SCIP_Real *lagrangian_costs_global; /**< lagrangian costs for the whole instance when all rows and columns are considered                 */
   SCIP_Real ub_greedy_local;          /**< bound computed by the greedy set cover algorithm for the restricted instance                     */
   SCIP_Real lb_lagrange_local;        /**< lower bound by lagrange relaxation for the restricted instance                                   */
   SCIP_Real lb_lagrange_global;       /**< lower bound by lagrange relaxation for the unrestricted instance                                 */
} SCP_Lagrange_Sol;

struct SCIP_HeurData
{
   SCP_Core core;                         /**< core (subcollection of columns) of the problem covering all rows */
   SCP_Instance inst;                     /**< reduced instance where some variables may be fixed and some rows be covered */
   SCP_Instance subinst;                  /**< reduced instance of 'inst', used during the three-phase */

   SCP_Lagrange_Sol mult_best_lb_total;   /**< lagrange multiplier that gives the best lower bound for the complete problem */
   SCP_Lagrange_Sol mult_best_lb_inst;    /**< best multiplier for instance 'inst' */
   SCP_Lagrange_Sol mult_best_lb_subinst; /**< best multiplier for instance 'subinst' */

   SCIP_Real best_ub;                     /**< best upper bound that could be obtained so far */
   SCIP_HASHTABLE *best_ub_sol;           /**< actual solution that gives the best upper bound */
   SCIP_Real best_ub_inst;                /**< best upper bound for the reduced instance 'inst' (including fixed costs) */
   SCIP_HASHTABLE *best_ub_inst_sol;      /**< actual solution for the instance 'inst' (including fixed variables) */
   SCIP_Real best_ub_subinst;             /**< best upper bound for the reduced instance 'subinst' (including fixed costs) */
   SCIP_HASHTABLE *best_ub_subinst_sol;   /**< actual solution for the instance 'subinst' (including fixed variables) */

   /* local data that is used by (most) local procedures */
   SCIP_VAR **vars;                       /**< used to iterate through the variables of a constraint */

   /* memory that is used locally by 'threephase' */
   SCP_Lagrange_Sol tp_mult_lb_subinst;   /**< lagrange multiplier */

   /* memory that is used locally by 'greedysetcover' */
   PQueue greedy_queue;                   /**< priority queue used by the greedy algorithm */
   int *greedy_colpos;                    /**< stores position of a variable within the priority queue */
   int *greedy_colmu;                     /**< value mu for each variable */
   SCIP_Real *greedy_colgamma;            /**< value gamma for each variable */
   SCIP_Real *greedy_colscore;            /**< score of each variable */
   SCP_Instance greedy_inst;              /**< instance that is used to hold covered rows */
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** initializes a priority queue */
static SCIP_RETCODE pqueue_init(SCIP *scip, PQueue *queue)
{
   queue->size = 0;
   queue->reserved = 32;

   SCIP_CALL( SCIPallocBufferArray(scip, &queue->keys, queue->reserved) );
   SCIP_CALL( SCIPallocBufferArray(scip, &queue->data, queue->reserved) );
   SCIP_CALL( SCIPallocBufferArray(scip, &queue->positions, queue->reserved) );

   return SCIP_OKAY;
}

/** removes all elements from a priority queue, but does not release any memory */
static SCIP_RETCODE pqueue_clear(SCIP *scip, PQueue *queue)
{
   queue->size = 0;
   return SCIP_OKAY;
}

/** releases all memory that is used by a priority queue */
static SCIP_RETCODE pqueue_destroy(SCIP *scip, PQueue *queue)
{
   if(queue->keys != NULL)
      SCIPfreeBufferArray(scip, &queue->keys);

   if(queue->data != NULL)
      SCIPfreeBufferArray(scip, &queue->data);

   if(queue->positions != NULL)
      SCIPfreeBufferArray(scip, &queue->positions);

   return SCIP_OKAY;
}

/** Inserts an element with key 'key' and value 'elem' into the queue.
 *  If 'position' is not NULL, the referenced memory location will always contain the internal position of the element
 * */
static SCIP_RETCODE pqueue_insert(SCIP *scip, PQueue *queue, SCIP_Real key, int elem, int *position)
{
   int pos;
   int parent = 0;

   if(queue->size == queue->reserved)
   {
      queue->reserved += 128;
      SCIP_CALL( SCIPreallocBufferArray(scip, &queue->keys, queue->reserved) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &queue->data, queue->reserved) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &queue->positions, queue->reserved) );
   }

   pos = queue->size++;
   queue->keys[pos] = key;
   queue->data[pos] = elem;

   if(pos > 0)
      parent = (pos - 1) / 2;

   while((pos > 0) && (key < queue->keys[parent]))
   {
      /* swap element with parent */
      queue->keys[pos] = queue->keys[parent];
      queue->data[pos] = queue->data[parent];
      queue->positions[pos] = queue->positions[parent];
      /*queue->data[parent] = elem;*/

      if(queue->positions[pos] != NULL)
         *queue->positions[pos] = pos;

      pos = (pos - 1) / 2;
      parent = (pos - 1) / 2;
   }

   queue->keys[pos] = key;
   queue->data[pos] = elem;
   queue->positions[pos] = position;

   if(position != NULL)
      *position = pos;

   return SCIP_OKAY;
}

/** Decreases the key to 'key' of the element that is currently at position 'pos' */
static SCIP_RETCODE pqueue_decrease_key(SCIP *scip, PQueue *queue, int pos, SCIP_Real key)
{
   int parent;
   int *positionptr;
   int elem;

   if(pos >= queue->size)
      return SCIP_OKAY;

   parent = (pos - 1) / 2;
   positionptr = queue->positions[pos];
   elem = queue->data[pos];

   while((pos > 0) && (key < queue->keys[parent]))
   {
      /* swap element with parent */
      queue->keys[pos] = queue->keys[parent];
      queue->data[pos] = queue->data[parent];
      queue->positions[pos] = queue->positions[parent];
      queue->data[parent] = elem;

      if(queue->positions[pos] != NULL)
         *queue->positions[pos] = pos;

      pos = (pos - 1) / 2;
      parent = (pos - 1) / 2;
   }

   queue->keys[pos] = key;
   queue->data[pos] = elem;
   queue->positions[pos] = positionptr;

   if(positionptr != NULL)
      *positionptr = pos;

   return SCIP_OKAY;
}

/** Increases the key to 'key' of the element that is currently at position 'pos' */
static SCIP_RETCODE pqueue_increase_key(SCIP *scip, PQueue *queue, int pos, SCIP_Real key)
{
   int elem;
   int *positionptr;
   int left;
   int right;

   left = 2 * pos + 1;
   right = 2 * pos + 2;

   if(pos >= queue->size)
      return SCIP_OKAY;

   elem = queue->data[pos];
   positionptr = queue->positions[pos];

   while(left < queue->size)
   {
      if(right < queue->size)
      {
         /* both children exist, so swap with smallest child */
         if(key <= queue->keys[left])
         {
            /* left child is not smaller than element */
            if(key <= queue->keys[right])
               break;
            else
            {
               /* swap with right child */
               queue->keys[pos] = queue->keys[right];
               queue->data[pos] = queue->data[right];
               queue->positions[pos] = queue->positions[right];

               if(queue->positions[pos] != NULL)
                  *(queue->positions[pos]) = pos;

               pos = right;
            }
         }
         else
         {
            /* left child is smaller than the element */
            if(queue->keys[left] <= queue->keys[right])
            {
               /* left child is smaller than right child and smaller than the element */
               queue->keys[pos] = queue->keys[left];
               queue->data[pos] = queue->data[left];
               queue->positions[pos] = queue->positions[left];

               if(queue->positions[pos] != NULL)
                  *(queue->positions[pos]) = pos;

               pos = left;
            }
            else
            {
               /* right child is smallest element */
               queue->keys[pos] = queue->keys[right];
               queue->data[pos] = queue->data[right];
               queue->positions[pos] = queue->positions[right];

               if(queue->positions[pos] != NULL)
                  *(queue->positions[pos]) = pos;

               pos = right;
            }
         }
      }
      else
      {
         /* only left child exists */
         if(key > queue->keys[left])
         {
            /* left child is smaller than element, so swap them */
            queue->keys[pos] = queue->keys[left];
            queue->data[pos] = queue->data[left];
            queue->positions[pos] = queue->positions[left];

            if(queue->positions[pos] != NULL)
               *(queue->positions[pos]) = pos;

            pos = left;
         }
         else
            break;
      }

      left = 2 * pos + 1;
      right = 2 * pos + 2;
   }

   queue->keys[pos] = key;
   queue->data[pos] = elem;
   queue->positions[pos] = positionptr;

   if(positionptr != NULL)
      *positionptr = pos;

   return SCIP_OKAY;
}

/** Returns the value of a minimum element in 'elem'. Sets 'elem' to -1 if the queue is empty */
static SCIP_RETCODE pqueue_get_min(SCIP *scip, PQueue *queue, int *elem)
{
   *elem = -1;

   if(queue->size > 0)
   {
      *elem = queue->data[0];

      queue->size--;

      /* copy last element to front and repair heap structure */
      if(queue->size > 0)
      {
         queue->keys[0] = queue->keys[queue->size];
         queue->data[0] = queue->data[queue->size];
         queue->positions[0] = queue->positions[queue->size];

         if(queue->positions[0] != NULL)
            *(queue->positions[0]) = 0;

         SCIP_CALL( pqueue_increase_key(scip, queue, 0, queue->keys[0]) );
      }
   }

   return SCIP_OKAY;
}

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

/** Allocates memory for a lagrange multiplier and a set covering solution */
static SCIP_RETCODE allocateMemoryForSolution(SCIP *scip, SCP_Core *core, SCP_Lagrange_Sol *mult)
{
   int i;

   SCIP_CALL( SCIPallocBufferArray(scip, &mult->u, core->nconstraints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->subgradient, core->nconstraints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->lagrangian_costs_local, core->nvariables) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mult->lagrangian_costs_global, core->nvariables) );
   SCIP_CALL( SCIPhashtableCreate(&mult->x_greedy_local, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );

   for(i = 0; i < core->nvariables; i++)
   {
      mult->lagrangian_costs_local[i] = 0.0;
      mult->lagrangian_costs_global[i] = 0.0;
   }

   for(i = 0; i < core->nconstraints; i++)
   {
      mult->u[i] = 0.0;
      mult->subgradient[i] = 0.0;
   }

   mult->lb_lagrange_global = SCIP_REAL_MIN;
   mult->lb_lagrange_local = SCIP_REAL_MIN;
   mult->ub_greedy_local = SCIP_REAL_MAX;

   return SCIP_OKAY;
}

/** checks if 'variable' is a core variable */
static SCIP_Bool isCoreVariable(SCP_Core *core, SCIP_VAR *variable)
{
   assert(core != NULL);
   assert(variable != NULL);

   return SCIPhashtableExists(core->corevariables, variable);
}

/** checks if 'variable' is fixed within instance 'inst' */
static SCIP_Bool isFixedVariable(SCP_Instance *inst, SCIP_VAR *variable)
{
   assert(inst != NULL);
   assert(variable != NULL);

   return SCIPhashtableExists(inst->varsfixed, variable);
}

/** checks if 'variable' is part of 'solution' */
static SCIP_Bool isVarInSolution(SCIP_HASHTABLE *solution, SCIP_VAR *variable)
{
   assert(solution != NULL);
   assert(variable != NULL);

   return SCIPhashtableExists(solution, variable);
}

/** fixes 'variable' within instance 'inst' */
static void fixVariable(SCP_Instance *inst, SCIP_VAR *variable)
{
   assert(inst != NULL);
   assert(variable != NULL);

   SCIPhashtableInsert(inst->varsfixed, variable);
}

/** adds variable 'variable' to 'solution' */
static void addVariableToSolution(SCIP_HASHTABLE *solution, SCIP_VAR *variable)
{
   assert(solution != NULL);
   assert(variable);

   /* no test whether the variable is already part of the solution */
   SCIPhashtableInsert(solution, variable);
}

/** checks if the row at position 'rowpos' is covered by fixed variables of 'inst' */
static SCIP_Bool isRowCovered(SCP_Instance *inst, int rowpos)
{
   assert(inst != NULL);
   return SCIPhashtableExists(inst->rowscovered, (void *) ((size_t) rowpos + 1));
}

/** marks the row at position 'rowpos' as covered within instance 'inst' */
static void markRowAsCovered(SCP_Instance *inst, int rowpos)
{
   assert(inst != NULL);
   SCIPhashtableInsert(inst->rowscovered, (void *) ((size_t) rowpos + 1));
}

/** returns the position of 'variable' within the array core->variables */
static SCIP_RETCODE getVarIndex(SCIP *scip, SCP_Core *core, SCIP_VAR *variable, int *pos)
{
   int varidx;

   assert(scip != NULL);
   assert(core != NULL);
   assert(variable != NULL);
   assert(pos != NULL);

   varidx = SCIPvarGetIndex(variable);
   *pos = (int) (size_t) SCIPhashmapGetImage(core->mapvariables, (void *) (((size_t) varidx) + 1));

   return SCIP_OKAY;
}

/** get all variables that are part of the constraint at position 'pos' and saves them into 'vars' */
static SCIP_RETCODE getConsVars(SCIP *scip, SCP_Core *core, int pos, SCIP_VAR **vars, int *nvars, SCIP_Bool *success)
{
   *success = FALSE;
   if(SCIPconsIsActive(core->constraints[pos]) == FALSE)
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[pos], nvars, success) );
   if(*success == FALSE)
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetConsVars(scip, core->constraints[pos], vars, core->maxconstraintvariables, success) );

   return SCIP_OKAY;
}

/** Releases all memory of a lagrange multiplier */
static SCIP_RETCODE freeMemoryForSolution(SCIP *scip, SCP_Lagrange_Sol *mult)
{
   SCIPfreeBufferArray(scip, &mult->u);
   SCIPfreeBufferArray(scip, &mult->subgradient);
   SCIPfreeBufferArray(scip, &mult->lagrangian_costs_local);
   SCIPfreeBufferArray(scip, &mult->lagrangian_costs_global);
   SCIPhashtableFree(&mult->x_greedy_local);
   return SCIP_OKAY;
}

/** Creates a set covering solution. Adds all fixed variables of 'inst' and all variables of 'source' to 'dest'.
 *  'destcosts' contains the total costs of the solution.
 */
static SCIP_RETCODE copySetCoverSolution(SCP_Core *core, SCP_Instance *inst, SCIP_HASHTABLE *dest, SCIP_Real *destcosts, SCIP_HASHTABLE *source)
{
   int i;
   SCIP_Real costs = 0.0;

   /* remove all entries from dest */
   SCIPhashtableRemoveAll(dest);

   /* iterate through all variables of the problem */
   for(i = 0; i < core->nvariables; i++)
   {
      if(isVarInSolution(source, core->variables[i]) || isFixedVariable(inst, core->variables[i]))
      {
         /* add variable 'i' if it is in the solution 'source' or fixed in instance 'inst' */
         costs += SCIPvarGetObj(core->variables[i]);
         addVariableToSolution(dest, core->variables[i]);
      }
   }

   /* save total costs if required */
   if(destcosts != NULL)
      *destcosts = costs;

   return SCIP_OKAY;
}

/** Copies all data of the lagrange multiplier 'source' to the lagrange multiplier 'dest' */
static SCIP_RETCODE copySolution(SCP_Core *core, SCP_Lagrange_Sol *dest, SCP_Lagrange_Sol *source)
{
   int i;

   SCIPhashtableRemoveAll(dest->x_greedy_local);

   memcpy(dest->lagrangian_costs_local, source->lagrangian_costs_local, sizeof(*dest->lagrangian_costs_local) * core->nvariables);
   memcpy(dest->lagrangian_costs_global, source->lagrangian_costs_global, sizeof(*dest->lagrangian_costs_global) * core->nvariables);
   memcpy(dest->u, source->u, sizeof(*dest->u) * core->nconstraints);
   memcpy(dest->subgradient, source->subgradient, sizeof(*dest->subgradient) * core->nconstraints);

   for(i = 0; i < core->nvariables; i++)
   {
      if(isVarInSolution(source->x_greedy_local, core->variables[i]))
         addVariableToSolution(dest->x_greedy_local, core->variables[i]);
   }

   dest->lb_lagrange_global = source->lb_lagrange_global;
   dest->lb_lagrange_local = source->lb_lagrange_local;
   dest->ub_greedy_local = source->ub_greedy_local;

   return SCIP_OKAY;
}

/** Initializes an instance where no variables are fixed and no rows are covered, yet */
static SCIP_RETCODE initInstance(SCIP *scip, SCP_Instance *inst)
{
   SCIP_CALL( SCIPhashtableCreate(&inst->varsfixed, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );
   SCIP_CALL( SCIPhashtableCreate(&inst->rowscovered, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyInt, hashKeyEqInt, hashKeyValInt, NULL) );
   inst->costsfixed = 0.0;
   return SCIP_OKAY;
}

/** Copies the fixed variables from 'source' to 'dest', but 'dest' does not automatically contain all rows covered by these variables */
static SCIP_RETCODE copyInstance(SCIP *scip, SCP_Core *core, SCP_Instance *dest, SCP_Instance *source)
{
   int i;

   SCIPhashtableRemoveAll(dest->varsfixed);
   SCIPhashtableRemoveAll(dest->rowscovered);
   dest->costsfixed = 0.0;

   for(i = 0; i < core->nvariables; i++)
   {
      if(isFixedVariable(source, core->variables[i]) == TRUE)
      {
         fixVariable(dest, core->variables[i]);
         dest->costsfixed += SCIPvarGetObj(core->variables[i]);
      }
   }

   return SCIP_OKAY;
}

/** Releases all memory used by an instance */
static SCIP_RETCODE freeInstance(SCIP *scip, SCP_Instance *inst)
{
   SCIPhashtableFree(&inst->varsfixed);
   SCIPhashtableFree(&inst->rowscovered);
   return SCIP_OKAY;
}

/** Initializes a tentative core: for each row the first few columns covering this row are added to the core */
static SCIP_RETCODE initTentativeCore(SCIP *scip, SCP_Core *core)
{
   SCIP_VAR **vars;
   SCIP_Bool success;
   int nvars;
   int i, j;

   assert(scip != NULL);
   assert(core != NULL);

   /* mapvariables is a mapping: SCIPvarGetIndex(core->variables[i]) -> i */
   SCIP_CALL( SCIPhashmapCreate(&core->mapvariables, SCIPblkmem(scip), SCIPcalcHashtableSize(SCIPgetNVars(scip))) );
   SCIP_CALL( SCIPhashtableCreate(&core->corevariables, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );

   core->rowsavailable = FALSE;
   core->rows = NULL;
   core->nrowvars = NULL;
   core->columnsavailable = FALSE;
   core->columns = NULL;
   core->nvariables = SCIPgetNVars(scip);
   core->variables = SCIPgetVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &core->delta, core->nvariables) );
   SCIP_CALL( SCIPallocBufferArray(scip, &core->delta_pos, core->nvariables) );
   SCIP_CALL( SCIPallocBufferArray(scip, &core->solgreedy, core->nvariables) );
   core->nsolgreedy = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &core->nvarconstraints, core->nvariables) );
   for(i = 0; i< core->nvariables; i++)
      core->nvarconstraints[i] = 0;

   /* construct mapping of variable-indices to array 'variables' */
   for(i = 0; i < core->nvariables; i++)
   {
      int varidx = SCIPvarGetIndex(core->variables[i]);
      SCIP_CALL( SCIPhashmapInsert(core->mapvariables, (void *) (((size_t) varidx) + 1), (void *) (size_t) i) );
   }

   core->nactiveconstraints = 0;
   core->maxconstraintvariables = 0;
   core->nconstraints = SCIPgetNConss(scip);
   core->constraints = SCIPgetConss(scip);
   core->ncorevariables = 0;
   core->listcorevariables = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->nvariables) );

   for(i = 0; i < core->nconstraints; i++)
   {
      if(SCIPconsIsActive(core->constraints[i]) == FALSE)
      {
         SCIPdebugMessage("constraint %i (%s) is inactive\n", i, SCIPconsGetName(core->constraints[i]));
         continue;
      }

      /* get all variables that are part of this constraint */
      SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[i], &nvars, &success) );
      if(success == FALSE)
      {
         SCIPdebugMessage("constraint %i (%s): can't get number of variables\n", i, SCIPconsGetName(core->constraints[i]));
         continue;
      }

      SCIP_CALL( SCIPgetConsVars(scip, core->constraints[i], vars, core->nvariables, &success) );
      if(success == FALSE)
      {
         SCIPdebugMessage("constraint %i (%s): can't get variables\n", i, SCIPconsGetName(core->constraints[i]));
         continue;
      }

      if(nvars > core->maxconstraintvariables)
         core->maxconstraintvariables = nvars;

      for(j = 0; j < nvars; j++)
      {
         int varpos;

         SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );
         if(j < SCP_CORE_TENT_SIZE)
         {
            /* add this variable to the core if it's not already in there */
            if(!isCoreVariable(core, core->variables[varpos]))
               SCIPhashtableInsert(core->corevariables, core->variables[varpos]);
         }

         /* increase the number of constraints this variable is part of */
         core->nvarconstraints[varpos]++;
      }

      core->nactiveconstraints++;
   }

   /* create list of core variables, so it is easy to traverse them */
   SCIP_CALL( SCIPallocBufferArray(scip, &core->listcorevariables, core->ncorevariables) );
   j = 0;
   core->ncorevariables = (int) SCIPhashtableGetNElements(core->corevariables);
   for(i = 0; i < core->nvariables; i++)
   {
      if(isCoreVariable(core, core->variables[i]))
         core->listcorevariables[j++] = i;
   }

   SCIPfreeBufferArray(scip, &vars);
   SCIPdebugMessage("%lli variables in the tentative core\n", SCIPhashtableGetNElements(core->corevariables));

   return SCIP_OKAY;
}

/** Adds all fixed variables of 'inst' to a set covering solution 'solution' */
static SCIP_RETCODE extendSolution(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCIP_HASHTABLE *solution)
{
   int i;

   for(i = 0; i < core->nvariables; i++)
   {
      if(!isFixedVariable(inst, core->variables[i]))
         continue;

      /* add fixed variable if it is not already part of the solution */
      if(!isVarInSolution(solution, core->variables[i]))
         addVariableToSolution(solution, core->variables[i]);
   }

   return SCIP_OKAY;
}

/** Constructs rows of all constraints, but only includes core variables */
static SCIP_RETCODE computeCoreRows(SCIP *scip, SCP_Core *core, SCIP_HEURDATA *heurdata)
{
   SCIP_Bool success;
   SCIP_VAR **vars;
   int i;
   int j;
   int ncorevars;
   int nvars;

   assert(scip != NULL);
   assert(core != NULL);
   assert(heurdata != NULL);

   /* don't compute again if already computed */
   if(core->rowsavailable)
      return SCIP_OKAY;

   assert(core->rows == NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );
   SCIP_CALL( SCIPallocBufferArray(scip, &core->rows, core->nconstraints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &core->nrowvars, core->nconstraints) );

   /* iterate through list of constraints */
   for(i = 0; i < core->nconstraints; i++)
   {
      core->rows[i] = NULL;
      core->nrowvars[i] = 0;

      SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
      if(success == FALSE)
         continue;

      /* count number of core variables this constraint contains */
      ncorevars = 0;
      for(j = 0; j < nvars; j++)
      {
         int varpos;

         SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

         if(isCoreVariable(core, core->variables[varpos]))
            ncorevars++;
      }

      if(ncorevars > 0)
      {
         /* create array of core variables of this constraint */
         core->nrowvars[i] = ncorevars;
         SCIP_CALL( SCIPallocBufferArray(scip, &core->rows[i], core->nrowvars[i]) );
         for(j = 0; j < nvars; j++)
         {
            int varpos;

            SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

            if(isCoreVariable(core, core->variables[varpos]))
            {
               ncorevars--;
               core->rows[i][ncorevars] = varpos;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &vars);
   core->rowsavailable = TRUE;
   return SCIP_OKAY;
}

/** Constructs columns of core variables to provide better access */
static SCIP_RETCODE computeCoreColumns(SCIP *scip, SCP_Core *core, SCIP_HEURDATA *heurdata)
{
   SCIP_Bool success;
   SCIP_VAR **vars;
   int i, j, k;
   int nvars;

   assert(scip != NULL);
   assert(core != NULL);

   /* don't compute columns again if already computed */
   if(core->columnsavailable)
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &core->columns, core->nvariables) );
   for(i = 0; i < core->nvariables; i++)
   {
      /* columns[i] = NULL for all non-core variables */
      core->columns[i] = NULL;

      /* check if variable is a core variable */
      if(SCIPhashtableExists(core->corevariables, core->variables[i]) == FALSE)
         continue;

      /* only allocate memory of it is part of any constraints at all (this should always be the case for core variables) */
      if(core->nvarconstraints[i] == 0)
         continue;

      SCIP_CALL( SCIPallocBufferArray(scip, &core->columns[i], core->nvarconstraints[i]) );

      for(j = 0; j < core->nvarconstraints[i]; j++)
         core->columns[i][j] = -1;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );

   for(i = 0; i < core->nconstraints; i++)
   {
      SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
      if(success == FALSE)
         continue;

      for(j = 0; j < nvars; j++)
      {
         int varpos;

         SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

         if(!isCoreVariable(core, core->variables[varpos]))
            continue;

         /* add this constraint to the column of the variable */
         for(k = 0; k < core->nvarconstraints[varpos]; k++)
         {
            /* find first position that is not used yet in column 'varpos' */
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
   SCIPfreeBufferArray(scip, &core->delta);
   SCIPfreeBufferArray(scip, &core->delta_pos);
   SCIPfreeBufferArray(scip, &core->solgreedy);

   if(core->columnsavailable)
   {
      int i;

      for(i = 0; i < core->nvariables; i++)
      {
         if(core->columns[i])
            SCIPfreeBufferArray(scip, &core->columns[i]);
      }

      SCIPfreeBufferArray(scip, &core->columns);
   }

   if(core->rowsavailable)
   {
      int i;

      for(i = 0; i < core->nconstraints; i++)
      {
         if(core->rows[i])
            SCIPfreeBufferArray(scip, &core->rows[i]);
      }

      SCIPfreeBufferArray(scip, &core->rows);
      SCIPfreeBufferArray(scip, &core->nrowvars);
   }

   if(core->listcorevariables != NULL)
      SCIPfreeBufferArray(scip, &core->listcorevariables);

   return SCIP_OKAY;
}

static SCIP_RETCODE redefineCore(SCIP *scip, SCIP_HEURDATA *heurdata)
{
   SCIP_Bool recomputecolumns = FALSE;
   SCIP_Bool recomputerows = FALSE;
   SCP_Core *core;
   int *delta_perm;
   int i, j;
   int nvars;
   SCIP_VAR **vars;

   /* assumption: delta values were already computed and are sorted in increasing order */

   core = &heurdata->core;
   SCIP_CALL( SCIPallocBufferArray(scip, &delta_perm, core->nvariables) );

   for(i = 0; i < core->nvariables; i++)
   {
      delta_perm[core->delta_pos[i]] = i;
   }

   /* remove data about core variables */
   SCIPhashtableRemoveAll(core->corevariables);
   if(core->columnsavailable)
   {
      recomputecolumns = TRUE;

      for(i = 0; i < core->nvariables; i++)
      {
         if(core->columns[i] != NULL)
            SCIPfreeBufferArray(scip, &core->columns[i]);
      }

      SCIPfreeBufferArray(scip, &core->columns);
      core->columns = NULL;
      core->columnsavailable = FALSE;
   }

   if(core->rowsavailable)
   {
      recomputerows = TRUE;

      for(i = 0; i < core->nconstraints; i++)
      {
         if(core->rows[i])
            SCIPfreeBufferArray(scip, &core->rows[i]);
      }

      SCIPfreeBufferArray(scip, &core->rows);
      SCIPfreeBufferArray(scip, &core->nrowvars);

      core->rows = NULL;
      core->nrowvars = NULL;
      core->rowsavailable = FALSE;
   }

   if(core->listcorevariables != NULL)
   {
      SCIPfreeBufferArray(scip, &core->listcorevariables);
      core->listcorevariables = NULL;
   }

   /* pick the first 'SCP_CORE_TENT_SIZE'*m columns with lowest delta values to be in the core */
   for(i = 0; i < core->nvariables; i++)
   {
      if(i >= SCP_CORE_TENT_SIZE * core->nactiveconstraints)
         break;

      SCIPhashtableInsert(core->corevariables, core->variables[core->delta_pos[i]]);
   }

   vars = heurdata->vars;
   /*SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );*/

   /* then add the first 'SCP_CORE_TENT_SIZE' columns covering each row in increasing order of their delta values */
   for(i = 0; i < core->nconstraints; i++)
   {
      SCIP_Bool success = FALSE;
      int cols[SCP_CORE_TENT_SIZE];
      SCIP_Real coldelta[SCP_CORE_TENT_SIZE];

      SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
      if(success == FALSE)
         continue;

      for(j = 0; j < SCP_CORE_TENT_SIZE; j++)
         coldelta[j] = SCIP_REAL_MAX;

      for(j = 0; j < nvars; j++)
      {
         int varpos;
         int k;
         SCIP_Real value;

         SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

         value = core->delta[delta_perm[varpos]];

         if(j < SCP_CORE_TENT_SIZE)
            k = j;
         else
            k = SCP_CORE_TENT_SIZE - 1;

         if((j < SCP_CORE_TENT_SIZE) || (coldelta[k] > value))
         {
            cols[k] = varpos;
            coldelta[k] = value;

            while((k > 0) && (coldelta[k] < coldelta[k - 1]))
            {
               int tmp1 = cols[k - 1];
               SCIP_Real tmp2 = coldelta[k - 1];

               cols[k - 1] = cols[k];
               coldelta[k - 1] = coldelta[k];
               cols[k] = tmp1;
               coldelta[k] = tmp2;
               k--;
            }
         }
      }

      for(j = 0; j < nvars; j++)
      {
         if(j >= SCP_CORE_TENT_SIZE)
            break;

         if(!isCoreVariable(core, core->variables[cols[j]]))
            SCIPhashtableInsert(core->corevariables, core->variables[cols[j]]);
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &core->listcorevariables, SCIPhashtableGetNElements(core->corevariables)) );
   core->ncorevariables = 0;
   for(i = 0; i < core->nvariables; i++)
   {
      if(isCoreVariable(core, core->variables[i]))
         core->listcorevariables[core->ncorevariables++] = i;
   }

   if(recomputecolumns)
      SCIP_CALL( computeCoreColumns(scip, core, heurdata) );

   if(recomputerows)
      SCIP_CALL( computeCoreRows(scip, core, heurdata) );

   SCIPfreeBufferArray(scip, &delta_perm);
   /*SCIPfreeBufferArray(scip, &vars);*/

   SCIPdebugMessage("%lli variables are in the refined core\n", SCIPhashtableGetNElements(core->corevariables));
   return SCIP_OKAY;
}

/* adds all indices of rows to inst->rowscovered for all rows that are covered by the variables in inst->varsfixed */
static SCIP_RETCODE markRowsCoveredByFixedVariables(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCIP_HEURDATA *heurdata)
{
   int i, j;

   SCIPhashtableRemoveAll(inst->rowscovered);

   if(core->columnsavailable == TRUE)
   {
      for(i = 0; i < core->ncorevariables; i++)
      {
         int corevar = core->listcorevariables[i];

         if(!isFixedVariable(inst, core->variables[corevar]))
            continue;

         for(j = 0; j < core->nvarconstraints[corevar]; j++)
         {
            int rowidx = core->columns[corevar][j];

            if(!isRowCovered(inst, rowidx))
               markRowAsCovered(inst, rowidx);
         }
      }
   }
   else
   {
      SCIP_VAR **vars;
      int nvars;

      /*SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );*/
      vars = heurdata->vars;

      for(i = 0; i < core->nconstraints; i++)
      {
         SCIP_Bool success = FALSE;

         SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
         if(success == FALSE)
            continue;

         for(j = 0; j < nvars; j++)
         {
            if(isFixedVariable(inst, vars[j]))
            {
               if(!isRowCovered(inst, i))
                  markRowAsCovered(inst, i);

               break;
            }
         }
      }

      /*SCIPfreeBufferArray(scip, &vars);*/
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE checkSetCover(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCIP_HASHTABLE *solution, SCIP_Bool *issetcover, SCIP_HEURDATA *heurdata)
{
   SCIP_Bool success;
   SCIP_VAR **vars;
   int nvars;
   int i, j;

   assert(scip != NULL);
   assert(core != NULL);
   assert(inst != NULL);
   assert(solution != NULL);
   assert(issetcover != NULL);
   assert(heurdata != NULL);

   *issetcover = TRUE;

   /*SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );*/
   vars = heurdata->vars;

   /* iterate through all constraints and check whether each of them contains a variable that is part of the cover */
   for(i = 0; i < core->nconstraints; i++)
   {
      SCIP_Bool rowcovered = FALSE;

      SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
      if(success == FALSE)
         continue;

      for(j = 0; j < nvars; j++)
      {
         int varpos;

         SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

         if(isVarInSolution(solution, core->variables[varpos]))
         {
            rowcovered = TRUE;
            break;
         }
      }

      if(rowcovered == FALSE)
      {
         SCIPdebugMessage("check set cover: row %i is not covered by any column\n", i);
         *issetcover = FALSE;
         break;
      }
   }

   /*SCIPfreeBufferArray(scip, &vars);*/

   return SCIP_OKAY;
}

static SCIP_RETCODE computeDelta(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCIP_Real *lagrangiancosts, SCIP_HASHTABLE *solution, SCIP_Real pi, SCIP_HEURDATA *heurdata)
{
   SCIP_Bool success;
   SCIP_Real delta_sum, delta_max;
   int *nvarcovering;
   SCIP_VAR **vars;
   int nvars;
   int i, j;
   int *delta_pos;
   SCIP_Real *delta;

   delta = core->delta;
   delta_pos = core->delta_pos;

   SCIP_CALL( SCIPallocBufferArray(scip, &nvarcovering, core->nconstraints) );
   vars = heurdata->vars;

   /* compute nvarcovering[i] = 'number of columns covering row i' */
   if(core->rowsavailable)
   {
      for(i = 0; i < core->nconstraints; i++)
      {
         nvarcovering[i] = 0;

         if(!SCIPconsIsActive(core->constraints[i]))
            continue;

         for(j = 0; j < core->nrowvars[i]; j++)
         {
            int varpos = core->rows[i][j];
            if(isVarInSolution(solution, core->variables[varpos]))
               nvarcovering[i]++;
         }
      }
   }
   else
   {
      for(i = 0; i < core->nconstraints; i++)
      {
         nvarcovering[i] = 0;

         SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
         if(success == FALSE)
            continue;

         for(j = 0; j < nvars; j++)
         {
            int varpos;

            SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

            if(isVarInSolution(solution, core->variables[varpos]))
               nvarcovering[i]++;
         }
      }
   }

   for(i = 0; i < core->nvariables; i++)
   {
      delta[i] = SCIP_REAL_MAX;
      delta_pos[i] = i;

      /* skip variables that are not part of the set covering solution */
      if(!isVarInSolution(solution, core->variables[i]))
         continue;

      delta[i] = lagrangiancosts[i];
      if(delta[i] < 0.0)
         delta[i] = 0.0;

      for(j = 0; j < core->nvarconstraints[i]; j++)
      {
         int rowpos = core->columns[i][j];

         if(isRowCovered(inst, rowpos))
            continue;

         delta[i] += lagrangiancosts[rowpos] * (nvarcovering[rowpos] - 1) / ((SCIP_Real) nvarcovering[rowpos]);
      }
   }

   SCIPsortRealInt(delta, delta_pos, core->nvariables);

   delta_sum = 0.0;
   delta_max = core->nactiveconstraints * pi;

   /* fix new variables of this instance */
   SCIPhashtableRemoveAll(inst->varsfixed);
   inst->costsfixed = 0.0;

   for(i = 0; i < core->nvariables; i++)
   {
      int varpos = delta_pos[i];

      if(!isVarInSolution(solution, core->variables[varpos]))
         break;

      inst->costsfixed += SCIPvarGetObj(core->variables[varpos]);
      delta_sum += delta[i];

      /* fix variable delta_pos[i] */
      fixVariable(inst, core->variables[varpos]);

      if(delta_sum >= delta_max)
         break;
   }

   /*SCIPfreeBufferArray(scip, &vars);*/
   SCIPfreeBufferArray(scip, &nvarcovering);
   return SCIP_OKAY;
}

static SCIP_RETCODE removeRedundantColumns(SCIP *scip, SCIP_HEURDATA *heurdata, SCIP_HASHTABLE *solution, SCIP_Real *solcosts)
{
   SCP_Core *core;
   SCIP_Real *costs;
   int *varpos;
   int *nvarcovering;
   SCIP_VAR **vars;
   int nvars;
   int i, j;
   SCIP_Bool success;
   int solsize;

   core = &heurdata->core;

   if(core->columnsavailable == FALSE)
   {
      SCIPdebugMessage("can only remove redundant columns if they are available in the core\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &costs, core->nvariables) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varpos, core->nvariables) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nvarcovering, core->nconstraints) );
   /*SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );*/
   vars = heurdata->vars;

   /* compute nvarcovering[i] = 'number of columns covering row i' */
   if(core->rowsavailable)
   {
      for(i = 0; i < core->nconstraints; i++)
      {
         nvarcovering[i] = 0;

         if(SCIPconsIsActive(core->constraints[i]) == FALSE)
            continue;

         if(core->nrowvars[i] == 0)
            continue;

         for(j = 0; j < core->nrowvars[i]; j++)
         {
            int vpos = core->rows[i][j];

            if(isVarInSolution(solution, core->variables[vpos]))
               nvarcovering[i]++;
         }
      }
   }
   else
   {
      for(i = 0; i < core->nconstraints; i++)
      {
         nvarcovering[i] = 0;

         SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
         if(success == FALSE)
            continue;

         for(j = 0; j < nvars; j++)
         {
            int vpos;

            SCIP_CALL( getVarIndex(scip, core, vars[j], &vpos) );

            if(isVarInSolution(solution, core->variables[vpos]))
               nvarcovering[i]++;
         }
      }
   }

   solsize = 0;
   for(i = 0; i < core->nvariables; i++)
   {
      if(isVarInSolution(solution, core->variables[i]))
      {
         costs[solsize] = -SCIPvarGetObj(core->variables[i]);
         varpos[solsize] = i;
         solsize++;
      }
   }

   SCIPsortRealInt(costs, varpos, solsize);

   if(core->columnsavailable == FALSE)
   {
      SCIPerrorMessage("columns need to be available in order to remove redundant columns\n");
      SCIPABORT();
   }

   for(i = 0; i < solsize; i++)
   {
      int vpos = varpos[i];
      for(j = 0; j < core->nvarconstraints[vpos]; j++)
      {
         if(nvarcovering[core->columns[vpos][j]] == 1)
            break;
      }

      if(j == core->nvarconstraints[vpos])
      {
         SCIPhashtableRemove(solution, core->variables[vpos]);
         *solcosts -= SCIPvarGetObj(core->variables[vpos]);
         /*SCIPdebugMessage("removing redundant column\n");*/

         for(j = 0; j < core->nvarconstraints[vpos]; j++)
         {
            nvarcovering[core->columns[vpos][j]]--;
         }
      }
   }

   SCIPfreeBufferArray(scip, &costs);
   SCIPfreeBufferArray(scip, &varpos);
   /*SCIPfreeBufferArray(scip, &vars);*/
   SCIPfreeBufferArray(scip, &nvarcovering);
   return SCIP_OKAY;
}

static SCIP_RETCODE greedySetCover(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *mult, SCIP_HEURDATA *heurdata)
{
   SCP_Instance *greedyinst;
   int i;
   int j;
   int nrowsuncovered = 0;
   SCIP_Bool success = FALSE;
   int nvars;

   PQueue *prioqueue;
   int *colpos;
   int *colmu;
   SCIP_Real *colgamma;
   SCIP_Real *colscore;
   SCIP_VAR **vars;
   int k;

   core->nsolgreedy = 0;
   mult->ub_greedy_local = 0.0;
   SCIPhashtableRemoveAll(mult->x_greedy_local);

   greedyinst = &heurdata->greedy_inst;
   SCIPhashtableRemoveAll(greedyinst->rowscovered);
   /*SCIP_CALL( SCIPhashtableCreate(&rowscovered, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyInt, hashKeyEqInt, hashKeyValInt, NULL) );*/

   for(i = 0; i < core->nconstraints; i++)
   {
      if(SCIPconsIsActive(core->constraints[i]) == FALSE)
         continue;

      /* this is actually necessary because there exist constraints where this fails, and we simply need to ignore them */
      SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[i], &nvars, &success) );
      if(success == FALSE)
         continue;

      if(isRowCovered(inst, i))
         markRowAsCovered(greedyinst, i);
      else
         nrowsuncovered++;
   }

   if(core->columnsavailable == FALSE)
   {
      SCIPerrorMessage("greedy algorithm requires core columns to be available\n");
      SCIPABORT();
   }

   /* compute scores and add them to the priority queue */
   colpos = heurdata->greedy_colpos;
   colmu = heurdata->greedy_colmu;
   colgamma = heurdata->greedy_colgamma;
   colscore = heurdata->greedy_colscore;
   vars = heurdata->vars;

   prioqueue = &heurdata->greedy_queue;
   SCIP_CALL( pqueue_clear(scip, prioqueue) );

   for(i = 0; i < core->nvariables; i++)
   {
      int varpos = i;
      int mu = 0;
      SCIP_Real gamma = SCIPvarGetObj(core->variables[varpos]);

      colmu[i] = 0;
      colgamma[i] = 0;
      colpos[i] = 0;
      colscore[i] = 0.0;

      if(!isCoreVariable(core, core->variables[varpos]))
          continue;
      else if(isFixedVariable(inst, core->variables[varpos]))
         continue;

      for(j = 0; j < core->nvarconstraints[varpos]; j++)
      {
         if(!isRowCovered(greedyinst, core->columns[varpos][j]))
         {
            gamma -= mult->u[core->columns[varpos][j]];
            mu++;
         }
      }

      /* skip columns that do not cover anything */
      if(mu > 0)
      {
         SCIP_Real score = (gamma > 0) ? (gamma / ((SCIP_Real) mu)) : (gamma * mu);
         colmu[i] = mu;
         colgamma[i] = gamma;
         colscore[i] = score;

         SCIP_CALL( pqueue_insert(scip, prioqueue, colscore[i], i, &(colpos[i])) );
      }
   }

   while(nrowsuncovered > 0)
   {
      int mincolumn;
      SCIP_CALL( pqueue_get_min(scip, prioqueue, &mincolumn) );

      /* add variable 'variables[mincolumn]' to the set cover */
      SCIPhashtableSafeInsert(mult->x_greedy_local, core->variables[mincolumn]);
      mult->ub_greedy_local += SCIPvarGetObj(core->variables[mincolumn]);
      core->solgreedy[core->nsolgreedy++] = mincolumn;

      /*SCIPdebugMessage("adding %i to the set cover, score = %f\n", mincolumn, colscore[mincolumn]);*/
      colmu[mincolumn] = 0;

      for(j = 0; j < core->nvarconstraints[mincolumn]; j++)
      {
         int columnpos = core->columns[mincolumn][j];

         if(!isRowCovered(greedyinst, columnpos))
         {
            markRowAsCovered(greedyinst, columnpos);
            nrowsuncovered--;

            /* update scores of columns covering this row */

            SCIP_CALL( getConsVars(scip, core, columnpos, vars, &nvars, &success) );
            if(success == FALSE)
               continue;

            /* for each core variable: subtract u[i] from the variable's costs */
            for(k = 0; k < nvars; k++)
            {
               SCIP_Real score;
               int varpos;

               SCIP_CALL( getVarIndex(scip, core, vars[k], &varpos) );

               /* skip non-core variables */
               if(!isCoreVariable(core, core->variables[varpos]))
                  continue;

               if(colmu[varpos] == 0)
                  continue;

               score = colscore[varpos];
               colscore[varpos] = SCIP_REAL_MAX;
               colmu[varpos]--;
               colgamma[varpos] += mult->u[columnpos];

               if(colmu[varpos] > 0)
               {

                  colscore[varpos] = (colgamma[varpos] > 0) ? (colgamma[varpos] / ((SCIP_Real) colmu[varpos])) : (colgamma[varpos] * colmu[varpos]);
               }

               if(score > colscore[varpos])
                  SCIP_CALL( pqueue_decrease_key(scip, prioqueue, colpos[varpos], colscore[varpos]) );
               else
                  SCIP_CALL( pqueue_increase_key(scip, prioqueue, colpos[varpos], colscore[varpos]) );
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** computes lagrangian costs for all columns, only considering rows that are uncovered by fixed variables in 'inst' */
static SCIP_RETCODE computeLocalLagrangianCosts(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *mult, SCIP_HEURDATA *heurdata)
{
   SCIP_VAR **vars;
   int nvars;
   int i, j;

   /* set all lagrangian costs to objective values */
   for(i = 0; i < core->nvariables; i++)
      mult->lagrangian_costs_local[i] = SCIPvarGetObj(core->variables[i]);

   vars = heurdata->vars;

   if(core->rowsavailable)
   {
      for(i = 0; i < core->nconstraints; i++)
      {
         if(SCIPconsIsActive(core->constraints[i]) == FALSE)
            continue;

         if(isRowCovered(inst, i))
            continue;

         for(j = 0; j < core->nrowvars[i]; j++)
         {
            int varpos = core->rows[i][j];
            mult->lagrangian_costs_local[varpos] -= mult->u[i];
         }
      }
   }
   else
   {
      for(i = 0; i < core->nconstraints; i++)
      {
         SCIP_Bool success;

         /* skip rows that are not part of the reduced instance */
         if(isRowCovered(inst, i))
            continue;

         SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
         if(success == FALSE)
            continue;

         /* for each core variable: subtract u[i] from the variable's costs */
         for(j = 0; j < nvars; j++)
         {
            int varpos;

            SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

            mult->lagrangian_costs_local[varpos] -= mult->u[i];
         }
      }
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE computeGlobalLagrangianCosts(SCIP *scip, SCP_Core *core, SCP_Lagrange_Sol *mult, SCIP_HEURDATA *heurdata)
{
   SCIP_VAR **vars;
   int nvars;
   int i, j;

   /* set all lagrangian costs to objective values */
   for(i = 0; i < core->nvariables; i++)
      mult->lagrangian_costs_global[i] = SCIPvarGetObj(core->variables[i]);

   /*SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );*/
   vars = heurdata->vars;

   for(i = 0; i < core->nconstraints; i++)
   {
      SCIP_Bool success = FALSE;

      SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
      if(success == FALSE)
         continue;

      /* for each core variable: subtract u[i] from the variable's costs */
      for(j = 0; j < nvars; j++)
      {
         int varpos;

         SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

         mult->lagrangian_costs_global[varpos] -= mult->u[i];
      }

      mult->lb_lagrange_global += mult->u[i];
   }

   for(i = 0; i < core->nvariables; i++)
   {
      if(mult->lagrangian_costs_global[i] < 0.0)
         mult->lb_lagrange_global += mult->lagrangian_costs_global[i];
   }

   /*SCIPfreeBufferArray(scip, &vars);*/

   return SCIP_OKAY;
}

/* computes an optimal solution to the lagrangian relaxation, see formulae (4), (5) in the paper */
static SCIP_RETCODE computeOptimalSolution(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *mult, SCIP_HEURDATA *heurdata)
{
   SCIP_VAR **vars;
   SCIP_Bool success;
   int nvars;
   int i, j;

   mult->lb_lagrange_local = 0.0;
   mult->lb_lagrange_global = 0.0;

   SCIP_CALL( computeLocalLagrangianCosts(scip, core, inst, mult, heurdata) );
   SCIP_CALL( computeGlobalLagrangianCosts(scip, core, mult, heurdata) );

   for(i = 0; i < core->ncorevariables; i++)
   {
      int varpos = core->listcorevariables[i];

      /*if(SCIPhashtableExists(core->corevariables, core->variables[varpos]) == FALSE)
         continue;*/

      if(mult->lagrangian_costs_local[varpos] < 0.0)
      {
         if(!isFixedVariable(inst, core->variables[varpos]))
            mult->lb_lagrange_local = mult->lb_lagrange_local + mult->lagrangian_costs_local[varpos];
      }
   }

   if(core->rowsavailable)
   {
      for(i = 0; i < core->nconstraints; i++)
      {
         mult->subgradient[i] = 0.0;

         if(SCIPconsIsActive(core->constraints[i]) == FALSE)
            continue;

         if(isRowCovered(inst, i))
            continue;

         if(core->nrowvars[i] == 0)
            continue;

         mult->subgradient[i] = 1.0;

         for(j = 0; j < core->nrowvars[i]; j++)
         {
            int varpos = core->rows[i][j];

            if(mult->lagrangian_costs_local[varpos] < 0.0)
               mult->subgradient[i] -= 1.0;
         }

         mult->lb_lagrange_local = mult->lb_lagrange_local + mult->u[i];
      }
   }
   else
   {
      /*SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );*/
      vars = heurdata->vars;

      for(i = 0; i < core->nconstraints; i++)
      {
         mult->subgradient[i] = 0.0;

         if(isRowCovered(inst, i))
            continue;

         SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
         if(success == FALSE)
            continue;

         mult->subgradient[i] = 1.0;

         for(j = 0; j < nvars; j++)
         {
            int varpos;

            SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

            if(!isCoreVariable(core, core->variables[varpos]))
               continue;

            if(mult->lagrangian_costs_local[varpos] < 0.0)
               mult->subgradient[i] -= 1.0;
         }

         mult->lb_lagrange_local = mult->lb_lagrange_local + mult->u[i];
      }

      /*SCIPfreeBufferArray(scip, &vars);*/
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE subgradientOptimization(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *best_mult_lb, SCIP_Real best_ub, SCIP_HEURDATA *heurdata)
{
   int max_iter;
   unsigned int seed;
   SCP_Lagrange_Sol last_mult, next_mult, tmp_mult;
   SCIP_Real norm, lambda = SCP_LAMBDA;
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
      if(!isRowCovered(inst, i))
         last_mult.u[i] = SCIPgetRandomReal(0.9, 1.1, &seed) * last_mult.u[i];
      else
         last_mult.u[i] = 0.0;
   }

   best_ub -= inst->costsfixed;

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
         SCIP_Real hk = last_mult.u[i] + lambda * (best_ub - last_mult.lb_lagrange_local) * last_mult.subgradient[i] / norm;
         next_mult.u[i] = 0.0;

         if(hk > 0.0)
            next_mult.u[i] = hk;
      }

      SCIP_CALL( computeOptimalSolution(scip, core, inst, &next_mult, heurdata) );

      /*SCIP_CALL( greedySetCover(scip, core, inst, &next_mult) );

      if(next_mult.ub_greedy_local + inst->costsfixed < heurdata->best_ub)
      {
         SCIP_CALL( copySetCoverSolution(core, inst, heurdata->best_ub_sol, &heurdata->best_ub, next_mult.x_greedy_local) );
         SCIP_CALL( removeRedundantColumns(scip, heurdata, heurdata->best_ub_sol, &heurdata->best_ub) );
         SCIPdebugMessage("new upper bound: %f\n", heurdata->best_ub);
      }*/

      if(next_mult.lb_lagrange_local > best_mult_lb->lb_lagrange_local)
         SCIP_CALL( copySolution(core, best_mult_lb, &next_mult) );

      if(next_mult.lb_lagrange_global > heurdata->mult_best_lb_total.lb_lagrange_global)
      {
         /*SCIP_CALL( greedySetCover(scip, core, inst, &next_mult) );
         SCIP_CALL( extendSolution(scip, core, inst, next_mult.x_greedy_local) );*/
         SCIP_CALL( copySolution(core, &heurdata->mult_best_lb_total, &next_mult) );
      }

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

static SCIP_RETCODE computeInitialLagrangeMultiplier(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCP_Lagrange_Sol *mult, SCIP_HEURDATA *heurdata)
{

   int i, j;

   if(core->columnsavailable == TRUE)
   {
      for(i = 0; i < core->nconstraints; i++)
         mult->u[i] = SCIP_REAL_MAX;

      for(i = 0; i < core->nvariables; i++)
      {
         int nuncovered = 0;

         if(!isCoreVariable(core, core->variables[i]))
            continue;

         if(isFixedVariable(inst, core->variables[i]))
            continue;

         /* count how many uncovered, active rows this column covers */
         for(j = 0; j < core->nvarconstraints[i]; j++)
         {
            int colpos = core->columns[i][j];

            /* skip inactive constraints */
            if(SCIPconsIsActive(core->constraints[colpos]) == FALSE)
               continue;

            /* skip covered rows */
            if(isRowCovered(inst, colpos))
               continue;

            nuncovered++;
         }

         /* if this column covers any rows, update their cost if necessary */
         if(nuncovered > 0)
         {
            SCIP_Real costs = SCIPvarGetObj(core->variables[i]) / ((SCIP_Real) nuncovered);
            int nvars;

            for(j = 0; j < core->nvarconstraints[i]; j++)
            {
               int colpos = core->columns[i][j];
               SCIP_Bool success = FALSE;

               if(isRowCovered(inst, colpos))
                  continue;

               if(SCIPconsIsActive(core->constraints[colpos]) == FALSE)
                  continue;

               SCIP_CALL( SCIPgetConsNVars(scip, core->constraints[colpos], &nvars, &success) );
               if(success == FALSE)
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

      /*SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );*/
      vars = heurdata->vars;
      SCIP_CALL( SCIPallocBufferArray(scip, &nuncoveredactive, core->nvariables) );

      for(i = 0; i < core->nvariables; i++)
         nuncoveredactive[i] = 0;

      for(i = 0; i < core->nconstraints; i++)
      {
         SCIP_Bool success = FALSE;

         if(isRowCovered(inst, i))
            continue;

         SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
         if(success == FALSE)
            continue;

         for(j = 0; j < nvars; j++)
         {
            int varpos;

            SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

            if(!isCoreVariable(core, vars[j]))
               continue;
            else
               nuncoveredactive[varpos]++;
         }
      }

      for(i = 0; i < core->nconstraints; i++)
      {
         SCIP_Bool found = FALSE;
         SCIP_Bool success = FALSE;

         if(isRowCovered(inst, i))
         {
            mult->u[i] = 0.0;
            continue;
         }

         SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
         if(success == FALSE)
            continue;

         for(j = 0; j < nvars; j++)
         {
            if(!isCoreVariable(core, vars[j]))
               continue;
            else
            {
               int varpos;
               SCIP_Real costs = SCIP_REAL_MAX;

               SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

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

      /*SCIPfreeBufferArray(scip, &vars);*/
      SCIPfreeBufferArray(scip, &nuncoveredactive);
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE exploreNeighborhood(SCIP *scip, SCP_Lagrange_Sol *startmult, SCIP_HEURDATA *heurdata)
{
   SCP_Core *core;
   SCP_Instance *subinst;
   SCP_Lagrange_Sol mult;
   SCIP_VAR **vars;
   SCIP_Bool success;
   SCIP_Real best_ub;
   SCIP_Real norm;
   SCIP_Real costs;
   int nvars;
   int iter;
   int i;
   int j;

   core = &heurdata->core;
   subinst = &heurdata->subinst;
   best_ub = heurdata->best_ub - subinst->costsfixed;

   /*SCIP_CALL( SCIPallocBufferArray(scip, &vars, core->maxconstraintvariables) );*/
   vars = heurdata->vars;
   SCIP_CALL( allocateMemoryForSolution(scip, &heurdata->core, &mult) );
   SCIP_CALL( copySolution(core, &mult, startmult) );

   for(iter = 0; iter < SCP_HEUR_MAX_ITER; iter++)
   {
      norm = 0.0;

      /* compute subgradient for 'mult' */
      if(core->rowsavailable)
      {
         for(i = 0; i < core->nconstraints; i++)
         {
            mult.subgradient[i] = 0.0;

            if(SCIPconsIsActive(core->constraints[i]) == FALSE)
               continue;

            if(isRowCovered(subinst, i))
               continue;

            if(core->nrowvars[i] == 0)
               continue;

            mult.subgradient[i] = 1.0;

            for(j = 0; j < core->nrowvars[i]; j++)
            {
               int varpos = core->rows[i][j];

               if(mult.lagrangian_costs_local[varpos] < 0.0)
                  mult.subgradient[i] -= 1.0;
            }
         }
      }
      else
      {
         for(i = 0; i < core->nconstraints; i++)
         {
            mult.subgradient[i] = 0.0;

            if(isRowCovered(subinst, i))
               continue;

            SCIP_CALL( getConsVars(scip, core, i, vars, &nvars, &success) );
            if(success == FALSE)
               continue;

            mult.subgradient[i] = 1.0;

            for(j = 0; j < nvars; j++)
            {
               int varpos;

               SCIP_CALL( getVarIndex(scip, core, vars[j], &varpos) );

               if(!isCoreVariable(core, core->variables[varpos]))
                  continue;

               if(mult.lagrangian_costs_local[varpos] < 0.0)
                  mult.subgradient[i] -= 1.0;
            }
         }
      }

      /* compute norm of subgradient */
      norm = 0.0;
      for(i = 0; i < core->nconstraints; i++)
         norm = norm + mult.subgradient[i] * mult.subgradient[i];

      /* Held-Karp update */
      for(i = 0; i < core->nconstraints; i++)
      {
         SCIP_Real hk = mult.u[i] + SCP_LAMBDA * (best_ub - mult.lb_lagrange_local) * mult.subgradient[i] / norm;
         mult.u[i] = 0.0;

         if(hk > 0.0)
            mult.u[i] = hk;
      }

      SCIP_CALL( computeOptimalSolution(scip, core, subinst, &mult, heurdata) );

      if(mult.lb_lagrange_local > heurdata->mult_best_lb_subinst.lb_lagrange_local)
         SCIP_CALL( copySolution(core, &heurdata->mult_best_lb_subinst, &mult) );

      if(mult.lb_lagrange_global > heurdata->mult_best_lb_total.lb_lagrange_global)
         SCIP_CALL( copySolution(core, &heurdata->mult_best_lb_total, &mult) );

      SCIP_CALL( greedySetCover(scip, core, subinst, &mult, heurdata) );
      SCIP_CALL( extendSolution(scip, core, subinst, mult.x_greedy_local) );

      costs = subinst->costsfixed + mult.ub_greedy_local;
      SCIP_CALL( removeRedundantColumns(scip, heurdata, mult.x_greedy_local, &costs) );

      if(costs < heurdata->best_ub_subinst)
         SCIP_CALL( copySetCoverSolution(core, subinst, heurdata->best_ub_subinst_sol, &heurdata->best_ub_subinst, mult.x_greedy_local) );
   }

   /*SCIPfreeBufferArray(scip, &vars);*/
   SCIP_CALL( freeMemoryForSolution(scip, &mult) );

   return SCIP_OKAY;
}

static SCIP_RETCODE reportSolution(SCIP *scip, SCP_Core *core, SCP_Instance *inst, SCIP_HASHTABLE *solution, SCIP_HEUR *heur)
{
   SCIP_VAR **solvars;
   SCIP_Real *solvals;
   SCIP_Real newcosts = 0.0;
   int nsolvars;
   int i;
   SCIP_SOL *newsol;
   /*SCIP_SOL *bestsol;*/
   SCIP_Bool success = FALSE;
   SCIP_Bool foundsol = FALSE;
   SCIP_CONS **cons;
   int ncons;

   SCIP_CALL( SCIPgetVarsData(scip, &solvars, &nsolvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nsolvars) );

   for(i = 0; i < nsolvars; i++)
   {
      if(SCIPhashtableExists(solution, solvars[i]) == TRUE)
         solvals[i] = 1.0;
      else if((isFixedVariable(inst, solvars[i]) == TRUE) && (SCIPisZero(scip, SCIPvarGetObj(solvars[i])) == FALSE))
         solvals[i] = 1.0;
      else
         solvals[i] = 0.0;

      if(solvals[i] > 0.0)
         newcosts += SCIPvarGetObj(solvars[i]);
   }

   /* test if best primal solution is better than this solution */
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nsolvars, solvars, solvals) );
   /*bestsol = SCIPgetBestSol(scip);
   if(bestsol)
   {
      if(SCIPsolGetOrigObj(bestsol) < newcosts)
      {
         SCIPdebugMessage("best primal solution is still better than the set cover solution\n");
         SCIPfreeSol(scip, &newsol);
         SCIPfreeBufferArray(scip, &solvals);
         return SCIP_OKAY;
      }
   }*/

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
            SCIP_Bool costszero = SCIPisZero(scip, SCIPvarGetObj(vars[j]));
            /* SCIP_Bool changevar = FALSE; */

            /*SCIPdebugMessage("before var %i: %f\n", j, SCIPgetSolVal(scip, newsol, vars[j]));*/
            if(costszero == TRUE)
            {
               if(costszero == TRUE)
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
            /*else
            {
               int norigvars = GCGmasterVarGetNOrigvars(vars[j]);
               SCIP_VAR **origvars = GCGmasterVarGetOrigvars(vars[j]);

               for(k = 0; k < norigvars; k++)
               {
                  if(SCIPisZero(scip, SCIPvarGetObj(origvars[k])))
                     SCIPdebugMessage("var %i, original var %i has zero costs\n", j, k);
               }
            }*/
            /*SCIPdebugMessage("after var %i: %f\n", j, SCIPgetSolVal(scip, newsol, vars[j]));*/
         }

         if(changed == FALSE)
            SCIPdebugMessage("could not find variable with zero costs\n");

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

static SCIP_RETCODE threePhase(SCIP *scip, SCIP_HEUR *heur, SCIP_HEURDATA *heurdata)
{
   SCP_Lagrange_Sol *mult_lb_subinst;
   int nfixvars;
   int i;
   SCIP_Bool success;
   SCP_Core *core;
   SCP_Instance *inst;
   SCP_Instance *subinst;

   core = &heurdata->core;
   inst = &heurdata->inst;
   subinst = &heurdata->subinst;
   mult_lb_subinst = &heurdata->tp_mult_lb_subinst;

   nfixvars = core->nconstraints / 200;
   if(nfixvars < 1)
      nfixvars = 1;

   /* we first create our own copy of the instance, as we need to mark variables as fixed until all variables are fixed */
   SCIP_CALL( copyInstance(scip, core, subinst, inst) );
   SCIP_CALL( markRowsCoveredByFixedVariables(scip, core, subinst, heurdata) );

   /* next, compute initial lagrange multipliers and find a first lower bound */
   SCIP_CALL( computeInitialLagrangeMultiplier(scip, core, subinst, mult_lb_subinst, heurdata) );

   /* computeOptimalSolution also computes the subgradient */
   SCIP_CALL( computeOptimalSolution(scip, core, subinst, mult_lb_subinst, heurdata) );
   SCIP_CALL( greedySetCover(scip, core, subinst, mult_lb_subinst, heurdata) );

   /* we now have a lower and upper bound in mult_lb_local for the instance 'inst' and take these as our starting values */
   SCIP_CALL( copySetCoverSolution(core, inst, heurdata->best_ub_inst_sol, &heurdata->best_ub_inst, mult_lb_subinst->x_greedy_local) );
   SCIP_CALL( copySetCoverSolution(core, subinst, heurdata->best_ub_subinst_sol, &heurdata->best_ub_subinst, mult_lb_subinst->x_greedy_local) );
   SCIP_CALL( copySolution(core, &heurdata->mult_best_lb_inst, mult_lb_subinst) );

   /* check whether 'best_ub_inst_sol' is a solution of the reduced instance 'subinst' */
   SCIP_CALL( checkSetCover(scip, core, subinst, heurdata->best_ub_inst_sol, &success, heurdata) );
   if(success == FALSE)
   {
      SCIPerrorMessage("three-phase: initial solution is not a valid set cover\n");
      SCIPABORT();
   }

   if(heurdata->best_ub_inst < heurdata->best_ub)
   {
      SCIP_CALL( copySetCoverSolution(core, inst, heurdata->best_ub_sol, &heurdata->best_ub, heurdata->best_ub_inst_sol) );
      SCIPdebugMessage("new upper bound: %f\n", heurdata->best_ub);
   }

   if(core->nactiveconstraints <= SCIPhashtableGetNElements(subinst->rowscovered))
   {
      SCIPdebugMessage("threephase: all rows are already covered\n");
   }

   /* stop if all rows are covered by fixed variables */
   while(core->nactiveconstraints > SCIPhashtableGetNElements(subinst->rowscovered))
   {
      SCIP_CALL( markRowsCoveredByFixedVariables(scip, core, subinst, heurdata) );
      SCIP_CALL( subgradientOptimization(scip, core, subinst, mult_lb_subinst, heurdata->best_ub_subinst, heurdata) );
      /*SCIP_CALL( greedySetCover(scip, core, subinst, mult_lb_subinst) );*/
      SCIP_CALL( exploreNeighborhood(scip, mult_lb_subinst, heurdata) );

      /*exit(0);*/
      /*if(mult_lb_subinst->ub_greedy_local + subinst->costsfixed < heurdata->best_ub_subinst)
      {
         SCIP_CALL( copySetCoverSolution(core, subinst, heurdata->best_ub_subinst_sol, &heurdata->best_ub_subinst, mult_lb_subinst->x_greedy_local) );

         if(heurdata->best_ub_subinst < heurdata->best_ub_inst)
         {
            SCIP_CALL( copySetCoverSolution(core, subinst, heurdata->best_ub_inst_sol, &heurdata->best_ub_inst, heurdata->best_ub_subinst_sol) );

            if(heurdata->best_ub_inst < heurdata->best_ub)
            {
               SCIP_CALL( copySetCoverSolution(core, inst, heurdata->best_ub_sol, &heurdata->best_ub, heurdata->best_ub_inst_sol) );

               SCIP_CALL( removeRedundantColumns(scip, heurdata, heurdata->best_ub_sol, &heurdata->best_ub) );
               SCIPdebugMessage("new upper bound: %f\n", heurdata->best_ub);
            }
         }
      }*/

      if(heurdata->best_ub_subinst < heurdata->best_ub_inst)
      {
         SCIP_CALL( copySetCoverSolution(core, subinst, heurdata->best_ub_inst_sol, &heurdata->best_ub_inst, heurdata->best_ub_subinst_sol) );

         if(heurdata->best_ub_inst < heurdata->best_ub)
         {
            SCIP_CALL( copySetCoverSolution(core, inst, heurdata->best_ub_sol, &heurdata->best_ub, heurdata->best_ub_inst_sol) );

            SCIP_CALL( removeRedundantColumns(scip, heurdata, heurdata->best_ub_sol, &heurdata->best_ub) );
            SCIPdebugMessage("new upper bound: %f\n", heurdata->best_ub);
         }
      }


      if(subinst->costsfixed + mult_lb_subinst->lb_lagrange_local >= heurdata->best_ub)
         break;

      for(i = 0; i < core->nsolgreedy; i++)
      {
         fixVariable(subinst, core->variables[core->solgreedy[i]]);
         subinst->costsfixed += SCIPvarGetObj(core->variables[core->solgreedy[i]]);

         /* fix at least max(1, nconstraints / 200) variables */
         if(i > core->nconstraints / 200)
            break;
      }
   }

   SCIP_CALL( checkSetCover(scip, core, inst, heurdata->best_ub_inst_sol, &success, heurdata) );
   if(success == FALSE)
      SCIPdebugMessage("three-phase: final solution is not a valid set cover\n");

   return SCIP_OKAY;
}

static SCIP_RETCODE setCoveringHeuristic(SCIP *scip, SCIP_HEUR *heur)
{
   SCIP_HEURDATA *heurdata;
   SCIP_Bool stopcft = FALSE;
   SCIP_Bool success;
   SCIP_Bool computecorecolumns = TRUE;
   SCIP_Bool computecorerows = TRUE;
   int niter = 0;
   int nitercore = 0;
   int coret = 10;
   SCIP_Real corelb = 0.0;
   SCP_Core *core;
   SCP_Instance *inst;
   SCIP_Real pi = SCP_PI_MIN;
   int niternoimp = 0; /* number of iterations for how long no improvement was made */

   heurdata = SCIPheurGetData(heur);

   /* basic setup, for each row find first five columns covering it */
   SCIP_CALL( initTentativeCore(scip, &heurdata->core) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->vars, heurdata->core.maxconstraintvariables) );

   SCIPdebugMessage("%i rows, %i columns\n", heurdata->core.nvariables, heurdata->core.nconstraints);

   if(computecorecolumns == TRUE)
      SCIP_CALL( computeCoreColumns(scip, &heurdata->core, heurdata) );

   if(computecorerows == TRUE)
      SCIP_CALL( computeCoreRows(scip, &heurdata->core, heurdata) );

   SCIP_CALL( pqueue_init(scip, &heurdata->greedy_queue) );

   /* set up basic instance. so far no variables are fixed */
   SCIP_CALL( initInstance(scip, &heurdata->inst) );
   SCIP_CALL( initInstance(scip, &heurdata->subinst) );

   core = &heurdata->core;
   inst = &heurdata->inst;

   SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->greedy_colpos, core->nvariables) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->greedy_colmu, core->nvariables) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->greedy_colgamma, core->nvariables) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->greedy_colscore, core->nvariables) );
   SCIP_CALL( initInstance(scip, &heurdata->greedy_inst) );

   SCIP_CALL( allocateMemoryForSolution(scip, &heurdata->core, &heurdata->mult_best_lb_inst) );
   SCIP_CALL( allocateMemoryForSolution(scip, &heurdata->core, &heurdata->mult_best_lb_subinst) );
   SCIP_CALL( allocateMemoryForSolution(scip, &heurdata->core, &heurdata->mult_best_lb_total) );
   SCIP_CALL( allocateMemoryForSolution(scip, &heurdata->core, &heurdata->tp_mult_lb_subinst) );

   SCIP_CALL( SCIPhashtableCreate(&heurdata->best_ub_inst_sol, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );
   SCIP_CALL( SCIPhashtableCreate(&heurdata->best_ub_subinst_sol, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );
   SCIP_CALL( SCIPhashtableCreate(&heurdata->best_ub_sol, SCIPblkmem(scip), SCIPcalcHashtableSize(10), hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );

   heurdata->mult_best_lb_total.lb_lagrange_global = 0;
   heurdata->best_ub = SCIP_REAL_MAX;

   while(stopcft == FALSE)
   {
      SCIP_Bool redefine = FALSE;
      SCIP_Real currentub = heurdata->best_ub;

      /* 1. derive the reduced subinstance by marking rows covered by fixed variables */
      SCIPhashtableRemoveAll(inst->rowscovered);
      markRowsCoveredByFixedVariables(scip, core, inst, heurdata);

      /* call three-phase as long as 'inst' contains uncovered rows */
      if(core->nactiveconstraints > SCIPhashtableGetNElements(inst->rowscovered))
      {
         SCIPdebugMessage("%lli variables are fixed, %lli rows are covered\n", SCIPhashtableGetNElements(inst->varsfixed), SCIPhashtableGetNElements(inst->rowscovered));

         /* 2. apply procedure three-phase to find an optimal lagrange multiplier */
         SCIP_CALL( threePhase(scip, heur, heurdata) );
      }
      else
      {
         /* compute new core when all of its variables are fixed */
         redefine = TRUE;
      }

      /* stop if maximum number of iterations is reached */
      if(niter++ == SCP_MAX_ITER)
         stopcft = TRUE;

      /* set pi to PI_MIN if the current best solution was found in this iteration */
      if(heurdata->best_ub < currentub)
      {
         niternoimp = 0;
         pi = SCP_PI_MIN;
      }
      else
      {
         /* increase pi if no better solution was found, i.e. fix more variables in order to cover more rows */
         pi = pi * SCP_PI_ALPHA;
         niternoimp++;
      }

      /* stop if no improvement during the last 'SCP_MAX_ITER_NO_IMP' iterations */
      if(niternoimp == SCP_MAX_ITER_NO_IMP)
         stopcft = TRUE;

      /* stop if UB <= beta * LB */
      if(heurdata->mult_best_lb_total.lb_lagrange_global * SCP_BETA >= heurdata->best_ub)
         stopcft = TRUE;

      /* redefine core if the current core was worked on for 'coret' iterations */
      if(nitercore++ == coret)
         redefine = TRUE;

      if(stopcft)
         break;

      /* TODO: adapt conditions under which the core is redefined */
      if(redefine == TRUE)
      {
         /* stop if the last core did not lead to any improvements */
         if(SCIPisEQ(scip, corelb, heurdata->mult_best_lb_total.lb_lagrange_global) == TRUE)
            stopcft = TRUE;
         else
         {
            SCIP_CALL( redefineCore(scip, heurdata) );
            SCIPhashtableRemoveAll(inst->varsfixed);
            inst->costsfixed = 0.0;
            pi = SCP_PI_MIN;
            nitercore = 0;
            corelb = heurdata->mult_best_lb_total.lb_lagrange_global;
         }
      }
      else
      {
         SCIP_CALL( markRowsCoveredByFixedVariables(scip, core, inst, heurdata) );
         SCIP_CALL( computeDelta(scip, core, inst, heurdata->mult_best_lb_total.lagrangian_costs_global, heurdata->best_ub_sol, pi, heurdata) );
      }

      SCIPdebugMessage("iteration %i: best lower bound: %f, best upper bound: %f\n", niter, heurdata->mult_best_lb_total.lb_lagrange_global, heurdata->best_ub);
   }

   SCIPhashtableRemoveAll(inst->varsfixed);

   SCIP_CALL( checkSetCover(scip, core, inst, heurdata->best_ub_sol, &success, heurdata) );
   if(success == FALSE)
      SCIPdebugMessage("final solution is not a valid set cover\n");
   else
      SCIPdebugMessage("final solution has costs %f\n", heurdata->best_ub);

   SCIP_CALL( reportSolution(scip, core, inst, heurdata->best_ub_sol, heur) );

   SCIP_CALL( pqueue_destroy(scip, &heurdata->greedy_queue) );
   SCIPfreeBufferArray(scip, &heurdata->greedy_colpos);
   SCIPfreeBufferArray(scip, &heurdata->greedy_colmu);
   SCIPfreeBufferArray(scip, &heurdata->greedy_colgamma);
   SCIPfreeBufferArray(scip, &heurdata->greedy_colscore);

   SCIPfreeBufferArray(scip, &heurdata->vars);

   SCIP_CALL( freeMemoryForSolution(scip, &heurdata->mult_best_lb_inst) );
   SCIP_CALL( freeMemoryForSolution(scip, &heurdata->mult_best_lb_subinst) );
   SCIP_CALL( freeMemoryForSolution(scip, &heurdata->mult_best_lb_total) );
   SCIP_CALL( freeMemoryForSolution(scip, &heurdata->tp_mult_lb_subinst) );
   SCIP_CALL( freeInstance(scip, &heurdata->greedy_inst) );

   SCIPhashtableFree(&heurdata->best_ub_sol);
   SCIPhashtableFree(&heurdata->best_ub_inst_sol);
   SCIPhashtableFree(&heurdata->best_ub_subinst_sol);

   SCIP_CALL( freeInstance(scip, &heurdata->inst) );
   SCIP_CALL( freeInstance(scip, &heurdata->subinst) );
   SCIP_CALL( freeCore(scip, &heurdata->core) );

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
static
SCIP_DECL_HEURFREE(heurFreeSetcover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

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
   /*return SCIP_OKAY;*/
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
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

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
