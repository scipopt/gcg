/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_roundbound.c
 * @brief  roundbound presolver: round fractional bounds on integer variables
 * @author Tobias Achterberg
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "gcg/presol_roundbound.h"


#define PRESOL_NAME            "roundbound"
#define PRESOL_DESC            "roundbound presolver: round fractional bounds on integers"
#define PRESOL_PRIORITY        +9000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_FAST /* timing of the presolver (fast, medium, or exhaustive) */

#ifdef FIXSIMPLEVALUE
#define MAXDNOM                 10000LL /**< maximal denominator for simple rational fixed values */
#endif


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyRoundbound)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( GCGincludePresolRoundbound(scip) );

   return SCIP_OKAY;
}


/** presolving execution method */
static
SCIP_DECL_PRESOLEXEC(presolExecRoundbound)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* get the problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* scan the variables for roundbound bound reductions
    * (loop backwards, since a variable fixing can change the current and the subsequent slots in the vars array)
    */
   for( v = nvars-1; v >= 0; --v )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      /* get variable's bounds */
      lb = SCIPvarGetLbGlobal(vars[v]);
      ub = SCIPvarGetUbGlobal(vars[v]);

      /* is variable integral? */
      if( SCIPvarGetType(vars[v]) != SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real newlb;
         SCIP_Real newub;

         /* round fractional bounds on integer variables */
         newlb = SCIPfeasCeil(scip, lb);
         newub = SCIPfeasFloor(scip, ub);

         /* check bounds on variable for infeasibility */
         if( newlb > newub + 0.5 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
               "problem infeasible: integral variable <%s> has bounds [%.17f,%.17f] rounded to [%.17f,%.17f]\n",
               SCIPvarGetName(vars[v]), lb, ub, newlb, newub);
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         /* round fractional bounds */
         if( !SCIPisFeasEQ(scip, lb, newlb) )
         {
            SCIPdebugMessage("rounding lower bound of integral variable <%s>: [%.17f,%.17f] -> [%.17f,%.17f]\n",
               SCIPvarGetName(vars[v]), lb, ub, newlb, ub);
            SCIP_CALL( SCIPchgVarLb(scip, vars[v], newlb) );
            (*nchgbds)++;
         }
         if( !SCIPisFeasEQ(scip, ub, newub) )
         {
            SCIPdebugMessage("rounding upper bound of integral variable <%s>: [%.17f,%.17f] -> [%.17f,%.17f]\n",
               SCIPvarGetName(vars[v]), newlb, ub, newlb, newub);
            SCIP_CALL( SCIPchgVarUb(scip, vars[v], newub) );
            (*nchgbds)++;
         }

      }
      else
      {
         /* check bounds on continuous variable for infeasibility */
         if( SCIPisFeasGT(scip, lb, ub) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
               "problem infeasible: continuous variable <%s> has bounds [%.17f,%.17f]\n",
               SCIPvarGetName(vars[v]), lb, ub);
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the roundbound presolver and includes it in SCIP */
SCIP_RETCODE GCGincludePresolRoundbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata = NULL;
   SCIP_PRESOL* presolptr;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presolptr, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecRoundbound, presoldata) );

   assert(presolptr != NULL);

   SCIP_CALL( SCIPsetPresolCopy(scip, presolptr, presolCopyRoundbound) );

   return SCIP_OKAY;
}
