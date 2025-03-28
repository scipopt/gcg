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

/**
 * @file   stat.c
 * @brief  Some printing methods for statistics
 * @author Alexander Gross
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "gcg/stat.h"
#include "gcg/scip_misc.h"
#include "gcg/pub_decomp.h"
#include "gcg/cons_decomp.h"
#include "gcg/struct_detector.h"
#include "gcg/pub_gcgvar.h"
#include "gcg/pricer_gcg.h"
#include "gcg/gcg.h"
#include "gcg/relax_gcg.h"


/** prints information about the best decomposition*/
SCIP_RETCODE GCGwriteDecompositionData(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_DECOMP* decomposition;
   SCIP* scip;
   GCG_DETECTOR* detector;
   GCG_DECTYPE type;
   const char* typeName;

   int i;
   int nblocks;
   int nlinkingconss;
   int nlinkingvars;
   int* nvarsinblocks;
   int* nconssinblocks;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);
   decomposition = GCGgetBestDecomp(gcg, TRUE);
   type = GCGdecompGetType(decomposition);
   typeName = GCGdecompGetStrType(type);

   detector = GCGdecompGetDetector(decomposition);

   nblocks = GCGdecompGetNBlocks(decomposition);

   nvarsinblocks = GCGdecompGetNSubscipvars(decomposition);
   nconssinblocks = GCGdecompGetNSubscipconss(decomposition);

   nlinkingvars = GCGdecompGetNLinkingvars(decomposition);
   nlinkingconss = GCGdecompGetNLinkingconss(decomposition);

   /* print information about decomposition type and number of blocks, vars, linking vars and cons */
   SCIPinfoMessage(scip, NULL, "Decomposition:\n");
   SCIPinfoMessage(scip, NULL, "Decomposition Type: %s \n", typeName);

   SCIPinfoMessage(scip, NULL, "Decomposition Detector: %s\n", detector == NULL ? "reader": detector->name);
   SCIPinfoMessage(scip, NULL, "Number of Blocks: %d \n", nblocks);
   SCIPinfoMessage(scip, NULL, "Number of LinkingVars: %d\n", nlinkingvars);
   SCIPinfoMessage(scip, NULL, "Number of LinkingCons: %d\n", nlinkingconss);

   /* print number of variables and constraints per block */
   SCIPinfoMessage(scip, NULL, "Block Information\n");
   SCIPinfoMessage(scip, NULL, "no.:\t\t#Vars\t\t#Constraints\n");
   for( i = 0; i < nblocks; i++ )
   {
      SCIPinfoMessage(scip, NULL, "%d:\t\t%d\t\t%d\n", i, nvarsinblocks[i], nconssinblocks[i]);
   }

   GCGdecompFree(gcg, &decomposition);

   return SCIP_OKAY;
}

/** prints additional solving statistics */
SCIP_RETCODE GCGwriteSolvingDetails(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_CLOCK* rootnodetime;
   SCIP* scip;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);
   rootnodetime = GCGgetRootNodeTime(gcg);
   SCIPinfoMessage(scip, NULL, "Solving Details    :\n");
   SCIPinfoMessage(scip, NULL, "  time in root node: %10.2f\n", SCIPgetClockTime(scip, rootnodetime));

   return SCIP_OKAY;
}

/** prints information about the creation of the Vars*/
SCIP_RETCODE GCGwriteVarCreationDetails(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* masterprob;
   SCIP_VAR** vars;
   SCIP_SOL* sol;
   SCIP_Real solvingtime;
   SCIP_Longint nnodes;

   int nvars, i, n;
   SCIP_Longint* createnodestat;
   int nodes[2];         /* Wurzel Knoten und nicht wurzelknoten  */
   SCIP_Longint createtimestat[10];
   int createiterstat[10];
   int m;

   assert(gcg != NULL);

   masterprob = GCGgetMasterprob(gcg);
   vars = SCIPgetVars(masterprob);
   nvars = SCIPgetNVars(masterprob);
   nnodes = SCIPgetNNodes(masterprob);
   sol = SCIPgetBestSol(masterprob);

   solvingtime = SCIPgetSolvingTime(masterprob);
   assert(nnodes < INT_MAX);
   SCIP_CALL( SCIPallocBufferArray(masterprob, &createnodestat, (int)nnodes) ); /* lld doesn't work here */

   SCIPinfoMessage(masterprob, NULL, "AddedVarDetails:\n");

   for( i = 0; i < 10; i++ )
   {
      createtimestat[i] = 0;
      createiterstat[i] = 0;
   }

   nodes[0] = 0;
   nodes[1] = 0;

   SCIPinfoMessage(masterprob, NULL, "VAR: name\tnode\ttime\titer\trootredcostcall\tredcost\tgap\tsolval\trootlpsolval\n");
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real redcost;
      SCIP_Real gap;
      SCIP_Longint  node;
      SCIP_Real time;
      SCIP_Longint iteration;
      SCIP_Longint rootredcostcall;
      SCIP_Real rootlpsolval;

      node = GCGgetCreationNode(vars[i]);
      time = GCGgetCreationTime(vars[i]);
      iteration = GCGgetIteration(vars[i]);
      redcost = GCGgetRedcost(vars[i]);
      gap = GCGgetVarGap(vars[i]);
      rootredcostcall = GCGgetRootRedcostCall(vars[i]);

      rootlpsolval = NAN;

#ifdef SCIP_STATISTIC
      rootlpsolval = SCIPgetSolVal(masterprob, GCGmasterGetRootLPSol(gcg), vars[i]);
#endif
      SCIPinfoMessage(masterprob, NULL, "VAR: <%s>\t%lld\t%f\t%lld\t%lld\t%f\t%f\t%f\t%f\n", SCIPvarGetName(vars[i]), node, time,
         iteration, rootredcostcall, redcost, gap, SCIPgetSolVal(masterprob, sol, vars[i]), rootlpsolval);

      if( SCIPisEQ(masterprob, SCIPgetSolVal(masterprob, sol, vars[i]), 0.0) )
      {
         continue;
      }
      else
      {
         SCIPdebugMessage("var <%s> has sol value %f (%lld, %f)\n", SCIPvarGetName(vars[i]),
            SCIPgetSolVal(masterprob, sol, vars[i]), node, time);
      }

      n = (int)(100 * time / solvingtime) % 10;
      m = (int)(100 * iteration / SCIPgetNLPIterations(masterprob)) % 10;
      createiterstat[n]++;
      createtimestat[m]++;

      if( node == 1 )
      {
         nodes[0]++;
      }
      else
      {
         nodes[1]++;
      }
   }

   SCIPinfoMessage(masterprob, NULL, "Root node:\tAdded Vars %d\n", nodes[0]);
   SCIPinfoMessage(masterprob, NULL, "Leftover nodes:\tAdded Vars %d\n", nodes[1]);

   for( i = 0; i < 10; i++ )
   {
      SCIPinfoMessage(masterprob, NULL, "Time %d-%d%%: Vars: %lld \n", 10 * i, 10 * (i + 1), createtimestat[i]);
   }

   for( i = 0; i < 10; i++ )
   {
      SCIPinfoMessage(masterprob, NULL, "Iter %d-%d%%: Vars: %d \n", 10 * i, 10 * (i + 1), createiterstat[i]);
   }

   SCIPfreeBufferArray(masterprob, &createnodestat);

   return SCIP_OKAY;
}
