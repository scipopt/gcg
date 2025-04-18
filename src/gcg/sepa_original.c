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

/**@file   sepa_original.c
 * @brief  original separator
 * @author Gerald Gamrath
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "scip/scip.h"
#include "gcg/sepa_original.h"
#include "gcg/gcg.h"
#include "gcg/relax_gcg.h"
#include "gcg/pricer_gcg.h"
#include "gcg/struct_gcg.h"


#define SEPA_NAME         "original"
#define SEPA_DESC         "separator for separating cuts in the original problem, called in the master"
#define SEPA_PRIORITY     1000

#define SEPA_FREQ         1
#define SEPA_MAXBOUNDDIST 1.0
#define SEPA_USESSUBSCIP  FALSE
#define SEPA_DELAY        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define STARTMAXCUTS 50       /**< maximal cuts used at the beginning */


/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   GCG*                  gcg;                /**< GCG data structure */
   SCIP_ROW**            mastercuts;         /**< cuts in the master problem */
   SCIP_ROW**            origcuts;           /**< cuts in the original problem */
   int                   ncuts;          /**< number of cuts in the original problem */
   int                   maxcuts;            /**< maximal number of allowed cuts */
   SCIP_Bool             enable;             /**< parameter returns if master separator is enabled */
   int                   separationsetting;  /**< parameter returns which parameter setting is used for separation */
   SCIP_HASHMAP*         origcutidxmap;      /**< origcuts mappend to index */
};


/*
 * Local methods
 */


/** allocates enough memory to hold more cuts */
static
SCIP_RETCODE ensureSizeCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data data structure */
   int                   size                /**< new size of cut arrays */
   )
{
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(sepadata->mastercuts != NULL);
   assert(sepadata->origcuts != NULL);
   assert(sepadata->ncuts <= sepadata->maxcuts);
   assert(sepadata->ncuts >= 0);

   if( sepadata->maxcuts < size )
   {
      int newmaxcuts = SCIPcalcMemGrowSize(scip, size);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(sepadata->mastercuts), sepadata->maxcuts, newmaxcuts) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(sepadata->origcuts), sepadata->maxcuts, newmaxcuts) );
      sepadata->maxcuts = newmaxcuts;
   }
   assert(sepadata->maxcuts >= size);

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

#define sepaCopyOriginal NULL
#define sepaInitOriginal NULL
#define sepaInitsolOriginal NULL
#define sepaExecsolOriginal NULL

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeOriginal)
{
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);

   SCIPfreeBlockMemoryArray(scip, &(sepadata->origcuts), sepadata->maxcuts);
   SCIPfreeBlockMemoryArray(scip, &(sepadata->mastercuts), sepadata->maxcuts);
   SCIPhashmapFree(&sepadata->origcutidxmap);

   sepadata->gcg->sepaorig = NULL;

   SCIPfreeBlockMemory(scip, &sepadata);

   return SCIP_OKAY;
}

/** deinitialization method of separator (called before transformed problem is freed) */
static
SCIP_DECL_SEPAEXIT(sepaExitOriginal)
{
   SCIP* origscip;
   SCIP_SEPADATA* sepadata;
   int i;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   origscip = GCGgetOrigprob(sepadata->gcg);
   assert(origscip != NULL);

   for( i = 0; i < sepadata->ncuts; i++ )
   {
      SCIP_CALL( SCIPreleaseRow(origscip, &(sepadata->origcuts[i])) );
   }

   sepadata->ncuts = 0;
   SCIPhashmapRemoveAll(sepadata->origcutidxmap);

   return SCIP_OKAY;
}

/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolOriginal)
{
   SCIP_SEPADATA* sepadata;
   int i;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   assert(GCGgetOrigprob(sepadata->gcg) != NULL);

   for( i = 0; i < sepadata->ncuts; i++ )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(sepadata->mastercuts[i])) );
   }

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpOriginal)
{
   SCIP* origscip;
   SCIP_Bool delayed;
   SCIP_Bool cutoff;
   SCIP_CUT** cuts;
   SCIP_ROW* mastercut;
   SCIP_VAR** rowvars;
   SCIP_SEPADATA* sepadata;
   SCIP_CUTPOOL* cutpool;
   GCG* gcg;

   SCIP_VAR** mastervars;
   SCIP_Real* mastervals;
   int nmastervars;

   char name[SCIP_MAXSTRLEN];

   int ncuts;
   int i;
   int j;
   SCIP_Bool feasible;

   SCIP_Bool isroot;

   assert(scip != NULL);
   assert(result != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   gcg = sepadata->gcg;

   origscip = GCGgetOrigprob(gcg);
   assert(origscip != NULL);

   SCIPdebugMessage("sepaExeclpOriginal\n");

   *result = SCIP_DIDNOTFIND;

   if( !sepadata->enable )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("master LP not solved to optimality, do no separation!\n");
      return SCIP_OKAY;
   }

   if( GCGgetNRelPricingprobs(gcg) < GCGgetNPricingprobs(gcg) )
   {
      SCIPdebugMessage("aggregated pricing problems, do no separation!\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* ensure to separate current sol */
   SCIP_CALL( GCGrelaxUpdateCurrentSol(gcg) );

   if( GCGrelaxIsOrigSolFeasible(gcg) )
   {
      SCIPdebugMessage("Current solution is feasible, no separation necessary!\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   isroot = SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip);

   /* set parameter setting for separation */
   SCIP_CALL( SCIPsetSeparating(origscip, (SCIP_PARAMSETTING) sepadata->separationsetting, TRUE) );

   SCIP_CALL( SCIPseparateSol(origscip, GCGrelaxGetCurrentOrigSol(gcg),
         isroot, TRUE, FALSE, &delayed, &cutoff) );

   if( delayed && !cutoff )
   {
      SCIPdebugMessage("call delayed separators\n");

      SCIP_CALL( SCIPseparateSol(origscip, GCGrelaxGetCurrentOrigSol(gcg),
            isroot, TRUE, TRUE, &delayed, &cutoff) );
   }

   cutpool = SCIPgetGlobalCutpool(origscip);
   SCIPdebugMessage("SCIPseparateSol() found %d cuts!\n", SCIPcutpoolGetNCuts(cutpool));

   /* if cut off is detected set result pointer and return SCIP_OKAY */
   if( cutoff )
   {
      *result = SCIP_CUTOFF;

      /* disable separating again */
      SCIP_CALL( SCIPsetSeparating(origscip, SCIP_PARAMSETTING_OFF, TRUE) );

      return SCIP_OKAY;
   }

   cuts = SCIPcutpoolGetCuts(cutpool);
   ncuts = SCIPcutpoolGetNCuts(cutpool);

   /* save cuts in the origcuts array in the separator data */
   SCIP_CALL( ensureSizeCuts(scip, sepadata, sepadata->ncuts + ncuts) );

   mastervars = SCIPgetVars(scip);
   nmastervars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &mastervals, nmastervars) );

   for( i = 0; i < ncuts; i++ )
   {
      SCIP_ROW* origcut;
      SCIP_COL** cols;
      int ncols;
      SCIP_Real* vals;
      SCIP_Real shift;

      origcut = SCIPcutGetRow(cuts[i]); /*lint !e679*/

      if( SCIPhashmapExists(sepadata->origcutidxmap, origcut) )
         continue;

      /* get columns and vals of the cut */
      ncols = SCIProwGetNNonz(origcut);
      cols = SCIProwGetCols(origcut);
      vals = SCIProwGetVals(origcut);

      /* get the variables corresponding to the columns in the cut */
      SCIP_CALL( SCIPallocBufferArray(scip, &rowvars, ncols) );
      for( j = 0; j < ncols; j++ )
      {
         rowvars[j] = SCIPcolGetVar(cols[j]);
      }

      /* transform the original variables to master variables */
      shift = GCGtransformOrigvalsToMastervals(gcg, rowvars, vals, ncols, mastervars, mastervals,
            nmastervars);

      /* create new cut in the master problem */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "mc_%s", SCIProwGetName(origcut));
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &mastercut, sepa, name,
            ( SCIPisInfinity(scip, -SCIProwGetLhs(origcut)) ?
              SCIProwGetLhs(origcut) : SCIProwGetLhs(origcut) - SCIProwGetConstant(origcut) - shift),
            ( SCIPisInfinity(scip, SCIProwGetRhs(origcut)) ?
              SCIProwGetRhs(origcut) : SCIProwGetRhs(origcut) - SCIProwGetConstant(origcut) - shift),
            SCIProwIsLocal(origcut), TRUE, FALSE) );

      /* add master variables to the cut */
      SCIP_CALL( SCIPaddVarsToRow(scip, mastercut, nmastervars, mastervars, mastervals) );

      /* add the cut to the master problem */
      SCIP_CALL( SCIPaddRow(scip, mastercut, FALSE, &feasible) );
      sepadata->mastercuts[sepadata->ncuts] = mastercut;
      sepadata->origcuts[sepadata->ncuts] = origcut;
      SCIP_CALL( SCIPhashmapInsertInt(sepadata->origcutidxmap, origcut, sepadata->ncuts) );
      SCIP_CALL( SCIPcaptureRow(origscip, origcut) );
      sepadata->ncuts++;

#ifdef SCIP_DEBUG
      SCIPdebugMessage("Cut %d (efficacious %d):\n", i, SCIPisCutEfficacious(scip, NULL, mastercut));
      SCIP_CALL( SCIPprintRow(scip, mastercut, NULL) );
      SCIPdebugMessage("\n\n");
#endif

      SCIPfreeBufferArray(scip, &rowvars);
   }

   if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   SCIPdebugMessage("%d cuts are in the original sepastore!\n", SCIPgetNCuts(origscip));
   SCIPdebugMessage("%d cuts are in the master sepastore!\n", SCIPgetNCuts(scip));

   /* disable separation for original problem again */
   SCIP_CALL( SCIPsetSeparating(origscip, SCIP_PARAMSETTING_OFF, TRUE) );

   SCIPfreeBufferArray(scip, &mastervals);

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the original separator and includes it in SCIP */
SCIP_RETCODE GCGincludeSepaOriginal(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP* masterprob = GCGgetMasterprob(gcg);

   /* create original separator data */
   SCIP_CALL( SCIPallocBlockMemory(masterprob, &sepadata) );
   sepadata->gcg = gcg;

   sepadata->maxcuts = SCIPcalcMemGrowSize(masterprob, STARTMAXCUTS);
   SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &(sepadata->origcuts), sepadata->maxcuts) ); /*lint !e506*/
   SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &(sepadata->mastercuts), sepadata->maxcuts) ); /*lint !e506*/
   SCIP_CALL( SCIPhashmapCreate(&sepadata->origcutidxmap, SCIPblkmem(masterprob), sepadata->maxcuts) );
   sepadata->ncuts = 0;

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(masterprob, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaCopyOriginal, sepaFreeOriginal, sepaInitOriginal, sepaExitOriginal,
         sepaInitsolOriginal, sepaExitsolOriginal,
         sepaExeclpOriginal, sepaExecsolOriginal,
         sepadata) );
   
   gcg->sepaorig = SCIPfindSepa(masterprob, SEPA_NAME);
   assert(gcg->sepaorig != NULL);

   SCIP_CALL( SCIPaddBoolParam(GCGgetOrigprob(gcg), "sepa/" SEPA_NAME "/enable", "enable original separator",
         &(sepadata->enable), FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(GCGgetOrigprob(gcg), "sepa/" SEPA_NAME "/paramsetting", "parameter returns which parameter setting is used for "
      "separation (default = 0, aggressive = 1, fast = 2", &(sepadata->separationsetting), FALSE, 0, 0, 2, NULL, NULL) );

   return SCIP_OKAY;
}


/** returns the array of original cuts in the original problem saved in the separator data */
SCIP_ROW** GCGsepaGetOriginalSepaOrigcuts(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;

   assert(gcg != NULL);

   sepa = GCGgetSepaorig(gcg);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   return sepadata->origcuts;
}


/** returns the number of cuts saved in the separator data */
int GCGsepaGetNOriginalSepaCuts(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;

   assert(gcg != NULL);

   sepa = GCGgetSepaorig(gcg);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   return sepadata->ncuts;
}

/** returns the array of original cuts in the master problem saved in the separator data */
SCIP_ROW** GCGsepaGetOriginalSepaMastercuts(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;

   assert(gcg != NULL);

   sepa = GCGgetSepaorig(gcg);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   return sepadata->mastercuts;
}

/** adds given original cut in both the original and master problem to master separator data */
SCIP_RETCODE GCGsepaAddOriginalSepaCuts(
   GCG*                 gcg,                /**< GCG data structure */
   SCIP_ROW*            origcut,            /**< pointer to orginal cut in the original problem */
   SCIP_ROW*            mastercut           /**< pointer to original cut in the master problem */
   )
{
   SCIP* scip;
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;

   assert(gcg != NULL);

   scip = GCGgetMasterprob(gcg);
   sepa = GCGgetSepaorig(gcg);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIP_CALL( ensureSizeCuts(scip, sepadata, sepadata->ncuts + 1) );

   sepadata->origcuts[sepadata->ncuts] = origcut;
   sepadata->mastercuts[sepadata->ncuts] = mastercut;
   SCIP_CALL( SCIPhashmapInsertInt(sepadata->origcutidxmap, origcut, sepadata->ncuts) );
   SCIP_CALL( SCIPcaptureRow(scip, sepadata->origcuts[sepadata->ncuts]) );
   SCIP_CALL( SCIPcaptureRow(scip, sepadata->mastercuts[sepadata->ncuts]) );

   ++(sepadata->ncuts);

   return SCIP_OKAY;
}

/** checks whether a given original cut in the original problem is already known */
SCIP_Bool GCGsepaOriginalSepaOrigcutExists(
   GCG*                 gcg,             /**< GCG data structure */
   SCIP_ROW*            origcut          /**< pointer to orginal cut in the original problem */
   )
{
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;

   assert(gcg != NULL);

   sepa = GCGgetSepaorig(gcg);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   return SCIPhashmapExists(sepadata->origcutidxmap, origcut);
}
