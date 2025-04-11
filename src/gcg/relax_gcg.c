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

/**@file    relax_gcg.c
 *
 * @brief   GCG relaxator
 * @author  Gerald Gamrath
 * @author  Martin Bergner
 * @author  Alexander Gross
 * @author  Michael Bastubbe
 * @author Erik Muehmer
 *
 * \bug
 * - The memory limit is not strictly enforced
 * - Dealing with timelimits is a working hack only
 * - CTRL-C handling is very flaky
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG

#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#include "gcg/struct_locks.h"
#endif
#include "gcg/type_locks.h"

#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "scip/misc.h"
#include "scip/clock.h"

#include "gcg/relax_gcg.h"
#include "gcg/gcg.h"
#include "gcg/struct_gcg.h"
#include "gcg/struct_branchgcg.h"

#include "gcg/cons_origbranch.h"
#include "gcg/cons_masterbranch.h"
#include "gcg/pricer_gcg.h"
#include "gcg/benders_gcg.h"
#include "gcg/masterplugins.h"
#include "gcg/bendersplugins.h"
#include "gcg/cons_decomp.h"
#include "gcg/scip_misc.h"
#include "gcg/solver_knapsack.h"

#include "gcg/params_visu.h"

#ifndef NO_AUT_LIB
#include "symmetry/pub_automorphism.h"
#include "symmetry/automorphism.h"
#endif

#define RELAX_NAME             "gcg"
#define RELAX_DESC             "relaxator for gcg project representing the master lp"
#define RELAX_PRIORITY         -1
#define RELAX_FREQ             1
#define RELAX_INCLUDESLP       TRUE

#define DEFAULT_DISCRETIZATION TRUE
#define DEFAULT_MIPDISCRETIZATION TRUE
#define DEFAULT_AGGREGATION TRUE
#define DEFAULT_DISPINFOS FALSE
#define DEFAULT_MODE GCG_DECMODE_DANTZIGWOLFE  /**< the decomposition mode that GCG will use. (0: Dantzig-Wolfe (default),
                                                    1: Benders' decomposition, 2: solve original problem) */
#define DEFAULT_BLISS TRUE
#define DEFAULT_BLISS_SEARCH_NODE_LIMIT 0
#define DEFAULT_BLISS_GENERATOR_LIMIT 100
#define DEFAULT_AGGREGATIONNCONSSLIMIT 300
#define DEFAULT_AGGREGATIONNVARSLIMIT 300

#define DELVARS

/*
 * Data structures
 */

/** relaxator data */

struct SCIP_RelaxData
{
   /* problems and convexity constraints */
   GCG*                  gcg;                /**< GCG data structure */
   SCIP**                pricingprobs;       /**< the array of pricing problems */
   int                   npricingprobs;      /**< the number of pricing problems */
   int                   maxpricingprobs;    /**< capacity of pricingprobs */
   int                   nrelpricingprobs;   /**< the number of relevant pricing problems */
   int*                  blockrepresentative;/**< number of the pricing problem, that represents the i-th problem */
   int*                  nblocksidentical;   /**< number of pricing blocks represented by the i-th pricing problem */
   SCIP_CONS**           convconss;          /**< array of convexity constraints, one for each block */
   int                   ntransvars;         /**< number of variables directly transferred to the master problem */
   int                   nlinkingvars;       /**< number of linking variables */
   int                   nvarlinkconss;      /**< number of constraints that ensure that copies of linking variables have the same value */
   SCIP_Real             pricingprobsmemused; /**< sum of memory used after problem creation stage of all pricing problems */

   /* constraint data */
   SCIP_CONS**           masterconss;        /**< array of constraints in the master problem */
   SCIP_CONS**           origmasterconss;    /**< array of constraints in the original problem that belong to the
                                              * master problem */
   SCIP_CONS**           linearmasterconss;  /**< array of linear constraints equivalent to the cons in
                                              * the original problem that belong to the master problem */
   SCIP_CONS**           varlinkconss;       /**< array of constraints ensuring linking vars equality */
   int*                  varlinkconsblock;   /**< array of constraints ensuring linking vars equality */
   int                   maxmasterconss;     /**< length of the array mastercons */
   int                   nmasterconss;       /**< number of constraints saved in mastercons */

   SCIP_SOL*             currentorigsol;     /**< current lp solution transformed into the original space */
   SCIP_Bool             origsolfeasible;    /**< is the current lp solution primal feasible in the original space? */
   SCIP_Longint          lastmasterlpiters;  /**< number of lp iterations when currentorigsol was updated the last time */
   SCIP_Longint          lastmasternode;     /**< number of current node when currentorigsol was updated the last time */
   SCIP_SOL*             lastmastersol;      /**< last feasible master solution that was added to the original problem */
   SCIP_CONS**           markedmasterconss;  /**< array of conss that are marked to be in the master */
   int                   nmarkedmasterconss; /**< number of elements in array of conss that are marked to be in the master */
   int                   maxmarkedmasterconss; /**< capacity of markedmasterconss */
   SCIP_Longint          lastsolvednodenr;   /**< node number of the node that was solved at the last call of the relaxator */

   /* branchrule data */
   GCG_BRANCHRULE**      branchrules;        /**< branching rules registered in the relaxator */
   int                   nbranchrules;       /**< number of branching rules registered in the relaxator */
   GCG_BRANCHRULE**      activebranchrules;  /**< branching rules that created extended master conss (cache) */
   GCG_BRANCHDATA**      activebranchdata;   /**< data represeting the branching decisions of the active nodes (cache) */
   GCG_EXTENDEDMASTERCONSDATA**   activebranchextendedmasterconss; /**< array of extended master conss that are active in the current node (cache) */
   int                   nactivebranchextendedmasterconss; /**< number of extended master conss that are active in the current node */
   int                   maxactivebranchextendedmasterconss; /**< capacity of activebranchextendedmasterconss */

   /* parameter data */
   SCIP_Bool             discretization;     /**< TRUE: use discretization approach; FALSE: use convexification approach */
   SCIP_Bool             mipdiscretization;  /**< TRUE: use discretization approach in MIPs; FALSE: use convexification approach in MIPs*/
   SCIP_Bool             aggregation;        /**< should identical blocks be aggregated (only for discretization approach)? */
   SCIP_Bool             masterissetpart;    /**< is the master a set partitioning problem? */
   SCIP_Bool             masterissetcover;   /**< is the master a set covering problem? */
   SCIP_Bool             dispinfos;          /**< should additional information be displayed? */
   GCG_DECMODE           mode;               /**< the decomposition mode for GCG. 0: Dantzig-Wolfe (default), 1: Benders' decomposition, 2: automatic */
   int                   origverblevel;      /**< the verbosity level of the original problem */
   SCIP_Bool             usesymmetrylib;     /**< should symmetry detection lib be used to check for identical blocks? */
   int                   searchnodelimit;    /**< bliss search node limit (requires patched bliss version) */
   int                   generatorlimit;     /**< bliss generator limit (requires patched bliss version) */
   int                   aggregationnconsslimit;          /**< if this limit on the number of constraints of a block is exceeded the aggregation information for this block is not calculated */
   int                   aggregationnvarslimit;           /**< if this limit on the number of variables of a block is exceeded the aggregation information for this block is not calculated */

   /* data for probing */
   SCIP_Bool             masterinprobing;    /**< is the master problem in probing mode? */
   SCIP_HEUR*            probingheur;        /**< heuristic that started probing in master problem, or NULL */
   SCIP_SOL*             storedorigsol;      /**< original solution that was stored before the probing */
   SCIP_Bool             storedfeasibility;  /**< is the stored original solution feasible? */

   /* structure information */
   GCG_DECOMP*           decomp;             /**< structure information */
   SCIP_Bool             relaxisinitialized; /**< indicates whether the relaxator is initialized */

   /* statistical information */
   SCIP_Longint          simplexiters;       /**< cumulative simplex iterations */
   SCIP_CLOCK*           rootnodetime;       /**< time in root node */

   /* visualization parameter */
   GCG_PARAMDATA*        paramsvisu;         /**< parameters for visualization */

   /* stashed limit settings */
   SCIP_Bool             limitsettingsstashed;   /**< are limit settings currently stashed? */
   SCIP_Longint          stashednodelimit;       /**< stashed node limit */
   SCIP_Longint          stashedstallnodelimit;  /**< stashed stalling node limit */
   SCIP_Real             stashedgaplimit;        /**< stashed gap limit */
   int                   stashedsollimit;        /**< stashed solution limit */
   SCIP_Real             stashedtimelimit;       /**< stashed time limit */

#ifdef _OPENMP
   GCG_LOCKS*            locks;                  /** OpenMP locks */
#endif
};


/*
 * Local methods
 */

/** add the activated branch node's extended master cons to the cache */
static
SCIP_RETCODE addActiveBranchExtendedmastercons(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< relaxator data data structure */
   GCG_BRANCHRULE*       branchrule,         /**< branchrule that was activated */
   GCG_BRANCHDATA*       branchdata          /**< branchdata that was activated */
   )
{
   GCG_EXTENDEDMASTERCONSDATA* extendedmasterconsdata;
   SCIP* origprob = GCGgetOrigprob(gcg);

   assert(origprob != NULL);
   assert(relaxdata != NULL);
   assert(branchrule != NULL);
   assert(branchdata != NULL);

   /* add only if branch creates an extended master cons */
   extendedmasterconsdata = NULL;
   if( branchrule->branchgetextendedmastercons == NULL )
      return SCIP_OKAY;
   SCIP_CALL( branchrule->branchgetextendedmastercons(gcg, branchdata, &extendedmasterconsdata) );
   if( extendedmasterconsdata == NULL )
      return SCIP_OKAY;

   if( relaxdata->maxactivebranchextendedmasterconss == 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(origprob, &(relaxdata->activebranchrules), 1) );
      SCIP_CALL( SCIPallocBlockMemoryArray(origprob, &(relaxdata->activebranchdata), 1) );
      SCIP_CALL( SCIPallocBlockMemoryArray(origprob, &(relaxdata->activebranchextendedmasterconss), 1) );
      relaxdata->maxactivebranchextendedmasterconss = 1;
   }
   else
   if( relaxdata->nactivebranchextendedmasterconss == relaxdata->maxactivebranchextendedmasterconss )
   {
      int newsize;
      newsize = SCIPcalcMemGrowSize(origprob, relaxdata->nactivebranchextendedmasterconss+1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(origprob, &(relaxdata->activebranchrules), relaxdata->maxactivebranchextendedmasterconss, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(origprob, &(relaxdata->activebranchdata), relaxdata->maxactivebranchextendedmasterconss, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(origprob, &(relaxdata->activebranchextendedmasterconss), relaxdata->maxactivebranchextendedmasterconss, newsize) );
      relaxdata->maxactivebranchextendedmasterconss = newsize;
   }
   assert(relaxdata->nactivebranchextendedmasterconss < relaxdata->maxactivebranchextendedmasterconss);
   assert(relaxdata->activebranchrules != NULL);
   assert(relaxdata->activebranchdata != NULL);
   assert(relaxdata->activebranchextendedmasterconss != NULL);

   relaxdata->activebranchrules[relaxdata->nactivebranchextendedmasterconss] = branchrule;
   relaxdata->activebranchdata[relaxdata->nactivebranchextendedmasterconss] = branchdata;
   relaxdata->activebranchextendedmasterconss[relaxdata->nactivebranchextendedmasterconss] = extendedmasterconsdata;
   relaxdata->nactivebranchextendedmasterconss++;

   return SCIP_OKAY;
}

/** drop the most recently added branch extended master cons data */
static
SCIP_RETCODE dropActiveBranchExtendedmastercons(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< relaxator data data structure */
   GCG_BRANCHRULE*       branchrule,         /**< branchrule that was deactivated */
   GCG_BRANCHDATA*       branchdata          /**< branchdata that was deactivated */
   )
{
   assert(gcg != NULL);
   assert(relaxdata != NULL);

   if( relaxdata->nactivebranchextendedmasterconss == 0 )
      return SCIP_OKAY;

   /* drop only if branch created an extended master cons */
   if( relaxdata->activebranchrules[relaxdata->nactivebranchextendedmasterconss-1] != branchrule
      || relaxdata->activebranchdata[relaxdata->nactivebranchextendedmasterconss-1] != branchdata )
      return SCIP_OKAY;

   assert(relaxdata->activebranchextendedmasterconss != NULL);
   relaxdata->nactivebranchextendedmasterconss--;
   relaxdata->activebranchextendedmasterconss[relaxdata->nactivebranchextendedmasterconss] = NULL;

   return SCIP_OKAY;
}

/** sets the number of the block, the given original variable belongs to */
static
SCIP_RETCODE setOriginalVarBlockNr(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< relaxator data data structure */
   SCIP_VAR*             var,                /**< variable to set the block number for */
   int                   newblock            /**< number of the block, the variable belongs to */
   )
{
   int blocknr;

   assert(gcg != NULL);
   assert(var != NULL);
   assert(newblock >= 0 || (GCGgetDecompositionMode(relaxdata->gcg) == GCG_DECMODE_BENDERS && newblock == -2));

   assert(SCIPvarIsOriginal(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(relaxdata != NULL);

   blocknr = GCGvarGetBlock(var);
   assert(GCGvarIsOriginal(var));

   assert(relaxdata->npricingprobs > 0);
   assert(newblock < relaxdata->npricingprobs);
   assert(blocknr >= -2 && blocknr < relaxdata->npricingprobs);

   /* var belongs to no block so far, just set the new block number */
   if( blocknr == -1 )
   {
      assert(newblock >= 0);
      GCGvarSetBlock(var, newblock);
   }
   /* if var already belongs to another block, it is a linking variable */
   else if( blocknr != newblock )
   {
      SCIP_CALL( GCGoriginalVarAddBlock(gcg, var, newblock, relaxdata->npricingprobs, relaxdata->mode) );
      assert(newblock == -2 || GCGisLinkingVarInBlock(var, newblock));
      assert(GCGoriginalVarIsLinking(var));
   }
   blocknr = GCGvarGetBlock(var);
   assert(blocknr == -2 || blocknr == newblock);

   return SCIP_OKAY;
}

/** marks the constraint to be transferred to the master problem */
static
SCIP_RETCODE markConsMaster(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< relaxator data data structure */
   SCIP_CONS*            cons                /**< constraint that is forced to be in the master */
   )
{
   SCIP* scip;
#ifndef NDEBUG
   int i;
#endif
   assert(gcg != NULL);
   assert(cons != NULL);
   assert(relaxdata != NULL);

   scip = GCGgetOrigprob(gcg);

   /* allocate array, if not yet done */
   if( relaxdata->markedmasterconss == NULL )
   {
      relaxdata->maxmarkedmasterconss = SCIPcalcMemGrowSize(scip, SCIPgetNConss(scip));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(relaxdata->markedmasterconss), relaxdata->maxmarkedmasterconss) );
      relaxdata->nmarkedmasterconss = 0;
   }
   assert(relaxdata->nmarkedmasterconss <= SCIPgetNConss(scip));

#ifndef NDEBUG
   /* check that constraints are not marked more than one time */
   for( i = 0; i < relaxdata->nmarkedmasterconss; i++ )
      assert(relaxdata->markedmasterconss[i] != cons);
#endif

   /* save constraint */
   relaxdata->markedmasterconss[relaxdata->nmarkedmasterconss] = cons;
   relaxdata->nmarkedmasterconss++;

   return SCIP_OKAY;
}


/** converts the structure to the GCG format by setting the appropriate blocks and master constraints */
static
SCIP_RETCODE convertStructToGCG(
   GCG*                  gcg,                /**< GCG data structure          */
   SCIP_RELAXDATA*       relaxdata,          /**< relaxator data structure     */
   GCG_DECOMP*           decomp              /**< decomp data structure        */
   )
{
   SCIP* origprob;
   int i;
   int j;
   int k;
   int v;
   int nblocks;
   int nvars;
   SCIP_VAR** origvars;
   SCIP_HASHMAP* transvar2origvar;
   SCIP_CONS** linkingconss;
   int nlinkingconss;
   SCIP_VAR** linkingvars;
   int nlinkingvars;
   SCIP_VAR*** subscipvars;
   int* nsubscipvars;
   SCIP_CONS*** subscipconss;
   int* nsubscipconss;

   assert(decomp != NULL);
   assert(relaxdata != NULL);
   assert(gcg != NULL);

   assert(GCGdecompGetLinkingconss(decomp) != NULL || GCGdecompGetNLinkingconss(decomp) == 0);
   assert(GCGdecompGetNSubscipvars(decomp) != NULL || GCGdecompGetSubscipvars(decomp) == NULL);


   SCIP_CALL( GCGdecompAddRemainingConss(gcg, decomp) );
   SCIP_CALL( GCGdecompCheckConsistency(gcg, decomp) );

   origprob = GCGgetOrigprob(gcg);

   origvars = SCIPgetVars(origprob);
   nvars = SCIPgetNVars(origprob);
   linkingconss = GCGdecompGetLinkingconss(decomp);
   nlinkingconss = GCGdecompGetNLinkingconss(decomp);
   linkingvars = GCGdecompGetLinkingvars(decomp);
   nlinkingvars = GCGdecompGetNLinkingvars(decomp);
   subscipvars = GCGdecompGetSubscipvars(decomp);
   nsubscipvars = GCGdecompGetNSubscipvars(decomp);

   subscipconss = GCGdecompGetSubscipconss(decomp);
   nsubscipconss = GCGdecompGetNSubscipconss(decomp);
   nblocks = GCGdecompGetNBlocks(decomp);

   SCIP_CALL( SCIPhashmapCreate(&transvar2origvar, SCIPblkmem(origprob), nvars) );
   relaxdata->npricingprobs = nblocks;
   SCIP_CALL( GCGcreateOrigVarsData(gcg) );

   SCIPdebugMessage("Copying structure with %d blocks, %d linking vars and %d linking constraints.\n", nblocks, nlinkingvars, nlinkingconss);

   /* set master constraints */
   for( i = 0; i < nlinkingconss; ++i )
   {
      assert(linkingconss[i] != NULL);
      /* SCIPdebugMessage("\tProcessing linking constraint %s.\n", SCIPconsGetName(linkingconss[i])); */
      if( SCIPconsIsActive(linkingconss[i]) )
      {
         SCIP_CALL( markConsMaster(gcg, relaxdata, linkingconss[i]) );
      }
   }

   /* prepare the map from transformed to original variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* transvar;

      SCIP_CALL( SCIPgetTransformedVar(origprob, origvars[i], &transvar) );
      assert(transvar != NULL);

      SCIP_CALL( SCIPhashmapInsert(transvar2origvar, transvar, origvars[i]) );
   }

   for( i = 0; i < nblocks; ++i )
   {
      /* SCIPdebugMessage("\tProcessing block %d (%d conss, %d vars).\n", i, nsubscipconss[i], nsubscipvars[i]); */
      assert((subscipvars[i] == NULL) == (nsubscipvars[i] == 0));
      for( j = 0; j < nsubscipvars[i]; ++j )
      {
         SCIP_VAR* relevantvar;
         assert(subscipvars[i][j] != NULL);
         relevantvar = SCIPvarGetProbvar(subscipvars[i][j]);

         /* If there is a corresponding original (untransformed) variable, assign it to the block */
         if( SCIPhashmapGetImage(transvar2origvar, subscipvars[i][j]) != NULL )
         {
            SCIP_VAR* origvar;

            origvar = (SCIP_VAR*) SCIPhashmapGetImage(transvar2origvar, subscipvars[i][j]);
            assert(SCIPvarGetData(origvar) != NULL);

            SCIP_CALL( setOriginalVarBlockNr(gcg, relaxdata, origvar, i) );
            SCIPdebugMessage("\t\tOriginal var %s (%p) in block %d\n", SCIPvarGetName(subscipvars[i][j]), (void*) subscipvars[i][j], i);
         }

         /* Assign the corresponding problem variable to the block */
         if( SCIPvarGetData(relevantvar) == NULL )
            SCIP_CALL( GCGorigVarCreateData(gcg, relevantvar) );
         SCIP_CALL( setOriginalVarBlockNr(gcg, relaxdata, relevantvar, i) );

         SCIPdebugMessage("\t\tTransformed var %s (%p) in block %d\n", SCIPvarGetName(relevantvar), (void*) relevantvar, i);

         assert(SCIPvarGetData(subscipvars[i][j]) != NULL || SCIPvarGetData(relevantvar) != NULL);
      }
   }
   SCIPdebugMessage("\tProcessing linking variables.\n");
   for( i = 0; i < nlinkingvars; ++i )
   {
      int nfound = 0;

      if( GCGoriginalVarIsLinking(linkingvars[i]) )
         continue;

      SCIPdebugMessage("\tDetecting constraint blocks of linking var %s\n", SCIPvarGetName(linkingvars[i]));
      /* HACK; @todo find out constraint blocks more intelligently */
      for( j = 0; j < nblocks; ++j )
      {
         int found = FALSE;
         for( k = 0; k < nsubscipconss[j]; ++k )
         {
            SCIP_VAR** curvars;
            int        ncurvars;
            if( SCIPconsIsDeleted(subscipconss[j][k]) )
               continue;
            ncurvars = GCGconsGetNVars(origprob, subscipconss[j][k]);
            curvars = NULL;
            if( ncurvars > 0 )
            {
               SCIP_CALL( SCIPallocBufferArray(origprob, &curvars, ncurvars) );
               SCIP_CALL( GCGconsGetVars(origprob, subscipconss[j][k], curvars, ncurvars) );

               for( v = 0; v < ncurvars; ++v )
               {
                  if( SCIPvarGetProbvar(curvars[v]) == linkingvars[i] || curvars[v] == linkingvars[i] )
                  {
                     SCIPdebugMessage("\t\t%s is in %d\n", SCIPvarGetName(SCIPvarGetProbvar(curvars[v])), j);
                     assert(SCIPvarGetData(linkingvars[i]) != NULL);
                     SCIP_CALL( setOriginalVarBlockNr(gcg, relaxdata, SCIPvarGetProbvar(linkingvars[i]), j) );
                     found = TRUE;
                     break;
                  }
               }

               SCIPfreeBufferArray(origprob, &curvars);
            }

            if( found )
            {
               nfound++;
               break;
            }

         }
      }

      /* if the linking variable is only in one block, then it would not have been flagged as a linking variable. In
       * the Benders' decomposition case, then linking variable needs to be flagged as linking so that it is added to
       * the master problem.
       */
      if( nfound == 1 && GCGgetDecompositionMode(gcg) == GCG_DECMODE_BENDERS )
      {
         SCIP_CALL( setOriginalVarBlockNr(gcg, relaxdata, SCIPvarGetProbvar(linkingvars[i]), -2) );
      }
   }

   SCIPhashmapFree(&transvar2origvar);
   return SCIP_OKAY;
}

/** ensures size of masterconss array */
static
SCIP_RETCODE ensureSizeMasterConss(
   GCG*                  gcg,
   SCIP_RELAXDATA*       relaxdata,
   int                   size
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   assert(relaxdata != NULL);
   assert(relaxdata->masterconss != NULL);

   scip = GCGgetOrigprob(gcg);

   if( relaxdata->maxmasterconss < size )
   {
      int newsize = SCIPcalcMemGrowSize(scip, size);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(relaxdata->masterconss), relaxdata->maxmasterconss, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(relaxdata->origmasterconss), relaxdata->maxmasterconss, newsize) );
      relaxdata->maxmasterconss = newsize;

   }
   assert(relaxdata->maxmasterconss >= size);

   return SCIP_OKAY;
}

/** ensures size of branchrules array: enlarges the array by 1 */
static
SCIP_RETCODE ensureSizeBranchrules(
   GCG*                  gcg,
   SCIP_RELAXDATA*       relaxdata
   )
{
   assert(gcg != NULL);
   assert(relaxdata != NULL);
   assert((relaxdata->branchrules == NULL) == (relaxdata->nbranchrules == 0));

   if( relaxdata->nbranchrules == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(GCGgetOrigprob(gcg), &(relaxdata->branchrules), 1) ); /*lint !e506*/
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(GCGgetOrigprob(gcg), &(relaxdata->branchrules), (size_t)relaxdata->nbranchrules+1) );
   }

   return SCIP_OKAY;
}


/** check whether the master problem has a set partitioning or set covering structure */
static
SCIP_RETCODE checkSetppcStructure(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata           /**< relaxator data structure */
   )
{
   SCIP* scip;
   SCIP_CONS** masterconss;
   int nmasterconss;
   int i;

   assert(relaxdata->decomp != NULL);

   scip = GCGgetOrigprob(gcg);
   masterconss = GCGdecompGetLinkingconss(relaxdata->decomp);
   nmasterconss = GCGdecompGetNLinkingconss(relaxdata->decomp);
   assert(nmasterconss >= 0);
   assert(masterconss != NULL || nmasterconss == 0);

   if( nmasterconss == 0 || relaxdata->nvarlinkconss > 0 )
   {
      relaxdata->masterissetcover = FALSE;
      relaxdata->masterissetpart = FALSE;
      return SCIP_OKAY;
   }

   relaxdata->masterissetcover = TRUE;
   relaxdata->masterissetpart = TRUE;

   for( i = 0; i < nmasterconss; ++i )
   {
      assert(masterconss != NULL);

      if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(masterconss[i])), "setppc") == 0 )
      {
         switch( SCIPgetTypeSetppc(scip, masterconss[i]) )
         {
         case SCIP_SETPPCTYPE_COVERING:
            relaxdata->masterissetpart = FALSE;
            break;
         case SCIP_SETPPCTYPE_PARTITIONING:
            relaxdata->masterissetcover = FALSE;
            break;
         case SCIP_SETPPCTYPE_PACKING:
            relaxdata->masterissetcover = FALSE;
            relaxdata->masterissetpart = FALSE;
            break;
         }
      }
      else if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(masterconss[i])), "logicor") == 0 )
      {
         relaxdata->masterissetpart = FALSE;
         break;
      }
      else if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(masterconss[i])), "linear") == 0 )
      {
         SCIP_SETPPCTYPE type;

         if( GCGgetConsIsSetppc(scip, masterconss[i], &type) )
         {
            switch( type )
            {
            case SCIP_SETPPCTYPE_COVERING:
               relaxdata->masterissetpart = FALSE;
               break;
            case SCIP_SETPPCTYPE_PARTITIONING:
               relaxdata->masterissetcover = FALSE;
               break;
            case SCIP_SETPPCTYPE_PACKING:
               relaxdata->masterissetcover = FALSE;
               relaxdata->masterissetpart = FALSE;
               break;
            }
         }
         else
         {
            relaxdata->masterissetcover = FALSE;
            relaxdata->masterissetpart = FALSE;
            break;
         }
      }
      else
      {
         relaxdata->masterissetcover = FALSE;
         relaxdata->masterissetpart = FALSE;
         break;
      }
   }

   if( relaxdata->masterissetcover )
   {
      assert(!relaxdata->masterissetpart);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Master problem is a set covering problem.\n");
   }
   if( relaxdata->masterissetpart )
   {
      assert(!relaxdata->masterissetcover);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Master problem is a set partitioning problem.\n");
   }

   return SCIP_OKAY;
}

/** checks whether there are identical pricing blocks */
static
SCIP_RETCODE checkIdenticalBlocks(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< the relaxator data data structure*/
   SCIP_HASHMAP**        hashorig2pricingvar /**< mapping from orig to pricingvar for each block */
   )
{
   PARTIALDECOMP_C* partialdec;
   SCIP* scip;

   int i;
   int j;
   int k;

   int nrelevant;

   SCIPdebugMessage("checking identical blocks \n");
   assert(gcg != NULL);
   assert(relaxdata != NULL);

   scip = GCGgetOrigprob(gcg);

   for( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      relaxdata->blockrepresentative[i] = i;
      relaxdata->nblocksidentical[i] = 1;
   }

   relaxdata->nrelpricingprobs = relaxdata->npricingprobs;
   nrelevant = 0;

   if(  ( !relaxdata->discretization || !relaxdata->aggregation ) )
   {
      SCIPdebugMessage("discretization is off, aggregation is off\n");
      return SCIP_OKAY;
   }

   assert( SCIPgetNConss(scip) == GCGconshdlrDecompGetNFormerDetectionConssForID(gcg, GCGdecompGetPartialdecID(relaxdata->decomp)) );
   SCIPdebugMessage( "nconss: %d; ndetectionconss: %d -> using partialdec information for identity test \n", SCIPgetNConss(scip), GCGconshdlrDecompGetNFormerDetectionConssForID(gcg, GCGdecompGetPartialdecID(relaxdata->decomp) ) );

   GCGconshdlrDecompGetPartialdecFromID(gcg, GCGdecompGetPartialdecID(relaxdata->decomp), &partialdec);

   if( !GCGconshdlrDecompPartialdecAggregationInformationCalculated(partialdec) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Calculating aggregation information.\n");
      GCGconshdlrDecompPartialdecCalcAggregationInformation(partialdec, TRUE);
   }

   nrelevant = GCGconshdlrDecompPartialdecGetNEquivalenceClasses(partialdec);
   assert(nrelevant > 0 || GCGconshdlrDecompPartialdecGetNBlocks(partialdec) == 0);
   for( j = 0; j < nrelevant; ++j )
   {
      int rb = GCGconshdlrDecompPartialdecGetReprBlockForEqClass(partialdec, j);
      const int* eqclassblocks = GCGconshdlrDecompPartialdecGetBlocksForEqClass(partialdec, j);
      int neqclassblocks = GCGconshdlrDecompPartialdecGetNBlocksForEqClass(partialdec, j);

      SCIPdebugMessage("Block %d is relevant!\n", rb);
      relaxdata->nblocksidentical[rb] = neqclassblocks;

      assert(eqclassblocks[0] == rb);
      for( i = 1; i < neqclassblocks; ++i )
      {
         int b = eqclassblocks[i];
         const int* repvarmap = GCGconshdlrDecompPartialdecGetRepVarMap(partialdec, j, i);

         // block b will be represented by block rb
         relaxdata->blockrepresentative[b] = rb;
         relaxdata->nblocksidentical[b] = 0;
         SCIPdebugMessage("Block %d is represented by block %d.\n", b, rb);

         for( k = 0; k < GCGconshdlrDecompPartialdecGetNVarsForBlock(partialdec, b); ++k )
         {
            int rvi = repvarmap[k];
            SCIP_VAR* origvar = GCGconshdlrDecompPartialdecGetOrigVarForBlock(partialdec, b, k);
            SCIP_VAR* repvar = GCGconshdlrDecompPartialdecGetOrigVarForBlock(partialdec, rb, rvi);
            SCIP_VAR* pricingvar = GCGoriginalVarGetPricingVar(repvar);

            assert(GCGvarIsPricing(pricingvar));
            assert(GCGvarIsOriginal(origvar));
            assert(GCGoriginalVarGetPricingVar(origvar) != NULL);
            GCGoriginalVarSetPricingVar(origvar, pricingvar);
            assert(GCGvarGetBlock(pricingvar) == rb);
            assert(b == rb || !GCGoriginalVarIsLinking(origvar));
            SCIP_CALL( GCGpricingVarAddOrigVar(relaxdata->pricingprobs[rb], pricingvar, origvar) );
            SCIPdebugMessage("Var <%s> is mapped to <%s> (<%s>).\n", SCIPvarGetName(origvar), SCIPvarGetName(repvar),
                  SCIPvarGetName(pricingvar));
         }
      }
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Matrix has %d blocks, using %d%s pricing problem%s.\n",
      relaxdata->npricingprobs, nrelevant, (relaxdata->npricingprobs == nrelevant ? "" : " aggregated"), (nrelevant == 1 ? "" : "s"));

   relaxdata->nrelpricingprobs = nrelevant;

   if( relaxdata->npricingprobs > nrelevant )
   {
      // this is a workaround (GCG cannot handle different bounds on aggregated variables, see checkAggregatedLocalBounds)
      SCIP_CALL( SCIPsetBoolParam(scip, "misc/allowstrongdualreds", FALSE) );
      assert(!SCIPallowStrongDualReds(scip));
   }

   return SCIP_OKAY;
}

/** sets the pricing problem parameters */
SCIP_RETCODE GCGsetPricingProblemParameters(
   GCG_DECTYPE           dectype,            /**< the dectype of the decomp */
   SCIP*                 pricingprob,        /**< SCIP data structure of the pricing problem */
   int                   clocktype,          /**< clocktype to use in the pricing problem */
   SCIP_Real             infinity,           /**< values larger than this are considered infinity in the pricing problem */
   SCIP_Real             epsilon,            /**< absolute values smaller than this are considered zero in the pricing problem */
   SCIP_Real             sumepsilon,         /**< absolute values of sums smaller than this are considered zero in the pricing problem */
   SCIP_Real             feastol,            /**< feasibility tolerance for constraints in the pricing problem */
   SCIP_Real             lpfeastolfactor,    /**< primal feasibility tolerance factor of LP solver in the pricing problem */
   SCIP_Real             dualfeastol,        /**< feasibility tolerance for reduced costs in LP solution in the pricing problem */
   SCIP_Bool             enableppcuts        /**< should ppcuts be stored for sepa_basis */
   )
{
   assert(pricingprob != NULL);

   if( dectype != GCG_DECTYPE_DIAGONAL )
   {
      /* disable conflict analysis */
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "conflict/useprop", FALSE) );
      SCIP_CALL( SCIPsetCharParam(pricingprob, "conflict/useinflp", 'o') );
      SCIP_CALL( SCIPsetCharParam(pricingprob, "conflict/useboundlp", 'o') );
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "conflict/usesb", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "conflict/usepseudo", FALSE) );

      /* reduce the effort spent for hash tables */
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "misc/usevartable", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "misc/useconstable", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "misc/usesmalltables", TRUE) );

      /* disable expensive presolving */
      /* @todo test whether this really helps, perhaps set presolving emphasis to fast? */
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "constraints/linear/presolpairwise", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "constraints/setppc/presolpairwise", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "constraints/logicor/presolpairwise", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "constraints/linear/presolusehashing", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "constraints/setppc/presolusehashing", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "constraints/logicor/presolusehashing", FALSE) );

      /* disable dual fixing presolver for the moment (propagator should be safe), because we want to avoid variables fixed to infinity */
      SCIP_CALL( SCIPsetIntParam(pricingprob, "propagating/dualfix/maxprerounds", 0) );
      SCIP_CALL( SCIPfixParam(pricingprob, "propagating/dualfix/maxprerounds") );


      /* disable solution storage ! */
      SCIP_CALL( SCIPsetIntParam(pricingprob, "limits/maxorigsol", 0) );
      SCIP_CALL( SCIPfixParam(pricingprob, "limits/maxorigsol") );

      /* @todo enable presolving and propagation of xor constraints if bug is fixed */

      /* disable presolving and propagation of xor constraints as work-around for a SCIP bug */
      SCIP_CALL( SCIPsetIntParam(pricingprob, "constraints/xor/maxprerounds", 0) );
      SCIP_CALL( SCIPsetIntParam(pricingprob, "constraints/xor/propfreq", -1) );

      /* jonas' stuff */
      if( enableppcuts )
      {
         int pscost;
         int prop;

         SCIP_CALL( SCIPgetIntParam(pricingprob, "branching/pscost/priority", &pscost) );
         SCIP_CALL( SCIPgetIntParam(pricingprob, "propagating/maxroundsroot", &prop) );
         SCIP_CALL( SCIPsetIntParam(pricingprob, "branching/pscost/priority", 11000) );
         SCIP_CALL( SCIPsetIntParam(pricingprob, "propagating/maxroundsroot", 0) );
         SCIP_CALL( SCIPsetPresolving(pricingprob, SCIP_PARAMSETTING_OFF, TRUE) );
      }
   }

   /* disable multiaggregation because of infinite values */
   SCIP_CALL( SCIPsetBoolParam(pricingprob, "presolving/donotmultaggr", TRUE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(pricingprob, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
   #if SCIP_VERSION > 210
      SCIP_CALL( SCIPsetBoolParam(pricingprob, "misc/printreason", FALSE) );
   #endif

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(pricingprob, "misc/catchctrlc", FALSE) );

   /* set clock type */
   SCIP_CALL( SCIPsetIntParam(pricingprob, "timing/clocktype", clocktype) );

   SCIP_CALL( SCIPsetBoolParam(pricingprob, "misc/calcintegral", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(pricingprob, "misc/finitesolutionstore", TRUE) );

   SCIP_CALL( SCIPsetRealParam(pricingprob, "numerics/infinity", infinity) );
   SCIP_CALL( SCIPsetRealParam(pricingprob, "numerics/epsilon", epsilon) );
   SCIP_CALL( SCIPsetRealParam(pricingprob, "numerics/sumepsilon", sumepsilon) );
   SCIP_CALL( SCIPsetRealParam(pricingprob, "numerics/feastol", feastol) );
   SCIP_CALL( SCIPsetRealParam(pricingprob, "numerics/lpfeastolfactor", lpfeastolfactor) );
   SCIP_CALL( SCIPsetRealParam(pricingprob, "numerics/dualfeastol", dualfeastol) );

   return SCIP_OKAY;
}


/** creates a variable in a pricing problem corresponding to the given original variable (belonging to exactly one block) */
static
SCIP_RETCODE createPricingVar(
   SCIP_RELAXDATA*       relaxdata,          /**< relaxator data data structure */
   SCIP_VAR*             origvar             /**< corresponding variable in the original program */
   )
{
   SCIP_VAR* var;
   int pricingprobnr;

   assert(relaxdata != NULL);
   assert(origvar != NULL);

   pricingprobnr = GCGvarGetBlock(origvar);
   assert(pricingprobnr >= 0);

   SCIP_CALL( GCGoriginalVarCreatePricingVar(relaxdata->pricingprobs[pricingprobnr], origvar, &var) );
   assert(var != NULL);

   GCGoriginalVarSetPricingVar(origvar, var);
   SCIP_CALL( SCIPaddVar(relaxdata->pricingprobs[pricingprobnr], var) );
   assert(GCGvarIsPricing(var));
   /* because the variable was added to the problem,
    * it is captured by SCIP and we can safely release it right now
    */
   SCIP_CALL( SCIPreleaseVar(relaxdata->pricingprobs[pricingprobnr], &var) );

   return SCIP_OKAY;
}

/** creates a variable in each of the pricing problems linked by given original variable */
static
SCIP_RETCODE createLinkingPricingVars(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< relaxator data data structure */
   SCIP_VAR*             origvar             /**< corresponding linking variable in the original program */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_CONS* linkcons;
#ifndef NDEBUG
   SCIP_CONS** linkconss;
   int nblocks;
#endif
   SCIP_VAR** pricingvars;
   int i;

   assert(origvar != NULL);
   assert(relaxdata != NULL);

   scip = GCGgetOrigprob(gcg);

   /* get variable data of the original variable */
   assert(GCGvarIsOriginal(origvar));
   assert(GCGoriginalVarIsLinking(origvar));
   pricingvars = GCGlinkingVarGetPricingVars(origvar);

#ifndef NDEBUG
   nblocks = GCGlinkingVarGetNBlocks(origvar);
   /* checks that GCGrelaxSetOriginalVarBlockNr() worked correctly */
   {
      int count;

      linkconss = GCGlinkingVarGetLinkingConss(origvar);
      /* the linking constraints could be NULL if the Benders' decomposition is used. */
      if( linkconss != NULL )
      {
         count = 0;
         for( i = 0; i < relaxdata->npricingprobs; i++ )
         {
            assert(linkconss[i] == NULL);

            if( pricingvars[i] != NULL )
               count++;
         }
         assert(nblocks == count);
      }
   }
#endif

   for( i = 0; i < relaxdata->npricingprobs; ++i )
   {
      if( pricingvars[i] == NULL )
         continue;

      SCIP_CALL( GCGlinkingVarCreatePricingVar(relaxdata->pricingprobs[i], i, origvar, &var) );

      GCGlinkingVarSetPricingVar(origvar, i, var);

      assert(GCGvarIsPricing(var));
      SCIP_CALL( SCIPaddVar(relaxdata->pricingprobs[i], var) );


      if( relaxdata->mode != GCG_DECMODE_BENDERS )
      {
         SCIP_CALL( GCGlinkingVarCreateMasterCons(relaxdata->gcg, i, origvar, &linkcons) );
         GCGlinkingVarSetLinkingCons(origvar, linkcons, i);
         SCIP_CALL( SCIPaddCons(GCGgetMasterprob(relaxdata->gcg), linkcons) );

         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &relaxdata->varlinkconss, relaxdata->nvarlinkconss, (size_t)relaxdata->nvarlinkconss+1) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &relaxdata->varlinkconsblock, relaxdata->nvarlinkconss, (size_t)relaxdata->nvarlinkconss+1) );

         relaxdata->varlinkconss[relaxdata->nvarlinkconss] = linkcons;
         relaxdata->varlinkconsblock[relaxdata->nvarlinkconss] = i;
         relaxdata->nvarlinkconss++;
      }

      /* because the variable was added to the problem,
       * it is captured by SCIP and we can safely release it right now
       */
      SCIP_CALL( SCIPreleaseVar(relaxdata->pricingprobs[i], &var) );
   }

#ifndef NDEBUG
   /* checks that createLinkingPricingVars() worked correctly */
   {
      int count;

      linkconss = GCGlinkingVarGetLinkingConss(origvar);
      /* the linking constraints could be NULL if the Benders' decomposition is used. */
      if( linkconss != NULL )
      {
         count = 0;
         for( i = 0; i < relaxdata->npricingprobs; i++ )
         {
            if( pricingvars[i] != NULL )
            {
               count++;
               assert(GCGvarIsPricing(pricingvars[i]));
               assert(relaxdata->mode == GCG_DECMODE_BENDERS || linkconss[i] != NULL);
            }
            else
               assert(linkconss[i] == NULL);
         }
         assert(nblocks == count);
      }
   }
#endif


   return SCIP_OKAY;
}

/** create pricing problem variables */
static
SCIP_RETCODE createPricingVariables(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< relaxator data data structure */
   SCIP_HASHMAP**        hashorig2pricingvar /**< hashmap mapping original variables to pricing variables */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;
   int nvars;
   int v;
   int i;
   int npricingprobs;
#ifndef NDEBUG
   SCIP_HASHMAP* hashorig2origvar;
#endif

   assert(gcg != NULL);
   assert(relaxdata != NULL);

   scip = GCGgetOrigprob(gcg);

   /* create pricing variables and map them to the original variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   npricingprobs = relaxdata->npricingprobs;

#ifndef NDEBUG
   SCIP_CALL( SCIPhashmapCreate(&hashorig2origvar, SCIPblkmem(scip), 10*SCIPgetNVars(scip)+1) );
#endif

   for( v = 0; v < nvars; v++ )
   {
      int blocknr;
      SCIP_VAR* probvar;

      assert(SCIPvarIsTransformed(vars[v]));

      probvar = SCIPvarGetProbvar(vars[v]);
      assert(SCIPvarIsTransformed(probvar));
      blocknr = GCGvarGetBlock(probvar);
      if( blocknr == -1 )
      {
         int tempblock;
         tempblock = (int) (size_t) SCIPhashmapGetImage(GCGdecompGetVartoblock(relaxdata->decomp), probvar)-1; /*lint !e507*/
         if( tempblock >= GCGdecompGetNBlocks(relaxdata->decomp) )
         {
            blocknr = -1;
         }
         else
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " changed block number to %d \n", tempblock );
            blocknr = tempblock; /*lint !e806*/
         }
      }

#ifndef NDEBUG
      SCIPdebugMessage("Creating map for (%p, %p) var %s:", (void*)(vars[v]), (void*)(probvar), SCIPvarGetName(probvar));
      assert( !SCIPhashmapExists(hashorig2origvar, probvar) );
      SCIP_CALL( SCIPhashmapInsert(hashorig2origvar, (void*)(probvar), (void*)(probvar)) );
#endif

      /* variable belongs to exactly one block --> create corresponding pricing variable*/
      if( blocknr >= 0 )
      {
         SCIPdebugPrintf("block %d", blocknr);

         assert(GCGoriginalVarGetPricingVar(probvar) == NULL);
         SCIP_CALL( createPricingVar(relaxdata, probvar) );
         assert(GCGoriginalVarGetPricingVar(probvar) != NULL);
         assert(hashorig2pricingvar != NULL);
         assert(hashorig2pricingvar[blocknr] != NULL);

         SCIPdebugPrintf("-> %p\n", (void*) GCGoriginalVarGetPricingVar(probvar));

         assert(!SCIPhashmapExists(hashorig2pricingvar[blocknr], probvar));
         SCIP_CALL( SCIPhashmapInsert(hashorig2pricingvar[blocknr], (void*)(probvar),
               (void*)(GCGoriginalVarGetPricingVar(probvar)) ));

         assert(GCGvarIsPricing((SCIP_VAR*) SCIPhashmapGetImage(hashorig2pricingvar[blocknr], probvar)));
      }
      /* variable is a linking variable --> create corresponding pricing variable in all linked blocks
       * and create corresponding linking constraints */
      else if( GCGoriginalVarIsLinking(probvar) )
      {
         SCIP_VAR** pricingvars;
         SCIPdebugPrintf("linking.\n");
         relaxdata->nlinkingvars++;
         SCIP_CALL( createLinkingPricingVars(gcg, relaxdata, probvar) );
         assert(GCGlinkingVarGetPricingVars(probvar) != NULL);


         pricingvars = GCGlinkingVarGetPricingVars(probvar);

         for( i = 0; i < npricingprobs; i++ )
         {
            if( pricingvars[i] != NULL )
            {
               assert(GCGvarIsPricing(pricingvars[i]));
               assert(hashorig2pricingvar != NULL);
               assert(hashorig2pricingvar[i] != NULL);
               assert(!SCIPhashmapExists(hashorig2pricingvar[i], probvar));
               SCIP_CALL( SCIPhashmapInsert(hashorig2pricingvar[i], (void*)(probvar),
                     (void*)(pricingvars[i])) );
               assert(GCGvarIsPricing((SCIP_VAR*) SCIPhashmapGetImage(hashorig2pricingvar[i], probvar)));
            }
         }
      }
      else
      {
         assert(GCGvarGetBlock(probvar) == -1);
         assert(GCGoriginalVarGetPricingVar(probvar) == NULL);
         SCIPdebugPrintf("master!\n");
         relaxdata->ntransvars++;
      }
      assert(SCIPhashmapExists(hashorig2origvar, probvar));
   }

#ifndef NDEBUG
   SCIPhashmapFree(&hashorig2origvar);
#endif

   return SCIP_OKAY;
}


/** displays statistics of the pricing problems */
static
SCIP_RETCODE displayPricingStatistics(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP**                pricingprobs,       /**< array of pricing problems */
   int                   npricingprobs,      /**< number of pricingproblems */
   int*                  blockrepresentative /**< array of representation information */
)
{
   SCIP* scip;
   char name[SCIP_MAXSTRLEN];
   int i;

   assert(gcg != NULL);
   assert(pricingprobs != NULL);
   assert(blockrepresentative != NULL);
   assert(npricingprobs > 0);

   scip = GCGgetOrigprob(gcg);

   for( i = 0; i < npricingprobs; i++ )
   {
      int nbin;
      int nint;
      int nimpl;
      int ncont;

      if( blockrepresentative[i] != i )
         continue;

      SCIP_CALL( SCIPgetVarsData(pricingprobs[i], NULL, NULL, &nbin, &nint, &nimpl, &ncont) );

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "pricing problem %d: %d conss, %d vars (%d bins, %d ints, %d impls and %d cont)\n", i,
         SCIPgetNConss(pricingprobs[i]), SCIPgetNVars(pricingprobs[i]), nbin, nint, nimpl, ncont);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricingprob_%d.lp", i);
      SCIP_CALL( SCIPwriteOrigProblem(pricingprobs[i], name, NULL, FALSE) );
   }

   return SCIP_OKAY;
}


/** allocates initial problem specific data */
static
SCIP_RETCODE initRelaxProblemdata(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata           /**< relaxatordata data structure */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   assert(relaxdata != NULL);

   scip = GCGgetOrigprob(gcg);

   /* initialize relaxator data */
   relaxdata->maxmasterconss = 16;
   assert(relaxdata->nmasterconss == 0);

   /* arrays of constraints belonging to the master problems */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(relaxdata->masterconss), relaxdata->maxmasterconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(relaxdata->origmasterconss), relaxdata->maxmasterconss) );

   if( relaxdata->npricingprobs > relaxdata->maxpricingprobs )
   {
      int i;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(relaxdata->pricingprobs), relaxdata->maxpricingprobs, relaxdata->npricingprobs) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(relaxdata->blockrepresentative), relaxdata->maxpricingprobs, relaxdata->npricingprobs) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(relaxdata->nblocksidentical), relaxdata->maxpricingprobs, relaxdata->npricingprobs) );

      /* array for saving convexity constraints belonging to one of the pricing problems */
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(relaxdata->convconss), relaxdata->maxpricingprobs, relaxdata->npricingprobs) );

      for( i = relaxdata->maxpricingprobs; i < relaxdata->npricingprobs; ++i)
      {
         relaxdata->pricingprobs[i] = NULL;
      }

       relaxdata->maxpricingprobs = relaxdata->npricingprobs;
   }

   return SCIP_OKAY;
}


/** creates the master problem with the specified name */
static
SCIP_RETCODE createMasterProblem(
   GCG*                  gcg,                /**< GCG data structure */
   const char*           name,               /**< name of the master problem */
   int                   clocktype,          /**< clocktype to use in the master problem */
   SCIP_Real             infinity,           /**< values larger than this are considered infinity in the master problem */
   SCIP_Real             epsilon,            /**< absolute values smaller than this are considered zero in the master problem */
   SCIP_Real             sumepsilon,         /**< absolute values of sums smaller than this are considered zero in the master problem */
   SCIP_Real             feastol,            /**< feasibility tolerance for constraints in the master problem */
   SCIP_Real             lpfeastolfactor,    /**< primal feasibility tolerance factor of LP solver in the master problem */
   SCIP_Real             dualfeastol,        /**< feasibility tolerance for reduced costs in LP solution in the master problem */
   GCG_DECMODE           mode                /**< the decomposition mode */
   )
{
   SCIP* masterprob = GCGgetMasterprob(gcg);
   assert(gcg != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPcreateProb(masterprob, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* set clocktype */
   SCIP_CALL( SCIPsetIntParam(masterprob, "timing/clocktype", clocktype) );

   /* set numerical tolerances */
   SCIP_CALL( SCIPsetRealParam(masterprob, "numerics/infinity", infinity) );
   SCIP_CALL( SCIPsetRealParam(masterprob, "numerics/epsilon", epsilon) );
   SCIP_CALL( SCIPsetRealParam(masterprob, "numerics/sumepsilon", sumepsilon) );
   SCIP_CALL( SCIPsetRealParam(masterprob, "numerics/feastol", feastol) );
   SCIP_CALL( SCIPsetRealParam(masterprob, "numerics/lpfeastolfactor", lpfeastolfactor) );
   SCIP_CALL( SCIPsetRealParam(masterprob, "numerics/dualfeastol", dualfeastol) );

   /* disable aggregation and multiaggregation of variables, as this might lead to issues with copied original variables */
   SCIP_CALL( SCIPsetBoolParam(masterprob, "presolving/donotaggr", TRUE) );
   SCIP_CALL( SCIPsetBoolParam(masterprob, "presolving/donotmultaggr", TRUE) );

   /* do not catch ctrl-c @todo: add this feature*/
   SCIP_CALL( SCIPsetBoolParam(masterprob, "misc/catchctrlc", FALSE) );

   /* the following settings are for decomposition, so if the original problem is solved directly, then these settings
    * are not required
    */
   if( mode == GCG_DECMODE_ORIGINAL )
   {
      return SCIP_OKAY;
   }

   if( mode == GCG_DECMODE_DANTZIGWOLFE )
      SCIP_CALL( SCIPactivatePricer(masterprob, SCIPfindPricer(masterprob, "gcg")) );

   /* do not modify the time limit after solving the master problem */
   SCIP_CALL( SCIPsetBoolParam(masterprob, "reoptimization/commontimelimit", FALSE) );

   /* for Benders' decomposition, some additional parameter settings are required for the master problem. */
   if( mode == GCG_DECMODE_BENDERS )
   {
      SCIP_CALL( SCIPsetSeparating(masterprob, SCIP_PARAMSETTING_OFF, TRUE) );
      SCIP_CALL( SCIPsetPresolving(masterprob, SCIP_PARAMSETTING_OFF, TRUE) );
      SCIP_CALL( SCIPsetIntParam(masterprob, "presolving/maxrestarts", 0) );
      SCIP_CALL( SCIPsetIntParam(masterprob, "propagating/maxroundsroot", 0) );
      SCIP_CALL( SCIPsetIntParam(masterprob, "heuristics/trysol/freq", 1) );
      SCIP_CALL( SCIPsetBoolParam(masterprob, "constraints/benders/active", TRUE) );
      SCIP_CALL( SCIPsetBoolParam(masterprob, "constraints/benderslp/active", TRUE) );
      SCIP_CALL( SCIPsetBoolParam(masterprob, "benders/gcg/lnscheck", FALSE) );
      SCIP_CALL( SCIPsetIntParam(masterprob, "presolving/maxrounds", 1) );
      SCIP_CALL( SCIPsetIntParam(masterprob, "constraints/benders/maxprerounds", 1) );

      /* the trysol heuristic must have a high priority to ensure the solutions found by the relaxator are added to the
       * original problem
       */
      SCIP_CALL( SCIPsetIntParam(GCGgetOrigprob(gcg), "heuristics/trysol/freq", 1) );

      /* disabling pricing problem aggregation */
      SCIP_CALL( SCIPsetBoolParam(GCGgetOrigprob(gcg), "relaxing/gcg/aggregation", FALSE) );
   }

   return SCIP_OKAY;
}


/** creates the pricing problem with the specified name */
static
SCIP_RETCODE createPricingProblem(
   SCIP_RELAXDATA*       relaxdata,          /**< the relaxator data data structure */
   SCIP**                pricingscip,        /**< Pricing GCG data structure */
   const char*           name,               /**< name of the pricing problem */
   int                   clocktype,          /**< clocktype to use in the pricing problem */
   SCIP_Real             infinity,           /**< values larger than this are considered infinity in the pricing problem */
   SCIP_Real             epsilon,            /**< absolute values smaller than this are considered zero in the pricing problem */
   SCIP_Real             sumepsilon,         /**< absolute values of sums smaller than this are considered zero in the pricing problem */
   SCIP_Real             feastol,            /**< feasibility tolerance for constraints in the pricing problem */
   SCIP_Real             lpfeastolfactor,    /**< primal feasibility tolerance factor of LP solver in the pricing problem */
   SCIP_Real             dualfeastol,        /**< feasibility tolerance for reduced costs in LP solution in the pricing problem */
   SCIP_Bool             enableppcuts        /**< should ppcuts be stored for sepa_basis */
   )
{
   assert(pricingscip != NULL);
   assert(name != NULL);
   assert(relaxdata->mode != GCG_DECMODE_ORIGINAL);

   if( *pricingscip == NULL )
   {
      SCIP_CALL( SCIPcreate(pricingscip) );
      SCIP_CALL( SCIPincludeDefaultPlugins(*pricingscip) );
      SCIP_CALL( GCGsetPricingProblemParameters(GCGdecompGetType(relaxdata->decomp), *pricingscip, clocktype, infinity, epsilon, sumepsilon, feastol, lpfeastolfactor, dualfeastol, enableppcuts) );
   }
   SCIP_CALL( SCIPcreateProb(*pricingscip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}


/** saves the coefficient of the masterconstraints in the original variable */
static
SCIP_RETCODE saveOriginalVarMastercoeffs(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_VAR**            origvars,           /**< original variables array */
   int                   norigvars,          /**< size of original variables array*/
   int                   nmasterconss,       /**< size of masterconns array */
   SCIP_CONS**           origmasterconss,    /**< orig master constraints array */
   SCIP_CONS**           masterconss         /**< master constraints */
   )
{
   SCIP* scip;
   int v;
   int i;

   assert(gcg != NULL);
   assert(origvars != NULL || norigvars == 0);
   assert(norigvars >= 0);
   assert(nmasterconss >= 0);
   assert(masterconss != NULL);
   assert(origmasterconss != NULL);

   scip = GCGgetOrigprob(gcg);

   /* for original variables, save the coefficients in the master problem */
   for( v = 0; v < norigvars; v++ )
   {
      SCIP_VAR* var;
      var = SCIPvarGetProbvar(origvars[v]); /*lint !e613*/
      assert(GCGvarIsOriginal(var));
      assert(GCGoriginalVarGetCoefs(var) == NULL);
      GCGoriginalVarSetNCoefs(var, 0);
   }

   /* save coefs */
   for( i = 0; i < nmasterconss; i++ )
   {
      SCIP_VAR** vars;
      SCIP_Real* vals;
      int nvars;

      nvars = GCGconsGetNVars(scip, origmasterconss[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
      GCGconsGetVars(scip, origmasterconss[i], vars, nvars);
      GCGconsGetVals(scip, origmasterconss[i], vals, nvars);
      for( v = 0; v < nvars; v++ )
      {
         SCIP_CALL( GCGoriginalVarAddCoef(gcg, vars[v], vals[v], masterconss[i]) );
      }
      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &vars);
   }

   return SCIP_OKAY;
}

/** creates the master problem constraints */
static
SCIP_RETCODE createMasterprobConss(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata           /**< the relaxator data data structure */
   )
{
   SCIP* scip;
   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_CONS* mastercons;
   int c;
   char name[SCIP_MAXSTRLEN];

   scip = GCGgetOrigprob(gcg);
   masterconss = GCGdecompGetLinkingconss(relaxdata->decomp);
   nmasterconss = GCGdecompGetNLinkingconss(relaxdata->decomp);

 //  assert(SCIPhashmapGetNElements(relaxdata->hashorig2origvar) == SCIPgetNVars(scip));
   for( c = 0; c < nmasterconss; ++c )
   {
      int nconsvars;
      int consvarssize;
      SCIP_VAR** consvars;
      SCIP_Real* consvals;
      SCIP_Bool* releasevars;
      int i;

      if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(masterconss[c])), "origbranch") == 0 )
         continue;

      /* in the Benders' decomposition mode, all variables from the linking constraints need to be added to the master
       * problem. Additionally, if the original problem is solved directly, then we must ensure that all variables are
       * added to the master problem.
       */
      if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_BENDERS || GCGgetDecompositionMode(gcg) == GCG_DECMODE_ORIGINAL )
      {
         nconsvars = GCGconsGetNVars(scip, masterconss[c]);
         consvarssize = nconsvars;

         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, consvarssize) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, consvarssize) );
         SCIP_CALL( SCIPallocClearBufferArray(scip, &releasevars, consvarssize) );

         SCIP_CALL( GCGconsGetVars(scip, masterconss[c], consvars, nconsvars) );
         SCIP_CALL( GCGconsGetVals(scip, masterconss[c], consvals, nconsvars) );

         for( i = 0; i < nconsvars; i++ )
         {
            SCIP_VAR* origvar;

            /* if the variable is a linking variable or is directly transferred to the master problem, then it is not
             * added to the constraint. This is because the linking variables and the transferred variables are added
             * later in GCGmasterCreateInitialMastervars().
             */
            while( i < nconsvars && (GCGoriginalVarIsLinking(consvars[i]) || GCGoriginalVarIsTransVar(consvars[i])) )
            {
               consvars[i] = consvars[nconsvars - 1];
               consvals[i] = consvals[nconsvars - 1];
               nconsvars--;
            }

            if( i >= nconsvars )
               break;

            /* assigning the origvar to the next variables that is not a linking variable */
            origvar = consvars[i];

            assert(GCGoriginalVarGetNMastervars(origvar) <= 1);

            /* if the original has already has a copy in the master problem, then this is used. Otherwise, the master
             * problem variable is created.
             */
            if( GCGoriginalVarGetNMastervars(origvar) > 0 )
            {
               consvars[i] = GCGoriginalVarGetMastervars(origvar)[0];
               releasevars[i] = FALSE;
            }
            else
            {
               SCIP_CALL( GCGcreateInitialMasterVar(gcg, consvars[i], &consvars[i]) );
               SCIP_CALL( SCIPaddVar(GCGgetMasterprob(gcg), consvars[i]) );

               SCIP_CALL( GCGoriginalVarAddMasterVar(gcg, origvar, consvars[i], 1.0) );

               releasevars[i] = TRUE;
            }

            assert(GCGoriginalVarGetNMastervars(origvar) <= 1);
         }
      }
      else
      {
         nconsvars = 0;
         consvars = NULL;
         consvals = NULL;
         releasevars = NULL;
      }

      /* create and add corresponding linear constraint in the master problem */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "m_%s", SCIPconsGetName(masterconss[c]));
      SCIP_CALL( SCIPcreateConsLinear(GCGgetMasterprob(relaxdata->gcg), &mastercons, name, nconsvars, consvars, consvals,
            GCGconsGetLhs(scip, masterconss[c]), GCGconsGetRhs(scip, masterconss[c]),
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(GCGgetMasterprob(relaxdata->gcg), mastercons) );
      SCIPdebugMessage("Copying %s to masterproblem\n", SCIPconsGetName(masterconss[c]));
      /* store the constraints in the arrays origmasterconss and masterconss in the problem data */
      SCIP_CALL( ensureSizeMasterConss(gcg, relaxdata, relaxdata->nmasterconss+1) );
      SCIP_CALL( SCIPcaptureCons(scip, masterconss[c]) );
      relaxdata->origmasterconss[relaxdata->nmasterconss] = masterconss[c];
      relaxdata->masterconss[relaxdata->nmasterconss] = mastercons;
      relaxdata->nmasterconss++;

      /* in the Benders' decomposition mode, the consvars and consvals arrays need to be freed */
      if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_BENDERS || GCGgetDecompositionMode(gcg) == GCG_DECMODE_ORIGINAL )
      {
         assert(releasevars != NULL);
         assert(consvars != NULL);
         for( i = 0; i < nconsvars; i++ )
         {
            if( releasevars[i] )
            {
               SCIP_CALL( SCIPreleaseVar(GCGgetMasterprob(relaxdata->gcg), &consvars[i]) );
            }
         }

         SCIPfreeBufferArray(scip, &releasevars);
         SCIPfreeBufferArray(scip, &consvals);
         SCIPfreeBufferArray(scip, &consvars);
      }
   }
   assert(relaxdata->nmasterconss == nmasterconss);
   SCIP_CALL( saveOriginalVarMastercoeffs(gcg, SCIPgetVars(scip), SCIPgetNVars(scip), relaxdata->nmasterconss, relaxdata->origmasterconss, relaxdata->masterconss) );

   return SCIP_OKAY;
}

/** creates the pricing problem constraints */
static
SCIP_RETCODE createPricingprobConss(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< the relaxator data data structure */
   SCIP_HASHMAP**        hashorig2pricingvar /**< hashmap mapping original to corresponding pricing variables */
   )
{
   SCIP* scip;
   SCIP_CONS*** subscipconss;
   int* nsubscipconss;
   SCIP_CONS* newcons;
   SCIP_HASHMAP* hashorig2pricingconstmp;
   int nblocks;
   int b;
   int c;
   char name[SCIP_MAXSTRLEN];
   SCIP_Bool success;

   assert(gcg != NULL);
   assert(relaxdata != NULL);

   scip = GCGgetOrigprob(gcg);
   subscipconss = GCGdecompGetSubscipconss(relaxdata->decomp);
   nsubscipconss = GCGdecompGetNSubscipconss(relaxdata->decomp);
   nblocks = GCGdecompGetNBlocks(relaxdata->decomp);

   SCIP_CALL( SCIPhashmapCreate(&hashorig2pricingconstmp, SCIPblkmem(scip), SCIPgetNConss(scip)) ); /*lint !e613*/

   for( b = 0; b < nblocks; ++b )
   {
      assert(hashorig2pricingvar != NULL);
      for( c = 0; c < nsubscipconss[b]; ++c )
      {
         SCIPdebugMessage("copying %s to pricing problem %d\n", SCIPconsGetName(subscipconss[b][c]), b);
         if( !SCIPconsIsActive(subscipconss[b][c]) )
         {
            SCIPdebugMessage("skipping, cons <%s> inactive\n", SCIPconsGetName(subscipconss[b][c]));
            continue;
         }
         SCIP_CALL( SCIPgetTransformedCons(scip, subscipconss[b][c], &subscipconss[b][c]) );
         assert(subscipconss[b][c] != NULL);

         /* copy the constraint */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "p%d_%s", b, SCIPconsGetName(subscipconss[b][c]));
         SCIP_CALL( SCIPgetConsCopy(scip, relaxdata->pricingprobs[b], subscipconss[b][c], &newcons, SCIPconsGetHdlr(subscipconss[b][c]),
               hashorig2pricingvar[b], hashorig2pricingconstmp, name,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, &success) );

         /* constraint was successfully copied */
         assert(success);

         SCIP_CALL( SCIPaddCons(relaxdata->pricingprobs[b], newcons) );


#ifndef NDEBUG
         {
            SCIP_VAR** curvars;
            int ncurvars;

            ncurvars = GCGconsGetNVars(relaxdata->pricingprobs[b], newcons);
            curvars = NULL;
            if( ncurvars > 0 )
            {
               int i;

               SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
               SCIP_CALL( GCGconsGetVars(relaxdata->pricingprobs[b], newcons, curvars, ncurvars) );

               for( i = 0; i < ncurvars; ++i )
               {
                  if( SCIPisFeasEQ( scip, SCIPvarGetLbGlobal(curvars[i]), SCIPvarGetUbGlobal(curvars[i]) ) && SCIPisFeasEQ( scip, SCIPvarGetUbGlobal(curvars[i]), 0. )  )
                     continue;

                  assert(GCGvarIsPricing(curvars[i]) || ( SCIPvarIsNegated(curvars[i]) && GCGvarIsPricing(SCIPvarGetNegatedVar(curvars[i]) ) ) );
               }

               SCIPfreeBufferArrayNull(scip, &curvars);
            }
         }
#endif
         SCIP_CALL( SCIPreleaseCons(relaxdata->pricingprobs[b], &newcons) );
      }
   }

   SCIPhashmapFree(&hashorig2pricingconstmp);

   return SCIP_OKAY;
}

/** creates the master problem and the pricing problems and copies the constraints into them */
static
SCIP_RETCODE createMaster(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata           /**< the relaxator data data structure */
   )
{
   SCIP* origprob;
   int npricingprobs;
   SCIP_HASHMAP** hashorig2pricingvar;
   SCIP_Bool enableppcuts;
   char name[SCIP_MAXSTRLEN];
   int clocktype;
   SCIP_Real infinity;
   SCIP_Real epsilon;
   SCIP_Real sumepsilon;
   SCIP_Real feastol;
   SCIP_Real lpfeastolfactor;
   SCIP_Real dualfeastol;
   int i;

   assert(gcg != NULL);
   assert(relaxdata != NULL);

   assert(relaxdata->decomp != NULL);

   origprob = GCGgetOrigprob(gcg);

   SCIP_CALL( convertStructToGCG(gcg, relaxdata, relaxdata->decomp) );

   /* if there are no pricing problems, then the original problem will be solved directly. */
   if( relaxdata->npricingprobs == 0 )
   {
      int origmode = relaxdata->mode;

      /* setting the mode to ORIGINAL */
      relaxdata->mode = GCG_DECMODE_ORIGINAL;
      SCIP_CALL( SCIPfixParam(origprob, "relaxing/gcg/mode") );

      if( origmode == GCG_DECMODE_DANTZIGWOLFE )
      {
         /* initialising the master problem */
         SCIP_CALL( SCIPsetIntParam(gcg->bendersmasterprob, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
         SCIP_CALL( SCIPsetBoolParam(gcg->bendersmasterprob, "display/relevantstats", FALSE) );

         /* disabling unnecessary display columns */
         SCIP_CALL( SCIPsetIntParam(origprob, "display/sumlpiterations/active", 0) );
         SCIP_CALL( SCIPsetIntParam(origprob, "display/lpiterations/active", 0) );
         SCIP_CALL( SCIPsetIntParam(origprob, "display/degeneracy/active", 0) );

         /* setting the total node limit to 1 for the original SCIP instance. This is because Benders' decomposition solves
          * the MIP within the relaxator of the root node. So no branching in the original problem is required.
          */
         SCIP_CALL( SCIPsetLongintParam(origprob, "limits/totalnodes", 1) );

         /* swapping the master problem with the original master problem */
         relaxdata->gcg->masterprob = gcg->bendersmasterprob;
      }

      SCIP_CALL( SCIPsetIntParam(GCGgetMasterprob(relaxdata->gcg), "constraints/components/maxprerounds", 0) );
      SCIP_CALL( SCIPsetBoolParam(origprob, "relaxing/gcg/discretization", FALSE) );
   }


   npricingprobs = relaxdata->npricingprobs;
   hashorig2pricingvar = NULL;

   if( npricingprobs > 0 )
   {
      /* create hashmaps for mapping from original to pricing variables */
      SCIP_CALL( SCIPallocBufferArray(origprob, &(hashorig2pricingvar), npricingprobs) );
   }

   SCIPdebugMessage("Creating master problem...\n");

   SCIP_CALL( initRelaxProblemdata(gcg, relaxdata) );

   /* get clocktype of the original SCIP instance in order to use the same clocktype in master and pricing problems */
   SCIP_CALL( SCIPgetIntParam(origprob, "timing/clocktype", &clocktype) );

   /* get numerical tolerances of the original SCIP instance in order to use the same numerical tolerances in master and pricing problems */
   SCIP_CALL( SCIPgetRealParam(origprob, "numerics/infinity", &infinity) );
   SCIP_CALL( SCIPgetRealParam(origprob, "numerics/epsilon", &epsilon) );
   SCIP_CALL( SCIPgetRealParam(origprob, "numerics/sumepsilon", &sumepsilon) );
   SCIP_CALL( SCIPgetRealParam(origprob, "numerics/feastol", &feastol) );
   SCIP_CALL( SCIPgetRealParam(origprob, "numerics/lpfeastolfactor", &lpfeastolfactor) );
   SCIP_CALL( SCIPgetRealParam(origprob, "numerics/dualfeastol", &dualfeastol) );

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "master_%s", SCIPgetProbName(origprob));
   SCIP_CALL( createMasterProblem(gcg, name, clocktype, infinity, epsilon, sumepsilon, feastol,
         lpfeastolfactor, dualfeastol, relaxdata->mode) );

   enableppcuts = FALSE;
   SCIP_CALL( SCIPgetBoolParam(origprob, "sepa/basis/enableppcuts", &enableppcuts) );

   /* create the pricing problems */
   for( i = 0; i < npricingprobs; i++ )
   {
      relaxdata->convconss[i] = NULL;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricing_block_%d", i);
      SCIP_CALL( createPricingProblem(relaxdata, &(relaxdata->pricingprobs[i]), name, clocktype, infinity, epsilon, sumepsilon,
            feastol, lpfeastolfactor, dualfeastol, enableppcuts) );
      SCIP_CALL( SCIPhashmapCreate(&(hashorig2pricingvar[i]), SCIPblkmem(origprob), SCIPgetNVars(origprob)) ); /*lint !e613*/

      /* disabling restarts from the tree size estimation */
      SCIP_CALL( SCIPsetCharParam(relaxdata->pricingprobs[i], "estimation/restarts/restartpolicy", 'n') );
   }

   SCIP_CALL( createPricingVariables(gcg, relaxdata, hashorig2pricingvar) );

   /* create master and pricing problem constraints
    * If the master problem is solved directly, then we can still call methods creating the pricing problems. These
    * methods check the number of pricing problems and number of blocks.  As such, if the original problem is solved
    * directly, then nothing will happen in these methods
    */
   SCIP_CALL( createMasterprobConss(gcg, relaxdata) );
   SCIP_CALL( createPricingprobConss(gcg, relaxdata, hashorig2pricingvar) );
   SCIP_CALL( GCGmasterCreateInitialMastervars(gcg) );

   /* check if the master problem is a set partitioning or set covering problem */
   SCIP_CALL( checkSetppcStructure(gcg, relaxdata) );

   /* check for identity of blocks */
   SCIP_CALL( checkIdenticalBlocks(gcg, relaxdata, hashorig2pricingvar) );

   /* the convexity constraints are only added in the Dantzig-Wolfe mode */
   if( relaxdata->mode == GCG_DECMODE_DANTZIGWOLFE )
   {
      for( i = 0; i < relaxdata->npricingprobs; i++ )
      {
         if( relaxdata->blockrepresentative[i] != i )
            continue;

         /* create the corresponding convexity constraint */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conv_block_%d", i);
         SCIP_CALL( SCIPcreateConsLinear(GCGgetMasterprob(relaxdata->gcg), &(relaxdata->convconss[i]), name, 0, NULL, NULL,
               relaxdata->nblocksidentical[i]*1.0, relaxdata->nblocksidentical[i]*1.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(GCGgetMasterprob(relaxdata->gcg), relaxdata->convconss[i]) );
      }
   }

   /* display statistics */
   if( relaxdata->dispinfos )
   {
      SCIP_CALL( displayPricingStatistics(gcg, relaxdata->pricingprobs, relaxdata->npricingprobs, relaxdata->blockrepresentative) );
      SCIP_CALL( SCIPwriteOrigProblem(GCGgetMasterprob(relaxdata->gcg), "masterprob.lp", "lp", FALSE) );
   }

   if( hashorig2pricingvar != NULL )
   {
      for( i = 0; i < npricingprobs; i++ )
         SCIPhashmapFree(&(hashorig2pricingvar[i]));

      SCIPfreeBufferArray(origprob, &(hashorig2pricingvar));
   }

   /* get used memory and save it for reference */
   for( i = 0; i < npricingprobs; ++i )
   {
      relaxdata->pricingprobsmemused += SCIPgetMemUsed(relaxdata->pricingprobs[i])/1048576.0;
   }

   return SCIP_OKAY;
}

/** combines the solutions from all (disjoint) problems to one solution */
static
SCIP_RETCODE combineSolutions(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_SOL**            newsol,             /**< pointer to store new solution */
   SCIP**                probs,              /**< array of (solved) subproblems */
   int                   nprobs              /**< number of subproblems */
   )
{
   SCIP* scip;
#ifdef SCIP_DEBUG
   int i;
#endif

   int v;
   int nvars;

   SCIP_VAR** vars;
   assert(gcg != NULL);
   assert(newsol != NULL);
   assert(probs != NULL);
   assert(nprobs > 0);

   scip = GCGgetOrigprob(gcg);

   SCIP_CALL( SCIPcreateSol(scip, newsol, NULL) );
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

#ifdef SCIP_DEBUG
   for( i = 0; i < nprobs; ++i )
   {
      if( probs[i] == NULL )
         continue;

      SCIPprintOrigProblem(probs[i], NULL, "lp", FALSE);
      SCIPprintSol(probs[i], SCIPgetBestSol(probs[i]), NULL, FALSE );
   }
#endif

   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* pricingvar;
      int block;

      pricingvar = GCGoriginalVarGetPricingVar(vars[v]);
      block = GCGvarGetBlock(pricingvar);
      assert(block >= 0);
      assert(block < nprobs);
      assert(probs[block] != NULL);

      /* @todo solval should be 0 before, anyway, check it with an assert */
      SCIP_CALL( SCIPincSolVal(scip, *newsol, vars[v], SCIPgetSolVal(probs[block], SCIPgetBestSol(probs[block]), pricingvar)) );
   }
   return SCIP_OKAY;
}

/** sets the pricing objective function to what is necessary */
static
SCIP_RETCODE setPricingObjsOriginal(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP**                probs,              /**< array of subproblems */
   int                   nprobs              /**< number of subproblems */
   )
{
   SCIP* scip;
   int v;
   int nvars;
   SCIP_VAR** vars;
   int i;

   assert(gcg != NULL);
   assert(probs != NULL);
   assert(nprobs > 0);

   scip = GCGgetOrigprob(gcg);
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   /* if the Benders' decomposition is used, then the transformed problem of the subproblems must be freed.
    * This is because the within the create subproblem stage, if the subproblem is an LP, then the SCIP instance is put
    * into probing mode.
    */
   if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_BENDERS )
   {
      for( i = 0; i < nprobs; i++ )
      {
         /* if the problem is not in SCIP_STAGE_PROBLEM, then the transformed problem must be freed. The subproblem
          * should also be in probing mode.
          */
         if( SCIPgetStage(probs[i]) != SCIP_STAGE_PROBLEM )
         {
            if( SCIPinProbing(probs[i]) )
            {
               SCIP_CALL( SCIPendProbing(probs[i]) );
            }

            SCIP_CALL( SCIPfreeTransform(probs[i]) );
         }
      }
   }

   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* pricingvar;
      SCIP_VAR* origvar;
      SCIP_Real objvalue;

      assert(GCGvarIsOriginal(vars[v]));
      origvar = SCIPvarGetProbvar(vars[v]);

      if( !GCGisPricingprobRelevant(gcg, GCGvarGetBlock(origvar)) )
         continue;

      pricingvar = GCGoriginalVarGetPricingVar(origvar);
      assert(pricingvar != NULL);

      objvalue = SCIPvarGetObj(origvar);
      /* SCIPinfoMessage(scip, NULL, "%s: %f block %d\n", SCIPvarGetName(origvar), SCIPvarGetObj(origvar),
         GCGvarGetBlock(origvar)); */
      SCIP_CALL( SCIPchgVarObj(probs[GCGvarGetBlock(pricingvar)], pricingvar, objvalue) );
   }
   return SCIP_OKAY;
}

/** solve a block problem when the decomposition is diagonal */
static
SCIP_RETCODE solveBlockProblem(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP*                 blockprob,          /**< the block problem that will be solved */
   SCIP_RELAXDATA*       relaxdata,          /**< relaxator data structure */
   SCIP_Real             timelimit,          /**< the original problem timelimit */
   int                   blocknum,           /**< the number of the block, -1 for the original problem */
   SCIP_RESULT*          result,             /**< result pointer to indicate success or failure */
   SCIP_Real*            objvalue            /**< the objective function value */
   )
{
   SCIP* scip;
   SCIP_Real blocktimelimit;
   SCIP_STATUS blockprobstatus = SCIP_STATUS_UNKNOWN;

   assert(gcg != NULL);
   assert(result != NULL);
   assert(objvalue != NULL);

   scip = GCGgetOrigprob(gcg);
   (*result) = SCIP_DIDNOTRUN;

#ifdef SCIP_DEBUG
   char name[SCIP_MAXSTRLEN];
#endif

   if( blockprob == NULL )
   {
      (*result) = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   if( GCGgetDecompositionMode(gcg) != GCG_DECMODE_ORIGINAL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Solving block %i.\n", blocknum+1);
   }

   SCIP_CALL( SCIPsetIntParam(blockprob, "display/verblevel", relaxdata->origverblevel) );

   /* give the pricing problem 2% more time then the original scip has left */
   if( SCIPgetStage(blockprob) > SCIP_STAGE_PROBLEM )
   {
      if( SCIPisInfinity(scip, timelimit) )
      {
         blocktimelimit = SCIPinfinity(blockprob);
      }
      else
      {
         blocktimelimit = (timelimit - SCIPgetSolvingTime(scip)) * 1.02 + SCIPgetSolvingTime(blockprob);
         blocktimelimit = MIN(SCIPinfinity(blockprob), blocktimelimit); /*lint !e666*/
      }
   }
   else
   {
      if( SCIPisInfinity(scip, timelimit) )
      {
         blocktimelimit = SCIPinfinity(blockprob);
      }
      else
      {
         blocktimelimit = (timelimit - SCIPgetSolvingTime(scip)) * 1.02;
         blocktimelimit = MIN(SCIPinfinity(blockprob), blocktimelimit); /*lint !e666*/
      }
   }

   if( blocktimelimit < 0 )
   {
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsetRealParam(blockprob, "limits/time", blocktimelimit) );

#ifdef SCIP_DEBUG
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "block_%i.lp", blocknum);
   SCIP_CALL( SCIPwriteOrigProblem(blockprob, name, "lp", FALSE) );
#endif

   if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_DANTZIGWOLFE || GCGgetDecompositionMode(gcg) == GCG_DECMODE_ORIGINAL )
   {
      /* try to solve with knapsack solver first */
      if( SCIPgetNConss(blockprob) == 1 )
      {
         SCIP_Real solval;
         GCG_PRICINGSTATUS status;
         int nsolvars;
         SCIP_VAR** solvars = NULL;
         SCIP_Real* solvals = NULL;

         SCIP_CALL( GCGsolverKnapsackSolveKnapsack(TRUE, blockprob, &solval, &status, &solvars, &solvals, &nsolvars) );

         if( solvars != NULL )
         {
            assert(solvals != NULL);

            if( status == GCG_PRICINGSTATUS_OPTIMAL )
            {
               SCIP_SOL* sol = NULL;
               SCIP_Bool stored = FALSE;
               SCIPcreateSol(blockprob, &sol, NULL);
               SCIPsetSolVals(blockprob, sol, nsolvars, solvars, solvals);
               SCIPaddSolFree(blockprob, &sol, &stored);
               assert(stored);
               if( stored )
               {
                  blockprobstatus = SCIP_STATUS_OPTIMAL;
                  (*objvalue) += solval;
               }
            }

            SCIPfreeBufferArray(blockprob, &solvals);
            SCIPfreeBufferArray(blockprob, &solvars);
         }

         if( status == GCG_PRICINGSTATUS_INFEASIBLE )
         {
            blockprobstatus = SCIP_STATUS_INFEASIBLE;
         }
      }

      if( blockprobstatus == SCIP_STATUS_UNKNOWN )
      {
         SCIP_CALL( SCIPsolve(blockprob) );
         blockprobstatus = SCIPgetStatus(blockprob);
      }
   }
   else
   {
      SCIP_BENDERS* benders;
      SCIP_Bool infeasible;

      assert(GCGgetDecompositionMode(gcg) == GCG_DECMODE_BENDERS);

      /* retrieving the Benders' decomposition */
      benders = SCIPfindBenders(GCGgetMasterprob(gcg), "gcg");

      /* since the diagonal blocks are being solved, this indicates that the subproblems are independent. As such, we
       * can declare this in the Benders' decomposition framework. This allows us to call
       * SCIPsolveBendersSubproblem() without setting up the problem
       */
      SCIPbendersSetSubproblemIsIndependent(benders, blocknum, TRUE);

      /* solving the Benders' decomposition subproblem */
      SCIP_CALL( SCIPsolveBendersSubproblem(GCGgetMasterprob(gcg), benders, NULL, blocknum, &infeasible,
            TRUE, NULL) );
   }

   switch( blockprobstatus )
   {
      case SCIP_STATUS_UNBOUNDED:
      case SCIP_STATUS_INFORUNBD:
      case SCIP_STATUS_INFEASIBLE:
         /* no other blocks should be solved. */
         *result = SCIP_CUTOFF;
         break;
      case SCIP_STATUS_BESTSOLLIMIT:
      case SCIP_STATUS_MEMLIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_NODELIMIT:
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_TIMELIMIT:
         /* no other blocks should be solved. */
         *result = SCIP_DIDNOTRUN;
         break;
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_OPTIMAL:
         (*result) = SCIP_SUCCESS;
         if( SCIPgetStage(blockprob) >= SCIP_STAGE_TRANSFORMED )
            (*objvalue) += SCIPgetDualbound(blockprob);
         break;
      default:
         break;
   } /*lint !e788*/

   return SCIP_OKAY;
}

/** frees the block problem */
static
SCIP_RETCODE freeBlockProblem(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP*                 blockprob,          /**< the block problem that will be solved */
   SCIP_RELAXDATA*       relaxdata,          /**< relaxator data structure */
   int                   blocknum            /**< the number of the block, -1 for the original problem */
   )
{
   assert(gcg != NULL);

   if( blockprob == NULL )
      return SCIP_OKAY;

   if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_DANTZIGWOLFE || GCGgetDecompositionMode(gcg) == GCG_DECMODE_ORIGINAL )
   {
      SCIP_CALL( SCIPfreeTransform(blockprob) );
   }
   else
   {
      SCIP_BENDERS* benders;

      assert(GCGgetDecompositionMode(gcg) == GCG_DECMODE_BENDERS);

      /* retrieving the Benders' decomposition */
      benders = SCIPfindBenders(GCGgetMasterprob(gcg), "gcg");

      /* freeing the Benders' decomposition subproblems */
      SCIP_CALL( SCIPfreeBendersSubproblem(GCGgetMasterprob(gcg), benders, blocknum) );
   }

   return SCIP_OKAY;
}

/** solves the blocks diagonal and individually */
static
SCIP_RETCODE solveDiagonalBlocks(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< relaxator data structure */
   SCIP_RESULT*          result,             /**< result pointer to indicate success or failure */
   SCIP_Real*            lowerbound          /**< lower bound pointer to return the lower bound */
   )
{
   SCIP* scip;
   int i;
   SCIP_Real timelimit;
   SCIP_Real objvalue;
   SCIP_SOL *newsol;
   SCIP_Bool isfeasible;
   SCIP_RESULT solveresult;

   scip = GCGgetOrigprob(gcg);

   /* set objective of pricing problems to original objective */
   if( GCGgetDecompositionMode(gcg) != GCG_DECMODE_ORIGINAL )
   {
      SCIP_CALL( setPricingObjsOriginal(gcg, relaxdata->pricingprobs, relaxdata->npricingprobs) );
   }

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   objvalue = 0.0;

   if( GCGgetDecompositionMode(gcg) != GCG_DECMODE_ORIGINAL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Block diagonal structure detected, solving blocks individually.\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "There is an objective function offset of %f.\n", SCIPgetTransObjoffset(scip));
   }

   /* if the original problem is solved directly, then we call  solveBlockProblem with the master problem */
   if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_ORIGINAL )
   {
      SCIP_CALL( solveBlockProblem(gcg, GCGgetMasterprob(gcg), relaxdata, timelimit, -1, &solveresult, &objvalue) );

      if( solveresult == SCIP_CUTOFF || solveresult == SCIP_DIDNOTRUN )
      {
         (*result) = solveresult;
         return SCIP_OKAY;
      }
   }
   else
   {
      /* solve pricing problems one after the other */
      for( i = 0; i < relaxdata->npricingprobs; ++i )
      {
         SCIP_CALL( solveBlockProblem(gcg, relaxdata->pricingprobs[i], relaxdata, timelimit, i, &solveresult, &objvalue) );

         if( solveresult == SCIP_CUTOFF || solveresult == SCIP_DIDNOTRUN )
         {
            (*result) = solveresult;
            return SCIP_OKAY;
         }
      }
   }

   /* get solution and glue it together */

   if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_ORIGINAL )
   {
      SCIP_CALL( GCGtransformMastersolToOrigsol(gcg, SCIPgetBestSol(GCGgetMasterprob(gcg)), &newsol, TRUE, NULL) );
   }
   else
   {
      SCIP_CALL( combineSolutions(gcg, &newsol, relaxdata->pricingprobs, relaxdata->npricingprobs) );
   }

   /* update lower bound pointer and add solution such that this node will be cut off automatically */
   if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE )
      *lowerbound = -objvalue;
   else
      *lowerbound = objvalue;

   SCIP_CALL( SCIPcheckSol(scip, newsol, TRUE, TRUE, TRUE, TRUE, TRUE, &isfeasible) );
   assert(isfeasible);

   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &isfeasible) );

   /** @todo maybe add a constraint here to indicate that it has been decomposed */

   /* if the original problem is solved directly, then we call freeBlockProblem with the master problem */
   if( GCGgetDecompositionMode(gcg) != GCG_DECMODE_ORIGINAL )
   {
      /* solve pricing problems one after the other */
      for( i = 0; i < relaxdata->npricingprobs; ++i )
      {
         SCIP_CALL( freeBlockProblem(gcg, relaxdata->pricingprobs[i], relaxdata, i) );
      }
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;

}

GCG_DECOMP* GCGgetStructDecomp(
   GCG*                  gcg
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->decomp;
}

/** sets the structure information */
static
SCIP_RETCODE GCGsetStructDecomp(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);
   assert(decomp != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if( relaxdata->decomp != NULL )
      SCIP_CALL( GCGdecompFree(gcg, &relaxdata->decomp ) );

   relaxdata->decomp = decomp;

   return SCIP_OKAY;
}

/** transforms the master problem **/
static
SCIP_RETCODE transformMaster(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAX*           relax               /**< relaxator data structure */
   )
{
   SCIP* scip;
   SCIP* masterprob;
   SCIP_VAR** vars;
   SCIP_CONS** oldconss;
   SCIP_RELAXDATA* relaxdata;
   int i;
   int nvars;

   assert(gcg != NULL);
   assert(relax != NULL);

   scip = GCGgetOrigprob(gcg);
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   masterprob = GCGgetMasterprob(relaxdata->gcg);
   assert(masterprob != NULL);
   SCIP_CALL( SCIPtransformProb(masterprob) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &oldconss, relaxdata->masterconss, relaxdata->nmasterconss) );

   /* transform the master constraints */
   SCIP_CALL( SCIPtransformConss(masterprob, relaxdata->nmasterconss,
                                 relaxdata->masterconss, relaxdata->masterconss) );
   for( i = 0; i < relaxdata->nmasterconss; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(masterprob, &(oldconss[i])) );
   }
   SCIPfreeBufferArray(scip, &oldconss);

   /* transform the convexity constraints */
   for( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      if( relaxdata->convconss[i] != NULL )
      {
         SCIP_CONS* oldcons = relaxdata->convconss[i];
         SCIP_CALL( SCIPreleaseCons(masterprob, &oldcons) );
         SCIP_CALL( SCIPtransformCons(masterprob, relaxdata->convconss[i], &(relaxdata->convconss[i])) );
      }
   }

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   /* transform the linking variable constraints */
   for( i = 0; i < nvars; ++i )
   {
      assert(GCGvarIsOriginal(vars[i]));

      if( GCGoriginalVarIsLinking(vars[i]) )
      {
         int j;
         SCIP_CONS** linkconss;
         linkconss = GCGlinkingVarGetLinkingConss(vars[i]);
         /* the linking constraints could be NULL if the Benders' decomposition is used. */
         if( linkconss != NULL )
         {
            for( j = 0; j < relaxdata->npricingprobs; ++j )
            {
               if( linkconss[j] != NULL )
               {
                  SCIP_CONS* tempcons;
                  SCIP_CALL( SCIPtransformCons(masterprob, linkconss[j], &(tempcons)) );
                  GCGlinkingVarSetLinkingCons(vars[i], tempcons, j);
               }
            }
         }
      }
   }
   for( i = 0; i < relaxdata->nvarlinkconss; ++i )
   {
      SCIP_CONS* transcons;

      SCIP_CALL( SCIPgetTransformedCons(masterprob, relaxdata->varlinkconss[i], &transcons) );
      assert(transcons != NULL);

      SCIP_CALL( SCIPreleaseCons(masterprob, &relaxdata->varlinkconss[i]) );
      relaxdata->varlinkconss[i] = transcons;
   }
   return SCIP_OKAY;
}

/** initializes and transforms relaxator data */
static
SCIP_RETCODE initRelaxator(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAX*           relax               /**< relaxator data structure */
   )
{
   SCIP* scip;
   SCIP_RELAXDATA* relaxdata;
   int permutationseed;
   int oxfordcomma;

   assert(gcg != NULL);
   assert(relax != NULL);

   scip = GCGgetOrigprob(gcg);
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* when the original problem should be solved directly, then a decomposition must be made with zero blocks */
   if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_ORIGINAL )
   {
      GCG_DECOMP* decomp;
      SCIP_RETCODE retcode;

      assert(relaxdata->decomp == NULL);

      retcode = GCGcreateBasicDecomp(gcg, &decomp, TRUE);
      assert(retcode == SCIP_OKAY);
      if( retcode != SCIP_OKAY )
      {
         SCIPerrorMessage("Could not add decomp to cons_decomp!\n");
         return SCIP_ERROR;
      }

      assert(decomp != NULL );

      GCGsetStructDecomp(gcg, decomp);
   }

   if( relaxdata->decomp == NULL )
   {
      relaxdata->decomp = GCGgetBestDecomp(gcg, TRUE);
      if( relaxdata->decomp == NULL )
      {
         PARTIALDECOMP_C* partialdec;
         SCIPwarningMessage(scip, "No complete decomposition available. Creating basic decomposition.\n");
         SCIP_CALL( GCGconshdlrDecompAddBasicPartialdec(gcg, TRUE, &partialdec) );
         SCIP_CALL( GCGconshdlrDecompSelectPartialdec(partialdec, TRUE) );

         relaxdata->decomp = GCGgetBestDecomp(gcg, FALSE);
         assert( relaxdata->decomp != NULL );
      }
   }

   oxfordcomma = 0;
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Chosen structure has %d blocks", GCGdecompGetNBlocks(relaxdata->decomp));
   /* every master-only variable internally also counts as linking, but should not be reported as linking variable */
   if ( GCGdecompGetNLinkingvars(relaxdata->decomp) - GCGdecompGetNMastervars(relaxdata->decomp) > 0)
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ", %d linking variables", GCGdecompGetNLinkingvars(relaxdata->decomp) - GCGdecompGetNMastervars(relaxdata->decomp));
      ++oxfordcomma;
   }
   if ( GCGdecompGetNMastervars(relaxdata->decomp) > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ", %d master-only (static) variables", GCGdecompGetNMastervars(relaxdata->decomp));
      ++oxfordcomma;
   }
   if ( oxfordcomma > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ",");
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " and %d linking constraints.\n", GCGdecompGetNLinkingconss(relaxdata->decomp));
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "This decomposition has a maxwhite score of %f.\n", GCGdecompGetMaxwhiteScore(relaxdata->decomp));

   /* permute the decomposition if the permutation seed is set */
   SCIP_CALL( SCIPgetIntParam(scip, "randomization/permutationseed", &permutationseed) );

   if( permutationseed > 0 )
   {
      SCIP_RANDNUMGEN* randnumgen;

      SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, (unsigned int) permutationseed, TRUE) );
      SCIP_CALL( GCGpermuteDecomp(gcg, relaxdata->decomp, randnumgen) );
      SCIPfreeRandom(scip, &randnumgen);
   }

   if( relaxdata->discretization && (SCIPgetNContVars(scip) > 0) )
   {
      if( relaxdata->mipdiscretization )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Warning: Discretization with continuous variables is only an experimental feature.\n");
      }
      else
      {
         SCIP_CALL( SCIPsetBoolParam(scip, "relaxing/gcg/discretization", FALSE) );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Warning: Discretization with continuous variables is disabled by parameter relaxing/gcg/mipdiscretization.\n");
      }
   }

   SCIP_CALL( createMaster(gcg, relaxdata) );

#ifdef _OPENMP
   if( relaxdata->mode == GCG_DECMODE_DANTZIGWOLFE && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_NORMAL )
   {
      int ompmaxthreads = omp_get_max_threads();
      int nthreads = GCGpricerGetMaxNThreads(gcg);
      if( nthreads > 0 )
         nthreads = MIN(nthreads, GCGgetNRelPricingprobs(gcg));
      else
         nthreads = MIN(ompmaxthreads, GCGgetNRelPricingprobs(gcg));
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Using up to %d (of %d) thread(s) to solve the pricing problems.\n", nthreads, ompmaxthreads);
   }
#endif

   /* for Benders' decomposition, the Benders' plugin must be activated */
   if( relaxdata->mode == GCG_DECMODE_BENDERS )
   {
      SCIP_CALL( SCIPactivateBenders(GCGgetMasterprob(gcg), SCIPfindBenders(GCGgetMasterprob(gcg), "gcg"),
            relaxdata->npricingprobs) );
   }

   relaxdata->lastsolvednodenr = -1;

   /* set objective limit in master problem if objective limit in original problem is finite */
   if( !SCIPisInfinity(scip, (int) SCIPgetObjsense(scip) * SCIPgetObjlimit(scip)) )
   {
      SCIP_CALL( SCIPsetObjlimit(GCGgetMasterprob(relaxdata->gcg), (int) SCIPgetObjsense(scip) * SCIPgetObjlimit(scip)) );
   }

   relaxdata->relaxisinitialized = TRUE;

   return SCIP_OKAY;
}

#ifdef _OPENMP
/** initializes all OpenMP locks */
static
void initLocks(
   GCG_LOCKS*            locks               /**< OpenMP locks */
   )
{
   assert(locks != NULL);
   GCG_INIT_LOCK(&locks->memorylock);
   GCG_INIT_LOCK(&locks->pricinglock);
   GCG_INIT_LOCK(&locks->pricinglimitslock);
   GCG_INIT_LOCK(&locks->pricestorelock);
   GCG_INIT_LOCK(&locks->printlock);
}

/** destroys all OpenMP locks */
static
void destroyLocks(
   GCG_LOCKS*            locks               /**< OpenMP locks */
   )
{
   assert(locks != NULL);
   GCG_DESTROY_LOCK(&locks->memorylock);
   GCG_DESTROY_LOCK(&locks->pricinglock);
   GCG_DESTROY_LOCK(&locks->pricinglimitslock);
   GCG_DESTROY_LOCK(&locks->pricestorelock);
   GCG_DESTROY_LOCK(&locks->printlock);
}
#endif

/** initializes relaxator data */
static
SCIP_RETCODE initRelaxdata(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_RELAXDATA*       relaxdata           /**< relaxdata data structure */
   )
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(relaxdata != NULL);

   relaxdata->gcg = gcg;
   relaxdata->decomp = NULL;

   relaxdata->nbranchrules = 0;
   relaxdata->branchrules = NULL;
   relaxdata->paramsvisu = NULL;

   relaxdata->blockrepresentative = NULL;
   relaxdata->convconss = NULL;
   relaxdata->lastsolvednodenr = 0;

   relaxdata->origmasterconss = NULL;
   relaxdata->masterconss = NULL;
   relaxdata->nmasterconss = 0;

   relaxdata->npricingprobs = -1;
   relaxdata->maxpricingprobs = 0;
   relaxdata->pricingprobs = NULL;
   relaxdata->nrelpricingprobs = 0;
   relaxdata->currentorigsol = NULL;
   relaxdata->storedorigsol = NULL;
   relaxdata->origsolfeasible = FALSE;
   relaxdata->storedfeasibility = FALSE;
   relaxdata->nblocksidentical = NULL;

   relaxdata->lastmastersol = NULL;
   relaxdata->lastmasterlpiters = 0;
   relaxdata->lastmasternode = -1;
   relaxdata->markedmasterconss = NULL;
   relaxdata->maxmarkedmasterconss = 0;
   relaxdata->masterinprobing = FALSE;
   relaxdata->probingheur = NULL;

   relaxdata->ntransvars = 0;
   relaxdata->nlinkingvars = 0;
   relaxdata->nvarlinkconss = 0;
   relaxdata->varlinkconss = NULL;
   relaxdata->varlinkconsblock = NULL;
   relaxdata->pricingprobsmemused = 0.0;

   relaxdata->relaxisinitialized = FALSE;
   relaxdata->simplexiters = 0;
   relaxdata->rootnodetime = NULL;

   relaxdata->limitsettingsstashed = FALSE;

   relaxdata->activebranchrules = NULL;
   relaxdata->activebranchdata = NULL;
   relaxdata->activebranchextendedmasterconss = NULL;
   relaxdata->nactivebranchextendedmasterconss = 0;
   relaxdata->maxactivebranchextendedmasterconss = 0;

   SCIP_CALL( GCGcreateParamsVisu(gcg, &(relaxdata->paramsvisu)) );
   assert(relaxdata->paramsvisu != NULL);

#ifdef _OPENMP
   SCIP_CALL( SCIPallocBlockMemory(origprob, &relaxdata->locks) );
   initLocks(relaxdata->locks);
#endif

return SCIP_OKAY;
}

/** resets relaxator data */
static
void resetRelaxdata(
   SCIP_RELAXDATA*       relaxdata           /**< relaxdata data structure */
   )
{
   assert(relaxdata != NULL);

   assert(relaxdata->decomp == NULL);

   relaxdata->lastsolvednodenr = 0;

   assert(relaxdata->origmasterconss == NULL);
   assert(relaxdata->masterconss == NULL);
   relaxdata->nmasterconss = 0;

   relaxdata->npricingprobs = -1;
   relaxdata->nrelpricingprobs = 0;
   assert(relaxdata->currentorigsol == NULL);
   assert(relaxdata->storedorigsol == NULL);
   relaxdata->origsolfeasible = FALSE;
   relaxdata->storedfeasibility = FALSE;

   relaxdata->lastmastersol = NULL;
   relaxdata->lastmasterlpiters = 0;
   relaxdata->lastmasternode = -1;
   assert(relaxdata->markedmasterconss == NULL);
   assert(relaxdata->maxmarkedmasterconss == 0);
   assert(relaxdata->masterinprobing == FALSE);
   assert(relaxdata->probingheur == NULL);

   relaxdata->ntransvars = 0;
   relaxdata->nlinkingvars = 0;
   relaxdata->nvarlinkconss = 0;
   assert(relaxdata->varlinkconss == NULL);
   assert(relaxdata->varlinkconsblock == NULL);
   relaxdata->pricingprobsmemused = 0.0;

   assert(relaxdata->relaxisinitialized == FALSE);
   relaxdata->simplexiters = 0;
   assert(relaxdata->rootnodetime == NULL);
}

/*
 * Callback methods of relaxator
 */

/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreeGcg)
{
   SCIP_RELAXDATA* relaxdata;
   int i;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* free pricing problems */
   if( relaxdata->pricingprobs != NULL )
   {
      for( i = relaxdata->maxpricingprobs - 1; i >= 0 ; i-- )
      {
         SCIP_CALL( SCIPfree(&(relaxdata->pricingprobs[i])) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->pricingprobs), relaxdata->maxpricingprobs);
      SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->blockrepresentative), relaxdata->maxpricingprobs);
      SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->nblocksidentical), relaxdata->maxpricingprobs);
      SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->convconss), relaxdata->npricingprobs);
   }

   /* free visualization parameters */
   if( relaxdata->paramsvisu != NULL )
   {
      GCGVisuFreeParams(relaxdata->gcg, relaxdata->paramsvisu);
   }

   /* free master problem */
   if( relaxdata->gcg->dwmasterprob != NULL )
   {
      SCIP_CALL( SCIPfree(&(relaxdata->gcg->dwmasterprob)) );
   }

   /* free the benders master problem */
   if( relaxdata->gcg->bendersmasterprob != NULL )
   {
      SCIP_CALL( SCIPfree(&(relaxdata->gcg->bendersmasterprob)) );
   }

   /* free used decomposition */
   if( relaxdata->decomp != NULL )
   {
      SCIP_CALL( GCGdecompFree(relaxdata->gcg, &relaxdata->decomp) );
   }

#ifdef _OPENMP
   /* free locks struct */
   if( relaxdata->locks != NULL )
   {
      destroyLocks(relaxdata->locks);
      SCIPfreeBlockMemory(scip, &relaxdata->locks);
   }
#endif

   relaxdata->gcg->masterprob = NULL;
   relaxdata->gcg->dwmasterprob = NULL;
   relaxdata->gcg->bendersmasterprob = NULL;
   relaxdata->gcg->relax = NULL;

   SCIPfreeMemory(scip, &relaxdata);
   return SCIP_OKAY;
}

/** deinitialization method of relaxator (called before transformed problem is freed) */

static
SCIP_DECL_RELAXEXIT(relaxExitGcg)
{
   SCIP_RELAXDATA* relaxdata;
   int i;

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* free pricing problems */
   for( i = relaxdata->npricingprobs - 1; i >= 0 ; i-- )
   {
      SCIP_CALL( SCIPfreeProb(relaxdata->pricingprobs[i]) );
   }

   if( relaxdata->decomp != NULL )
   {
      SCIP_CALL( GCGdecompFree(relaxdata->gcg, &relaxdata->decomp) );
      relaxdata->decomp = NULL;
   }

   /* free array for branchrules*/
   if( relaxdata->nbranchrules > 0 )
   {
      for( i = 0; i < relaxdata->nbranchrules; i++ )
      {
         SCIPfreeMemory(scip, &(relaxdata->branchrules[i]));
      }
      SCIPfreeMemoryArray(scip, &(relaxdata->branchrules));
   }


   relaxdata->nbranchrules = 0;
   relaxdata->relaxisinitialized = FALSE;
   relaxdata->limitsettingsstashed = FALSE;

   return SCIP_OKAY;
}


/** initialize the relaxator and master problem for solving the original problem by Dantzig-Wolfe reformulation and
 * Benders' decomposition
 */
static
SCIP_RETCODE initializeMasterProblemSolve(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP_RELAX*           relax               /**< the relaxator */
)
{
   SCIP* scip;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);
   assert(relax != NULL);

   scip = GCGgetOrigprob(gcg);
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if( !SCIPisTransformed(GCGgetMasterprob(relaxdata->gcg)) )
   {
      /* set integral objective status in the extended problem, if possible */
      if( SCIPisObjIntegral(scip) && relaxdata->discretization && SCIPgetNContVars(scip) == 0
          && relaxdata->mode == GCG_DECMODE_DANTZIGWOLFE )
      {
         SCIP_CALL( SCIPsetObjIntegral(GCGgetMasterprob(relaxdata->gcg)) );
      }
      SCIP_CALL( transformMaster(gcg, relax) );
      /* transform the decomposition */
      // SCIP_CALL( GCGdecompTransform(scip, relaxdata->decomp) );
      SCIP_CALL( GCGconsOrigbranchAddRootCons(gcg) );
      assert(relaxdata->decomp != NULL);
   }

   return SCIP_OKAY;
}


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
static
SCIP_DECL_RELAXINITSOL(relaxInitsolGcg)
{
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(GCGgetMasterprob(relaxdata->gcg) != NULL);

   resetRelaxdata(relaxdata);
   SCIP_CALL( SCIPcreateClock(scip, &(relaxdata->rootnodetime)) );

   /* set active masterprob */
   switch( GCGgetDecompositionMode(relaxdata->gcg) )
   {
      case GCG_DECMODE_ORIGINAL:
      case GCG_DECMODE_BENDERS:
      {
         relaxdata->gcg->masterprob = relaxdata->gcg->bendersmasterprob;
         break;
      }
      case GCG_DECMODE_DANTZIGWOLFE:
      {
         relaxdata->gcg->masterprob = relaxdata->gcg->dwmasterprob;
         break;
      }
      default:
      {
         SCIPerrorMessage("Unknown decomposition mode.");
         return SCIP_ERROR;
      }
   }

   /* alternative verbosity levels are used for the Benders' decomposition and original mode compared to the Dantzig-Wolfe
    * decomposition mode.
    */
   if( GCGgetDecompositionMode(relaxdata->gcg) == GCG_DECMODE_BENDERS || GCGgetDecompositionMode(relaxdata->gcg) == GCG_DECMODE_ORIGINAL )
   {
      /* first getting the verbosity level for the original problem before setting it to none. While the verbosity level
       * was collected previously, the user may have changed this in the mean time.
       */
      SCIP_CALL( SCIPgetIntParam(scip, "display/verblevel", &relaxdata->origverblevel) );

      /* deactivating display columns */
      SCIP_CALL( SCIPsetIntParam(scip, "display/sumlpiterations/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/lpiterations/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/degeneracy/active", 0) );

      /* setting the total node limit to 1 for the original SCIP instance. This is because Benders' decomposition solves
       * the MIP within the relaxator of the root node. So no branching in the original problem is required.
       */
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/totalnodes", 1LL) );
   }

   /* fixing the GCG mode parameter. This ensure that the user does not change this during the solution process. If the
    * mode parameter were to change, the behaviour is unknown.
    */
   SCIP_CALL( SCIPfixParam(scip, "relaxing/gcg/mode") );

   /* Informing the user of the decomposition technique that is being used to solve the original problem */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "\n");
   if( relaxdata->mode == GCG_DECMODE_DANTZIGWOLFE )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "A Dantzig-Wolfe reformulation is applied to solve the original problem.\n");
   }
   else if( relaxdata->mode == GCG_DECMODE_BENDERS )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "A Benders' decomposition is applied to solve the original problem.\n");
   }
   else if( relaxdata->mode == GCG_DECMODE_ORIGINAL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "No reformulation will be performed. Solving the original model.\n");
   }

   if( !SCIPisStopped(scip) )
      SCIP_CALL( initRelaxator(relaxdata->gcg, relax) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolGcg)
{
   SCIP_RELAXDATA* relaxdata;
   int i;

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->markedmasterconss), relaxdata->maxmarkedmasterconss);
   relaxdata->maxmarkedmasterconss = 0;

   /* free arrays for constraints */
   for( i = 0; i < relaxdata->nmasterconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &relaxdata->origmasterconss[i]) );
      SCIP_CALL( SCIPreleaseCons(GCGgetMasterprob(relaxdata->gcg), &relaxdata->masterconss[i]) );
   }
   for( i = 0; i < relaxdata->npricingprobs; i++ )
   {
      if( relaxdata->convconss[i] != NULL )
         SCIP_CALL( SCIPreleaseCons(GCGgetMasterprob(relaxdata->gcg), &relaxdata->convconss[i]) );
   }
   for( i = 0; i < relaxdata->nvarlinkconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(GCGgetMasterprob(relaxdata->gcg), &relaxdata->varlinkconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->varlinkconss), relaxdata->nvarlinkconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->varlinkconsblock), relaxdata->nvarlinkconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->origmasterconss), relaxdata->maxmasterconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->masterconss), relaxdata->maxmasterconss);

   /* free master problem */
   if( GCGgetMasterprob(relaxdata->gcg) != NULL )
   {
      SCIP_CALL( SCIPfreeProb(GCGgetMasterprob(relaxdata->gcg)) );
   }

   /* free solutions */
   if( relaxdata->currentorigsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &relaxdata->currentorigsol) );
      relaxdata->currentorigsol = NULL;
   }
   if( relaxdata->storedorigsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &relaxdata->storedorigsol) );
      relaxdata->storedorigsol = NULL;
   }

   if( relaxdata->decomp != NULL )
   {
      SCIP_CALL( GCGdecompFree(relaxdata->gcg, &relaxdata->decomp) );
      relaxdata->decomp = NULL;
   }

   SCIP_CALL( GCGfreeOrigVarsData(relaxdata->gcg) );

   /* free root node clock */
   if( relaxdata->rootnodetime != NULL )
   {
      SCIP_CALL( SCIPfreeClock(scip, &(relaxdata->rootnodetime)) );
      relaxdata->rootnodetime = NULL;
   }

   if( relaxdata->maxactivebranchextendedmasterconss > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->activebranchrules), relaxdata->maxactivebranchextendedmasterconss);
      SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->activebranchdata), relaxdata->maxactivebranchextendedmasterconss);
      SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->activebranchextendedmasterconss), relaxdata->maxactivebranchextendedmasterconss);
      relaxdata->maxactivebranchextendedmasterconss = 0;
      relaxdata->nactivebranchextendedmasterconss = 0;
   }

   relaxdata->relaxisinitialized = FALSE;

   return SCIP_OKAY;
}


/** sets (time and gap) limits in the master problem based on limits of the original problem */
static
SCIP_RETCODE setMasterLimits(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP*                 masterprob,         /**< the master problem SCIP instance */
   SCIP_Real             origtimelimit,      /**< time limit fo the original problem */
   SCIP_Real             origgaplimit        /**< gap limit of the original problem */
   )
{
   SCIP* scip;
   SCIP_Real mastertimelimit;
   
   scip = GCGgetOrigprob(gcg);
   mastertimelimit = SCIPinfinity(scip);

   if( masterprob == NULL )
      masterprob = GCGgetMasterprob(gcg);

   if( !SCIPisInfinity(scip, origtimelimit) )
   {
      /* give the master 0.5 seconds more time than the original scip has left */
      mastertimelimit = (origtimelimit - SCIPgetSolvingTime(scip)) + 0.5 + SCIPgetSolvingTime(masterprob);
      assert(origtimelimit - SCIPgetSolvingTime(scip) > 0 || SCIPisStopped(scip));

      SCIPdebugMessage("  time limit for master: %f, left: %f, left for original problem: %f\n",
            mastertimelimit,
            mastertimelimit - SCIPgetSolvingTime(masterprob),
            origtimelimit - SCIPgetSolvingTime(scip));
   }
   SCIP_CALL( SCIPsetRealParam(masterprob, "limits/time", mastertimelimit) );

   /* set gap limit */
   SCIP_CALL( SCIPsetRealParam(masterprob, "limits/gap", origgaplimit) );
   return SCIP_OKAY;
}


/** method to solve the master problem that is used by Dantzig-Wolfe and Benders' decomposition */
static
SCIP_RETCODE solveMasterProblem(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP*                 masterprob,         /**< the master problem SCIP instance */
   SCIP_RELAXDATA*       relaxdata,          /**< the relaxator data */
   SCIP_Longint          nodelimit,          /**< the number of nodes the will be solved in this master problem */
   SCIP_Real*            lowerbound,         /**< the lowerbound computed by the relaxator for the current node */
   SCIP_RESULT*          result              /**< the result of the relaxation call */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   assert(masterprob != NULL);
   assert(relaxdata != NULL);

   scip = GCGgetOrigprob(gcg);

   /* update the number of the last solved node */
   relaxdata->lastsolvednodenr = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   /* increase the node limit for the master problem by 1 */
   SCIP_CALL( SCIPsetLongintParam(masterprob, "limits/nodes", nodelimit) );

   /* loop to solve the master problem, this is a workaround and does not fix any problem */
   do
   {
      SCIP_Real timelimit;
      SCIP_Real memorylimit;
      SCIP_Real gaplimit;

      /* set memorylimit for master */
      SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
      if( !SCIPisInfinity(scip, memorylimit) )
         memorylimit -= SCIPgetMemUsed(scip)/1048576.0;

      SCIP_CALL( SCIPsetRealParam(masterprob, "limits/memory", memorylimit) );

      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      SCIP_CALL( SCIPgetRealParam(scip, "limits/gap", &gaplimit) );

      SCIP_CALL( setMasterLimits(gcg, masterprob, timelimit, gaplimit) );

      /* if we have a blockdetection, see whether the node is block diagonal. Additionally, the solveDiagonalBlocks can
       * be called when the original problem is solved directly.
       */
      if( GCGdecompGetType(relaxdata->decomp) == GCG_DECTYPE_DIAGONAL || relaxdata->mode == GCG_DECMODE_ORIGINAL )
      {
         SCIP_CALL( solveDiagonalBlocks(gcg, relaxdata, result, lowerbound) );
         if( *result == SCIP_SUCCESS || *result == SCIP_CUTOFF )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
      }
      /* We are solving the masterproblem regularly */
      else
      {
         SCIP_CALL( SCIPsolve(masterprob) );
      }


      if( SCIPgetStatus(masterprob) != SCIP_STATUS_TIMELIMIT )
      {
         break;
      }

      if( !SCIPisInfinity(scip, timelimit) && !SCIPisStopped(scip) )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "time for master problem was too short, extending time.\n");
   }
   while( !SCIPisStopped(scip) );

   if( SCIPgetStatus(masterprob) == SCIP_STATUS_TIMELIMIT && SCIPisStopped(scip) )
   {
      if( SCIPgetCurrentNode(masterprob) == NULL || !GCGmasterIsCurrentSolValid(gcg) || !SCIPisGT(scip, SCIPgetLocalDualbound(masterprob), SCIPgetLocalLowerbound(scip)) )
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }
   }

   /* set the lower bound pointer */
   if( SCIPgetStage(masterprob) == SCIP_STAGE_SOLVING && GCGmasterIsCurrentSolValid(gcg) )
   {
      *lowerbound = SCIPgetLocalDualbound(masterprob);
      if( SCIPisInfinity(scip, *lowerbound) )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }
   else
   {
      SCIPdebugMessage("  stage: %d\n", SCIPgetStage(masterprob));
      assert(SCIPgetStatus(masterprob) == SCIP_STATUS_TIMELIMIT || SCIPgetBestSol(masterprob) != NULL || SCIPgetStatus(masterprob) == SCIP_STATUS_INFEASIBLE || SCIPgetStatus(masterprob) == SCIP_STATUS_UNKNOWN);
      if( SCIPgetStatus(masterprob) == SCIP_STATUS_OPTIMAL && GCGmasterIsCurrentSolValid(gcg) )
      {
         *lowerbound = SCIPgetSolOrigObj(masterprob, SCIPgetBestSol(masterprob));
         assert(!SCIPisInfinity(scip, *lowerbound));
      }
      else if( SCIPgetStatus(masterprob) == SCIP_STATUS_INFEASIBLE || SCIPgetStatus(masterprob) == SCIP_STATUS_TIMELIMIT || !GCGmasterIsCurrentSolValid(gcg) )
      {
         SCIP_Real tilim;
         SCIP_CALL( SCIPgetRealParam(masterprob, "limits/time", &tilim) );
         if( tilim-SCIPgetSolvingTime(masterprob) < 0 )
         {
            *result = SCIP_DIDNOTRUN;
            return SCIP_OKAY;
         }
         *lowerbound = SCIPinfinity(scip);
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else if( SCIPgetStatus(masterprob) == SCIP_STATUS_UNKNOWN )
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }
      else
      {
         SCIPwarningMessage(scip, "Stage <%d> is not handled!\n", SCIPgetStage(masterprob));
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }
   }

   SCIPdebugMessage("  update lower bound (value = %g).\n", *lowerbound);

   /* NOTE: All other points when result is set, the function is exited immediately. Ensure that this is checked for
    * future changes to this function
    */
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** execution method of the relaxator for Dantzig-Wolfe reformulation */
static
SCIP_RETCODE relaxExecGcgDantzigWolfe(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP_RELAX*           relax,              /**< the relaxator */
   SCIP_Real*            lowerbound,         /**< the lowerbound computed by the relaxator for the current node */
   SCIP_RESULT*          result              /**< the result of the relaxation call */
   )
{
   SCIP* scip;
   SCIP* masterprob;
   SCIP_RELAXDATA* relaxdata;
   SCIP_Longint oldnnodes;
   SCIP_Longint nodelimit;
   SCIP_Bool stored;

   assert(gcg != NULL);
   assert(relax != NULL);
   assert(result != NULL);
   assert(GCGgetDecompositionMode(gcg) == GCG_DECMODE_DANTZIGWOLFE);

   scip = GCGgetOrigprob(gcg);
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   *result = SCIP_DIDNOTRUN;

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);

   assert(relaxdata->gcg->masterprob == relaxdata->gcg->dwmasterprob);

   /* solve the next node in the master problem */
   SCIPdebugMessage("Solving node %"SCIP_LONGINT_FORMAT"'s relaxation.\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   /* only solve the relaxation if it was not yet solved at the current node */
   if( SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) != relaxdata->lastsolvednodenr )
   {
      SCIP_CONS* activeorigcons;

      /* start root node time clock */
      if( SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) )
      {
         SCIP_CALL( SCIPstartClock(scip, relaxdata->rootnodetime) );
         SCIPdebugMessage("  root node time clock started.\n");
      }

      /* increase the node limit for the master problem by 1 */
      SCIP_CALL( SCIPgetLongintParam(masterprob, "limits/nodes", &oldnnodes) );

      nodelimit = (SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) ? 1 : oldnnodes + 1);
      /* solving the master problem */
      SCIP_CALL( solveMasterProblem(gcg, masterprob, relaxdata, nodelimit, lowerbound, result) );
      assert(*result == SCIP_CUTOFF || !SCIPisInfinity(scip, *lowerbound));

      if( relaxdata->currentorigsol != NULL )
      {
         SCIP_CALL( SCIPtrySol(scip, relaxdata->currentorigsol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );
      }

      /* if a new primal solution was found in the master problem, transfer it to the original problem */
      if( SCIPgetBestSol(GCGgetMasterprob(gcg)) != NULL && relaxdata->lastmastersol != SCIPgetBestSol(GCGgetMasterprob(gcg)) && GCGmasterIsCurrentSolValid(gcg) )
      {
         SCIP_SOL* newsol;

         relaxdata->lastmastersol = SCIPgetBestSol(GCGgetMasterprob(gcg));

         SCIP_CALL( GCGtransformMastersolToOrigsol(gcg, relaxdata->lastmastersol, &newsol, TRUE, NULL) );
   #ifdef SCIP_DEBUG
         SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, TRUE, &stored) );
   #else
         SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );
   #endif
   #ifndef NDEBUG
         /* only check failed solution if best master solution is valid */
         if( !stored && GCGmasterIsBestsolValid(gcg) )
         {
            SCIP_CALL( SCIPcheckSolOrig(scip, newsol, &stored, TRUE, TRUE) );
         }
   #endif
         /** @bug The solution doesn't have to be accepted, numerics might bite us, so the transformation might fail.
          *  A remedy could be: Round the values or propagate changes or call a heuristic to fix it.
          *  SCIP rejects a solution if it is equal to a known one
          */
         SCIP_CALL( SCIPfreeSol(scip, &newsol) );

         if( stored )
            SCIPdebugMessage("  updated current best primal feasible solution.\n");
      }

      activeorigcons = GCGconsOrigbranchGetActiveCons(gcg);
      if( GCGconsOrigbranchGetBranchrule(activeorigcons) != NULL )
      {
         SCIP_CALL( GCGrelaxBranchMasterSolved(gcg, GCGconsOrigbranchGetBranchrule(activeorigcons ),
               GCGconsOrigbranchGetBranchdata(activeorigcons), *lowerbound) );
      }

      /* stop root node clock */
      if( SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip))
      {
         SCIP_CALL( SCIPstopClock(scip, relaxdata->rootnodetime) );
         SCIPdebugMessage("  root node time clock stopped at %6.2fs.\n", SCIPgetClockTime(scip, relaxdata->rootnodetime));
      }
   }
   else
   {
      SCIPdebugMessage("Problem has been already solved at this node\n");
   }

   if( SCIPgetStatus(masterprob) == SCIP_STATUS_OPTIMAL )
      *result = SCIP_CUTOFF;

   return SCIP_OKAY;
}


/** method to solve the master problem for Benders' decomposition and when solving the original problem directly. */
static
SCIP_RETCODE solveMasterProblemAndEvaluate(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP_RELAX*           relax,              /**< the relaxator */
   SCIP_Real*            lowerbound,         /**< the lowerbound computed by the relaxator for the current node */
   SCIP_RESULT*          result              /**< the result of the relaxation call */
   )
{
   SCIP* scip;
   SCIP* masterprob;
   SCIP_RELAXDATA* relaxdata;
   SCIP_Longint nodelimit;

   assert(gcg != NULL);
   assert(relax != NULL);
   assert(result != NULL);

   scip = GCGgetOrigprob(gcg);
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   *result = SCIP_DIDNOTRUN;

   masterprob = GCGgetMasterprob(relaxdata->gcg);
   assert(masterprob != NULL);

   assert(relaxdata->gcg->masterprob == relaxdata->gcg->bendersmasterprob);

   /* solve the next node in the master problem */
   SCIPdebugMessage("Solving node %"SCIP_LONGINT_FORMAT"'s relaxation.\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   /* prior to performing the decomposition the original problem verbosity is changed to NONE. This avoids output from
    * the original problem before the decomposition output. Once the decomposition has been performed, then the
    * verbosity level of the original problem is returned to the original verbosity level.
    */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", relaxdata->origverblevel) );
   SCIP_CALL( SCIPsetIntParam(masterprob, "display/verblevel", relaxdata->origverblevel) );

   /* getting the node limit from the original problem. This is because the master problem is solved to optimality in
    * the execution of the relaxator.
    */
   SCIP_CALL( SCIPgetLongintParam(scip, "limits/nodes", &nodelimit) );

   /* solving the master problem */
   SCIP_CALL( solveMasterProblem(gcg, masterprob, relaxdata, nodelimit, lowerbound, result) );

   /* if the master problem has been detected as infeasible, then the result must be set to SCIP_CUTOFF. */
   if( SCIPgetStatus(masterprob) == SCIP_STATUS_INFEASIBLE )
      (*result) = SCIP_CUTOFF;

   /* if the master problem has been solved to optimality, the we cutoff the root node. This informs that original
    * problem that no further processing is required.
    */
   if( SCIPgetStatus(masterprob) == SCIP_STATUS_OPTIMAL )
   {
      (*result) = SCIP_CUTOFF;
   }

   /* if there is no primal solution for the original problem, then the master solution is transferred */
   if( SCIPgetBestSol(GCGgetMasterprob(gcg)) != NULL && relaxdata->lastmastersol != SCIPgetBestSol(GCGgetMasterprob(relaxdata->gcg)) )
   {
      SCIP_SOL* newsol;
      SCIP_Bool stored;

      relaxdata->lastmastersol = SCIPgetBestSol(GCGgetMasterprob(relaxdata->gcg));

      SCIP_CALL( GCGtransformMastersolToOrigsol(gcg, SCIPgetBestSol(GCGgetMasterprob(gcg)), &newsol, TRUE, NULL) );
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, TRUE, &stored) );
#else
      SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );
#endif
      /* only check failed solution if best master solution is valid */
      if( !stored && GCGmasterIsBestsolValid(gcg) )
      {
         SCIP_CALL( SCIPcheckSolOrig(scip, newsol, &stored, TRUE, TRUE) );
      }
      /** @bug The solution doesn't have to be accepted, numerics might bite us, so the transformation might fail.
       *  A remedy could be: Round the values or propagate changes or call a heuristic to fix it.
       */
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );

      if( stored )
         SCIPdebugMessage("  updated current best primal feasible solution.\n");
   }

   /* set the lower bound pointer */
   if( GCGmasterIsCurrentSolValid(gcg)
      && (SCIPgetStage(masterprob) == SCIP_STAGE_SOLVED || SCIPgetStage(masterprob) == SCIP_STAGE_SOLVING) )
   {
      *lowerbound = SCIPgetDualbound(masterprob);
   }

   /* if the time, memory or node limit is hit in the Original or Benders mode, then we need to interrupt the solve.
    * This is required because the original problem is not solved in either of these modes, so it is not certain that
    * the original SCIP will also exceed the limit (definitely not for the node limit).
    */
   if( SCIPgetStatus(masterprob) == SCIP_STATUS_TIMELIMIT || SCIPgetStatus(masterprob) == SCIP_STATUS_NODELIMIT
      || SCIPgetStatus(masterprob) == SCIP_STATUS_MEMLIMIT )
   {
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   /* if the result pointer is DIDNOTRUN, this implies that the master problem was interrupted during solving. Since
    * Benders' decomposition uses a one-tree approach, then the user limits must be adhered to. This means, the if a
    * limit is exceeded, this is still a success for the solving.
    */
   if( (*result) == SCIP_DIDNOTRUN )
      (*result) = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** execution method of the relaxator for Benders' decomposition */
static
SCIP_RETCODE relaxExecGcgBendersDecomposition(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP_RELAX*           relax,              /**< the relaxator */
   SCIP_Real*            lowerbound,         /**< the lowerbound computed by the relaxator for the current node */
   SCIP_RESULT*          result              /**< the result of the relaxation call */
   )
{
   assert(gcg != NULL);
   assert(relax != NULL);
   assert(result != NULL);

   SCIP_CALL( solveMasterProblemAndEvaluate(gcg, relax, lowerbound, result) );

   return SCIP_OKAY;
}

/** execution method of the relaxator when the original problem is solved directly */
static
SCIP_RETCODE relaxExecGcgOriginalProblem(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP_RELAX*           relax,              /**< the relaxator */
   SCIP_Real*            lowerbound,         /**< the lowerbound computed by the relaxator for the current node */
   SCIP_RESULT*          result              /**< the result of the relaxation call */
   )
{
   assert(gcg != NULL);
   assert(relax != NULL);
   assert(result != NULL);

   SCIP_CALL( solveMasterProblemAndEvaluate(gcg, relax, lowerbound, result) );

   return SCIP_OKAY;
}


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecGcg)
{
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);
   assert(relax != NULL);
   assert(result != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   SCIP_CALL( initializeMasterProblemSolve(relaxdata->gcg, relax) );

   if( !SCIPisLPConstructed(scip) && !SCIPinProbing(scip) ) {
      SCIP_Bool cutoff;
      /* construct the LP in the original problem */
      SCIP_CALL(SCIPconstructLP(scip, &cutoff));
      assert(!cutoff);
      SCIP_CALL(SCIPflushLP(scip));
   }

   /* selecting the solving algorithm based upon the decomposition mode selected by the user, or whether the original
    * problem should be solved directly
    */
   if( GCGgetDecompositionMode(relaxdata->gcg) == GCG_DECMODE_ORIGINAL )
   {
      SCIP_CALL( GCGrestoreLimitSettings(relaxdata->gcg) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "There are no pricing problems in the decomposition. The original problem will be solved directly.\n");
      SCIP_CALL( relaxExecGcgOriginalProblem(relaxdata->gcg, relax, lowerbound, result) );
   }
   else if( relaxdata->mode == GCG_DECMODE_DANTZIGWOLFE )
   {
      SCIP_CALL( relaxExecGcgDantzigWolfe(relaxdata->gcg, relax, lowerbound, result) );
   }
   else if( relaxdata->mode == GCG_DECMODE_BENDERS )
   {
      SCIP_CALL( GCGrestoreLimitSettings(relaxdata->gcg) );
      SCIP_CALL( relaxExecGcgBendersDecomposition(relaxdata->gcg, relax, lowerbound, result) );
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "Sorry, the automatic selection is not currently available\n");
   }

   assert(*result == SCIP_CUTOFF || !SCIPisInfinity(scip, *lowerbound));

   return SCIP_OKAY;
}

#define relaxCopyGcg NULL
#define relaxInitGcg NULL


/*
 * relaxator specific interface methods
 */

/** creates the GCG relaxator and includes it in SCIP */
SCIP_RETCODE GCGincludeRelaxGcg(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP* origprob;

   origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

#ifdef WITH_BLISS
   {
      char name[SCIP_MAXSTRLEN];
      GCGgetBlissName(name, SCIP_MAXSTRLEN);
      SCIP_CALL( SCIPincludeExternalCodeInformation(origprob, name, "A Tool for Computing Automorphism Groups of Graphs by T. Junttila and P. Kaski (http://www.tcs.hut.fi/Software/bliss/)") );
   }
#endif

#ifdef WITH_NAUTY
   {
      char name[SCIP_MAXSTRLEN];
      GCGgetNautyName(name, SCIP_MAXSTRLEN);
      SCIP_CALL( SCIPincludeExternalCodeInformation(origprob, name, "A Tool for Computing Automorphism Groups of Graphs by B.D. McKay and A. Piperno (https://pallini.di.uniroma1.it/)") );
   }
#endif

#ifdef WITH_CLIQUER
      SCIP_CALL( SCIPincludeExternalCodeInformation(origprob, "Cliquer", "A set of C routines for finding cliques in an arbitrary weighted graph by S. Niskanen and P. Ostergard (https://users.aalto.fi/~pat/cliquer.html)") );
#endif

   /* create GCG relaxator data */
   SCIP_CALL( SCIPallocMemory(origprob, &relaxdata) );
   SCIP_CALL( initRelaxdata(gcg, relaxdata) );

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelax(origprob, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, relaxCopyGcg, relaxFreeGcg, relaxInitGcg,
         relaxExitGcg, relaxInitsolGcg, relaxExitsolGcg, relaxExecGcg, relaxdata) );
   relax = SCIPfindRelax(origprob, RELAX_NAME);
   assert(relax != NULL);
   gcg->relax = relax;

   /* inform the main scip, that no LPs should be solved */
   SCIP_CALL( SCIPsetIntParam(origprob, "lp/solvefreq", 0) );

   /* Disable restarts */
   SCIP_CALL( SCIPsetIntParam(origprob, "presolving/maxrestarts", 0) );
   SCIP_CALL( SCIPsetBoolParam(origprob, "misc/calcintegral", FALSE) );

   /* initialize the scip data structure for the master problem. The master problem is initialized as the Dantzig-Wolfe
    * master problem. The alternate master problem is initialized as the Benders' decomposition master problem.
    */
   SCIP_CALL( SCIPcreate(&(gcg->dwmasterprob)) );
   gcg->masterprob = gcg->dwmasterprob;
   SCIP_CALL( GCGincludePricerGcg(relaxdata->gcg) );
   SCIP_CALL( GCGincludeMasterPlugins(gcg) );
   SCIP_CALL( SCIPsetMessagehdlr(gcg->masterprob, SCIPgetMessagehdlr(origprob)) );

   /* getting the verbosity level of the original problem */
   SCIP_CALL( SCIPgetIntParam(origprob, "display/verblevel", &relaxdata->origverblevel) );

   /* disable display output in the master problem */
   SCIP_CALL( SCIPsetIntParam(gcg->masterprob, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );

   /* set parameters in master problem */
   SCIP_CALL( SCIPsetIntParam(gcg->masterprob, "pricing/maxvars", INT_MAX) );
   SCIP_CALL( SCIPsetIntParam(gcg->masterprob, "pricing/maxvarsroot", INT_MAX) );
   SCIP_CALL( SCIPsetRealParam(gcg->masterprob, "pricing/abortfac", 1.0) );
   SCIP_CALL( SCIPsetIntParam(gcg->masterprob, "lp/disablecutoff", 1) );
#ifdef DELVARS
   /* set paramteters to allow deletion of variables */
   SCIP_CALL( SCIPsetBoolParam(gcg->masterprob, "pricing/delvars", TRUE) );
   SCIP_CALL( SCIPsetBoolParam(gcg->masterprob, "pricing/delvarsroot", TRUE) );
   SCIP_CALL( SCIPsetBoolParam(gcg->masterprob, "lp/cleanupcols", TRUE) );
   SCIP_CALL( SCIPsetBoolParam(gcg->masterprob, "lp/cleanupcolsroot", TRUE) );
#endif

   /* initializing the alternate master problem. The alternate master problem is initially the Benders' decomposition
    * master problem
    */
   SCIP_CALL( SCIPcreate(&(gcg->bendersmasterprob)) );
   SCIP_CALL( GCGincludeBendersGcg(gcg) );
   SCIP_CALL( GCGincludeBendersPlugins(gcg) );
   SCIP_CALL( SCIPsetMessagehdlr(gcg->bendersmasterprob, SCIPgetMessagehdlr(origprob)) );

   SCIP_CALL( SCIPsetIntParam(gcg->bendersmasterprob, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
   SCIP_CALL( SCIPsetBoolParam(gcg->bendersmasterprob, "display/relevantstats", FALSE) );

   /* add GCG relaxator parameters */
   SCIP_CALL( SCIPaddBoolParam(origprob, "relaxing/gcg/discretization",
         "should discretization (TRUE) or convexification (FALSE) approach be used?",
         &(relaxdata->discretization), FALSE, DEFAULT_DISCRETIZATION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(origprob, "relaxing/gcg/mipdiscretization",
         "should discretization (TRUE) or convexification (FALSE) approach be used in mixed-integer programs?",
         &(relaxdata->mipdiscretization), FALSE, DEFAULT_MIPDISCRETIZATION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(origprob, "relaxing/gcg/aggregation/enabled",
         "should identical blocks be aggregated (only for discretization approach)?",
         &(relaxdata->aggregation), FALSE, DEFAULT_AGGREGATION, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(origprob, "relaxing/gcg/aggregation/limitnconssperblock",
         "Limits the number of constraints of a block (aggregation information for block is not calculated when exceeded)",
         &(relaxdata->aggregationnconsslimit), FALSE, DEFAULT_AGGREGATIONNCONSSLIMIT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(origprob, "relaxing/gcg/aggregation/limitnvarsperblock",
         "Limits the number of variables of a block (aggregation information for block is not calculated when exceeded)",
         &(relaxdata->aggregationnvarslimit), FALSE, DEFAULT_AGGREGATIONNVARSLIMIT, 0, INT_MAX, NULL, NULL) );
#ifndef NO_AUT_LIB
   SCIP_CALL( SCIPaddBoolParam(origprob, "relaxing/gcg/aggregation/usesymmetrylib",
         "should a symmetry detection library be used to check for identical blocks?",
         &(relaxdata->usesymmetrylib), FALSE, DEFAULT_BLISS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(origprob, "relaxing/gcg/aggregation/searchnodelimit",
         "search node limit (0: unlimited), requires patched bliss version",
         &(relaxdata->searchnodelimit), TRUE, (int)DEFAULT_BLISS_SEARCH_NODE_LIMIT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(origprob, "relaxing/gcg/aggregation/generatorlimit",
         "generator limit (0: unlimited), requires patched bliss version or version >= 0.76",
         &(relaxdata->generatorlimit), TRUE, (int)DEFAULT_BLISS_GENERATOR_LIMIT, 0, INT_MAX, NULL, NULL) );
#else
   relaxdata->usesymmetrylib = FALSE;
   relaxdata->searchnodelimit = 0;
   relaxdata->generatorlimit = 0;
#endif
   SCIP_CALL( SCIPaddBoolParam(origprob, "relaxing/gcg/dispinfos",
         "should additional information about the blocks be displayed?",
         &(relaxdata->dispinfos), FALSE, DEFAULT_DISPINFOS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(origprob, "relaxing/gcg/mode",
            "the decomposition mode that GCG will use. (0: Dantzig-Wolfe (default), 1: Benders' decomposition, "
            "2: no decomposition will be performed)",
            (int*)&(relaxdata->mode), FALSE, (int)DEFAULT_MODE, 0, 2, NULL, NULL) );

   return SCIP_OKAY;
}


/*
 * relaxator specific interface methods for coordination of branching rules
 */

/** includes a branching rule into the relaxator data */
SCIP_RETCODE GCGrelaxIncludeBranchrule(
   GCG*                  gcg,                /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule for which callback methods are saved */
   GCG_BRANCHRULE**      gcgbranchrule,      /**< pointer to store created GCG branch rule (can be NULL) */
   GCG_DECL_BRANCHACTIVEMASTER((*branchactivemaster)),/**<  activation method for branchrule */
   GCG_DECL_BRANCHDEACTIVEMASTER((*branchdeactivemaster)),/**<  deactivation method for branchrule */
   GCG_DECL_BRANCHPROPMASTER((*branchpropmaster)),/**<  propagation method for branchrule */
   GCG_DECL_BRANCHMASTERSOLVED((*branchmastersolved)),/**<  master solved method for branchrule */
   GCG_DECL_BRANCHDATADELETE((*branchdatadelete)),/**<  branchdata deletion method for branchrule */
   GCG_DECL_BRANCHNEWCOL ((*branchnewcol)),  /**< new column handler method of branching rule */
   GCG_DECL_BRANCHGETEXTENDEDMASTERCONS ((*branchgetextendedmastercons)), /**< extended master cons getter of branching rule */
   GCG_DECL_BRANCHGETEXTENDEDMASTERCONSCOEFF((*branchgetextendedmasterconscoeff)) /**< column coefficient calculation method for extended master conss */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int pos;

   assert(gcg != NULL);
   assert(branchrule != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   SCIP_CALL( ensureSizeBranchrules(gcg, relaxdata) );

   pos = relaxdata->nbranchrules;

   /* store callback functions */
   SCIP_CALL( SCIPallocMemory(origprob, &(relaxdata->branchrules[pos])) ); /*lint !e866*/
   relaxdata->branchrules[pos]->branchrule = branchrule;
   relaxdata->branchrules[pos]->branchactivemaster = branchactivemaster;
   relaxdata->branchrules[pos]->branchdeactivemaster = branchdeactivemaster;
   relaxdata->branchrules[pos]->branchpropmaster = branchpropmaster;
   relaxdata->branchrules[pos]->branchmastersolved = branchmastersolved;
   relaxdata->branchrules[pos]->branchdatadelete = branchdatadelete;
   relaxdata->branchrules[pos]->branchnewcol = branchnewcol;
   relaxdata->branchrules[pos]->branchgetextendedmastercons = branchgetextendedmastercons;
   relaxdata->branchrules[pos]->branchgetextendedmasterconscoeff = branchgetextendedmasterconscoeff;
   relaxdata->nbranchrules++;

   if( gcgbranchrule != NULL )
      *gcgbranchrule = relaxdata->branchrules[pos];

   return SCIP_OKAY;
}

/** perform activation method of the given branchrule for the given branchdata */
SCIP_RETCODE GCGrelaxBranchActiveMaster(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata          /**< data representing the branching decision */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;

   assert(gcg != NULL);
   assert(branchrule != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         /* call activation method of branching rule */
         if( relaxdata->branchrules[i]->branchactivemaster != NULL )
            SCIP_CALL( relaxdata->branchrules[i]->branchactivemaster(gcg, branchdata) );

         SCIP_CALL( addActiveBranchExtendedmastercons(gcg, relaxdata, relaxdata->branchrules[i], branchdata) );

         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

   return SCIP_OKAY;
}

/** perform deactivation method of the given branchrule for the given branchdata */
SCIP_RETCODE GCGrelaxBranchDeactiveMaster(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata          /**< data representing the branching decision */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;

   assert(gcg != NULL);
   assert(branchrule != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         /* call deactivation method of branching rule */
         if( relaxdata->branchrules[i]->branchdeactivemaster != NULL )
            SCIP_CALL( relaxdata->branchrules[i]->branchdeactivemaster(gcg, branchdata) );

         SCIP_CALL( dropActiveBranchExtendedmastercons(gcg, relaxdata, relaxdata->branchrules[i], branchdata) );

         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

   return SCIP_OKAY;
}

/** perform propagation method of the given branchrule for the given branchdata */
SCIP_RETCODE GCGrelaxBranchPropMaster(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation call */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;

   assert(gcg != NULL);
   assert(branchrule != NULL);
   assert(result != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   *result = SCIP_DIDNOTRUN;

   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         /* call propagation method of branching rule*/
         if( relaxdata->branchrules[i]->branchpropmaster != NULL )
            SCIP_CALL( relaxdata->branchrules[i]->branchpropmaster(gcg, branchdata, result) );

         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

   return SCIP_OKAY;
}

/** frees branching data created by the given branchrule */
SCIP_RETCODE GCGrelaxBranchDataDelete(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA**      branchdata,         /**< data representing the branching decision */
   SCIP_Bool             origbranch,         /**< true iff an origbranch triggered this call */
   SCIP_Bool             force               /**< branch data must be deleted if true */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;

   assert(gcg != NULL);
   assert(branchrule != NULL);
   assert(branchdata != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         /* call branchrule data deletion method of the branching rule */
         if( relaxdata->branchrules[i]->branchdatadelete != NULL )
            SCIP_CALL( relaxdata->branchrules[i]->branchdatadelete(gcg, branchdata, origbranch, force) );
         else
         {
            if( *branchdata != NULL )
            {
               SCIPfreeMemory(GCGgetMasterprob(gcg), branchdata);
            }
         }
         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

   return SCIP_OKAY;
}

/** notifies the branching rule that a new mastervariable was created while this node was active */
SCIP_RETCODE GCGrelaxBranchNewCol(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   SCIP_VAR*             mastervar           /**< new mastervariable that was created */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;

   assert(gcg != NULL);
   assert(branchrule != NULL);
   assert(branchdata != NULL);
   assert(mastervar != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         SCIP_CALL( GCGrelaxBranchNewColWithGCGBranchrule(gcg, relaxdata->branchrules[i], branchdata, mastervar) );

         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

   return SCIP_OKAY;
}

/** notifies the branching rule that a new mastervariable was created while this node was active */
SCIP_RETCODE GCGrelaxBranchNewColWithGCGBranchrule(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE*       branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   SCIP_VAR*             mastervar           /**< new mastervariable that was created */
   )
{
   assert(gcg != NULL);
   assert(branchrule != NULL);
   assert(branchdata != NULL);
   assert(mastervar != NULL);

   /* call new mastervariable handler method of branching rule*/
   if( branchrule->branchnewcol != NULL )
      SCIP_CALL( branchrule->branchnewcol(gcg, branchdata, mastervar) );

   return SCIP_OKAY;
}

/** gets the extendedmasterconsdata created by this branching rule, if any */
SCIP_RETCODE GCGrelaxBranchGetExtendedMasterCons(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   GCG_EXTENDEDMASTERCONSDATA**   extendedmasterconsdata       /**< the extendedmasterconsdata to grab */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;

   assert(gcg != NULL);
   assert(branchrule != NULL);
   assert(extendedmasterconsdata == NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         /* call extended master cons getter method of branching rule */
         if( relaxdata->branchrules[i]->branchgetextendedmastercons != NULL )
         {
            SCIP_CALL( relaxdata->branchrules[i]->branchgetextendedmastercons(gcg, branchdata, extendedmasterconsdata) );
            assert(*extendedmasterconsdata != NULL);
         }

         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

   return SCIP_OKAY;
}

/** get extended master conss of all active nodes */
SCIP_RETCODE GCGrelaxBranchGetAllActiveExtendedMasterConss(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_BRANCHRULE***     branchrules,        /**< branching rules that created extended master conss */
   GCG_BRANCHDATA***     branchdata,         /**< data represeting the branching decisions of the active nodes */
   GCG_EXTENDEDMASTERCONSDATA***  extendedmasterconsdata,      /**< array of extended master conss generated by branching in all currently active nodes */
   int*                  nextendedmasterconss         /**< number of currently active branching rules that created extended master conss */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   *branchrules = relaxdata->activebranchrules;
   *branchdata = relaxdata->activebranchdata;
   *extendedmasterconsdata = relaxdata->activebranchextendedmasterconss;
   *nextendedmasterconss = relaxdata->nactivebranchextendedmasterconss;

   return SCIP_OKAY;
}

/** perform method of the given branchrule that is called after the master LP is solved */
SCIP_RETCODE GCGrelaxBranchMasterSolved(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule that did the branching */
   GCG_BRANCHDATA*       branchdata,         /**< data representing the branching decision */
   SCIP_Real             newlowerbound       /**< the new local lowerbound */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   int i;

   assert(gcg != NULL);
   assert(branchrule != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* search for the branching rule in the branchrules array */
   for( i = 0; i < relaxdata->nbranchrules; i++ )
   {
      if( branchrule == relaxdata->branchrules[i]->branchrule )
      {
         /* call master problem solved method of the branching rule */
         if( relaxdata->branchrules[i]->branchmastersolved != NULL )
            SCIP_CALL( relaxdata->branchrules[i]->branchmastersolved(gcg, branchdata, newlowerbound) );

         break;
      }
   }

   assert(i < relaxdata->nbranchrules);

   return SCIP_OKAY;
}

/** transforms a constraint of the original problem into the master variable space
 *  and stores information about the constraints in the variable */
SCIP_RETCODE GCGrelaxTransOrigToMasterCons(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_CONS*            cons,               /**< the constraint that should be transformed */
   SCIP_CONS**           transcons           /**< pointer to store the transformed constraint */
   )
{
   SCIP* scip;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_CONS* mastercons;
   char name[SCIP_MAXSTRLEN];

   SCIP_VAR** mastervars;
   int nmastervars;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;
   int v;
   int i;
   int j;

   assert(gcg != NULL);
   assert(cons != NULL);

   scip = GCGgetOrigprob(gcg);
   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* create and add corresponding linear constraint in the master problem */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "m_%s", SCIPconsGetName(cons));
   SCIP_CALL( SCIPcreateConsLinear(GCGgetMasterprob(relaxdata->gcg), &mastercons, name, 0, NULL, NULL,
         GCGconsGetLhs(scip, cons), GCGconsGetRhs(scip, cons),
         TRUE, TRUE, TRUE, TRUE, TRUE, SCIPconsIsLocal(cons), TRUE, FALSE, FALSE,
         SCIPconsIsStickingAtNode(cons)) );

   /* now compute coefficients of the master variables in the master constraint */
   mastervars = SCIPgetVars(GCGgetMasterprob(gcg));
   nmastervars = SCIPgetNVars(GCGgetMasterprob(gcg));

   consvars = SCIPgetVarsLinear(scip, cons);
   nconsvars = SCIPgetNVarsLinear(scip, cons);
   consvals = SCIPgetValsLinear(scip, cons);


   /* add coefs of the original variables in the constraint to their variable data */
   for( v = 0; v < nconsvars; v++ )
   {
      SCIP_CALL( GCGoriginalVarAddCoef(gcg, consvars[v], consvals[v], mastercons) );
   }

   /* add master variables to the corresponding master constraint */
   for( v = 0; v < nmastervars; v++ )
   {
      SCIP_VAR** origvars;
      SCIP_Real* origvals;
      int norigvars;
      SCIP_Real coef = 0.0;

      origvars = GCGmasterVarGetOrigvars(mastervars[v]);
      norigvars = GCGmasterVarGetNOrigvars(mastervars[v]);
      origvals = GCGmasterVarGetOrigvals(mastervars[v]);

      for( i = 0; i < norigvars; i++ )
         for( j = 0; j < nconsvars; j++ )
            if( consvars[j] == origvars[i] )
               coef += consvals[j] * origvals[i];

      if( !SCIPisFeasZero(scip, coef) )
      {
         SCIP_CALL( SCIPaddCoefLinear(GCGgetMasterprob(gcg), mastercons, mastervars[v], coef) );
      }
   }

   /* store the constraints in the arrays origmasterconss and masterconss in the problem data */
   SCIP_CALL( ensureSizeMasterConss(gcg, relaxdata, relaxdata->nmasterconss+1) );
   SCIP_CALL( SCIPcaptureCons(scip, cons) );
   relaxdata->origmasterconss[relaxdata->nmasterconss] = cons;
   relaxdata->masterconss[relaxdata->nmasterconss] = mastercons;

   SCIP_CALL( GCGmasterAddMasterconsToHashmap(gcg, relaxdata->masterconss[relaxdata->nmasterconss],
         relaxdata->nmasterconss) );

   relaxdata->nmasterconss++;

   *transcons = mastercons;

   return SCIP_OKAY;
}

/** returns the pricing problem of the given number */
SCIP* GCGgetPricingprob(
   GCG*                  gcg,                /**< GCG data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->pricingprobs[pricingprobnr];
}

/** returns the number of relevant pricing problems */
int GCGgetNRelPricingprobs(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   assert(relaxdata->nrelpricingprobs >= -1);
   return relaxdata->nrelpricingprobs;
}

/** returns the number of pricing problems */
int GCGgetNPricingprobs(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   assert(relaxdata->npricingprobs >= -1);
   return relaxdata->npricingprobs;
}

/** returns TRUE iff the pricing problem of the given number is relevant, that means is not identical to
 *  another and represented by it */
SCIP_Bool GCGisPricingprobRelevant(
   GCG*                  gcg,                /**< GCG data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return (relaxdata->blockrepresentative[pricingprobnr] == pricingprobnr);

}

/**
 *  for a given block, return the block by which it is represented
 */
int GCGgetBlockRepresentative(
   GCG*                  gcg,                /**< GCG data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   if( pricingprobnr == -1 )
      return -1;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   assert(pricingprobnr >= 0);
   assert(pricingprobnr < relaxdata->npricingprobs);
   assert(relaxdata->nblocksidentical[pricingprobnr] >= 0);
   assert((relaxdata->blockrepresentative[pricingprobnr] == pricingprobnr)
      == (relaxdata->nblocksidentical[pricingprobnr] > 0));

   return relaxdata->blockrepresentative[pricingprobnr];
}

/** returns the number of blocks in the original formulation, that are represented by
 *  the pricingprob with the given number */
int GCGgetNIdenticalBlocks(
   GCG*                  gcg,                /**< GCG data structure */
   int                   pricingprobnr       /**< number of the pricing problem */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);
   assert(pricingprobnr >= 0);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(pricingprobnr <= relaxdata->npricingprobs);
   assert(relaxdata->nblocksidentical[pricingprobnr] >= 0);
   assert((relaxdata->blockrepresentative[pricingprobnr] == pricingprobnr)
      == (relaxdata->nblocksidentical[pricingprobnr] > 0));

   return relaxdata->nblocksidentical[pricingprobnr];

}

/** returns the number of constraints in the master problem */
int GCGgetNMasterConss(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->nmasterconss;
}

/** returns the contraints in the master problem */
SCIP_CONS** GCGgetMasterConss(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->masterconss;
}

/** returns the linking constraints in the original problem that correspond to the constraints in the master problem */
SCIP_CONS** GCGgetOrigMasterConss(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->origmasterconss;
}

/** returns the convexity constraint for the given block */
SCIP_CONS* GCGgetConvCons(
   GCG*                  gcg,                /**< GCG data structure */
   int                   blocknr             /**< the number of the block for which we
                                              *   need the convexity constraint */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);
   assert(blocknr >= 0);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(blocknr < relaxdata->npricingprobs);

   return relaxdata->convconss[blocknr];
}

/** returns the visualization parameters */
GCG_PARAMDATA* GCGgetParamsVisu(
   GCG*                  gcg                /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   GCG_PARAMDATA* paramdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->paramsvisu != NULL);

   paramdata = relaxdata->paramsvisu;
   assert(paramdata != NULL);

   return paramdata;
}

/** returns the current solution for the original problem */
SCIP_SOL* GCGrelaxGetCurrentOrigSol(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->currentorigsol;
}

/** returns whether the current solution is primal feasible in the original problem */
SCIP_Bool GCGrelaxIsOrigSolFeasible(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->origsolfeasible;
}

/** returns whether the master problem is a set covering problem */
SCIP_Bool GCGisMasterSetCovering(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->masterissetcover;
}

/** returns whether the master problem is a set partitioning problem */
SCIP_Bool GCGisMasterSetPartitioning(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->masterissetpart;
}

/** start probing mode on both the original and master problems
 *
 *  @note This mode is intended for working on the original variables but using the master LP;
 *        it currently only supports bound changes on the original variables,
 *        but no additional rows
 */
SCIP_RETCODE GCGrelaxStartProbing(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_HEUR*            probingheur         /**< heuristic that started probing mode, or NULL */
   )
{
   SCIP* scip;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP* masterprob;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);
   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if( relaxdata->masterinprobing )
   {
      SCIPerrorMessage("already in GCG probing mode\n");
      return SCIP_INVALIDCALL;
   }

   masterprob = GCGgetMasterprob(relaxdata->gcg);
   assert(masterprob != NULL);

   /* start probing in both the original and the master problem */
   SCIP_CALL( SCIPstartProbing(scip) );
   SCIP_CALL( SCIPstartProbing(masterprob) );

   relaxdata->masterinprobing = TRUE;
   relaxdata->probingheur = probingheur;

   /* remember the current original solution */
   assert(relaxdata->storedorigsol == NULL);
   if( relaxdata->currentorigsol != NULL )
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &relaxdata->storedorigsol, relaxdata->currentorigsol) );
      relaxdata->storedfeasibility = relaxdata->origsolfeasible;
   }

   return SCIP_OKAY;
}

/** returns the  heuristic that started probing in the master problem, or NULL */
SCIP_HEUR* GCGrelaxGetProbingheur(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->probingheur;
}

/** add a new probing node the original problem together with an original branching constraint
 *
 *  @note A corresponding probing node must be added to the master problem right before solving the probing LP
 */
SCIP_RETCODE GCGrelaxNewProbingnodeOrig(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_CONS* probingcons;
   SCIP_NODE* probingnode;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);
   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if( !relaxdata->masterinprobing )
   {
      SCIPerrorMessage("not in GCG probing mode\n");
      return SCIP_INVALIDCALL;
   }

   if( SCIPgetProbingDepth(scip) != SCIPgetProbingDepth(GCGgetMasterprob(gcg)) )
   {
      SCIPerrorMessage("original and master problem not at same probing depth\n");
      return SCIP_INVALIDCALL;
   }

   /* add a probing node in the original problem together with an original branching constraint */
   SCIP_CALL( SCIPnewProbingNode(scip) );
   probingnode = SCIPgetCurrentNode(scip);
   SCIP_CALL( GCGcreateConsOrigbranch(gcg, &probingcons, "probingcons", probingnode,
      GCGconsOrigbranchGetActiveCons(gcg), NULL, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, probingnode, probingcons, NULL) );
   SCIP_CALL( SCIPreleaseCons(scip, &probingcons) );


   return SCIP_OKAY;
}

/** add a new probing node the master problem together with a master branching constraint
 *  which ensures that bound changes are transferred to master and pricing problems
 *
 *  @note A corresponding probing node must have been added to the original problem beforehand;
 *        furthermore, this method must be called after bound changes to the original problem have been made
 */
SCIP_RETCODE GCGrelaxNewProbingnodeMaster(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP* masterprob;
   SCIP_CONS* probingcons;
   SCIP_NODE* probingnode;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);
   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if( !relaxdata->masterinprobing )
   {
      SCIPerrorMessage("not in GCG probing mode\n");
      return SCIP_INVALIDCALL;
   }

   masterprob = GCGgetMasterprob(relaxdata->gcg);
   assert(masterprob != NULL);

   if( SCIPgetProbingDepth(scip) != SCIPgetProbingDepth(masterprob) + 1 )
   {
      SCIPerrorMessage("master probing node must be created after original probing node\n");
      return SCIP_INVALIDCALL;
   }

   /* add a probing node in the master problem together with a master branching constraint */
   SCIP_CALL( SCIPnewProbingNode(masterprob) );
   probingnode = SCIPgetCurrentNode(masterprob);
   assert(GCGconsMasterbranchGetActiveCons(gcg) != NULL);
   SCIP_CALL( GCGcreateConsMasterbranch(gcg, &probingcons, "mprobingcons", probingnode,
      GCGconsMasterbranchGetActiveCons(gcg), NULL, NULL, NULL, 0, 0) );
   SCIP_CALL( SCIPaddConsNode(masterprob, probingnode, probingcons, NULL) );
   SCIP_CALL( SCIPreleaseCons(masterprob, &probingcons) );

   return SCIP_OKAY;
}

/** add a new probing node the master problem together with a master branching constraint
 *  which ensures that bound changes are transferred to master and pricing problems as well as additional
 *  constraints
 *
 *  @note A corresponding probing node must have been added to the original problem beforehand;
 *        furthermore, this method must be called after bound changes to the original problem have been made
 */
SCIP_RETCODE GCGrelaxNewProbingnodeMasterCons(
   GCG*                  gcg,                 /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< pointer to the branching rule */
   GCG_BRANCHDATA*       branchdata,         /**< branching data */
   SCIP_CONS**           origbranchconss,    /**< original constraints enforcing the branching decision */
   int                   norigbranchconss,   /**< number of original constraints */
   int                   maxorigbranchconss  /**< capacity of origbranchconss */
   )
{
   SCIP* scip;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP* masterprob;
   SCIP_CONS* probingcons;
   SCIP_NODE* probingnode;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);
   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if( !relaxdata->masterinprobing )
   {
      SCIPerrorMessage("not in GCG probing mode\n");
      return SCIP_INVALIDCALL;
   }

   masterprob = GCGgetMasterprob(relaxdata->gcg);
   assert(masterprob != NULL);

   if( SCIPgetProbingDepth(scip) != SCIPgetProbingDepth(masterprob) + 1 )
   {
      SCIPerrorMessage("master probing node must be created after original probing node\n");
      return SCIP_INVALIDCALL;
   }

   /* add a probing node in the master problem together with a master branching constraint */
   SCIP_CALL( SCIPnewProbingNode(masterprob) );
   probingnode = SCIPgetCurrentNode(masterprob);
   assert(GCGconsMasterbranchGetActiveCons(gcg) != NULL);
   SCIP_CALL( GCGcreateConsMasterbranch(relaxdata->gcg, &probingcons, "mprobingcons", probingnode,
      GCGconsMasterbranchGetActiveCons(gcg), branchrule, branchdata, origbranchconss, norigbranchconss, maxorigbranchconss) );
   SCIP_CALL( SCIPaddConsNode(masterprob, probingnode, probingcons, NULL) );
   SCIP_CALL( SCIPreleaseCons(masterprob, &probingcons) );

   return SCIP_OKAY;
}

/** add probing nodes to both the original and master problem;
 *  furthermore, add origbranch and masterbranch constraints to transfer branching decisions
 *  from the original to the master problem
 */
SCIP_RETCODE GCGrelaxBacktrackProbing(
   GCG*                  gcg,                /**< GCG data structure */
   int                   probingdepth        /**< probing depth of the node in the probing path that should be reactivated */
   )
{
   SCIP* scip;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP* masterprob;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);
   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if( !relaxdata->masterinprobing )
   {
      SCIPerrorMessage("not in GCG probing mode\n");
      return SCIP_INVALIDCALL;
   }

   masterprob = GCGgetMasterprob(relaxdata->gcg);
   assert(masterprob != NULL);

   SCIP_CALL( SCIPbacktrackProbing(scip, probingdepth) );
   SCIP_CALL( SCIPbacktrackProbing(masterprob, probingdepth) );

   return SCIP_OKAY;
}

/** solve the master probing LP with or without pricing */
static
SCIP_RETCODE performProbing(
   GCG*                  gcg,                /**< GCG data structure */
   int                   maxlpiterations,    /**< maximum number of lp iterations allowed */
   int                   maxpricerounds,     /**< maximum number of pricing rounds allowed */
   SCIP_Bool             usepricing,         /**< should the LP be solved with or without pricing? */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of performed LP iterations (or NULL) */
   int*                  npricerounds,       /**< pointer to store the number of performed pricing rounds (or NULL) */
   SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP* masterprob;
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_Longint oldnlpiters;
   int oldpricerounds;
   SCIP_Longint nodelimit;

   assert(gcg != NULL);

   /* get the relaxator */
   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   /* get the relaxator data */
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if( !relaxdata->masterinprobing )
   {
      SCIPerrorMessage("not in GCG probing mode\n");
      return SCIP_INVALIDCALL;
   }

   /* get master problem */
   masterprob = GCGgetMasterprob(relaxdata->gcg);
   assert(masterprob != NULL);

   /* increase node limit for the master problem by 1 */
   SCIP_CALL( SCIPgetLongintParam(masterprob, "limits/nodes", &nodelimit) );
   SCIP_CALL( SCIPsetLongintParam(masterprob, "limits/nodes", nodelimit + 1) );

   /* propagate probing bound changes to the master problem */
   SCIP_CALL( SCIPpropagateProbing(masterprob, -1, cutoff, NULL) );
   assert(!(*cutoff));

   /* remember LP iterations and pricing rounds before LP solving */
   oldnlpiters = SCIPgetNLPIterations(masterprob);
   oldpricerounds = SCIPgetNPriceRounds(masterprob);

   *lpobjvalue = 0.0;
   *lpsolved = FALSE;

   /* solve the probing LP */
   if( usepricing )
   {
      /* LP iterations are unlimited when probing LP is solved with pricing */
      assert(maxlpiterations == -1);
      SCIP_CALL( SCIPsolveProbingLPWithPricing(masterprob, FALSE, TRUE, maxpricerounds, lperror, NULL) );
   }
   else
   {
      assert(maxpricerounds == 0);
      SCIP_CALL( SCIPsolveProbingLP(masterprob, maxlpiterations, lperror, NULL) );
   }
   lpsolstat = SCIPgetLPSolstat(masterprob);

   /* reset the node limit */
   SCIP_CALL( SCIPsetLongintParam(masterprob, "limits/nodes", nodelimit) );

   /* calculate number of LP iterations and pricing rounds performed */
   if( nlpiterations != NULL )
      *nlpiterations = SCIPgetNLPIterations(masterprob) - oldnlpiters;
   if( npricerounds != NULL )
      *npricerounds = SCIPgetNPriceRounds(masterprob) - oldpricerounds;

   if( !(*lperror) )
   {
      /* get LP solution status, objective value */
      *cutoff = *cutoff || (lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT || lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         SCIPdebugMessage("lpobjval = %g\n", SCIPgetLPObjval(masterprob));
         *lpobjvalue = SCIPgetLPObjval(masterprob);
         *lpsolved = TRUE;
         SCIP_CALL( GCGrelaxUpdateCurrentSol(gcg) );
      }
   }
   else
   {
      SCIPdebugMessage("something went wrong, an lp error occurred\n");
   }

   return SCIP_OKAY;
}


/** solve the master probing LP without pricing */
SCIP_RETCODE GCGrelaxPerformProbing(
   GCG*                  gcg,                /**< GCG data structure */
   int                   maxlpiterations,    /**< maximum number of lp iterations allowed */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of performed LP iterations (or NULL) */
   SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   )
{
   SCIP_CALL( performProbing(gcg, maxlpiterations, 0, FALSE, nlpiterations,
         NULL, lpobjvalue, lpsolved, lperror, cutoff) );

   return SCIP_OKAY;
}


/** solve the master probing LP with pricing */
SCIP_RETCODE GCGrelaxPerformProbingWithPricing(
   GCG*                  gcg,                /**< GCG data structure */
   int                   maxpricerounds,     /**< maximum number of pricing rounds allowed */
   SCIP_Longint*         nlpiterations,      /**< pointer to store the number of performed LP iterations (or NULL) */
   int*                  npricerounds,       /**< pointer to store the number of performed pricing rounds (or NULL) */
   SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   )
{
   SCIP_CALL( performProbing(gcg, -1, maxpricerounds, TRUE, nlpiterations,
         npricerounds, lpobjvalue, lpsolved, lperror, cutoff) );

   return SCIP_OKAY;
}

/** end probing mode in both the original and master problems */
SCIP_RETCODE GCGrelaxEndProbing(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP* masterprob;

   SCIP_VAR** vars;
   int nvars;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);
   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if( !relaxdata->masterinprobing )
   {
      SCIPerrorMessage("not in GCG probing mode\n");
      return SCIP_INVALIDCALL;
   }

   masterprob = GCGgetMasterprob(relaxdata->gcg);
   assert(masterprob != NULL);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(vars != NULL);
   assert(nvars >= 0);

   SCIP_CALL( SCIPendProbing(masterprob) );
   SCIP_CALL( SCIPendProbing(scip) );

   relaxdata->masterinprobing = FALSE;
   relaxdata->probingheur = NULL;

   /* if a new primal solution was found in the master problem, transfer it to the original problem
    * @todo: this is probably not necessary anymore since it is done by an event handler
    */
   if( SCIPgetBestSol(masterprob) != NULL && relaxdata->lastmastersol != SCIPgetBestSol(masterprob) )
   {
      SCIP_SOL* newsol;
      SCIP_Bool stored;

      relaxdata->lastmastersol = SCIPgetBestSol(masterprob);

      SCIP_CALL( GCGtransformMastersolToOrigsol(gcg, relaxdata->lastmastersol, &newsol, TRUE, NULL) );

      SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );
      if( !stored )
      {
         SCIP_CALL( SCIPcheckSolOrig(scip, newsol, &stored, TRUE, TRUE) );
      }
      assert(stored);
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );

      SCIPdebugMessage("probing finished in master problem\n");
   }

   /* restore old relaxation solution and branching candidates */
   if( relaxdata->currentorigsol != NULL )
   {
      SCIPdebugMessage("Freeing previous solution origsol\n");
      SCIP_CALL( SCIPfreeSol(scip, &(relaxdata->currentorigsol)) );
   }
   SCIPclearExternBranchCands(scip);

   if( relaxdata->storedorigsol != NULL )
   {
      int i;

      SCIP_CALL( SCIPcreateSol(scip, &relaxdata->currentorigsol, NULL) );
      SCIP_CALL( SCIPsetRelaxSolValsSol(scip, relax, relaxdata->storedorigsol, RELAX_INCLUDESLP) );

      for( i = 0; i < nvars; i++ )
      {
         SCIP_VAR* var;
         SCIP_Real solval;

         var = vars[i];
         solval = SCIPgetSolVal(scip, relaxdata->storedorigsol, var);

         SCIP_CALL( SCIPsetSolVal(scip, relaxdata->currentorigsol, var, solval) );

         if( SCIPvarGetType(var) <= SCIP_VARTYPE_INTEGER && !SCIPisFeasIntegral(scip, solval) )
         {
            /* this was an assertion but I think it is ok to fail as the old solution is restored
             * and probing may happen directly after branching */
            if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
               SCIP_CALL( SCIPaddExternBranchCand(scip, var, solval - SCIPfloor(scip, solval), solval) );
         }
      }
      assert(SCIPisFeasEQ(scip, SCIPgetRelaxSolObj(scip), SCIPgetSolTransObj(scip, relaxdata->currentorigsol)));

      SCIP_CALL( SCIPfreeSol(scip, &relaxdata->storedorigsol) );

      relaxdata->origsolfeasible = relaxdata->storedfeasibility;
   }

   /** @todo solve master problem again */

   return SCIP_OKAY;
}


/** checks whether a variable should be added as an external branching candidate, if so it is added */
static
SCIP_RETCODE checkAndAddExternalBranchingCandidate(
   GCG*                  gcg,                /**< the GCG data structure */
   SCIP_VAR*             var                 /**< the variable to check whether to add as a branching candidate */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   assert(var != NULL);

   scip = GCGgetOrigprob(gcg);

   if( SCIPvarGetType(var) <= SCIP_VARTYPE_INTEGER && !SCIPisFeasIntegral(scip, SCIPgetRelaxSolVal(scip, var)) )
   {
      if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
      {
         SCIPdebugMessage("lblocal = %g, ublocal = %g\n", SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         SCIPdebugMessage("var = %s, vartype = %d, val = %g\n", SCIPvarGetName(var), SCIPvarGetType(var),
            SCIPgetRelaxSolVal(scip, var));
      }

      assert(!SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

      SCIP_CALL( SCIPaddExternBranchCand(scip, var, SCIPgetRelaxSolVal(scip, var) -
            SCIPfloor(scip, SCIPgetRelaxSolVal(scip, var)), SCIPgetRelaxSolVal(scip, var)) );
   }

   return SCIP_OKAY;
}


/* frees current currentorigsol, sets origsolfeasible to FALSE and clears external branching candidates */
static
SCIP_RETCODE freeCurrentOrigSol(
   GCG*  gcg,                                /**< the GCG data structure */
   SCIP_RELAXDATA* relaxdata                 /**< relaxdata of the relaxator */
   )
{
   SCIP* scip;

   scip = GCGgetOrigprob(gcg);
   relaxdata->origsolfeasible = FALSE;
   /* free previous solution and clear branching candidates */
   if( relaxdata->currentorigsol != NULL )
   {
      SCIPdebugMessage("Freeing previous solution origsol\n");
      SCIP_CALL( SCIPfreeSol(scip, &(relaxdata->currentorigsol)) );
   }

   if( SCIPgetStage(GCGgetMasterprob(relaxdata->gcg)) == SCIP_STAGE_SOLVING )
   {
      SCIPclearExternBranchCands(scip);
   }
   return SCIP_OKAY;
}


/** transforms the current solution of the master problem into the original problem's space
 *  and saves this solution as currentsol in the relaxator's data
 */
SCIP_RETCODE GCGrelaxUpdateCurrentSol(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP* masterprob;
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_VAR** origvars;
   int norigvars;
   SCIP_Bool stored;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);
   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   origvars = SCIPgetVars(scip);
   norigvars = SCIPgetNVars(scip);
   assert(origvars != NULL);

   /* retrieving the master problem */
   masterprob = GCGgetMasterprob(gcg);

   /* if the master problem has not been solved, don't try to update the solution */
   if( SCIPgetStage(masterprob) == SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( freeCurrentOrigSol(gcg, relaxdata) );
      return SCIP_OKAY;
   }

   if( SCIPgetStage(masterprob) == SCIP_STAGE_SOLVED || SCIPgetLPSolstat(masterprob) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_SOL* mastersol;
      SCIP_Longint currentnode;

      currentnode = SCIPgetCurrentNode(masterprob) == NULL ? -1 : SCIPnodeGetNumber(SCIPgetCurrentNode(masterprob));

      /* create new solution */
      if( SCIPgetStage(masterprob) == SCIP_STAGE_SOLVING )
      {
         SCIPdebugMessage("Masterproblem still solving, mastersol = NULL\n");
         mastersol = NULL;

         if( relaxdata->lastmasternode == currentnode && relaxdata->lastmasterlpiters >= SCIPgetNLPIterations(masterprob) )
         {
            SCIPdebugMessage("no new lp iterations\n");
            return SCIP_OKAY;
         }
      }
      else if( SCIPgetStage(masterprob) == SCIP_STAGE_SOLVED )
      {
         mastersol = SCIPgetBestSol(masterprob);
         if( mastersol == NULL )
         {
            SCIP_CALL( freeCurrentOrigSol(gcg, relaxdata) );
            SCIPdebugMessage("Masterproblem solved, no master sol present\n");
            return SCIP_OKAY;
         }
         SCIPdebugMessage("Masterproblem solved, mastersol = %p\n", (void*) mastersol);
      }
      else
      {
         SCIPdebugMessage("stage in master not solving and not solved!\n");
         return SCIP_OKAY;
      }

      /* free previous solution and clear branching candidates */
      SCIP_CALL( freeCurrentOrigSol(gcg, relaxdata) );

      relaxdata->lastmasterlpiters = SCIPgetNLPIterations(masterprob);
      relaxdata->lastmasternode = currentnode;

      if( !SCIPisInfinity(scip, SCIPgetSolOrigObj(masterprob, mastersol)) && GCGmasterIsSolValid(gcg, mastersol) )
      {
         int i;
         int j;
         SCIP_Bool violatesvarbnds;

         /* transform the master solution to the original variable space */
         SCIP_CALL( GCGtransformMastersolToOrigsol(gcg, mastersol, &(relaxdata->currentorigsol), FALSE, &violatesvarbnds) );
         assert(!violatesvarbnds || !GCGmasterIsSolValid(gcg, mastersol));

         /* store the solution as relaxation solution */
         SCIP_CALL( SCIPsetRelaxSolValsSol(scip, relax, relaxdata->currentorigsol, RELAX_INCLUDESLP) );
         assert(SCIPisEQ(scip, SCIPgetRelaxSolObj(scip), SCIPgetSolTransObj(scip, relaxdata->currentorigsol)));

         if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_BENDERS )
            SCIP_CALL( SCIPtrySol(scip, relaxdata->currentorigsol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );
         else
            SCIP_CALL( SCIPcheckSolOrig(scip, relaxdata->currentorigsol, &stored, FALSE, TRUE) );

         SCIPdebugMessage("updated current original LP solution, %s feasible in the original problem!\n",
            (stored ? "" : "not"));

         if( stored )
            relaxdata->origsolfeasible = TRUE;

         /* in the case of Benders decomposition, only the master variables can be added as branching candidates */
         if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_BENDERS )
         {
            SCIP_VAR** mastervars;
            SCIP_VAR** masterorigvars;
            int nmastervars;
            int nmasterorigvars;

            /* get variables of the master problem and their solution values */
            SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

            /* looping over all master variables to get the original variable for branching candidates */
            for( i = 0; i < nmastervars; i++ )
            {
               masterorigvars = GCGmasterVarGetOrigvars(mastervars[i]);
               nmasterorigvars = GCGmasterVarGetNOrigvars(mastervars[i]);

               for( j = 0; j < nmasterorigvars; j++ )
                  SCIP_CALL( checkAndAddExternalBranchingCandidate(gcg, masterorigvars[j]) );
            }
         }
         else
         {
            assert(GCGgetDecompositionMode(gcg) == GCG_DECMODE_DANTZIGWOLFE );
            /* store branching candidates */
            for( i = 0; i < norigvars; i++ )
               SCIP_CALL( checkAndAddExternalBranchingCandidate(gcg, origvars[i]) );
         }

         SCIPdebugMessage("updated relaxation branching candidates\n");
      }
   }
   else
   {
      SCIP_CALL( freeCurrentOrigSol(gcg, relaxdata) );
   }

   return SCIP_OKAY;
}


/** gets the total memory used after problem creation stage for all pricingproblems */
SCIP_Real GCGgetPricingprobsMemUsed(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   int p;
   SCIP_Real memused;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   memused = 0.0;

   /* @todo replace the computation by relaxdata->pricingprobsmemused if we can assure that the memory
    * used by the pricing problems is constant */

   /* compute memory that is used by all pricing problems */
   for( p = 0; p < relaxdata->npricingprobs; ++p )
   {
      memused += SCIPgetMemUsed(relaxdata->pricingprobs[p])/1048576.0;
   }

   return memused;
}

/** returns whether the relaxator has been initialized */
SCIP_Bool GCGrelaxIsInitialized(
   GCG*                  gcg                 /**< GCG data structure */
   )
{

   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->relaxisinitialized;
}

/** returns the average degeneracy */
SCIP_Real GCGgetDegeneracy(
   GCG*                  gcg                 /**< GCG data structure */
   )
{

   SCIP_Real degeneracy;
   SCIP* masterprob;
   SCIP* origprob;

   assert(gcg != NULL);

   origprob = GCGgetOrigprob(gcg);
   masterprob  = GCGgetMasterprob(gcg);

   degeneracy = 0.0;
   if( masterprob != NULL )
   {
      degeneracy = GCGmasterGetDegeneracy(gcg);
      if( SCIPisInfinity(masterprob, degeneracy) )
         degeneracy = SCIPinfinity(origprob);
   }
   return degeneracy;
}

/** return linking constraints for variables */
SCIP_CONS** GCGgetVarLinkingconss(
   GCG*                  gcg                 /**< GCG data structure */
  )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->varlinkconss;
}

/** return blocks of linking constraints for variables */
int* GCGgetVarLinkingconssBlock(
   GCG*                  gcg                 /**< GCG data structure */
  )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->varlinkconsblock;
}

/** return number of linking constraints for variables */
int GCGgetNVarLinkingconss(
   GCG*                  gcg                 /**< GCG data structure */
  )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->nvarlinkconss;
}

/** return number of linking variables */
int GCGgetNLinkingvars(
   GCG*                  gcg                 /**< GCG data structure */
  )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->nlinkingvars;
}

/** return number of variables directly transferred to the master problem */
int GCGgetNTransvars(
   GCG*                  gcg                 /**< GCG data structure */
  )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->ntransvars;
}

/** returns the relaxation solution from the Benders' decomposition */
SCIP_SOL* GCGgetBendersRelaxationSol(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;
   SCIP_BENDERS* benders;

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   benders = SCIPfindBenders(GCGgetMasterprob(relaxdata->gcg), "gcg");
   assert(benders != NULL);

   return GCGbendersGetRelaxSol(benders);
}

/** returns the decomposition mode */
GCG_DECMODE GCGgetDecompositionMode(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return (GCG_DECMODE)relaxdata->mode;
}

/** return root node time clock */
SCIP_CLOCK* GCGgetRootNodeTime(
   GCG*                 gcg               /**< GCG data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->rootnodetime;
}

SCIP_RETCODE GCGtransformProb(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   switch( SCIPgetStage(origprob) )
   {
   case SCIP_STAGE_INIT:
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPconshdlrDecompRepairConsNames(gcg) );
      SCIP_CALL( SCIPtransformProb(origprob) );
      break;

   case SCIP_STAGE_TRANSFORMED:
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "problem is already transformed\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }

   return SCIP_OKAY;
}

SCIP_RETCODE GCGpresolve(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   switch( SCIPgetStage(origprob) )
   {
   case SCIP_STAGE_INIT:
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( GCGtransformProb(gcg) );
      assert(SCIPgetStage(origprob) == SCIP_STAGE_TRANSFORMED);

      /*lint -fallthrough*/

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
      if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_DANTZIGWOLFE )
      {
         SCIP_CALL( GCGstashLimitSettings(gcg) );
      }
      SCIP_CALL( SCIPpresolve(origprob) );
      SCIP_CALL( GCGconshdlrDecompTranslateOrigPartialdecs(gcg) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "problem is already presolved\n");
      break;

   case SCIP_STAGE_SOLVED:
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "problem is already solved\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }

   return SCIP_OKAY;
}

SCIP_RETCODE GCGdetect(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RESULT result;
   SCIP* origprob = GCGgetOrigprob(gcg);

   switch( SCIPgetStage(origprob) )
   {
   case SCIP_STAGE_INIT:
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_DIALOG, NULL, "no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( GCGtransformProb(gcg) );
      assert(SCIPgetStage(origprob) == SCIP_STAGE_TRANSFORMED);

      /*lint -fallthrough*/

   case SCIP_STAGE_TRANSFORMED:
      if( GCGdetectionTookPlace(gcg, TRUE) )
      {
         SCIPverbMessage(origprob, SCIP_VERBLEVEL_DIALOG, NULL, "The detection for the original problem took place already.\n");
      }
      else
      {
         SCIPverbMessage(origprob, SCIP_VERBLEVEL_DIALOG, NULL, "starting detection\n");
         SCIP_CALL( GCGdetectStructure(gcg, &result) );
      }
      break;
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
      if( GCGdetectionTookPlace(gcg, FALSE) )
      {
         SCIPverbMessage(origprob, SCIP_VERBLEVEL_DIALOG, NULL, "The detection for the presolved problem took place already.\n");
      }
      else
      {
         SCIPverbMessage(origprob, SCIP_VERBLEVEL_DIALOG, NULL, "starting detection\n");
         SCIP_CALL( GCGdetectStructure(gcg, &result) );
      }
      break;
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }

   return SCIP_OKAY;
}

SCIP_RETCODE GCGsolve(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_RESULT result;
   int presolrounds;
   SCIP_Bool exit = FALSE;
   SCIP* origprob = GCGgetOrigprob(gcg);

   presolrounds = -1;

   assert(GCGconshdlrDecompCheckConsistency(gcg) );

   while( !exit )
   {
      switch( SCIPgetStage(origprob) )
      {
      case SCIP_STAGE_INIT:
         SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "No problem exists\n");
         exit = TRUE;
         break;

      case SCIP_STAGE_PROBLEM:
         SCIP_CALL( GCGtransformProb(gcg) );
         assert(SCIPgetStage(origprob) == SCIP_STAGE_TRANSFORMED);

         /*lint -fallthrough*/

      case SCIP_STAGE_TRANSFORMED:
      case SCIP_STAGE_PRESOLVING:
         if( GCGconshdlrDecompOrigPartialdecExists(gcg) )
         {
            SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "there is an original decomposition and problem is not presolved yet -> disable presolving and start optimizing (rerun with presolve command before detect command for detecting in presolved problem)  \n");
            SCIP_CALL( SCIPgetIntParam(origprob, "presolving/maxrounds", &presolrounds) );
            SCIP_CALL( SCIPsetIntParam(origprob, "presolving/maxrounds", 0) );
         }
         SCIP_CALL( GCGpresolve(gcg) );
         assert(SCIPgetStage(origprob) > SCIP_STAGE_PRESOLVING);

         break;

      case SCIP_STAGE_PRESOLVED:
         assert(GCGconshdlrDecompCheckConsistency(gcg) );

         if( !GCGdetectionTookPlace(gcg, TRUE) && !GCGdetectionTookPlace(gcg, FALSE) && GCGconshdlrDecompGetNFinishedPartialdecsTransformed(gcg) == 0 )
         {
            SCIP_CALL( GCGdetectStructure(gcg, &result) );
            if( result == SCIP_DIDNOTFIND )
            {
               GCG_DECOMP* bestdecomp;
               bestdecomp = GCGgetBestDecomp(gcg, TRUE);
               assert(bestdecomp == NULL && (GCGdetectionTookPlace(gcg, TRUE) || GCGdetectionTookPlace(gcg, FALSE)));
               GCGdecompFree(gcg, &bestdecomp);
               SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "No decomposition exists or could be detected. Solution process started with original problem...\n");
            }
         }
         else if( !GCGdetectionTookPlace(gcg, TRUE) && !GCGdetectionTookPlace(gcg, FALSE) && GCGconshdlrDecompGetNFinishedPartialdecsTransformed(gcg) > 0 )
         {
   #ifndef NDEBUG
            GCG_DECOMP* bestdecomp;
            bestdecomp = GCGgetBestDecomp(gcg, TRUE);
            assert(bestdecomp != NULL);
            GCGdecompFree(gcg, &bestdecomp);
   #endif
            SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "Preexisting decomposition found. Solution process started...\n");
         }
         else if( GCGconshdlrDecompGetNFinishedPartialdecsTransformed(gcg) == 0 )
         {
            SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "No decomposition exists or could be detected. Solution process started with original problem...\n");
         }
         assert(GCGconshdlrDecompCheckConsistency(gcg));
         assert(SCIPgetNConss(origprob) == SCIPgetNActiveConss(origprob));

         /*lint -fallthrough*/
      case SCIP_STAGE_SOLVING:
         if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_DANTZIGWOLFE && SCIPgetNNodes(origprob) == 0 )
         {
            SCIP_CALL( GCGstashLimitSettings(gcg) );
         }
         SCIP_CALL( SCIPsolve(origprob) );
         exit = TRUE;
         break;

      case SCIP_STAGE_SOLVED:
         SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, "Problem is already solved\n");
         exit = TRUE;
         break;

      case SCIP_STAGE_TRANSFORMING:
      case SCIP_STAGE_INITPRESOLVE:
      case SCIP_STAGE_EXITPRESOLVE:
      case SCIP_STAGE_INITSOLVE:
      case SCIP_STAGE_EXITSOLVE:
      case SCIP_STAGE_FREETRANS:
      case SCIP_STAGE_FREE:
      default:
         SCIPerrorMessage("invalid SCIP stage <%d>\n", SCIPgetStage(origprob));
         return SCIP_INVALIDCALL;
      }
   }

   if( presolrounds != -1 )
   {
      SCIP_CALL( SCIPsetIntParam(origprob, "presolving/maxrounds", presolrounds) );
   }

   return SCIP_OKAY;
}

SCIP_Real GCGgetDualbound(
   GCG*                 gcg               /**< GCG data structure */
   )
{
   SCIP* origprob;
   SCIP* masterprob;
   SCIP_Real dualbound;

   assert(gcg != NULL);

   origprob = GCGgetOrigprob(gcg);

   /* get master problem */
   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);

   dualbound = SCIPgetDualbound(origprob);

   /* @todo find a better way to do this */
   if( SCIPgetStage(masterprob) >= SCIP_STAGE_SOLVING )
   {
      SCIP_Real masterdualbound;

      masterdualbound = SCIPgetDualbound(masterprob);
      masterdualbound = SCIPretransformObj(origprob, masterdualbound);
      dualbound = MAX(dualbound, masterdualbound);
   }

   return dualbound;
}

SCIP_Real GCGgetPrimalbound(
   GCG*                 gcg               /**< GCG data structure */
   )
{
   SCIP* origprob;
   SCIP* masterprob;
   SCIP_Real primalbound;

   assert(gcg != NULL);

   origprob = GCGgetOrigprob(gcg);

   /* get master problem */
   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);

   primalbound = SCIPgetPrimalbound(origprob);

   /* @todo find a better way to do this */
   if( SCIPgetStage(masterprob) >= SCIP_STAGE_SOLVING && GCGmasterIsBestsolValid(gcg) )
   {
      SCIP_Real masterprimalbound;
      masterprimalbound = SCIPgetPrimalbound(masterprob);
      masterprimalbound = SCIPretransformObj(origprob, masterprimalbound);

      primalbound = MIN(primalbound, masterprimalbound);
   }

   return primalbound;
}

SCIP_Real GCGgetGap(
   GCG*                 gcg               /**< GCG data structure */
   )
{
   SCIP* origprob;
   SCIP_Real dualbound;
   SCIP_Real primalbound;
   SCIP_Real gap;

   assert(gcg != NULL);

   origprob = GCGgetOrigprob(gcg);
   primalbound = GCGgetPrimalbound(gcg);
   dualbound = GCGgetDualbound(gcg);

   /* this is the gap calculation from SCIPgetGap() */
   if( SCIPisEQ(origprob, primalbound, dualbound) )
      gap = 0.0;
   else if( SCIPisZero(origprob, dualbound )
      || SCIPisZero(origprob, primalbound)
      || SCIPisInfinity(origprob, REALABS(primalbound))
      || SCIPisInfinity(origprob, REALABS(dualbound))
      || primalbound * dualbound < 0.0 )
      gap = SCIPinfinity(origprob);
   else
   {
      SCIP_Real absdual = REALABS(dualbound);
      SCIP_Real absprimal = REALABS(primalbound);

      gap = REALABS((primalbound - dualbound)/MIN(absdual, absprimal));
   }

   return gap;
}

SCIP_RETCODE GCGinitializeMasterProblemSolve(
   GCG*                  gcg
   )
{
   SCIP_RELAX* relax;

   assert(gcg != NULL);

   relax = GCGgetRelax(gcg);
   assert(relax != NULL);
   assert(SCIPgetStage(GCGgetOrigprob(gcg)) >= SCIP_STAGE_TRANSFORMED);
   return initializeMasterProblemSolve(gcg, relax);
}

SCIP_RETCODE GCGstashLimitSettings(
   GCG*                 gcg
   )
{
   SCIP* scip;
   SCIP* masterprob;
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;

   scip = GCGgetOrigprob(gcg);
   masterprob = GCGgetMasterprob(gcg);
   relax = GCGgetRelax(gcg);
   assert(relax != NULL);
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if( relaxdata->limitsettingsstashed == FALSE )
   {
      relaxdata->limitsettingsstashed = TRUE;

      SCIP_CALL( SCIPgetLongintParam(scip, "limits/nodes", &relaxdata->stashednodelimit) );
      SCIP_CALL( SCIPgetLongintParam(scip, "limits/stallnodes", &relaxdata->stashedstallnodelimit) );
      SCIP_CALL( SCIPgetRealParam(scip, "limits/gap", &relaxdata->stashedgaplimit) );
      SCIP_CALL( SCIPgetIntParam(scip, "limits/solutions", &relaxdata->stashedsollimit) );
      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &relaxdata->stashedtimelimit) );

      SCIPresetParam(scip, "limits/nodes");
      SCIPresetParam(scip, "limits/stallnodes");
      SCIPresetParam(scip, "limits/gap");
      SCIPresetParam(scip, "limits/solutions");
      SCIPresetParam(scip, "limits/time");
      SCIPresetParam(masterprob, "limits/gap");
      SCIPresetParam(masterprob, "limits/time");
   }

   return SCIP_OKAY;
}

SCIP_RETCODE GCGrestoreLimitSettings(
   GCG*                 gcg
   )
{
   SCIP* scip;
   SCIP* masterprob;
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;

   scip = GCGgetOrigprob(gcg);
   masterprob = GCGgetMasterprob(gcg);
   relax = GCGgetRelax(gcg);
   assert(relax != NULL);
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if( relaxdata->limitsettingsstashed == TRUE )
   {
      relaxdata->limitsettingsstashed = FALSE;

      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", relaxdata->stashednodelimit) );
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/stallnodes", relaxdata->stashedstallnodelimit) );
      SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", relaxdata->stashedgaplimit) );
      SCIP_CALL( SCIPsetIntParam(scip, "limits/solutions", relaxdata->stashedsollimit) );
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", relaxdata->stashedtimelimit) );

      SCIP_CALL( setMasterLimits(gcg, masterprob, relaxdata->stashedtimelimit, relaxdata->stashedgaplimit) );
   }

   return SCIP_OKAY;
}

#ifdef _OPENMP
GCG_LOCKS* GCGgetLocks(
   GCG*                  gcg                /**< the GCG data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax = GCGgetRelax(gcg);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->locks;
}
#endif

/** returns the GCG data structure */
GCG_EXPORT
GCG* GCGrelaxGetGcg(
   SCIP*                 origprob            /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   assert(origprob != NULL);

   relax = SCIPfindRelax(origprob, RELAX_NAME);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->gcg;
}
