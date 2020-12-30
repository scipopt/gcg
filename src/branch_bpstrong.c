/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2020 Operations Research, RWTH Aachen University       */
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

/**@file    branch_bpstrong.c
 * @ingroup BRANCHINGRULES
 * @brief   generic branch and price strong branching as described in
 *          Pecin, D., Pessoa, A., Poggi, M., Uchoa, E. Improved branch-cut-and-price for capacitated vehicle routing.
 *          In: Math. Prog. Comp. 9:61-100. Springer (2017).
 * @author  Oliver Gaul
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>

#include "branch_bpstrong.h"
#include "type_branchgcg.h"
#include "gcg.h"
#include "branch_orig.h"

#include <string.h>

#include "gcg.h"
#include "branch_relpsprob.h"
#include "cons_integralorig.h"
#include "cons_masterbranch.h"
#include "cons_origbranch.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "type_branchgcg.h"

#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"

#define BRANCHRULE_NAME "bpstrong"                              /**< name of branching rule */
#define BRANCHRULE_DESC "strong branching for branch-and-price" /**< short description of branching rule */
#define BRANCHRULE_PRIORITY -99999                              /**< priority of this branching rule */
#define BRANCHRULE_MAXDEPTH 0                                  /**< maximal depth level of the branching rule */
#define BRANCHRULE_MAXBOUNDDIST 1.0                             /**< maximal relative distance from current node's      \
                                                                     dual bound to primal bound compared to best node's \
                                                                     dual bound for applying branching */
#define DEFAULT_ENFORCEBYCONS  FALSE
#define DEFAULT_MOSTFRAC       FALSE
#define DEFAULT_USEPSEUDO      TRUE
#define DEFAULT_USEPSSTRONG    FALSE

#define DEFAULT_USESTRONG      FALSE
#define DEFAULT_STRONGLITE     FALSE
#define DEFAULT_STRONGTRAIN    FALSE
#define DEFAULT_IMMEDIATEINF   TRUE

#define DEFAULT_REEVALAGE           1
#define DEFAULT_MINCOLGENCANDS      4
#define DEFAULT_MINPHASE0OUTCANDS   10
#define DEFAULT_MAXPHASE0OUTCANDS   50
#define DEFAULT_PHASE1GAPWEIGHT     0.25
#define DEFAULT_MINPHASE1OUTCANDS   3
#define DEFAULT_MAXPHASE1OUTCANDS   20
#define DEFAULT_PHASE2GAPWEIGHT     1
#define DEFAULT_HISTWEIGHT          0.5   


/** branching rule data */
struct SCIP_BranchruleData
{
   int                   lastcand;              /**< last evaluated candidate of last branching rule execution */
   int                   nvars;                 /**< the number of vars currently in the hashmap */
   int                   maxvars;               /**< the maximal number of vars that were in the hashmap at the same time */
   SCIP_HASHMAP*         varhashmap;            /**< hashmap mapping variables to their last result in strong branching */
   SCIP_Real             *score;                /**< the variables' last scores */
   int                   *uniqueblockflags;     /**< flags assigned by assignUniqueBlockFlags() */
   SCIP_Real             *strongbranchscore;    /**< the variables' last score from strong branching with column generation */
   SCIP_Bool             *sbscoreisrecent;      /**< was the score saved in strongbranchscore computed in a parent of the current node *
                                                 *   where all node on the path to the parent were created for domainreduction due to infeasibility? */
   int                   *lastevalnode;         /**< the last node at which the variables were evaluated */

   SCIP_Bool             mostfrac;              /**< should branching be performed on the most fractional variable instead of the first variable? */
   SCIP_Bool             usepseudocosts;        /**< should pseudocosts be used to determine the variable on which the branching is performed? */

   SCIP_Bool             usestronglite;         /**< should strong branching use column generation during variable evaluation? */
   SCIP_Bool             usestrongtrain;        /**< should strong branching run as precise as possible (to generate more valuable training data)? */
   SCIP_Bool             immediateinf;          /**< should infeasibility detected during strong branching be handled immediately, or only if the variable is selected? */
   int                   reevalage;             /**< how many times can bounds be changed due to infeasibility during strong branching until an already evaluated variable needs to be reevaluated? */
   int                   mincolgencands;        /**< minimum number of variables for phase 2 to be executed, otherwise the best candidate from phase 1 will be chosen */
   int                   minphasezerooutcands;  /**< minimum number of output candidates from phase 0 */
   int                   maxphasezerooutcands;  /**< maximum number of output candidates from phase 0 */
   SCIP_Real             phaseonegapweight;     /**< how much impact should the nodegap have on the number of precisely evaluated candidates in phase 1? */
   int                   minphaseoneoutcands;   /**< minimum number of output candidates from phase 1 */
   int                   maxphaseoneoutcands;   /**< maximum number of output candidates from phase 1 */
   SCIP_Real             phasetwogapweight;     /**< how much impact should the nodegap have on the number of precisely evaluated candidates in phase 2? */
   SCIP_Real             histweight;            /**< how many candidates should be chosen based on historical strong branching scores as opposed to current heuristic scores in phase 0 (e.g. 0.5 = 50%)? */
};

/* needed for compare_function (for now)*/
SCIP_BRANCHRULEDATA* this_branchruledata;

/* calculates the number of needed candidates based on the min and max number of candidates as well as the node gap */ 
static
int calculateNCands(
   SCIP* scip,
   int min,             /**< minimum number of candidates */
   int max,             /**< maximum number of candidates */
   SCIP_Real nodegap,   /**< node gap in current focus node */
   SCIP_Real gapweight  /**< how much influence the node gap should have */
)
{
   int dif = max-min;

   assert( min>=1 );

   return min + (int) SCIPceil(scip, MIN(dif, dif * nodegap * gapweight + dif * (1-gapweight)));
}

/* return  1: integer variables belonging to a unique block with fractional value
 * return  0: variables that belong to no block but were directly transferred to the 
 *            master problem and which have a fractional value in the current solution
 * return -1: neither
 */
static
int assignUniqueBlockFlags(
   SCIP* scip,
   SCIP_VAR* branchcand
)
{
   assert(GCGvarIsOriginal(branchcand));

   for (int iter = 0; iter <= 1; iter++)
   {
      /* continue if variable belongs to a block in second iteration*/
      if (iter == 0)
      {
         /* variable belongs to no block */
         if (GCGvarGetBlock(branchcand) == -1)
            continue;

         /* block is not unique (non-linking variables) */
         if (!GCGoriginalVarIsLinking(branchcand) && GCGgetNIdenticalBlocks(scip, GCGvarGetBlock(branchcand)) != 1)
            continue;

         /* check that blocks of linking variable are unique */
         if (GCGoriginalVarIsLinking(branchcand))
         {
            int nvarblocks;
            int *varblocks;
            SCIP_Bool unique;
            int j;

            nvarblocks = GCGlinkingVarGetNBlocks(branchcand);
            SCIP_CALL(SCIPallocBufferArray(scip, &varblocks, nvarblocks));
            SCIP_CALL(GCGlinkingVarGetBlocks(branchcand, nvarblocks, varblocks));

            unique = TRUE;
            for (j = 0; j < nvarblocks; ++j)
               if (GCGgetNIdenticalBlocks(scip, varblocks[j]) != 1)
                  unique = FALSE;

            SCIPfreeBufferArray(scip, &varblocks);

            if (!unique)
               continue;
         }
         /* candidate is valid in first iteration */
         return 1;
         
      }
      else /* iter == 1 */
      {
         if (GCGvarGetBlock(branchcand) != -1)
            return -1;

         /* candidate is valid in second iteration */
         return 0;
      }
   }
   return -1;
}

/** adds branching candidates to branchruledata to collect infos about it */
static
SCIP_RETCODE addBranchcandsToData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_VAR**            branchcands,        /**< branching candidates */
   int                   npriobranchcands    /**< number of priority branching candidates */
   )
{

   SCIP_BRANCHRULEDATA* branchruledata;
   int i;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( branchruledata->nvars == 0 )
   { 
      assert(branchruledata->varhashmap != NULL);

      /* create arrays */
      branchruledata->maxvars = SCIPcalcMemGrowSize(scip, npriobranchcands);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->score, branchruledata->maxvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->uniqueblockflags, branchruledata->maxvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->strongbranchscore, branchruledata->maxvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->sbscoreisrecent, branchruledata->maxvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->lastevalnode, branchruledata->maxvars) );
      branchruledata->nvars = npriobranchcands;

      /* store each variable in hashmap and initialize array entries */
      for( i = 0; i < npriobranchcands; ++i )
      {
         SCIP_CALL( SCIPhashmapInsertInt(branchruledata->varhashmap, branchcands[i], i) );
         branchruledata->score[i] = -1;
         branchruledata->strongbranchscore[i] = -1;
         branchruledata->sbscoreisrecent[i] = FALSE;
         branchruledata->lastevalnode[i] = -1;
         branchruledata->uniqueblockflags[i] = -2;
      }
   }
   else  /* possibly new variables need to be added */
   {

      /* if var is not in hashmap, insert it */
      for( i = 0; i < npriobranchcands; i++ )
      {
         SCIP_VAR* var;
         int nvars;

         var = branchcands[i];
         assert(var != NULL);
         nvars = branchruledata->nvars;

         /* if variable is not in hashmap insert it, initialize its array entries, and increase array sizes */
         if( !SCIPhashmapExists(branchruledata->varhashmap, var) )
         {
            int newsize = SCIPcalcMemGrowSize(scip, nvars + 1);
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &branchruledata->score, branchruledata->maxvars,
               newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &branchruledata->strongbranchscore, branchruledata->maxvars,
               newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &branchruledata->sbscoreisrecent, branchruledata->maxvars,
               newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &branchruledata->lastevalnode, branchruledata->maxvars,
               newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &branchruledata->uniqueblockflags, branchruledata->maxvars,
               newsize) );
            branchruledata->maxvars = newsize;

            SCIP_CALL( SCIPhashmapInsertInt(branchruledata->varhashmap, var, nvars) );
            branchruledata->score[nvars] = -1;
            branchruledata->strongbranchscore[nvars] = -1;
            branchruledata->sbscoreisrecent[nvars] = FALSE;
            branchruledata->lastevalnode[nvars] = -1;
            branchruledata->uniqueblockflags[nvars] = -2;

            assert(SCIPhashmapExists(branchruledata->varhashmap, var)
               && SCIPhashmapGetImageInt(branchruledata->varhashmap, var) == nvars); /*lint !e507*/

            ++(branchruledata->nvars);
         }
      }
   }

   return SCIP_OKAY;
}

/* compare two indices corresponding to entries in branchruledata->score/uniqueblockflags */
static int score_compare_function(const void *index1, const void *index2)
{
   return (this_branchruledata->score[*(int *)index1] > this_branchruledata->score[*(int *)index2] ||
           this_branchruledata->uniqueblockflags[*(int *)index1] > this_branchruledata->uniqueblockflags[*(int *)index2]) ? -1 : 1;
}

/* compare two indices corresponding to entries in branchruledata->strongbranchscore/uniqueblockflags */
static int hist_compare_function(const void *index1, const void *index2)
{
   return (this_branchruledata->strongbranchscore[*(int *)index1] > this_branchruledata->strongbranchscore[*(int *)index2] ||
           this_branchruledata->uniqueblockflags[*(int *)index1] > this_branchruledata->uniqueblockflags[*(int *)index2]) ? -1 : 1;
}

/* compare two indices based on descending numerical order */
static int geq_compare_function(const void *index1, const void *index2)
{
   return *(int *)index1 < *(int *)index2 ? -1 : 1;
}

/* executes strong branching on one variable, with or without pricing */
static
SCIP_RETCODE executeStrongBranching(
    SCIP *scip,           /* SCIP data structure */
    SCIP_BRANCHRULE*      branchrule,         /* pointer to the branching rule */
    SCIP_VAR *branchvar,        /* variable to get strong branching values for */
    SCIP_Real solval,     /* value of the variable in the current solution */
    SCIP_Bool pricing,    /* should pricing be applied? */
    int maxpricingrounds, /* maximal number of pricing rounds, -1 for no limit */
    SCIP_Real *up,      /* stores dual bound after branching column up */
    SCIP_Real *down,         /* stores dual bound after branching column down */
    SCIP_Bool *upvalid, /* stores whether the upbranch was solved properly */
    SCIP_Bool *downvalid,    /* stores whether the downbranch was solved properly */
    SCIP_Bool *upinf,   /*stores whether the upbranch is infeasible */
    SCIP_Bool *downinf  /*stores whether the downbranch is infeasible */
)
{
   /* get bound values */
   char name[SCIP_MAXSTRLEN];
   SCIP* masterscip;

   SCIP_Real downub;
   SCIP_Real uplb;

   SCIP_Bool cutoff;
   SCIP_Bool lperror;
   SCIP_Bool lpsolved;

   downub = SCIP_INVALID;
   uplb = SCIP_INVALID;
   *downvalid = FALSE;
   *upvalid = FALSE;
   *downinf = FALSE;
   *upinf = FALSE;

    /* get master problem */
   masterscip = GCGgetMasterprob(scip);
   assert(masterscip != NULL);

   assert(scip != NULL);

   /* create child nodes with x <= floor(x'), and x >= ceil(x') */
   downub = SCIPfeasFloor(scip, solval);
   uplb = downub + 1.0;

   //SCIPdebugMessage("Probing on var %s with value %g in current solution\n", SCIPvarGetName(branchvar), solval);

   /* probe for each child node */
   for( int cnode = 0; cnode <= 1; cnode++ )
   {
      if ((cnode == 0 && downub != SCIP_INVALID) || (cnode == 1 && uplb != SCIP_INVALID))
      {
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchvar),
                           cnode == 0? "<=":">=", cnode == 0? downub : uplb);

         /* start probing */
         SCIP_CALL( GCGrelaxStartProbing(scip, NULL) );
         SCIP_CALL( GCGrelaxNewProbingnodeOrig(scip) );

         cutoff = FALSE;
         lperror = FALSE;
         lpsolved = FALSE;

         if( cnode == 0 )
         {
            SCIP_CALL( SCIPchgVarUbProbing(scip, branchvar, downub) );
         }
         else
         {
            SCIP_CALL( SCIPchgVarLbProbing(scip, branchvar, uplb) );
         }

         /* propagate the new b&b-node */
         SCIP_CALL(SCIPpropagateProbing(scip, -1, &cutoff, NULL));

         /* solve the LP with or without pricing */
         if( !cutoff )
         {
            SCIP_CALL( GCGrelaxNewProbingnodeMaster(scip) );
            if (pricing)
            {
               SCIP_CALL( GCGrelaxPerformProbingWithPricing(scip, -1, NULL, NULL,
                        cnode == 0? down : up, &lpsolved, &lperror, &cutoff) );
            }
            else
            {
               SCIP_CALL( GCGrelaxPerformProbing(scip, -1, NULL,
                        cnode == 0? down : up, &lpsolved, &lperror, &cutoff) );
            }
         }

         if( cnode == 0 )
         {
            *downvalid = lpsolved;
            *downinf = cutoff && pricing;
         }
         else
         {
            *upvalid = lpsolved;
            *upinf = cutoff && pricing;
         }

         //SCIPdebugMessage("probing results in cutoff/lpsolved/lpobj: %s / %s / %g\n",
         //      cutoff?"cutoff":"no cutoff", lpsolved?"lpsolved":"lp not solved", cnode == 0? *down : *up);
         SCIP_CALL( GCGrelaxEndProbing(scip) );
      }
   }
   return SCIP_OKAY;
}

/* Returns true iff the the second node is a k-successor of the first to the first number corresponding node
 * (i.e. iff there are at most k edges between them)
 */
static
SCIP_Bool isKAncestor(
    SCIP* scip, 
    int ancestornodenr,          /**< number of the supposed ancestor */
    SCIP_NODE *successornode,    /**< the supposed successor */
    int k                        /**< maximal allowed distance between the nodes */
)
{
   SCIP_NODE* curnode;
   curnode = successornode;

   for( int i = 0; i<=k && SCIPnodeGetNumber(curnode) >= ancestornodenr; i++ )
   {
      if( SCIPnodeGetNumber(curnode) == ancestornodenr )
         return TRUE;

      if( SCIPnodeGetNumber(curnode) == 1)
         break;

      curnode = SCIPnodeGetParent(curnode);
   }

   return FALSE;
}

/* Evaluates the given variable based on a score function of choice. Higher scores are given to better
 * variables.
 */
static
SCIP_Real score_function(
    SCIP *scip,
    SCIP_BRANCHRULE*      branchrule,         /* pointer to the branching rule */
    SCIP_VAR *var,            /* var to be scored */
    SCIP_Real solval,         /* the var's current solution value */
    SCIP_Bool useheuristic,   /* should heuristics be used instead of strong branching? */
    SCIP_Bool usehistorical,  /* should historical data from phase 2 be used as heuristic? */
    SCIP_Bool usecolgen,      /* should column generation be used during strong branching? */
    SCIP_Real *score,         /* stores the computed score */
    SCIP_Bool *upinf,         /* stores whether the upbranch is infeasible */
    SCIP_Bool *downinf        /* stores whether the downbranch is infeasible */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* define score functions and calculate score for all variables for sorting dependent on used heuristic */
   // phase 0
   if( useheuristic)
   {
      if( usehistorical )
      {
         int hashindex;

         assert(SCIPhashmapExists(branchruledata->varhashmap, var));
         hashindex = SCIPhashmapGetImageInt(branchruledata->varhashmap, var);

         return branchruledata->strongbranchscore[hashindex];
      }
      else if( branchruledata->usepseudocosts )
      {
         *score = SCIPgetVarPseudocostScore(scip, var, solval);
      }
      else /* no parameter for fractional variable selection? */
      {
         if( !branchruledata->mostfrac )
            return 1;
            
         *score = solval - SCIPfloor(scip, solval);
         *score = MIN(*score, 1.0 - *score);
      }
   }
   else
   //phase 1 & 2
   {
      SCIP* masterscip;

      int hashindex;
      int currentnodenr;
      
      SCIP_Real down;
      SCIP_Real up;
      SCIP_Real downgain;
      SCIP_Real upgain;
      SCIP_Bool upvalid;
      SCIP_Bool downvalid;
      SCIP_Real lpobjval;

      /* get master problem */
      masterscip = GCGgetMasterprob(scip);
      assert(masterscip != NULL);

      assert(SCIPhashmapExists(branchruledata->varhashmap, var));
      hashindex = SCIPhashmapGetImageInt(branchruledata->varhashmap, var);
      currentnodenr = SCIPnodeGetNumber(SCIPgetFocusNode(scip));
   
      if( !usecolgen 
          || !branchruledata->sbscoreisrecent[hashindex]
          || !isKAncestor(scip, branchruledata->lastevalnode[hashindex], SCIPgetFocusNode(scip), branchruledata->reevalage) )
      {
         up = -SCIPinfinity(scip);
         down = -SCIPinfinity(scip);

         lpobjval = SCIPgetLPObjval(masterscip);

         // usecolgen is True for phase 1 and False for phase 2
         SCIP_CALL( executeStrongBranching(scip, branchrule, var, solval, usecolgen, -1, &up, &down, &upvalid, &downvalid, upinf, downinf) );

         //TODO handle this better
         down = downvalid? down : upvalid? up : 0;
         up = upvalid? up : down;

         downgain = down - lpobjval;
         upgain = up - lpobjval;

         *score = SCIPgetBranchScore(scip, var, downgain, upgain);

         if( usecolgen && upvalid && downvalid && !*upinf && !*downinf )
         {
            branchruledata->strongbranchscore[hashindex] = *score;
            branchruledata->sbscoreisrecent[hashindex] = TRUE;
            branchruledata->lastevalnode[hashindex] = currentnodenr;
         }
         //SCIPdebugMessage("Variable %s has downgain %f and upgain %f\n", SCIPvarGetName(var), downgain, upgain);
      }
      else
      {
         *score = branchruledata->strongbranchscore[hashindex];
      }
   }

   return SCIP_OKAY;
}

/** branching method for relaxation solutions */
static
SCIP_RETCODE branchExtern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< pointer to the branching rule */
   SCIP_VAR**            branchvar,          /**< pointer to store the selected variable pointer */
   SCIP_Bool*            bestupinf,          /**< pointer to store whether strong branching detected infeasibility in the upbranch */
   SCIP_Bool*            bestdowninf,        /**< pointer to store whether strong branching detected infeasibility in the downbranch */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching call */
   )
{
   SCIP* masterscip;
   SCIP_BRANCHRULEDATA* branchruledata;

   /* branching candidates */
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   int nbranchcands;
   int npriobranchcands;

   SCIP_VAR* checkvar;

   int ncands;
   int nneededcands;

   SCIP_Real nodegap;
   SCIP_Real upperbound;
   SCIP_Real nodelowerbound;

   int hashindex;

   /* values for choosing the variable to branch on */
   SCIP_Real solval;

   SCIP_Real maxscore;
   SCIP_Real score;

   /* infeasibility results during strong branching */
   SCIP_Bool upinf;
   SCIP_Bool downinf;

   /* storing best candidates */
   int *indices;
   int nvalidcands;

   /* storing best candidates based on historical strong branching scores */
   int *histindices;
   int nvalidhistcands;
   int nneededhistcands;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPisRelaxSolValid(scip));

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   /* get master problem */
   masterscip = GCGgetMasterprob(scip);
   assert(masterscip != NULL);

   /* get the branching candidates */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &branchcands, &branchcandssol, NULL, &nbranchcands,
         &npriobranchcands, NULL, NULL, NULL) );

   *branchvar = NULL;
   solval = 0.0;

   maxscore = -1.0;

   upinf = FALSE;
   downinf = FALSE;
   *bestupinf = FALSE;
   *bestdowninf = FALSE;

   upperbound = SCIPgetUpperbound(scip);
   nodelowerbound = SCIPnodeGetLowerbound( SCIPgetFocusNode(scip) );
   nodegap = ((upperbound >= 0) == (nodelowerbound >= 0))? 
             MIN(ABS((upperbound-nodelowerbound)/MIN(ABS(upperbound), ABS(nodelowerbound))), 1) : 1;
   assert(0<=nodegap && nodegap<=1);

   /* number of candidates we evaluate precisely should be based on the likely relevance of this branching decision via the nodegap */
   nneededcands =  calculateNCands(scip, branchruledata->minphasezerooutcands, branchruledata->maxphasezerooutcands, nodegap, branchruledata->phaseonegapweight);
   nvalidcands = 0;
   nvalidhistcands = 0;

   /* insert branchcands into hashmap */
   SCIP_CALL( addBranchcandsToData(scip, branchrule, branchcands, npriobranchcands) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &indices, npriobranchcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &histindices, npriobranchcands) );
   indices[0] = 0;

   /* iter = 0: integer variables belonging to a unique block with fractional value,
    * iter = 1: we did not find enough variables to branch on so far, so we look for integer variables that belong to no block
    * but were directly transferred to the master problem and which have a fractional value in the current solution
    */
   for( int iter = 0; iter <= 1 && nvalidcands < nneededcands; iter++ )
   {
      for( int i = 0; i < npriobranchcands; i++ )
      {
         hashindex = SCIPhashmapGetImageInt(branchruledata->varhashmap, branchcands[i]);     

         if (iter == 0)
         {  
            if( branchruledata->uniqueblockflags[hashindex] < -1 )
            {
               branchruledata->uniqueblockflags[hashindex] = assignUniqueBlockFlags(scip, branchcands[i]);
            }
            
            if( branchruledata->uniqueblockflags[hashindex] == 1 )
            {
               indices[nvalidcands] = i;
               nvalidcands++;

               if( branchruledata->strongbranchscore[hashindex] != -1)
               {
                  histindices[nvalidhistcands] = i;
                  nvalidhistcands++;
               }
            }
         }
         else /* iter == 1 */
         {
            if( branchruledata->uniqueblockflags[hashindex] == 0 )
            {
               indices[nvalidcands] = i;
               nvalidcands++;
               if( branchruledata->strongbranchscore[hashindex] != -1)
               {
                  histindices[nvalidhistcands] = i;
                  nvalidhistcands++;
               }
            }
         }
      }
   }

   /* the number of candidates we select based on historical strong branching scores needs to depend on the number of
   * candidates for which we have historical scores, otherwise some candidates would be selected simply because they
   * have been scored before
   */
   nneededhistcands = SCIPfloor(scip, MIN((SCIP_Real)nvalidhistcands/(SCIP_Real)(nvalidcands+nvalidhistcands), branchruledata->histweight) * nvalidcands);

   qsort(histindices, nvalidhistcands, sizeof(int), hist_compare_function);
   qsort(histindices, nneededhistcands, sizeof(int), geq_compare_function);

   /* if we are performing strong branching, go through the three phases:
    * - phase 0: select a first selection (50 to 10, based on |T S(v)|) of candidates based on some traditional variable selection
    *            heuristic, some (half) of the variables are new, and some are selected based on previous calls
    * - phase 1: find 20 to 3 (based on |T S(v)|) best candidates by evaluating the Master LP, w/o column and cut generation
    * - phase 2: select the best of the candidates from phase 1 by solving the Master LP with column and cut generation.
    * else, select the candidate that perfroms best on the given heuristic (i.e. phase 0 with only one output candidate)
    */
   for( int phase = 0; phase<=2; phase++ )
   {
      switch( phase )
      {
         case 0:
            ncands = nvalidcands;
            break;

         case 1:
            nneededcands = calculateNCands(scip, branchruledata->minphaseoneoutcands, branchruledata->maxphaseoneoutcands,
                                           nodegap, branchruledata->phasetwogapweight);

            /* skip phase 2 if we are in lite mode,
             * or if the number of available candidates is lower than the min amount for phase 2
             */
            if( branchruledata->usestronglite 
               || nneededcands < branchruledata->mincolgencands
               || ncands < branchruledata->mincolgencands )
               nneededcands = 1;

            break;

         case 2:
            nneededcands = 1;
            break;

      }

      if( nneededcands >= ncands )
      {
         continue;
      }

      /* compute scores */
      for( int i = 0, c=branchruledata->lastcand; i < ncands; i++, c++ )
      {
         c = c % ncands;
         
         assert(GCGvarIsOriginal(branchcands[indices[c]]));
         //checkvar = branchcands[indices[c]];

         /* select the variable as new best candidate (if it is) if we look for only one candidate,
          * or remember its score if we look for multiple
          */
         SCIP_CALL( score_function(scip, branchrule, branchcands[indices[c]], branchcandssol[indices[c]], phase == 0, FALSE, phase == 2 && !branchruledata->usestronglite, &score, &upinf, &downinf) );
         
         /* variable pointers sometimes change during probing in strong branching */
         //if( checkvar != branchcands[indices[c]] )
         //{
         //   if(phase==2)
         //      SCIP_CALL( SCIPgetExternBranchCands(scip, &branchcands, &branchcandssol, NULL, NULL,
         //         NULL, NULL, NULL, NULL) );
         
         //}

         /* handle infeasibility detected during strong branching */
         if( phase == 2 && !branchruledata->usestronglite && branchruledata->immediateinf && (upinf || downinf) )
         {
            if( upinf && downinf )
            {
               for( int k=0; k<branchruledata->maxvars; k++ )
               {
                  branchruledata->sbscoreisrecent[k] = FALSE;
               }
               *result = SCIP_CUTOFF;

               SCIPfreeBufferArray(scip, &indices);
               SCIPfreeBufferArray(scip, &histindices);

               *bestupinf = TRUE;
               *bestdowninf = TRUE;

               SCIPdebugMessage("Original branching rule detected current node to be infeasible!\n");
               return SCIP_OKAY;
            }

            branchruledata->lastcand = c;
            indices[0] = indices[c];
            *bestupinf = upinf;
            *bestdowninf = downinf;
            break;
         }

         if( nneededcands == 1 )
         {
            if( score > maxscore )
            {
               indices[0] = indices[c];
               maxscore = score;
               *bestupinf = upinf;
               *bestdowninf = downinf;
            }
         }
         else
         {
            hashindex = SCIPhashmapGetImageInt(branchruledata->varhashmap, branchcands[c]);
            branchruledata->score[hashindex] = score;
         } 
      }

      if( nneededcands > 1 )
      {
         qsort(indices, ncands, sizeof(int), score_compare_function);
         ncands = MIN(ncands, nneededcands);

         if( phase == 0 && nneededhistcands )
         {
            /* swap out the worst performing "new" candidates with the best performing historical candidates */
            int *indicescopy;
            int pos;

            SCIP_CALL( SCIPallocBufferArray(scip, &indicescopy, ncands) );
            pos = nneededhistcands;

            for( int i = 0; i<ncands; i++ )
            {
               indicescopy[i] = indices[i];
            }

            for( int i = 0; i<nneededhistcands; i++ )
            {
               indices[i] = histindices[i];
            }

            /* concatenate the two arrays, while avoiding duplicates */
            for( int i = 0; i<ncands && pos<=ncands; i++ )
            {
               for( int j = 0; j<=nneededhistcands; j++ )
               {
                  if( j == nneededhistcands )
                  {
                     indices[pos] = indicescopy[i];
                     pos++;
                  }
                  else if( indices[j] == indicescopy[i] )
                     break;
               }
               
            }
            
            SCIPfreeBufferArray(scip, &indicescopy);
         }

      }
      else
      {
         break;
      }     
   }

   *branchvar = branchcands[indices[0]];
   solval = SCIPgetRelaxSolVal(scip, branchcands[indices[0]]);

   /* free memory */
   SCIPfreeBufferArray(scip, &indices);
   SCIPfreeBufferArray(scip, &histindices);

   if( *branchvar == NULL )
   {
      SCIPdebugMessage("Strong branching could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(*branchvar != NULL);
   //TODO
   assert(!(*bestupinf && *bestdowninf));
   
   SCIPdebugMessage("Strong branching selected variable %s with solval %f%s\n", SCIPvarGetName(*branchvar), solval, (*bestupinf || *bestdowninf)? ", which is infeasible in one direction" : "");
   
   if( !*bestupinf && !*bestdowninf )
   {
      for( int i=0; i<branchruledata->maxvars; i++ )
      {
         branchruledata->sbscoreisrecent[i] = FALSE;
      }
   }

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/*
 * Callback methods
 */
#define branchDeactiveMasterBPStrong NULL
#define branchPropMasterBPStrong NULL
#define branchActiveMasterBPStrong NULL
#define branchMasterSolvedBPStrong NULL
#define branchDataDeleteBPStrong NULL

#define branchExeclpBPStrong NULL
#define branchExecextBPStrong NULL
#define branchExecpsBPStrong NULL

static
SCIP_DECL_BRANCHFREE(branchFreeBPStrong)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_HASHMAP* varhashmap;

   branchruledata = SCIPbranchruleGetData(branchrule);
   varhashmap = branchruledata->varhashmap;

   if( branchruledata->varhashmap != NULL )
   {
      SCIPhashmapFree(&varhashmap);
   }

   SCIPfreeBlockMemoryArray(scip, &branchruledata->score, branchruledata->maxvars);
   SCIPfreeBlockMemoryArray(scip, &branchruledata->strongbranchscore, branchruledata->maxvars);
   SCIPfreeBlockMemoryArray(scip, &branchruledata->uniqueblockflags, branchruledata->maxvars);
   SCIPfreeBlockMemoryArray(scip, &branchruledata->lastevalnode, branchruledata->maxvars);

   SCIPfreeBlockMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitBPStrong)
{
   SCIP* origprob;
   SCIP_BRANCHRULEDATA* branchruledata;

   origprob = GCGmasterGetOrigprob(scip);
   assert(branchrule != NULL);
   assert(origprob != NULL);

   SCIPdebugMessage("Init BPStrong branching rule\n");

   SCIP_CALL( GCGrelaxIncludeBranchrule( origprob, branchrule, branchActiveMasterBPStrong,
         branchDeactiveMasterBPStrong, branchPropMasterBPStrong, branchMasterSolvedBPStrong, branchDataDeleteBPStrong) );

   branchruledata = SCIPbranchruleGetData(branchrule);
   branchruledata->lastcand = 0;
   branchruledata->nvars = 0;
   branchruledata->maxvars = 0;

   this_branchruledata = branchruledata;
   SCIP_CALL( SCIPhashmapCreate(&(branchruledata->varhashmap), SCIPblkmem(scip), SCIPgetNVars(scip)) );

   return SCIP_OKAY;
}

/** creates the b&p strong-branching branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleBPStrong(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP* origscip;
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   SCIPdebugMessage("Include BPStrong branching rule\n");

   /* get original problem */
   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   /* alloc branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );
   
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
            BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );
   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitBPStrong) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeBPStrong) );

   /* add branching rule parameters */
   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/bp_strong/stronglite",
         "should strong branching use column generation during variable evaluation?",
         &branchruledata->usestronglite, FALSE, DEFAULT_STRONGLITE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/bp_strong/strongtraining",
         "should strong branching run as precise as possible (to generate more valuable training data)?",
         &branchruledata->usestrongtrain, FALSE, DEFAULT_STRONGTRAIN, NULL, NULL) );
      
   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/bp_strong/immediateinf",
         "should infeasibility detected during strong branching be handled immediately, or only if the variable is selected?",
         &branchruledata->immediateinf, FALSE, DEFAULT_IMMEDIATEINF, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "branching/bp_strong/reevalage",
         "how many times can bounds be changed due to infeasibility during strong branching until an already evaluated variable needs to be reevaluated?",
         &branchruledata->reevalage, FALSE, DEFAULT_REEVALAGE, 0, 100, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "branching/bp_strong/mincolgencands",
         "minimum number of variables for phase 2 to be executed, otherwise the best candidate from phase 1 will be chosen",
         &branchruledata->mincolgencands, FALSE, DEFAULT_MINCOLGENCANDS, 0, 100000, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "branching/bp_strong/minphase0outcands",
         "minimum number of output candidates from phase 0",
         &branchruledata->minphasezerooutcands, FALSE, DEFAULT_MINPHASE0OUTCANDS, 1, 100000, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "branching/bp_strong/maxphase0outcands",
         "maximum number of output candidates from phase 0",
         &branchruledata->maxphasezerooutcands, FALSE, DEFAULT_MAXPHASE0OUTCANDS, 1, 100000, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origscip, "branching/bp_strong/phase1gapweight",
         "how much impact should the nodegap have on the number of precisely evaluated candidates in phase 1?",
         &branchruledata->phaseonegapweight, FALSE, DEFAULT_PHASE1GAPWEIGHT, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "branching/bp_strong/minphase1outcands",
         "minimum number of output candidates from phase 1",
         &branchruledata->minphaseoneoutcands, FALSE, DEFAULT_MINPHASE1OUTCANDS, 1, 100000, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "branching/bp_strong/maxphase1outcands",
         "maximum number of output candidates from phase 1",
         &branchruledata->maxphaseoneoutcands, FALSE, DEFAULT_MAXPHASE1OUTCANDS, 1, 100000, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origscip, "branching/bp_strong/phase2gapweight",
         "how much impact should the nodegap have on the number of precisely evaluated candidates in phase 2?",
         &branchruledata->phasetwogapweight, FALSE, DEFAULT_PHASE2GAPWEIGHT, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origscip, "branching/bp_strong/histweight",
         "how many candidates should be chosen based on historical strong branching scores as opposed to current heuristic scores in phase 0 (e.g. 0.5 = 50%)?",
         &branchruledata->histweight, FALSE, DEFAULT_HISTWEIGHT, 0, 1, NULL, NULL) );


   /* notify cons_integralorig about the branching rule */
   SCIP_CALL( GCGconsIntegralorigAddBranchrule(scip, branchrule) );

   return SCIP_OKAY;
}

SCIP_RETCODE
GCGbranchSelectCandidateStrongBranchingOrig(
   SCIP* scip,                      /**< SCIP data structure */
   SCIP_BRANCHRULE *origbranchrule, /**< pointer storing original branching rule */
   SCIP_VAR **branchvar,            /**< pointer to store output var pointer */
   SCIP_Bool *upinf,                /**< pointer to store whether strong branching detected infeasibility in the upbranch */
   SCIP_Bool *downinf,              /**< pointer to store whether strong branching detected infeasibility in the downbranch */
   SCIP_RESULT *result              /**< pointer to store result */
)
{
   SCIP_BRANCHRULEDATA *origbranchruledata;
   SCIP_BRANCHRULEDATA *branchruledata;
   SCIP_BRANCHRULE *branchrule;
   SCIP* masterscip;

   masterscip = GCGgetMasterprob(scip);
   branchrule = SCIPfindBranchrule(masterscip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   origbranchruledata = SCIPbranchruleGetData(origbranchrule);

   branchruledata->usepseudocosts = origbranchruledata->usepseudocosts;
   branchruledata->mostfrac = origbranchruledata->mostfrac;

   branchExtern(scip, branchrule, branchvar, upinf, downinf, result);

   return SCIP_OKAY;
}