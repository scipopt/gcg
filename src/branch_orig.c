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

/**@file   branch_orig.c
 * @brief  branching rule for the original problem in GCG
 * @author Gerald Gamrath
 * @author Marcel Schmickerath
 * @author Christian Puchert
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>
#include <string.h>

#include "branch_orig.h"
#include "gcg.h"
#include "branch_relpsprob.h"
#include "cons_integralorig.h"
#include "cons_masterbranch.h"
#include "cons_origbranch.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "type_branchgcg.h"

#include "scip/cons_linear.h"


#define BRANCHRULE_NAME          "orig"
#define BRANCHRULE_DESC          "branching for the original program in generic column generation"
#define BRANCHRULE_PRIORITY      100
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define DEFAULT_ENFORCEBYCONS  FALSE
#define DEFAULT_MOSTFRAC       FALSE
#define DEFAULT_USEPSEUDO      TRUE
#define DEFAULT_USEPSSTRONG    FALSE

#define DEFAULT_USESTRONG      FALSE
#define DEFAULT_STRONGLITE     FALSE
#define DEFAULT_STRONGTRAIN    FALSE
#define DEFAULT_IMMEDIATEINF   TRUE

#define DEFAULT_REEVALAGE      1
#define DEFAULT_MINCOLGENCANDS 4
#define DEFAULT_PHASE0OUTCANDS 40
#define DEFAULT_PHASE1OUTCANDS 20
#define DEFAULT_GAPWEIGHT      1   


/** branching rule data */
struct SCIP_BranchruleData
{
   int                   lastcand;           /**< last evaluated candidate of last branching rule execution */
   int                   nvars;              /**< the number of vars currently in the hashmap */
   int                   maxvars;            /**< the maximal number of vars that were in the hashmap at the same time */
   SCIP_HASHMAP*         varhashmap;         /**< hashmap mapping variables to their last result in strong branching */
   SCIP_Real             *score;             /**< the variables' last scores */
   int                   *uniqueblockflags;  /**< flags assigned by assignUniqueBlockFlags() */
   SCIP_Real             *strongbranchscore; /**< the variables' last score from strong branching with column generation */
   SCIP_Bool             *sbscoreisrecent;   /**< was the score saved in strongbranchscore computed in a parent of the current node *
                                              *   where all node on the path to the parent were created for domainreduction due to infeasibility? */
   int                   *lastevalnode;      /**< the last node at which the variables were evaluated */

   SCIP_Bool             enforcebycons;      /**< should bounds on variables be enforced by constraints(TRUE) or by bounds(FALSE) */
   SCIP_Bool             mostfrac;           /**< should branching be performed on the most fractional variable instead of the first variable? */
   SCIP_Bool             usepseudocosts;     /**< should pseudocosts be used to determine the variable on which the branching is performed? */
   SCIP_Bool             usepsstrong;        /**< should strong branching with propagation be used to determine the variable on which the branching is performed? */
   SCIP_Bool             usestrong;          /**< should strong branching be used to determine the variable on which the branching is performed? */

   SCIP_Bool             usestronglite;      /**< should strong branching use column generation during variable evaluation? */
   SCIP_Bool             usestrongtrain;     /**< should strong branching run as precise as possible (to generate more valuable training data)? */
   SCIP_Bool             immediateinf;       /**< should infeasibility detected during strong branching be handled immediately, or only if the variable is selected? */
   int                   reevalage;          /**< how many times can bounds be changed due to infeasibility during strong branching until an already evaluated variable needs to be reevaluated? */
   int                   mincolgencands;     /**< minimum number of variables for phase 2 to be executed, otherwise the best candidate from phase 1 will be chosen */
   int                   phasezerooutcands;  /**< maximum number of output candidates from phase 0 */
   int                   phaseoneoutcands;   /**< maximum number of output candidates from phase 1 */
   SCIP_Real             gapweight;          /**< how much impact should the nodegap have on the number of precisely evaluated candidates? */
};

/** branching data for branching decisions */
struct GCG_BranchData
{
   SCIP_VAR*             origvar;            /**< original variable on which the branching is done */
   GCG_BOUNDTYPE         boundtype;          /**< type of the new bound of original variable */
   SCIP_Real             newbound;           /**< new lower/upper bound of the original variable */
   SCIP_Real             oldbound;           /**< old lower/upper bound of the pricing variable */
   SCIP_Real             oldvalue;           /**< old value of the original variable */
   SCIP_Real             olddualbound;       /**< dual bound before the branching was performed */
   SCIP_CONS*            cons;               /**< constraint that enforces the branching restriction in the original
                                              *   problem, or NULL if this is done by variable bounds */
};

/* needed for compare_function (for now)*/
SCIP_BRANCHRULEDATA* this_branchruledata;

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
   int                   npriobranchcands        /**< number of priority branching candidates */
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
         SCIP_CALL( SCIPhashmapInsert(branchruledata->varhashmap, branchcands[i], (void*) (size_t)i) );
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

            SCIP_CALL( SCIPhashmapInsert(branchruledata->varhashmap, var, (void*) (size_t)nvars) );
            branchruledata->score[nvars] = -1;
            branchruledata->strongbranchscore[nvars] = -1;
            branchruledata->sbscoreisrecent[nvars] = FALSE;
            branchruledata->lastevalnode[nvars] = -1;
            branchruledata->uniqueblockflags[nvars] = -2;

            assert(SCIPhashmapExists(branchruledata->varhashmap, var)
               && (int)(size_t) SCIPhashmapGetImage(branchruledata->varhashmap, var) == nvars); /*lint !e507*/

            ++(branchruledata->nvars);
         }
      }
   }

   return SCIP_OKAY;
}

/** branches on a integer variable x
 *  if solution value x' is fractional, two child nodes will be created
 *  (x <= floor(x'), x >= ceil(x')),
 *  if solution value is integral, the bounds of x are finite, then two child nodes will be created
 *  (x <= x", x >= x"+1 with x" = floor((lb + ub)/2)),
 *  otherwise (up to) three child nodes will be created
 *  (x <= x'-1, x == x', x >= x'+1)
 *  if solution value is equal to one of the bounds and the other bound is infinite, only two child nodes
 *  will be created (the third one would be infeasible anyway)
 */
static
SCIP_RETCODE branchVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< pointer to the original variable branching rule */
   SCIP_VAR*             branchvar,          /**< variable to branch on */
   SCIP_Real             solval,             /**< value of the variable in the current solution */
   SCIP_Bool             upinf,              /**< have we seen during strong branching that the upbranch is infeasible? */
   SCIP_Bool             downinf             /**< have we seen during strong branching that the downbranch is infeasible? */
   )
{
   /* data for b&b child creation */
   SCIP* masterscip;
   SCIP_Real downub;
   SCIP_Real fixval;
   SCIP_Real uplb;

   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(branchvar != NULL);

   /* get master problem */
   masterscip = GCGgetMasterprob(scip);
   assert(masterscip != NULL);

   downub = SCIP_INVALID;
   fixval = SCIP_INVALID;
   uplb = SCIP_INVALID;

   if( SCIPisFeasIntegral(scip, solval) )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbLocal(branchvar);
      ub = SCIPvarGetUbLocal(branchvar);

      /* if there was no explicit value given for branching, the variable has a finite domain and the current LP/pseudo
       * solution is one of the bounds, we branch in the center of the domain */
      if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
      {
         SCIP_Real center;

         /* create child nodes with x <= x", and x >= x"+1 with x" = floor((lb + ub)/2);
          * if x" is integral, make the interval smaller in the child in which the current solution x'
          * is still feasible
          */
         center = (ub + lb) / 2.0;
         if( solval <= center )
         {
            downub = SCIPfeasFloor(scip, center);
            uplb = downub + 1.0;
         }
         else
         {
            uplb = SCIPfeasCeil(scip, center);
            downub = uplb - 1.0;
         }
      }
      else
      {
         /* create child nodes with x <= x'-1, x = x', and x >= x'+1 */
         assert(SCIPisEQ(scip, SCIPfeasCeil(scip, solval), SCIPfeasFloor(scip, solval)));

         fixval = solval;

         /* create child node with x <= x'-1, if this would be feasible */
         if( SCIPisFeasGE(scip, fixval-1.0, lb) )
            downub = fixval - 1.0;

         /* create child node with x >= x'+1, if this would be feasible */
         if( SCIPisFeasLE(scip, fixval+1.0, ub) )
            uplb = fixval + 1.0;
      }
      SCIPdebugMessage("integral branch on variable <%s> with value %g, priority %d (current lower bound: %g)\n",
         SCIPvarGetName(branchvar), solval, SCIPvarGetBranchPriority(branchvar), SCIPgetLocalLowerbound(GCGgetMasterprob(scip)));
   }
   else
   {
      /* create child nodes with x <= floor(x'), and x >= ceil(x') */
      downub = SCIPfeasFloor(scip, solval);
      uplb = downub + 1.0;
      assert( SCIPisEQ(scip, SCIPfeasCeil(scip, solval), uplb) );
      //SCIPdebugMessage("fractional branch on variable <%s> with value %g, root value %g, priority %d (current lower bound: %g)\n",
         //SCIPvarGetName(branchvar), solval, SCIPvarGetRootSol(branchvar), SCIPvarGetBranchPriority(branchvar), SCIPgetLocalLowerbound(GCGgetMasterprob(scip)));
   }

   //SCIPdebugMessage("Branching on var %s with value %g in current solution\n", SCIPvarGetName(branchvar), solval);


   if( uplb != SCIP_INVALID && !upinf ) /*lint !e777*/
   {
      SCIP_CONS* cons;
      SCIP_NODE* child;
      SCIP_CONS** origbranchconss;
      GCG_BRANCHDATA* branchdata;
      char name[SCIP_MAXSTRLEN];

      int norigbranchconss;
      int maxorigbranchconss;

      origbranchconss = NULL;
      norigbranchconss = 0;
      maxorigbranchconss = 0;

      /* create child node x >= uplb */
      SCIP_CALL( SCIPcreateChild(masterscip, &child, 0.0, SCIPgetLocalTransEstimate(masterscip)) );

      SCIP_CALL( SCIPallocBlockMemory(scip, &branchdata) );

      branchdata->origvar = branchvar;
      branchdata->oldvalue = solval;
      branchdata->olddualbound = SCIPgetLocalLowerbound(masterscip);
      branchdata->boundtype = GCG_BOUNDTYPE_LOWER;
      branchdata->newbound = uplb;
      branchdata->oldbound = SCIPvarGetLbLocal(branchvar);

      SCIPdebugMessage(" -> creating child: <%s> >= %g\n",
         SCIPvarGetName(branchvar), uplb);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchdata->origvar),
         ">=", branchdata->newbound);

      if( branchruledata->enforcebycons )
      {
         /* enforce new bounds by linear constraints */
         SCIP_CONS* consup;

         SCIPdebugMessage("enforced by cons\n");

         norigbranchconss = 1;
         maxorigbranchconss = SCIPcalcMemGrowSize(scip, 1);
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &origbranchconss, maxorigbranchconss) );

         /* create corresponding constraints */
         SCIP_CALL( SCIPcreateConsLinear(scip, &consup, name, 0, NULL, NULL,
            SCIPceil(scip, solval), SCIPinfinity(scip),
            TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );

         SCIP_CALL( SCIPaddCoefLinear(scip, consup, branchvar, 1.0) );

         origbranchconss[0] = consup;
         branchdata->cons = consup;
      }
      else
         branchdata->cons = NULL;

      /* create and add the masterbranch constraint */
      SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &cons, name, child,
         GCGconsMasterbranchGetActiveCons(masterscip), branchrule, branchdata, origbranchconss, norigbranchconss,
         maxorigbranchconss) );
      SCIP_CALL( SCIPaddConsNode(masterscip, child, cons, NULL) );
   }

   if( downub != SCIP_INVALID && !downinf ) /*lint !e777*/
   {
      SCIP_CONS* cons;
      SCIP_NODE* child;
      SCIP_CONS** origbranchconss;
      GCG_BRANCHDATA* branchdata;
      char name[SCIP_MAXSTRLEN];

      int norigbranchconss;
      int maxorigbranchconss;

      origbranchconss = NULL;
      norigbranchconss = 0;
      maxorigbranchconss = 0;

      /* create child node x <= downub */
      SCIP_CALL( SCIPcreateChild(masterscip, &child, 0.0, SCIPgetLocalTransEstimate(masterscip)) );

      SCIP_CALL( SCIPallocBlockMemory(scip, &branchdata) );

      branchdata->origvar = branchvar;
      branchdata->oldvalue = solval;
      branchdata->olddualbound = SCIPgetLocalLowerbound(masterscip);
      branchdata->boundtype = GCG_BOUNDTYPE_UPPER;
      branchdata->newbound = downub;
      branchdata->oldbound = SCIPvarGetUbLocal(branchvar);

      SCIPdebugMessage(" -> creating child: <%s> <= %g\n",
         SCIPvarGetName(branchvar), downub);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchdata->origvar),
         "<=", branchdata->newbound);

      /* enforce branching decision by a constraint rather than by bound changes */
      if( branchruledata->enforcebycons )
      {
         /* enforce new bounds by linear constraints */
         SCIP_CONS* consdown;

         norigbranchconss = 1;
         maxorigbranchconss = SCIPcalcMemGrowSize(scip, 1);
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &origbranchconss, maxorigbranchconss) );

         /* create corresponding constraints */
         SCIP_CALL( SCIPcreateConsLinear(scip, &consdown, name, 0, NULL, NULL,
            -1.0 * SCIPinfinity(scip), SCIPfloor(scip, solval),
            TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, consdown, branchvar, 1.0) );

         origbranchconss[0] = consdown;
         branchdata->cons = consdown;
      }
      else
         branchdata->cons = NULL;

      /* create and add the masterbranch constraint */
      SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &cons, name, child,
         GCGconsMasterbranchGetActiveCons(masterscip), branchrule, branchdata, origbranchconss, norigbranchconss,
         maxorigbranchconss) );
      SCIP_CALL( SCIPaddConsNode(masterscip, child, cons, NULL) );
   }

   if( fixval != SCIP_INVALID ) /*lint !e777*/
   {
      SCIP_CONS* cons;
      SCIP_NODE* child;
      SCIP_CONS** origbranchconss;
      GCG_BRANCHDATA* branchdata;
      char name[SCIP_MAXSTRLEN];

      int norigbranchconss;
      int maxorigbranchconss;

      origbranchconss = NULL;
      norigbranchconss = 0;
      maxorigbranchconss = 0;

      /* create child node x = fixval */
      SCIP_CALL( SCIPcreateChild(masterscip, &child, 0.0, SCIPgetLocalTransEstimate(masterscip)) );

      SCIP_CALL( SCIPallocBlockMemory(scip, &branchdata) );

      branchdata->origvar = branchvar;
      branchdata->oldvalue = solval;
      branchdata->olddualbound = SCIPgetLocalLowerbound(masterscip);
      branchdata->boundtype = GCG_BOUNDTYPE_FIXED;
      branchdata->newbound = fixval;
      branchdata->oldbound = SCIPvarGetUbLocal(branchvar);

      SCIPdebugMessage(" -> creating child: <%s> == %g\n",
         SCIPvarGetName(branchvar), fixval);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchdata->origvar),
         "==", branchdata->newbound);

      /* enforce branching decision by a constraint rather than by bound changes */
      if( branchruledata->enforcebycons )
      {
         /* enforce new bounds by linear constraints */
         SCIP_CONS* consfix;

         norigbranchconss = 1;
         maxorigbranchconss = SCIPcalcMemGrowSize(scip, 1);
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &origbranchconss, maxorigbranchconss) );

         /* create corresponding constraints */
         SCIP_CALL( SCIPcreateConsLinear(scip, &consfix, name, 0, NULL, NULL,
            fixval, fixval, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, consfix, branchvar, 1.0) );

         origbranchconss[0] = consfix;
         branchdata->cons = consfix;
      }
      else
         branchdata->cons = NULL;

      /* create and add the masterbranch constraint */
      SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &cons, name, child,
         GCGconsMasterbranchGetActiveCons(masterscip), branchrule, branchdata, origbranchconss, norigbranchconss,
         maxorigbranchconss) );
      SCIP_CALL( SCIPaddConsNode(masterscip, child, cons, NULL) );
   }

   return SCIP_OKAY;
}

/* compare two indices corresponding to entries in branchruledata->score/uniqueblockflags */
static int compare_function(const void *index1, const void *index2)
{
   return (this_branchruledata->score[*(int *)index1] > this_branchruledata->score[*(int *)index2] ||
           this_branchruledata->uniqueblockflags[*(int *)index1] > this_branchruledata->uniqueblockflags[*(int *)index2]) ? -1 : 1;
}


/* executes strong branching on one variable, with or without pricing */
static SCIP_RETCODE executeStrongBranching(
    SCIP *scip,           /* SCIP data structure */
    SCIP_BRANCHRULE*      branchrule,         /* pointer to the original variable branching rule */
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
static SCIP_Real score_function(
    SCIP *scip,
    SCIP_BRANCHRULE*      branchrule,         /* pointer to the original variable branching rule */
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
         hashindex = (int)(size_t) SCIPhashmapGetImage(branchruledata->varhashmap, var);

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
      hashindex = (int)(size_t) SCIPhashmapGetImage(branchruledata->varhashmap, var);
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

#if FALSE
/** branching method for relaxation solutions */
static
SCIP_RETCODE branchExtern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< pointer to the original variable branching rule */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching call */
   )
{
   SCIP* masterscip;
   int i;

   /* parameter data */
   SCIP_Bool mostfrac;
   SCIP_Bool usepseudocosts;
   SCIP_Bool usepsstrong;

   /* branching candidates */
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   int nbranchcands;
   int npriobranchcands;

   /* values for choosing the variable to branch on */
   SCIP_VAR* branchvar;
   SCIP_Real solval;
   SCIP_Real maxfrac;
   SCIP_Real frac;
   SCIP_Real maxpsscore;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPisRelaxSolValid(scip));

   *result = SCIP_DIDNOTRUN;

   /* get master problem */
   masterscip = GCGgetMasterprob(scip);
   assert(masterscip != NULL);

   /* get values of parameters */
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/mostfrac", &mostfrac) );
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/usepseudocosts", &usepseudocosts) );
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/usepsstrong", &usepsstrong) );

   /* get the branching candidates */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &branchcands, &branchcandssol, NULL, &nbranchcands,
         &npriobranchcands, NULL, NULL, NULL) );

   branchvar = NULL;
   solval = 0.0;

   maxfrac = 0.0;
   maxpsscore = -1.0;

   if( usepsstrong )
   {
      SCIP_CALL( SCIPgetRelpsprobBranchVar(masterscip, branchcands, branchcandssol, npriobranchcands,
            npriobranchcands, result, &branchvar) );
      assert(branchvar != NULL || *result == SCIP_CUTOFF);
      assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF);

      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;

      solval = SCIPgetRelaxSolVal(scip, branchvar);
   }

   /* branch on an integer variable belonging to a unique block with fractional value */
   if( branchvar == NULL )
      for( i = 0; i < npriobranchcands; i++ )
      {
         
         assert(GCGvarIsOriginal(branchcands[i]));

         /* variable belongs to no block */
         if( GCGvarGetBlock(branchcands[i]) == -1 )
            continue;

         /* block is not unique (non-linking variables) */
         if( !GCGoriginalVarIsLinking(branchcands[i]) && GCGgetNIdenticalBlocks(scip, GCGvarGetBlock(branchcands[i])) != 1 )
            continue;

         /* check that blocks of linking variable are unique */
         if( GCGoriginalVarIsLinking(branchcands[i]) )
         {
            int nvarblocks;
            int* varblocks;
            SCIP_Bool unique;
            int j;

            nvarblocks = GCGlinkingVarGetNBlocks(branchcands[i]);
            SCIP_CALL( SCIPallocBufferArray(scip, &varblocks, nvarblocks) );
            SCIP_CALL( GCGlinkingVarGetBlocks(branchcands[i], nvarblocks, varblocks) );

            unique = TRUE;
            for( j = 0; j < nvarblocks; ++j )
               if( GCGgetNIdenticalBlocks(scip, varblocks[j]) != 1 )
                  unique = FALSE;

            SCIPfreeBufferArray(scip, &varblocks);

            if( !unique )
               continue;
         }
         SCIPdebugMessage("Looking at variable %s\n", SCIPvarGetName(branchcands[i]));
         /* use pseudocost variable selection rule */
         if( usepseudocosts )
         {
            /* select the variable, if its pseudocost are higher than the ones of the currently saved variable */
            if( SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]) > maxpsscore )
            {
               branchvar = branchcands[i];
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               maxpsscore = SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]);
            }
         }
         /* use most fractional variable selection rule */
         else
         {
            /* compute the fractionality */
            frac = branchcandssol[i] - SCIPfloor(scip, branchcandssol[i]);
            frac = MIN(frac, 1.0 - frac);
            assert(frac > 0);

            /* fractionality is higher than that of the current highest fractionality */
            if( frac >= maxfrac )
            {
               SCIPdebugMessage("Var %s has fractional value in current solution: %f\n", SCIPvarGetName(branchcands[i]), branchcandssol[i]);
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               branchvar = branchcands[i];
               /* if we do not look for the most fractional variable, but for the first fractional variable,
                * we can stop here since we found a variable to branch on */
               if( !mostfrac )
                  break;
            }
         }
      }

   /* we did not find a variable to branch on so far, so we look for an integer variable that belongs to no block
    * but was directly transferred to the master problem and which has fractional value in the current solution */
   if( branchvar == NULL )
   {
      for( i = 0; i < npriobranchcands; i++ )
      {
         assert(GCGvarIsOriginal(branchcands[i]));

         /* continue if variable belongs to a block */
         if( GCGvarGetBlock(branchcands[i]) != -1 )
            continue;

         SCIPdebugMessage("Looking at variable %s\n", SCIPvarGetName(branchcands[i]));

         /* use pseudocost variable selection rule */
         if( usepseudocosts )
         {
            if( SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]) > maxpsscore )
            {
               branchvar = branchcands[i];
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               maxpsscore = SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]);
            }
         }
         /* use most fractional variable selection rule */
         else
         {
            /* compute fractionality */
            frac = branchcandssol[i] - SCIPfloor(scip, branchcandssol[i]);
            frac = MIN(frac, 1.0 - frac);
            assert(frac > 0);

            if( frac >= maxfrac )
            {
               SCIPdebugMessage("Var %s has fractional value in current solution: %f\n",
                  SCIPvarGetName(branchcands[i]), branchcandssol[i]);
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               branchvar = branchcands[i];
               /* if we do not look for the most fractional variable, but for the first fractional variable,
                * we stop here since we found a variable to branch on */
               if( !mostfrac )
                  break;
            }
         }
      }
   }

   if( branchvar == NULL )
   {
      SCIPdebugMessage("Original branching rule could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(branchvar != NULL);

   SCIPdebugMessage("Original branching rule selected variable %s with solval %f\n", SCIPvarGetName(branchvar), solval);
   SCIP_CALL( branchVar(scip, branchrule, branchvar, solval, FALSE, FALSE) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}
#else
/** branching method for relaxation solutions */
static
SCIP_RETCODE branchExtern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< pointer to the original variable branching rule */
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

   int ncands;
   int nneededcands;

   int c;

   /* values for choosing the variable to branch on */
   SCIP_VAR* branchvar;
   SCIP_Real solval;

   SCIP_Real maxscore;
   SCIP_Real score;

   /* infeasibility results during strong branching */
   SCIP_Bool upinf;
   SCIP_Bool downinf;
   SCIP_Bool bestupinf;
   SCIP_Bool bestdowninf;

   int *indices;
   int nvalidcands;

   int hashindex;

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

   branchvar = NULL;
   solval = 0.0;

   maxscore = -1.0;

   upinf = FALSE;
   downinf = FALSE;
   bestupinf = FALSE;
   bestdowninf = FALSE;

   if( branchruledata->usepsstrong )
   {
      SCIP_CALL( SCIPgetRelpsprobBranchVar(masterscip, branchcands, branchcandssol, npriobranchcands,
            npriobranchcands, result, &branchvar) );
      assert(branchvar != NULL || *result == SCIP_CUTOFF);
      assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF);

      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;

      solval = SCIPgetRelaxSolVal(scip, branchvar);
   }

   //TODO
   nneededcands = branchruledata->usestrong? branchruledata->phasezerooutcands : 1;
   nvalidcands = 0;

   SCIPdebugMessage("Current Nodenr: %lld\n", SCIPnodeGetNumber(SCIPgetFocusNode(scip)));

   /* insert branchcands into hashmap */
   SCIP_CALL( addBranchcandsToData(scip, branchrule, branchcands, npriobranchcands) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &indices, branchruledata->usestrong? npriobranchcands : 1) );
   indices[0] = 0;

   /* iter = 0: integer variables belonging to a unique block with fractional value,
    * iter = 1: we did not find enough variables to branch on so far, so we look for integer variables that belong to no block
    * but were directly transferred to the master problem and which have a fractional value in the current solution
    */
   for( int iter = 0; iter <= 1 && nvalidcands < nneededcands; iter++ )
   {
      for( int i = 0; i < npriobranchcands; i++ )
      {
         hashindex = (int)(size_t) SCIPhashmapGetImage(branchruledata->varhashmap, branchcands[i]);     

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
            }
         }
         else /* iter == 1 */
         {
            if( branchruledata->uniqueblockflags[hashindex] == 0 )
            {
               indices[nvalidcands] = i;
               nvalidcands++;
            }
         }
      }
   }

   /* if we are performing strong branching, go through the three phases:
    * - phase 0: select a first selection (50 to 10, based on |T S(v)|) of candidates based on some traditional variable selection
    *            heuristic, some (half) of the variables are new, and some are selected based on previous calls
    * - phase 1: find 20 to 3 (based on |T S(v)|) best candidates by evaluating the Master LP, w/o column and cut generation
    * - phase 2: select the best of the candidates from phase 1 by solving the Master LP with column and cut generation.
    * else, select the candidate that perfroms best on the given heuristic (i.e. phase 0 with only one output candidate)
    */
   for( int phase = 0; phase<=0 || (branchruledata->usestrong && phase<=2); phase++ )
   {
      switch( phase )
      {
         case 0:
            ncands = nvalidcands;
            break;

         case 1:
            //TODO
            nneededcands = branchruledata->phaseoneoutcands;

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

         /* select the variable as new best candidate (if it is) if we look for only one candidate,
          * or add remember its score if we look for multiple
          */
         SCIP_CALL( score_function(scip, branchrule, branchcands[indices[c]], branchcandssol[indices[c]], phase == 0, FALSE, phase == 2 && !branchruledata->usestronglite, &score, &upinf, &downinf) );
         
         /* handle infeasibility detected during strong branching */
         if( phase == 2 && !branchruledata->usestronglite && branchruledata->immediateinf && (upinf || downinf) )
         {
            if( upinf && downinf )
            {
               //TODO actually handle this further up...
               for( int i=0; i<branchruledata->maxvars; i++ )
               {
                  branchruledata->sbscoreisrecent[i] = FALSE;
               }
               *result = SCIP_CUTOFF;
               SCIPdebugMessage("Original branching rule detected current node to be infeasible!\n");
               return SCIP_OKAY;
            }

            branchruledata->lastcand = c;
            indices[0] = indices[c];
            bestupinf = upinf;
            bestdowninf = downinf;
            *result = SCIP_REDUCEDDOM;
            break;
         }

         if( nneededcands == 1 )
         {
            if( score > maxscore )
            {
               indices[0] = indices[c];
               maxscore = score;
               bestupinf = upinf;
               bestdowninf = downinf;
               /* if we do not look for the most fractional variable, but for the first fractional variable,
               * we can stop here since we found a variable to branch on */
               if( !branchruledata->mostfrac && !branchruledata->usepseudocosts && !branchruledata->usestrong )
                  break;
            }
         }
         else
         {
            hashindex = (int)(size_t) SCIPhashmapGetImage(branchruledata->varhashmap, branchcands[c]);
            branchruledata->score[hashindex] = score;
         }
         SCIPdebugMessage("Looked at variable %s with current score: %f\n",
                                 SCIPvarGetName(branchcands[indices[c]]), score);   
      }
      if( nneededcands > 1 )
      {
         qsort(indices, ncands, sizeof(int), compare_function);
         ncands = MIN(ncands, nneededcands);
      }
      else
      {
         break;
      }     
   }

   branchvar = branchcands[indices[0]];
   solval = SCIPgetRelaxSolVal(scip, branchcands[indices[0]]);

   /* free memory */
   SCIPfreeBufferArray(scip, &indices);

   if( branchvar == NULL )
   {
      SCIPdebugMessage("Original branching rule could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(branchvar != NULL);
   assert(!(bestupinf && bestdowninf));
   
   SCIPdebugMessage("Original branching rule selected variable %s with solval %f%s\n", SCIPvarGetName(branchvar), solval, (bestupinf || bestdowninf)? ", which is infeasible in one direction" : "");
   
   if( !bestupinf && !bestdowninf )
   {
      for( int i=0; i<branchruledata->maxvars; i++ )
      {
         branchruledata->sbscoreisrecent[i] = FALSE;
      }
   }

   SCIP_CALL( branchVar(scip, branchrule, branchvar, solval, bestupinf, bestdowninf) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}
#endif

/*
 * Callback methods for enforcing branching constraints
 */

#define branchDeactiveMasterOrig NULL
#define branchPropMasterOrig NULL

/** callback activation method */
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterOrig)
{
   SCIP* origscip;
   SCIP_CONS* mastercons;

   assert(scip != NULL);
   assert(branchdata != NULL);

   /* branching restrictions are enforced by variable bounds, this is done automatically, so we can abort here */
   if( branchdata->cons == NULL )
      return SCIP_OKAY;

   assert(branchdata->origvar != NULL);

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   SCIPdebugMessage("branchActiveMasterOrig: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == GCG_BOUNDTYPE_LOWER ? ">=" : branchdata->boundtype == GCG_BOUNDTYPE_UPPER ? "<=" : "==" ), branchdata->newbound);

   /* transform constraint to the master variable space */
   SCIP_CALL( GCGrelaxTransOrigToMasterCons(origscip, branchdata->cons, &mastercons) );
   assert(mastercons != NULL);

   /* add constraint to the master problem */
   SCIP_CALL( SCIPaddConsNode(scip, SCIPgetCurrentNode(scip), mastercons, NULL) );

   /* constraint was added locally to the node where it is needed, so we do not need to care about this
    * at the next activation of the node and can set the constraint pointer to NULL */
   SCIP_CALL( SCIPreleaseCons(scip, &branchdata->cons) );
   branchdata->cons = NULL;

   return SCIP_OKAY;
}

/** callback solved method */
static
GCG_DECL_BRANCHMASTERSOLVED(branchMasterSolvedOrig)
{
   assert(scip != NULL);
   assert(GCGisOriginal(scip));
   assert(branchdata != NULL);
   assert(branchdata->origvar != NULL);

   SCIPdebugMessage("branchMasterSolvedOrig: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == GCG_BOUNDTYPE_LOWER ? ">=" : branchdata->boundtype == GCG_BOUNDTYPE_UPPER ? "<=" : "==" ), branchdata->newbound);

   if( !SCIPisInfinity(scip, newlowerbound) && SCIPgetStage(GCGgetMasterprob(scip)) == SCIP_STAGE_SOLVING
      && SCIPisRelaxSolValid(GCGgetMasterprob(scip)) )
   {
      SCIP_CALL( SCIPupdateVarPseudocost(scip, branchdata->origvar,
            SCIPgetRelaxSolVal(scip, branchdata->origvar) - branchdata->oldvalue,
            newlowerbound - branchdata->olddualbound, 1.0) );
   }

   return SCIP_OKAY;
}

/** callback deletion method for branching data */
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteOrig)
{
   assert(scip != NULL);
   assert(branchdata != NULL);

   if( *branchdata == NULL )
      return SCIP_OKAY;

   SCIPdebugMessage("branchDataDeleteOrig: %s %s %f\n", SCIPvarGetName((*branchdata)->origvar),
      ( (*branchdata)->boundtype == GCG_BOUNDTYPE_LOWER ? ">=" : (*branchdata)->boundtype == GCG_BOUNDTYPE_UPPER ? "<=" : "==" ),
      (*branchdata)->newbound);

   /* release constraint */
   if( (*branchdata)->cons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*branchdata)->cons) );
      (*branchdata)->cons = NULL;
   }

   SCIPfreeBlockMemoryNull(scip, branchdata);
   *branchdata = NULL;

   return SCIP_OKAY;
}

static
SCIP_DECL_BRANCHFREE(branchFreeOrig)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_HASHMAP* varhashmap;

   branchruledata = SCIPbranchruleGetData(branchrule);
   varhashmap = branchruledata->varhashmap;

   SCIPfreeBlockMemoryArray(scip, &branchruledata->score, branchruledata->maxvars);
   SCIPfreeBlockMemoryArray(scip, &branchruledata->strongbranchscore, branchruledata->maxvars);
   SCIPfreeBlockMemoryArray(scip, &branchruledata->uniqueblockflags, branchruledata->maxvars);
   SCIPfreeBlockMemoryArray(scip, &branchruledata->lastevalnode, branchruledata->maxvars);

   if( branchruledata->varhashmap != NULL )
   {
      SCIPhashmapFree(&varhashmap);
   }

   SCIPfreeMemory(scip, branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpOrig)
{  /*lint --e{715}*/
   SCIP* origscip;

   //SCIPdebugMessage("Execlp method of orig branching\n");

   /* get original problem */
   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   if( GCGcurrentNodeIsGeneric(scip) )
   {
      SCIPdebugMessage("Not executing orig branching, node was branched by generic branchrule\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* if the transferred master solution is feasible, the current node is solved to optimality and can be pruned */
   if( GCGrelaxIsOrigSolFeasible(origscip) )
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMessage("solution was feasible, node can be cut off!");
   }

   if( SCIPgetNExternBranchCands(origscip) > 0 )
   {
      assert(SCIPisRelaxSolValid(origscip));
      SCIP_CALL( branchExtern(origscip, branchrule, result) );
   }

   return SCIP_OKAY;
}

/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextOrig)
{  /*lint --e{715}*/
   SCIP* origscip;

   SCIPdebugMessage("Execext method of orig branching\n");

   /* get original problem */
   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   if( GCGcurrentNodeIsGeneric(scip) )
   {
      SCIPdebugMessage("Not executing orig branching, node was branched by generic branchrule\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* if the transferred master solution is feasible, the current node is solved to optimality and can be pruned */
   if( GCGrelaxIsOrigSolFeasible(origscip) )
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMessage("solution was feasible, node can be cut off!");
   }
   SCIP_CALL( branchExtern(origscip, branchrule, result) );

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitOrig)
{
   SCIP* origprob;
   SCIP_BRANCHRULEDATA* branchruledata;

   origprob = GCGmasterGetOrigprob(scip);
   assert(branchrule != NULL);
   assert(origprob != NULL);

   SCIPdebugMessage("Init orig branching rule\n");

   SCIP_CALL( GCGrelaxIncludeBranchrule( origprob, branchrule, branchActiveMasterOrig,
         branchDeactiveMasterOrig, branchPropMasterOrig, branchMasterSolvedOrig, branchDataDeleteOrig) );

   branchruledata = SCIPbranchruleGetData(branchrule);
   branchruledata->lastcand = 0;
   branchruledata->nvars = 0;
   branchruledata->maxvars = 0;
   this_branchruledata = branchruledata;
   SCIP_CALL( SCIPhashmapCreate(&(branchruledata->varhashmap), SCIPblkmem(scip), SCIPgetNVars(scip)) );

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsOrig)
{  /*lint --e{715}*/
   int i;
   SCIP* origscip;

   /* branching candidates */
   SCIP_VAR** branchcands;
   int nbranchcands;
   int npriobranchcands;

   /* values for choosing the variable to branch on */
   SCIP_VAR* branchvar;
   SCIP_Real solval;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("EXECPS orig branching rule\n");

   /* get original problem */
   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   if( GCGcurrentNodeIsGeneric(scip) )
   {
      SCIPdebugMessage("Not executing orig branching, node was branched by generic branchrule\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Execps method of orig branching\n");

   *result = SCIP_DIDNOTRUN;
   if( SCIPgetStage(scip) > SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* get the branching candidates */
   SCIP_CALL( SCIPgetPseudoBranchCands(origscip, &branchcands, &nbranchcands, &npriobranchcands) );

   branchvar = NULL;
   solval = 0.0;

   /* branch on an integer variable belonging to a unique block with fractional value */
   for( i = 0; i < npriobranchcands; i++ )
   {
      assert(GCGvarIsOriginal(branchcands[i]));

      /* variable belongs to no block or the block is not unique */
      if( GCGvarGetBlock(branchcands[i]) <= -1 || GCGgetNIdenticalBlocks(origscip, GCGvarGetBlock(branchcands[i])) != 1 )
         continue;

      branchvar = branchcands[i];
      lb = SCIPvarGetLbLocal(branchvar);
      ub = SCIPvarGetUbLocal(branchvar);

      assert(ub - lb > 0.8);

      /*  if the bounds of the branching variable x are finite, then the solution value
       *  is floor((lb + ub)/2)) + 0.5,
       *  otherwise the solution value is set to a finite bound
       *  if no finite bound exists, the solution value is set to 0.
       */
      if( !SCIPisInfinity(origscip, ub) && !SCIPisInfinity(origscip, -lb) )
         solval =  SCIPfeasFloor(scip, (ub + lb) / 2.0) + 0.5;
      else if( !SCIPisInfinity(origscip, -lb) )
         solval = lb;
      else if( !SCIPisInfinity(origscip, ub) )
         solval = ub;
      else
         solval = 0.0;

      break;
   }

   /* we did not find a variable to branch on so far, so we look for an unfixed linking variable or an integer variable
    * that belongs to no block but was directly transferred to the master problem
    */
   if( branchvar == NULL )
   {
      for( i = 0; i < npriobranchcands; i++ )
      {
         assert(GCGvarIsOriginal(branchcands[i]));

         /* continue if variable belongs to a block */
         if( GCGvarGetBlock(branchcands[i]) > -1 )
            continue;

         /* check that blocks of linking variable are unique */
         if( GCGoriginalVarIsLinking(branchcands[i]) )
         {
            int nvarblocks;
            int* varblocks;
            SCIP_Bool unique;
            int j;

            nvarblocks = GCGlinkingVarGetNBlocks(branchcands[i]);
            SCIP_CALL( SCIPallocBufferArray(origscip, &varblocks, nvarblocks) );
            SCIP_CALL( GCGlinkingVarGetBlocks(branchcands[i], nvarblocks, varblocks) );

            unique = TRUE;
            for( j = 0; j < nvarblocks; ++j )
               if( GCGgetNIdenticalBlocks(origscip, varblocks[j]) != 1 )
                  unique = FALSE;

            SCIPfreeBufferArray(origscip, &varblocks);

            if( !unique )
               continue;
         }

         branchvar = branchcands[i];
         lb = SCIPvarGetLbLocal(branchvar);
         ub = SCIPvarGetUbLocal(branchvar);

         assert(ub - lb > 0.8);

         /*  if the bounds of the branching variable x are finite, then the solution value
          *  is floor((lb + ub)/2)) + 0.5,
          *  otherwise the solution value is set to a finite bound
          *  if no finite bound exists, the solution value is set to 0.
          */
         if( !SCIPisInfinity(origscip, ub) && !SCIPisInfinity(origscip, -lb) )
            solval =  SCIPfeasFloor(scip, (ub + lb) / 2.0) + 0.5;
         else if( !SCIPisInfinity(origscip, -lb) )
            solval = lb;
         else if( !SCIPisInfinity(origscip, ub) )
            solval = ub;
         else
            solval = 0.0;

         break;
      }
   }

   if( branchvar == NULL )
   {
      SCIPdebugMessage("Original branching rule could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(branchvar != NULL);

   SCIP_CALL( branchVar(origscip, branchrule, branchvar, solval, FALSE, FALSE) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}



/*
 * branching specific interface methods
 */

/** creates the branching on original variable branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleOrig(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP* origscip;
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   SCIPdebugMessage("Include orig branching rule\n");

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
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitOrig) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpOrig) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextOrig) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsOrig) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeOrig) );

   /* add original variable branching rule parameters */

   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/orig/enforcebycons",
         "should bounds on variables be enforced by constraints(TRUE) or by bounds(FALSE)",
         &branchruledata->enforcebycons, FALSE, DEFAULT_ENFORCEBYCONS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/orig/mostfrac",
         "should branching be performed on the most fractional variable instead of the first variable?",
         &branchruledata->mostfrac, FALSE, DEFAULT_MOSTFRAC, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/orig/usepseudocosts",
         "should pseudocosts be used to determine the variable on which the branching is performed?",
         &branchruledata->usepseudocosts, FALSE, DEFAULT_USEPSEUDO, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/orig/usepsstrong",
         "should strong branching with propagation be used to determine the variable on which the branching is performed?",
         &branchruledata->usepsstrong, FALSE, DEFAULT_USEPSSTRONG, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/orig/usestrong",
         "should strong branching be used to determine the variable on which the branching is performed?",
         &branchruledata->usestrong, FALSE, DEFAULT_USESTRONG, NULL, NULL) );

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

   SCIP_CALL( SCIPaddIntParam(origscip, "branching/bp_strong/phase0outcands",
         "maximum number of output candidates from phase 0",
         &branchruledata->phasezerooutcands, FALSE, DEFAULT_PHASE0OUTCANDS, 1, 100000, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "branching/bp_strong/phase1outcands",
         "maximum number of output candidates from phase 1",
         &branchruledata->phaseoneoutcands, FALSE, DEFAULT_PHASE1OUTCANDS, 1, 100000, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(origscip, "branching/bp_strong/gapweight",
         "how much impact should the nodegap have on the number of precisely evaluated candidates?",
         &branchruledata->gapweight, FALSE, DEFAULT_GAPWEIGHT, 0, 1, NULL, NULL) );


   /* notify cons_integralorig about the original variable branching rule */
   SCIP_CALL( GCGconsIntegralorigAddBranchrule(scip, branchrule) );

   return SCIP_OKAY;
}

/** get the original variable on which the branching was performed */
SCIP_VAR* GCGbranchOrigGetOrigvar(
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   )
{
   assert(branchdata != NULL);

   return branchdata->origvar;
}

/** get the type of the new bound which resulted of the performed branching */
GCG_BOUNDTYPE GCGbranchOrigGetBoundtype(
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   )
{
   assert(branchdata != NULL);

   return branchdata->boundtype;
}

/** get the new bound which resulted of the performed branching */
SCIP_Real GCGbranchOrigGetNewbound(
   GCG_BRANCHDATA*       branchdata          /**< branching data */
   )
{
   assert(branchdata != NULL);

   return branchdata->newbound;
}

/** updates extern branching candidates before branching */
SCIP_RETCODE GCGbranchOrigUpdateExternBranchcands(
   SCIP*                 scip               /**< SCIP data structure */
)
{
   SCIP_VAR** origvars;
   int norigvars;
   int i;

   assert(GCGisOriginal(scip));

   origvars = SCIPgetVars(scip);
   norigvars = SCIPgetNVars(scip);
   assert(origvars != NULL);

   SCIPclearExternBranchCands(scip);

   /* store branching candidates */
   for( i = 0; i < norigvars; i++ )
   {
      if( SCIPvarGetType(origvars[i]) <= SCIP_VARTYPE_INTEGER && !SCIPisFeasIntegral(scip, SCIPgetRelaxSolVal(scip, origvars[i])) )
      {
         assert(!SCIPisEQ(scip, SCIPvarGetLbLocal(origvars[i]), SCIPvarGetUbLocal(origvars[i])));

         SCIP_CALL( SCIPaddExternBranchCand(scip, origvars[i], SCIPgetRelaxSolVal(scip,
            origvars[i]) - SCIPfloor(scip, SCIPgetRelaxSolVal(scip, origvars[i])),
            SCIPgetRelaxSolVal(scip, origvars[i])) );
      }
   }
   SCIPdebugMessage("updated relaxation branching candidates\n");

   return SCIP_OKAY;
}
