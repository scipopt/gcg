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

/**@file    branch_compbnd.h
 *
 * @brief   component bound branching rule
 * @author  Til Mohr
 *
 * This is an implementation of the component bound branching rule based on the papers:
 *
 * F. Vanderbeck.
 * On Dantzig-Wolfe decomposition in integer programming and ways to perform branching in a branch-and-price algorithm.
 * Oper. Res. 48(1):111-128 (2000)
 *
 * J. Desrosiers, M. L¨ubbecke, G. Desaulniers,
 * J. B. Gauthier (Juin 2024). Branch-and-Price, Technical report,
 * Les Cahiers du GERAD G–2024–36, GERAD, HEC Montr´eal, Canada.
 *
 * Vanderbeck, François, and Laurence A. Wolsey. "An exact algorithm for IP column generation."
 * Operations research letters 19.4 (1996): 151-159.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

// will write information about the component bound sequences to a file (slow!)
//#define COMPBND_STATS

#include "gcg/branch_compbnd.h"
#include "gcg/cons_integralorig.h"
#include "gcg/cons_masterbranch.h"
#include "gcg/gcg.h"
#include "gcg/pricer_gcg.h"
#include "gcg/relax_gcg.h"
#include "gcg/pub_extendedmasterconsdata.h"
#include "gcg/extendedmasterconsdata.h"
#include "gcg/pub_gcgvar.h"
#include "gcg/type_extendedmasterconsdata.h"
#include "gcg/type_branchgcg.h"
#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"

#define BRANCHRULE_NAME            "compbnd"                      /**< name of branching rule */
#define BRANCHRULE_DESC            "component bound branching"    /**< short description of branching rule */
#define BRANCHRULE_PRIORITY        -200                           /**< priority of this branching rule */
#define BRANCHRULE_MAXDEPTH        -1                             /**< maximal depth level of the branching rule */
#define BRANCHRULE_MAXBOUNDDIST    1.0                            /**< maximal relative distance from current node's
                                                                   dual bound to primal bound compared to best node's
                                                                   dual bound for applying branching */

/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   GCG*                  gcg;                            /**< GCG data structure */
   GCG_BRANCHRULE*       gcgbranchrule;                  /**< GCG branchrule structure */
   SCIP_Bool             useMaxRangeMidrangeHeuristic;   /** should the MaxRangeMidrangeHeuristic be used */
   SCIP_Bool             useMostDistinctMedianHeuristic; /** should the MostDistinctMedianHeuristic be used */
#ifdef COMPBND_STATS
   FILE*                 statsfile;          /**< file to write statistics to */
#endif
};

/** branching data */
struct GCG_BranchData
{
   GCG_BRANCH_TYPE       branchtype;         /**< type of branching, up or down */
   int                   constant;           /**< constant value of the branching constraint in the master problem - either lhs or rhs, depending on branchtype */
   GCG_COMPBND*          compBndSeq;         /**< component bound sequence which induce the current branching constraint */
   int                   compBndSeqSize;     /**< size of the component bound sequence compBndSeq */
   int                   blocknr;            /**< id of the pricing problem (or block) to which this branching constraint belongs */
   int                   nvars;              /**< number of master variables the last time the node has been visited - neccessary to later include newly generated master variables */
   GCG_EXTENDEDMASTERCONSDATA*    mastercons;         /**< master constraint along with its corresponding inferred pricing modifications */
};



/** define the signature of functions to choose the component bound */
#define CHOOSE_COMPBND(x) SCIP_RETCODE x (GCG* gcg, SCIP_VAR** satisfyingMastervars, int satisfyingMastervarsSize, SCIP_VAR** indexSet, int indexSetSize, GCG_COMPBND* compBndSeq, int compBndSeqSize, int blocknr, SCIP_VAR** selectedOrigvar, SCIP_Real* value)

/** define the signature of functions that chose one component bound sequence out of many */
#define CHOOSE_COMPBNDSEQ(x) SCIP_RETCODE x (GCG* gcg, GCG_COMPBND** compBndSeqList, int* compBndSeqSizeList, int compBndListSize, int blocknr, GCG_COMPBND** selectedCompBndSeq, int* selectedCompBndSeqSize)

/*
 * Local methods
 */

/* floor and ceil returning int */
#define FLOOR(scip, x) (SCIPconvertRealToInt((scip), SCIPfloor((scip), (x))))
#define CEIL(scip, x) (SCIPconvertRealToInt((scip), SCIPceil((scip), (x))))

/** set components of componentbound sequence to NULL */
static
void freeComponentBoundSequence(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_COMPBND**         compBndSeq,         /**< component bound sequence to set to NULL */
   int*                  compBndSeqSize      /**< size of the component bound sequence compBndSeq */
   )
{
   assert(scip != NULL);
   assert(GCGisMaster(scip));

   SCIPfreeBlockMemoryArrayNull(scip, compBndSeq, *compBndSeqSize);
   assert(*compBndSeq == NULL);
   *compBndSeqSize = 0;
}

/** initialize branchdata at the node */
static
SCIP_RETCODE initNodeBranchdata(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_BRANCHDATA**      nodebranchdata,     /**< branching data to set */
   GCG_BRANCH_TYPE       branchtype,         /**< type of branch to generate */
   int                   constant,           /**< constant value of the branching constraint in the master problem - either lhs or rhs, depending on branchtype */
   GCG_COMPBND*          compBndSeq,         /**< component bound sequence which induce the current branching constraint */
   int                   compBndSeqSize,     /**< size of the component bound sequence compBndSeq */
   int                   blocknr             /**< block we are branching in */
   )
{
   assert(GCGisMaster(scip));

   SCIP_CALL( SCIPallocBlockMemory(scip, nodebranchdata) );

   (*nodebranchdata)->branchtype = branchtype;
   (*nodebranchdata)->blocknr = blocknr;
   (*nodebranchdata)->constant = constant;
   (*nodebranchdata)->compBndSeqSize = compBndSeqSize;
   (*nodebranchdata)->nvars = 0;
   (*nodebranchdata)->mastercons = NULL;

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &((*nodebranchdata)->compBndSeq), compBndSeq, compBndSeqSize));

   return SCIP_OKAY;
}

/** simplify a component bound sequence by strengthening bounds of the same sense on the same component */
static
SCIP_RETCODE simplifyComponentBoundSequence(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_COMPBND**         compBndSeq,         /**< Component Bound Sequence defining the nodes */
   int*                  compBndSeqSize      /**< size of compBndSeq */
   )
{
   GCG_COMPBND* newCompBndSeq;
   int newCompBndSeqSize;
   SCIP_Bool* alreadyAdded;
   int i;
   int j;

   assert(GCGisMaster(scip));

   newCompBndSeq = NULL;
   newCompBndSeqSize = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &alreadyAdded, *compBndSeqSize) );
   for( i = 0; i < *compBndSeqSize; ++i )
   {
      alreadyAdded[i] = FALSE;
   }


   for( i=0; i<*compBndSeqSize; ++i )
   {
      if( alreadyAdded[i] )
         continue;

      SCIP_VAR* component = (*compBndSeq)[i].component;
      GCG_COMPBND_SENSE sense = (*compBndSeq)[i].sense;
      SCIP_Real bound = (*compBndSeq)[i].bound;

      for( j=i+1; j<*compBndSeqSize; ++j )
      {
         if( component != (*compBndSeq)[j].component )
            continue;
         if( (*compBndSeq)[i].sense != (*compBndSeq)[j].sense )
            continue;

         alreadyAdded[j] = TRUE;

         if( (*compBndSeq)[i].sense == GCG_COMPBND_SENSE_LE )
         {
            bound = MIN((*compBndSeq)[j].bound, bound);
         }
         else if( (*compBndSeq)[i].sense == GCG_COMPBND_SENSE_GE )
         {
            bound = MAX((*compBndSeq)[j].bound, bound);
         }
      }

      alreadyAdded[i] = TRUE;

      if( newCompBndSeqSize == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newCompBndSeq, 1) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &newCompBndSeq, newCompBndSeqSize, newCompBndSeqSize+1) );
      }

      newCompBndSeq[newCompBndSeqSize].component = component;
      newCompBndSeq[newCompBndSeqSize].sense = sense;
      newCompBndSeq[newCompBndSeqSize].bound = bound;
      newCompBndSeqSize += 1;
   }

   assert(1 <= newCompBndSeqSize && newCompBndSeqSize <= *compBndSeqSize);
   SCIPdebugMessage("Simplified compBndSeq from %d to %d\n", *compBndSeqSize, newCompBndSeqSize);
   for (i = 0; i < newCompBndSeqSize; ++i)
   {
      SCIPdebugMessage("compBndSeq[%d]: %s %s %d\n", i, SCIPvarGetName(newCompBndSeq[i].component),
                        newCompBndSeq[i].sense == GCG_COMPBND_SENSE_LE ? "<=" : ">=",
                        newCompBndSeq[i].bound);
   }

   SCIPfreeBlockMemoryArray(scip, &alreadyAdded, *compBndSeqSize);
   assert(alreadyAdded == NULL);

   // free old compBndSeq and replace it with newB
   freeComponentBoundSequence(scip, compBndSeq, compBndSeqSize);
   assert(*compBndSeq == NULL);
   *compBndSeq = newCompBndSeq;
   *compBndSeqSize = newCompBndSeqSize;

   return SCIP_OKAY;
}

/** computes the generator of mastervar for the entry in origvar
 * @return entry of the generator corresponding to origvar */
static
SCIP_Real getGeneratorEntryCol(
   SCIP_VAR**            solvars,            /**< column solution variables */
   SCIP_Real*            solvals,            /**< column solution values */
   int                   nsolvars,           /**< number of column solution variables */
   SCIP_VAR*             origvar             /**< corresponding origvar */
   )
{
   int i;
   SCIP_VAR* pricingvar;
   SCIP_VAR* cmp;

   assert(origvar != NULL);
   assert(GCGvarIsOriginal(origvar));

   pricingvar = GCGoriginalVarGetPricingVar(origvar);

   for( i = 0; i < nsolvars; ++i )
   {
      if( GCGvarIsOriginal(solvars[i]) )
         cmp = origvar;
      else
         cmp = pricingvar;

      if( solvars[i] == cmp )
      {
         return solvals[i];
      }
   }

   return 0.0;
}

/** computes the generator of mastervar for the entry in origvar
 * @return entry of the generator corresponding to origvar */
static
SCIP_Real getGeneratorEntry(
   SCIP_VAR*             mastervar,          /**< current mastervariable */
   SCIP_VAR*             origvar             /**< corresponding origvar */
   )
{
   SCIP_Real entry = GCGmasterVarGetOrigval(mastervar, origvar);

   return entry != SCIP_INVALID ? entry : 0.;
}

/** whether a master variable is in the component bound sequence or not */
static
SCIP_Bool isColInCompBndSeq(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            solvars,            /**< column solution variables */
   SCIP_Real*            solvals,            /**< column solution values */
   int                   nsolvars,           /**< number of column solution variables */
   GCG_COMPBND*          compBndSeq,         /**< component bound sequence to check */
   int                   compBndSeqSize      /**< size of compBndSeq */
   )
{
   int i;
   SCIP_Real generatorentry;

   assert(compBndSeqSize > 0);
   assert(compBndSeq != NULL);

   for( i = 0; i < compBndSeqSize; ++i )
   {
      generatorentry = getGeneratorEntryCol(solvars, solvals, nsolvars, compBndSeq[i].component);
      if ( (compBndSeq[i].sense == GCG_COMPBND_SENSE_GE && !SCIPisGE(scip, generatorentry, compBndSeq[i].bound)) ||
         (compBndSeq[i].sense == GCG_COMPBND_SENSE_LE && !SCIPisLE(scip, generatorentry, compBndSeq[i].bound)) )
      {
         return FALSE;
      }
   }

   return TRUE;
}

/** whether a master variable is in compBndSeq or not */
static
SCIP_Bool isMasterVarInCompBndSeq(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             mastervar,          /**< master variable to check */
   GCG_COMPBND*          compBndSeq,         /**< component bound sequence to check */
   int                   compBndSeqSize,     /**< size of compBndSeq */
   int                   blocknr             /**< block we are branching in */
   )
{
   int i;
   SCIP_Real generatorentry;

   assert(mastervar != NULL);
   assert(compBndSeqSize > 0);
   assert(compBndSeq != NULL);

   if( !GCGisMasterVarInBlock(mastervar, blocknr) )
      return FALSE;

   for( i = 0; i < compBndSeqSize; ++i )
   {
      generatorentry = getGeneratorEntry(mastervar, compBndSeq[i].component);
      if ( (compBndSeq[i].sense == GCG_COMPBND_SENSE_GE && !SCIPisGE(scip, generatorentry, compBndSeq[i].bound)) ||
         (compBndSeq[i].sense == GCG_COMPBND_SENSE_LE && !SCIPisLE(scip, generatorentry, compBndSeq[i].bound)) )
      {
         return FALSE;
      }
   }

   return TRUE;
}

/** adds a variable to a branching constraint */
static
SCIP_RETCODE addVarToMasterbranch(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_VAR*             mastervar,          /**< the variable to add */
   GCG_BRANCHDATA*       branchdata,         /**< branching data structure where the variable should be added */
   SCIP_Bool*            added               /**< whether the variable was added */
)
{
   SCIP_Real coef;
   SCIP_VAR** origvars;
   SCIP_Real* origvals;
   int norigvars;
   SCIP* masterprob = GCGgetMasterprob(gcg);

   assert(masterprob != NULL);
   assert(mastervar != NULL);
   assert(branchdata != NULL);
   assert(added != NULL);

   origvars = GCGmasterVarGetOrigvars(mastervar);
   norigvars = GCGmasterVarGetNOrigvars(mastervar);
   origvals = GCGmasterVarGetOrigvals(mastervar);

   *added = FALSE;

   SCIP_CALL( GCGextendedmasterconsGetCoeff(gcg, branchdata->mastercons, origvars, origvals, norigvars, GCGvarGetBlock(mastervar), &coef) );

   if( !SCIPisZero(masterprob, coef) )
   {
      SCIPdebugMessage("Adding variable %s to branching constraint: %s %d\n", SCIPvarGetName(mastervar), branchdata->branchtype == GCG_BRANCH_DOWN ? "<=" : ">=", branchdata->constant);
      SCIP_CALL( GCGextendedmasterconsAddMasterVar(gcg, branchdata->mastercons, mastervar, coef) );
      *added = TRUE;
   }

   return SCIP_OKAY;
}

/** callback new column method */
static
GCG_DECL_BRANCHGETEXTENDEDMASTERCONSCOEFF(branchGetExtendedmasterconsCoeffCompBnd)
{
   SCIP* masterprob;

   masterprob = GCGgetMasterprob(gcg);

   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));

   assert(branchdata != NULL);
   assert(branchdata->mastercons == extendedmasterconsdata);

   *coef = 0.0;

   if( probnr == -1 || probnr != branchdata->blocknr )
      return SCIP_OKAY;

   if( !isColInCompBndSeq(masterprob, solvars, solvals, nsolvars, branchdata->compBndSeq, branchdata->compBndSeqSize) )
      return SCIP_OKAY;

   *coef = 1.0;

   return SCIP_OKAY;
}

/** creates the constraint for branching directly on a master variable */
static
SCIP_RETCODE createBranchingCons(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branch rule structure */
   SCIP_NODE*            node,               /**< node to add constraint */
   GCG_BRANCHDATA*       branchdata          /**< branching data structure */
)
{
   SCIP_Longint uuid;

   SCIP* pricingscip;
   SCIP* masterprob;
   SCIP_BRANCHRULEDATA* branchruledata;

   SCIP_VAR** mastervars;
   int nmastervars;
   int i;
   SCIP_Bool added = FALSE;

   SCIP_CONS* branchcons;

   SCIP_VAR* coefvar;
   SCIP_VAR** additionalvars;
   int nadditionalvars;
   SCIP_CONS** additionalcons;
   int nadditionalcons;
   GCG_PRICINGMODIFICATION** pricingmods; // will always only contain one element in this branching rule

   char name[SCIP_MAXSTRLEN];

   masterprob = GCGgetMasterprob(gcg);

   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));
   assert(node != NULL);
   assert(branchdata != NULL);
   assert(branchdata->mastercons == NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   uuid = SCIPnodeGetNumber(node);

   mastervars = SCIPgetVars(masterprob);
   nmastervars = SCIPgetNVars(masterprob);
   assert(nmastervars >= 0);
   assert(nmastervars == 0 || mastervars != NULL);

   /*  create constraint for child */
   if (branchdata->branchtype == GCG_BRANCH_DOWN)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "compbnd_%d_child(%d_LE_%d)", uuid, branchdata->compBndSeqSize, branchdata->constant);
      SCIP_CALL( SCIPcreateConsLinear(masterprob, &branchcons, name, 0, NULL, NULL,
         -SCIPinfinity(masterprob), (SCIP_Real) branchdata->constant, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
   }
   else
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "compbnd_%d_child(%d_GE_%d)", uuid, branchdata->compBndSeqSize, branchdata->constant);
      SCIP_CALL( SCIPcreateConsLinear(masterprob, &branchcons, name, 0, NULL, NULL,
         (SCIP_Real) branchdata->constant, SCIPinfinity(masterprob), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
   }

   SCIP_CALL( SCIPaddConsNode(masterprob, node, branchcons, NULL) );

   pricingscip = GCGgetPricingprob(gcg, branchdata->blocknr);
   assert(pricingscip != NULL);

   assert(branchdata->mastercons == NULL);

   if( branchdata->branchtype == GCG_BRANCH_DOWN )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "compbnd_%d_down", uuid);
   }
   else
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "compbnd_%d_up", uuid);
   }

   char pricingvarname[SCIP_MAXSTRLEN];
   (void) SCIPsnprintf(pricingvarname, SCIP_MAXSTRLEN, "g(%s)", name);

   // create the g_x variable
   SCIP_CALL( GCGcreateInferredPricingVar(pricingscip, &coefvar, pricingvarname, 0.0, 1.0, TRUE, 0.0, SCIP_VARTYPE_BINARY, branchdata->blocknr) );

   // create the y_j variables

   /* @todo: when the pricing problem only contains binary variables, simplifications apply:
    *        we may not even need extra y_j variables to enforce the correct value of the y variable
    */

   nadditionalvars = branchdata->compBndSeqSize;
   SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &additionalvars, nadditionalvars) );
   for (i = 0; i < branchdata->compBndSeqSize; ++i)
   {
      (void) SCIPsnprintf(pricingvarname, SCIP_MAXSTRLEN, "y(%s,%s_%s_%d)", name, SCIPvarGetName(branchdata->compBndSeq[i].component), branchdata->compBndSeq[i].sense == GCG_COMPBND_SENSE_GE ? "GE" : "LE", branchdata->compBndSeq[i].bound);
      SCIP_CALL( GCGcreateInferredPricingVar(pricingscip, &additionalvars[i], pricingvarname, 0.0, 1.0, FALSE, 0.0, SCIP_VARTYPE_BINARY, branchdata->blocknr) );
   }

   // create the pricing constraints
   if( branchdata->branchtype == GCG_BRANCH_DOWN )
   {
      nadditionalcons = branchdata->compBndSeqSize + 1;
      SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &additionalcons, nadditionalcons) );

      /* g_x >= 1 + sum_{j=1}^{n} y_j - n
         * 1 - n <= g_x - sum_{j=1}^{n} y_j */
      char consname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "c(g(%s))", name);
      SCIP_CALL( SCIPcreateConsLinear(pricingscip, &additionalcons[0], consname, 0, NULL, NULL, (SCIP_Real) (1 - branchdata->compBndSeqSize), SCIPinfinity(pricingscip),
            TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[0], coefvar, 1.0) );
      for( i = 0; i < branchdata->compBndSeqSize; ++i)
      {
         SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[0], additionalvars[i], -1.0) );
      }

      for( i = 0; i < branchdata->compBndSeqSize; ++i)
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "c(y(%s,%s_%s_%d))", name, SCIPvarGetName(branchdata->compBndSeq[i].component), branchdata->compBndSeq[i].sense == GCG_COMPBND_SENSE_GE ? "GE" : "LE", branchdata->compBndSeq[i].bound);

         if( branchdata->compBndSeq[i].sense == GCG_COMPBND_SENSE_LE )
         {
            /* y_j >= ((bound + 1) - x_j) / ((bound + 1) - l_j)
               * ((bound + 1) - l_j) * y_j >= (bound + 1) - x_j
               * (bound + 1) <= x_j + ((bound + 1) - l_j) * y_j */
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->compBndSeq[i].component);
            SCIP_Real bound = (SCIP_Real) (branchdata->compBndSeq[i].bound + 1);
            SCIP_Real lowerbound = SCIPvarGetLbGlobal(branchdata->compBndSeq[i].component);
            assert(SCIPisPositive(pricingscip, bound - lowerbound));
            SCIP_CALL( SCIPcreateConsVarbound(pricingscip, &additionalcons[i+1], consname, pricing_var, additionalvars[i], bound - lowerbound, bound, SCIPinfinity(pricingscip),
               TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         }
         else
         {
            /* y_j >= (x_j - (bound - 1)) / (u_j - (bound - 1))
               * (u_j - (bound - 1)) * y_j >= x_j - (bound - 1)
               * -(bound - 1) <= (u_j - (bound - 1)) * y_j - x_j
               * x_j + ((bound - 1) - u_j) * y_j <= (bound - 1) */
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->compBndSeq[i].component);
            SCIP_Real bound = (SCIP_Real) (branchdata->compBndSeq[i].bound - 1);
            SCIP_Real upperbound = SCIPvarGetUbGlobal(branchdata->compBndSeq[i].component);
            assert(SCIPisPositive(pricingscip, upperbound - bound));
            SCIP_CALL( SCIPcreateConsVarbound(pricingscip, &additionalcons[i+1], consname, pricing_var, additionalvars[i], bound - upperbound, -SCIPinfinity(pricingscip), bound,
               TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         }
      }
   }
   else
   {
      nadditionalcons = branchdata->compBndSeqSize * 2;
      SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &additionalcons, nadditionalcons) );

      /* g_x <= y_j
         * 0 <= y_j - g_y*/
      char consname[SCIP_MAXSTRLEN];
      for( i = 0; i < branchdata->compBndSeqSize; ++i )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "c(g(%s,%s))", name, SCIPvarGetName(additionalvars[i]));
         SCIP_CALL( SCIPcreateConsVarbound(pricingscip, &additionalcons[i], consname, additionalvars[i], coefvar, -1.0, 0.0, SCIPinfinity(pricingscip),
            TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      }

      for( i = 0; i < branchdata->compBndSeqSize; ++i )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "c1(y(%s,%s_%s_%d))", name, SCIPvarGetName(branchdata->compBndSeq[i].component), branchdata->compBndSeq[i].sense == GCG_COMPBND_SENSE_GE ? "GE" : "LE", branchdata->compBndSeq[i].bound);

         if( branchdata->compBndSeq[i].sense == GCG_COMPBND_SENSE_LE )
         {
            /* y_j <= (u_j - x_j) / (u_j - bound)
               * (u_j - bound) * y_j <= u_j - x_j
               * x_j + (u_j - bound) * y_j <= u_j */
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->compBndSeq[i].component);
            SCIP_Real bound = (SCIP_Real) branchdata->compBndSeq[i].bound;
            SCIP_Real upperbound = SCIPvarGetUbGlobal(branchdata->compBndSeq[i].component);
            assert(SCIPisPositive(pricingscip, upperbound - bound));
            SCIP_CALL( SCIPcreateConsVarbound(pricingscip, &additionalcons[i+branchdata->compBndSeqSize], consname, pricing_var, additionalvars[i], upperbound - bound, -SCIPinfinity(pricingscip), upperbound,
               TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         }
         else
         {
            /* y_j <= (x_j - l_j) / (bound - l_j)
               * (bound - l_j) * y_j <= x_j - l_j
               * (bound - l_j) * y_j - x_j <= -l_j
               * l_j <= x_j + (l_j - bound) * y_j */
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->compBndSeq[i].component);
            SCIP_Real bound = (SCIP_Real) branchdata->compBndSeq[i].bound;
            SCIP_Real lowerbound = SCIPvarGetLbGlobal(branchdata->compBndSeq[i].component);
            assert(SCIPisPositive(pricingscip, bound - lowerbound));
            SCIP_CALL( SCIPcreateConsVarbound(pricingscip, &additionalcons[i+branchdata->compBndSeqSize], consname, pricing_var, additionalvars[i], lowerbound - bound, lowerbound, SCIPinfinity(pricingscip),
               TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         }
      }
   }

   // create the pricing modification
   SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &pricingmods, 1) );
   assert(pricingmods != NULL);
   SCIP_CALL( GCGpricingmodificationCreate(
      gcg,
      &(pricingmods[0]),
      branchdata->blocknr,
      coefvar,
      additionalvars,
      nadditionalvars,
      additionalcons,
      nadditionalcons
   ) );

   // create the master constraint
   SCIP_CALL( GCGbranchCreateExtendedmastercons(gcg, branchruledata->gcgbranchrule, &branchdata->mastercons, branchcons, pricingmods, 1, branchdata) );

   // add the variables to the constraint
   for( i = 0; i < nmastervars; i++ )
   {
      added = FALSE;
      SCIP_CALL( addVarToMasterbranch(gcg, mastervars[i], branchdata, &added) );
   }

   return SCIP_OKAY;
}

/** for given component bound sequence create the 2 child nodes */
static
SCIP_RETCODE createChildNodesCompBndSeq(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   GCG_COMPBND*          compBndSeq,         /**< Component Bound Sequence defining the nodes */
   int                   compBndSeqSize,     /**< size of compBndSeq */
   int                   blocknr             /**< number of the block */
   )
{
   SCIP* origprob;
   SCIP* masterprob;
   int i;
   SCIP_Real constantSum;
   int nmastervars;
   SCIP_VAR** mastervars;
   GCG_BRANCHDATA* downBranchData;
   GCG_BRANCHDATA* upBranchData;
   char downChildname[SCIP_MAXSTRLEN];
   char upChildname[SCIP_MAXSTRLEN];
   SCIP_NODE* downChild;
   SCIP_CONS* downChildcons;
   SCIP_NODE* upChild;
   SCIP_CONS* upChildcons;

   assert(gcg != NULL);
   assert(compBndSeqSize > 0);
   assert(compBndSeq != NULL);

   origprob = GCGgetOrigprob(gcg);
   masterprob = GCGgetMasterprob(gcg);

   SCIPdebugMessage("Component bound branching rule Node creation for blocknr %d with %d identical blocks \n", blocknr, GCGgetNIdenticalBlocks(gcg, blocknr));

   /*  get variable data of the master problem */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(nmastervars >= 0);
   /* determine the constant value of the master constraint, i.e. the rhs for the down branch, and lhs for the up branch */
   constantSum = 0.0;
   for( i = 0; i < nmastervars; ++i )
   {
      if(
         GCGisMasterVarInBlock(
            mastervars[i],
            GCGgetBlockRepresentative(gcg, blocknr)
         )
         && isMasterVarInCompBndSeq(origprob, mastervars[i], compBndSeq, compBndSeqSize, blocknr)
      )
      {
         constantSum += SCIPgetSolVal(masterprob, NULL, mastervars[i]);
      }
   }
   // sanity check: the sum must be fractional, otherwise something went wrong during separation, i.e. compBndSeq is incorrect
   assert(!SCIPisFeasIntegral(origprob, constantSum));

   /* create two nodes */
   SCIPdebugMessage("Component bound branching rule: creating 2 nodes\n");
   SCIP_CALL( initNodeBranchdata(masterprob, &downBranchData, GCG_BRANCH_DOWN, FLOOR(origprob, constantSum), compBndSeq, compBndSeqSize, blocknr) );
   SCIP_CALL( initNodeBranchdata(masterprob, &upBranchData, GCG_BRANCH_UP, CEIL(origprob, constantSum), compBndSeq, compBndSeqSize, blocknr) );
   freeComponentBoundSequence(masterprob, &compBndSeq, &compBndSeqSize);
   assert(downBranchData != NULL);
   assert(upBranchData != NULL);
   assert(compBndSeq == NULL);
   assert(compBndSeqSize == 0);

   /* define names for origbranch constraints */
   (void) SCIPsnprintf(downChildname, SCIP_MAXSTRLEN, "node(%d) (last comp=%s %s %d) <= %d", blocknr,
      SCIPvarGetName(downBranchData->compBndSeq[downBranchData->compBndSeqSize-1].component),
      downBranchData->compBndSeq[downBranchData->compBndSeqSize-1].sense == GCG_COMPBND_SENSE_GE ? ">=": "<=",
      downBranchData->compBndSeq[downBranchData->compBndSeqSize-1].bound,
      downBranchData->constant);
   (void) SCIPsnprintf(upChildname, SCIP_MAXSTRLEN, "node(%d) (last comp=%s %s %d) >= %d", blocknr,
      SCIPvarGetName(upBranchData->compBndSeq[upBranchData->compBndSeqSize-1].component),
      upBranchData->compBndSeq[upBranchData->compBndSeqSize-1].sense == GCG_COMPBND_SENSE_GE ? ">=": "<=",
      upBranchData->compBndSeq[upBranchData->compBndSeqSize-1].bound,
      upBranchData->constant);

   /* add the nodes to the tree */
   SCIP_CALL( SCIPcreateChild(masterprob, &downChild, 0.0, SCIPgetLocalTransEstimate(masterprob)) );
   SCIP_CALL( GCGcreateConsMasterbranch(gcg, &downChildcons, downChildname, downChild,
      GCGconsMasterbranchGetActiveCons(gcg), branchrule, downBranchData, NULL, 0, 0) );
   SCIP_CALL( SCIPaddConsNode(masterprob, downChild, downChildcons, NULL) );
   SCIP_CALL( createBranchingCons(gcg, branchrule, downChild, downBranchData) );

   SCIP_CALL( SCIPcreateChild(masterprob, &upChild, 0.0, SCIPgetLocalTransEstimate(masterprob)) );
   SCIP_CALL( GCGcreateConsMasterbranch(gcg, &upChildcons, upChildname, upChild,
      GCGconsMasterbranchGetActiveCons(gcg), branchrule, upBranchData, NULL, 0, 0) );
   SCIP_CALL( SCIPaddConsNode(masterprob, upChild, upChildcons, NULL) );
   SCIP_CALL( createBranchingCons(gcg, branchrule, upChild, upBranchData) );

   /*  release constraints */
   SCIP_CALL( SCIPreleaseCons(masterprob, &upChildcons) );
   SCIP_CALL( SCIPreleaseCons(masterprob, &downChildcons) );

   return SCIP_OKAY;
}

/** determine fractionality of a set of mastervariables */
static
SCIP_Real calcFractionality(
   SCIP*                 masterscip,              /**< SCIP data structure */
   SCIP_VAR**            mastervars,              /**< mastervariables */
   int                   nmastervars              /**< number of mastervariables */
   )
{
   int i;
   SCIP_Real fractionality;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));

   assert((mastervars != NULL) || (nmastervars == 0));

   if( nmastervars == 0 )
      return 0.0;

   fractionality = 0.0;

   for( i = 0; i < nmastervars; ++i )
   {
      SCIP_Real solval = SCIPgetSolVal(masterscip, NULL, mastervars[i]);
      fractionality += SCIPfrac(masterscip, solval);
   }

   return fractionality;
}

/** method for initializing the indexSet for a given set of master variables
  *
  * the index set denotes the set of all integral original variables contained in the provided master variables
  */
static
SCIP_RETCODE initIndexSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           indexSet,           /**< set to initialize */
   int*                  indexSetSize,       /**< size of the index set */
   SCIP_VAR**            mastervars,         /**< mastervariables currently satisfying the component bound sequence in the specified block */
   int                   mastervarsSize      /**< size of mastervars */
   )
{
   int i;
   int j;
   SCIP_VAR** origvars;
   int norigvars;

   assert( scip != NULL);
   assert(GCGisMaster(scip));
   assert( indexSet != NULL);
   assert( indexSetSize != NULL);

   *indexSet = NULL;
   *indexSetSize = 0;

   for( i = 0; i < mastervarsSize; ++i )
   {
      origvars = GCGmasterVarGetOrigvars(mastervars[i]);
      norigvars = GCGmasterVarGetNOrigvars(mastervars[i]);

      if( *indexSetSize == 0 && norigvars > 0 )
      {
         for( j = 0; j < norigvars; ++j )
         {
            if( SCIPvarGetType(origvars[j]) > SCIP_VARTYPE_INTEGER )
               continue;

            if( *indexSetSize == 0 )
            {
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, indexSet, 1) );
            }
            else
            {
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, indexSet, *indexSetSize, *indexSetSize + 1) );
            }

            (*indexSet)[*indexSetSize] = origvars[j];
            *indexSetSize += 1;
         }
      }
      else
      {
         for( j = 0; j < norigvars; ++j )
         {
            int k;
            int oldsize = *indexSetSize;

            if( SCIPvarGetType(origvars[j]) > SCIP_VARTYPE_INTEGER )
               continue;

            for( k = 0; k < oldsize; ++k )
            {
               /*  if variable already in set */
               if( (*indexSet)[k] == origvars[j] )
                  break;

               if( k == oldsize-1 )
               {
                  /*  add variable to the end */
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, indexSet, *indexSetSize, *indexSetSize + 1) );
                  (*indexSet)[*indexSetSize] = origvars[j];
                  *indexSetSize += 1;
               }
            }
         }
      }
   }

   assert(*indexSetSize > 0 && *indexSet != NULL);

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** get the lower bound of a original variable from the component bound sequence, or the lb */
static
SCIP_Real getLowerBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar,            /**< original variable */
   GCG_COMPBND*          compBndSeq,         /**< component bound sequence */
   int                   compBndSeqSize      /**< size of compBndSeq */
   )
{
   int i;
   SCIP_Real lb;

   assert(scip != NULL);
   assert(origvar != NULL);

   lb = SCIPvarGetLbGlobal(origvar);

   for( i = 0; i < compBndSeqSize; ++i )
   {
      if(
         compBndSeq[i].sense == GCG_COMPBND_SENSE_GE
         && origvar == compBndSeq[i].component
         && SCIPisLT(scip, lb, compBndSeq[i].bound)
      )
      {
         lb = (SCIP_Real) compBndSeq[i].bound;
      }
   }

   return lb;
}

/** get the upper bound of a original variable from the component bound sequence, or the ub */
static
SCIP_Real getUpperBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar,            /**< original variable */
   GCG_COMPBND*          compBndSeq,         /**< component bound sequence */
   int                   compBndSeqSize      /**< size of compBndSeq */
   )
{
   int i;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(origvar != NULL);

   ub = SCIPvarGetUbGlobal(origvar);

   for( i = 0; i < compBndSeqSize; ++i )
   {
      if(
         compBndSeq[i].sense == GCG_COMPBND_SENSE_LE
         && origvar == compBndSeq[i].component
         && SCIPisLT(scip, compBndSeq[i].bound, ub)
      )
      {
         ub = (SCIP_Real) compBndSeq[i].bound;
      }
   }

   return ub;
}
#endif

/** separation algorithm determining multiple component bound sequences one could branch on */
static
SCIP_RETCODE separation_helper(
   GCG*                  gcg,                     /**< GCG data structure */
   SCIP_VAR**            satisfyingMastervars,    /**< mastervariables currently satisfying the component bound sequence in the specified block */
   int                   satisfyingMastervarsSize,/**< size of mastervars */
   GCG_COMPBND*          compBndSeq,              /**< current Component Bound Sequence defining the nodes */
   int                   compBndSeqSize,          /**< size of compBndSeq */
   int                   blocknr,                 /**< number of the block */
   CHOOSE_COMPBND((*chooseCompBnd)),              /**< method for choosing the variable to branch on */
   GCG_COMPBND***        foundCompBndSeqList,     /**< list of component bound sequences */
   int**                 foundCompBndSeqSizeList, /**< sizes of the component bound sequences */
   int*                  foundCompBndSeqListSize, /**< size of Blist */
   SCIP_RESULT*          result                   /**< pointer to store the result of the branching call */
   )
{
   SCIP_Real fractionality;
   int i;
   SCIP* masterprob;

   masterprob = GCGgetMasterprob(gcg);

   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));

   // Sanity check: All variables in satisfyingMastervars must satisfy the component bound sequence
   if( compBndSeqSize > 0 ) {
      for( i = 0; i < satisfyingMastervarsSize; ++i )
      {
         assert(isMasterVarInCompBndSeq(masterprob, satisfyingMastervars[i], compBndSeq, compBndSeqSize, blocknr));
      }
   }

   fractionality = calcFractionality(masterprob, satisfyingMastervars, satisfyingMastervarsSize);
   assert(fractionality >= 0.0);

   if( SCIPisEQ(masterprob, fractionality, 0.0) ) // should never happen - the branching scheme is sound and complete
   {
      // all variables are integral, nothing to do
      SCIPdebugMessage("All variables are integral, nothing to do\n");

      *result = SCIP_DIDNOTFIND;

      return SCIP_ERROR;
   }

   if( compBndSeqSize > 0 && !SCIPisFeasIntegral(masterprob, fractionality) )
   {
      /* now we branch on <= floor(fractionality) and >= ceil(fractionality)
       * this is handled by the createChildNodesCompBnd method
       */
      SCIPdebugMessage("fractionality is fractional, we need to branch\n");

      *result = SCIP_BRANCHED;

      // add current compBndSeq to Blist
      if( *foundCompBndSeqListSize == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, foundCompBndSeqList, 1) );
         SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, foundCompBndSeqSizeList, 1) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, foundCompBndSeqList, *foundCompBndSeqListSize, *foundCompBndSeqListSize + 1) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, foundCompBndSeqSizeList, *foundCompBndSeqListSize, *foundCompBndSeqListSize + 1) );
      }
      (*foundCompBndSeqList)[*foundCompBndSeqListSize] = compBndSeq;
      (*foundCompBndSeqSizeList)[*foundCompBndSeqListSize] = compBndSeqSize;
      *foundCompBndSeqListSize += 1;

      return SCIP_OKAY;
   }

   // the fractionality is integral, we need to impose an additional bound
   SCIP_VAR** indexSet = NULL;
   int indexSetSize = 0;
   SCIP_CALL( initIndexSet(masterprob, &indexSet, &indexSetSize, satisfyingMastervars, satisfyingMastervarsSize) );

   SCIP_VAR* selectedOrigvar = NULL;
   SCIP_Real value = 0.0;
   SCIP_CALL( chooseCompBnd(gcg, satisfyingMastervars, satisfyingMastervarsSize, indexSet, indexSetSize, compBndSeq, compBndSeqSize, blocknr, &selectedOrigvar, &value) );
   assert(selectedOrigvar != NULL);
   assert(SCIPvarGetLbGlobal(selectedOrigvar) < value && value < SCIPvarGetUbGlobal(selectedOrigvar));

   SCIPfreeBlockMemoryArray(masterprob, &indexSet, indexSetSize);
   indexSetSize = 0;

   // create two component bound options, xj <= floor(value) and xj >= floor(value)+1
   GCG_COMPBND* lowCompBndSeq;
   GCG_COMPBND* highCompBndSeq;
   int newCompBndSeqSize = compBndSeqSize + 1;
   // allocate memory for lowCompBndSeq and highCompBndSeq
   SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &lowCompBndSeq, newCompBndSeqSize) );

   SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &highCompBndSeq, newCompBndSeqSize) );
   // copy the old component bounds to lowCompBndSeq and highCompBndSeq
   for(i = 0; i < compBndSeqSize; ++i)
   {
      lowCompBndSeq[i] = compBndSeq[i];
      highCompBndSeq[i] = compBndSeq[i];
   }
   // add the new component bounds to lowCompBndSeq and highCompBndSeq
   lowCompBndSeq[compBndSeqSize] = (GCG_COMPBND){selectedOrigvar, GCG_COMPBND_SENSE_LE, FLOOR(masterprob, value)};
   highCompBndSeq[compBndSeqSize] = (GCG_COMPBND){selectedOrigvar, GCG_COMPBND_SENSE_GE, FLOOR(masterprob, value) + 1};

   freeComponentBoundSequence(masterprob, &compBndSeq, &compBndSeqSize);
   assert(compBndSeq == NULL);

   SCIPdebugMessage("lowCompBndSeq and highCompBndSeq after adding the new component bounds\n");
   for (i = 0; i < newCompBndSeqSize; ++i)
   {
      SCIPdebugMessage("lowCompBndSeq[%d]: %s %s %d\n", i, SCIPvarGetName(lowCompBndSeq[i].component),
                        lowCompBndSeq[i].sense == GCG_COMPBND_SENSE_LE ? "<=" : ">=",
                        lowCompBndSeq[i].bound);
      SCIPdebugMessage("highCompBndSeq[%d]: %s %s %d\n", i, SCIPvarGetName(highCompBndSeq[i].component),
                        highCompBndSeq[i].sense == GCG_COMPBND_SENSE_LE ? "<=" : ">=",
                        highCompBndSeq[i].bound);
   }

   // assign the master variables to lowMastervars and highMastervars, depending on whether they satisfy the new component bound sequences
   SCIP_VAR** lowMastervars = NULL;
   SCIP_VAR** highMastervars = NULL;
   int lowMastervarsSize = 0;
   int highMastervarsSize = 0;
   for( i=0; i<satisfyingMastervarsSize; ++i )
   {
#ifndef NDEBUG
      SCIP_Bool satisfiesLowCompBndSeq = FALSE;
      SCIP_Bool satisfiesHighCompBndSeq = FALSE;
#endif
      if( isMasterVarInCompBndSeq(masterprob, satisfyingMastervars[i], lowCompBndSeq, newCompBndSeqSize, blocknr) )
      {
         // increase the size of lowMastervars by 1, and add the current variable to lowMastervars
         if( lowMastervarsSize == 0 )
            SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &lowMastervars, lowMastervarsSize+1) );
         else
            SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &lowMastervars, lowMastervarsSize, lowMastervarsSize+1) );

         lowMastervars[lowMastervarsSize] = satisfyingMastervars[i];
         lowMastervarsSize += 1;
#ifndef NDEBUG
         satisfiesLowCompBndSeq = TRUE;
#endif
      }
      if( isMasterVarInCompBndSeq(masterprob, satisfyingMastervars[i], highCompBndSeq, newCompBndSeqSize, blocknr) )
      {
         // increase the size of highMastervars by 1, and add the current variable to highMastervars
         if( highMastervarsSize == 0 )
            SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &highMastervars, highMastervarsSize+1) );
         else
            SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &highMastervars, highMastervarsSize, highMastervarsSize+1) );

         highMastervars[highMastervarsSize] = satisfyingMastervars[i];
         highMastervarsSize += 1;
#ifndef NDEBUG
         satisfiesHighCompBndSeq = TRUE;
#endif
      }
#ifndef NDEBUG
      assert(satisfiesLowCompBndSeq + satisfiesHighCompBndSeq <= 1);
#endif
   }
   SCIPdebugMessage("lowMastervarsSize: %d, highMastervarsSize: %d\n", lowMastervarsSize, highMastervarsSize);
   assert(lowMastervarsSize > 0);
   assert(highMastervarsSize > 0);

#ifndef NDEBUG
   // determine the fractionality of lowMastervars and highMastervars
   SCIP_Real fractionality1 = calcFractionality(masterprob, lowMastervars, lowMastervarsSize);
   SCIP_Real fractionality2 = calcFractionality(masterprob, highMastervars, highMastervarsSize);

   assert(SCIPisPositive(masterprob, fractionality1));
   assert(SCIPisPositive(masterprob, fractionality2));
   assert(SCIPisLE(masterprob, fractionality / 2, fractionality1) || SCIPisLE(masterprob, fractionality / 2, fractionality2));
   assert(
      (SCIPisFeasIntegral(masterprob, fractionality1) && SCIPisFeasIntegral(masterprob, fractionality2)) ||
      (!SCIPisFeasIntegral(masterprob, fractionality1) && !SCIPisFeasIntegral(masterprob, fractionality2))
   );
#endif

   // recursive call
   SCIP_CALL( separation_helper(gcg, lowMastervars, lowMastervarsSize, lowCompBndSeq, newCompBndSeqSize, blocknr, chooseCompBnd, foundCompBndSeqList, foundCompBndSeqSizeList, foundCompBndSeqListSize, result) );

   SCIPfreeBlockMemoryArray(masterprob, &lowMastervars, lowMastervarsSize);
   lowMastervarsSize = 0;
   assert(lowMastervars == NULL);

   SCIP_CALL( separation_helper(gcg, highMastervars, highMastervarsSize, highCompBndSeq, newCompBndSeqSize, blocknr, chooseCompBnd, foundCompBndSeqList, foundCompBndSeqSizeList, foundCompBndSeqListSize, result) );

   SCIPfreeBlockMemoryArray(masterprob, &highMastervars, highMastervarsSize);
   highMastervarsSize = 0;
   assert(highMastervars == NULL);

   return SCIP_OKAY;
}

/** collects all currently mastervariables in the given block that satisfy the given component bound sequence */
static
SCIP_RETCODE findSatisfyingMastervars(
   SCIP*                 masterscip,              /**< SCIP data structure */
   SCIP_VAR***           satisfyingMastervars,    /**< mastervariables in the given block satisfying the compBndSeq */
   int*                  satisfyingMastervarsSize,/**< size of mastervars */
   int                   blocknr,                 /**< number of the block */
   GCG_COMPBND*          compBndSeq,              /**< Component Bound Sequence defining the nodes */
   int                   compBndSeqSize           /**< size of compBndSeq */
   )
{
   SCIP_VAR** allMastervars;
   int allMastervarsSize;
   int i;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(satisfyingMastervars != NULL);
   assert(satisfyingMastervarsSize != NULL);

   allMastervars = NULL;
   allMastervarsSize = 0;

   assert(masterscip != NULL);

   SCIP_CALL( SCIPgetVarsData(masterscip, &allMastervars, &allMastervarsSize, NULL, NULL, NULL, NULL) );

   if(*satisfyingMastervarsSize > 0)
   {
      assert(*satisfyingMastervars != NULL);
      SCIPfreeBlockMemoryArray(masterscip, satisfyingMastervars, *satisfyingMastervarsSize);
      *satisfyingMastervarsSize = 0;
   }

   assert(*satisfyingMastervarsSize == 0 && *satisfyingMastervars == NULL);

   for( i = 0; i < allMastervarsSize; ++i )
   {
      if( !GCGisMasterVarInBlock(allMastervars[i], blocknr) )
         continue;

      if( compBndSeqSize > 0 && !isMasterVarInCompBndSeq(masterscip, allMastervars[i], compBndSeq, compBndSeqSize, blocknr) )
            continue;

      if( *satisfyingMastervarsSize == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, satisfyingMastervars, 1) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, satisfyingMastervars, *satisfyingMastervarsSize, *satisfyingMastervarsSize+1) );
      }
      (*satisfyingMastervars)[*satisfyingMastervarsSize] = allMastervars[i];
      (*satisfyingMastervarsSize) += 1;
   }

   return SCIP_OKAY;
}

/* @todo: Vanderbeck (2000) suggests other (also optimal?) algorithms to find a good (best?) component bound sequence; we should check that */

/** choose to separate on original variable where max-min is maximal */
static
CHOOSE_COMPBND(chooseMaxMin)
{
   // Task: Find fractional mastervarMin, mastervarMax, indexed by x1, x2, s.t. for some j, xj1 < xj2
   SCIP_VAR* currentMastervar;
   SCIP_Real currentSolutionValue; // of the currMastervar
   SCIP_VAR* currentOrigvar;
   SCIP_Real currentMin;
   SCIP_Real currentMax;
   SCIP_Real currentDiff;
   SCIP_Real largestDiffMin;
   SCIP_Real largestDiffMax;
   SCIP_Real largestDiff;
   int i;
   int j;
   SCIP_Real generatorentry;
   SCIP_Bool found = FALSE;
   SCIP* masterprob;

#ifndef NDEBUG
   SCIP_Real givenMin;
   SCIP_Real givenMax;
#endif

   masterprob = GCGgetMasterprob(gcg);
   largestDiff = -SCIPinfinity(masterprob);

   for( j=0; j<indexSetSize; ++j )
   {
      currentOrigvar = indexSet[j];
#ifndef NDEBUG
      givenMin = getLowerBound(masterprob, currentOrigvar, compBndSeq, compBndSeqSize);
      givenMax = getUpperBound(masterprob, currentOrigvar, compBndSeq, compBndSeqSize);
      assert(SCIPisFeasLE(masterprob, givenMin, givenMax));
#endif
      currentMin = SCIPinfinity(masterprob);
      currentMax = -SCIPinfinity(masterprob);
      for( i=0; i<satisfyingMastervarsSize; ++i )
      {
         currentMastervar = satisfyingMastervars[i];

         // Current solution value of the master variable must be fractional and > 0
         currentSolutionValue = SCIPgetSolVal(masterprob, NULL, currentMastervar);
         if ( !SCIPisFeasPositive(masterprob, currentSolutionValue) || SCIPisFeasIntegral(masterprob, currentSolutionValue) ) {
            // skip
            continue;
         }

         generatorentry = getGeneratorEntry(currentMastervar, currentOrigvar);

         // first, check if we can update the min and max values
         if( SCIPisInfinity(masterprob, currentMin) )
         {
            assert(SCIPisInfinity(masterprob, -currentMax));
            currentMin = generatorentry;
            currentMax = generatorentry;
         }
         else if( SCIPisLT(masterprob, generatorentry, currentMin) )
            currentMin = generatorentry;
         else if( SCIPisLT(masterprob, currentMax, generatorentry) )
            currentMax = generatorentry;

         // now could be that min<max
         if( SCIPisLT(masterprob, currentMin, currentMax) )
         {
            found = TRUE;
            currentDiff = currentMax - currentMin;
            if( SCIPisGT(masterprob, currentDiff, largestDiff) )
            {
               largestDiffMin = currentMin;
               largestDiffMax = currentMax;
               largestDiff = currentDiff;
               *selectedOrigvar = currentOrigvar;
            }
         }
      }
   }

   if( !found ) // should always be able to find - the branching scheme is sound and complete
   {
      return SCIP_ERROR;
   }

   assert(largestDiffMin < largestDiffMax);
   assert(largestDiff > 0);
   assert(SCIPisFeasLE(masterprob, SCIPvarGetLbGlobal(*selectedOrigvar), largestDiffMin));
   assert(SCIPisFeasLE(masterprob, largestDiffMax, SCIPvarGetUbGlobal(*selectedOrigvar)));
   SCIPdebugMessage("Found %d for origvar %s, j=%d, min=%f, max=%f\n", found, SCIPvarGetName(*selectedOrigvar), j, largestDiffMin, largestDiffMax);

   // j and origivar are still set to the correct values after execution of the loop
   *value = (largestDiffMin+largestDiffMax)/2;

   return SCIP_OKAY;
}

/** choose to separate on original variable which has the largest number of distinct numbers, and find the median value */
static
CHOOSE_COMPBND(chooseMostDistinctValuesMedian)
{
   SCIP* masterprob;
   int selectedNumDistinctValues = 0;
   SCIP_VAR* currentOrigvar;
   int currentNumDistinctValues;
   SCIP_Real* distinctValues = NULL;
   int distinctValuesSize = 0;
   SCIP_Bool inserted;
   SCIP_Bool duplicate;
   int i;
   int j;
   int k;
   int l;

   masterprob = GCGgetMasterprob(gcg);
   assert(*selectedOrigvar == NULL);

   SCIP_Real generatorentry;

   for( j=0; j<indexSetSize; ++j )
   {
      currentOrigvar = indexSet[j];
      currentNumDistinctValues = 0;

      for( i=0; i<satisfyingMastervarsSize; ++i )
      {
         SCIP_VAR* current_mastervar = satisfyingMastervars[i];

         // Current solution value of the master variable must be fractional and > 0
         SCIP_Real current_solution_value = SCIPgetSolVal(masterprob, NULL, current_mastervar);
         if ( !SCIPisFeasPositive(masterprob, current_solution_value) || SCIPisFeasIntegral(masterprob, current_solution_value) ) {
            // skip
            continue;
         }

         generatorentry = getGeneratorEntry(current_mastervar, currentOrigvar);

         // insert sorted, if the value is not already in the array
         if( distinctValuesSize == 0 )
         {
            assert(distinctValues == NULL);
            assert(currentNumDistinctValues == 0);
            distinctValuesSize = SCIPcalcMemGrowSize(masterprob, 1);
            SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &distinctValues, distinctValuesSize) );
            distinctValues[0] = generatorentry;
            currentNumDistinctValues = 1;
         }
         else
         {
            // current_distinctValues is sorted ascendingly
            inserted = FALSE;
            duplicate = FALSE;
            for( k=0; k<currentNumDistinctValues; ++k )
            {
               if( SCIPisEQ(masterprob, generatorentry, distinctValues[k]) )
               {
                  duplicate = TRUE;
                  break;
               }

               if( SCIPisLT(masterprob, generatorentry, distinctValues[k]) )
               {
                  // insert at position k
                  if( currentNumDistinctValues == distinctValuesSize )
                  {
                     int new_current_distinctValuesSize = SCIPcalcMemGrowSize(masterprob, currentNumDistinctValues + 1);
                     SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &distinctValues, distinctValuesSize, new_current_distinctValuesSize) );
                     distinctValuesSize = new_current_distinctValuesSize;
                  }

                  for( l = currentNumDistinctValues; l > k; --l )
                  {
                     distinctValues[l] = distinctValues[l-1];
                  }

                  distinctValues[k] = generatorentry;
                  currentNumDistinctValues += 1;

                  inserted = TRUE;
                  break;
               }
            }

            if( !inserted && !duplicate )
            {
               // insert at the end
               if( currentNumDistinctValues == distinctValuesSize )
               {
                  int new_current_distinctValuesSize = SCIPcalcMemGrowSize(masterprob, currentNumDistinctValues + 1);
                  SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &distinctValues, distinctValuesSize, new_current_distinctValuesSize) );
                  distinctValuesSize = new_current_distinctValuesSize;
               }

               distinctValues[currentNumDistinctValues] = generatorentry;
               currentNumDistinctValues += 1;
            }
         }
      }

      if( currentNumDistinctValues > 1 && currentNumDistinctValues > selectedNumDistinctValues )
      {
         assert(currentNumDistinctValues > 1);
         *selectedOrigvar = currentOrigvar;
         selectedNumDistinctValues = currentNumDistinctValues;

         // calculate median
         if( currentNumDistinctValues % 2 == 0 )
         {
            *value = 0.5 * (distinctValues[currentNumDistinctValues / 2] + distinctValues[currentNumDistinctValues / 2 - 1]);
         }
         else
         {
            *value = distinctValues[currentNumDistinctValues / 2];
         }
      }
   }

   assert(*selectedOrigvar != NULL);

   SCIPfreeBlockMemoryArrayNull(masterprob, &distinctValues, distinctValuesSize);

   return SCIP_OKAY;
}

/** calculate sum of component bound sequence */
static
SCIP_Real calcSum(
   SCIP*                 masterscip,              /**< SCIP data structure */
   GCG_COMPBND*          compBndSeq,              /**< Component Bound Sequence defining the nodes */
   int                   compBndSeqSize,                   /**< size of compBndSeq */
   int                   blocknr                  /**< number of the block */
   )
{
   int i;
   SCIP_Real sum = 0.0;
   SCIP_VAR** satisfyingMastervars = NULL;
   int satisfyingMastervarsSize = 0;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));

   assert(satisfyingMastervars == NULL && satisfyingMastervarsSize == 0);
   findSatisfyingMastervars(masterscip, &satisfyingMastervars, &satisfyingMastervarsSize, blocknr, compBndSeq, compBndSeqSize);
   assert((satisfyingMastervarsSize > 0) == (satisfyingMastervars != NULL));

   sum = 0;
   for( i = 0; i < satisfyingMastervarsSize; ++i )
   {
      sum += SCIPgetSolVal(masterscip, NULL, satisfyingMastervars[i]);
   }

   SCIPfreeBlockMemoryArrayNull(masterscip, &satisfyingMastervars, satisfyingMastervarsSize);

   return sum;
}

/** choose the component bound of which the sum is closest to K/2, where K is the number of identical subproblems in the current block */
static
CHOOSE_COMPBNDSEQ(chooseClosestToKHalf)
{
   SCIP* masterprob;
   int i;
   SCIP_Real sum;
   SCIP_Real kHalf;

   int bestIndex = -1;
   SCIP_Real best;
   
   masterprob = GCGgetMasterprob(gcg);

   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));
   assert(compBndSeqList != NULL);
   assert(compBndSeqSizeList != NULL);
   assert(compBndListSize > 0);

   best = SCIPinfinity(masterprob);
   kHalf = GCGgetNIdenticalBlocks(gcg, blocknr) * 0.5;

   for( i = 0; i < compBndListSize; ++i )
   {
      sum = ABS(calcSum(masterprob, compBndSeqList[i], compBndSeqSizeList[i], blocknr) - kHalf);

      if( sum < best )
      {
         best = sum;
         bestIndex = i;
      }
   }

   assert(bestIndex >= 0);

   *selectedCompBndSeq = compBndSeqList[bestIndex];
   *selectedCompBndSeqSize = compBndSeqSizeList[bestIndex];

   return SCIP_OKAY;
}

/** choose the component bound of which the sum is most fractional */
static
CHOOSE_COMPBNDSEQ(chooseMostFractional)
{
   int i;
   SCIP_Real sum;
   SCIP_Real fractionality;
   SCIP* masterprob;

   int bestIndex = -1;
   SCIP_Real best = 0.0;

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));
   assert(compBndSeqList != NULL);
   assert(compBndSeqSizeList != NULL);
   assert(compBndListSize > 0);

   for( i = 0; i < compBndListSize; ++i )
   {
      sum = calcSum(masterprob, compBndSeqList[i], compBndSeqSizeList[i], blocknr);
      fractionality = SCIPfrac(masterprob, sum);

      fractionality = MIN(fractionality, 1.0 - fractionality);

      if( best < fractionality )
      {
         best = fractionality;
         bestIndex = i;
      }
   }

   assert(bestIndex >= 0);

   *selectedCompBndSeq = compBndSeqList[bestIndex];
   *selectedCompBndSeqSize = compBndSeqSizeList[bestIndex];

   return SCIP_OKAY;
}

/** choose the component bound sequence with the smallest size */
static
CHOOSE_COMPBNDSEQ(chooseSmallestSize)
{
   int i;
   int smallestSize = INT_MAX;
   SCIP* masterprob;

   GCG_COMPBND** newCompBndSeqList = NULL;
   int* newCompBndSeqSizeList = NULL;
   int newCompBndSeqListSize = 0;

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);
   assert(GCGisMaster(masterprob));
   assert(compBndSeqList != NULL);
   assert(compBndSeqSizeList != NULL);
   assert(compBndListSize > 0);

   // step 1: find smallest size
   for( i = 0; i < compBndListSize; ++i )
   {
      if( compBndSeqSizeList[i] < smallestSize )
      {
         smallestSize = compBndSeqSizeList[i];
      }
   }
   assert(smallestSize > 0);

   // step 2: create new Blist with smallest size
   for( i = 0; i < compBndListSize; ++i )
   {
      if( compBndSeqSizeList[i] == smallestSize )
      {
         if( newCompBndSeqListSize == 0 )
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &newCompBndSeqList, 1) );
            SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &newCompBndSeqSizeList, 1) );
         }
         else
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &newCompBndSeqList, newCompBndSeqListSize, newCompBndSeqListSize + 1) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &newCompBndSeqSizeList, newCompBndSeqListSize, newCompBndSeqListSize + 1) );
         }

         newCompBndSeqList[newCompBndSeqListSize] = compBndSeqList[i];
         newCompBndSeqSizeList[newCompBndSeqListSize] = compBndSeqSizeList[i];
         newCompBndSeqListSize += 1;
      }
   }

   assert(newCompBndSeqListSize > 0);

   // step 3: refine the selection
   //SCIP_CALL( chooseClosestToKHalf(masterscip, newCompBndSeqList, newCompBndSeqSizeList, newCompBndSeqListSize, blocknr, selectedCompBndSeq, selectedCompBndSeqSize) );
   SCIP_CALL( chooseMostFractional(gcg, newCompBndSeqList, newCompBndSeqSizeList, newCompBndSeqListSize, blocknr, selectedCompBndSeq, selectedCompBndSeqSize) );

   // step 4: free memory
   SCIPfreeBlockMemoryArray(masterprob, &newCompBndSeqList, newCompBndSeqListSize);
   SCIPfreeBlockMemoryArray(masterprob, &newCompBndSeqSizeList, newCompBndSeqListSize);

   return SCIP_OKAY;
}

/** separation algorithm determining a component bound sequence to branch on */
static
SCIP_RETCODE separation(
   GCG*                  gcg,                     /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,              /**< branching rule */
   GCG_COMPBND**         compBndSeq,              /**< Component Bound Sequence defining the nodes */
   int*                  compBndSeqSize,                   /**< size of compBndSeq */
   int                   blocknr,                 /**< number of the block */
   SCIP_RESULT*          result                   /**< pointer to store the result of the branching call */
   )
{
   // a list of component bound sequences (also lists)
   GCG_COMPBND** compBndSeqList;
   int* compBndSeqSizeList;
   int compBndSeqListSize;
   SCIP* masterprob;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** satisfyingMastervars;
   int satisfyingMastervarsSize;
   int i;

   compBndSeqList = NULL;
   compBndSeqSizeList = NULL;
   compBndSeqListSize = 0;

   satisfyingMastervars = NULL;
   satisfyingMastervarsSize = 0;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   masterprob = GCGgetMasterprob(gcg);

   findSatisfyingMastervars(masterprob, &satisfyingMastervars, &satisfyingMastervarsSize, blocknr, NULL, 0);
#ifndef NDEBUG
   int satisfyingMastervarsSizeCopy = satisfyingMastervarsSize;
   int BlistsizeCopy = compBndSeqListSize;
#endif

   // call the separation algorithm with all chooseVars functions
   assert(branchruledata->useMaxRangeMidrangeHeuristic || branchruledata->useMostDistinctMedianHeuristic);

   if( branchruledata->useMaxRangeMidrangeHeuristic )
   {
      SCIP_CALL( separation_helper(gcg, satisfyingMastervars, satisfyingMastervarsSize, NULL, 0, blocknr, chooseMaxMin, &compBndSeqList, &compBndSeqSizeList, &compBndSeqListSize, result) );
#ifndef NDEBUG
      assert(satisfyingMastervarsSize == satisfyingMastervarsSizeCopy);
      assert(compBndSeqListSize > BlistsizeCopy);
      BlistsizeCopy = compBndSeqListSize;
#endif
   }

   if( branchruledata->useMostDistinctMedianHeuristic )
   {
      SCIP_CALL( separation_helper(gcg, satisfyingMastervars, satisfyingMastervarsSize, NULL, 0, blocknr, chooseMostDistinctValuesMedian, &compBndSeqList, &compBndSeqSizeList, &compBndSeqListSize, result) );
#ifndef NDEBUG
      assert(satisfyingMastervarsSize == satisfyingMastervarsSizeCopy);
      assert(compBndSeqListSize > BlistsizeCopy);
      BlistsizeCopy = compBndSeqListSize;
#endif
   }

   SCIPfreeBlockMemoryArray(masterprob, &satisfyingMastervars, satisfyingMastervarsSize);
   satisfyingMastervarsSize = 0;
   assert(satisfyingMastervars == NULL);

   assert(compBndSeqListSize > 0);
   assert(compBndSeqList != NULL);
   assert(compBndSeqSizeList != NULL);

   SCIP_CALL( chooseSmallestSize(gcg, compBndSeqList, compBndSeqSizeList, compBndSeqListSize, blocknr, compBndSeq, compBndSeqSize) );

   assert(*compBndSeqSize > 0);
   assert(*compBndSeq != NULL);

   // free others, as well as Blist
   for( i = 0; i < compBndSeqListSize; ++i )
   {
      if( compBndSeqList[i] != *compBndSeq )
      {
         freeComponentBoundSequence(masterprob, &compBndSeqList[i], &compBndSeqSizeList[i]);
      }
   }

   SCIPfreeBlockMemoryArray(masterprob, &compBndSeqList, compBndSeqListSize);
   SCIPfreeBlockMemoryArray(masterprob, &compBndSeqSizeList, compBndSeqListSize);

   return SCIP_OKAY;
}

/** prepares information for using the generic branching scheme */
static
SCIP_RETCODE initBranch(
   GCG*                  gcg,                     /**< GCG data structure */
   SCIP_BRANCHRULE*      branchrule,              /**< branching rule */
   SCIP_RESULT*          result                   /**< pointer to store the result of the branching call */
   )
{
#ifdef COMPBND_STATS
   SCIP_BRANCHRULEDATA* branchruledata;
#endif
   SCIP* masterprob;
   SCIP_VAR** branchcands;
   int nbranchcands;
   GCG_COMPBND* compBndSeq;
   int compBndSeqSize;
   int blocknr;
   SCIP_Bool foundBlock;
   SCIP_VAR** satisfyingMastervars;
   int satisfyingMastervarsSize;
   SCIP_Real fractionality;

   blocknr = -2;
   compBndSeq = NULL;
   compBndSeqSize = 0;

   masterprob = GCGgetMasterprob(gcg);

   assert(masterprob != NULL);

   SCIPdebugMessage("get information for component bound branching\n");

   SCIP_CALL( SCIPgetLPBranchCands(masterprob, &branchcands, NULL, NULL, &nbranchcands, NULL, NULL) );

#ifdef COMPBND_STATS
   branchruledata = SCIPbranchruleGetData(branchrule);
#endif

   /* in case original problem contains continuous variables, there are no branching cands */
   assert(nbranchcands > 0);

   foundBlock = FALSE;
   satisfyingMastervarsSize = 0;
   satisfyingMastervars = NULL;

   /* 1. Determine in what block we are branching. We select the first available block,
    *     i.e. the first block that contains a branching candidate, starting from the master block.
    */
   for( blocknr = 0; blocknr < GCGgetNPricingprobs(gcg); blocknr++ )
   {
      if( !GCGisPricingprobRelevant(gcg, blocknr) )
         continue;

      SCIPdebugMessage("\nTrying to branch in block %d:\n", blocknr);

      /* 2. Check whether the fractionality of all master variables in this block is not 0. */
      findSatisfyingMastervars(masterprob, &satisfyingMastervars, &satisfyingMastervarsSize, blocknr, NULL, 0);
      assert((satisfyingMastervarsSize > 0) == (satisfyingMastervars != NULL));
      fractionality = calcFractionality(masterprob, satisfyingMastervars, satisfyingMastervarsSize);
      assert(SCIPisFeasIntegral(masterprob, fractionality));
      SCIPfreeBlockMemoryArray(masterprob, &satisfyingMastervars, satisfyingMastervarsSize);
      satisfyingMastervarsSize = 0;
      assert(satisfyingMastervars == NULL);
      if( SCIPisZero(masterprob, fractionality) )
      {
         SCIPdebugMessage("No fractional integer variables in block %d\n", blocknr);
         continue;
      }

      /* 3. Call to separation algorithm to find a suitable compBndSeq to branch on in the current block. */
      SCIP_CALL( separation(gcg, branchrule, &compBndSeq, &compBndSeqSize, blocknr, result) );

      if( *result == SCIP_BRANCHED )
      {
         assert(compBndSeqSize > 0);
         assert(compBndSeq != NULL);
         SCIPdebugMessage("Branching in block %d\n", blocknr);
         foundBlock = TRUE;
         break;
      }
      else
      {
         if( compBndSeqSize > 0 )
         {
            freeComponentBoundSequence(masterprob, &compBndSeq, &compBndSeqSize);
            compBndSeqSize = 0;
         }
      }
   }

   if( !foundBlock )
   {
      assert(compBndSeqSize == 0);
      assert(compBndSeq == NULL);
      SCIPdebugMessage("No block found to branch on\n");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Component bound sequence to branch on:\n");
   SCIPdebugMessage("compBndSeqSize: %d\n", compBndSeqSize);
   for (int i = 0; i < compBndSeqSize; ++i)
   {
      SCIPdebugMessage("compBndSeq[%d]: %s %s %d\n", i, SCIPvarGetName((compBndSeq)[i].component),
                        (compBndSeq)[i].sense == GCG_COMPBND_SENSE_LE ? "<=" : ">=",
                        (compBndSeq)[i].bound);
   }

   /* 4. Remove component bounds that are strengthened by others */
   SCIP_CALL( simplifyComponentBoundSequence(masterprob, &compBndSeq, &compBndSeqSize) );
   assert(compBndSeqSize > 0);
   assert(compBndSeq != NULL);

#ifdef COMPBND_STATS
   int depth = SCIPgetDepth(masterprob);
   SCIP_Real sum = calcSum(masterprob, compBndSeq, compBndSeqSize, blocknr);
   int K = GCGgetNIdenticalBlocks(origprob, blocknr);
   fprintf(branchruledata->statsfile, "%d,%d,%f,%d\n", depth, compBndSeqSize, sum, K);
#endif

   /* 5. Create the child nodes. */
   SCIP_CALL( createChildNodesCompBndSeq(gcg, branchrule, compBndSeq, compBndSeqSize, blocknr) );

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeCompBnd)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);

#ifdef COMPBND_STATS
   if( branchruledata->statsfile != NULL )
   {
      fclose(branchruledata->statsfile);
      branchruledata->statsfile = NULL;
   }
#endif

   SCIPfreeBlockMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitCompBnd NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolCompBnd NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolCompBnd NULL


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpCompBnd)
{  /*lint --e{715}*/
   SCIP* origprob;
   SCIP_Bool discretization;
   SCIP_BRANCHRULEDATA* branchruledata;
   
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   origprob = GCGgetOrigprob(branchruledata->gcg);
   assert(origprob != NULL);

   SCIPdebugMessage("Execlp method of component bound branching\n");

   *result = SCIP_DIDNOTRUN;

   /* the branching scheme only works for the discretization approach
    * actually, this branching scheme _enforces_ a discretization concept
    */
   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( !discretization )
   {
      SCIPdebugMessage("Component bound branching only for discretization approach\n");
      return SCIP_OKAY;
   }

   if( GCGisMasterSetCovering(branchruledata->gcg) || GCGisMasterSetPartitioning(branchruledata->gcg) )
   {
      SCIPdebugMessage("Component bound branching executed on a set covering or set partitioning problem\n");
   }

   if( GCGrelaxIsOrigSolFeasible(branchruledata->gcg) )
   {
      SCIPdebugMessage("node cut off, since origsol was feasible, solval = %f\n",
         SCIPgetSolOrigObj(origprob, GCGrelaxGetCurrentOrigSol(branchruledata->gcg)));

      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   *result = SCIP_BRANCHED;

   SCIP_CALL( initBranch(branchruledata->gcg, branchrule, result) );

   return SCIP_OKAY;
}


static
SCIP_DECL_BRANCHEXECEXT(branchExecextCompBnd)
{  /*lint --e{715}*/
   SCIPdebugMessage("Execext method of component bound branching\n");

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}


/** copy method for branchrule plugins (called when SCIP copies plugins) */
#define branchCopyCompBnd NULL


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsCompBnd)
{  /*lint --e{715}*/

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of component bound branching\n");

   if( SCIPisStopped(scip) )
   {
      SCIPwarningMessage(scip, "No branching could be created, solving process cannot be restarted...\n" );

      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   else
   {
      SCIPerrorMessage("This method is not implemented, aborting since we cannot recover!\n");
      SCIPdialogMessage(scip, NULL, "Due to numerical issues, the problem could not be solved.\n");
      SCIPdialogMessage(scip, NULL, "You can try to disable discretization and aggregation and resolve the problem.\n");

      *result = SCIP_DIDNOTRUN;
      return SCIP_ERROR;
   }
}

/*
 * GCG specific branching rule callbacks
 */

/** activation method for branchrule, called when a node in the master problem is activated,
 *  should perform changes to the current node's problem due to the branchdata
 */
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterCompBnd)
{
   assert(gcg != NULL);

   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   SCIPdebugMessage("branchActiveMasterCompBnd: Block %d, Ssize %d\n", branchdata->blocknr, branchdata->compBndSeqSize);

   return SCIP_OKAY;
}


/** deactivation method for branchrule, called when a node in the master problem is deactivated,
 *  should undo changes to the current node's problem due to the branchdata
 */
static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterCompBnd)
{
   assert(gcg != NULL);

   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   SCIPdebugMessage("branchDeactiveMasterCompBnd: Block %d, Ssize %d\n", branchdata->blocknr, branchdata->compBndSeqSize);

   return SCIP_OKAY;
}

/** propagation method for branchrule, called when a node in the master problem is propagated,
 *  should perform propagation at the current node due to the branchdata
 */
static
GCG_DECL_BRANCHPROPMASTER(branchPropMasterCompBnd)
{
   assert(gcg != NULL);
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);
   assert(branchdata->compBndSeq != NULL);

   *result = SCIP_DIDNOTFIND;
   return SCIP_OKAY;
}

/** method for branchrule, called when the master LP is solved at one node,
 *  can store pseudocosts for the branching decisions
 */
#define branchMasterSolvedCompBnd NULL

/** frees branching data of an origbranch constraint (called when the origbranch constraint is deleted) */
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteCompBnd)
{
   SCIP* masterscip;

   assert(gcg != NULL);
   assert(branchdata != NULL);

   if( origbranch && !force )
      return SCIP_OKAY;

   masterscip = GCGgetMasterprob(gcg);
   assert(masterscip != NULL);

   SCIPdebugMessage("branchDataDeleteCompBnd: Block %d, Ssize %d\n", (*branchdata)->blocknr, (*branchdata)->compBndSeqSize);

   assert(*branchdata != NULL);

   /* release constraint that enforces the branching decision */
   assert((*branchdata)->mastercons != NULL);
   SCIP_CALL( GCGbranchFreeExtendedmasterconsBranchData(gcg, &(*branchdata)->mastercons) );
   assert((*branchdata)->mastercons == NULL);

   assert((*branchdata)->compBndSeq != NULL && (*branchdata)->compBndSeqSize > 0);
   SCIPfreeBlockMemoryArray(masterscip, &((*branchdata)->compBndSeq), (*branchdata)->compBndSeqSize);
   assert((*branchdata)->compBndSeq == NULL);
   (*branchdata)->compBndSeqSize = 0;

   SCIPfreeBlockMemory(masterscip, branchdata);
   assert(*branchdata == NULL);

   return SCIP_OKAY;
}

/** callback new column method */
static
GCG_DECL_BRANCHNEWCOL(branchNewColCompBnd)
{
   assert(gcg != NULL);

   assert(mastervar != NULL);
   assert(GCGvarIsMaster(mastervar));
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   SCIP_Bool added = FALSE;
   SCIP_CALL( addVarToMasterbranch(gcg, mastervar, branchdata, &added) );

   return SCIP_OKAY;
}

static
GCG_DECL_BRANCHGETEXTENDEDMASTERCONS(branchGetExtendedmasterconsCompBnd)
{
   assert(gcg != NULL);

   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   *extendedmasterconsdata = branchdata->mastercons;

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitCompBnd)
{
#ifdef COMPBND_STATS
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_Bool stabilization_tree;

   assert(branchrule != NULL);
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   SCIP_CALL( SCIPgetBoolParam(origprob, "pricing/masterpricer/stabilizationtree", &stabilization_tree) );
#endif

   SCIPdebugMessage("Init method of component bound branching\n");

#ifdef COMPBND_STATS
   /* determine the filename */
   // "stats_compbnd" + (if useMaxRangeMidrangeHeuristic && useMostDistinctMedianHeuristic ? "" : (if useMaxRangeMidrangeHeuristic ? "-mrm" : "-mdm")) + (if stabilization_tree ? "+dvs" : "") + ".csv"

   char filename[SCIP_MAXSTRLEN];
   (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "stats_compbnd%s%s.csv",
      branchruledata->useMaxRangeMidrangeHeuristic && branchruledata->useMostDistinctMedianHeuristic ? "" :
      (branchruledata->useMaxRangeMidrangeHeuristic ? "-mrm" : "-mdm"),
      stabilization_tree ? "+dvs" : "");

   /* create the file with header, if not exists */
   branchruledata->statsfile = fopen(filename, "r");
   if( branchruledata->statsfile == NULL )
   {
      branchruledata->statsfile = fopen(filename, "w");
      assert(branchruledata->statsfile != NULL);
      fprintf(branchruledata->statsfile, "depth,num_bounds,sum,K\n");
   }
   else
   {
      fclose(branchruledata->statsfile);
      branchruledata->statsfile = fopen(filename, "a");
      assert(branchruledata->statsfile != NULL);
   }
#endif

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the compbnd branching rule and includes it in SCIP */
SCIP_RETCODE GCGincludeBranchruleCompBnd(
   GCG*                  gcg                 /**< GCG data structure */
)
{
   SCIP* masterprob;
   SCIP* origprob;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   SCIP_UNUSED(chooseClosestToKHalf);

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);
   origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   /* create compbnd branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(masterprob, &branchruledata) );
   branchruledata->gcg = gcg;
   branchruledata->gcgbranchrule = NULL;

   SCIPdebugMessage("Include method of component bound branching\n");

   /* include branching rule */
   SCIP_CALL( GCGrelaxIncludeBranchrule(branchruledata->gcg, &branchrule, &branchruledata->gcgbranchrule, BRANCHRULE_NAME,
         BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata, branchActiveMasterCompBnd,
         branchDeactiveMasterCompBnd, branchPropMasterCompBnd, branchMasterSolvedCompBnd, branchDataDeleteCompBnd,
         branchNewColCompBnd, branchGetExtendedmasterconsCompBnd, branchGetExtendedmasterconsCoeffCompBnd) );
   SCIP_CALL( SCIPsetBranchruleInit(masterprob, branchrule, branchInitCompBnd) );
   SCIP_CALL( SCIPsetBranchruleFree(masterprob, branchrule, branchFreeCompBnd) );
   SCIP_CALL( SCIPsetBranchruleExecLp(masterprob, branchrule, branchExeclpCompBnd) );
   SCIP_CALL( SCIPsetBranchruleExecExt(masterprob, branchrule, branchExecextCompBnd) );
   SCIP_CALL( SCIPsetBranchruleExecPs(masterprob, branchrule, branchExecpsCompBnd) );

   /* add compbnd branching rule parameters */
   SCIP_CALL( SCIPaddBoolParam(origprob, "branching/compbnd/useMRMH",
      "should the max range midrange heuristic be used?", &branchruledata->useMaxRangeMidrangeHeuristic, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(origprob, "branching/compbnd/useMDMH",
      "should the most distinct median heuristic be used?", &branchruledata->useMostDistinctMedianHeuristic, FALSE, TRUE, NULL, NULL) );

   branchrule = SCIPfindBranchrule(masterprob, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   return SCIP_OKAY;
}

/** returns true when the branch rule is the component bound branching rule */
SCIP_Bool GCGisBranchruleCompBnd(
   SCIP_BRANCHRULE*      branchrule          /**< branchrule to check */
)
{
   return (branchrule != NULL) && (strcmp(BRANCHRULE_NAME, SCIPbranchruleGetName(branchrule)) == 0);
}
