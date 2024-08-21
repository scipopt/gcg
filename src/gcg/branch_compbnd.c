/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file    branch_compbnd.h
 *
 * @brief   component bound branching rule
 * @author  Til Mohr
 *
 * This is an implementation of the component bound branching rule based on the papers:
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
//#define STATS

#include <gcg/pricer_gcg.h>
#include <scip/pub_message.h>
#include <scip/pub_var.h>
#include <scip/type_branch.h>
#include <scip/type_misc.h>
#include <assert.h>
#include <scip/def.h>
#include <scip/cons_linear.h>
#include <scip/cons_varbound.h>
#include <scip/scip.h>
#include <scip/type_result.h>
#include <scip/type_retcode.h>
#include <scip/type_scip.h>
#include <scip/type_var.h>
#include <string.h>
#include <scip/struct_var.h>
#include <scip/pub_lp.h>
#include <scip/pub_misc_linear.h>
#include <scip/scip_mem.h>
#include <scip/scip_numerics.h>
#include <scip/scip_var.h>

#include "branch_compbnd.h"
#include "cons_integralorig.h"
#include "cons_masterbranch.h"
#include "gcg.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "mastercutdata.h"
#include "pub_gcgvar.h"
#include "type_mastercutdata.h"
#include "type_branchgcg.h"

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
   SCIP_Bool             useMaxRangeMidrangeHeuristic; /** should the MaxRangeMidrangeHeuristic be used */
   SCIP_Bool             useMostDistinctMedianHeuristic; /** should the MostDistinctMedianHeuristic be used */
#ifdef STATS
   FILE*                 statsfile;          /**< file to write statistics to */
#endif
};

/** branching data */
struct GCG_BranchData
{
   GCG_BRANCH_TYPE       branchtype;         /**< type of branching */
   int                   constant;           /**< constant value of the branching constraint in the master problem - either lhs or rhs, depending on branchtype */
   GCG_COMPBND*          B;                  /**< component bound sequence which induce the current branching constraint */
   int                   Bsize;              /**< size of the component bound sequence B */
   int                   blocknr;            /**< id of the pricing problem (or block) to which this branching constraint belongs */
   int                   nvars;              /**< number of master variables the last time the node has been visited - neccessary to later include newly generated master variables */
   GCG_MASTERCUTDATA*    mastercons;         /**< master constraint along with its corresponding inferred pricing modifications */
};



/** define the signature of functions to choose the component bound */
#define CHOOSE_COMPBND(x) SCIP_RETCODE x (SCIP* masterscip, SCIP_VAR** X, int Xsize, SCIP_VAR** indexSet, int indexSetSize, GCG_COMPBND* B, int Bsize, int blocknr, SCIP_VAR** selected_origvar, SCIP_Real* value)

/** define the signature of functions that chose one component bound sequence out of many */
#define CHOOSE_COMPBNDSEQ(x) SCIP_RETCODE x (SCIP* masterscip, GCG_COMPBND** Blist, int* Bsizes, int Blistsize, int blocknr, GCG_COMPBND** selected_B, int* selected_Bsize)

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
   GCG_COMPBND**         B,                  /**< component bound sequence to set to NULL */
   int*                  Bsize               /**< size of the component bound sequence B */
   )
{
   assert(scip != NULL);
   assert(GCGisMaster(scip));

   SCIPfreeBlockMemoryArray(scip, B, *Bsize);
   assert(*B == NULL);
   *Bsize = 0;
}

/** initialize branchdata at the node */
static
SCIP_RETCODE initNodeBranchdata(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_BRANCHDATA**      nodebranchdata,     /**< branching data to set */
   GCG_BRANCH_TYPE       branchtype,         /**< type of branch to generate */
   int                   constant,           /**< constant value of the branching constraint in the master problem - either lhs or rhs, depending on branchtype */
   GCG_COMPBND*          B,                  /**< component bound sequence which induce the current branching constraint */
   int                   Bsize,              /**< size of the component bound sequence B */
   int                   blocknr             /**< block we are branching in */
   )
{
   assert(GCGisMaster(scip));

   SCIP_CALL( SCIPallocBlockMemory(scip, nodebranchdata) );

   (*nodebranchdata)->branchtype = branchtype;
   (*nodebranchdata)->blocknr = blocknr;
   (*nodebranchdata)->constant = constant;
   (*nodebranchdata)->Bsize = Bsize;
   (*nodebranchdata)->nvars = 0;
   (*nodebranchdata)->mastercons = NULL;

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &((*nodebranchdata)->B), B, Bsize));

   return SCIP_OKAY;
}

/** simplify B by strengthening bounds of the same sense on the same component */
static
SCIP_RETCODE simplifyComponentBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_COMPBND**         B,                  /**< Component Bound Sequence defining the nodes */
   int*                  Bsize               /**< size of B */
   )
{
   GCG_COMPBND* newB;
   int newBsize;
   SCIP_Bool* already_added;
   int i;
   int j;

   assert(GCGisMaster(scip));

   newB = NULL;
   newBsize = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &already_added, *Bsize) );
   for( i = 0; i < *Bsize; ++i )
   {
      already_added[i] = FALSE;
   }


   for( i=0; i<*Bsize; ++i )
   {
      if( already_added[i] )
         continue;

      SCIP_VAR* component = (*B)[i].component;
      GCG_COMPBND_SENSE sense = (*B)[i].sense;
      SCIP_Real bound = (*B)[i].bound;

      for( j=i+1; j<*Bsize; ++j )
      {
         if( SCIPvarCompare(component, (*B)[j].component) != 0 )
            continue;
         if( (*B)[i].sense != (*B)[j].sense )
            continue;

         already_added[j] = TRUE;

         if( (*B)[i].sense == GCG_COMPBND_SENSE_LE )
         {
            bound = MIN((*B)[j].bound, bound);
         }
         else if( (*B)[i].sense == GCG_COMPBND_SENSE_GE )
         {
            bound = MAX((*B)[j].bound, bound);
         }
      }

      already_added[i] = TRUE;

      if( newBsize == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newB, 1) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &newB, newBsize, newBsize+1) );
      }

      newB[newBsize].component = component;
      newB[newBsize].sense = sense;
      newB[newBsize].bound = bound;
      newBsize += 1;
   }

   assert(1 <= newBsize && newBsize <= *Bsize);
   SCIPdebugMessage("Simplified B from %d to %d\n", *Bsize, newBsize);
   for (i = 0; i < newBsize; ++i)
   {
      SCIPdebugMessage("B[%d]: %s %s %d\n", i, SCIPvarGetName(newB[i].component),
                        newB[i].sense == GCG_COMPBND_SENSE_LE ? "<=" : ">=",
                        newB[i].bound);
   }

   SCIPfreeBlockMemoryArray(scip, &already_added, *Bsize);
   assert(already_added == NULL);

   // free old B and replace it with newB
   freeComponentBoundSequence(scip, B, Bsize);
   assert(*B == NULL);
   *B = newB;
   *Bsize = newBsize;

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

      if( SCIPvarCompare(solvars[i], cmp) == 0 )
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
   SCIP_VAR** origvars;
   SCIP_Real* origvals;
   int norigvars;

   assert(mastervar != NULL);
   assert(origvar != NULL);

   origvars = GCGmasterVarGetOrigvars(mastervar);
   norigvars = GCGmasterVarGetNOrigvars(mastervar);
   origvals = GCGmasterVarGetOrigvals(mastervar);

   return getGeneratorEntryCol(origvars, origvals, norigvars, origvar);
}

/** whether a master variable is in B or not */
static
SCIP_Bool isColInB(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            solvars,            /**< column solution variables */
   SCIP_Real*            solvals,            /**< column solution values */
   int                   nsolvars,           /**< number of column solution variables */
   GCG_COMPBND*          B,                  /**< component bound sequence to check */
   int                   Bsize,              /**< size of B */
   int                   blocknr             /**< block we are branching in */
   )
{
   int i;
   SCIP_Real generatorentry;

   assert(Bsize > 0);
   assert(B != NULL);

   for( i = 0; i < Bsize; ++i )
   {
      generatorentry = getGeneratorEntryCol(solvars, solvals, nsolvars, B[i].component);
      if ( (B[i].sense == GCG_COMPBND_SENSE_GE && !SCIPisGE(scip, generatorentry, B[i].bound)) ||
         (B[i].sense == GCG_COMPBND_SENSE_LE && !SCIPisLE(scip, generatorentry, B[i].bound)) )
      {
         return FALSE;
      }
   }

   return TRUE;
}

/** whether a master variable is in B or not */
static
SCIP_Bool isMasterVarInB(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             mastervar,          /**< master variable to check */
   GCG_COMPBND*          B,                  /**< component bound sequence to check */
   int                   Bsize,              /**< size of B */
   int                   blocknr             /**< block we are branching in */
   )
{
   int i;
   SCIP_Real generatorentry;

   assert(mastervar != NULL);
   assert(Bsize > 0);
   assert(B != NULL);

   if( !GCGisMasterVarInBlock(mastervar, blocknr) )
      return FALSE;

   for( i = 0; i < Bsize; ++i )
   {
      generatorentry = getGeneratorEntry(mastervar, B[i].component);
      if ( (B[i].sense == GCG_COMPBND_SENSE_GE && !SCIPisGE(scip, generatorentry, B[i].bound)) ||
         (B[i].sense == GCG_COMPBND_SENSE_LE && !SCIPisLE(scip, generatorentry, B[i].bound)) )
      {
         return FALSE;
      }
   }

   return TRUE;
}

/** adds a variable to a branching constraint */
static
SCIP_RETCODE addVarToMasterbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             mastervar,          /**< the variable to add */
   GCG_BRANCHDATA*       branchdata,         /**< branching data structure where the variable should be added */
   SCIP_Bool*            added               /**< whether the variable was added */
)
{
   SCIP_Real coef;
   SCIP_VAR** origvars;
   SCIP_Real* origvals;
   int norigvars;

   assert(scip != NULL);
   assert(mastervar != NULL);
   assert(branchdata != NULL);
   assert(added != NULL);

   origvars = GCGmasterVarGetOrigvars(mastervar);
   norigvars = GCGmasterVarGetNOrigvars(mastervar);
   origvals = GCGmasterVarGetOrigvals(mastervar);

   *added = FALSE;

   coef = GCGmastercutGetCoeff(scip, branchdata->mastercons, origvars, origvals, norigvars, GCGvarGetBlock(mastervar));

   if( !SCIPisZero(scip, coef) )
   {
      SCIPdebugMessage("Adding variable %s to branching constraint: %s %d\n", SCIPvarGetName(mastervar), branchdata->branchtype == GCG_BRANCH_DOWN ? "<=" : ">=", branchdata->constant);
      SCIP_CALL( GCGmastercutAddMasterVar(scip, branchdata->mastercons, mastervar, coef) );
      *added = TRUE;
   }

   return SCIP_OKAY;
}

/** callback new column method */
static
GCG_DECL_MASTERCUTGETCOEFF(mastercutGetCoeffCompBnd)
{
   GCG_BRANCHDATA* branchdata;

   assert(scip != NULL);
   assert(GCGisMaster(scip));

   branchdata = (GCG_BRANCHDATA*) GCGmastercutGetData(mastercutdata);

   assert(branchdata != NULL);
   assert(branchdata->mastercons == mastercutdata);

   *coef = 0.0;

   if( probnr == -1 || probnr != branchdata->blocknr )
      return SCIP_OKAY;

   if( !isColInB(scip, solvars, solvals, nsolvars, branchdata->B, branchdata->Bsize, branchdata->blocknr) )
      return SCIP_OKAY;

   *coef = 1.0;

   return SCIP_OKAY;
}

/** creates the constraint for branching directly on a master variable */
static
SCIP_RETCODE createBranchingCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to add constraint */
   GCG_BRANCHDATA*       branchdata          /**< branching data structure */
)
{
   SCIP_Longint uuid;

   SCIP* origscip;
   SCIP* pricingscip;

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
   GCG_PRICINGMODIFICATION* pricingmods; // will always only contain one element in this branching rule

   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(node != NULL);
   assert(branchdata != NULL);
   assert(branchdata->mastercons == NULL);

   uuid = SCIPnodeGetNumber(node);

   mastervars = SCIPgetVars(scip);
   nmastervars = SCIPgetNVars(scip);
   assert(nmastervars >= 0);
   assert(nmastervars == 0 || mastervars != NULL);

   /*  create constraint for child */
   if (branchdata->branchtype == GCG_BRANCH_DOWN)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "compbnd_%d_child(%d_LE_%d)", uuid, branchdata->Bsize, branchdata->constant);
      SCIP_CALL( SCIPcreateConsLinear(scip, &branchcons, name, 0, NULL, NULL,
         -SCIPinfinity(scip), (SCIP_Real) branchdata->constant, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
   }
   else
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "compbnd_%d_child(%d_GE_%d)", uuid, branchdata->Bsize, branchdata->constant);
      SCIP_CALL( SCIPcreateConsLinear(scip, &branchcons, name, 0, NULL, NULL,
         (SCIP_Real) branchdata->constant, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
   }

   SCIP_CALL( SCIPaddConsNode(scip, node, branchcons, NULL) );

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);
   pricingscip = GCGgetPricingprob(origscip, branchdata->blocknr);
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
   SCIP_CALL( GCGcreateInferredPricingVar(pricingscip, &coefvar, pricingvarname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, branchdata->blocknr) );

   // create the y_j variables
   nadditionalvars = branchdata->Bsize;
   SCIP_CALL( SCIPallocBlockMemoryArray(pricingscip, &additionalvars, nadditionalvars) );
   for (i = 0; i < branchdata->Bsize; ++i)
   {
      (void) SCIPsnprintf(pricingvarname, SCIP_MAXSTRLEN, "y(%s,%s_%s_%d)", name, SCIPvarGetName(branchdata->B[i].component), branchdata->B[i].sense == GCG_COMPBND_SENSE_GE ? "GE" : "LE", branchdata->B[i].bound);
      SCIP_CALL( GCGcreateInferredPricingVar(pricingscip, &additionalvars[i], pricingvarname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, branchdata->blocknr) );
   }

   // create the pricing constraints
   if( branchdata->branchtype == GCG_BRANCH_DOWN )
   {
      nadditionalcons = branchdata->Bsize + 1;
      SCIP_CALL( SCIPallocBlockMemoryArray(pricingscip, &additionalcons, nadditionalcons) );

      /* g_x >= 1 + sum_{j=1}^{n} y_j - n
         * 1 - n <= g_x - sum_{j=1}^{n} y_j */
      char consname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "c(g(%s))", name);
      SCIP_CALL( SCIPcreateConsLinear(pricingscip, &additionalcons[0], consname, 0, NULL, NULL, (SCIP_Real) (1 - branchdata->Bsize), SCIPinfinity(pricingscip),
            TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[0], coefvar, 1.0) );
      for( i = 0; i < branchdata->Bsize; ++i)
      {
         SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[0], additionalvars[i], -1.0) );
      }

      for( i = 0; i < branchdata->Bsize; ++i)
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "c(y(%s,%s_%s_%d))", name, SCIPvarGetName(branchdata->B[i].component), branchdata->B[i].sense == GCG_COMPBND_SENSE_GE ? "GE" : "LE", branchdata->B[i].bound);

         if( branchdata->B[i].sense == GCG_COMPBND_SENSE_LE )
         {
            /* y_j >= ((bound + 1) - x_j) / ((bound + 1) - l_j)
               * ((bound + 1) - l_j) * y_j >= (bound + 1) - x_j
               * (bound + 1) <= x_j + ((bound + 1) - l_j) * y_j */
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->B[i].component);
            SCIP_Real bound = (SCIP_Real) (branchdata->B[i].bound + 1);
            SCIP_Real lowerbound = SCIPvarGetLbGlobal(branchdata->B[i].component);
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
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->B[i].component);
            SCIP_Real bound = (SCIP_Real) (branchdata->B[i].bound - 1);
            SCIP_Real upperbound = SCIPvarGetUbGlobal(branchdata->B[i].component);
            assert(SCIPisPositive(pricingscip, upperbound - bound));
            SCIP_CALL( SCIPcreateConsVarbound(pricingscip, &additionalcons[i+1], consname, pricing_var, additionalvars[i], bound - upperbound, -SCIPinfinity(pricingscip), bound,
               TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         }
      }
   }
   else
   {
      nadditionalcons = branchdata->Bsize * 2;
      SCIP_CALL( SCIPallocBlockMemoryArray(pricingscip, &additionalcons, nadditionalcons) );

      /* g_x <= y_j
         * 0 <= y_j - g_y*/
      char consname[SCIP_MAXSTRLEN];
      for( i = 0; i < branchdata->Bsize; ++i )
      {
         SCIP_CALL( SCIPcreateConsVarbound(pricingscip, &additionalcons[i], consname, additionalvars[i], coefvar, -1.0, 0.0, SCIPinfinity(pricingscip),
            TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      }

      for( i = 0; i < branchdata->Bsize; ++i )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "c1(y(%s,%s_%s_%d))", name, SCIPvarGetName(branchdata->B[i].component), branchdata->B[i].sense == GCG_COMPBND_SENSE_GE ? "GE" : "LE", branchdata->B[i].bound);

         if( branchdata->B[i].sense == GCG_COMPBND_SENSE_LE )
         {
            /* y_j <= (u_j - x_j) / (u_j - bound)
               * (u_j - bound) * y_j <= u_j - x_j
               * x_j + (u_j - bound) * y_j <= u_j */
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->B[i].component);
            SCIP_Real bound = (SCIP_Real) branchdata->B[i].bound;
            SCIP_Real upperbound = SCIPvarGetUbGlobal(branchdata->B[i].component);
            assert(SCIPisPositive(pricingscip, upperbound - bound));
            SCIP_CALL( SCIPcreateConsVarbound(pricingscip, &additionalcons[i+branchdata->Bsize], consname, pricing_var, additionalvars[i], upperbound - bound, -SCIPinfinity(pricingscip), upperbound,
               TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         }
         else
         {
            /* y_j <= (x_j - l_j) / (bound - l_j)
               * (bound - l_j) * y_j <= x_j - l_j
               * (bound - l_j) * y_j - x_j <= -l_j
               * l_j <= x_j + (l_j - bound) * y_j */
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->B[i].component);
            SCIP_Real bound = (SCIP_Real) branchdata->B[i].bound;
            SCIP_Real lowerbound = SCIPvarGetLbGlobal(branchdata->B[i].component);
            assert(SCIPisPositive(pricingscip, bound - lowerbound));
            SCIP_CALL( SCIPcreateConsVarbound(pricingscip, &additionalcons[i+branchdata->Bsize], consname, pricing_var, additionalvars[i], lowerbound - bound, lowerbound, SCIPinfinity(pricingscip),
               TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         }
      }
   }

   // create the pricing modification
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &pricingmods, 1) );
   assert(pricingmods != NULL);
   SCIP_CALL( GCGpricingmodificationCreate(
      scip,
      &pricingmods[0],
      branchdata->blocknr,
      coefvar,
      additionalvars,
      nadditionalvars,
      additionalcons,
      nadditionalcons
   ) );

   // create the master constraint
   SCIP_CALL( GCGmastercutCreateFromCons(scip, &branchdata->mastercons, branchcons, pricingmods, 1, branchdata, mastercutGetCoeffCompBnd) );

   // add the variables to the constraint
   for( i = 0; i < nmastervars; i++ )
   {
      added = FALSE;
      SCIP_CALL( addVarToMasterbranch(scip, mastervars[i], branchdata, &added) );
   }

   return SCIP_OKAY;
}

/** for given component bound sequence B create the 2 child nodes */
static
SCIP_RETCODE createChildNodesCompBnd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   GCG_COMPBND*          B,                  /**< Component Bound Sequence defining the nodes */
   int                   Bsize,              /**< size of B */
   int                   blocknr,            /**< number of the block */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching call */
   )
{
   SCIP*  masterscip;
   int i;
   int identicalBlocks;
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

   assert(scip != NULL);
   assert(Bsize > 0);
   assert(B != NULL);

   identicalBlocks = GCGgetNIdenticalBlocks(scip, blocknr);
   SCIPdebugMessage("Component bound branching rule Node creation for blocknr %d with %d identical blocks \n", blocknr, identicalBlocks);

   /*  get variable data of the master problem */
   masterscip = GCGgetMasterprob(scip);
   assert(masterscip != NULL);
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(nmastervars >= 0);
   /* determine the constant value of the master constraint, i.e. the rhs for the down branch, and lhs for the up branch */
   constantSum = 0.0;
   for( i = 0; i < nmastervars; ++i )
   {
      if(
         GCGisMasterVarInBlock(
            mastervars[i],
            GCGgetBlockRepresentative(scip, blocknr)
         )
         && isMasterVarInB(scip, mastervars[i], B, Bsize, blocknr)
      )
      {
         constantSum += SCIPgetSolVal(masterscip, NULL, mastervars[i]);
      }
   }
   // sanity check: the sum must be fractional, otherwise something went wrong during separation, i.e. B is incorrect
   assert(!SCIPisFeasIntegral(scip, constantSum));

   /* create two nodes */
   SCIPdebugMessage("Component bound branching rule: creating 2 nodes\n");
   SCIP_CALL( initNodeBranchdata(masterscip, &downBranchData, GCG_BRANCH_DOWN, FLOOR(scip, constantSum), B, Bsize, blocknr) );
   SCIP_CALL( initNodeBranchdata(masterscip, &upBranchData, GCG_BRANCH_UP, CEIL(scip, constantSum), B, Bsize, blocknr) );
   freeComponentBoundSequence(masterscip, &B, &Bsize);
   assert(downBranchData != NULL);
   assert(upBranchData != NULL);
   assert(B == NULL);
   assert(Bsize == 0);

   /* define names for origbranch constraints */
   (void) SCIPsnprintf(downChildname, SCIP_MAXSTRLEN, "node(%d) (last comp=%s %s %d) <= %d", blocknr,
      SCIPvarGetName(downBranchData->B[downBranchData->Bsize-1].component),
      downBranchData->B[downBranchData->Bsize-1].sense == GCG_COMPBND_SENSE_GE ? ">=": "<=",
      downBranchData->B[downBranchData->Bsize-1].bound,
      downBranchData->constant);
   (void) SCIPsnprintf(upChildname, SCIP_MAXSTRLEN, "node(%d) (last comp=%s %s %d) >= %d", blocknr,
      SCIPvarGetName(upBranchData->B[upBranchData->Bsize-1].component),
      upBranchData->B[upBranchData->Bsize-1].sense == GCG_COMPBND_SENSE_GE ? ">=": "<=",
      upBranchData->B[upBranchData->Bsize-1].bound,
      upBranchData->constant);


   /* add the nodes to the tree */
   SCIP_CALL( SCIPcreateChild(masterscip, &downChild, 0.0, SCIPgetLocalTransEstimate(masterscip)) );
   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &downChildcons, downChildname, downChild,
      GCGconsMasterbranchGetActiveCons(masterscip), branchrule, downBranchData, NULL, 0, 0) );
   SCIP_CALL( SCIPaddConsNode(masterscip, downChild, downChildcons, NULL) );
   SCIP_CALL( createBranchingCons(masterscip, downChild, downBranchData) );

   SCIP_CALL( SCIPcreateChild(masterscip, &upChild, 0.0, SCIPgetLocalTransEstimate(masterscip)) );
   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &upChildcons, upChildname, upChild,
      GCGconsMasterbranchGetActiveCons(masterscip), branchrule, upBranchData, NULL, 0, 0) );
   SCIP_CALL( SCIPaddConsNode(masterscip, upChild, upChildcons, NULL) );
   SCIP_CALL( createBranchingCons(masterscip, upChild, upBranchData) );

   /*  release constraints */
   SCIP_CALL( SCIPreleaseCons(masterscip, &upChildcons) );
   SCIP_CALL( SCIPreleaseCons(masterscip, &downChildcons) );


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
      // We do this check to avoid numerical issues
      if( !SCIPisFeasIntegral(masterscip, solval) )
      {
         fractionality += solval - SCIPfloor(masterscip, solval);
      }
   }

   return fractionality;
}

/** method for initializing the set of respected indices */
static
SCIP_RETCODE initIndexSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           indexSet,           /**< set to initialize */
   int*                  indexSetSize,       /**< size of the index set */
   SCIP_VAR**            X,                  /**< mastervariables currently satisfying the component bound sequence in the specified block */
   int                   Xsize               /**< size of X */
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

   for( i = 0; i < Xsize; ++i )
   {
      origvars = GCGmasterVarGetOrigvars(X[i]);
      norigvars = GCGmasterVarGetNOrigvars(X[i]);

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
               /*  if variable already in union */
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
   GCG_COMPBND*          B,                  /**< component bound sequence */
   int                   Bsize               /**< size of B */
   )
{
   int i;
   SCIP_Real lb;

   assert(scip != NULL);
   assert(origvar != NULL);

   lb = SCIPvarGetLbGlobal(origvar);

   for( i = 0; i < Bsize; ++i )
   {
      if(
         B[i].sense == GCG_COMPBND_SENSE_GE
         && SCIPvarCompare(origvar, B[i].component) == 0
         && SCIPisLT(scip, lb, B[i].bound)
      )
      {
         lb = (SCIP_Real) B[i].bound;
      }
   }

   return lb;
}

/** get the upper bound of a original variable from the component bound sequence, or the ub */
static
SCIP_Real getUpperBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar,            /**< original variable */
   GCG_COMPBND*          B,                  /**< component bound sequence */
   int                   Bsize               /**< size of B */
   )
{
   int i;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(origvar != NULL);

   ub = SCIPvarGetUbGlobal(origvar);

   for( i = 0; i < Bsize; ++i )
   {
      if(
         B[i].sense == GCG_COMPBND_SENSE_LE
         && SCIPvarCompare(origvar, B[i].component) == 0
         && SCIPisLT(scip, B[i].bound, ub)
      )
      {
         ub = (SCIP_Real) B[i].bound;
      }
   }

   return ub;
}
#endif

/** separation algorithm determining a component bound sequence to branch on */
static
SCIP_RETCODE _separation(
   SCIP*                 masterscip,              /**< SCIP data structure */
   SCIP_VAR**            X,                       /**< mastervariables currently satisfying the component bound sequence in the specified block */
   int                   Xsize,                   /**< size of X */
   GCG_COMPBND*          B,                       /**< current Component Bound Sequence defining the nodes */
   int                   Bsize,                   /**< size of B */
   int                   blocknr,                 /**< number of the block */
   CHOOSE_COMPBND((*chooseCompbnd)),              /**< method for choosing the variable to branch on */
   GCG_COMPBND***        Blist,                   /**< list of component bound sequences */
   int**                 Bsizes,                  /**< sizes of the component bound sequences */
   int*                  Blistsize,               /**< size of Blist */
   SCIP_RESULT*          result                   /**< pointer to store the result of the branching call */
   )
{
   SCIP_Real fractionality;
   int i;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));

   // Sanity check: All variables in X must satisfy the component bound sequence
   if( Bsize > 0 ) {
      for( i = 0; i < Xsize; ++i )
      {
         assert(isMasterVarInB(masterscip, X[i], B, Bsize, blocknr));
      }
   }

   fractionality = calcFractionality(masterscip, X, Xsize);
   assert(fractionality >= 0.0);

   if( SCIPisEQ(masterscip, fractionality, 0.0) ) // should never happen - the branching scheme is sound and complete
   {
      // all variables are integral, nothing to do
      SCIPdebugMessage("All variables are integral, nothing to do\n");

      *result = SCIP_DIDNOTFIND;

      return SCIP_ERROR;
   }

   if( Bsize > 0 && !SCIPisFeasIntegral(masterscip, fractionality) )
   {
      /* now we branch on <= floor(fractionality) and >= ceil(fractionality)
       * this is handled by the createChildNodesCompBnd method
       */
      SCIPdebugMessage("fractionality is fractional, we need to branch\n");

      *result = SCIP_BRANCHED;

      // add current B to Blist
      if( *Blistsize == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, Blist, 1) );
         SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, Bsizes, 1) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, Blist, *Blistsize, *Blistsize + 1) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, Bsizes, *Blistsize, *Blistsize + 1) );
      }
      (*Blist)[*Blistsize] = B;
      (*Bsizes)[*Blistsize] = Bsize;
      *Blistsize += 1;

      return SCIP_OKAY;
   }

   // the fractionality is integral, we need to impose an additional bound
   SCIP_VAR** indexSet = NULL;
   int indexSetSize = 0;
   SCIP_CALL( initIndexSet(masterscip, &indexSet, &indexSetSize, X, Xsize) );

   SCIP_VAR* selected_origvar = NULL;
   SCIP_Real value = 0.0;
   SCIP_CALL( chooseCompbnd(masterscip, X, Xsize, indexSet, indexSetSize, B, Bsize, blocknr, &selected_origvar, &value) );
   assert(selected_origvar != NULL);
   assert(SCIPvarGetLbGlobal(selected_origvar) < value && value < SCIPvarGetUbGlobal(selected_origvar));

   SCIPfreeBlockMemoryArray(masterscip, &indexSet, indexSetSize);
   indexSetSize = 0;

   // create two component bound options, xj <= floor(value) and xj >= floor(value)+1
   GCG_COMPBND* B1;
   GCG_COMPBND* B2;
   int new_Bsize = Bsize + 1;
   // allocate memory for B1 and B2
   SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &B1, new_Bsize) );

   SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &B2, new_Bsize) );
   // copy the old component bounds to B1 and B2
   for(i = 0; i < Bsize; ++i)
   {
      B1[i] = B[i];
      B2[i] = B[i];
   }
   // add the new component bounds to B1 and B2
   B1[Bsize] = (GCG_COMPBND){selected_origvar, GCG_COMPBND_SENSE_LE, FLOOR(masterscip, value)};
   B2[Bsize] = (GCG_COMPBND){selected_origvar, GCG_COMPBND_SENSE_GE, FLOOR(masterscip, value) + 1};

   freeComponentBoundSequence(masterscip, &B, &Bsize);
   assert(B == NULL);

   SCIPdebugMessage("B1 and B2 after adding the new component bounds\n");
   for (i = 0; i < new_Bsize; ++i)
   {
      SCIPdebugMessage("B1[%d]: %s %s %d\n", i, SCIPvarGetName(B1[i].component),
                        B1[i].sense == GCG_COMPBND_SENSE_LE ? "<=" : ">=",
                        B1[i].bound);
      SCIPdebugMessage("B2[%d]: %s %s %d\n", i, SCIPvarGetName(B2[i].component),
                        B2[i].sense == GCG_COMPBND_SENSE_LE ? "<=" : ">=",
                        B2[i].bound);
   }

   // assign the master variables to X1 and X2, depending on whether they satisfy the new component bound sequences
   SCIP_VAR** X1 = NULL;
   SCIP_VAR** X2 = NULL;
   int X1size = 0;
   int X2size = 0;
   int x;
   for( x=0; x<Xsize; ++x )
   {
#ifndef NDEBUG
      SCIP_Bool inB1 = FALSE;
      SCIP_Bool inB2 = FALSE;
#endif
      if( isMasterVarInB(masterscip, X[x], B1, new_Bsize, blocknr) )
      {
         // increase the size of X1 by 1, and add the current variable to X1
         if( X1size == 0 )
            SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &X1, X1size+1) );
         else
            SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, &X1, X1size, X1size+1) );

         X1[X1size] = X[x];
         X1size += 1;
#ifndef NDEBUG
         inB1 = TRUE;
#endif
      }
      if( isMasterVarInB(masterscip, X[x], B2, new_Bsize, blocknr) )
      {
         // increase the size of X2 by 1, and add the current variable to X2
         if( X2size == 0 )
            SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &X2, X2size+1) );
         else
            SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, &X2, X2size, X2size+1) );

         X2[X2size] = X[x];
         X2size += 1;
#ifndef NDEBUG
         inB2 = TRUE;
#endif
      }
#ifndef NDEBUG
      assert(inB1 + inB2 <= 1);
#endif
   }
   SCIPdebugMessage("X1size: %d, X2size: %d\n", X1size, X2size);
   assert(X1size > 0);
   assert(X2size > 0);

#ifndef NDEBUG
   // determine the fractionality of B1 and B2
   SCIP_Real fractionality1 = calcFractionality(masterscip, X1, X1size);
   SCIP_Real fractionality2 = calcFractionality(masterscip, X2, X2size);

   assert(SCIPisPositive(masterscip, fractionality1));
   assert(SCIPisPositive(masterscip, fractionality2));
   assert(SCIPisLE(masterscip, fractionality / 2, fractionality1) || SCIPisLE(masterscip, fractionality / 2, fractionality2));
   assert(
      (SCIPisFeasIntegral(masterscip, fractionality1) && SCIPisFeasIntegral(masterscip, fractionality2)) ||
      (!SCIPisFeasIntegral(masterscip, fractionality1) && !SCIPisFeasIntegral(masterscip, fractionality2))
   );
#endif

   // recursive call
   SCIP_CALL( _separation(masterscip, X1, X1size, B1, new_Bsize, blocknr, chooseCompbnd, Blist, Bsizes, Blistsize, result) );

   SCIPfreeBlockMemoryArray(masterscip, &X1, X1size);
   X1size = 0;
   assert(X1 == NULL);

   SCIP_CALL( _separation(masterscip, X2, X2size, B2, new_Bsize, blocknr, chooseCompbnd, Blist, Bsizes, Blistsize, result) );

   SCIPfreeBlockMemoryArray(masterscip, &X2, X2size);
   X2size = 0;
   assert(X2 == NULL);

   return SCIP_OKAY;
}

/** create initial set X */
static
SCIP_RETCODE createInitialSetX(
   SCIP*                 masterscip,              /**< SCIP data structure */
   SCIP_VAR***           X,                       /**< mastervariables */
   int*                  Xsize,                   /**< size of X */
   int                   blocknr,                 /**< number of the block */
   GCG_COMPBND*          B,                       /**< Component Bound Sequence defining the nodes */
   int                   Bsize                    /**< size of B */
   )
{
   SCIP_VAR** mastervars;
   int nmastervars;
   int i;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(X != NULL);
   assert(Xsize != NULL);

   mastervars = NULL;
   nmastervars = 0;

   assert(masterscip != NULL);

   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   if(*Xsize > 0)
   {
      assert(*X != NULL);
      SCIPfreeBlockMemoryArray(masterscip, X, *Xsize);
      *Xsize = 0;
   }

   assert(*Xsize == 0 && *X == NULL);

   for( i = 0; i < nmastervars; ++i )
   {
      if( !GCGisMasterVarInBlock(mastervars[i], blocknr) )
         continue;

      if( Bsize > 0 && !isMasterVarInB(masterscip, mastervars[i], B, Bsize, blocknr) )
            continue;

      if( *Xsize == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, X, 1) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, X, *Xsize, *Xsize+1) );
      }
      (*X)[*Xsize] = mastervars[i];
      (*Xsize) += 1;
   }

   return SCIP_OKAY;
}

/** choose to separate on original variable where max-min is maximal */
static
CHOOSE_COMPBND(chooseMaxMin)
{
   // Task: Find fractional mastervarMin, mastervarMax, indexed by x1, x2, s.t. for some j, xj1 < xj2
   SCIP_VAR* current_mastervar;
   SCIP_Real current_solution_value; // of the currMastervar
   SCIP_VAR* current_origvar;
   SCIP_Real current_min;
   SCIP_Real current_max;
   SCIP_Real current_diff;
   SCIP_Real largest_diff_min;
   SCIP_Real largest_diff_max;
   SCIP_Real largest_diff = -SCIPinfinity(masterscip);
   int i;
   int j;
   SCIP_Real generatorentry;
   SCIP_Bool found = FALSE;

#ifndef NDEBUG
   SCIP_Real given_min;
   SCIP_Real given_max;
#endif

   for( j=0; j<indexSetSize; ++j )
   {
      current_origvar = indexSet[j];
#ifndef NDEBUG
      given_min = getLowerBound(masterscip, current_origvar, B, Bsize);
      given_max = getUpperBound(masterscip, current_origvar, B, Bsize);
      assert(SCIPisFeasLE(masterscip, given_min, given_max));
#endif
      current_min = SCIPinfinity(masterscip);
      current_max = -SCIPinfinity(masterscip);
      for( i=0; i<Xsize; ++i )
      {
         current_mastervar = X[i];

         // Current solution value of the master variable must be fractional and > 0
         current_solution_value = SCIPgetSolVal(masterscip, NULL, current_mastervar);
         if ( !SCIPisFeasPositive(masterscip, current_solution_value) || SCIPisFeasIntegral(masterscip, current_solution_value) ) {
            // skip
            continue;
         }

         generatorentry = getGeneratorEntry(current_mastervar, current_origvar);

         // first, check if we can update the min and max values
         if( SCIPisInfinity(masterscip, current_min) )
         {
            assert(SCIPisInfinity(masterscip, -current_max));
            current_min = generatorentry;
            current_max = generatorentry;
         }
         else if( SCIPisLT(masterscip, generatorentry, current_min) )
            current_min = generatorentry;
         else if( SCIPisLT(masterscip, current_max, generatorentry) )
            current_max = generatorentry;

         // now could be that min<max
         if( SCIPisLT(masterscip, current_min, current_max) )
         {
            found = TRUE;
            current_diff = current_max - current_min;
            if( SCIPisGT(masterscip, current_diff, largest_diff) )
            {
               largest_diff_min = current_min;
               largest_diff_max = current_max;
               largest_diff = current_diff;
               *selected_origvar = current_origvar;
            }
         }
      }
   }

   if( !found ) // should always be able to find - the branching scheme is sound and complete
   {
      return SCIP_ERROR;
   }

   assert(largest_diff_min < largest_diff_max);
   assert(largest_diff > 0);
   assert(SCIPisFeasLE(masterscip, SCIPvarGetLbGlobal(*selected_origvar), largest_diff_min));
   assert(SCIPisFeasLE(masterscip, largest_diff_max, SCIPvarGetUbGlobal(*selected_origvar)));
   SCIPdebugMessage("Found %d for origvar %s, j=%d, min=%f, max=%f\n", found, SCIPvarGetName(*selected_origvar), j, largest_diff_min, largest_diff_max);

   // j and origivar are still set to the correct values after execution of the loop
   *value = (largest_diff_min+largest_diff_max)/2;

   return SCIP_OKAY;
}

/** choose to separate on original variable which has the largest number of distinct numbers, and find the median value */
static
CHOOSE_COMPBND(chooseMostDistinctValuesMedian)
{
   int selected_numDistinctValues = 0;

   SCIP_VAR* current_origvar;
   int current_numDistinctValues;
   SCIP_Real* distinctValues = NULL;
   int distinctValuesSize = 0;

   SCIP_Bool inserted;
   SCIP_Bool duplicate;

   int i;
   int j;
   int k;
   int l;

   assert(*selected_origvar == NULL);

   SCIP_Real generatorentry;

   for( j=0; j<indexSetSize; ++j )
   {
      current_origvar = indexSet[j];
      current_numDistinctValues = 0;

      for( i=0; i<Xsize; ++i )
      {
         SCIP_VAR* current_mastervar = X[i];

         // Current solution value of the master variable must be fractional and > 0
         SCIP_Real current_solution_value = SCIPgetSolVal(masterscip, NULL, current_mastervar);
         if ( !SCIPisFeasPositive(masterscip, current_solution_value) || SCIPisFeasIntegral(masterscip, current_solution_value) ) {
            // skip
            continue;
         }

         generatorentry = getGeneratorEntry(current_mastervar, current_origvar);

         // insert sorted, if the value is not already in the array
         if( distinctValuesSize == 0 )
         {
            assert(distinctValues == NULL);
            assert(current_numDistinctValues == 0);
            distinctValuesSize = SCIPcalcMemGrowSize(masterscip, 1);
            SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &distinctValues, distinctValuesSize) );
            distinctValues[0] = generatorentry;
            current_numDistinctValues = 1;
         }
         else
         {
            // current_distinctValues is sorted ascendingly
            inserted = FALSE;
            duplicate = FALSE;
            for( k=0; k<current_numDistinctValues; ++k )
            {
               if( SCIPisEQ(masterscip, generatorentry, distinctValues[k]) )
               {
                  duplicate = TRUE;
                  break;
               }

               if( SCIPisLT(masterscip, generatorentry, distinctValues[k]) )
               {
                  // insert at position k
                  if( current_numDistinctValues == distinctValuesSize )
                  {
                     int new_current_distinctValuesSize = SCIPcalcMemGrowSize(masterscip, current_numDistinctValues + 1);
                     SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, &distinctValues, distinctValuesSize, new_current_distinctValuesSize) );
                     distinctValuesSize = new_current_distinctValuesSize;
                  }

                  for( l = current_numDistinctValues; l > k; --l )
                  {
                     distinctValues[l] = distinctValues[l-1];
                  }

                  distinctValues[k] = generatorentry;
                  current_numDistinctValues += 1;

                  inserted = TRUE;
                  break;
               }
            }

            if( !inserted && !duplicate )
            {
               // insert at the end
               if( current_numDistinctValues == distinctValuesSize )
               {
                  int new_current_distinctValuesSize = SCIPcalcMemGrowSize(masterscip, current_numDistinctValues + 1);
                  SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, &distinctValues, distinctValuesSize, new_current_distinctValuesSize) );
                  distinctValuesSize = new_current_distinctValuesSize;
               }

               distinctValues[current_numDistinctValues] = generatorentry;
               current_numDistinctValues += 1;
            }
         }
      }

      if( current_numDistinctValues > 1 && current_numDistinctValues > selected_numDistinctValues )
      {
         assert(current_numDistinctValues > 1);
         *selected_origvar = current_origvar;
         selected_numDistinctValues = current_numDistinctValues;

         // calculate median
         if( current_numDistinctValues % 2 == 0 )
         {
            *value = 0.5 * (distinctValues[current_numDistinctValues / 2] + distinctValues[current_numDistinctValues / 2 - 1]);
         }
         else
         {
            *value = distinctValues[current_numDistinctValues / 2];
         }
      }
   }

   assert(*selected_origvar != NULL);

   SCIPfreeBlockMemoryArrayNull(masterscip, &distinctValues, distinctValuesSize);

   return SCIP_OKAY;
}

/** calculate sum of component bound sequence */
static
SCIP_Real calcSum(
   SCIP*                 masterscip,              /**< SCIP data structure */
   GCG_COMPBND*          B,                       /**< Component Bound Sequence defining the nodes */
   int                   Bsize,                   /**< size of B */
   int                   blocknr                  /**< number of the block */
   )
{
   int i;
   SCIP_Real sum = 0.0;
   SCIP_VAR** X = NULL;
   int Xsize = 0;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));

   assert(X == NULL && Xsize == 0);
   createInitialSetX(masterscip, &X, &Xsize, blocknr, B, Bsize);
   assert((Xsize > 0) == (X != NULL));

   sum = 0;
   for( i = 0; i < Xsize; ++i )
   {
      sum += SCIPgetSolVal(masterscip, NULL, X[i]);
   }

   SCIPfreeBlockMemoryArrayNull(masterscip, &X, Xsize);

   return sum;
}

/** choose the component bound of which the sum is closest to K/2, where K is the number of identical subproblems in the current block */
static
CHOOSE_COMPBNDSEQ(chooseClosestToKHalf)
{
   SCIP* origscip;
   int i;
   SCIP_Real sum;
   SCIP_Real K_half;

   int best_index = -1;
   SCIP_Real best = SCIPinfinity(masterscip);

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(Blist != NULL);
   assert(Bsizes != NULL);
   assert(Blistsize > 0);

   origscip = GCGmasterGetOrigprob(masterscip);
   assert(origscip != NULL);

   K_half = GCGgetNIdenticalBlocks(origscip, blocknr) * 0.5;

   for( i = 0; i < Blistsize; ++i )
   {
      sum = ABS(calcSum(masterscip, Blist[i], Bsizes[i], blocknr) - K_half);

      if( sum < best )
      {
         best = sum;
         best_index = i;
      }
   }

   assert(best_index >= 0);

   *selected_B = Blist[best_index];
   *selected_Bsize = Bsizes[best_index];

   return SCIP_OKAY;
}

/** choose the component bound of which the sum is most fractional */
static
CHOOSE_COMPBNDSEQ(chooseMostFractional)
{
   int i;
   SCIP_Real sum;
   SCIP_Real fractionality;

   int best_index = -1;
   SCIP_Real best = 0.0;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(Blist != NULL);
   assert(Bsizes != NULL);
   assert(Blistsize > 0);

   for( i = 0; i < Blistsize; ++i )
   {
      sum = calcSum(masterscip, Blist[i], Bsizes[i], blocknr);
      fractionality = SCIPfrac(masterscip, sum);

      fractionality = MIN(fractionality, 1.0 - fractionality);

      if( best < fractionality )
      {
         best = fractionality;
         best_index = i;
      }
   }

   assert(best_index >= 0);

   *selected_B = Blist[best_index];
   *selected_Bsize = Bsizes[best_index];

   return SCIP_OKAY;
}

/** choose the component bound sequence with the smallest size */
static
CHOOSE_COMPBNDSEQ(chooseSmallestSize)
{
   int i;
   int smallest_size = INT_MAX;

   GCG_COMPBND** newBlist = NULL;
   int* newBsizes = NULL;
   int newBlistsize = 0;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(Blist != NULL);
   assert(Bsizes != NULL);
   assert(Blistsize > 0);

   // step 1: find smallest size
   for( i = 0; i < Blistsize; ++i )
   {
      if( Bsizes[i] < smallest_size )
      {
         smallest_size = Bsizes[i];
      }
   }
   assert(smallest_size > 0);

   // step 2: create new Blist with smallest size
   for( i = 0; i < Blistsize; ++i )
   {
      if( Bsizes[i] == smallest_size )
      {
         if( newBlistsize == 0 )
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &newBlist, 1) );
            SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &newBsizes, 1) );
         }
         else
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, &newBlist, newBlistsize, newBlistsize + 1) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, &newBsizes, newBlistsize, newBlistsize + 1) );
         }

         newBlist[newBlistsize] = Blist[i];
         newBsizes[newBlistsize] = Bsizes[i];
         newBlistsize += 1;
      }
   }

   assert(newBlistsize > 0);

   // step 3: refine the selection
   //SCIP_CALL( chooseClosestToKHalf(masterscip, newBlist, newBsizes, newBlistsize, blocknr, selected_B, selected_Bsize) );
   SCIP_CALL( chooseMostFractional(masterscip, newBlist, newBsizes, newBlistsize, blocknr, selected_B, selected_Bsize) );

   // step 4: free memory
   SCIPfreeBlockMemoryArray(masterscip, &newBlist, newBlistsize);
   SCIPfreeBlockMemoryArray(masterscip, &newBsizes, newBlistsize);

   return SCIP_OKAY;
}

/** separation algorithm determining a component bound sequence to branch on */
static
SCIP_RETCODE separation(
   SCIP*                 masterscip,              /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,              /**< branching rule */
   GCG_COMPBND**         B,                       /**< Component Bound Sequence defining the nodes */
   int*                  Bsize,                   /**< size of B */
   int                   blocknr,                 /**< number of the block */
   SCIP_RESULT*          result                   /**< pointer to store the result of the branching call */
   )
{
   // a list of component bound sequences (also lists)
   GCG_COMPBND** Blist;
   int* Bsizes;
   int Blistsize;

   SCIP_BRANCHRULEDATA* branchruledata;
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIP_VAR** X;
   int Xsize;

   int i;

   Blist = NULL;
   Bsizes = NULL;
   Blistsize = 0;

   X = NULL;
   Xsize = 0;

   createInitialSetX(masterscip, &X, &Xsize, blocknr, NULL, 0);
#ifndef NDEBUG
   int XsizeCopy = Xsize;
   int BlistsizeCopy = Blistsize;
#endif

   // call the separation algorithm with all chooseVars functions
   assert(branchruledata->useMaxRangeMidrangeHeuristic || branchruledata->useMostDistinctMedianHeuristic);

   if( branchruledata->useMaxRangeMidrangeHeuristic )
   {
      SCIP_CALL( _separation(masterscip, X, Xsize, NULL, 0, blocknr, chooseMaxMin, &Blist, &Bsizes, &Blistsize, result) );
#ifndef NDEBUG
      assert(Xsize == XsizeCopy);
      assert(Blistsize > BlistsizeCopy);
      BlistsizeCopy = Blistsize;
#endif
   }

   if( branchruledata->useMostDistinctMedianHeuristic )
   {
      SCIP_CALL( _separation(masterscip, X, Xsize, NULL, 0, blocknr, chooseMostDistinctValuesMedian, &Blist, &Bsizes, &Blistsize, result) );
#ifndef NDEBUG
      assert(Xsize == XsizeCopy);
      assert(Blistsize > BlistsizeCopy);
      BlistsizeCopy = Blistsize;
#endif
   }

   SCIPfreeBlockMemoryArray(masterscip, &X, Xsize);
   Xsize = 0;
   assert(X == NULL);

   assert(Blistsize > 0);
   assert(Blist != NULL);
   assert(Bsizes != NULL);

   SCIP_CALL( chooseSmallestSize(masterscip, Blist, Bsizes, Blistsize, blocknr, B, Bsize) );

   assert(*Bsize > 0);
   assert(*B != NULL);

   // free others, as well as Blist
   for( i = 0; i < Blistsize; ++i )
   {
      if( Blist[i] != *B )
      {
         freeComponentBoundSequence(masterscip, &Blist[i], &Bsizes[i]);
      }
   }

   SCIPfreeBlockMemoryArray(masterscip, &Blist, Blistsize);
   SCIPfreeBlockMemoryArray(masterscip, &Bsizes, Blistsize);

   return SCIP_OKAY;
}

/** prepares information for using the generic branching scheme */
SCIP_RETCODE GCGbranchCompBndInitbranch(
   SCIP*                 masterscip,              /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,              /**< branching rule */
   SCIP_RESULT*          result                   /**< pointer to store the result of the branching call */
   )
{
   SCIP* origscip;
#ifdef STATS
   SCIP_BRANCHRULEDATA* branchruledata;
#endif
   SCIP_VAR** branchcands;
   int nbranchcands;
   GCG_COMPBND* B;
   int Bsize;
   int blocknr;
   SCIP_Bool foundBlock;
   SCIP_VAR** X;
   int Xsize;
   SCIP_Real fractionality;

   blocknr = -2;
   B = NULL;
   Bsize = 0;

   assert(masterscip != NULL);

   SCIPdebugMessage("get information for component bound branching\n");

   origscip = GCGmasterGetOrigprob(masterscip);

   assert(origscip != NULL);
   SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL, NULL) );

#ifdef STATS
   branchruledata = SCIPbranchruleGetData(branchrule);
#endif

   /* in case original problem contains continuous variables, there are no branching cands */
   assert(nbranchcands > 0);

   foundBlock = FALSE;
   Xsize = 0;
   X = NULL;

   /* 1. Determine in what block we are branching. We select the first available block,
    *     i.e. the first block that contains a branching candidate, starting from the master block.
    */
   for( blocknr = 0; blocknr < GCGgetNPricingprobs(origscip); blocknr++ )
   {
      if( !GCGisPricingprobRelevant(origscip, blocknr) )
         continue;

      SCIPdebugMessage("\nTrying to branch in block %d:\n", blocknr);

      /* 2. Check whether the fractionality of all master variables in this block is not 0. */
      createInitialSetX(masterscip, &X, &Xsize, blocknr, NULL, 0);
      assert((Xsize > 0) == (X != NULL));
      fractionality = calcFractionality(masterscip, X, Xsize);
      assert(SCIPisFeasIntegral(masterscip, fractionality));
      SCIPfreeBlockMemoryArray(masterscip, &X, Xsize);
      Xsize = 0;
      assert(X == NULL);
      if( SCIPisZero(masterscip, fractionality) )
      {
         SCIPdebugMessage("No fractional integer variables in block %d\n", blocknr);
         continue;
      }

      /* 3. Call to separation algorithm to find a suitable B to branch on in the current block. */
      SCIP_CALL( separation(masterscip, branchrule, &B, &Bsize, blocknr, result) );

      if( *result == SCIP_BRANCHED )
      {
         assert(Bsize > 0);
         assert(B != NULL);
         SCIPdebugMessage("Branching in block %d\n", blocknr);
         foundBlock = TRUE;
         break;
      }
      else
      {
         if( Bsize > 0 )
         {
            freeComponentBoundSequence(masterscip, &B, &Bsize);
            Bsize = 0;
         }
      }
   }

   if( !foundBlock )
   {
      assert(Bsize == 0);
      assert(B == NULL);
      SCIPdebugMessage("No block found to branch on\n");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Component bound sequence to branch on:\n");
   SCIPdebugMessage("Bsize: %d\n", Bsize);
   for (int i = 0; i < Bsize; ++i)
   {
      SCIPdebugMessage("B[%d]: %s %s %d\n", i, SCIPvarGetName((B)[i].component),
                        (B)[i].sense == GCG_COMPBND_SENSE_LE ? "<=" : ">=",
                        (B)[i].bound);
   }

   /* 4. Remove component bounds that are strengthened by others */
   SCIP_CALL( simplifyComponentBounds(masterscip, &B, &Bsize) );
   assert(Bsize > 0);
   assert(B != NULL);

#ifdef STATS
   int depth = SCIPgetDepth(masterscip);
   SCIP_Real sum = calcSum(masterscip, B, Bsize, blocknr);
   int K = GCGgetNIdenticalBlocks(origscip, blocknr);
   fprintf(branchruledata->statsfile, "%d,%d,%f,%d\n", depth, Bsize, sum, K);
#endif

   /* 5. Create the child nodes. */
   SCIP_CALL( createChildNodesCompBnd(origscip, branchrule, B, Bsize, blocknr, result) );

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

#ifdef STATS
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
   SCIP* origscip;
   SCIP_Bool discretization;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   SCIPdebugMessage("Execlp method of component bound branching\n");

   *result = SCIP_DIDNOTRUN;

   /* the branching scheme only works for the discretization approach */
   SCIP_CALL( SCIPgetBoolParam(origscip, "relaxing/gcg/discretization", &discretization) );
   if( !discretization )
   {
      SCIPdebugMessage("Component bound branching only for discretization approach\n");
      return SCIP_OKAY;
   }

   if( GCGisMasterSetCovering(origscip) || GCGisMasterSetPartitioning(origscip) )
   {
      SCIPdebugMessage("Component bound branching executed on a set covering or set partitioning problem\n");
   }

   if( GCGrelaxIsOrigSolFeasible(origscip) )
   {
      SCIPdebugMessage("node cut off, since origsol was feasible, solval = %f\n",
         SCIPgetSolOrigObj(origscip, GCGrelaxGetCurrentOrigSol(origscip)));

      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   *result = SCIP_BRANCHED;

   SCIP_CALL( GCGbranchCompBndInitbranch(scip, branchrule, result) );

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
   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   SCIPdebugMessage("branchActiveMasterCompBnd: Block %d, Ssize %d\n", branchdata->blocknr, branchdata->Bsize);

   return SCIP_OKAY;
}


/** deactivation method for branchrule, called when a node in the master problem is deactivated,
 *  should undo changes to the current node's problem due to the branchdata
 */
static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterCompBnd)
{
   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   SCIPdebugMessage("branchDeactiveMasterCompBnd: Block %d, Ssize %d\n", branchdata->blocknr, branchdata->Bsize);

   return SCIP_OKAY;
}

/** propagation method for branchrule, called when a node in the master problem is propagated,
 *  should perform propagation at the current node due to the branchdata
 */
static
GCG_DECL_BRANCHPROPMASTER(branchPropMasterCompBnd)
{
   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);
   assert(branchdata->B != NULL);

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

   assert(scip != NULL);
   assert(GCGisOriginal(scip));
   assert(branchdata != NULL);

   masterscip = GCGgetMasterprob(scip);
   assert(masterscip != NULL);

   SCIPdebugMessage("branchDataDeleteCompBnd: Block %d, Ssize %d\n", (*branchdata)->blocknr, (*branchdata)->Bsize);

   assert(*branchdata != NULL);

   /* release constraint that enforces the branching decision */
   assert((*branchdata)->mastercons != NULL);
   SCIP_CALL( GCGmastercutFree(masterscip, &(*branchdata)->mastercons) );
   assert((*branchdata)->mastercons == NULL);

   assert((*branchdata)->B != NULL && (*branchdata)->Bsize > 0);
   SCIPfreeBlockMemoryArray(masterscip, &((*branchdata)->B), (*branchdata)->Bsize);
   assert((*branchdata)->B == NULL);
   (*branchdata)->Bsize = 0;

   SCIPfreeBlockMemory(masterscip, branchdata);
   assert(*branchdata == NULL);

   return SCIP_OKAY;
}

/** callback new column method */
static
GCG_DECL_BRANCHNEWCOL(branchNewColCompBnd)
{
   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(mastervar != NULL);
   assert(GCGvarIsMaster(mastervar));
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   SCIP_Bool added = FALSE;
   SCIP_CALL( addVarToMasterbranch(scip, mastervar, branchdata, &added) );

   return SCIP_OKAY;
}

static
GCG_DECL_BRANCHGETMASTERCUT(branchGetMastercutCompBnd)
{
   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   *mastercutdata = branchdata->mastercons;

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitCompBnd)
{
   SCIP* origscip;
#ifdef STATS
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_Bool stabilization_tree;
#endif

   origscip = GCGmasterGetOrigprob(scip);
   assert(branchrule != NULL);
   assert(origscip != NULL);

#ifdef STATS
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   SCIP_CALL( SCIPgetBoolParam(origscip, "pricing/masterpricer/stabilizationtree", &stabilization_tree) );
#endif

   SCIPdebugMessage("Init method of component bound branching\n");

   SCIP_CALL( GCGrelaxIncludeBranchrule(origscip, branchrule, branchActiveMasterCompBnd,
         branchDeactiveMasterCompBnd, branchPropMasterCompBnd, branchMasterSolvedCompBnd, branchDataDeleteCompBnd,
         branchNewColCompBnd, branchGetMastercutCompBnd) );

#ifdef STATS
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
SCIP_RETCODE SCIPincludeBranchruleCompBnd(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP* origscip;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* get original problem */
   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   /* create compbnd branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );

   SCIPdebugMessage("Include method of component bound branching\n");

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
         BRANCHRULE_MAXBOUNDDIST,
         branchCopyCompBnd, branchFreeCompBnd, branchInitCompBnd, branchExitCompBnd, branchInitsolCompBnd, branchExitsolCompBnd,
         branchExeclpCompBnd, branchExecextCompBnd, branchExecpsCompBnd,
         branchruledata) );

   /* add compbnd branching rule parameters */
   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/compbnd/useMRMH",
      "should the max range midrange heuristic be used?", &branchruledata->useMaxRangeMidrangeHeuristic, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(origscip, "branching/compbnd/useMDMH",
      "should the most distinct median heuristic be used?", &branchruledata->useMostDistinctMedianHeuristic, FALSE, TRUE, NULL, NULL) );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   SCIP_CALL( GCGconsIntegralorigAddBranchrule(scip, branchrule) );

   return SCIP_OKAY;
}

/** returns true when the branch rule is the component bound branching rule */
SCIP_Bool GCGisBranchruleCompBnd(
   SCIP_BRANCHRULE*      branchrule          /**< branchrule to check */
)
{
   return (branchrule != NULL) && (strcmp(BRANCHRULE_NAME, SCIPbranchruleGetName(branchrule)) == 0);
}
