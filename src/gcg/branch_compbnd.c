/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       */
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

/**@file   branch_compbnd.c
 *
 * @brief  component bound branching rule
 * @author Til Mohr
 * @author Marcel Schmickerath
 * @author Martin Bergner
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#define SCIP_DEBUG
#include <scip/pub_message.h>
#include <scip/pub_var.h>
#include <scip/type_misc.h>
#include <assert.h>
#include <scip/def.h>
#include <scip/cons_linear.h>
#include <scip/scip.h>
#include <scip/type_result.h>
#include <scip/type_retcode.h>
#include <scip/type_scip.h>
#include <scip/type_var.h>
#include <string.h>
#include <scip/struct_var.h>

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
#define BRANCHRULE_PRIORITY        -10                              /**< priority of this branching rule */
#define BRANCHRULE_MAXDEPTH        -1                             /**< maximal depth level of the branching rule */
#define BRANCHRULE_MAXBOUNDDIST    1.0                            /**< maximal relative distance from current node's
                                                                   dual bound to primal bound compared to best node's
                                                                   dual bound for applying branching */

/*
 * Data structures
 */

/** branching data */
struct GCG_BranchData
{
   GCG_BRANCH_TYPE       branchtype;         /**< type of branching */
   SCIP_Real             constant;           /**< constant value of the branching constraint in the master problem - either lhs or rhs, depending on branchtype */
   GCG_COMPBND*          B;                  /**< component bound sequence which induce the current branching constraint */
   int                   Bsize;              /**< size of the component bound sequence B */
   int                   blocknr;            /**< id of the pricing problem (or block) to which this branching constraint belongs */
   int                   nvars;              /**< number of master variables the last time the node has been visited - neccessary to later include newly generated master variables */
   GCG_MASTERCUTDATA*    mastercons;         /**< master constraint along with its corresponding inferred pricing modifications */
};

/*
 * Local methods
 */

/** initialize branchdata at the node */
static
SCIP_RETCODE initNodeBranchdata(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_BRANCHDATA**      nodebranchdata,     /**< branching data to set */
   GCG_BRANCH_TYPE       branchtype,         /**< type of branch to generate */
   SCIP_Real             constant,           /**< constant value of the branching constraint in the master problem - either lhs or rhs, depending on branchtype */
   GCG_COMPBND*          B,                  /**< component bound sequence which induce the current branching constraint */
   int                   Bsize,              /**< size of the component bound sequence B */
   int                   blocknr             /**< block we are branching in */
   )
{
   SCIP_CALL( SCIPallocBlockMemory(scip, nodebranchdata) );

   (*nodebranchdata)->branchtype = branchtype;
   (*nodebranchdata)->blocknr = blocknr;
   (*nodebranchdata)->mastercons = NULL;
   (*nodebranchdata)->constant = constant;
   (*nodebranchdata)->B = B;
   (*nodebranchdata)->Bsize = Bsize;
   (*nodebranchdata)->nvars = 0;
   (*nodebranchdata)->mastercons = NULL;

   return SCIP_OKAY;
}

/** initialize B based on B&B-tree ancestors */
static
SCIP_RETCODE initComponentBoundsFromAncestors(
   SCIP*                 masterscip,              /**< SCIP data structure */
   GCG_COMPBND**         B,                       /**< Component Bound Sequence defining the nodes */
   int*                  Bsize,                   /**< size of B */
   int                   blocknr                  /**< number of the block */
   )
{
   SCIP_CONS* parentcons;
   GCG_BRANCHDATA* parentdata;

   parentcons = GCGconsMasterbranchGetActiveCons(masterscip);

   while(parentcons != NULL) {
      if(GCGconsMasterbranchGetBranchrule(parentcons) == NULL
         || strcmp(SCIPbranchruleGetName(GCGconsMasterbranchGetBranchrule(parentcons)), BRANCHRULE_NAME) != 0)
      {
         parentcons = GCGconsMasterbranchGetParentcons(parentcons);
         continue;
      }

      parentdata = GCGconsMasterbranchGetBranchdata(parentcons);
      assert(parentdata != NULL);

      if(blocknr != parentdata->blocknr)
      {
         parentcons = GCGconsMasterbranchGetParentcons(parentcons);
         continue;
      }

      // resize B to the size of the parent's B
      if(*Bsize == 0)
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, B, parentdata->Bsize) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, B, *Bsize, parentdata->Bsize) );
      }

      // copy the parent's B into B
      for(int i = 0; i < parentdata->Bsize; ++i)
      {
         (*B)[i] = parentdata->B[i];
      }

      *Bsize = parentdata->Bsize;

      break;
   }

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
         if( SCIPvarGetName(component) != SCIPvarGetName((*B)[j].component) )
            continue;
         if( (*B)[i].sense != (*B)[j].sense )
            continue;

         already_added[j] = TRUE;

         if( (*B)[i].sense == GCG_COMPBND_SENSE_LE )
         {
            bound = MIN((*B)[j].bound, (*B)[i].bound);
         } else if ( (*B)[i].sense == GCG_COMPBND_SENSE_GE )
         {
            bound = MAX((*B)[j].bound, (*B)[i].bound);
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

   SCIPdebugMessage("Simplified B from %d to %d\n", *Bsize, newBsize);
   assert(1 <= newBsize && newBsize <= *Bsize);

   SCIPfreeBlockMemoryArrayNull(scip, &already_added, *Bsize);

   // free old B and replace it with newB
   SCIPfreeBlockMemoryArrayNull(scip, B, *Bsize);
   *B = newB;
   *Bsize = newBsize;

   return SCIP_OKAY;
}

/** computes the generator of mastervar for the entry in origvar
 * @return entry of the generator corresponding to origvar */
static
SCIP_Real getGeneratorEntry(
   SCIP_VAR*             mastervar,          /**< current mastervariable */
   SCIP_VAR*             origvar             /**< corresponding origvar */
   )
{
   int i;
   SCIP_VAR** origvars;
   SCIP_Real* origvals;
   int norigvars;

   assert(mastervar != NULL);
   assert(origvar != NULL);

   origvars = GCGmasterVarGetOrigvars(mastervar);
   norigvars = GCGmasterVarGetNOrigvars(mastervar);
   origvals = GCGmasterVarGetOrigvals(mastervar);

   for( i = 0; i < norigvars; ++i )
   {
      if( origvars[i] == origvar )
      {
         return origvals[i];
      }
   }

   return 0.0;
}

/** whether a master variable is in B or not */
static
SCIP_Bool isMasterVarInB(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             mastervar,          /**< master variable to check */
   GCG_COMPBND*          B,                  /**< component bound sequence to check */
   int                   Bsize               /**< size of B */
   )
{
   int i;

   assert(mastervar != NULL);
   assert(Bsize > 0);
   assert(B != NULL);

   for( i = 0; i < Bsize; ++i )
   {
      SCIP_Real generatorentry = getGeneratorEntry(mastervar, B[i].component);
      if ( (B[i].sense == GCG_COMPBND_SENSE_GE && SCIPisLT(scip, generatorentry, SCIPfloor(scip, B[i].bound) + 1.0)) ||
           (B[i].sense == GCG_COMPBND_SENSE_LE && SCIPisGT(scip, generatorentry, SCIPfloor(scip, B[i].bound))) )
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
   SCIP_CONS* branchingcons;
#ifdef SCIP_DEBUG
   SCIP_VAR** mastervars;
   int nmastervars;
   int i;
   SCIP_Real constantSum;
#endif

   assert(scip != NULL);
   assert(mastervar != NULL);
   assert(branchdata != NULL);
   assert(added != NULL);

   *added = FALSE;

   if( GCGvarGetBlock(mastervar) == -1 || !GCGisMasterVarInBlock(mastervar, branchdata->blocknr) )
      return SCIP_OKAY;

   if( !isMasterVarInB(scip, mastervar, branchdata->B, branchdata->Bsize) )
      return SCIP_OKAY;

   branchingcons = NULL;
   SCIP_CALL( GCGmastercutGetCons(branchdata->mastercons, &branchingcons) );
   assert(branchingcons != NULL);
   SCIPdebugMessage("Adding variable %s to branching constraint: %s%.2f\n", SCIPvarGetName(mastervar), branchdata->branchtype == GCG_BRANCH_DOWN ? "<=" : ">=", branchdata->constant);
   SCIP_CALL( SCIPaddCoefLinear(scip, branchingcons, mastervar, 1.0) );
   *added = TRUE;

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(nmastervars >= 0);
   /* determine the constant value of the master constraint, i.e. the rhs for the down branch, and lhs for the up branch */
   constantSum = 0;
   for( i = 0; i < nmastervars; ++i )
   {
      if( GCGisMasterVarInBlock(mastervars[i], branchdata->blocknr) && isMasterVarInB(scip, mastervars[i], branchdata->B, branchdata->Bsize) )
      {
         constantSum += SCIPgetSolVal(scip, NULL, mastervars[i]);
      }
   }
   SCIPdebugMessage("Constant sum: %g\n", constantSum);
#endif

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
   // get the dual value of the master constraint
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
   GCG_PRICINGMODIFICATION* pricingmod;
   GCG_PRICINGMODIFICATION** pricingmods; // will always only contain one element in this branching rule

   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(node != NULL);
   assert(branchdata != NULL);

   assert(branchdata->mastercons == NULL);

   mastervars = SCIPgetVars(scip);
   nmastervars = SCIPgetNVars(scip);
   assert(nmastervars >= 0);
   assert(nmastervars == 0 || mastervars != NULL);

   /*  create constraint for child */
   if (branchdata->branchtype == GCG_BRANCH_DOWN)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "child(%d <= %g)", branchdata->Bsize, branchdata->constant);
      SCIP_CALL( SCIPcreateConsLinear(scip, &branchcons, name, 0, NULL, NULL,
         -SCIPinfinity(scip), branchdata->constant, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
   }
   else
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "child(%d >= %g)", branchdata->Bsize, branchdata->constant);
      SCIP_CALL( SCIPcreateConsLinear(scip, &branchcons, name, 0, NULL, NULL,
         branchdata->constant, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
   }

   SCIP_CALL( SCIPaddConsNode(scip, node, branchcons, NULL) );

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);
   pricingscip = GCGgetPricingprob(origscip, branchdata->blocknr);
   assert(pricingscip != NULL);

   assert(branchdata->mastercons == NULL);

   if( branchdata->branchtype == GCG_BRANCH_DOWN )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "down(%s%s%.2f,%d)<=%.2f",
            SCIPvarGetName(branchdata->B[0].component), branchdata->B[0].sense == GCG_COMPBND_SENSE_GE ? ">=" : "<=",
                           SCIPfloor(scip, branchdata->B[i].bound) + (branchdata->B[0].sense == GCG_COMPBND_SENSE_GE ? 1.0 : 0.0),
                           branchdata->Bsize, branchdata->constant);
   }
   else
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "up(%s%s%.2f,%d)>=%.2f",
            SCIPvarGetName(branchdata->B[0].component), branchdata->B[0].sense == GCG_COMPBND_SENSE_GE ? ">=" : "<=",
                           SCIPfloor(scip, branchdata->B[i].bound) + (branchdata->B[0].sense == GCG_COMPBND_SENSE_GE ? 1.0 : 0.0),
                           branchdata->Bsize, branchdata->constant);
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
      (void) SCIPsnprintf(pricingvarname, SCIP_MAXSTRLEN, "y(%s,%s,%f.2)", SCIPvarGetName(branchdata->B[i].component), branchdata->B[i].sense == GCG_COMPBND_SENSE_GE ? ">=" : "<=",
                           SCIPfloor(scip, branchdata->B[i].bound) + (branchdata->B[i].sense == GCG_COMPBND_SENSE_GE ? 1.0 : 0.0));
      SCIP_CALL( GCGcreateInferredPricingVar(pricingscip, &additionalvars[i], pricingvarname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, branchdata->blocknr) );
   }

   // create the pricing constraints
   if( branchdata->branchtype == GCG_BRANCH_DOWN )
   {
      /* g_x >= 1 + sum_{j=1}^{n} y_j - n
         * 1 - n <= g_x - sum_{j=1}^{n} y_j */
      nadditionalcons = branchdata->Bsize + 1;
      SCIP_CALL( SCIPallocBlockMemoryArray(pricingscip, &additionalcons, nadditionalcons) );
      char consname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "c(%s)", pricingvarname);
      SCIP_CALL( SCIPcreateConsLinear(pricingscip, &additionalcons[0], consname, 0, NULL, NULL, 1.0 - branchdata->Bsize, SCIPinfinity(pricingscip),
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) );
      SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[0], coefvar, 1.0) );
      for( i = 0; i < branchdata->Bsize; ++i)
      {
         SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[0], additionalvars[i], -1.0) );
      }

      for( i = 0; i < branchdata->Bsize; ++i)
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "c(%s)(%s,%s,%f.2)", name, SCIPvarGetName(branchdata->B[i].component), branchdata->B[i].sense == GCG_COMPBND_SENSE_GE ? ">=" : "<=",
                              SCIPfloor(scip, branchdata->B[i].bound) + (branchdata->B[i].sense == GCG_COMPBND_SENSE_GE ? 1.0 : 0.0));

         if( branchdata->B[i].sense == GCG_COMPBND_SENSE_LE )
         {
            /* y_j >= (floor(bound) + 1 - x_j) / (floor(bound) + 1 - l_j)
               * (floor(bound) + 1 - l_j) * y_j >= floor(bound) + 1 - x_j
               * floor(bound) + 1 <= (floor(bound) + 1 - l_j) * y_j + x_j */
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->B[i].component);
            SCIP_Real bound = SCIPfloor(pricingscip, branchdata->B[i].bound);
            SCIP_Real lowerbound = SCIPvarGetLbLocal(pricing_var);
            assert(SCIPisPositive(pricingscip, bound + 1.0 - lowerbound));
            SCIP_CALL( SCIPcreateConsLinear(pricingscip, &additionalcons[i+1], consname, 0, NULL, NULL, bound + 1.0, SCIPinfinity(pricingscip),
               TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) );
            SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[i+1], additionalvars[i], bound + 1.0 - lowerbound) );
            SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[i+1], pricing_var, 1.0) );
         }
         else
         {
            /* y_j >= (x_j - floor(bound)) / (u_j - floor(bound))
               * (u_j - floor(bound)) * y_j >= x_j - floor(bound)
               * -floor(bound) <= (u_j - floor(bound)) * y_j - x_j */
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->B[i].component);
            SCIP_Real bound = SCIPfloor(pricingscip, branchdata->B[i].bound);
            SCIP_Real upperbound = SCIPvarGetUbLocal(pricing_var);
            assert(SCIPisPositive(pricingscip, upperbound - bound));
            SCIP_CALL( SCIPcreateConsLinear(pricingscip, &additionalcons[i+1], consname, 0, NULL, NULL, -bound, SCIPinfinity(pricingscip),
               TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) );
            SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[i+1], additionalvars[i], upperbound - bound) );
            SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[i+1], pricing_var, -1.0) );
         }
      }
   }
   else
   {
      /* g_x <= y_j
         * g_x - y_j <= 0 */
      nadditionalcons = branchdata->Bsize * 2;
      SCIP_CALL( SCIPallocBlockMemoryArray(pricingscip, &additionalcons, nadditionalcons) );
      char consname[SCIP_MAXSTRLEN];
      for( i = 0; i < branchdata->Bsize; ++i )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "c0(%s)(%s,%s,%f.2)", name, SCIPvarGetName(branchdata->B[i].component), branchdata->B[i].sense == GCG_COMPBND_SENSE_GE ? ">=" : "<=",
                              SCIPfloor(scip, branchdata->B[i].bound) + (branchdata->B[i].sense == GCG_COMPBND_SENSE_GE ? 1.0 : 0.0));
         SCIP_CALL( SCIPcreateConsLinear(pricingscip, &additionalcons[i], consname, 0, NULL, NULL, -SCIPinfinity(pricingscip), 0.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) );
         SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[i], coefvar, 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[i], additionalvars[i], -1.0) );
      }

      for( i = 0; i < branchdata->Bsize; ++i )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "c1(%s)(%s,%s,%f.2)", name, SCIPvarGetName(branchdata->B[i].component), branchdata->B[i].sense == GCG_COMPBND_SENSE_GE ? ">=" : "<=",
                              SCIPfloor(scip, branchdata->B[i].bound) + (branchdata->B[i].sense == GCG_COMPBND_SENSE_GE ? 1.0 : 0.0));

         if( branchdata->B[i].sense == GCG_COMPBND_SENSE_LE )
         {
            /* y_j <= (u_j - x_j) / (u_j - floor(bound))
               * (u_j - floor(bound)) * y_j <= u_j - x_j
               * (u_j - floor(bound)) * y_j + x_j <= u_j */
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->B[i].component);
            SCIP_Real bound = SCIPfloor(pricingscip, branchdata->B[i].bound);
            SCIP_Real upperbound = SCIPvarGetUbLocal(pricing_var);
            assert(SCIPisPositive(pricingscip, upperbound - bound));
            SCIP_CALL( SCIPcreateConsLinear(pricingscip, &additionalcons[i+branchdata->Bsize], consname, 0, NULL, NULL, -SCIPinfinity(pricingscip), upperbound,
               TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) );
            SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[i+branchdata->Bsize], additionalvars[i], upperbound - bound) );
            SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[i+branchdata->Bsize], pricing_var, -1.0) );
         }
         else
         {
            /* y_j <= (x_j - l_j) / (floor(bound) + 1 - l_j)
               * (floor(bound) + 1 - l_j) * y_j <= x_j - l_j
               * (floor(bound) + 1 - l_j) * y_j - x_j <= -l_j */
            SCIP_VAR* pricing_var = GCGoriginalVarGetPricingVar(branchdata->B[i].component);
            SCIP_Real bound = SCIPfloor(pricingscip, branchdata->B[i].bound);
            SCIP_Real lowerbound = SCIPvarGetLbLocal(pricing_var);
            assert(SCIPisPositive(pricingscip, bound + 1.0 - lowerbound));
            SCIP_CALL( SCIPcreateConsLinear(pricingscip, &additionalcons[i+branchdata->Bsize], consname, 0, NULL, NULL, -SCIPinfinity(pricingscip), -lowerbound,
               TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) );
            SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[i+branchdata->Bsize], additionalvars[i], bound + 1.0 - lowerbound) );
            SCIP_CALL( SCIPaddCoefLinear(pricingscip, additionalcons[i+branchdata->Bsize], pricing_var, -1.0) );
         }
      }
   }

   // create the pricing modification
   pricingmod = NULL;
   SCIP_CALL( GCGpricingmodificationCreate(
      scip,
      &pricingmod,
      branchdata->blocknr,
      coefvar,
      additionalvars,
      nadditionalvars,
      additionalcons,
      nadditionalcons
   ) );
   assert(pricingmod != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &pricingmods, 1) );
   assert(pricingmods != NULL);
   pricingmods[0] = pricingmod;

   // create the master constraint
   SCIP_CALL( GCGmastercutCreateFromCons(scip, &branchdata->mastercons, branchcons, pricingmods, 1) );

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
   constantSum = 0;
   for( i = 0; i < nmastervars; ++i )
   {
      if( GCGisMasterVarInBlock(mastervars[i], blocknr) && isMasterVarInB(scip, mastervars[i], B, Bsize) )
      {
         constantSum += SCIPgetSolVal(masterscip, NULL, mastervars[i]);
      }
   }
   // sanity check: the sum must be fractional, otherwise something went wrong during separation, i.e. B is incorrect
   assert(!SCIPisFeasIntegral(scip, constantSum));

   /* create two nodes */
   SCIPdebugMessage("Component bound branching rule: creating 2 nodes\n");
   SCIP_CALL( initNodeBranchdata(scip, &downBranchData, GCG_BRANCH_DOWN, SCIPfloor(scip, constantSum), B, Bsize, blocknr) );
   SCIP_CALL( initNodeBranchdata(scip, &upBranchData, GCG_BRANCH_UP, SCIPceil(scip, constantSum), B, Bsize, blocknr) );
   assert(downBranchData != NULL);
   assert(upBranchData != NULL);

   /* define names for origbranch constraints */
   (void) SCIPsnprintf(downChildname, SCIP_MAXSTRLEN, "node(%d) (last comp=%s %s %g) <= %g", blocknr,
      SCIPvarGetName(downBranchData->B[downBranchData->Bsize-1].component),
      downBranchData->B[downBranchData->Bsize-1].sense == GCG_COMPBND_SENSE_GE? ">=": "<=",
      SCIPfloor(scip, downBranchData->B[downBranchData->Bsize-1].bound) + (downBranchData->B[downBranchData->Bsize-1].sense == GCG_COMPBND_SENSE_GE ? 1.0 : 0.0),
      downBranchData->constant);
   (void) SCIPsnprintf(upChildname, SCIP_MAXSTRLEN, "node(%d) (last comp=%s %s %g) >= %g", blocknr,
      SCIPvarGetName(upBranchData->B[upBranchData->Bsize-1].component),
      upBranchData->B[upBranchData->Bsize-1].sense == GCG_COMPBND_SENSE_GE? ">=": "<=",
      SCIPfloor(scip, upBranchData->B[upBranchData->Bsize-1].bound) + (upBranchData->B[upBranchData->Bsize-1].sense == GCG_COMPBND_SENSE_GE ? 1.0 : 0.0),
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
   assert(mastervars != NULL);
   assert(nmastervars > 0);

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
   SCIP_VAR**            X,                  /**< array of mastervars */
   int                   Xsize,              /**< number of mastervars */
   SCIP_VAR***           indexSet,           /**< set to initialize */
   int*                  indexSetSize        /**< size of the index set */
   )
{
   int i;

   assert( scip != NULL);
   assert( X != NULL);
   assert( Xsize > 0);
   assert( indexSet != NULL);
   assert( indexSetSize != NULL);

   *indexSet = NULL;
   *indexSetSize = 0;


   for( i = 0; i < Xsize; ++i )
   {
      int j;
      SCIP_VAR** origvars = GCGmasterVarGetOrigvars(X[i]);
      int norigvars = GCGmasterVarGetNOrigvars(X[i]);

      if( *indexSetSize == 0 && norigvars > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, indexSet, norigvars) );
         for( j = 0; j < norigvars; ++j )
         {

            (*indexSet)[*indexSetSize] = origvars[j];
            ++(*indexSetSize);
         }
      }
      else
      {
         for( j = 0; j < norigvars; ++j )
         {
            int k;
            int oldsize = *indexSetSize;
            SCIP_Bool found = FALSE;

            // Make sure to not include duplicates in the set
            for( k = 0; k < oldsize; ++k )
            {
               /*  if variable already in union */
               if( (*indexSet)[k] == origvars[j] )
               {
                  found = TRUE;
                  break;
               }
            }
            if( !found )
            {
               /*  add variable to the end */
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, indexSet, *indexSetSize, *indexSetSize+1) );
               (*indexSet)[*indexSetSize] = origvars[j];
               ++(*indexSetSize);
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** recursive helper for the separation algorithm determining a component bound sequence to branch on */
static
SCIP_RETCODE _separation(
   SCIP*                 masterscip,              /**< SCIP data structure */
   SCIP_VAR**            X,                       /**< mastervariables currently satisfying the component bound sequence in the specified block */
   int                   Xsize,                   /**< size of X */
   GCG_COMPBND**         B,                       /**< Component Bound Sequence defining the nodes */
   int*                  Bsize,                   /**< size of B */
   SCIP_RESULT*          result                   /**< pointer to store the result of the branching call */
   )
{
   SCIP_Real fractionality;

   // Sanity check: All variables in X must satisfy the component bound sequence
   if( *Bsize > 0 ) {
      for( int i = 0; i < Xsize; ++i )
      {
         assert(isMasterVarInB(masterscip, X[i], *B, *Bsize));
      }
   }

   fractionality = calcFractionality(masterscip, X, Xsize);
   assert(fractionality >= 0.0);

   if( SCIPisEQ(masterscip, fractionality, 0.0) )
   {
      // all variables are integral, nothing to do
      SCIPdebugMessage("All variables are integral, nothing to do\n");
      return SCIP_OKAY;
   }

   if( !SCIPisFeasIntegral(masterscip, fractionality) )
   {
      /* now we branch on <= floor(fractionality) and >= ceil(fractionality)
       * this is handled by the createChildNodesCompBnd method
       */
       SCIPdebugMessage("fractionality is fractional, we need to branch\n");
       return SCIP_OKAY;
   }

   // the fractionality is integral, we need to impose an additional bound
   // Task: Find fractional mastervarMin, mastervarMax, indexed by x1, x2, s.t. for some j, xj1 < xj2
   SCIP_VAR* currMastervar;
   SCIP_Real currSolutionValue; // of the currMastervar
   SCIP_VAR* origvar;
   SCIP_VAR** indexSet = NULL;
   int indexSetSize = 0;
   SCIP_Real min;
   SCIP_Real max;
   int i;
   int j;
   SCIP_Real generatorentry;
   SCIP_Bool found = FALSE;

   SCIP_CALL( initIndexSet(masterscip, X, Xsize, &indexSet, &indexSetSize) );

   for( j=0; j<indexSetSize; ++j )
   {
      min = SCIPinfinity(masterscip);
      max = -SCIPinfinity(masterscip);
      SCIPdebugMessage("j=%d\n", j);
      SCIPdebugMessage("INIT min=%f, max=%f\n", min, max);
      origvar = indexSet[j];
      for( i=0; i<Xsize; ++i )
      {
         currMastervar = X[i];
         // Current solution value of the master variable must be fractional and > 0
         currSolutionValue = SCIPgetSolVal(masterscip, NULL, currMastervar);
         if ( SCIPisFeasIntegral(masterscip, currSolutionValue) || !SCIPisPositive(masterscip, currSolutionValue) ) {
            // skip
            continue;
         }

         generatorentry = getGeneratorEntry(currMastervar, origvar);
         // first, check if we can update the min and max values
         if( SCIPisLT(masterscip, generatorentry, min) )
            min = generatorentry;
         else if( SCIPisLT(masterscip, max, generatorentry) )
            max = generatorentry;

         // now could be that min<max
         if( SCIPisLT(masterscip, min, max) )
         {
            found=TRUE;
         }
      }
      if (found)
         break;
   }
   SCIPdebugMessage("Found %d for origvar %s, j=%d, min=%f, max=%f\n", found, SCIPvarGetName(origvar), j, min, max);
   assert(found); // Must exist, otherwise we have a bug in the code
   assert(min<max);

   // j and origivar are still set to the correct values after execution of the loop
   SCIP_Real value = (min+max)/2;

   // create two component bound options, xj <= floor(value) and xj >= floor(value)+1
   GCG_COMPBND* B1;
   GCG_COMPBND* B2;
   int new_Bsize = *Bsize + 1;
   if( *Bsize == 0 )
   {
      // allocate memory for B1 and B2
      SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &B1, new_Bsize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &B2, new_Bsize) );
   } else
   {
      // copy the current component bound sequence B (of type GCG_COMPBND**) into B1 and B2, increasing the size by 1
      SCIP_CALL( SCIPduplicateBlockMemoryArray(masterscip, &B1, *B, new_Bsize) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(masterscip, &B2, *B, new_Bsize) );
   }
   // add the new component bounds to B1 and B2
   B1[*Bsize] = (GCG_COMPBND){origvar, GCG_COMPBND_SENSE_LE, value};
   B2[*Bsize] = (GCG_COMPBND){origvar, GCG_COMPBND_SENSE_GE, value};

   SCIPdebugMessage("B1 and B2 after adding the new component bounds\n");
   for (i = 0; i < new_Bsize; ++i)
   {
      SCIPdebugMessage("B1[%d]: %s %s %f\n", i, SCIPvarGetName(B1[i].component),
                        B1[i].sense == GCG_COMPBND_SENSE_LE ? "<=" : ">=",
                        SCIPfloor(masterscip, B1[i].bound) + (B1[i].sense == GCG_COMPBND_SENSE_GE ? 1.0 : 0.0));
      SCIPdebugMessage("B2[%d]: %s %s %f\n", i, SCIPvarGetName(B2[i].component),
                        B2[i].sense == GCG_COMPBND_SENSE_LE ? "<=" : ">=",
                        SCIPfloor(masterscip, B2[i].bound) + (B2[i].sense == GCG_COMPBND_SENSE_GE ? 1.0 : 0.0));
   }

   // assign the master variables to X1 and X2, depending on whether they satisfy the new component bound sequences
   SCIP_VAR** X1 = NULL;
   SCIP_VAR** X2 = NULL;
   int X1size = 0;
   int X2size = 0;
   int x;
   for( x=0; x<Xsize; ++x )
   {
#ifdef SCIP_DEBUG
      SCIP_Bool inB1 = FALSE;
      SCIP_Bool inB2 = FALSE;
#endif
      if( isMasterVarInB(masterscip, X[x], B1, new_Bsize) )
      {
         // increase the size of X1 by 1, and add the current variable to X1
         SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, &X1, X1size, X1size+1) );
         X1[X1size] = X[x];
         X1size += 1;
#ifdef SCIP_DEBUG
         inB1 = TRUE;
#endif
      }
      if( isMasterVarInB(masterscip, X[x], B2, new_Bsize) )
      {
         // increase the size of X2 by 1, and add the current variable to X2
         SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, &X2, X2size, X2size+1) );
         X2[X2size] = X[x];
         X2size += 1;
#ifdef SCIP_DEBUG
         inB2 = TRUE;
#endif
      }
#ifdef SCIP_DEBUG
      assert(inB1 != inB2); // one or the other must hold
#endif
   }
   SCIPdebugMessage("X1size: %d, X2size: %d\n", X1size, X2size);
   assert(X1size > 0 && X2size > 0);

   // determine the fractionality of B1 and B2
   SCIP_Real fractionality1 = calcFractionality(masterscip, X1, X1size);
   SCIP_Real fractionality2 = calcFractionality(masterscip, X2, X2size);

   assert(!SCIPisGT(masterscip, fractionality1, fractionality / 2) || !SCIPisGT(masterscip, fractionality2, fractionality / 2));
   assert(
      (SCIPisFeasIntegral(masterscip, fractionality1) && SCIPisFeasIntegral(masterscip, fractionality2)) ||
      (!SCIPisFeasIntegral(masterscip, fractionality1) && !SCIPisFeasIntegral(masterscip, fractionality2))
   );

   if( SCIPisFeasIntegral(masterscip, fractionality1) && SCIPisFeasIntegral(masterscip, fractionality2) )
   {
      // both are integral, recursively branch on the B with smallest fractionality
      if( SCIPisLT(masterscip, fractionality1, fractionality2) )
      {
         // free B2 early to save memory
         SCIPfreeBlockMemoryArray(masterscip, &B2, new_Bsize);

         // recursive call
         SCIP_CALL( _separation(masterscip, X1, X1size, &B1, &new_Bsize, result) );

         // copy the new component bound sequence B1 into B, and free B1
         SCIPreallocBlockMemoryArray(masterscip, B, *Bsize, new_Bsize);
         for( i = 0; i < new_Bsize; ++i )
         {
            (*B)[i] = B1[i];
         }
         *Bsize = new_Bsize;
         SCIPfreeBlockMemoryArray(masterscip, &B1, new_Bsize);
      }
      else
      {
         // free B1 early to save memory
         SCIPfreeBlockMemoryArray(masterscip, &B1, new_Bsize);

         // recursive call
         SCIP_CALL( _separation(masterscip, X2, X2size, &B2, &new_Bsize, result) );

         // copy the new component bound sequence B2 into B, and free B2
         SCIPreallocBlockMemoryArray(masterscip, B, *Bsize, new_Bsize);
         for( i = 0; i < new_Bsize; ++i )
         {
            (*B)[i] = B2[i];
         }
         *Bsize = new_Bsize;
         SCIPfreeBlockMemoryArray(masterscip, &B2, new_Bsize);
      }
   } else
   {
      // both are fractional, select the one with greater Xsize, and return
      if( X1size > X2size )
      {
         // copy the new component bound sequence B1 into B, and free B1 and B2
         SCIPreallocBlockMemoryArray(masterscip, B, *Bsize, new_Bsize);
         for( i = 0; i < new_Bsize; ++i )
         {
            (*B)[i] = B1[i];
         }
         *Bsize = new_Bsize;
         SCIPfreeBlockMemoryArray(masterscip, &B1, new_Bsize);
         SCIPfreeBlockMemoryArray(masterscip, &B2, new_Bsize);
      }
      else
      {
         // copy the new component bound sequence B2 into B, and free B1 and B2
         SCIPreallocBlockMemoryArray(masterscip, B, *Bsize, new_Bsize);
         for( i = 0; i < new_Bsize; ++i )
         {
            (*B)[i] = B2[i];
         }
         *Bsize = new_Bsize;
         SCIPfreeBlockMemoryArray(masterscip, &B1, new_Bsize);
         SCIPfreeBlockMemoryArray(masterscip, &B2, new_Bsize);
      }
   }

   return SCIP_OKAY;
}

/** separation algorithm determining a component bound sequence to branch on
 *
 * @param[out] B the component bound sequence to branch on
 * @param[out] Bsize the size of B
 */
static
SCIP_RETCODE separation(
   SCIP*                 masterscip,              /**< SCIP data structure */
   GCG_COMPBND**         B,                       /**< Component Bound Sequence defining the nodes */
   int*                  Bsize,                   /**< size of B */
   int                   blocknr,                 /**< number of the block */
   SCIP_RESULT*          result                   /**< pointer to store the result of the branching call */
   )
{
   SCIP* origscip;
   SCIP_VAR** mastervars;
   int nmastervars;
   SCIP_VAR** X;
   int Xsize;

   assert(*Bsize == 0 || *B != NULL);

   origscip = GCGmasterGetOrigprob(masterscip);
   mastervars = NULL;
   nmastervars = 0;
   X = NULL;
   Xsize = 0;

   assert(masterscip != NULL);

   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   /* 1. Set X to include all mastervariables in the selected block */
   for( int i = 0; i < nmastervars; ++i )
   {
      if( !GCGisMasterVarInBlock(mastervars[i], blocknr) )
         continue;

      if( *Bsize > 0 && !isMasterVarInB(masterscip, mastervars[i], *B, *Bsize) )
         continue;

      if( Xsize == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(origscip, &X, 1) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(origscip, &X, Xsize, Xsize + 1) );
      }
      X[Xsize] = mastervars[i];
      Xsize++;
   }

   /* 2. Call the recursive separation algorithm */
   SCIP_CALL( _separation(masterscip, X, Xsize, B, Bsize, result) );
   SCIPdebugMessage("Component bound sequence to branch on:\n");
   SCIPdebugMessage("Bsize: %d\n", *Bsize);
   for (int i = 0; i < *Bsize; ++i)
   {
      SCIPdebugMessage("B[%d]: %s %s %f\n", i, SCIPvarGetName((*B)[i].component),
                        (*B)[i].sense == GCG_COMPBND_SENSE_LE ? "<=" : ">=",
                        SCIPfloor(masterscip, (*B)[i].bound) + ((*B)[i].sense == GCG_COMPBND_SENSE_GE ? 1.0 : 0.0));
   }

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
   SCIP_VAR** branchcands;
   SCIP_VAR** allorigvars;
   SCIP_VAR** mastervars;
   int nmastervars;
   int nbranchcands;
   SCIP_VAR* mastervar;
   GCG_COMPBND* B;
   int Bsize;
   int blocknr;
   int i;
   int j;
   int allnorigvars;

   blocknr = -2;
   B = NULL;
   Bsize = 0;

   assert(masterscip != NULL);

   SCIPdebugMessage("get information for component bound branching\n");

   origscip = GCGmasterGetOrigprob(masterscip);

   assert(origscip != NULL);
   SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL, NULL) );

   SCIP_CALL( SCIPgetVarsData(origscip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   /* in case original problem contains continuous variables, there are no branching cands */
   assert(nbranchcands > 0 || SCIPgetNContVars(origscip) > 0);
   mastervar = NULL;

   /* 1. Determine in what block we are branching. We select the first available block,
    *     i.e. the first block that contains a branching candidate, starting from the master block.
    */
   /* in case of continuous origvar look for "fractional" blocks using the representation (currentorigsol) in the original problem */
   if(SCIPgetNContVars(origscip) > 0)
   {
      int norigvars;
      SCIP_VAR** origvars;
      SCIP_VAR* origvar;

      norigvars = SCIPgetNVars(origscip);
      origvars = SCIPgetVars(origscip);

      nbranchcands = SCIPgetNVars(masterscip);
      branchcands = SCIPgetVars(masterscip);

      assert(nbranchcands > 0);

      for( i = 0; i < norigvars; ++i )
      {
         origvar = origvars[i];

         if( SCIPvarGetType(origvar) > SCIP_VARTYPE_INTEGER )
            continue;

         if( SCIPisIntegral(origscip, SCIPgetSolVal(origscip, GCGrelaxGetCurrentOrigSol(origscip), origvar)) )
            continue;

         blocknr = GCGgetBlockRepresentative(origscip, GCGvarGetBlock(origvar));

         SCIPdebugMessage("Variable %s belonging to block %d with representative %d is not integral!\n", SCIPvarGetName(origvar), GCGvarGetBlock(origvar), blocknr);

         if( blocknr == -1 )
         {
            assert(GCGoriginalVarGetNMastervars(origvar) == 1);
            mastervar = GCGoriginalVarGetMastervars(origvar)[0];
            break;
         }

         break;
      }
   } else
   {
      /* loop over all branching candidates */
      for( i = 0; i < nbranchcands; ++i )
      {
         mastervar = branchcands[i];
         assert(GCGvarIsMaster(mastervar));

         /* if we have a master variable, we branch on it */
         if( GCGvarGetBlock(mastervar) == -1 )
         {
            assert(!GCGmasterVarIsArtificial(mastervar));
            blocknr = -1;
            break;
         }

         /* else, check if the candidate is in an procing block */
         for( j = 0; j < GCGgetNPricingprobs(origscip); ++j )
         {
            if( GCGisMasterVarInBlock(mastervar, j) )
            {
               blocknr = j;
               break;
            }
         }
      }
   }

   if( blocknr < -1 )
   {
      SCIPdebugMessage("Generic branching rule could not find variables to branch on!\n");
      SCIP_Bool rays;
      SCIP_CALL( GCGpricerExistRays(masterscip, &rays) );
      if( rays )
         SCIPwarningMessage(masterscip, "Generic branching is not compatible with unbounded problems!\n");
      return SCIP_ERROR;
   }

   // static mastervariables should be handled by the staticvar branchrule
   assert(blocknr >= -1);

   SCIPdebugMessage("branching in block %d \n", blocknr);

   /* 2. Check B&B-tree ancestors for previous compbnd branching in the node */
   SCIP_CALL( initComponentBoundsFromAncestors(masterscip, &B, &Bsize, blocknr) );

   /* 3. Call to separation algorithm to find a suitable B to branch on in the current block.*/
   SCIP_CALL( separation(masterscip, &B, &Bsize, blocknr, result) );
   assert(Bsize > 0);
   assert(B != NULL);

   /* 4. Remove component bounds that are strengthened by others */
   SCIP_CALL( simplifyComponentBounds(masterscip, &B, &Bsize) );
   assert(Bsize > 0);
   assert(B != NULL);

   /* 5. Create the child nodes. */
   SCIP_CALL( createChildNodesCompBnd(origscip, branchrule, B, Bsize, blocknr, result) );

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeCompBnd NULL


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
   assert(scip != NULL);
   assert(branchdata != NULL);

   if( *branchdata == NULL )
   {
      SCIPdebugMessage("branchDataDeleteCompBnd: cannot delete empty branchdata\n");

      return SCIP_OKAY;
   }

   /* release constraint that enforces the branching decision */
   if( (*branchdata)->mastercons != NULL )
   {
      SCIP_CALL( GCGmastercutFree(scip, &(*branchdata)->mastercons) );
   }

   if( (*branchdata)->B != NULL && (*branchdata)->Bsize > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &((*branchdata)->B), (*branchdata)->Bsize);
      (*branchdata)->B = NULL;
      (*branchdata)->Bsize = 0;
   }

   SCIPfreeBlockMemoryNull(scip, branchdata);
   *branchdata = NULL;

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

/** callback dual variable update method */
static
GCG_DECL_BRANCHUPDATEDUAL(branchUpdateDualCompBnd)
{
   SCIP* origscip;
   SCIP* pricingprob;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);
   assert(branchdata->mastercons->npricingmodifications == 1);
   assert(branchdata->mastercons->pricingmodifications[0] != NULL);
   assert(branchdata->mastercons->pricingmodifications[0]->coefvar != NULL);

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);
   pricingprob = GCGgetPricingprob(origscip, branchdata->blocknr);
   assert(pricingprob != NULL);

   assert((dual <= 0 && branchdata->branchtype == GCG_BRANCH_DOWN) || (dual >= 0 && branchdata->branchtype == GCG_BRANCH_UP));

   SCIP_CALL( SCIPchgVarObj(pricingprob, branchdata->mastercons->pricingmodifications[0]->coefvar, -dual) );

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

   origscip = GCGmasterGetOrigprob(scip);
   assert(branchrule != NULL);
   assert(origscip != NULL);

   SCIPdebugMessage("Init method of component bound branching\n");

   SCIP_CALL( GCGrelaxIncludeBranchrule(origscip, branchrule, branchActiveMasterCompBnd,
         branchDeactiveMasterCompBnd, branchPropMasterCompBnd, branchMasterSolvedCompBnd, branchDataDeleteCompBnd,
         branchNewColCompBnd, branchUpdateDualCompBnd, branchGetMastercutCompBnd) );

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
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create compbnd branching rule data */
   branchruledata = NULL;

   SCIPdebugMessage("Include method of component bound branching\n");

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
         BRANCHRULE_MAXBOUNDDIST,
         branchCopyCompBnd, branchFreeCompBnd, branchInitCompBnd, branchExitCompBnd, branchInitsolCompBnd, branchExitsolCompBnd,
         branchExeclpCompBnd, branchExecextCompBnd, branchExecpsCompBnd,
         branchruledata) );

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
