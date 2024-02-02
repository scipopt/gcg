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

/**@file    branch_compbnd.c
 *
 * @brief   branching rule based on vanderbeck's component bound branching
 * @author  Til Mohr
 * @author Marcel Schmickerath
 * @author Martin Bergner
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#define SCIP_DEBUG
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
#include "type_branchgcg.h"

#define BRANCHRULE_NAME            "compbnd"                      /**< name of branching rule */
#define BRANCHRULE_DESC            "component bound branching"    /**< short description of branching rule */
#define BRANCHRULE_PRIORITY        100                              /**< priority of this branching rule */
#define BRANCHRULE_MAXDEPTH        -1                             /**< maximal depth level of the branching rule */
#define BRANCHRULE_MAXBOUNDDIST    1.0                            /**< maximal relative distance from current node's
                                                                   dual bound to primal bound compared to best node's
                                                                   dual bound for applying branching */

#define EVENTHDLR_NAME         "compbndbranchvaradd"
#define EVENTHDLR_DESC         "event handler for adding a new generated mastervar into the right branching constraints by using component bound branching"

/*
 * Data structures
 */

typedef enum {
   GCG_BRANCH_DOWN = 0,
   GCG_BRANCH_UP = 1
} GCG_BRANCH_TYPE;

typedef enum {
   GCG_COMPSENSE_GE = 0,
   GCG_COMPSENSE_LE = 1
} GCG_COMPSENSE;

/** component bound structure */
struct ComponentBound
{
   SCIP_VAR*             component;          /**< variable to which this bound belongs */
   GCG_COMPSENSE         sense;              /**< sense of the bound */
   SCIP_Real             bound;              /**< bound value */
};
typedef struct ComponentBound GCG_COMPBND;

/** branching data */
struct GCG_BranchData
{
   GCG_BRANCH_TYPE       branchtype;         /**< type of branching */
   SCIP_CONS*            mastercons;         /**< constraint enforcing the branching restriction in the master problem */
   SCIP_Real             constant;           /**< constant value of the branching constraint in the master problem - either lhs or rhs, depending on branchtype */
   GCG_COMPBND*          B;                  /**< component bound sequence which induce the current branching constraint */
   int                   Bsize;              /**< size of the component bound sequence B */
   int                   consblocknr;        /**< id of the pricing problem (or block) to which this branching constraint belongs */
   int                   nvars;              /**< number of master variables the last time the node has been visited - neccessary to later include newly generated master variables */
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
   (*nodebranchdata)->consblocknr = blocknr;
   (*nodebranchdata)->mastercons = NULL;
   (*nodebranchdata)->constant = constant;
   (*nodebranchdata)->B = B;
   (*nodebranchdata)->Bsize = Bsize;
   (*nodebranchdata)->nvars = 0;

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
      if ( (B[i].sense == GCG_COMPSENSE_GE && SCIPisLT(scip, generatorentry, B[i].bound)) ||
           (B[i].sense == GCG_COMPSENSE_LE && SCIPisGT(scip, generatorentry, B[i].bound)) )
      {
         return FALSE;
      }
   }

   return TRUE;
}

/** creates the constraint for branching directly on a master variable */
static
SCIP_RETCODE createBranchingCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to add constraint */
   GCG_BRANCHDATA*       branchdata          /**< branching data structure */
)
{
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(node != NULL);
   assert(branchdata != NULL);

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "child(%d, %g)", branchdata->Bsize, branchdata->constant);

   assert(branchdata->mastercons == NULL);

   /*  create constraint for child */
   if (branchdata->branchtype == GCG_BRANCH_DOWN)
   {
      SCIP_CALL( SCIPcreateConsLinear(scip, &(branchdata->mastercons), name, 0, NULL, NULL,
         -SCIPinfinity(scip), branchdata->constant, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
   }
   else
   {
      SCIP_CALL( SCIPcreateConsLinear(scip, &(branchdata->mastercons), name, 0, NULL, NULL,
         branchdata->constant, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
   }

   SCIP_CALL( SCIPaddConsNode(scip, node, branchdata->mastercons, NULL) );

   return SCIP_OKAY;
}

 /** for given component bound sequence S, create |S|+1 Vanderbeck branching nodes */
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
   int nbranchcands;
   SCIP_VAR** mastervars;
   SCIP_VAR** branchcands;
   GCG_BRANCHDATA* downBranchData;
   GCG_BRANCHDATA* upBranchData;

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

   SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL, NULL) );

   /* determine the constant value of the master constraint, i.e. the rhs for the down branch, and lhs for the up branch */
   constantSum = 0;
   for( i = 0; i < nmastervars; ++i )
   {
      if( GCGisMasterVarInBlock(mastervars[i], blocknr) && isMasterVarInB(scip, mastervars[i], B, Bsize) )
      {
         constantSum += SCIPgetSolVal(masterscip, NULL, mastervars[i]);
      }
   }
   // sanity check: the sum must be fractional, otherwise something went wrong during seperation, i.e. B is incorrect
   assert(!SCIPisFeasIntegral(scip, constantSum));

   /* create two nodes */
   SCIPdebugMessage("Component bound branching rule: creating 2 nodes\n");
   SCIP_CALL( initNodeBranchdata(scip, &downBranchData, GCG_BRANCH_DOWN, SCIPfloor(scip, constantSum), B, Bsize, blocknr) );
   SCIP_CALL( initNodeBranchdata(scip, &upBranchData, GCG_BRANCH_UP, SCIPceil(scip, constantSum), B, Bsize, blocknr) );
   assert(downBranchData != NULL);
   assert(upBranchData != NULL);

   /* add the nodes to the tree */
   char downChildname[SCIP_MAXSTRLEN];
   char upChildname[SCIP_MAXSTRLEN];
   SCIP_NODE* downChild;
   SCIP_CONS* downChildcons;
   SCIP_NODE* upChild;
   SCIP_CONS* upChildcons;

   /* define names for origbranch constraints */
   (void) SCIPsnprintf(downChildname, SCIP_MAXSTRLEN, "node(%d) (last comp=%s %s %g) >= %g", blocknr,
      SCIPvarGetName(downBranchData->B[downBranchData->Bsize-1].component),
      downBranchData->B[downBranchData->Bsize-1].sense == GCG_COMPSENSE_GE? ">=": "<=",
      downBranchData->B[downBranchData->Bsize-1].bound,
      downBranchData->constant);
   (void) SCIPsnprintf(upChildname, SCIP_MAXSTRLEN, "node(%d) (last comp=%s %s %g) >= %g", blocknr,
      SCIPvarGetName(upBranchData->B[upBranchData->Bsize-1].component),
      upBranchData->B[upBranchData->Bsize-1].sense == GCG_COMPSENSE_GE? ">=": "<=",
      upBranchData->B[upBranchData->Bsize-1].bound,
      upBranchData->constant);

   SCIP_CALL( SCIPcreateChild(masterscip, &downChild, 0.0, SCIPgetLocalTransEstimate(masterscip)) );
   SCIP_CALL( SCIPcreateChild(masterscip, &upChild, 0.0, SCIPgetLocalTransEstimate(masterscip)) );

   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &downChildcons, downChildname, downChild,
      GCGconsMasterbranchGetActiveCons(masterscip), branchrule, downBranchData, NULL, 0, 0) );
   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &upChildcons, upChildname, upChild,
      GCGconsMasterbranchGetActiveCons(masterscip), branchrule, upBranchData, NULL, 0, 0) );

   SCIP_CALL( SCIPaddConsNode(masterscip, downChild, downChildcons, NULL) );
   SCIP_CALL( SCIPaddConsNode(masterscip, upChild, upChildcons, NULL) );

   SCIP_CALL( createBranchingCons(masterscip, downChild, downBranchData) );
   SCIP_CALL( createBranchingCons(masterscip, upChild, upBranchData) );

   /*  release constraints */
   SCIP_CALL( SCIPreleaseCons(masterscip, &downChildcons) );
   SCIP_CALL( SCIPreleaseCons(masterscip, &upChildcons) );


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

            // Make sure to not include dublicates in the set
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
               (*indexSet)[*indexSetSize-1] = origvars[j];
               ++(*indexSetSize);
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** recursive helper for the seperation algorithm determining a component bound sequence to branch on */
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
      for( i=0; i<Xsize; ++i)
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
   SCIP_Real value = SCIPfloor(masterscip, (min+max)/2);

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
   B1[*Bsize] = (GCG_COMPBND){origvar, GCG_COMPSENSE_LE, value};
   B2[*Bsize] = (GCG_COMPBND){origvar, GCG_COMPSENSE_GE, value+1};

   SCIPdebugMessage("B1 and B2 after adding the new component bounds\n");
   for (i = 0; i < new_Bsize; ++i)
   {
      SCIPdebugMessage("B1[%d]: %s %s %f\n", i, SCIPvarGetName(B1[i].component), B1[i].sense == GCG_COMPSENSE_LE ? "<=" : ">=", B1[i].bound);
      SCIPdebugMessage("B2[%d]: %s %s %f\n", i, SCIPvarGetName(B2[i].component), B2[i].sense == GCG_COMPSENSE_LE ? "<=" : ">=", B2[i].bound);
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
         assert(!inB1);
#endif
      }
   }
   SCIPdebugMessage("X1size: %d, X2size: %d\n", X1size, X2size);
   assert(X1size > 0 && X2size > 0);

   // determine the fractionality of B1 and B2
   SCIP_Real fractionality1 = calcFractionality(masterscip, X1, X1size);
   SCIP_Real fractionality2 = calcFractionality(masterscip, X2, X2size);
   // skipping assertions due to numerical issues
#ifdef SCIP_DISABLED_CODE
   assert(SCIPisLT(masterscip, fractionality1, fractionality / 2));
   assert(SCIPisLT(masterscip, fractionality2, fractionality / 2));
#endif
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

/** seperation algorithm determining a component bound sequence to branch on
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
      if( GCGisMasterVarInBlock(mastervars[i], blocknr) )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(origscip, &X, Xsize, Xsize + 1) );
         X[Xsize] = mastervars[i];
         Xsize++;
      }
   }

   /* 2. Call the recursive separation algorithm */
   SCIP_CALL( _separation(masterscip, X, Xsize, B, Bsize, result) );
   SCIPdebugMessage("Component bound sequence to branch on:\n");
   SCIPdebugMessage("Bsize: %d\n", *Bsize);
   for (int i = 0; i < *Bsize; ++i)
   {
      SCIPdebugMessage("B[%d]: %s %s %f\n", i, SCIPvarGetName((*B)[i].component), (*B)[i].sense == GCG_COMPSENSE_LE ? "<=" : ">=", (*B)[i].bound);
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

   /* a special case; branch on copy of an origvar directly */
#ifdef SCIP_DISABLED_CODE
   // TODO-TIL: Generic Branching handles this special case, but I have no idea why and for what purpose.
   if( blocknr == -1 )
   {
      assert(!GCGmasterVarIsLinking(mastervar));
      SCIPdebugMessage("branching on master variable\n");
      SCIP_CALL( branchDirectlyOnMastervar(origscip, mastervar, branchrule) );
      return SCIP_OKAY;
   }
#endif

   SCIPdebugMessage("branching in block %d \n", blocknr);

   /* 2. Call to seperation algorithm to find a suitable B to branch on in the current block.*/
   SCIP_CALL( separation(masterscip, &B, &Bsize, blocknr, result) );
   assert(Bsize > 0);
   assert(B != NULL);

   /* 3. Create the child nodes. */
   SCIP_CALL( createChildNodesCompBnd(origscip, branchrule, B, Bsize, blocknr, result) );


   SCIPdebugMessage("free B\n");
   SCIPfreeBlockMemoryArrayNull(origscip, &B, Bsize);

   return SCIP_OKAY;
}

/*
 * Callback methods of event handler
 *
 * The event handler is necessary to react on variable additions to the master problem.
 * If the new master variable satisfies all the component bound constraints, it must be added to the branching constraint.
 */

/* define not used callback as NULL*/
#define branchFreeCompBnd NULL
#define branchExitCompBnd NULL
#define branchInitsolCompBnd NULL
#define branchExitsolCompBnd NULL

/** adds a variable to a branching constraint */
static
SCIP_RETCODE addVarToMasterbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             mastervar,          /**< the variable to add */
   GCG_BRANCHDATA*       branchdata,         /**< branching data structure where the variable should be added */
   SCIP_Bool*            added               /**< whether the variable was added */
)
{
   assert(scip != NULL);
   assert(mastervar != NULL);
   assert(branchdata != NULL);
   assert(added != NULL);

   *added = FALSE;

   if( GCGvarGetBlock(mastervar) == -1 || branchdata->consblocknr == -3 || !GCGisMasterVarInBlock(mastervar, branchdata->consblocknr) )
      return SCIP_OKAY;

   if( isMasterVarInB(scip, mastervar, branchdata->B, branchdata->Bsize) )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, branchdata->mastercons, mastervar, 1.0) );
      *added = TRUE;
   }

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolCompBndbranchvaradd)
{  /*lint --e{715}*/
    assert(scip != NULL);
    assert(eventhdlr != NULL);
    assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

    /* notify SCIP that your event handler wants to react on the event type */
    SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, NULL) );

    return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolCompBndbranchvaradd)
{  /*lint --e{715}*/
    assert(scip != NULL);
    assert(eventhdlr != NULL);
    assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

    /* notify SCIP that your event handler wants to drop the event type */
    SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, -1) );

    return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecCompBndbranchvaradd)
{  /*lint --e{715}*/
   SCIP* origscip;
   SCIP_CONS* masterbranchcons;
   SCIP_CONS* parentcons;
   SCIP_VAR* mastervar;
   SCIP_VAR** allorigvars;
   SCIP_VAR** mastervars;
   GCG_BRANCHDATA* branchdata;
   int allnorigvars;
   int nmastervars;

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(scip != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

   mastervar = SCIPeventGetVar(event);
   if( !GCGvarIsMaster(mastervar) )
      return SCIP_OKAY;

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   masterbranchcons = GCGconsMasterbranchGetActiveCons(scip);
   assert(masterbranchcons != NULL);

   /* if branch rule is not component bound, abort */
   if( !GCGisBranchruleCompBnd(GCGconsMasterbranchGetBranchrule(masterbranchcons)) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetVarsData(origscip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   parentcons = masterbranchcons;
   branchdata = GCGconsMasterbranchGetBranchdata(parentcons);


   if( GCGvarIsMaster(mastervar) && GCGconsMasterbranchGetBranchrule(parentcons) != NULL )
   {
      SCIP_Bool added = FALSE;
      while( parentcons != NULL && branchdata != NULL
            && branchdata->B != NULL && branchdata->Bsize > 0 )
      {
         if( GCGconsMasterbranchGetBranchrule(parentcons) == NULL || strcmp(SCIPbranchruleGetName(GCGconsMasterbranchGetBranchrule(parentcons)), "generic") != 0 )
            break;

         assert(branchdata != NULL);


         if( (branchdata->consblocknr != GCGvarGetBlock(mastervar) && GCGvarGetBlock(mastervar) != -1 )
               || (GCGvarGetBlock(mastervar) == -1 && !GCGmasterVarIsLinking(mastervar)) )
         {
            parentcons = GCGconsMasterbranchGetParentcons(parentcons);

            if( parentcons != NULL )
               branchdata = GCGconsMasterbranchGetBranchdata(parentcons);

            continue;
         }

         SCIP_CALL( addVarToMasterbranch(scip, mastervar, branchdata, &added) );

         parentcons = GCGconsMasterbranchGetParentcons(parentcons);
         branchdata = GCGconsMasterbranchGetBranchdata(parentcons);
      }
   }

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
   SCIP_VAR** mastervars;
   int nmastervars;
   int i;
   int nvarsadded;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   SCIPdebugMessage("branchActiveMasterCompBnd: Block %d, Ssize %d\n", branchdata->consblocknr, branchdata->Bsize);

   if( branchdata->nvars >= SCIPgetNVars(scip) )
      return SCIP_OKAY;

   nmastervars = SCIPgetNVars(scip);
   mastervars = SCIPgetVars(scip);
   nvarsadded = 0;

   for( i = branchdata->nvars; i < nmastervars; ++i )
   {
      SCIP_Bool added = FALSE;
      assert(mastervars[i] != NULL);
      assert(GCGvarIsMaster(mastervars[i]));

      SCIP_CALL( addVarToMasterbranch(scip, mastervars[i], branchdata, &added) );
      if( added )
         ++nvarsadded;

   }
   SCIPdebugMessage("%d/%d vars added with contant=%g\n", nvarsadded, nmastervars-branchdata->nvars, branchdata->constant);
   branchdata->nvars = nmastervars;

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

   SCIPdebugMessage("branchDeactiveMasterCompBnd: Block %d, Ssize %d\n", branchdata->consblocknr, branchdata->Bsize);

   /* set number of variables since last call */
   branchdata->nvars = SCIPgetNVars(scip);
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

   if( (*branchdata)->mastercons != NULL )
   {
      SCIPdebugMessage("branchDataDeleteCompBnd: child blocknr %d, %s\n", (*branchdata)->consblocknr,
         SCIPconsGetName((*branchdata)->mastercons) );
   }
   else
   {
      SCIPdebugMessage("branchDataDeleteCompBnd: child blocknr %d, empty mastercons\n", (*branchdata)->consblocknr);
   }

   /* release constraint that enforces the branching decision */
   if( (*branchdata)->mastercons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(GCGgetMasterprob(scip), &(*branchdata)->mastercons) );
      (*branchdata)->mastercons = NULL;
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


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitCompBnd)
{
   SCIP* origscip;

   origscip = GCGmasterGetOrigprob(scip);
   assert(branchrule != NULL);
   assert(origscip != NULL);

   SCIPdebugMessage("Init method of Vanderbecks generic branching\n");

   SCIP_CALL( GCGrelaxIncludeBranchrule(origscip, branchrule, branchActiveMasterCompBnd,
         branchDeactiveMasterCompBnd, branchPropMasterCompBnd, branchMasterSolvedCompBnd, branchDataDeleteCompBnd) );

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

   /* include event handler for adding generated mastervars to the branching constraints */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, eventInitsolCompBndbranchvaradd, eventExitsolCompBndbranchvaradd,
         NULL, eventExecCompBndbranchvaradd,
         NULL) );

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
