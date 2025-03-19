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

/**@file    misc.c
 * @brief   miscellaneous methods
 * @author  Gerald Gamrath
 * @author  Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/gcg.h"
#include "gcg/relax_gcg.h"
#include "gcg/pricer_gcg.h"
#include "gcg/benders_gcg.h"
#include "gcg/pub_gcgvar.h"
#include "gcg/cons_decomp.h"
#include "gcg/gcgsort.h"
#include "gcg/stat.h"
#include <string.h>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#else
#include <unistd.h>
#endif

/** comparefunction for lexicographical sort */
static
GCG_DECL_SORTPTRCOMP(mastervarcomp)
{
   SCIP* origprob = (SCIP *) userdata; /* TODO: continue here */
   SCIP_VAR* mastervar1;
   SCIP_VAR* mastervar2;
   SCIP_VAR** origvars1;
   SCIP_VAR** origvars2;
   SCIP_Real* origvals1;
   SCIP_Real* origvals2;
   int norigvars1;
   int norigvars2;
   int i;
   int j;

   mastervar1 = (SCIP_VAR*) elem1;
   mastervar2 = (SCIP_VAR*) elem2;

   assert(mastervar1 != NULL);
   assert(mastervar2 != NULL);

   if( GCGvarGetBlock(mastervar1) < 0 )
   {
      SCIPdebugMessage("linkingvar or directy transferred var\n");
   }
   if( GCGvarGetBlock(mastervar2) < 0 )
   {
      SCIPdebugMessage("linkingvar or directy transferred var\n");
   }

   origvars1 = GCGmasterVarGetOrigvars(mastervar1);
   norigvars1 = GCGmasterVarGetNOrigvars(mastervar1);
   origvals1 = GCGmasterVarGetOrigvals(mastervar1);
   origvars2 = GCGmasterVarGetOrigvars(mastervar2);
   norigvars2 = GCGmasterVarGetNOrigvars(mastervar2);
   origvals2 = GCGmasterVarGetOrigvals(mastervar2);

   for( i = 0, j = 0; i < norigvars1 || j < norigvars2; )
   {
      if( i < norigvars1 && SCIPvarGetType(origvars1[i]) > SCIP_VARTYPE_INTEGER )
         ++i;
      else if( j < norigvars2 && SCIPvarGetType(origvars2[j]) > SCIP_VARTYPE_INTEGER )
         ++j;
      else if( i < norigvars1 && (j >= norigvars2 || SCIPvarGetProbindex(origvars1[i]) < SCIPvarGetProbindex(origvars2[j])) )
      {
         if( SCIPisFeasGT(origprob, origvals1[i], 0.0) )
            return -1;
         ++i;
      }
      else if( j < norigvars2 && (i >= norigvars1 || SCIPvarGetProbindex(origvars1[i]) > SCIPvarGetProbindex(origvars2[j])) )
      {
         if( SCIPisFeasGT(origprob, origvals2[j], 0.0) )
            return 1;
         ++j;
      }
      else if( i < norigvars1 && j < norigvars2 )
      {
         if( SCIPisFeasGT(origprob, origvals1[i], origvals2[j]) )
            return -1;
         else if( SCIPisFeasGT(origprob, origvals2[j], origvals1[i]) )
            return 1;
         ++i;
         ++j;
      }

   }

   return 0;
}

/* transforms given solution of the master problem into solution of the original problem
 *  @todo think about types of epsilons used in this method
 *  @returns SCIP return code
 */
SCIP_RETCODE GCGtransformMastersolToOrigsol(
   GCG*                  gcg,                /* GCG data structure */
   SCIP_SOL*             mastersol,          /* solution of the master problem, or NULL for current LP solution */
   SCIP_SOL**            origsol,            /* pointer to store the new created original problem's solution */
   SCIP_Bool             ignorelocalvarbnds, /* check global or local varbounds */
   SCIP_Bool*            violatesvarbnds     /* pointer to variable to store whether the solution violates orig varbnds (can be NULL) */
   )
{
   SCIP* origprob;
   SCIP* masterprob;
   int npricingprobs;
   int* blocknrs;
   SCIP_Real* blockvalue;
   SCIP_Real increaseval;
   SCIP_VAR** mastervars;
   SCIP_VAR** mastervarsall;
   SCIP_Real* mastervals;
   int nmastervarsall;
   int nmastervars;
   SCIP_VAR** vars;
   int nvars;
   SCIP_Real feastol;
   SCIP_Bool discretization;
   int i;
   int j;

   assert(gcg != NULL);
   assert(origsol != NULL);

   origprob = GCGgetOrigprob(gcg);
   masterprob = GCGgetMasterprob(gcg);
   npricingprobs = GCGgetNPricingprobs(gcg);

   if( violatesvarbnds != NULL )
      *violatesvarbnds = FALSE;

   assert( !SCIPisInfinity(origprob, SCIPgetSolOrigObj(masterprob, mastersol)) );

   if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_BENDERS )
   {
      SCIP_SOL* relaxsol;

      relaxsol = GCGgetBendersRelaxationSol(gcg);

      SCIP_CALL( SCIPcreateSolCopy(origprob, origsol, relaxsol) );
      SCIP_CALL( SCIPunlinkSol(origprob, *origsol) );

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPcreateSol(origprob, origsol, GCGrelaxGetProbingheur(gcg)) );

   if( GCGgetDecompositionMode(gcg) != GCG_DECMODE_ORIGINAL && !GCGmasterIsSolValid(gcg, mastersol) )
   {
      if( violatesvarbnds != NULL )
         *violatesvarbnds = TRUE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(origprob, &blockvalue, npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(origprob, &blocknrs, npricingprobs) );

   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );

   if( discretization && (SCIPgetNContVars(origprob) > 0) )
   {
      SCIP_VAR** fixedvars;
      SCIP_Real* fixedvals;
      int nfixedvars;

      /* get variables of the master problem and their solution values */
      SCIP_CALL( SCIPgetVarsData(masterprob, &mastervarsall, &nmastervarsall, NULL, NULL, NULL, NULL) );

      /* getting the fixed variables */
      fixedvars = SCIPgetFixedVars(masterprob);
      nfixedvars = SCIPgetNFixedVars(masterprob);

      assert(mastervarsall != NULL || nmastervarsall == 0);
      assert(nmastervarsall >= 0);

      SCIP_CALL( SCIPallocBufferArray(origprob, &mastervars, nmastervarsall + nfixedvars) );

      SCIP_CALL( SCIPallocBufferArray(origprob, &mastervals, nmastervarsall + nfixedvars) );
      SCIP_CALL( SCIPallocBufferArray(origprob, &fixedvals, nfixedvars) );

      SCIP_CALL( SCIPgetSolVals(masterprob, mastersol, nmastervarsall, mastervarsall, mastervals) );
      SCIP_CALL( SCIPgetSolVals(masterprob, mastersol, nfixedvars, fixedvars, fixedvals) );

      nmastervars = 0;

      for( i = 0; i < nmastervarsall; ++i )
      {
         SCIP_Real solval;

         assert( i >= nmastervars );
         solval = mastervals[i];

         if( !SCIPisZero(origprob, solval) )
         {
            mastervars[nmastervars] = mastervarsall[i];
            mastervals[nmastervars] = solval;

            ++nmastervars;
         }
      }

      /* adding the fixed variables to the mastervars array */
      for( i = 0; i < nfixedvars; i++ )
      {
         SCIP_Real solval;

         solval = fixedvals[i];

         if( !SCIPisZero(origprob, solval) )
         {
            mastervars[nmastervars] = fixedvars[i];
            mastervals[nmastervars] = solval;

            ++nmastervars;
         }
      }

      SCIPfreeBufferArray(origprob, &fixedvals);

      /* sort mastervariables lexicographically */
      GCGsortPtrPtr((void**)mastervars, (void**) mastervals, mastervarcomp, origprob, nmastervars );
   }
   else
   {
      SCIP_VAR** fixedvars;
      SCIP_Real* fixedvals;
      int nfixedvars;

      /* get variables of the master problem and their solution values */
      SCIP_CALL( SCIPgetVarsData(masterprob, &mastervarsall, &nmastervarsall, NULL, NULL, NULL, NULL) );

      /* getting the fixed variables */
      fixedvars = SCIPgetFixedVars(masterprob);
      nfixedvars = SCIPgetNFixedVars(masterprob);

      assert(mastervarsall != NULL);
      assert(nmastervarsall >= 0);

      SCIP_CALL( SCIPallocBufferArray(origprob, &mastervars, nmastervarsall + nfixedvars) );
      SCIP_CALL( SCIPallocBufferArray(origprob, &mastervals, nmastervarsall + nfixedvars) );
      SCIP_CALL( SCIPallocBufferArray(origprob, &fixedvals, nfixedvars) );

      SCIP_CALL( SCIPgetSolVals(masterprob, mastersol, nmastervarsall, mastervarsall, mastervals) );
      SCIP_CALL( SCIPgetSolVals(masterprob, mastersol, nfixedvars, fixedvars, fixedvals) );

      nmastervars = 0;

      for( i = 0; i < nmastervarsall; ++i )
      {
         SCIP_Real solval;

         assert( i >= nmastervars );
         solval = mastervals[i];

         if( !SCIPisZero(origprob, solval) )
         {
            mastervars[nmastervars] = mastervarsall[i];
            mastervals[nmastervars] = solval;

            ++nmastervars;
         }
      }

      /* adding the fixed variables to the mastervars array */
      for( i = 0; i < nfixedvars; i++ )
      {
         SCIP_Real solval;

         solval = fixedvals[i];

         if( !SCIPisZero(origprob, solval) )
         {
            mastervars[nmastervars] = fixedvars[i];
            mastervals[nmastervars] = solval;

            ++nmastervars;
         }
      }

      SCIPfreeBufferArray(origprob, &fixedvals);
   }

   /* initialize the block values for the pricing problems */
   for( i = 0; i < npricingprobs; i++ )
   {
      blockvalue[i] = 0.0;
      blocknrs[i] = 0;
   }

   /* loop over all given master variables */
   for( i = 0; i < nmastervars; i++ )
   {
      SCIP_VAR** origvars;
      int norigvars;
      SCIP_Real* origvals;
      SCIP_Bool isray;
      int blocknr;

      if( SCIPisZero(masterprob, mastervals[i]) )
         continue;

      origvars = GCGmasterVarGetOrigvars(mastervars[i]);
      norigvars = GCGmasterVarGetNOrigvars(mastervars[i]);
      origvals = GCGmasterVarGetOrigvals(mastervars[i]);
      blocknr = GCGvarGetBlock(mastervars[i]);
      isray = GCGmasterVarIsRay(mastervars[i]);

      assert(GCGvarIsMaster(mastervars[i]));

      /** @todo handle infinite master solution values */
      assert(!SCIPisInfinity(origprob, mastervals[i]));

      /* first of all, handle variables representing rays */
      if( isray )
      {
         assert(blocknr >= 0);
         /* we also want to take into account variables representing rays, that have a small value (between normal and feas eps),
          * so we do no feas comparison here */
         if( SCIPisPositive(masterprob, mastervals[i]) )
         {
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < norigvars; j++ )
            {
               if( SCIPisZero(origprob, origvals[j]) )
                  continue;

               assert(!SCIPisZero(origprob, origvals[j]));

               /* the original variable is a linking variable: just transfer the solution value of the direct copy (this is done later) */
               if( GCGoriginalVarIsLinking(origvars[j]) )
                  continue;

               SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origvars[j]), origvals[j] * mastervals[i], SCIPvarGetName(mastervars[i]));
               /* increase the corresponding value */
               SCIP_CALL( SCIPincSolVal(origprob, *origsol, origvars[j], origvals[j] * mastervals[i]) );
            }
         }
         mastervals[i] = 0.0;
         continue;
      }

      /* variable was directly transferred to the master problem (only in linking conss or linking variable) */
      /** @todo this may be the wrong place for this case, handle it before the while loop
       * and remove the similar case in the next while loop */
      if( blocknr == -1 )
      {
         assert(norigvars == 1);
         assert(origvals[0] == 1.0);

         /* increase the corresponding value */
         SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origvars[0]), origvals[0] * mastervals[i],  SCIPvarGetName(mastervars[i]));
         SCIP_CALL( SCIPincSolVal(origprob, *origsol, origvars[0], origvals[0] * mastervals[i]) );
         mastervals[i] = 0.0;
         continue;
      }
      if( blocknr == -2 )
      {
         assert(norigvars == 0);

         mastervals[i] = 0.0;
         continue;
      }

      /* handle the variables with value >= 1 to get integral values in original solution */
      while( SCIPisFeasGE(masterprob, mastervals[i], 1.0) )
      {
         assert(blocknr >= 0);
         /* loop over all original variables contained in the current master variable */
         for( j = 0; j < norigvars; j++ )
         {
            SCIP_VAR* pricingvar;
            int norigpricingvars;
            SCIP_VAR** origpricingvars;
            if( SCIPisZero(origprob, origvals[j]) )
               continue;
            assert(!SCIPisZero(origprob, origvals[j]));

            /* the original variable is a linking variable: just transfer the solution value of the direct copy (this is done above) */
            if( GCGoriginalVarIsLinking(origvars[j]) )
               continue;

            pricingvar = GCGoriginalVarGetPricingVar(origvars[j]);
            assert(GCGvarIsPricing(pricingvar));

            norigpricingvars = GCGpricingVarGetNOrigvars(pricingvar);
            origpricingvars = GCGpricingVarGetOrigvars(pricingvar);

            /* just in case a variable has a value higher than the number of blocks, it represents */
            if( norigpricingvars <= blocknrs[blocknr] )
            {
               SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[norigpricingvars-1]), mastervals[i] * origvals[j], SCIPvarGetName(mastervars[i]));
               /* increase the corresponding value */
               SCIP_CALL( SCIPincSolVal(origprob, *origsol, origpricingvars[norigpricingvars-1], mastervals[i] * origvals[j]) );
               mastervals[i] = 1.0;
            }
            /* this should be default */
            else
            {
               SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[blocknrs[blocknr]]), origvals[j], SCIPvarGetName(mastervars[i]) );
               /* increase the corresponding value */
               SCIP_CALL( SCIPincSolVal(origprob, *origsol, origpricingvars[blocknrs[blocknr]], origvals[j]) );
            }
         }
         mastervals[i] = mastervals[i] - 1.0;
         blocknrs[blocknr]++;
      }
      assert(!SCIPisFeasNegative(masterprob, mastervals[i]));
   }

   /* TODO: Change order of mastervars */

   /* loop over all given master variables */
   for( i = 0; i < nmastervars; i++ )
   {
      SCIP_VAR** origvars;
      int norigvars;
      SCIP_Real* origvals;
      int blocknr;

      origvars = GCGmasterVarGetOrigvars(mastervars[i]);
      norigvars = GCGmasterVarGetNOrigvars(mastervars[i]);
      origvals = GCGmasterVarGetOrigvals(mastervars[i]);
      blocknr = GCGvarGetBlock(mastervars[i]);

      if( SCIPisFeasZero(masterprob, mastervals[i]) )
      {
         continue;
      }
      assert(SCIPisFeasGE(masterprob, mastervals[i], 0.0) && SCIPisFeasLT(masterprob, mastervals[i], 1.0));

      while( SCIPisFeasPositive(masterprob, mastervals[i]) )
      {
         assert(blocknr >= 0);
         assert(GCGvarIsMaster(mastervars[i]));
         assert(!GCGmasterVarIsRay(mastervars[i]));

         increaseval = MIN(mastervals[i], 1.0 - blockvalue[blocknr]);
         /* loop over all original variables contained in the current master variable */
         for( j = 0; j < norigvars; j++ )
         {
            SCIP_VAR* pricingvar;
            int norigpricingvars;
            SCIP_VAR** origpricingvars;

            if( SCIPisZero(origprob, origvals[j]) )
               continue;

            /* the original variable is a linking variable: just transfer the solution value of the direct copy (this is done above) */
            if( GCGoriginalVarIsLinking(origvars[j]) )
               continue;

            pricingvar = GCGoriginalVarGetPricingVar(origvars[j]);
            assert(GCGvarIsPricing(pricingvar));

            norigpricingvars = GCGpricingVarGetNOrigvars(pricingvar);
            origpricingvars = GCGpricingVarGetOrigvars(pricingvar);

            if( norigpricingvars <= blocknrs[blocknr] )
            {
               increaseval = mastervals[i];

               SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[norigpricingvars-1]), origvals[j] * increaseval, SCIPvarGetName(mastervars[i]) );
               /* increase the corresponding value */
               SCIP_CALL( SCIPincSolVal(origprob, *origsol, origpricingvars[norigpricingvars-1], origvals[j] * increaseval) );
            }
            else
            {
               /* increase the corresponding value */
               SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[blocknrs[blocknr]]), origvals[j] * increaseval, SCIPvarGetName(mastervars[i]) );
               SCIP_CALL( SCIPincSolVal(origprob, *origsol, origpricingvars[blocknrs[blocknr]], origvals[j] * increaseval) );
            }
         }

         mastervals[i] = mastervals[i] - increaseval;
         if( SCIPisFeasZero(masterprob, mastervals[i]) )
         {
            mastervals[i] = 0.0;
         }
         blockvalue[blocknr] += increaseval;

         /* if the value assigned to the block is equal to 1, this block is full and we take the next block */
         if( SCIPisFeasGE(masterprob, blockvalue[blocknr], 1.0) )
         {
            blockvalue[blocknr] = 0.0;
            blocknrs[blocknr]++;
         }
      }
   }

   SCIPfreeBufferArray(origprob, &mastervals);
   SCIPfreeBufferArray(origprob, &mastervars);

   SCIPfreeBufferArray(origprob, &blocknrs);
   SCIPfreeBufferArray(origprob, &blockvalue);

   /* if the solution violates one of its bounds by more than feastol
    * and less than 10*feastol, round it and print a warning
    */
   SCIP_CALL( SCIPgetVarsData(origprob, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetRealParam(origprob, "numerics/feastol", &feastol) );

   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real solval;
      SCIP_Real lb;
      SCIP_Real ub;

      solval = SCIPgetSolVal(origprob, *origsol, vars[i]);
      lb = ignorelocalvarbnds ? SCIPvarGetLbGlobal(vars[i]) : SCIPvarGetLbLocal(vars[i]);
      ub = ignorelocalvarbnds ? SCIPvarGetUbGlobal(vars[i]) : SCIPvarGetUbLocal(vars[i]);

      if( SCIPisFeasGT(origprob, solval, ub) && EPSEQ(solval, ub, 10 * feastol) )
      {
         SCIP_CALL( SCIPsetSolVal(origprob, *origsol, vars[i], ub) );
         SCIPwarningMessage(origprob, "Variable %s rounded from %g to %g in relaxation solution\n",
            SCIPvarGetName(vars[i]), solval, ub);
      }
      else if( SCIPisFeasLT(origprob, solval, lb) && EPSEQ(solval, lb, 10 * feastol) )
      {
         SCIP_CALL( SCIPsetSolVal(origprob, *origsol, vars[i], lb) );
         SCIPwarningMessage(origprob, "Variable %s rounded from %g to %g in relaxation solution\n",
            SCIPvarGetName(vars[i]), solval, lb);
      }
      if( !SCIPisFeasGE(origprob, SCIPgetSolVal(origprob, *origsol, vars[i]), lb) || !SCIPisFeasLE(origprob, SCIPgetSolVal(origprob, *origsol, vars[i]), ub) )
      {
         if( violatesvarbnds != NULL )
            *violatesvarbnds = TRUE;
      }
   }

   return SCIP_OKAY;
}


/* transforms given values of the given original variables into values of the given master variables
 * @returns the sum of the values of the corresponding master variables that are fixed */
SCIP_Real GCGtransformOrigvalsToMastervals(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_VAR**            origvars,           /**< array with (subset of the) original variables */
   SCIP_Real*            origvals,           /**< array with values (coefs) for the given original variables */
   int                   norigvars,          /**< number of given original variables */
   SCIP_VAR**            mastervars,         /**< array of (all present) master variables */
   SCIP_Real*            mastervals,         /**< array to store the values of the master variables */
   int                   nmastervars         /**< number of master variables */
   )
{
   SCIP* origprob;
   int i;
   int j;
   int k;
   SCIP_Real sum;

   assert(gcg != NULL);
   assert(origvars != NULL);
   assert(origvals != NULL);
   assert(mastervars != NULL);
   assert(mastervals != NULL);
   assert(nmastervars >= 0);

   sum = 0.0;

   origprob = GCGgetOrigprob(gcg);

   /* set all values to 0 initially */
   for( i = 0; i < nmastervars; i++ )
      mastervals[i] = 0.0;

   /* iterate over all original variables */
   for( i = 0; i < norigvars; i++ )
   {
      SCIP_VAR** varmastervars;
      SCIP_Real* varmastervals;
      int blocknr;

      assert(GCGvarIsOriginal(origvars[i]));
      varmastervars = GCGoriginalVarGetMastervars(origvars[i]);
      varmastervals = GCGoriginalVarGetMastervals(origvars[i]);
      blocknr = GCGvarGetBlock(origvars[i]);

      /* variable belongs to no block (or is a linking variable), so it was transferred directly to the master problem,
       * hence, we transfer the value directly to the corresponding master variable
       */
      if( blocknr < 0 )
      {
         assert(blocknr == -1 || blocknr == -2);
         assert(SCIPvarIsOriginal(varmastervars[0]));
         assert(SCIPvarGetTransVar(varmastervars[0]) != NULL);

         for( k = 0; k < nmastervars; k++ )
         {
            if( mastervars[k] == SCIPvarGetTransVar(varmastervars[0]) )
            {
               mastervals[k] += (varmastervals[0] * origvals[i]);
               break;
            }
         }

         if ( k >= nmastervars )
         {
            // inactive variable, check whether it is fixed
            if ( SCIPisFeasEQ(origprob, SCIPvarGetLbGlobal(varmastervars[0]), SCIPvarGetUbGlobal(varmastervars[0])) )
            {
               // variable is fixed in the master problem
               sum += (SCIPvarGetLbGlobal(varmastervars[0]) * varmastervals[0] * origvals[i]);
            }
            else
            {
#ifdef SCIP_DEBUG
               SCIP* masterprob = GCGgetMasterprob(gcg);
               SCIP_VAR** vars = SCIPgetVars(masterprob);

               SCIPdebugMessage("OrigVar %s [%f,%f]\n", SCIPvarGetName(origvars[i]), SCIPvarGetLbGlobal(origvars[i]),
                     SCIPvarGetUbGlobal(origvars[i]));
               SCIPdebugMessage("MasterVar %s [%f,%f]\n", SCIPvarGetName(varmastervars[0]),
                     SCIPvarGetLbGlobal(varmastervars[0]), SCIPvarGetUbGlobal(varmastervars[0]));
#endif
               assert(FALSE);
            }
         }
      }
      /* variable belongs to exactly one block, so we have to look at all master variables and increase their values
       * if they contain the original variable
       */
      else
      {
         SCIP_VAR* pricingvar;
         SCIP_VAR* origvar;
         SCIP_VAR** curmastervars;
         SCIP_Real* curmastervals;
         int ncurmastervars;

         pricingvar = GCGoriginalVarGetPricingVar(origvars[i]);
         assert(GCGvarIsPricing(pricingvar));

         origvar = GCGpricingVarGetOriginalVar(pricingvar);
         assert(GCGvarIsOriginal(origvar));
         curmastervars = GCGoriginalVarGetMastervars(origvar);
         curmastervals = GCGoriginalVarGetMastervals(origvar);
         ncurmastervars = GCGoriginalVarGetNMastervars(origvar);

         for( j = 0; j < ncurmastervars; j++ )
         {
            assert(SCIPvarIsTransformed(curmastervars[j]));
            for( k = 0; k < nmastervars; k++ )
               if( mastervars[k] == curmastervars[j] )
               {
                  mastervals[k] += (curmastervals[j] * origvals[i]);
                  break;
               }
            assert(k < nmastervars);
         }
      }

   }
   return sum;
}

/* checks whether the scip is the original scip instance
 * @returns whether the scip is the original scip instance */
SCIP_Bool GCGisOriginal(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   assert(scip != NULL);
   return SCIPfindRelax(scip, "gcg") != NULL;
}

/* checks whether the scip is the master problem scip
 * @returns whether the scip is the master problem scip */
SCIP_Bool GCGisMaster(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   assert(scip != NULL);
   return SCIPfindPricer(scip, "gcg") != NULL || SCIPfindBenders(scip, "gcg") != NULL;
}

/* print out GCG statistics
 * @returns SCIP return code */
SCIP_RETCODE GCGprintStatistics(
   GCG*                  gcg,                /* GCG data structure */
   FILE*                 file                /* output file or NULL for standard output */
)
{
   SCIP* origprob;
   SCIP* masterprob;
   assert(gcg != NULL);

   origprob = GCGgetOrigprob(gcg);
   masterprob = GCGgetMasterprob(gcg);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterprob), file, "\nMaster Program statistics:\n");
   SCIP_CALL( SCIPprintStatistics(masterprob, file) );
   if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_DANTZIGWOLFE
      && SCIPgetStage(masterprob) > SCIP_STAGE_PRESOLVED )
   {
      GCGpricerPrintPricingStatistics(gcg, file);
      SCIP_CALL( GCGwriteSolvingDetails(gcg) );
   }

   if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_DANTZIGWOLFE )
   {
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(origprob), file, "\nOriginal Program statistics:\n");
      SCIP_CALL( SCIPprintStatistics(origprob, file) );
   }
   else
   {
      assert(GCGgetDecompositionMode(gcg) == GCG_DECMODE_BENDERS
         || GCGgetDecompositionMode(gcg) == GCG_DECMODE_ORIGINAL);
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(origprob), file, "\nOriginal Program Solution statistics:\n");
      SCIPprintSolutionStatistics(origprob, file);
   }
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterprob), file, "\n");
   if( GCGgetDecompositionMode(gcg) == GCG_DECMODE_DANTZIGWOLFE
      && SCIPgetStage(origprob) >= SCIP_STAGE_SOLVING )
   {
      SCIP_CALL( GCGmasterPrintSimplexIters(gcg, file) );
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterprob), file, "\n");
   }
   if( GCGgetDecompositionMode(gcg) != GCG_DECMODE_ORIGINAL )
   {
      SCIP_CALL( GCGconshdlrDecompPrintDetectorStatistics(gcg, file) );
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterprob), file, "\n");
   }
   if( SCIPgetStage(origprob) >= SCIP_STAGE_PRESOLVING && GCGgetNPricingprobs(gcg) > 0 )
   {
      GCG_DECOMP* decomp = GCGgetStructDecomp(gcg);
      if( decomp != NULL )
      {
         SCIP_CALL( GCGprintDecompStatistics(gcg, file, decomp) );
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterprob), file, "\n");
      }
   }
   if( SCIPgetStage(origprob) >= SCIP_STAGE_PRESOLVING && GCGgetDecompositionMode(gcg) != GCG_DECMODE_ORIGINAL )
   {
      SCIP_CALL( GCGconshdlrDecompPrintScoreStatistics(gcg, file) );
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterprob), file, "\n");
   }
   return SCIP_OKAY;
}

/* print name of current instance to given output
 * @returns SCIP return code */
SCIP_RETCODE GCGprintInstanceName(
   GCG*                  gcg,                /* GCG data structure */
   FILE*                 file                /* output file or NULL for standard output */
)
{
   char problemname[SCIP_MAXSTRLEN];
   char* outputname;
   SCIP* origprob = GCGgetOrigprob(gcg);
   (void) SCIPsnprintf(problemname, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(origprob));
   SCIPsplitFilename(problemname, NULL, &outputname, NULL, NULL);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(origprob), file, "filename: %s \n", outputname );
   return SCIP_OKAY;
}


/* print out complete detection statistics
 * @returns SCIP return code */
SCIP_RETCODE GCGprintCompleteDetectionStatistics(
   GCG*                  gcg,                /* GCG data structure */
   FILE*                 file                /* output file or NULL for standard output */
)
{
   SCIP* origprob;
   assert(gcg != NULL);

   origprob = GCGgetOrigprob(gcg);

   if( !GCGdetectionTookPlace(gcg, TRUE) && !GCGdetectionTookPlace(gcg, FALSE) )
   {
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(origprob), file, "\nDetection did not take place so far\n");
      return SCIP_OKAY;
   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(origprob), file, "\nStart writing complete detection information:\n");

   SCIP_CALL( GCGprintInstanceName(gcg, file) );

   GCGprintBlockcandidateInformation(gcg, file);

   GCGprintCompleteDetectionTime(gcg, file);

   GCGprintPartitionInformation(gcg, file);

   GCGprintDecompInformation(gcg, file);

   return SCIP_OKAY;
}



/* Checks whether the constraint belongs to GCG or not
 *  @returns whether the constraint belongs to GCG or not */
SCIP_Bool GCGisConsGCGCons(
   SCIP_CONS*            cons                /* constraint to check */
   )
{
   SCIP_CONSHDLR* conshdlr;
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   if( strcmp("origbranch", SCIPconshdlrGetName(conshdlr)) == 0 )
      return TRUE;
   else if( strcmp("masterbranch", SCIPconshdlrGetName(conshdlr)) == 0 )
      return TRUE;

   return FALSE;
}
