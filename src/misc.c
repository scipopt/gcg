/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* #define SCIP_DEBUG */
/* #define CHECKCONSISTENCY */
/**@file    misc.c
 * @brief   miscellaneous methods
 * @author  Gerald Gamrath
 * @author  Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "relax_gcg.h"
#include "pub_gcgvar.h"

/** transforms given solution of the master problem into solution of the original problem
 *  TODO: think about types of epsilons used in this method*/
SCIP_RETCODE GCGrelaxTransformMastersolToOrigsol(
   SCIP*                 scip,               /** SCIP data structure */
   SCIP_SOL*             mastersol,          /** solution of the master problem, or NULL for current LP solution */
   SCIP_SOL**            origsol             /** pointer to store the new created original problem's solution */
   )
{
   SCIP* masterprob;
   int npricingprobs;
   int* blocknrs;
   SCIP_Real* blockvalue;
   SCIP_Real increaseval;
   SCIP_VAR** mastervars;
   SCIP_Real* mastervals;
   int nmastervars;
   int i;
   int j;

   assert(scip != NULL);
   assert(origsol != NULL);

   masterprob = GCGrelaxGetMasterprob(scip);
   npricingprobs = GCGrelaxGetNPricingprobs(scip);

   assert( !SCIPisInfinity(scip, SCIPgetSolOrigObj(masterprob, mastersol)) );

   SCIP_CALL( SCIPcreateSol(scip, origsol, NULL) );

   SCIP_CALL( SCIPallocBufferArray(scip, &blockvalue, npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocknrs, npricingprobs) );

   /* get variables of the master problem and their solution values */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &mastervals, nmastervars) );
   SCIP_CALL( SCIPgetSolVals(masterprob, mastersol, nmastervars, mastervars, mastervals) );

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

      origvars = GCGmasterVarGetOrigvars(mastervars[i]);
      norigvars = GCGmasterVarGetNOrigvars(mastervars[i]);
      origvals = GCGmasterVarGetOrigvals(mastervars[i]);
      blocknr = GCGvarGetBlock(mastervars[i]);
      isray = GCGmasterVarIsRay(mastervars[i]);

      assert(GCGvarIsMaster(mastervars[i]));
      assert(!SCIPisFeasNegative(scip, mastervals[i]));

      /* TODO: handle infinite master solution values */
      assert(!SCIPisInfinity(scip, mastervals[i]));

      /* first of all, handle variables representing rays */
      if( isray )
      {
         assert(blocknr >= 0);
         /* we also want to take into account variables representing rays, that have a small value (between normal and feas eps),
          * so we do no feas comparison here */
         if( SCIPisPositive(scip, mastervals[i]) )
         {
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < norigvars; j++ )
            {
               if(SCIPisZero(scip, origvals[j]))
                  break;

               assert(!SCIPisZero(scip, origvals[j]));

               /* the original variable is a linking variable: just transfer the solution value of the direct copy (this is done later) */
               if( GCGvarIsLinking(origvars[j]) )
                  continue;

//               SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origvars[j]), origvals[j] * mastervals[i], SCIPvarGetName(mastervars[i]));
               /* increase the corresponding value */
               SCIP_CALL( SCIPincSolVal(scip, *origsol, origvars[j], origvals[j] * mastervals[i]) );
            }
         }
         mastervals[i] = 0.0;
         continue;
      }

      /* handle the variables with value >= 1 to get integral values in original solution */
      while( SCIPisFeasGE(scip, mastervals[i], 1.0) )
      {
         /* variable was directly transferred to the master problem (only in linking conss or linking variable) */
         /* TODO: this may be the wrong place for this case, handle it before the while loop
          * and remove the similar case in the next while loop */
         if( blocknr == -1 )
         {
            assert(norigvars == 1);
            assert(origvals[0] == 1.0);

            /* increase the corresponding value */
//            SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origvars[0]), origvals[0] * mastervals[i],  SCIPvarGetName(mastervars[i]));
            SCIP_CALL( SCIPincSolVal(scip, *origsol, origvars[0], origvals[0] * mastervals[i]) );
            mastervals[i] = 0.0;
         }
         else
         {
            assert(blocknr >= 0);
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < norigvars; j++ )
            {
               SCIP_VAR* pricingvar;
               int norigpricingvars;
               SCIP_VAR** origpricingvars;
               if(SCIPisZero(scip, origvals[j]))
                  break;
               assert(!SCIPisZero(scip, origvals[j]));

               /* the original variable is a linking variable: just transfer the solution value of the direct copy (this is done above) */
               if( GCGvarIsLinking(origvars[j]) )
                  continue;

               pricingvar = GCGoriginalVarGetPricingVar(origvars[j]);
               assert(GCGvarIsPricing(pricingvar));

               norigpricingvars = GCGpricingVarGetNOrigvars(pricingvar);
               origpricingvars = GCGpricingVarGetOrigvars(pricingvar);

               /* just in case a variable has a value higher than the number of blocks, it represents */
               if( norigpricingvars <= blocknrs[blocknr] )
               {
//                  SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[norigpricingvars-1]), mastervals[i] * origvals[j], SCIPvarGetName(mastervars[i]));
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, origpricingvars[norigpricingvars-1], mastervals[i] * origvals[j]) );
                  mastervals[i] = 1.0;
               }
               /* this should be default */
               else
               {
//                  SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[blocknrs[blocknr]]), origvals[j], SCIPvarGetName(mastervars[i]) );
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, origpricingvars[blocknrs[blocknr]], origvals[j]) );
               }
            }
            mastervals[i] = mastervals[i] - 1.0;
            blocknrs[blocknr]++;
         }
      }
   }

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

      if( SCIPisFeasZero(scip, mastervals[i]) )
      {
         continue;
      }
      assert(SCIPisFeasGE(scip, mastervals[i], 0.0) && SCIPisFeasLT(scip, mastervals[i], 1.0));

      while( SCIPisFeasPositive(scip, mastervals[i]) )
      {
         assert(GCGvarIsMaster(mastervars[i]));
         assert(!GCGmasterVarIsRay(mastervars[i]));

         if( blocknr == -1 )
         {
            assert(norigvars == 1);
            assert(origvals[0] == 1.0);

//            SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origvars[0]), origvals[0] * mastervals[i], SCIPvarGetName(mastervars[i]) );
            /* increase the corresponding value */
            SCIP_CALL( SCIPincSolVal(scip, *origsol, origvars[0], origvals[0] * mastervals[i]) );
            mastervals[i] = 0.0;
         }
         else
         {
            increaseval = MIN(mastervals[i], 1.0 - blockvalue[blocknr]);
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < norigvars; j++ )
            {
               SCIP_VAR* pricingvar;
               int norigpricingvars;
               SCIP_VAR** origpricingvars;

               if( SCIPisZero(scip, origvals[j]) )
                  continue;

               /* the original variable is a linking variable: just transfer the solution value of the direct copy (this is done above) */
               if( GCGvarIsLinking(origvars[j]) )
                  continue;

               pricingvar = GCGoriginalVarGetPricingVar(origvars[j]);
               assert(GCGvarIsPricing(pricingvar));

               norigpricingvars = GCGpricingVarGetNOrigvars(pricingvar);
               origpricingvars = GCGpricingVarGetOrigvars(pricingvar);

               if( norigpricingvars <= blocknrs[blocknr] )
               {
                  increaseval = mastervals[i];

//                  SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[norigpricingvars-1]), origvals[j] * increaseval, SCIPvarGetName(mastervars[i]) );
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, origpricingvars[norigpricingvars-1], origvals[j] * increaseval) );
               }
               else
               {
                  /* increase the corresponding value */
//                  SCIPdebugMessage("Increasing value of %s by %f because of %s\n", SCIPvarGetName(origpricingvars[blocknrs[blocknr]]), origvals[j] * increaseval, SCIPvarGetName(mastervars[i]) );
                  SCIP_CALL( SCIPincSolVal(scip, *origsol, origpricingvars[blocknrs[blocknr]], origvals[j] * increaseval) );
               }
            }

            mastervals[i] = mastervals[i] - increaseval;
            if( SCIPisFeasZero(scip, mastervals[i]) )
            {
               mastervals[i] = 0.0;
            }
            blockvalue[blocknr] += increaseval;

            /* if the value assigned to the block is equal to 1, this block is full and we take the next block */
            if( SCIPisFeasGE(scip, blockvalue[blocknr], 1.0) )
            {
               blockvalue[blocknr] = 0.0;
               blocknrs[blocknr]++;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &mastervals);
   SCIPfreeBufferArray(scip, &blocknrs);
   SCIPfreeBufferArray(scip, &blockvalue);

   return SCIP_OKAY;
}


/** transforms given values of the given original variables into values of the given master variables */
void GCGrelaxTransformOrigvalsToMastervals(
   SCIP*                 scip,               /** SCIP data structure */
   SCIP_VAR**            origvars,           /** array with (subset of the) original variables */
   SCIP_Real*            origvals,           /** array with values (coefs) for the given original variables */
   int                   norigvars,          /** number of given original variables */
   SCIP_VAR**            mastervars,         /** array of (all present) master variables */
   SCIP_Real*            mastervals,         /** array to store the values of the master variables */
   int                   nmastervars         /** number of master variables */
   )
{
   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(origvars != NULL);
   assert(origvals != NULL);
   assert(mastervars != NULL);
   assert(mastervals != NULL);
   assert(nmastervars >= 0);

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
       * hence, we transfer the value directly to the corresponding master variabe
       */
      if( blocknr < 0 )
      {
         assert(blocknr == -1 || blocknr == -2);
         for( k = 0; k < nmastervars; k++ )
         {
            assert(!SCIPvarIsTransformedOrigvar(mastervars[k]));
            if( mastervars[k] == varmastervars[0])
            {
               assert(!SCIPvarIsTransformedOrigvar(varmastervars[0]));
               mastervals[k] += (varmastervals[0] * origvals[i]);
               break;
            }
         }
         assert(k < nmastervars);
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
}

/**  prints the given variable: name, type (original, master or pricing) block number,
 * and the list of all variables related to the given variable
 */
void GCGrelaxPrintVar(
   SCIP_VAR*             var                 /**< variable that should be printed */
   )
{
   int i;
   int blocknr;
   assert(GCGvarIsOriginal(var) || GCGvarIsMaster(var) || GCGvarIsPricing(var));

   blocknr = GCGvarGetBlock(var);

   if( GCGvarIsOriginal(var) )
   {
      SCIP_VAR** mastervars;
      SCIP_Real* mastervals;
      int  nmastervars;

      if( GCGvarIsLinking(var) )
      {
         SCIP_VAR** pricingvars;
         int nblocks;
         int j;
         pricingvars = GCGlinkingVarGetPricingVars(var);
         nblocks = GCGlinkingVarGetNBlocks(var);
         printf("Variable %s (linking): %d block%s (", SCIPvarGetName(var), nblocks, nblocks == 1 ? "":"s" );
         /*lint --e{440}*/
         for( i = 0, j = 0; j < nblocks; ++i)
         {
            if( pricingvars[i] != NULL )
            {
               printf("%d ", i);
               ++j;
            }
         }
         printf(")\n");
      }
      else
      {
         printf("Variable %s (original): block %d\n", SCIPvarGetName(var), blocknr);
      }

      mastervars = GCGoriginalVarGetMastervars(var);
      mastervals = GCGoriginalVarGetMastervals(var);
      nmastervars = GCGoriginalVarGetNMastervars(var);
      printf("mastervars:");
      for( i = 0; i < nmastervars-1; i++ )
      {
         printf("%s (%g), ", SCIPvarGetName(mastervars[i]), mastervals[i]);
      }
      printf("%s (%g)\n", SCIPvarGetName(mastervars[nmastervars-1]), mastervals[nmastervars-1]);
   }
   else if( GCGvarIsPricing(var) )
   {
      SCIP_VAR** origvars;
      int  norigvars;

      origvars = GCGpricingVarGetOrigvars(var);
      norigvars = GCGpricingVarGetNOrigvars(var);

      printf("Variable %s (pricing): block %d\n", SCIPvarGetName(var), blocknr);
      printf("origvars:");
      for( i = 0; i < norigvars-1; i++ )
      {
         printf("%s, ", SCIPvarGetName(origvars[i]));
      }
      printf("%s\n", SCIPvarGetName(origvars[norigvars-1]));
   }
   else if( GCGvarIsMaster(var) )
   {
      SCIP_VAR** origvars;
      int  norigvars;
      SCIP_Real* origvals;

      origvars = GCGmasterVarGetOrigvars(var);
      norigvars = GCGmasterVarGetNOrigvars(var);
      origvals = GCGmasterVarGetOrigvals(var);
      printf("Variable %s (master): block %d\n", SCIPvarGetName(var), blocknr);
      printf("origvars:");
      for( i = 0; i < norigvars-1; i++ )
      {
         printf("%s (%g), ", SCIPvarGetName(origvars[i]), origvals[i]);
      }
      printf("%s (%g)\n", SCIPvarGetName(origvars[norigvars-1]), origvals[norigvars-1]);
   }
}
