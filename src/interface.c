/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    interface.c
 * @brief   scip interface for gcg
 * @author  Martin Bergner
 */
#include "interface.h"
#include "relax_gcg.h"

/**
 * Set the blocks for each variable and the master constraints without the need
 * to write blk files.
 *
 * @param scip Pointer to the scip instance
 * @param blocksPerVar an array sorted like the SCIP var array, partition for each variable
 * @param nblocks The number of partitions
 * @param masterConstraints the id of the masterconstraints in the constraint array
 * @param nMasterConstraints The number of masterconstraints
 * @return
 */
SCIP_RETCODE GCGsetBlocksForProblem(
   SCIP *scip, 
   int* blocksPerVar, 
   int nblocks, 
   int* masterConstraints, 
   int nMasterConstraints
   ) 
{
   int nvars;
   SCIP_VAR** vars;
   int i;
   int nconss;
   SCIP_CONS** conss;

   assert(scip!= 0);
   assert(blocksPerVar != 0);

   GCGrelaxSetNPricingprobs(scip, nblocks);
   SCIP_CALL( GCGrelaxCreateOrigVarsData(scip) );
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);
   for(i = 0; i < nvars; i++)
   {
      assert(vars[i] != 0);
      assert(blocksPerVar[i] >= 0 );
      assert(blocksPerVar[i] <= nblocks);
      if(blocksPerVar[i] >= nblocks) {
         continue;
      }
      assert(blocksPerVar[i] < nblocks);
      SCIP_CALL(GCGrelaxSetOriginalVarBlockNr(vars[i], blocksPerVar[i]));
   }

   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   for(i = 0; i < nMasterConstraints; i++)
   {
      assert(masterConstraints[i] >= 0);
      assert(masterConstraints[i] < nconss);
      SCIP_CALL(GCGrelaxMarkConsMaster(scip, conss[masterConstraints[i]]));
   }

   return SCIP_OKAY;
}
