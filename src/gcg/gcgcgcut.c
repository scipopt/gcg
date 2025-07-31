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

/**@file   gcgcgcut.c
 * @brief  methods for working with Chvatal-Gomory cuts
 * @author Chantal Reinartz Groba
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "gcg/gcg.h"
#include "gcg/pub_mastersepacut.h"
#include "gcg/struct_mastersepacut.h"
#include "gcg/struct_gcgcol.h"

/** additional data for subset row cuts */
struct GCG_SeparatorMasterCutData
{
   int*                    conssindices;           /**< indices of constraints used to create cut */
   SCIP_Real*              weights;                /**< weights used to create cut */
   int                     nconssindices;          /**< number of constraints used to create cut */

};

/**< creates a Chvatal-Gomory cut */
SCIP_RETCODE GCGcreateChvatalGomoryCut(
   GCG*                       gcg,                   /**< GCG data structure */
   GCG_SEPARATORMASTERCUT**   mastersepacut,         /**< pointer to store separator master cut */
   GCG_SEPA*                  sepa,                  /**< separator creating this cut */
   GCG_VARHISTORY*            varhistory,            /**< variables history of Chvatal-Gomory cut */
   SCIP_Real*                 weights,               /**< weights which were used to create the cut */
   int*                       indices,               /**< indices of constraints used to create the cut */
   int                        n                      /**< number of constraints used to create the cut */
   )
{
   SCIP* masterscip;
   GCG_SEPARATORMASTERCUTDATA* data;

   assert(gcg != NULL);
   assert(mastersepacut != NULL);
   assert(n >= 0);

   masterscip = GCGgetMasterprob(gcg);

   SCIP_CALL( SCIPallocBlockMemory(masterscip, &data) );
   data->nconssindices = n;
   data->weights = NULL;
   data->conssindices = NULL;

   if( n > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(masterscip, &(data->weights), weights, n) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(masterscip, &(data->conssindices), indices, n) );
   }

   SCIP_CALL( GCGcreateMasterSepaCut(gcg, mastersepacut, sepa, varhistory, data) );
   return SCIP_OKAY;
}

/** frees data of Chvatal-Gomory cut */
SCIP_RETCODE GCGfreeChvatalGomoryCutData(
   GCG*                             gcg,           /**< GCG data structure */
   GCG_SEPARATORMASTERCUTDATA**     data           /**< pointer to data of subset row cut */
   )
{
   SCIP* masterscip;
   if( *data == NULL )
      return SCIP_OKAY;

   masterscip = GCGgetMasterprob(gcg);
    
   if( (*data)->nconssindices > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(masterscip, &((*data)->weights), (*data)->nconssindices);
      SCIPfreeBlockMemoryArrayNull(masterscip, &((*data)->conssindices), (*data)->nconssindices);
   }

   assert((*data)->weights == NULL);
   assert((*data)->conssindices == NULL);

   SCIPfreeBlockMemory(masterscip, data);
   *data = NULL;

   return SCIP_OKAY;
}


/**< returns the number of weights of Chvatal-Gomory cut */
int GCGchvatalGomoryCutGetNWeights(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   )
{
   GCG_SEPARATORMASTERCUTDATA* data;

   assert(mastersepacut != NULL);

   data = GCGmastersepacutGetData(mastersepacut);
   assert(data != NULL);

   return data->nconssindices;
}

/**< returns the weights of Chvatal-Gomory cut */
SCIP_Real* GCGchvatalGomoryCutGetWeights(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   )
{
   GCG_SEPARATORMASTERCUTDATA* data;

   assert(mastersepacut != NULL);

   data = GCGmastersepacutGetData(mastersepacut);
   assert(data != NULL);

   return data->weights;
}

/**< returns the constraint indices of Chvatal-Gomory cut */
int* GCGchvatalGomoryCutGetConssIndices(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   )
{
   GCG_SEPARATORMASTERCUTDATA* data;

   assert(mastersepacut != NULL);

   data = GCGmastersepacutGetData(mastersepacut);
   assert(data != NULL);

   return data->conssindices;
}

/** computes the coefficient of a column for a Chvatal-Gomory cut */
SCIP_RETCODE GCGchvatalGomoryCutGetColumnCoefficient(
   GCG*                          gcg,        /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*   mastercutdata,    /**< master separator cut */
   GCG_COL*                      gcgcol,     /**< gcg column */
   SCIP_Real*                    coeff       /**< pointer to store the coefficient */
   )
{
   SCIP* scip;
   SCIP_Real*     weights;
   SCIP_Real*     mastercoeffs;
   int*           conssindices;
   int            n;
   int            i;
   GCG_SEPARATORMASTERCUT* cut;

   assert(gcg != NULL);
   assert(gcgcol != NULL);
   assert(GCGcolGetInitializedCoefs(gcgcol));

   scip = GCGgetMasterprob(gcg);
   cut = GCGextendedmasterconsGetSepamastercut(mastercutdata);
   mastercoeffs = GCGcolGetMastercoefs(gcgcol);
   weights = GCGchvatalGomoryCutGetWeights(cut);
   conssindices = GCGchvatalGomoryCutGetConssIndices(cut);
   n = GCGchvatalGomoryCutGetNWeights(cut);
   assert(mastercoeffs != NULL);
   assert(weights != NULL);
   assert(conssindices != NULL);

   /* use the coefficients of the master constraints to compute coefficient for cut */
   *coeff = 0.0;
   for( i = 0; i < n; i++ )
   {
      SCIPdebugMessage("w[%i]: %f, i[%i]: %i --> %f\n", i, weights[i], i, conssindices[i], mastercoeffs[conssindices[i]]);
      *coeff += weights[i] * mastercoeffs[conssindices[i]];
   }

   *coeff = SCIPfeasFloor(scip, *coeff);
   SCIPdebugMessage("column coefficient: %f\n", *coeff);
   return SCIP_OKAY;
}

/** computes the coefficient of a master variable for a Chvatal-Gomory cut */
SCIP_RETCODE GCGchvatalGomoryCutGetVariableCoefficient(
   GCG*                          gcg,              /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*   mastercutdata,    /**< master separator cut */
   SCIP_VAR**                    vars,             /**< pricing variables which define the master variable */
   SCIP_Real*                    vals,             /**< values of the pricing variables which define the master variables */
   int                           nvars,            /**< number of pricing variables which define the master variable */
   int                           probnr,           /**< index of the pricing problem which generated the master variable */
   SCIP_Real*                    coef              /**< pointer to store the coefficient */
   )
{
   SCIP*                      scip;
   SCIP*                      pricingscip;
   GCG_PRICINGMODIFICATION*   pricingmod;
   SCIP_CONS**                pricingconss;
   SCIP_VAR**                 pricingconsvars;
   SCIP_Bool                  success;
   SCIP_Bool                  found;
   SCIP_Real*                 pricingconscoeffs;
   int                        npricingconsvars;
   int                        i;
   int                        pos;

   assert(gcg != NULL);
   assert(mastercutdata != NULL);

   scip = GCGgetMasterprob(gcg);
   *coef = 0.0;
   pricingmod = GCGextendedmasterconsGetPricingModification(gcg, mastercutdata, probnr);

   /* no pricing modification for this problem: coefficient is zero */
   if( pricingmod == NULL )
   {
      SCIPdebugMessage("no pricing modification for pp%i --> variable coefficient 0\n", probnr);
      return SCIP_OKAY;
   }

   pricingscip = GCGgetPricingprob(gcg, probnr);
   pricingconss = GCGpricingmodificationGetAdditionalConss(pricingmod);

   /* get all the pricing variables and their coefficients in the constraint */
   SCIP_CALL( SCIPgetConsNVars(pricingscip, pricingconss[0], &npricingconsvars, &success) );
   assert(success);
   SCIPallocBufferArray(scip, &pricingconsvars, npricingconsvars);
   SCIPallocBufferArray(scip, &pricingconscoeffs, npricingconsvars);
   SCIP_CALL( SCIPgetConsVars(pricingscip, pricingconss[0], pricingconsvars, npricingconsvars, &success) );
   assert(success);
   SCIP_CALL( SCIPgetConsVals(pricingscip, pricingconss[0], pricingconscoeffs, npricingconsvars, &success) );
   assert(success);

   /* compute w^TAx using the pricing constraint */
   for( i = 0; i < npricingconsvars; i++ )
   {
      if( GCGvarIsInferredPricing(pricingconsvars[i]) )
         continue;

      found = SCIPsortedvecFindPtr((void**) vars, SCIPvarComp, (void*) pricingconsvars[i], nvars, &pos);

      if( !found )
         continue;

      *coef += pricingconscoeffs[i] * vals[pos];
   }

   /* finally, we round down w^TAx */
   SCIPdebugMessage("variable coefficient %f\n", *coef);
   *coef = SCIPfeasFloor(scip, *coef);

   /* clean-up */
   SCIPfreeBufferArray(scip, &pricingconsvars);
   SCIPfreeBufferArray(scip, &pricingconscoeffs);

   return SCIP_OKAY;
}

/** adapts the objectives of all the necessary pricing problems such that they consider the Chvatal-Gomory cut */
SCIP_RETCODE GCGchvatalGomorySetPricingObjectives(
   GCG*                          gcg,     /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*   cut,     /**< separator master cut */
   SCIP_Real                     dual     /**< the dual value of the separator master cut */
   )
{
   SCIP*                      pricingproblem;
   SCIP_VAR*                  coeffvar;
   GCG_PRICINGMODIFICATION**  pricingmodifications;
   int                        npricingmodifications;
   int                        pricingblocknr;
   int                        i;

   assert(gcg != NULL);
   assert(cut != NULL);

   /* get all the pricing modifications associated with this master cut */
   npricingmodifications = GCGextendedmasterconsGetNPricingModifications(cut);
   pricingmodifications = GCGextendedmasterconsGetPricingModifications(cut);

   /* set the objective value of each coefficient variable y to -dual of the cut it is associated with */
   for( i = 0; i < npricingmodifications; i++ )
   {
      pricingblocknr = GCGpricingmodificationGetBlock(pricingmodifications[i]);
      pricingproblem = GCGgetPricingprob(gcg, pricingblocknr);
      coeffvar = GCGpricingmodificationGetCoefVar(pricingmodifications[i]);

      if( dual >= 0.0 ) // @todo: theoretically, dual should always be non-positive: 'correct' it to zero
         SCIP_CALL( SCIPchgVarObj(pricingproblem, coeffvar, 0.0) );
      else
         SCIP_CALL( SCIPchgVarObj(pricingproblem, coeffvar, -1.0 * dual) );
   }

   return SCIP_OKAY;
}

/** adapts a GCG column such that it respects the pricing modification imposed by the Chvatal-Gomory cut */
SCIP_RETCODE GCGchvatalGomoryAdjustGCGColumn(
   GCG*                          gcg,        /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*   cut,        /**< separator master cut */
   GCG_COL*                      gcgcol      /**< gcg column */
   )
{
   GCG_PRICINGMODIFICATION* pricemod;

   assert(gcg != NULL);
   assert(cut != NULL);

   if( !GCGextendedmasterconsIsActive(cut) )
      return SCIP_OKAY;

   pricemod = GCGextendedmasterconsGetPricingModification(gcg, cut, GCGcolGetProbNr(gcgcol));
   if( pricemod != NULL )
   {
      SCIP_VAR* coefvar;
      SCIP_Real coefvarval;
      SCIP_Bool append = FALSE;
      SCIP_Bool insert = FALSE;
      int pos;

      coefvar = GCGpricingmodificationGetCoefVar(pricemod);
      assert(coefvar != NULL);

      /* we compute the value of y */
      if( GCGcolGetInitializedCoefs(gcgcol) )
         SCIP_CALL( GCGchvatalGomoryCutGetColumnCoefficient(gcg, cut, gcgcol,
                                                            &coefvarval) );
      else
         SCIP_CALL(
            GCGchvatalGomoryCutGetVariableCoefficient(gcg, cut,
               gcgcol->inferredpricingvars, gcgcol->inferredpricingvals, gcgcol->ninferredpricingvars, gcgcol->probnr,
               &coefvarval) );

      /* if the computed coefficient is zero, we do not need to modify the column */
      if( coefvarval == 0.0 )
         return SCIP_OKAY;

      /* 1. variable already in column: replace value (this indicates that this was not the violating constraint)
       * 2. variable not yet in column:
       *    a. variable can be appended and variable order (based on index) remains correct
       *    b. variable has to be inserted to maintain the correct order */
      if( SCIPvarCompare(gcgcol->inferredpricingvars[gcgcol->ninferredpricingvars - 1], coefvar) == -1 )
         append = TRUE; // variable can simply be appended (2.a)
      else
      {
         /* check if variable already in column */
         SCIP_VAR** vars;
         SCIP_Bool found;
         int nvars;

         vars = gcgcol->inferredpricingvars;
         nvars = gcgcol->ninferredpricingvars;

         /* search position variable should be at */
         found = SCIPsortedvecFindPtr((void**) vars, SCIPvarComp, (void*) coefvar, nvars, &pos);
         if( found )
            gcgcol->inferredpricingvals[pos] = coefvarval; // variable already in column (1)
         else
            insert = TRUE; // variable needs to be inserted in column at position pos (2.b)
      }

      /* cases 2.a and 2.b */
      if( !SCIPisZero(gcgcol->pricingprob, coefvarval) && (append || insert))
      {
         /* ensure column has enough space to include variable */
         if( gcgcol->maxvars < gcgcol->ninferredpricingvars + 1 )
         {
            int newmaxvars = SCIPcalcMemGrowSize(gcgcol->pricingprob, gcgcol->ninferredpricingvars + 1);
            SCIP_CALL( SCIPreallocBlockMemoryArray(gcgcol->pricingprob, &(gcgcol->inferredpricingvars), gcgcol->maxvars,
                                                   newmaxvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(gcgcol->pricingprob, &(gcgcol->inferredpricingvals), gcgcol->maxvars,
                                                   newmaxvars) );
            gcgcol->maxvars = newmaxvars;
         }

         if( append )
         {
            /* variable can simply be appended to array and order remains correct */
            gcgcol->inferredpricingvars[gcgcol->ninferredpricingvars] = coefvar;
            gcgcol->inferredpricingvals[gcgcol->ninferredpricingvars] = coefvarval;
            SCIPcaptureVar(gcgcol->pricingprob, gcgcol->inferredpricingvars[gcgcol->ninferredpricingvars]);
         }
         else if( insert )
         {
            int i = 0;

            /* we have to move all the variables (& and their values) stored behind pos */
            for( i = gcgcol->ninferredpricingvars; i > pos ; i-- )
            {
               gcgcol->inferredpricingvars[i] = gcgcol->inferredpricingvars[i - 1];
               gcgcol->inferredpricingvals[i] = gcgcol->inferredpricingvals[i - 1];
            }

            /* add variable at correct position */
            gcgcol->inferredpricingvars[pos] = coefvar;
            gcgcol->inferredpricingvals[pos] = coefvarval;
            SCIPcaptureVar(gcgcol->pricingprob, gcgcol->inferredpricingvars[pos]);
         }

         (gcgcol->ninferredpricingvars)++;
      }
   }

   return SCIP_OKAY;
}
