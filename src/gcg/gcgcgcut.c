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
struct GCG_MastersepacutData
{
   int*                    conssindices;           /**< indices of constraints used to create cut */
   SCIP_Real*              weights;                /**< weights used to create cut */
   int                     nconssindices;          /**< number of constraints used to create cut */

};

/**< creates a Chvatal-Gomory cut */
SCIP_RETCODE GCGcreateChvatalGomoryCut(
   GCG*                       gcg,                   /**< GCG data structure */
   GCG_MASTERSEPACUT**   mastersepacut,         /**< pointer to store separator master cut */
   GCG_SEPA*                  sepa,                  /**< separator creating this cut */
   GCG_VARHISTORY*            varhistory,            /**< variables history of Chvatal-Gomory cut */
   SCIP_Real*                 weights,               /**< weights which were used to create the cut */
   int*                       indices,               /**< indices of constraints used to create the cut */
   int                        n                      /**< number of constraints used to create the cut */
   )
{
   SCIP* masterscip;
   GCG_MASTERSEPACUTDATA* data;

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

   SCIP_CALL( GCGcreateMastersepacut(gcg, mastersepacut, sepa, varhistory, data) );
   return SCIP_OKAY;
}

/** frees data of Chvatal-Gomory cut */
SCIP_RETCODE GCGfreeChvatalGomoryCutData(
   GCG*                             gcg,           /**< GCG data structure */
   GCG_MASTERSEPACUTDATA**     data           /**< pointer to data of subset row cut */
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
   GCG_MASTERSEPACUT*      mastersepacut     /**< separator master cut */
   )
{
   GCG_MASTERSEPACUTDATA* data;

   assert(mastersepacut != NULL);

   data = GCGmastersepacutGetData(mastersepacut);
   assert(data != NULL);

   return data->nconssindices;
}

/**< returns the weights of Chvatal-Gomory cut */
SCIP_Real* GCGchvatalGomoryCutGetWeights(
   GCG_MASTERSEPACUT*      mastersepacut     /**< separator master cut */
   )
{
   GCG_MASTERSEPACUTDATA* data;

   assert(mastersepacut != NULL);

   data = GCGmastersepacutGetData(mastersepacut);
   assert(data != NULL);

   return data->weights;
}

/**< returns the constraint indices of Chvatal-Gomory cut */
int* GCGchvatalGomoryCutGetConssIndices(
   GCG_MASTERSEPACUT*      mastersepacut     /**< separator master cut */
   )
{
   GCG_MASTERSEPACUTDATA* data;

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
   GCG_MASTERSEPACUT* cut;

   assert(gcg != NULL);
   assert(gcgcol != NULL);
   assert(GCGcolGetInitializedCoefs(gcgcol));

   scip = GCGgetMasterprob(gcg);
   cut = GCGextendedmasterconsGetMastersepacut(mastercutdata);
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
