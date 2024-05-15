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

/**@file    sepa_xyz.c
 *
 * @brief   xyz separator for master problem (put your description here)
 * @author  Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>

#include "mastercutdata.h"
#include "gcg.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "sepa_subsetrow.h"
#include "struct_mastercutdata.h"
#include "struct_sepagcg.h"
#include "type_mastercutdata.h"
#include "type_sepagcg.h"
#include "event_sepacuts.h"
#include "scip/cons_linear.h"
#include "pub_gcgcol.h"

#define SEPA_NAME              "subsetrow"
#define SEPA_DESC              "subsetrow separator"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    100
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define STARTMAXCUTS                50
#define DEFAULT_RANDSEED       71

/*
 * Data structures
 */

/* TODO: fill in the necessary separator data */

/** separator data */
struct SCIP_SepaData
{
   int                     sepaidx;       /**< index of the separator in relaxdata->separators */
   int                     ngeneratedcut; /**< counts the number of cuts generated */
   SCIP_RANDNUMGEN*        randnumgen;    /**< random number generator (for strategy RANDOM) */
   SCIP_Bool               enable;        /**< is this separator enabled?*/
};


/*
 * Callback methods of separator
 */

/* TODO: Implement all necessary separator methods. The methods with an #if 0 ... #else #define ... are optional */


#define sepaCopySubsetrow NULL
#define sepaExecsolSubsetrow NULL

static
SCIP_DECL_SEPAEXIT(sepaExitSubsetrow)
{
   SCIP_SEPADATA* sepadata;
   SCIPdebugMessage("exit sgcg sepa subsetrow\n");


   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeSubsetrow)
{
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);
   assert(sepa != NULL);

   SCIPdebugMessage("free gcg sepa subsetrow sepa data\n");
   sepadata = SCIPsepaGetData(sepa);
   SCIPfreeRandom(scip, &(sepadata->randnumgen));
   SCIPfreeBlockMemory(scip, &sepadata);

   return SCIP_OKAY;
}


/** randomly selects n different constraints from the master problem */
static
SCIP_RETCODE selectRandomRows(
   SCIP_SEPADATA* sepadata,                   /**< sepa data */
   int            nmasterconss,               /**< number of constraints in the master problem */
   int*           selectedmasterconssidx,     /**< pointer to store the indices of the selected constraints */
   int            n                           /**< number of constraints to be selected */
   )
{
   int i;
   int j;

   assert(n > 0 && n < nmasterconss);

   /* randomly select n of indices out of [0, ..., nmasterconss - 1] */
   i = 0;
   while( i < n )
   {
      selectedmasterconssidx[i] = SCIPrandomGetInt(sepadata->randnumgen, 0, nmasterconss-1);

      /* if we get an element that we already had, we will draw again */
      for( j = 0; j < i; j++ )
      {
         if( selectedmasterconssidx[i] == selectedmasterconssidx[j] )
         {
            --i;
            break;
         }
      }
      ++i;
   }

   for( i = 0; i < n; i++ )
   {
      SCIPdebugMessage("select index %i\n", selectedmasterconssidx[i]);
      assert(0 <= selectedmasterconssidx[i] && selectedmasterconssidx[i] < nmasterconss);

   }

   return SCIP_OKAY;
}

/** create a new row for master problem and fill it with the variables (and their coefficients) stored in the map */
static
SCIP_RETCODE createSubsetRowCut(
   SCIP*          masterscip,
   SCIP_ROW**     ssrc,
   SCIP_HASHMAP*  mapmastervarxcoeff,
   SCIP_Real      rhs_ssrc,
   SCIP_SEPA*     sepa
   )
{
   int nssrcvars;
   int i;

   assert(GCGisMaster(masterscip));
   assert(sepa != NULL);

   SCIPdebugMessage("create subset row cut\n");
   /* create 'empty' subset row cut of form -inf <= ... <= rhs_ssrc
    * - local, removable, modifiable */
   SCIP_CALL( SCIPcreateEmptyRowSepa(masterscip, &(*ssrc), sepa, "subsetrow", -SCIPinfinity(masterscip), floor(rhs_ssrc), TRUE, TRUE, TRUE) );
   assert(ssrc != NULL);
   assert(*ssrc != NULL);

   /* fill the row with master variables and their coefficients */
   nssrcvars = SCIPhashmapGetNEntries(mapmastervarxcoeff);
   for( i = 0; i < nssrcvars; i++ )
   {
      SCIP_HASHMAPENTRY* entry;
      SCIP_VAR* mastervar;
      SCIP_Real ssrc_coeff;

      entry = SCIPhashmapGetEntry(mapmastervarxcoeff, i);
      if( entry == NULL )
      {
         continue;
      }
      mastervar = (SCIP_VAR*) SCIPhashmapEntryGetOrigin(entry);
      ssrc_coeff = SCIPhashmapEntryGetImageReal(entry);
      ssrc_coeff = floor(ssrc_coeff);
      if( !SCIPisZero(masterscip, ssrc_coeff) )
      {
         SCIPdebugMessage("Add master variable %s with coefficient %f to subset row cut\n", SCIPvarGetName(mastervar), ssrc_coeff);
         SCIP_CALL( SCIPaddVarToRow(masterscip, *ssrc, mastervar, ssrc_coeff) );
      }
   }
   return SCIP_OKAY;
}

/** computes the rhs (w^Tb) and the coefficient for each variable (w^Ta_p) in the cut (still non-rounded) */
static
SCIP_RETCODE computeSubsetRowCoefficientsAndRHS(
   SCIP*          masterscip,             /**< SCIP data structure of master problem */
   SCIP_CONS**    masterconss,            /**< constraints of the master problem */
   int*           selectedconssidx,       /**< indices of the selected constraints */
   int            nselectedconss,         /**< number of selected constraints */
   int            k,                      /**< defines the weights: 1/k or -1/k */
   SCIP_Real**    weights,                /**< pointer to store weights of the selected constraints */
   SCIP_HASHMAP*  mapmastervarxcoeff,     /**< map to store the computed subsetrow cut coefficient for each master variable */
   SCIP_Real*     rhs_ssrc                /**< pointer to store rhs of subsetrow cut */
   )
{
   SCIP_CONS* mastercons;
   SCIP_VAR** masterconsvars;
   SCIP_Real* masterconscoeffs;
   SCIP_Real  rhs_mastercons;
   SCIP_Real  lhs_mastercons;
   SCIP_Real  coeff_ssrc;
   SCIP_Bool  success;
   int        nmasterconsvars;
   int        i;
   int        j;

   assert(masterscip != NULL);
   assert(mapmastervarxcoeff != NULL);
   assert(k != 0);

   SCIPdebugMessage("compute the rhs and coefficients of the subsetrow cut\n");
   *rhs_ssrc = 0;
   for( i = 0; i < nselectedconss; i++ )
   {
      /* lhs <= ax <= rhs */
      mastercons = masterconss[selectedconssidx[i]];
      lhs_mastercons = SCIPconsGetLhs(masterscip, mastercons, &success);
      assert(success);
      rhs_mastercons = SCIPconsGetRhs(masterscip, mastercons, &success);
      assert(success);

      /* we want all constraints to have form ax <= b to simplify the summation
       * - if constraint has form ax <= b: weight is  1/k
       * - if constraint has form b <= ax: weight is -1/k (to transform it to -ax <= -b) */
      if( SCIPisInfinity(masterscip, -lhs_mastercons ) ) // lhs = -inf --> constraint has form: ax <= b
      {
         (*weights)[i] = 1.0 / k;
         (*rhs_ssrc) += (*weights)[i] * rhs_mastercons;
         SCIPdebugMessage("master constraint %s (ax <= %f) with weight %f\n", SCIPconsGetName(mastercons), rhs_mastercons, (*weights)[i]);
      }
      else if( SCIPisInfinity(masterscip, rhs_mastercons ) ) // lhs != -inf & rhs = inf --> constraint has form: b <= ax
      {
         (*weights)[i] = -1.0 / k;
         (*rhs_ssrc) += (*weights)[i] * lhs_mastercons;
         SCIPdebugMessage("master constraint %s (%f <= ax) with weight %f\n", SCIPconsGetName(mastercons), lhs_mastercons, (*weights)[i]);
      }
      else // lhs != -inf & rhs != inf --> constraint has form: b_l <= ax <= b_r
      {
         /* we use the right side: ax <= b_r */
         (*weights)[i] = 1.0 / k;
         (*rhs_ssrc) += (*weights)[i] * rhs_mastercons;
         SCIPdebugMessage("master constraint %s (%f <= ax <= %f) with weight\n", SCIPconsGetName(mastercons), lhs_mastercons, rhs_mastercons);
      }
      SCIPdebugMessage("rhs of subset row cut: %f\n", *rhs_ssrc);

      /* get all variables and their corresponding coefficients in the master constraint */
      masterconsvars = NULL;
      masterconscoeffs = NULL;
      SCIP_CALL( SCIPgetConsNVars(masterscip, mastercons, &nmasterconsvars, &success) );
      assert(success);
      if( nmasterconsvars == 0 )
      {
         SCIPdebugMessage("constraint has no variables\n");
         continue;
      }

      SCIPallocBufferArray(masterscip, &masterconsvars, nmasterconsvars);
      SCIPallocBufferArray(masterscip, &masterconscoeffs, nmasterconsvars);
      SCIP_CALL( SCIPgetConsVars(masterscip, mastercons, masterconsvars, nmasterconsvars, &success) );
      assert(success);
      SCIP_CALL( SCIPgetConsVals(masterscip, mastercons, masterconscoeffs, nmasterconsvars, &success) );
      assert(success);
      /* for each variable add its weighted coefficient in this constraint to its coefficient for the subsetrow cut */
      for( j = 0; j < nmasterconsvars; j++ )
      {
         coeff_ssrc = SCIPhashmapGetImageReal(mapmastervarxcoeff, masterconsvars[j]);
         if( coeff_ssrc == SCIP_INVALID )
         {
            SCIPhashmapSetImageReal(mapmastervarxcoeff, masterconsvars[j], (*weights)[i] * masterconscoeffs[j]);
         }
         else
         {
            SCIPhashmapSetImageReal(mapmastervarxcoeff, masterconsvars[j], coeff_ssrc + (*weights)[i] * masterconscoeffs[j]);
         }
      }

      SCIPfreeBufferArray(masterscip, &masterconsvars);
      SCIPfreeBufferArray(masterscip, &masterconscoeffs);
   }

   return SCIP_OKAY;
}

/** computes the (non-rounded) coefficients for the pricing variables used in the pricing constraints */
SCIP_RETCODE computePricingConssCoefficients(
   SCIP*             origscip,               /**< SCIP data structure of original problem */
   SCIP_CONS**       originalconss,          /**< constraints in the original problem */
   int*              selectedconssidx,       /**< indices of selected constraints */
   int               nselectedconss,         /**< number of selected constraints */
   SCIP_Real*        weights,                /**< weights for selected constraints */
   SCIP_HASHMAP*     mappricingvarxcoeff     /**< map to store the coefficient for each pricing variable */
   )
{
   SCIP_CONS* origcons;
   SCIP_VAR** origconsvars;
   SCIP_VAR* pricingvar;
   SCIP_Real* origconscoeffs;
   SCIP_Real coeff_pricing;
   SCIP_Bool success;
   int norigconsvars;
   int i;
   int j;

   assert(GCGisOriginal(origscip));
   assert(nselectedconss > 0);

   SCIPdebugMessage("compute the coefficients of the pricing variables\n");
   /* - compute w_1 * A_1j + ... + w_m * A_mj for each pricing variable x_j */
   for( i = 0; i < nselectedconss; i++ )
   {
      /* get all variables and their corresponding coefficients in the original constraint */
      origcons = originalconss[selectedconssidx[i]];
      SCIP_CALL( SCIPgetConsNVars(origscip, origcons, &norigconsvars, &success) );
      assert(success);
      if( norigconsvars == 0 )
      {
         SCIPdebugMessage("constraint has no variables\n");
         continue;
      }

      origconsvars = NULL;
      origconscoeffs = NULL;
      SCIPallocBufferArray(origscip, &origconsvars, norigconsvars);
      SCIPallocBufferArray(origscip, &origconscoeffs, norigconsvars);
      SCIP_CALL( SCIPgetConsVars(origscip, origcons, origconsvars, norigconsvars, &success) );
      assert(success);
      SCIP_CALL( SCIPgetConsVals(origscip, origcons, origconscoeffs, norigconsvars, &success) );
      assert(success);
      SCIPdebugMessage("number of orig vars %i\n", norigconsvars);
      for( j = 0; j < norigconsvars; j++ )
      {
         /* use the pricing variable corresponding to the original variable as key in map */
         //SCIPdebugMessage("orig var %s belongs to block %i\n", SCIPvarGetName(origconsvars[j]), GCGvarGetBlock(origconsvars[j]));
         pricingvar = GCGoriginalVarGetPricingVar(origconsvars[j]);
         if( pricingvar == NULL )
            continue;
         assert(pricingvar != NULL);
         coeff_pricing = SCIPhashmapGetImageReal(mappricingvarxcoeff, pricingvar);
         if( coeff_pricing == SCIP_INVALID )
         {
            //SCIPdebugMessage("pricing var %s: add %f * %f\n", SCIPvarGetName(pricingvar), weights[i], origconscoeffs[j]);
            SCIPhashmapSetImageReal(mappricingvarxcoeff, pricingvar, weights[i] * origconscoeffs[j]);
         }
         else
         {
            //SCIPdebugMessage("pricing var %s: add %f * %f to %f\n", SCIPvarGetName(pricingvar), weights[i], origconscoeffs[j], coeff_pricing);
            SCIPhashmapSetImageReal(mappricingvarxcoeff, pricingvar, coeff_pricing + weights[i] * origconscoeffs[j]);
         }
      }

      SCIPfreeBufferArray(origscip, &origconsvars);
      SCIPfreeBufferArray(origscip, &origconscoeffs);
   }
   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpSubsetrow)
{
   SCIP*                      origscip;            // original problem
   SCIP_SEPADATA*             sepadata;
   SCIP_ROW*                  ssrc;                // stores the subsetrow cut
   SCIP_HASHMAP*              mapmastervarxcoeff;  // maps master variable to its coefficient in the subset row
   SCIP_HASHMAP*              mappricingvarxcoeff; // maps pricing variable to its coefficient in its pricing constraint
   SCIP_CONS**                masterconss;         // a              : constraints in master problem
   SCIP_CONS**                originalconss;       // A              : constraints in original problem
   GCG_PRICINGMODIFICATION**  pricingmodifications;// pricing modifications associated with subsetrow cut
   SCIP_Bool                  success;
   SCIP_Real*                 weights;             // w in {-1/k, 1/k}^n: vector of weights for selected constraints
   SCIP_Real                  rhs_ssrc;            // stores the rhs of the subsetrow cut
   int                        npricingmodifcations;// number of pricing modifications associated with subsetrow cut
   int*                       selectedconssidx;    // indices of constraints used to construct subset row (weights are non-zero)
   int                        nmasterconss;        // m              : number of constraints in master problem
   int                        npricingproblems;    // K              : number of pricing problems
   int                        strategy;            // 0 = RANDOM, 1 = ZERO-HALF
   int                        nvars;
   int                        n;                   // n = |S| > 0    : number of constraints used to construct subset row
   int                        k;                   // k > 0          : defines the possible weights
   int                        i;
   int                        j;

   assert(scip != NULL);
   assert(result != NULL);
   assert(sepa != NULL);
   assert(GCGisMaster(scip));

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   *result = SCIP_DIDNOTFIND;
   SCIPdebugMessage("Start subsetrow cut exec\n");
   if( !sepadata->enable )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("master LP not solved to optimality, do no separation!\n");
      return SCIP_OKAY;
   }

   if( GCGgetNRelPricingprobs(origscip) < GCGgetNPricingprobs(origscip) )
   {
      SCIPdebugMessage("aggregated pricing problems, do no separation!\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* ensure to separate current sol */
   SCIP_CALL( GCGrelaxUpdateCurrentSol(origscip) );

   n = 3;
   k = 2;

   /* get info of master problem */
   originalconss = GCGgetOrigMasterConss(origscip);
   masterconss = GCGgetMasterConss(origscip);
   nmasterconss = GCGgetNMasterConss(origscip);
   if( n >= nmasterconss )
   {
      SCIPdebugMessage("Not enough constraints to build subsetrow: n = %i >= nmasterconss = %i!\n", n, nmasterconss);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( SCIPgetCurrentNode(scip) != SCIPgetRootNode(scip) || sepadata->ngeneratedcut > 0 )
   {
      SCIPdebugMessage("Only run on root\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }


   for( i = 0; i < 1; i ++ )
   {
      /* select n constraints and store their indices in selectedconssidx */
      SCIPdebugMessage("Select %i out of %i constraints\n", n, nmasterconss);
      selectedconssidx = NULL;
      SCIPallocBufferArray(scip, &selectedconssidx, n);
      SCIP_CALL( selectRandomRows(sepadata, nmasterconss, selectedconssidx, n) );

      /* determine the master variables, their coefficients and rhs for subset row */
      weights = NULL;
      mapmastervarxcoeff = NULL;
      SCIP_CALL( SCIPgetConsNVars(scip, masterconss[selectedconssidx[0]], &nvars, &success) );
      assert(success);
      SCIPallocBufferArray(scip, &weights, n);
      SCIPhashmapCreate(&mapmastervarxcoeff, SCIPblkmem(scip), nvars);
      SCIP_CALL( computeSubsetRowCoefficientsAndRHS(scip, masterconss, selectedconssidx, n, k,
                                                    &weights, mapmastervarxcoeff, &rhs_ssrc) );


      /* create the subsetrow cut */
      ssrc = NULL;
      SCIP_CALL( createSubsetRowCut(scip, &ssrc, mapmastervarxcoeff, rhs_ssrc, sepa) );
      assert(ssrc != NULL);


      /* determine the pricing variables and their coefficients for the pricing constraints */
      SCIP_CALL( SCIPgetConsNVars(origscip, originalconss[selectedconssidx[0]], &nvars, &success) );
      assert(success);
      SCIPhashmapCreate(&mappricingvarxcoeff, SCIPblkmem(scip), nvars);
      SCIP_CALL( computePricingConssCoefficients(origscip, originalconss, selectedconssidx, n, weights,
                                                      mappricingvarxcoeff) );


      /* create the pricing modification for every pricing problem */
      npricingmodifcations = 0;
      pricingmodifications = NULL;
      npricingproblems = GCGgetNPricingprobs(origscip);

      for( j = 0; j < npricingproblems; j++ )
      {
         SCIP*                      pricingproblem;
         GCG_PRICINGMODIFICATION*   pricingmodification = NULL;
         SCIP_CONS**                pricingconss = NULL;
         SCIP_VAR**                 pricingvars;
         SCIP_VAR*                  coeffvar; // y
         SCIP_Real                  pricingcoeff;
         int                        npricingvars;
         int                        l;

         pricingproblem = GCGgetPricingprob(origscip, j);
         if( pricingproblem == NULL )
            continue;

         npricingvars = SCIPgetNVars(pricingproblem);
         pricingvars = SCIPgetVars(pricingproblem);

         /* we add at most one constraint to each pricing problem */
         SCIP_CALL( SCIPallocBlockMemoryArray(pricingproblem, &pricingconss, 1) ); // freed via GCGpricingmodificationFree

         /* create (and capture) 'empty' pricing constraint: 0 <= ... <= 1 - EPSILON
          * - in case the dual of the ssrc is 0: lhs=0 ensures that y=0 */
         SCIP_CALL( SCIPcreateConsBasicLinear(pricingproblem, &pricingconss[0], "inferred_pricing",
                                              0, NULL, NULL, 0.0,
                                              1 - 0.0001) ); // released via GCGpricingmodificationFree -SCIPinfinity(pricingproblem)

         /* fill constraint such that 0 <= w^TAx <= 1 - EPSILON */
         for( l = 0; l < npricingvars; l++ )
         {
            if( !GCGvarIsPricing(pricingvars[l]) )
               continue;

            pricingcoeff = SCIPhashmapGetImageReal(mappricingvarxcoeff, pricingvars[l]);
            if( pricingcoeff == SCIP_INVALID )
               continue;

            if( !SCIPisZero(pricingproblem, pricingcoeff) )
            {
               SCIPdebugMessage("add pricing var %s with coefficient %f to constraint for pricing problem %i\n",
                                SCIPvarGetName(pricingvars[l]), pricingcoeff, j);
               SCIP_CALL( SCIPaddCoefLinear(pricingproblem, pricingconss[0], pricingvars[l], pricingcoeff) );
            }
         }

         /* if no variables were actually added, the constraint is useless and can be released */
         SCIP_CALL( SCIPgetConsNVars(pricingproblem, pricingconss[0], &nvars, &success) );
         if( nvars == 0  )
         {
            SCIPdebugMessage("constraint was empty --> release\n");
            SCIP_CALL( SCIPreleaseCons(pricingproblem, &pricingconss[0]) );
            SCIPfreeBlockMemoryArray(pricingproblem, &pricingconss, 1);
            continue;
         }

         /* create (and capture) y: 0 <= y <= inf */
         SCIPdebugMessage("create new (inferred) pricing variable y for pricing problem %i\n", j);
         coeffvar = NULL;
         SCIP_CALL( GCGcreateInferredPricingVar(pricingproblem, &coeffvar, "subsetrow_coef", 0.0,
                                                SCIPinfinity(pricingproblem), -1.0, SCIP_VARTYPE_INTEGER, j) ); // released in GCGpricingmodificationFree


         /* add y to constraint such that: 0 <= w^TAx - y <= 1 - EPSILON */
         SCIP_CALL( SCIPaddCoefLinear(pricingproblem, pricingconss[0], coeffvar, -1.0) );

         /* create pricing modifications containing y as the coeffvar and a single constraint */
         SCIPdebugMessage("create pricing modification for pricing problem %i\n", j);
         SCIP_CALL( GCGpricingmodificationCreate(scip, &pricingmodification, j, coeffvar, NULL,
                                                 0, pricingconss, 1) ); // released in GCGpricingmodificationFree

         /* ensure we have enough memory for all the pricing modifications */
         if( npricingmodifcations == 0 )
         {
            SCIPallocBlockMemoryArray(scip, &pricingmodifications, 1); // freed in GCGmastercutFree
         }
         else
         {
            SCIPreallocBlockMemoryArray(scip, &pricingmodifications, npricingmodifcations, npricingmodifcations + 1); // freed in GCGmastercutFree
         }

         pricingmodifications[npricingmodifcations] = pricingmodification;
         npricingmodifcations++;
      }
      SCIPdebugMessage("number of pricing mods: %i\n", npricingmodifcations);


      /* create a mastercutdata containing y, ssrc and the pricing modifications */
      SCIPdebugMessage("create mastercut for subsetrow cut\n");
      GCG_MASTERCUTDATA* mastercutdata = NULL;
      SCIP_CALL( GCGmastercutCreateFromRow(scip, &mastercutdata, ssrc, pricingmodifications, npricingmodifcations) ); // freed in GCGmastercutFree
      sepadata->ngeneratedcut++;

      /* add the subsetrow cut to the sepa store of the master problem and the generated cuts*/
      SCIP_CALL( SCIPaddRow(scip, ssrc, TRUE, &success) );
      SCIP_CALL( GCGaddCutToGeneratedCutsSepa(scip, mastercutdata, sepadata->sepaidx) );

      /* free used data structure */
      SCIPhashmapFree(&mapmastervarxcoeff);
      SCIPhashmapFree(&mappricingvarxcoeff);
      SCIPfreeBufferArray(scip, &selectedconssidx);
      SCIPfreeBufferArray(scip, &weights);
      SCIPdebugMessage("End subsetrow cut exec\n");
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of GCG separator
 */

/** get all the cuts generated by this master separator */
static
GCG_DECL_SEPAGETCOLCOEFFICIENTS(gcgsepaGetColCoefficientSubsetrow)
{
   SCIP_VARDATA* vardata;
   SCIP_VAR** solvars;
   SCIP_Real* solvals;
   int nsolvars;
   int i;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(sepa != NULL);
   assert(gcgcol != NULL);
   assert(ncoeffs >= 0);


   if( ncoeffs == 0 )
      return SCIP_OKAY;

   for( i = 0; i < ncoeffs; i++ )
   {
      (*coeffs)[i] = 0.0;
   }

   nsolvars = GCGcolGetNVars(gcgcol);
   solvars = GCGcolGetVars(gcgcol);
   solvals = GCGcolGetVals(gcgcol);

   for( i = 0; i < nsolvars; i++ )
   {
      /* ignore pricing variables which are not relatd to mastercuts*/
      if( !GCGvarIsInferredPricing(solvars[i]) )
         continue;

      vardata = SCIPvarGetData(solvars[i]);
      assert(vardata != NULL);

      /* ignore inferred variables which are not related to active cuts generated by separators */
      if( vardata->data.inferredpricingvardata.index == -1 || vardata->data.inferredpricingvardata.mastercutdata->type != GCG_MASTERCUTTYPE_ROW)
         continue;

      /* coefficient in mastercut for this column is y */
      SCIPdebugMessage("index %i\n", vardata->data.inferredpricingvardata.index);
      (*coeffs)[vardata->data.inferredpricingvardata.index] = solvals[i];
   }

   for( i = 0; i < ncoeffs; i++ )
   {
      SCIPdebugMessage("coeffs[%i] = %f\n", i, (*coeffs)[i]);
   }

   return SCIP_OKAY;
}

/** method for adding new master variable to cut */
static
GCG_DECL_SEPAGETVARCOEFFICIENT(gcgsepaGetVarCoefficientSubsetrow)
{
   SCIP* origscip;
   SCIP* pricingscip;
   GCG_PRICINGMODIFICATION** pricingmods;
   SCIP_CONS** pricingconss;
   SCIP_VAR** pricingconsvars;
   SCIP_Bool success;
   SCIP_Real* pricingconscoeffs;
   int npricingconsvars;
   int npricingmods;
   int i;
   int j;
   int k;


   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(sepa != NULL);
   assert(cut != NULL);

   *coef = 0.0;

   origscip = GCGgetOriginalprob(scip);
   npricingmods = GCGmastercutGetNPricingModifications(cut);
   pricingmods = GCGmastercutGetPricingModifications(cut);

   for( i = 0; i < npricingmods; i++ )
   {
      /* find the pricing modification that corresponds to the pricing problem which create this column/variable */
      if( probnr == GCGpricingmodificationGetBlock(pricingmods[i]) )
      {
         pricingscip = GCGgetPricingprob(origscip, probnr);
         pricingconss = GCGpricingmodificationGetAdditionalConss(pricingmods[i]);

         SCIP_CALL( SCIPgetConsNVars(pricingscip, pricingconss[0], &npricingconsvars, &success) );
         assert(success);
         SCIPallocBufferArray(scip, &pricingconsvars, npricingconsvars);
         SCIPallocBufferArray(scip, &pricingconscoeffs, npricingconsvars);
         SCIP_CALL( SCIPgetConsVars(pricingscip, pricingconss[0], pricingconsvars, npricingconsvars, &success) );
         assert(success);
         SCIP_CALL( SCIPgetConsVals(pricingscip, pricingconss[0], pricingconscoeffs, npricingconsvars, &success) );
         assert(success);

         /* compute floor(w^TAx) */
         for( j = 0; j < npricingconsvars; j++ )
         {
            if( GCGvarIsInferredPricing(pricingconsvars[j]) )
               continue;

            for( k = 0; k < nvars; k++ )
            {
               if( pricingconsvars[j] == vars[k] )
                  *coef += pricingconscoeffs[j] * vals[k];
            }
         }

         SCIPdebugMessage("compute coef %f\n", *coef);
         *coef = floor(*coef);
         SCIPfreeBufferArray(scip, &pricingconsvars);
         SCIPfreeBufferArray(scip, &pricingconscoeffs);
         return SCIP_OKAY;
      }
   }

   return SCIP_ERROR;
}

/** get all the cuts generated by this master separator */
static
GCG_DECL_SEPASETOBJECTIVE(gcgsepaSetObjectiveSubsetrow)
{
   SCIP* pricingproblem;
   SCIP* origscip;
   SCIP_VAR* coeffvar;
   GCG_PRICINGMODIFICATION** pricingmodifications;
   int npricingmodifications;
   int pricingblocknr;
   int i;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(sepa != NULL);
   assert(cut != NULL);

   origscip = GCGgetOriginalprob(scip);
   assert(origscip != NULL);

   npricingmodifications = GCGmastercutGetNPricingModifications(cut);
   pricingmodifications = GCGmastercutGetPricingModifications(cut);

   for( i = 0; i < npricingmodifications; i++ )
   {
      pricingblocknr = GCGpricingmodificationGetBlock(pricingmodifications[i]);
      pricingproblem = GCGgetPricingprob(origscip, pricingblocknr);
      coeffvar = GCGpricingmodificationGetCoefVar(pricingmodifications[i]);
      SCIP_CALL( SCIPchgVarObj(pricingproblem, coeffvar, -1.0 * dual) );
      assert(SCIPvarGetProbindex(coeffvar) != -1);
      (*realdualvalues)[pricingblocknr][SCIPvarGetProbindex(coeffvar)] = -1.0 * dual;

   }
   SCIPdebugMessage("dual is %f\n", dual);

   return SCIP_OKAY;
}


/*
 * MASTER separator specific interface methods
 */


/** initialization method of separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitSubsetrow)
{
   SCIP* origscip;
   SCIP_SEPADATA* sepadata;
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(GCGisMaster(scip));

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPdebugMessage("initialize gcg sepa subsetrow\n");
   /* creates the subsetrow gcg separator and includes in the relaxator data of the original problem */
   sepadata->sepaidx = GCGrelaxIncludeSeparator(origscip, sepa, gcgsepaGetVarCoefficientSubsetrow,
                                                gcgsepaGetColCoefficientSubsetrow, gcgsepaSetObjectiveSubsetrow);


   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the scip sepa of the subsetrow separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaSubsetrow(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create subsetrow separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   sepa = NULL;
   sepadata->sepaidx = 0;
   sepadata->ngeneratedcut = 0;
   sepadata->randnumgen = NULL;
   sepadata->enable = FALSE;
   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &(sepadata->randnumgen),
                               SCIPinitializeRandomSeed(scip, DEFAULT_RANDSEED), TRUE) );

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
                                   SEPA_USESSUBSCIP, SEPA_DELAY,
                                   sepaExeclpSubsetrow, sepaExecsolSubsetrow,
                                   sepadata) );
   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeSubsetrow) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitSubsetrow) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitSubsetrow) );

   SCIP_CALL( SCIPaddBoolParam(GCGmasterGetOrigprob(scip), "sepa/" SEPA_NAME "/enable", "enable subsetrow separator",
      &(sepadata->enable), FALSE, TRUE, NULL, NULL) );


   return SCIP_OKAY;
}