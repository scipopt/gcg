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

/**@file    sepa_subsetrow.c
 * @brief   subset row separator for master problem
 * @author  Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#include <assert.h>
#include <scip/cons_linear.h>

#include "gcg/mastersepacut.h"
#include "gcg/pricer_gcg.h"
#include "gcg/relax_gcg.h"
#include "gcg/scip_misc.h"
#include "gcg/sepa_subsetrow.h"
#include "gcg/zerohalf_selector.h"
#include "gcg/struct_sepagcg.h"
#include "gcg/event_mastersepacut.h"

#define SEPA_NAME           "subsetrow"
#define SEPA_DESC "subsetrow separator"
#define SEPA_PRIORITY               100
#define SEPA_FREQ                     1 /**< default frequency: 1 --> callback is executed for subproblems at every level*/
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                 TRUE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_RANDSEED             71 /**< initial seed for RandomNumber Generator which picks the constraint indices*/
#define DEFAULT_MAXROUNDS             1 /**< maximal number of subset row separation rounds per non-root node */
#define DEFAULT_MAXROUNDSROOT         1 /**< maximal number of subset row separation calls in the root node */
#define DEFAULT_MAXSEPACUTS         100 /**< maximal number of subset row cuts separated per call in non-root nodes */
#define DEFAULT_MAXSEPACUTSROOT    1000 /**< maximal number of subset row cuts separated per call in root node */
#define DEFAULT_MAXCUTCANDS        2000 /**< maximal number of subset row cuts in total */
#define DEFAULT_ONLYROOT          FALSE /**< only apply separator in root node */
#define DEFAULT_STRATEGY              0 /**< strategy which is used to determine which rows to consider for cut computation */
#define DEFAULT_N                     3 /**< number of rows used to create a new cut */
#define DEFAULT_K                     2 /**< inverse of weight used for cut generation */
#define DEFAULT_ONLYAGGREGATED     TRUE /**< run only on aggregated problems */
#define MAXAGGRLEN(nvars)           ((int)(0.1*(nvars) + 1000 )) // OG:+ 1000

/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   GCG*                    gcg;                 /**< GCG data structure */
   GCG_SEPA*               sepa;                /**< gcg master separator instance */
   SCIP_RANDNUMGEN*        randnumgen;          /**< random number generator (for strategy RANDOM) */
   SCIP_Bool               enable;              /**< is this separator enabled? */
   SCIP_Bool               onlyroot;            /**< indicates if separator should only be applied at root node */
   int                     ngeneratedcut;       /**< counts the total number of cuts generated */
   int                     maxrounds;           /**< maximal number of separation calls per non-root node */
   int                     maxroundsroot;       /**< maximal number of separation calls for root node */
   int                     maxsepacutsroot;     /**< number of cuts generated per separation call of root node */
   int                     maxsepacuts;         /**< number of cuts generated per separation call at non-root node */
   int                     maxcutcands;         /**< maximal number of cuts generated in total */
   int                     strategy;            /**< RANDOM (0), KOSTER-ET-A (1) */
   int                     n;                   /**< n = |S| > 0    : number of constraints used to construct cut */
   int                     k;                   /**< k > 0          : defines the weights 1/k */
   SCIP_Bool               onlyaggregated;      /**< indicates if separator should only run on aggregated problems */
   SCIP_CONS**             origmasterconss;     /**< array of (filtered) orig master conss (is only allocated if not equal to data of relaxator) */
   SCIP_CONS**             masterconss;         /**< array of (filtered) master conss (is only allocated if not equal to data of relaxator) */
   int                     nmasterconss;        /**< number of (filtered) master conss */
   SCIP_Bool               filteredmasterconss; /**< do origmasterconss and masterconss point to allocated arrays? */
};


/** randomly selects n different constraints from the master problem */
static
SCIP_RETCODE selectRandomRows(
   SCIP_RANDNUMGEN* randnumgen,                 /**< random number generator */
   int              nmasterconss,               /**< number of constraints in the master problem */
   int*             selectedmasterconssidx,     /**< pointer to store the indices of the selected constraints */
   int              n                           /**< number of constraints to be selected */
   )
{
   int i;
   int j;

   assert(n > 0 && n < nmasterconss);

   /* randomly select n indices out of [0, ..., nmasterconss - 1] */
   i = 0;
   while( i < n )
   {
      selectedmasterconssidx[i] = SCIPrandomGetInt(randnumgen, 0, nmasterconss-1);

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
#ifdef SCIP_DEBUG
   /*for( i = 0; i < n; i++ )
   {
      SCIPdebugMessage("select index %i\n", selectedmasterconssidx[i]);
      assert(0 <= selectedmasterconssidx[i] && selectedmasterconssidx[i] < nmasterconss);
   }*/
#endif

   return SCIP_OKAY;
}

/** create a new row for master problem and fill it with the variables (+ their coefficients) and right-hand side*/
static
SCIP_RETCODE createSubsetRowCut(
   SCIP*          masterscip,                   /**< SCIP data structure (master problem) */
   SCIP_ROW**     ssrc,                         /**< pointer to store subset row cut */
   SCIP_HASHMAP*  mapmastervarxcoeffs,          /**< maps master variables to their coefficient in subset-row cut */
   SCIP_Real      rhs_ssrc,                     /**< right hand side of subset row cut */
   SCIP_SEPA*     sepa                          /**< separator which creates subset row cut */
   )
{
   SCIP_SEPADATA* sepadata;
   char name[SCIP_MAXSTRLEN]; // name of the ssrc
   int nentries;
   int i;

   assert(GCGisMaster(masterscip));
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);

   /* create 'empty' subset row cut of form -inf <= ... <= rhs_ssrc
    * - local, removable, modifiable */
   rhs_ssrc = SCIPfeasFloor(masterscip, rhs_ssrc);
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ssrc_%i", sepadata->ngeneratedcut);
   SCIP_CALL( SCIPcreateEmptyRowSepa(masterscip, &(*ssrc), sepa, name, -SCIPinfinity(masterscip), rhs_ssrc, TRUE, TRUE, TRUE) );
   assert(ssrc != NULL);
   assert(*ssrc != NULL);

   /* fill the row with master variables and their coefficients */
   nentries = SCIPhashmapGetNEntries(mapmastervarxcoeffs);
   SCIP_CALL( SCIPcacheRowExtensions(masterscip, *ssrc) );
   for( i = 0; i < nentries; i++ )
   {
      SCIP_HASHMAPENTRY* entry;
      SCIP_Real varcoeff;
      SCIP_VAR* mastervar;

      entry = SCIPhashmapGetEntry(mapmastervarxcoeffs, i);
      if( entry == NULL )
         continue;

      mastervar = (SCIP_VAR*) SCIPhashmapEntryGetOrigin(entry);
      varcoeff = SCIPhashmapEntryGetImageReal(entry);
      varcoeff = SCIPfeasFloor(masterscip, varcoeff);
      if( !SCIPisZero(masterscip, varcoeff) )
         SCIP_CALL( SCIPaddVarToRow(masterscip, *ssrc, mastervar, varcoeff) );
   }
   SCIP_CALL( SCIPflushRowExtensions(masterscip, *ssrc) );

   return SCIP_OKAY;
}

/** computes the rhs (w^Tb) and the coefficient for each variable (w^Ta_p) in the cut (still non-rounded) */
static
SCIP_RETCODE computeSubsetRowCoefficientsAndRHS(
   SCIP*          masterscip,                 /**< SCIP data structure (master problem) */
   SCIP_SEPA*     sepa,                       /**< SCIP separator */
   SCIP_CONS**    masterconss,                /**< constraints of the master problem */
   int*           selectedconssidx,           /**< indices of the selected constraints */
   int            nselectedconss,             /**< number of selected constraints */
   SCIP_Real**    weights,                    /**< pointer to store weights of the selected constraints */
   SCIP_Real*     rhs_ssrc,                   /**< pointer to store rhs of subset row cut */
   SCIP_HASHMAP*  mapmastervarxcoeff          /**< maps master variable to its coefficient in subset-row cut */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_Bool success;
   int nmasterconsvars;
   int i;
   int j;

   assert(masterscip != NULL);
   assert(sepa != NULL);
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   *rhs_ssrc = 0;

   for( i = 0; i < nselectedconss; i++ )
   {
      SCIP_CONS* mastercons;
      SCIP_VAR** masterconsvars = NULL;
      SCIP_Real* masterconscoeffs = NULL;
      SCIP_Real  rhs_mastercons;
      SCIP_Real  lhs_mastercons;

      /* lhs <= ax <= rhs */
      SCIPdebugMessage("select constraint: %i\n", selectedconssidx[i]);
      mastercons = masterconss[selectedconssidx[i]];

      SCIPdebugPrintCons(masterscip, mastercons, NULL);
      lhs_mastercons = SCIPconsGetLhs(masterscip, mastercons, &success);
      assert(success);
      rhs_mastercons = SCIPconsGetRhs(masterscip, mastercons, &success);
      assert(success);

      /* we want all constraints to have form ax <= b to simplify the summation
       * - if constraint has form ax <= b: weight is  1/k
       * - if constraint has form b <= ax: weight is -1/k (to transform it to -ax <= -b) */
      if( SCIPisInfinity(masterscip, -lhs_mastercons ) ) // lhs = -inf --> constraint has form: ax <= b
      {
         (*weights)[i] = 1.0 / sepadata->k;
         (*rhs_ssrc) += (*weights)[i] * rhs_mastercons;
         SCIPdebugMessage("master constraint %s (ax <= %f) with weight %f\n", SCIPconsGetName(mastercons), rhs_mastercons, (*weights)[i]);
      }
      else if( SCIPisInfinity(masterscip, rhs_mastercons ) ) // lhs != -inf & rhs = inf --> constraint has form: b <= ax
      {
         (*weights)[i] = -1.0 / sepadata->k;
         (*rhs_ssrc) += (*weights)[i] * lhs_mastercons;
         SCIPdebugMessage("master constraint %s (%f <= ax) with weight %f\n", SCIPconsGetName(mastercons), lhs_mastercons, (*weights)[i]);
      }
      else // lhs != -inf & rhs != inf --> constraint has form: b_l <= ax <= b_r
      {
         /* we use the right side: ax <= b_r */
         (*weights)[i] = 1.0 / sepadata->k;
         (*rhs_ssrc) += (*weights)[i] * rhs_mastercons;
         SCIPdebugMessage("master constraint %s (%f <= ax <= %f) with weight %f\n", SCIPconsGetName(mastercons), lhs_mastercons, rhs_mastercons, (*weights)[i]);
      }

      /* get all variables and their corresponding coefficients in the master constraint */
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

      /* for each variable: add its weighted coefficient in this constraint to its coefficient for the subset row cut */
      for( j = 0; j < nmasterconsvars; j++ )
      {
         SCIP_Real varcoeff;

         varcoeff = SCIPhashmapGetImageReal(mapmastervarxcoeff, masterconsvars[j]);
         if( varcoeff == SCIP_INVALID )
            varcoeff = (*weights)[i] * masterconscoeffs[j];
         else
            varcoeff += (*weights)[i] * masterconscoeffs[j];

         SCIPhashmapSetImageReal(mapmastervarxcoeff, masterconsvars[j], varcoeff);
      }

      SCIPfreeBufferArrayNull(masterscip, &masterconsvars);
      SCIPfreeBufferArrayNull(masterscip, &masterconscoeffs);
   }

   return SCIP_OKAY;
}

/** computes the (non-rounded) coefficients for the pricing variables used in the pricing constraints */
static
SCIP_RETCODE computePricingConssCoefficients(
   GCG*              gcg,                    /**< GCG data structure */
   SCIP_CONS**       originalconss,          /**< constraints in the original problem */
   int*              selectedconssidx,       /**< indices of selected constraints */
   int               nselectedconss,         /**< number of selected constraints */
   SCIP_Real*        weights,                /**< weights for selected constraints */
   SCIP_HASHMAP*     mappricingvarxcoeff     /**< map to store the coefficient for each pricing variable */
   )
{
   SCIP* origscip;
   SCIP_CONS* origcons;
   SCIP_VAR** origconsvars;
   SCIP_VAR*  pricingvar;
   SCIP_Real* origconscoeffs;
   SCIP_Real  coeff_pricing;
   int norigconsvars;
   int i;
   int j;

   assert(nselectedconss > 0);

   origscip = GCGgetOrigprob(gcg);

   SCIPdebugMessage("compute the coefficients of the pricing variables\n");
   /* - compute w_1 * A_1j + ... + w_m * A_mj for each pricing variable x_j */
   for( i = 0; i < nselectedconss; i++ )
   {
      /* get all variables and their corresponding coefficients in the original constraint */
      origcons = originalconss[selectedconssidx[i]];
      SCIPdebugPrintCons(origscip, origcons, NULL);
      norigconsvars = GCGconsGetNVars(origscip, origcons);

      if( norigconsvars == 0 )
      {
         SCIPdebugMessage("constraint has no variables\n");
         continue;
      }

      origconsvars = NULL;
      origconscoeffs = NULL;
      SCIPallocBufferArray(origscip, &origconsvars, norigconsvars);
      SCIPallocBufferArray(origscip, &origconscoeffs, norigconsvars);
      SCIP_CALL( GCGconsGetVars(origscip, origcons, origconsvars, norigconsvars) );
      SCIP_CALL( GCGconsGetVals(origscip, origcons, origconscoeffs, norigconsvars) );

      for( j = 0; j < norigconsvars; j++ )
      {
         assert(GCGvarIsOriginal(origconsvars[j]));

         /* use the pricing variable corresponding to the original variable as key in map */
         pricingvar = GCGoriginalVarGetPricingVar(origconsvars[j]);
         if( pricingvar == NULL )
         {
            SCIPdebugMessage("original variable %s does not have a corresponding pricing var!\n", SCIPvarGetName(origconsvars[j]));
            continue;
         }
         assert(pricingvar != NULL);
         assert(GCGvarGetBlock(pricingvar) >= 0 && GCGvarGetBlock(origconsvars[j]) >= 0);

         /* in case we have aggregated pricing problems:
          * - multiple original variables in the same constraint are mapped to the same pricing variable
          * - we only need to consider one of those original variables
          * --> only consider the original variables associated with the relevant pricing problem block */
         if( GCGvarGetBlock(pricingvar) != GCGvarGetBlock(origconsvars[j]))
         {
            assert(!GCGisPricingprobRelevant(gcg, GCGvarGetBlock(origconsvars[j])));
            continue;
         }
         assert(GCGisPricingprobRelevant(gcg, GCGvarGetBlock(origconsvars[j])));

         coeff_pricing = SCIPhashmapGetImageReal(mappricingvarxcoeff, pricingvar);
         if( coeff_pricing == SCIP_INVALID )
            SCIP_CALL( SCIPhashmapSetImageReal(mappricingvarxcoeff, pricingvar, weights[i] * origconscoeffs[j]) );
         else
            SCIP_CALL( SCIPhashmapSetImageReal(mappricingvarxcoeff, pricingvar, coeff_pricing + weights[i] * origconscoeffs[j]) );
      }

      SCIPfreeBufferArrayNull(origscip, &origconsvars);
      SCIPfreeBufferArrayNull(origscip, &origconscoeffs);
   }
   return SCIP_OKAY;
}

/** select the indices of constraints to use for construction of cuts randomly */
static
SCIP_RETCODE selectConstraintsRandom(
   SCIP*                   masterscip,       /**< SCIP data structure (master problem) */
   GCG_CUTINDICES***       cutindices,       /**< pointer to store array of constraint indices for cuts */
   int*                    ncutindices,      /**< pointer to store the number of cuts for which constraint indices were found */
   int                     maxcuts,          /**< maximal number of cuts for which indices should be computed */
   SCIP_SEPADATA*          sepadata          /**< data of separator */
   )
{
   int nconss;
   int nmasterconss;
   SCIP_RANDNUMGEN* randnumgen;

   nconss = sepadata->n;
   nmasterconss = sepadata->nmasterconss;
   randnumgen = sepadata->randnumgen;

   assert(masterscip != NULL);
   assert(randnumgen != NULL);

   int i;
   *ncutindices = 0;

   for( i = 0; i < maxcuts; i++ )
   {
      GCG_CUTINDICES* cutindex;
      int* selectedindices;

      SCIPallocBufferArray(masterscip, &selectedindices, nconss);
      SCIP_CALL( selectRandomRows(randnumgen, nmasterconss, selectedindices, nconss) );

      SCIP_CALL( GCGcreateCutIndicesFromArray(masterscip, &cutindex, nconss, selectedindices) );
      (*cutindices)[(*ncutindices)] = cutindex;
      (*ncutindices)++;
   }

   assert(*ncutindices == maxcuts);

   return SCIP_OKAY;
}

/** select the indices of constraints to use for construction of cuts using strategy devised for zero-half-cuts
 * @to-do: finished implementing alternative strategy */
static
SCIP_RETCODE selectConstraintsKosterEtAl(
   GCG*                    gcg,           /**< GCG data structure */
   GCG_CUTINDICES***       cutindices,    /**< pointer to store array of constraint indices for cuts */
   int*                    ncutindices,   /**< pointer to store the number of cuts for which constraint indices were found */
   int                     ncalls,        /**< number of time separator has bin called at current node */
   SCIP_Bool               allowlocal,    /**< local cuts allowed? */
   int                     depth,         /**< depth of current node */
   int                     maxcuts,       /**< maximum number of cuts to be generated */
   SCIP_SOL*               sol,           /**< current best solution, if available (else NULL) */
   SCIP_SEPADATA*          sepadata       /**< data of separator */
   )
{
   GCG_ZEROHALFDATA zhdata;

   zhdata.maxroundsroot = 20;
   zhdata.maxrounds = 5;
   zhdata.maxslack = 0.0;
   zhdata.maxslackroot = 0.0;
   zhdata.minviol = 0.1;
   zhdata.dynamiccuts = TRUE;
   zhdata.maxrowdensity = 0.05;
   zhdata.densityoffset = 100;
   zhdata.infeasible = FALSE;
   zhdata.nreductions = 0;
   zhdata.nmasterconss = sepadata->nmasterconss;
   zhdata.origmasterconss = sepadata->origmasterconss;
   zhdata.masterconss = sepadata->masterconss;
   *ncutindices = 0;

   SCIP_CALL( GCGselectConstraintsZeroHalf(gcg, sol, allowlocal, depth, &zhdata, ncalls, maxcuts, cutindices, ncutindices) );

   return SCIP_OKAY;
}

/** creates a subset row cut */
static
SCIP_RETCODE createCut(
   SCIP*             masterscip,    /**< SCIP data structure (master problem) */
   GCG_CUTINDICES*   cutindex,      /**< indices of the master constraints to use for the construction of subset row cut */
   SCIP_SEPA*        sepa,          /**< subset row separator */
   int               nmastervars,   /**< number of variables in the master problem*/
   SCIP_CONS**       masterconss,   /**< constraints of the master problem */
   SCIP_Real**       weights,       /**< sign (1, -1) with which to multiply each selected constraint before adding to cut*/
   SCIP_ROW**        ssrc,          /**< pointer to store subset row cut */
   SCIP_HASHMAP*     mapmastervarxcoeff /**< temporary map to store coefficients */
   )
{
   
   SCIP_Real      rhs_ssrc;

   assert(masterscip != NULL);
   /* determine the master variables, their coefficients and rhs for subset row (non-rounded) */
   SCIP_CALL( computeSubsetRowCoefficientsAndRHS(masterscip, sepa, masterconss, cutindex->indices, cutindex->nindices,
                                                weights, &rhs_ssrc, mapmastervarxcoeff) );

   /* create the subset row cut */
   SCIP_CALL( createSubsetRowCut(masterscip, ssrc, mapmastervarxcoeff, rhs_ssrc, sepa) );

   SCIPhashmapRemoveAll(mapmastervarxcoeff);

   return SCIP_OKAY;
}

/** creates the master cut data for the subset row cut */
static
GCG_EXTENDEDMASTERCONSDATA* createMastercutData(
   GCG*           gcg,                    /**< GCG data structure */
   SCIP_ROW*      ssrc,                   /**< subset row cut */
   GCG_MASTERSEPACUT* mastercut,          /**< master sepa cut */
   int            npricingproblems,       /**< number of pricing problems */
   SCIP_SEPADATA* sepadata,               /**< subset row separator data */
   SCIP_HASHMAP*  mappricingvarxcoeff     /**< map storing the coefficients for the pricing constraints for subset row cut */
   )
{
   SCIP* masterscip;
   GCG_EXTENDEDMASTERCONSDATA* mastercutdata = NULL;         // master cut data for subset row cut
   GCG_PRICINGMODIFICATION** pricingmodifications = NULL;  // pricing modifications associated with cut
   GCG_PRICINGMODIFICATION** pricingmodificationsbuffer = NULL;
   SCIP_Bool success;
   char name[SCIP_MAXSTRLEN];         // name of the subset row cut
   int npricingmodifications = 0;    // counter for number of pricing modifications
   int j;

   masterscip = GCGgetMasterprob(gcg);

   SCIPallocBufferArray(masterscip, &pricingmodificationsbuffer, GCGgetNRelPricingprobs(gcg));

   /* create the pricing modification for every (relevant) pricing problem */
   for( j = 0; j < npricingproblems; j++ )
   {
      SCIP* pricingproblem;
      /* we add at most one constraint to each pricing problem */
      SCIP_CONS* pricingcons;
      SCIP_CONS** pricingconss = NULL;
      SCIP_VAR** pricingvars;
      SCIP_VAR* coeffvar = NULL; // y
      SCIP_Real pricingcoeff;
      int npricingvars;
      int nvars;
      int l;

      /* in case of aggregated pricing problems, we skip the non-representative ones */
      pricingproblem = GCGgetPricingprob(gcg, j);
      if( pricingproblem == NULL || !GCGisPricingprobRelevant(gcg, j) )
         continue;

      assert(GCGisPricingprobRelevant(gcg, j));
      assert(SCIPgetObjsense(pricingproblem) == SCIP_OBJSENSE_MINIMIZE);

      npricingvars = SCIPgetNVars(pricingproblem);
      pricingvars = SCIPgetVars(pricingproblem);

      /* create (and capture) 'empty' pricing constraint: -inf <= ... <= 1 - EPSILON */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pp%i_cons_ssrc_%i", j, sepadata->ngeneratedcut);
      SCIPcreateConsBasicLinear(pricingproblem, &pricingcons, name, 0, NULL, NULL, -SCIPinfinity(pricingproblem),
                                1 - 0.0001); // released via GCGpricingmodificationFree

      /* fill constraint such that -inf <= w^TAx <= 1 - EPSILON */
      for( l = 0; l < npricingvars; l++ )
      {
         assert(GCGvarIsPricing(pricingvars[l]));
         pricingcoeff = SCIPhashmapGetImageReal(mappricingvarxcoeff, pricingvars[l]);
         if( pricingcoeff == SCIP_INVALID )
            continue;

         if( !SCIPisZero(pricingproblem, pricingcoeff) ) // @todo: == 0.0??
            SCIPaddCoefLinear(pricingproblem, pricingcons, pricingvars[l], pricingcoeff);
      }

      /* if no variables were actually added, the constraint is useless and can be released */
      SCIPgetConsNVars(pricingproblem, pricingcons, &nvars, &success);
      if( nvars == 0  )
      {
         SCIPdebugMessage("constraint was empty --> release\n");
         SCIPreleaseCons(pricingproblem, &pricingcons);
         continue;
      }

      /* create (and capture) y: -inf <= y <= inf (integer) */
      SCIPdebugMessage("create new (inferred) pricing variable y for pricing problem %i\n", j);
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pp%i_y_ssrc_%i", j, sepadata->ngeneratedcut);
      GCGcreateInferredPricingVar(pricingproblem, &coeffvar, name, -SCIPinfinity(pricingproblem),
                                  SCIPinfinity(pricingproblem), TRUE, -1.0, SCIP_VARTYPE_INTEGER, j) ; // released in GCGpricingmodificationFree
      assert(coeffvar != NULL);

      /* add y to constraint such that: -inf <= w^TAx - y <= 1 - EPSILON  <=> w^TAx - 1 + EPSILON <= y */
      SCIPaddCoefLinear(pricingproblem, pricingcons, coeffvar, -1.0);
      SCIPdebugPrintCons(pricingproblem, pricingcons, NULL);

      assert(GCGgetNRelPricingprobs(gcg) > npricingmodifications);

      /* create pricing modifications containing y as the coeffvar and the single constraint we created */
      SCIPallocBlockMemoryArray(masterscip, &pricingconss, 1); // freed via GCGpricingmodificationFree
      pricingconss[0] = pricingcons;
      GCGpricingmodificationCreate(gcg, &pricingmodificationsbuffer[npricingmodifications], j, coeffvar, NULL, 0, pricingconss, 1); // released in GCGpricingmodificationFree

      npricingmodifications++;
   }
   SCIPdebugMessage("number of pricing mods: %i\n", npricingmodifications);

   if( npricingmodifications > 0 )
   {
      SCIPduplicateBlockMemoryArray(masterscip, &pricingmodifications, pricingmodificationsbuffer, npricingmodifications);
   }
   SCIPfreeBufferArray(masterscip, &pricingmodificationsbuffer);

   /* create master cut data containing y, ssrc and the pricing modifications */
   GCGextendedmasterconsCreateForSepamastercut(gcg, &mastercutdata, mastercut, ssrc, pricingmodifications, npricingmodifications); // freed in GCGextendedmasterconsFree
   sepadata->ngeneratedcut++;

   return mastercutdata;
}

/*
 * Callback methods of separator
 */

#define sepaCopySubsetrow NULL

static
SCIP_DECL_SEPAEXIT(sepaExitSubsetrow)
{
   SCIP_SEPADATA* sepadata;
   assert(sepa != NULL);
   sepadata = SCIPsepaGetData(sepa);

   if( sepadata->origmasterconss != NULL && sepadata->filteredmasterconss )
   {
      SCIPfreeBlockMemoryArray(scip, &sepadata->origmasterconss, sepadata->nmasterconss);
      SCIPfreeBlockMemoryArray(scip, &sepadata->masterconss, sepadata->nmasterconss);
   }

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeSubsetrow)
{
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);
   assert(sepa != NULL);

   SCIPdebugMessage("free separator data for subset row separator\n");
   sepadata = SCIPsepaGetData(sepa);
   SCIPfreeRandom(scip, &(sepadata->randnumgen));
   SCIPfreeBlockMemory(scip, &sepadata);

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpSubsetrow)
{
   GCG* gcg;
   SCIP_SEPADATA* sepadata;
   SCIP_HASHMAP* mappricingvarxcoeff = NULL;    // maps pricing variable to its coefficient in its pricing constraint
   GCG_CUTINDICES** cutindices = NULL;          // Ss             : for every cut contains the S defining it
   SCIP_Bool success;
   SCIP_Bool isroot;
   int ncutindices;                  // |Ss|           : number of S-sets created
   int npricingproblems;             // |K|            : number of pricing problems
   int nmastervars;                  // |X'|           : number of variables in (reduced) master problem
   int ncontmastervars;              // number of continuous variables in (reduced) master problem
   int maxcuts;                      // number of cuts to be generated in a separation round
   int ncalls;                       // number of times separator was called in this node
   int curfound;                     // number of generated subset-row cuts
   int i;
   GCG_MASTERSEPACUT* mastercut;
   SCIP_HASHMAP* mapmastervarxcoeff = NULL;
   SCIP* origprob;

   assert(scip != NULL);
   assert(result != NULL);
   assert(sepa != NULL);
   assert(GCGisMaster(scip));

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   gcg = sepadata->gcg;
   origprob = GCGgetOrigprob(gcg);

   isroot = SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);
   *result = SCIP_DIDNOTFIND;

   if( !sepadata->enable )
   {
      SCIPdebugMessage("subset row separator is not enabled.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("master LP not solved to optimality, do no separation!\n");
      return SCIP_OKAY;
   }

   if( ncalls >= sepadata->maxrounds )
   {
      SCIPdebugMessage("exceeded max rounds for this node.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( (!isroot && sepadata->onlyroot) || (!isroot && !allowlocal))
   {
      SCIPdebugMessage("subset row separator is only configured to run on root node.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( SCIPnodeGetType(SCIPgetCurrentNode(scip)) == SCIP_NODETYPE_PROBINGNODE ||
       SCIPnodeGetType(SCIPgetCurrentNode(scip)) == SCIP_NODETYPE_REFOCUSNODE )
   {
      SCIPdebugMessage("subset row separator does not run on probing or refocus nodes.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( sepadata->ngeneratedcut >= sepadata->maxcutcands )
   {
      SCIPdebugMessage("already generated the maximal number of cuts.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( sepadata->onlyaggregated && GCGgetNRelPricingprobs(gcg) == GCGgetNPricingprobs(gcg) )
   {
      SCIPdebugMessage("subset row separator is skipped since no pricing problem is aggregated.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* get info of master problem */
   npricingproblems = GCGgetNPricingprobs(gcg);
   nmastervars = SCIPgetNVars(scip);
   ncontmastervars = SCIPgetNContVars(scip);

   if( sepadata->origmasterconss == NULL )
   {
      if( sepadata->onlyaggregated )
      {
         SCIP_CONS** origmasterconss;
         SCIP_CONS** masterconss;
         int nmasterconss;
         SCIP_VAR** vars;
         int nvars;
         int j;
         SCIP_Bool add;

         nmasterconss = GCGgetNMasterConss(gcg);
         origmasterconss = GCGgetOrigMasterConss(gcg);
         masterconss = GCGgetMasterConss(gcg);

         sepadata->filteredmasterconss = TRUE;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sepadata->origmasterconss, nmasterconss) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sepadata->masterconss, nmasterconss) );

         /* ignore constraints that contain variables of not aggregated pricing problems */
         sepadata->nmasterconss = 0;
         for( i = 0; i < nmasterconss; ++i )
         {
            SCIP_VAR* pricingvar;
            add = TRUE;
            nvars = GCGconsGetNVars(origprob, origmasterconss[i]);
            SCIPallocBufferArray(scip, &vars, nvars);
            GCGconsGetVars(origprob, origmasterconss[i], vars, nvars);
            for( j = 0; j < nvars; ++j )
            {
               assert(GCGvarIsOriginal(vars[j]));
               pricingvar = GCGoriginalVarGetPricingVar(vars[j]);
               if( pricingvar != NULL && GCGpricingVarGetNOrigvars(pricingvar) <= 1)
               {
                  add = FALSE;
                  break;
               }
            }
            SCIPfreeBufferArray(scip, &vars);

            if( add )
            {
               sepadata->origmasterconss[sepadata->nmasterconss] = origmasterconss[i];
               sepadata->masterconss[sepadata->nmasterconss] = masterconss[i];
               sepadata->nmasterconss++;
            }
         }
      }
      else
      {
         sepadata->nmasterconss = GCGgetNMasterConss(gcg);
         sepadata->origmasterconss = GCGgetOrigMasterConss(gcg);
         sepadata->masterconss = GCGgetMasterConss(gcg);
      }
   }

   if( sepadata->n >= sepadata->nmasterconss )
   {
      SCIPdebugMessage("not enough constraints to build subset row cut: n = %i >= number of master constraints = %i!\n", sepadata->n, sepadata->nmasterconss);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* we only apply the separator on pure integer master problems */
   if( ncontmastervars > 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* determine the number of cuts to generate based on node type */
   if( isroot )
      maxcuts = sepadata->maxsepacutsroot;
   else
      maxcuts = sepadata->maxsepacuts;


   /* select which constraints to use for new subset row cuts */
   SCIP_CALL( SCIPallocBufferArray(scip, &cutindices, maxcuts) );
   SCIP_CALL( SCIPhashmapCreate(&mappricingvarxcoeff, SCIPblkmem(scip), sepadata->n) );
   ncutindices = 0;
   if( sepadata->strategy == 0 )
      SCIP_CALL( selectConstraintsRandom(scip, &cutindices, &ncutindices, maxcuts, sepadata) );
   else if( sepadata->strategy == 1 )
      SCIP_CALL( selectConstraintsKosterEtAl(gcg, &cutindices, &ncutindices, ncalls, allowlocal, depth, maxcuts, NULL, sepadata) );
   else
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   curfound = 0;
   SCIP_CALL( SCIPhashmapCreate(&mapmastervarxcoeff, SCIPblkmem(scip), sepadata->n) );
   for( i = ncutindices - 1; i >= 0; i-- )
   {
      GCG_EXTENDEDMASTERCONSDATA*   mastercutdata = NULL;
      SCIP_ROW*            ssrc = NULL;
      SCIP_Real*           weights = NULL;

      SCIPallocBufferArray(scip, &weights, cutindices[i]->nindices);

      // create the subset row cut based on the selected indices
      SCIP_CALL( createCut(scip, cutindices[i], sepa, nmastervars, sepadata->masterconss, &weights, &ssrc, mapmastervarxcoeff) );

      // row is empty --> useless
      if( ssrc == NULL || SCIProwGetNNonz(ssrc) == 0 )
      {
         SCIPdebugMessage("created an empty row: release row\n");
         SCIPfreeBufferArrayNull(scip, &weights);
         SCIP_CALL( GCGfreeCutIndices(scip, &cutindices[i]) );
         if( ssrc != NULL )
            SCIP_CALL( SCIPreleaseRow(scip, &ssrc) );

         continue;
      }
#ifdef SCIP_DEBUG
      SCIPprintRow(scip, ssrc, NULL);
#endif
      // determine the pricing variables and their coefficients for the pricing constraints
      SCIP_CALL( computePricingConssCoefficients(gcg, sepadata->origmasterconss, cutindices[i]->indices, cutindices[i]->nindices,
                                                 weights, mappricingvarxcoeff) );

      // add cut to separation store and corresponding master separator cut to sepacut event handler
      SCIP_CALL( GCGcreateChvatalGomoryCut(gcg, &mastercut, sepadata->sepa, NULL, weights, cutindices[i]->indices, cutindices[i]->nindices) );
      assert(mastercut != NULL);
      mastercutdata = createMastercutData(gcg, ssrc, mastercut, npricingproblems, sepadata, mappricingvarxcoeff);
      SCIP_CALL( SCIPaddRow(scip, ssrc, FALSE, &success) );
      SCIP_CALL( GCGaddGeneratedMastersepacut(gcg, mastercutdata) );
      curfound++;

      // cleanup
      SCIPfreeBufferArrayNull(scip, &weights);
      SCIP_CALL( GCGfreeCutIndices(scip, &cutindices[i]) );
      SCIP_CALL( SCIPhashmapRemoveAll(mappricingvarxcoeff) );
   }

   SCIPhashmapFree(&mapmastervarxcoeff);

   SCIPfreeBufferArrayNull(scip, &cutindices);
   SCIPhashmapFree(&mappricingvarxcoeff);

   if( curfound > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** compute cut coefficient for a column */
static
GCG_DECL_SEPAGETCOLCOEFFICIENT(sepaGetColCoefficientSubsetrow)
{
   assert(strcmp(SEPA_NAME, SCIPsepaGetName(GCGsepaGetScipSeparator(sepa))) == 0);
   if( gcgcol != NULL )
      return GCGchvatalGomoryCutGetColumnCoefficient(gcg, cut, gcgcol, coef);
   else
      return GCGchvatalGomoryCutGetVariableCoefficient(gcg, cut, vars, vals, nvars, probnr, coef);
}

/** modifies the objective values of the pricing variables affected by the master cut */
static
GCG_DECL_SEPASETOBJECTIVE(sepaSetObjectiveSubsetrow)
{
   assert(strcmp(SEPA_NAME, SCIPsepaGetName(GCGsepaGetScipSeparator(sepa))) == 0);
   return GCGchvatalGomorySetPricingObjectives(gcg, cut, dual);
}

/** modifies outdated column to respect cut */
static
GCG_DECL_SEPAADJUSTCOL(sepaAdjustColSubsetrow)
{
   assert(strcmp(SEPA_NAME, SCIPsepaGetName(GCGsepaGetScipSeparator(sepa))) == 0);
   return GCGchvatalGomoryAdjustGCGColumn(gcg, cut, gcgcol);
}

/** deletes mastersepacutdata */
static
GCG_DECL_SEPAMASTERCUTDELETE(sepaMastercutDeleteSubsetrow)
{
   assert(strcmp(SEPA_NAME, SCIPsepaGetName(GCGsepaGetScipSeparator(sepa))) == 0);
   return GCGfreeChvatalGomoryCutData(gcg, data);
}


/** initialization method of separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitSubsetrow)
{
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);
   assert(sepa != NULL);
   assert(GCGisMaster(scip));

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   sepadata->ngeneratedcut = 0;

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the scip separator of the subset row separator and includes it in master SCIP*/
SCIP_RETCODE SCIPincludeSepaSubsetrow(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* scip;
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;
   GCG_SEPA* gcgsepa;
   SCIP* origscip;

   scip = GCGgetMasterprob(gcg);

   /* create subsetrow separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   sepadata->gcg = gcg;
   sepadata->ngeneratedcut = 0;
   sepadata->randnumgen = NULL;
   sepadata->enable = FALSE;
   sepadata->origmasterconss = NULL;
   sepadata->masterconss = NULL;
   sepadata->nmasterconss = 0;
   sepadata->filteredmasterconss = FALSE;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &(sepadata->randnumgen),
                               SCIPinitializeRandomSeed(scip, DEFAULT_RANDSEED), TRUE) );

   SCIP_CALL( GCGrelaxIncludeSepa(gcg, &sepa, &gcgsepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
      SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpSubsetrow, NULL, sepadata, sepaAdjustColSubsetrow, sepaGetColCoefficientSubsetrow,
      sepaSetObjectiveSubsetrow, sepaMastercutDeleteSubsetrow) );

   assert(sepa != NULL);
   assert(gcgsepa != NULL);
   sepadata->sepa = gcgsepa;

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeSubsetrow) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitSubsetrow) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitSubsetrow) );

   origscip = GCGgetOrigprob(gcg);
   assert(origscip != NULL);

   /* define setting parameters */
   SCIP_CALL( SCIPaddBoolParam(origscip, "sepa/" SEPA_NAME "/enable", "enable subsetrow separator",
      &(sepadata->enable), FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "sepa/" SEPA_NAME "/maxrounds", "maximal number of subsetrow separation rounds per node (-1: unlimited)",
      &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "sepa/" SEPA_NAME "/maxroundsroot", "maximal number of subsetrow separation rounds in the root node (-1: unlimited)",
      &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "sepa/" SEPA_NAME "/maxsepacuts", "maximal number of subsetrow cuts separated per separation round",
      &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "sepa/" SEPA_NAME "/maxsepacutsroot", "maximal number of subsetrow cuts separated per separation round in the root node",
      &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "sepa/" SEPA_NAME "/maxcutcands", "maximal number of total subsetrow cuts considered",
      &sepadata->maxcutcands, FALSE, DEFAULT_MAXCUTCANDS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origscip, "sepa/" SEPA_NAME "/onlyroot", "apply subsetrow separator only on root",
      &(sepadata->onlyroot), FALSE, DEFAULT_ONLYROOT, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "sepa/" SEPA_NAME "/strategy", "RANDOM (0)",
      &sepadata->strategy, FALSE, DEFAULT_STRATEGY, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip,"sepa/" SEPA_NAME "/n", "number of rows used to create a new subset row cut ",
      &(sepadata->n), FALSE, DEFAULT_N, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(origscip, "sepa/" SEPA_NAME "/l", "weight used to create new subset row cut",
      &(sepadata->k), FALSE, DEFAULT_K, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(origscip, "sepa/" SEPA_NAME "/onlyaggregated", "apply subsetrow separator only on aggregated problems",
      &(sepadata->onlyaggregated), FALSE, DEFAULT_ONLYAGGREGATED, NULL, NULL) );

   return SCIP_OKAY;
}


