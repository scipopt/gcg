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
#define SEPA_PRIORITY                 1000
#define SEPA_FREQ                    1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define STARTMAXCUTS                50
#define DEFAULT_RANDSEED       71

#define DEFAULT_MAXROUNDS             5 /**< maximal number of zerohalf separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        20 /**< maximal number of zerohalf separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS         500 /**< maximal number of zerohalf cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT      10 /**< maximal number of zerohalf cuts separated per separation round in root node */
#define DEFAULT_MAXCUTCANDS        2000 /**< maximal number of zerohalf cuts considered per separation round */
#define DEFAULT_INITSEED         0x5EED /**< default initial seed used for random tie-breaking in cut selection */
#define DEFAULT_ONLYROOT           TRUE /**< only apply separator in root node */
#define DEFAULT_STRATEGY              0 /**< strategy which is used to determine which rows to considers cut computation */
#define DEFAULT_N                     3 /**< number of rows used to create a new cut */
#define DEFAULT_K                     2 /**< inverse of weight used for cut generation*/

/*
 * Data structures
 */

/* TODO: fill in the necessary separator data */

/** separator data */
struct SCIP_SepaData
{
   SCIP_RANDNUMGEN*        randnumgen;          /**< random number generator (for strategy RANDOM) */
   SCIP_Bool               enable;              /**< is this separator enabled? */
   SCIP_Bool               onlyroot;            /**< indicates if separator should only be applied at root node */
   int                     sepaidx;             /**< index of the separator in relaxdata->separators */
   int                     ngeneratedcut;       /**< counts the number of cuts generated */
   int                     maxrounds;           /**< maximal number of separation rounds per node */
   int                     maxsepacutsroot;     /**< number of cuts generated per separation round of root node */
   int                     maxroundsroot;       /**< maximal number of separation rounds for root node */
   int                     maxsepacuts;         /**< number of cuts generated per separation round */
   int                     initseed;            /**< */
   int                     maxcutcands;         /**< */
   int                     strategy;            /**< RANDOM (0), KOSTER-AL (1) */
   int                     n;                   /**< k > 0          : defines the possible weights 1/k */
   int                     k;                   /**< n = |S| > 0    : number of constraints used to construct subset row */
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
   SCIP_SEPADATA* sepadata;
   char name[SCIP_MAXSTRLEN]; // name of the ssrc
   int nssrcvars;
   int i;

   assert(GCGisMaster(masterscip));
   assert(sepa != NULL);
   sepadata = SCIPsepaGetData(sepa);

   SCIPdebugMessage("create subset row cut\n");
   /* create 'empty' subset row cut of form -inf <= ... <= rhs_ssrc
    * - local, removable, modifiable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ssrc_%i", sepadata->ngeneratedcut);
   SCIP_CALL( SCIPcreateEmptyRowSepa(masterscip, &(*ssrc), sepa, name, -SCIPinfinity(masterscip), SCIPfeasFloor(masterscip, rhs_ssrc), TRUE, TRUE, FALSE) );
   assert(ssrc != NULL);
   assert(*ssrc != NULL);

   /* fill the row with master variables and their coefficients */
   SCIP_CALL( SCIPcacheRowExtensions(masterscip, *ssrc) );
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
      assert(GCGvarIsMaster(mastervar));
      ssrc_coeff = SCIPhashmapEntryGetImageReal(entry);
      ssrc_coeff = SCIPfeasFloor(masterscip, ssrc_coeff);
      if( !SCIPisZero(masterscip, ssrc_coeff) )
      {
         SCIPdebugMessage("Add master variable %s with coefficient %f to subset row cut\n", SCIPvarGetName(mastervar), ssrc_coeff);
         SCIP_CALL( SCIPaddVarToRow(masterscip, *ssrc, mastervar, ssrc_coeff) );
      }
   }
   SCIP_CALL( SCIPflushRowExtensions(masterscip, *ssrc) );

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
   SCIP_Real  coeff_ssrc;
   SCIP_Bool  success;
   int        nmasterconsvars;
   int        i;
   int        j;

   assert(masterscip != NULL);
   assert(mapmastervarxcoeff != NULL);
   assert(k != 0);

   *rhs_ssrc = 0;
   for( i = 0; i < nselectedconss; i++ )
   {
      SCIP_CONS* mastercons = NULL;
      SCIP_VAR** masterconsvars = NULL;
      SCIP_Real* masterconscoeffs = NULL;
      SCIP_Real  rhs_mastercons;
      SCIP_Real  lhs_mastercons;

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
   SCIPdebugMessage("rhs of subset row cut: %f\n", *rhs_ssrc);

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
      SCIPdebugMessage("Process original constraint %s\n", SCIPconsGetName(origcons));
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
         pricingvar = GCGoriginalVarGetPricingVar(origconsvars[j]);
         if( pricingvar == NULL )
         {
            SCIPdebugMessage("Original vars %s did not have a corresponding pricing var!\n", SCIPvarGetName(origconsvars[j]));
            continue;
         }

         assert(pricingvar != NULL);
         coeff_pricing = SCIPhashmapGetImageReal(mappricingvarxcoeff, pricingvar);
         if( coeff_pricing == SCIP_INVALID )
         {
            SCIPhashmapSetImageReal(mappricingvarxcoeff, pricingvar, weights[i] * origconscoeffs[j]);
         }
         else
         {
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
   SCIP_HASHMAP*              mapmastervarxcoeff = NULL;  // maps master variable to its coefficient in the subset row
   SCIP_HASHMAP*              mappricingvarxcoeff = NULL; // maps pricing variable to its coefficient in its pricing constraint
   SCIP_CONS**                masterconss;         // a              : constraints in master problem
   SCIP_CONS**                originalconss;       // A              : constraints in original problem
   SCIP_Bool                  success;
   SCIP_Real*                 weights = NULL;      // w in {-1/k, 1/k}^n: vector of weights for selected constraints
   int*                       selectedconssidx = NULL;    // indices of constraints used to construct subset row (weights are non-zero)
   int                        nmasterconss;        // m              : number of constraints in master problem
   int                        npricingproblems;    // K              : number of pricing problems
   int                        i;
   int                        j;
   int                        nmastervars;

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

   /* get info of master problem */
   originalconss = GCGgetOrigMasterConss(origscip);
   masterconss = GCGgetMasterConss(origscip);
   nmasterconss = GCGgetNMasterConss(origscip);
   npricingproblems = GCGgetNPricingprobs(origscip);
   nmastervars = SCIPgetNVars(scip);

   if( sepadata->n >= nmasterconss )
   {
      SCIPdebugMessage("Not enough constraints to build subsetrow: n = %i >= nmasterconss = %i!\n", sepadata->n, nmasterconss);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( (SCIPgetCurrentNode(scip) != SCIPgetRootNode(scip) && sepadata->onlyroot) || sepadata->ngeneratedcut > 0 )
   {
      SCIPdebugMessage("only run on root node\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* alloc memory for structures used to store
    * - constraint indices
    * - weights
    * - coefficients for pricing variables
    * - coefficients for master variables */
   SCIPallocBufferArray(scip, &selectedconssidx, sepadata->n);
   SCIPallocBufferArray(scip, &weights, sepadata->n);
   SCIPhashmapCreate(&mapmastervarxcoeff, SCIPblkmem(scip), nmastervars);
   SCIPhashmapCreate(&mappricingvarxcoeff, SCIPblkmem(scip), nmastervars);

   for( i = 0; i < sepadata->maxsepacuts; i ++ )
   {
      GCG_MASTERCUTDATA*         mastercutdata = NULL;         // mastercutdata for cut
      GCG_PRICINGMODIFICATION**  pricingmodifications = NULL;  // pricing modifications associated with cut
      SCIP_ROW*                  ssrc = NULL;                  // cut
      SCIP_Real                  rhs_ssrc;                     // rhs of the cut
      int                        npricingmodifications = 0;    // number of pricing modifications associated with cut
      char                       name[SCIP_MAXSTRLEN];

      /* select n constraints and store their indices in selectedconssidx */
      SCIPdebugMessage("Select %i out of %i constraints\n", sepadata->n, nmasterconss);
      if( sepadata->strategy == 0 )
         SCIP_CALL( selectRandomRows(sepadata, nmasterconss, selectedconssidx, sepadata->n) );

      /* determine the master variables, their coefficients and rhs for subset row */
      SCIP_CALL( computeSubsetRowCoefficientsAndRHS(scip, masterconss, selectedconssidx, sepadata->n,
                                                    sepadata->k, &weights, mapmastervarxcoeff, &rhs_ssrc) );

      /* create the subsetrow cut */
      SCIP_CALL( createSubsetRowCut(scip, &ssrc, mapmastervarxcoeff, rhs_ssrc, sepa) );
      assert(ssrc != NULL);

      /* row is empty --> useless */
      if( SCIProwGetNNonz(ssrc) == 0 )
      {
         SCIPdebugMessage("created an empty row: release row\n");
         SCIPreleaseRow(scip, &ssrc);
         SCIP_CALL( SCIPhashmapRemoveAll(mapmastervarxcoeff) );
         continue;
      }

      /* determine the pricing variables and their coefficients for the pricing constraints */
      SCIP_CALL( computePricingConssCoefficients(origscip, originalconss, selectedconssidx, sepadata->n,
                                                 weights, mappricingvarxcoeff) );

      /* create the pricing modification for every pricing problem */
      for( j = 0; j < npricingproblems; j++ )
      {
         SCIP*                      pricingproblem;
         GCG_PRICINGMODIFICATION*   pricingmodification = NULL;
         SCIP_CONS**                pricingconss = NULL;
         SCIP_VAR**                 pricingvars;
         SCIP_VAR*                  coeffvar = NULL; // y
         SCIP_Real                  pricingcoeff;
         int                        npricingvars;
         int                        nvars;
         int                        l;

         pricingproblem = GCGgetPricingprob(origscip, j);
         assert(pricingproblem != NULL);
         assert(GCGisPricingprobRelevant(origscip, j));
         assert(SCIPgetObjsense(pricingproblem) == SCIP_OBJSENSE_MINIMIZE);

         npricingvars = SCIPgetNVars(pricingproblem);
         pricingvars = SCIPgetVars(pricingproblem);

         /* we add at most one constraint to each pricing problem */
         SCIP_CALL( SCIPallocBlockMemoryArray(pricingproblem, &pricingconss, 1) ); // freed via GCGpricingmodificationFree

         /* create (and capture) 'empty' pricing constraint: -inf <= ... <= 1 - EPSILON */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pp%i_cons_ssrc_%i", j, sepadata->ngeneratedcut);
         SCIP_CALL( SCIPcreateConsBasicLinear(pricingproblem, &pricingconss[0], name,
                                              0, NULL, NULL, -SCIPinfinity(pricingproblem),
                                              1 - 0.0001) ); // released via GCGpricingmodificationFree

         /* fill constraint such that -inf <= w^TAx <= 1 - EPSILON */
         for( l = 0; l < npricingvars; l++ )
         {
            assert(GCGvarIsPricing(pricingvars[l]));

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
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pp%i_y_ssrc_%i", j, sepadata->ngeneratedcut);
         SCIP_CALL( GCGcreateInferredPricingVar(pricingproblem, &coeffvar, name, 0.0,
                                                SCIPinfinity(pricingproblem), -1.0, SCIP_VARTYPE_INTEGER, j) ); // released in GCGpricingmodificationFree
         assert(coeffvar != NULL);

         /* add y to constraint such that: -inf <= w^TAx - y <= 1 - EPSILON */
         SCIP_CALL( SCIPaddCoefLinear(pricingproblem, pricingconss[0], coeffvar, -1.0) );

         /* create pricing modifications containing y as the coeffvar and a single constraint */
         SCIP_CALL( GCGpricingmodificationCreate(scip, &pricingmodification, j, coeffvar, NULL,
                                                 0, pricingconss, 1) ); // released in GCGpricingmodificationFree

         /* ensure we have enough memory for all the pricing modifications */
         if( npricingmodifications == 0 )
         {
            SCIPallocBlockMemoryArray(scip, &pricingmodifications, 1); // freed in GCGmastercutFree
         }
         else
         {
            SCIPreallocBlockMemoryArray(scip, &pricingmodifications, npricingmodifications, npricingmodifications + 1); // freed in GCGmastercutFree
         }

         pricingmodifications[npricingmodifications] = pricingmodification;
         npricingmodifications++;
      }
      SCIPdebugMessage("number of pricing mods: %i\n", npricingmodifications);

      /* create a mastercutdata containing y, ssrc and the pricing modifications */
      SCIP_CALL( GCGmastercutCreateFromRow(scip, &mastercutdata, ssrc, pricingmodifications, npricingmodifications) ); // freed in GCGmastercutFree
      sepadata->ngeneratedcut++;

      /* add the subsetrow cut to the sepa store of the master problem and the generated cuts*/
      SCIP_CALL( SCIPaddRow(scip, ssrc, FALSE, &success) );
      SCIP_CALL( GCGaddCutToGeneratedCutsSepa(scip, mastercutdata, weights, selectedconssidx, sepadata->n, sepadata->sepaidx) );

      /* empty the hashmaps containing the coeffcients for variables */
      SCIP_CALL( SCIPhashmapRemoveAll(mapmastervarxcoeff) );
      SCIP_CALL( SCIPhashmapRemoveAll(mappricingvarxcoeff) );
   }

   /* free data structures */
   SCIPhashmapFree(&mapmastervarxcoeff);
   SCIPhashmapFree(&mappricingvarxcoeff);
   SCIPfreeBufferArray(scip, &selectedconssidx);
   SCIPfreeBufferArray(scip, &weights);

   return SCIP_OKAY;
}


/*
 * Callback methods of GCG separator
 */

/** get all the cuts generated by this master separator */
static
GCG_DECL_SEPAGETCOLCOEFFICIENTS(gcgsepaGetColCoefficientSubsetrow)
{
   SCIP_Real* weights;
   SCIP_Real* mastercoeffs;
   SCIP_ROW* mastercutrow;
   GCG_MASTERCUTDATA* mastercutdata;
   int* conssindices;
   int n;
   int i;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(sepa != NULL);
   assert(sepa->separator != NULL);
   assert(gcgcol != NULL);

   mastercoeffs = GCGcolGetMastercoefs(gcgcol);
   weights = GCGsepamastercutGetWeights(cut);
   conssindices = GCGsepamastercutGetConssIndices(cut);
   n = GCGsepamastercutGetNWeights(cut);

   assert(mastercoeffs != NULL);
   assert(weights != NULL);
   assert(conssindices != NULL);

   /* use the coefficients of the master constraints to compute coefficient for cut */
   *coeff = 0.0;
   for( i = 0; i < n; i++ )
   {
      *coeff += weights[i] * mastercoeffs[conssindices[i]];
   }

#ifdef SCIP_DEBUG
   mastercutdata = GCGsepamastercutGetMastercutData(cut);
   assert(mastercutdata != NULL);
   SCIP_CALL( GCGmastercutGetRow(mastercutdata, &mastercutrow) );
   assert(mastercutrow != NULL);
   SCIPdebugMessage("column coeff for cut %s: %f\n", SCIProwGetName(mastercutrow), *coeff);
#endif
   *coeff = SCIPfeasFloor(scip, *coeff);

   return SCIP_OKAY;
}

/** method for adding new master variable to cut */
static
GCG_DECL_SEPAGETVARCOEFFICIENT(gcgsepaGetVarCoefficientSubsetrow)
{
   SCIP* origscip;
   SCIP* pricingscip;
   GCG_PRICINGMODIFICATION* pricingmod;
   SCIP_CONS** pricingconss;
   SCIP_VAR** pricingconsvars;
   SCIP_Bool success;
   SCIP_Real* pricingconscoeffs;
   int npricingconsvars;
   int j;
   int k;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(sepa != NULL);
   assert(cut != NULL);

   origscip = GCGgetOriginalprob(scip);
   pricingmod = GCGmastercutGetPricingModification(scip, cut, probnr);
   pricingscip = GCGgetPricingprob(origscip, probnr);
   pricingconss = GCGpricingmodificationGetAdditionalConss(pricingmod);

   SCIP_CALL( SCIPgetConsNVars(pricingscip, pricingconss[0], &npricingconsvars, &success) );
   assert(success);
   SCIPallocBufferArray(scip, &pricingconsvars, npricingconsvars);
   SCIPallocBufferArray(scip, &pricingconscoeffs, npricingconsvars);
   SCIP_CALL( SCIPgetConsVars(pricingscip, pricingconss[0], pricingconsvars, npricingconsvars, &success) );
   assert(success);
   SCIP_CALL( SCIPgetConsVals(pricingscip, pricingconss[0], pricingconscoeffs, npricingconsvars, &success) );
   assert(success);

   /* compute floor(w^TAx) */
   *coef = 0.0;
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

   SCIPdebugMessage("variable coeff %f\n", *coef);
   *coef = SCIPfeasFloor(scip, *coef);

   SCIPfreeBufferArray(scip, &pricingconsvars);
   SCIPfreeBufferArray(scip, &pricingconscoeffs);

   return SCIP_OKAY;

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

   /* set the objective value of each coefficient variable y to dual of the cut it is associated with */
   for( i = 0; i < npricingmodifications; i++ )
   {
      pricingblocknr = GCGpricingmodificationGetBlock(pricingmodifications[i]);
      pricingproblem = GCGgetPricingprob(origscip, pricingblocknr);
      coeffvar = GCGpricingmodificationGetCoefVar(pricingmodifications[i]);
      assert(SCIPvarGetProbindex(coeffvar) != -1);
      if( SCIPisZero(pricingproblem, dual) )
      {
         SCIP_CALL( SCIPchgVarObj(pricingproblem, coeffvar, 0.0) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarObj(pricingproblem, coeffvar, -1.0 * dual) );
         SCIPdebugMessage("%s's objective: %f\n", SCIPvarGetName(coeffvar), SCIPvarGetObj(coeffvar));
      }

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

   SCIP_CALL( SCIPaddIntParam(scip,
                              "sepa/" SEPA_NAME "/maxrounds",
      "maximal number of subsetrow separation rounds per node (-1: unlimited)",
      &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
                              "sepa/" SEPA_NAME "/maxroundsroot",
      "maximal number of subsetrow separation rounds in the root node (-1: unlimited)",
      &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
                              "sepa/" SEPA_NAME "/maxsepacuts",
      "maximal number of subsetrow cuts separated per separation round",
      &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
                              "sepa/" SEPA_NAME "/initseed",
      "initial seed used for the row selection in cut selection",
      &sepadata->initseed, FALSE, DEFAULT_INITSEED, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
                              "sepa/" SEPA_NAME "/maxsepacutsroot",
      "maximal number of zerohalf cuts separated per separation round in the root node",
      &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
                              "sepa/" SEPA_NAME "/maxcutcands",
      "maximal number of zerohalf cuts considered per separation round",
      &sepadata->maxcutcands, FALSE, DEFAULT_MAXCUTCANDS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(GCGmasterGetOrigprob(scip), "sepa/" SEPA_NAME "/onlyroot", "apply subsetrow separator only on root",
      &(sepadata->onlyroot), FALSE, DEFAULT_ONLYROOT, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
                              "sepa/" SEPA_NAME "/strategy",
      "RANDOM (0), ZERO-HALF (1)",
      &sepadata->strategy, FALSE, DEFAULT_STRATEGY, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
                              "sepa/" SEPA_NAME "/n",
      "number of rows used to create a new subset row cut ",
      &sepadata->n, FALSE, DEFAULT_N, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
                              "sepa/" SEPA_NAME "/k",
      "weight used to create new subset row cut",
      &sepadata->k, FALSE, DEFAULT_K, 1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}