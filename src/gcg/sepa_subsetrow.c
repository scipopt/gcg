/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file    sepa_subsetrow.c
 * @brief   subset row separator for master problem
 * @author  Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#include <assert.h>
#include <scip/cutsel_hybrid.h>
#include <scip/cutsel_dynamic.h>
#include <scip/cons_linear.h>

#include "gcg.h"
#include "mastercutdata.h"
#include "mastersepacut.h"
#include "pricer_gcg.h"
#include "pub_gcgcol.h"
#include "relax_gcg.h"
#include "scip_misc.h"
#include "sepa_subsetrow.h"
#include "struct_gcgcol.h"
#include "struct_sepagcg.h"
#include "type_sepagcg.h"
#include "type_mastersepacut.h"

#define SEPA_NAME           "subsetrow"
#define SEPA_DESC "subsetrow separator"
#define SEPA_PRIORITY               500
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define STARTMAXCUTS                 50
#define DEFAULT_RANDSEED             71
#define DEFAULT_MAXROUNDS             1 /**< maximal number of subset row separation rounds per non-root node */
#define DEFAULT_MAXROUNDSROOT         1 /**< maximal number of subset row separation calls in the root node */
#define DEFAULT_MAXSEPACUTS         200 /**< maximal number of subset row cuts separated per call in non-root nodes */
#define DEFAULT_MAXSEPACUTSROOT    2000 /**< maximal number of subset row cuts separated per call in root node */
#define DEFAULT_MAXCUTCANDS        4000 /**< maximal number of subset row cuts in total */
#define DEFAULT_ONLYROOT          FALSE /**< only apply separator in root node */
#define DEFAULT_STRATEGY              0 /**< strategy which is used to determine which rows to consider for cut computation */
#define DEFAULT_N                     3 /**< number of rows used to create a new cut */
#define DEFAULT_K                     2 /**< inverse of weight used for cut generation */

/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_RANDNUMGEN*        randnumgen;          /**< random number generator (for strategy RANDOM) */
   SCIP_Bool               enable;              /**< is this separator enabled? */
   SCIP_Bool               onlyroot;            /**< indicates if separator should only be applied at root node */
   int                     ngeneratedcut;       /**< counts the total number of cuts generated */
   int                     maxrounds;           /**< maximal number of separation calls per non-root node */
   int                     maxroundsroot;       /**< maximal number of separation calls for root node */
   int                     maxsepacutsroot;     /**< number of cuts generated per separation call of root node */
   int                     maxsepacuts;         /**< number of cuts generated per separation call at non-root node */
   int                     maxcutcands;         /**< maximal number of cuts generated in total */
   int                     strategy;            /**< RANDOM (0), */
   int                     n;                   /**< k > 0          : defines the possible weights 1/k */
   int                     k;                   /**< n = |S| > 0    : number of constraints used to construct subset row */
   GCG_SEPA*               sepa;                /**< gcg master separator instance */
};


/*
 * Callback methods of separator
 */

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

   SCIPdebugMessage("free separator data for subset row separator\n");
   sepadata = SCIPsepaGetData(sepa);
   SCIPfreeRandom(scip, &(sepadata->randnumgen));
   SCIPfreeBlockMemory(scip, &sepadata);

   return SCIP_OKAY;
}

/**< create and add a subset row cut the generated cuts */
static
SCIP_RETCODE addSubsetRowCutToGeneratedCuts(
   SCIP*                masterscip,       /**< SCIP data structure */
   GCG_MASTERCUTDATA*   mastercutdata,    /**< master cut data */
   SCIP_Real*           weights,          /**< weights used to create the cut */
   int*                 conssindices,     /**< indices of constraints used to create the cut */
   int                  n,                /**< number of constraints used to create the cut */
   GCG_SEPA*            sepa              /**< separator which generated the cut */
)
{
   GCG_MASTERSEPACUT* mastersepacut;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(mastercutdata != NULL);

   /* create a subset row cut */
   SCIP_CALL( GCGcreateSubsetRowCut(masterscip, &mastersepacut, sepa, mastercutdata, NULL, weights, conssindices, n) );
   assert(mastersepacut != NULL);

   /* register it with the event handler managing active master separator cuts */
   SCIP_CALL( GCGaddCutToGeneratedCuts(masterscip, mastersepacut) );

   return SCIP_OKAY;
}


/** randomly selects n different constraints from the master problem */
static
SCIP_RETCODE selectRandomRows(
   SCIP_SEPADATA* sepadata,                   /**< separator data */
   int            nmasterconss,               /**< number of constraints in the master problem */
   int*           selectedmasterconssidx,     /**< pointer to store the indices of the selected constraints */
   int            n                           /**< number of constraints to be selected */
   )
{
   int i;
   int j;

   assert(n > 0 && n < nmasterconss);

   /* randomly select n indices out of [0, ..., nmasterconss - 1] */
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
#ifdef SCIP_DEBUG
   for( i = 0; i < n; i++ )
   {
      SCIPdebugMessage("select index %i\n", selectedmasterconssidx[i]);
      assert(0 <= selectedmasterconssidx[i] && selectedmasterconssidx[i] < nmasterconss);
   }
#endif

   return SCIP_OKAY;
}

/** create a new row for master problem and fill it with the variables (+ their coefficients) and right-hand side*/
static
SCIP_RETCODE createSubsetRowCut_alt(
   SCIP*          masterscip,                   /**< SCIP data structure (master problem) */
   SCIP_ROW**     ssrc,                         /**< pointer to store subset row cut */
   SCIP_Real*     ssrccoeffs,                   /**< coefficients of variables in subset row cut */
   int*           nonzerossrccoeffsindices,     /**< indices of variables with non-zero coefficient in subset row cut */
   int            nnonzerossrccoeffsindices,    /**< number of variables with non-zero coefficient in subset row cut */
   SCIP_Real      rhs_ssrc,                     /**< right hand side of subset row cut */
   SCIP_SEPA*     sepa                          /**< separator which creates subset row cut */
)
{
   SCIP_SEPADATA* sepadata;
   SCIP_VAR**     mastervars;
   char           name[SCIP_MAXSTRLEN]; // name of the ssrc
   int            i;
   int            nmastervar;

   assert(GCGisMaster(masterscip));
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   nmastervar = SCIPgetNVars(masterscip);
   mastervars = SCIPgetVars(masterscip);

   SCIPdebugMessage("create subset row cut\n");
   /* create 'empty' subset row cut of form -inf <= ... <= rhs_ssrc
    * - local, non-removable, modifiable */
   rhs_ssrc = rhs_ssrc / sepadata->k;
   rhs_ssrc = SCIPfeasFloor(masterscip, rhs_ssrc);
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ssrc_%i", sepadata->ngeneratedcut);
   SCIP_CALL( SCIPcreateEmptyRowSepa(masterscip, &(*ssrc), sepa, name, -SCIPinfinity(masterscip), rhs_ssrc, TRUE, TRUE, TRUE) );
   assert(ssrc != NULL);
   assert(*ssrc != NULL);

   /* fill the row with master variables and their coefficients */
   SCIP_CALL( SCIPcacheRowExtensions(masterscip, *ssrc) );
   for( i = 0; i < nnonzerossrccoeffsindices; i++ )
   {
      SCIP_Real varcoeff;
      SCIP_VAR* mastervar;
      int       varindex;

      varindex = nonzerossrccoeffsindices[i];
      assert(varindex <= nmastervar);

      varcoeff = ssrccoeffs[varindex] / sepadata->k;
      varcoeff = SCIPfeasFloor(masterscip, varcoeff);
      ssrccoeffs[varindex] = 0.0; // have to reset it to zero, re-allocating a clean buffer array does not set them back to zero????? vllt clear benutzen??

      if( varcoeff == 0.0 )
         continue;

      mastervar = mastervars[varindex];
      assert(SCIPvarGetProbindex(mastervar) == varindex);
      SCIPdebugMessage("add var %s with coefficient %f to subset row cut %s\n", SCIPvarGetName(mastervar), varcoeff, name);
      SCIP_CALL( SCIPaddVarToRow(masterscip, *ssrc, mastervar, varcoeff) );
   }
   SCIP_CALL( SCIPflushRowExtensions(masterscip, *ssrc) );

   return SCIP_OKAY;
}

/** computes the rhs (w^Tb) and the coefficient for each variable (w^Ta_p) in the cut (still non-rounded) */
static
SCIP_RETCODE computeSubsetRowCoefficientsAndRHS_alt(
   SCIP*          masterscip,                 /**< SCIP data structure of master problem */
   SCIP_CONS**    masterconss,                /**< constraints of the master problem */
   int*           selectedconssidx,           /**< indices of the selected constraints */
   int            nselectedconss,             /**< number of selected constraints */
   SCIP_Real**    weights,                    /**< pointer to store weights of the selected constraints */
   SCIP_Real*     rhs_ssrc,                   /**< pointer to store rhs of subset row cut */
   SCIP_Real**    mastervarcoeffs,            /**< pointer to store subset row cut variable coefficients */
   int**          nonzeromastervardindices,   /**< pointer to store indices of master variables with non-zero coefficients */
   int*           nnonzeromastervardindices   /**< pointer to store number of master variables with non-zero coefficients */
)
{
   SCIP_Bool  success;
   int        nmastervars;
   int        nmasterconsvars;
   int        i;
   int        j;

   assert(masterscip != NULL);

   nmastervars = SCIPgetNVars(masterscip);
   *nnonzeromastervardindices = 0;
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
         (*weights)[i] = 1.0;
         (*rhs_ssrc) += (*weights)[i] * rhs_mastercons;
         SCIPdebugMessage("master constraint %s (ax <= %f) with weight %f\n", SCIPconsGetName(mastercons), rhs_mastercons, (*weights)[i]);
      }
      else if( SCIPisInfinity(masterscip, rhs_mastercons ) ) // lhs != -inf & rhs = inf --> constraint has form: b <= ax
      {
         (*weights)[i] = -1.0;
         (*rhs_ssrc) += (*weights)[i] * lhs_mastercons;
         SCIPdebugMessage("master constraint %s (%f <= ax) with weight %f\n", SCIPconsGetName(mastercons), lhs_mastercons, (*weights)[i]);
      }
      else // lhs != -inf & rhs != inf --> constraint has form: b_l <= ax <= b_r
      {
         /* we use the right side: ax <= b_r */
         (*weights)[i] = 1.0;
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
         int varindex;

         assert(GCGvarIsMaster(masterconsvars[j]));

         varindex = SCIPvarGetProbindex(masterconsvars[j]);
         assert(varindex <= nmastervars);

         if( masterconscoeffs[j] != 0.0 )
         {
            SCIP_Real varcoeff;

            varcoeff = (*mastervarcoeffs)[varindex];

            /* new non-zero variable coefficient */
            if( varcoeff == 0.0 )
            {
               (*nonzeromastervardindices)[*nnonzeromastervardindices] = varindex;
               (*nnonzeromastervardindices)++;
            }

            varcoeff += (*weights)[i] * masterconscoeffs[j];
            (*mastervarcoeffs)[varindex] = varcoeff;
            SCIPdebugMessage("increase %s-coefficient to %f\n", SCIPvarGetName(masterconsvars[j]), varcoeff);
         }
      }

      SCIPfreeBufferArray(masterscip, &masterconsvars);
      SCIPfreeBufferArray(masterscip, &masterconscoeffs);
   }

   SCIPdebugMessage("right-hand side of subset row cut: %f\n", *rhs_ssrc);

   return SCIP_OKAY;
}

/** computes the (non-rounded) coefficients for the pricing variables used in the pricing constraints */
static
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
   SCIP_VAR*  pricingvar;
   SCIP_Real* origconscoeffs;
   SCIP_Real  coeff_pricing;
   int        norigconsvars;
   int        i;
   int        j;

   assert(GCGisOriginal(origscip));
   assert(nselectedconss > 0);

   SCIPdebugMessage("compute the coefficients of the pricing variables\n");
   /* - compute w_1 * A_1j + ... + w_m * A_mj for each pricing variable x_j */
   for( i = 0; i < nselectedconss; i++ )
   {
      /* get all variables and their corresponding coefficients in the original constraint */
      origcons = originalconss[selectedconssidx[i]];
      SCIPdebugMessage("Process original constraint %s\n", SCIPconsGetName(origcons));
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
            assert(!GCGisPricingprobRelevant(origscip, GCGvarGetBlock(origconsvars[j])));
            continue;
         }
         assert(GCGisPricingprobRelevant(origscip, GCGvarGetBlock(origconsvars[j])));
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
   SCIP*             origscip;                     // original problem
   SCIP_SEPADATA*    sepadata;
   SCIP_HASHMAP*     mappricingvarxcoeff = NULL;   // maps pricing variable to its coefficient in its pricing constraint
   SCIP_CONS**       masterconss;                  // a              : constraints in master problem
   SCIP_CONS**       originalconss;                // A              : constraints in original problem
   SCIP_Bool         success;
   SCIP_Real*        weights = NULL;               // w in {-1/k, 1/k}^n: vector of weights for selected constraints
   int*              selectedconssidx = NULL;      // indices of constraints used to construct subset row (weights are non-zero)
   int               nmasterconss;                 // m              : number of constraints in master problem
   int               npricingproblems;             // K              : number of pricing problems
   int               nmastervars;
   int               maxcuts;                      // number of cuts to be generated in a separation round
   int               i;
   int               j;
   SCIP_ROW**        cuts = NULL;
   int               ncuts = 0;

   assert(scip != NULL);
   assert(result != NULL);
   assert(sepa != NULL);
   assert(GCGisMaster(scip));

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIP_Bool isroot;
   isroot = SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip);

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

   if( SCIPsepaGetNCallsAtNode(sepa) >= sepadata->maxrounds )
   {
      SCIPdebugMessage("exceeded max rounds for this node.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( (isroot && sepadata->onlyroot) || (isroot && !allowlocal))
   {
      SCIPdebugMessage("subset row separator is only configured to run on root node.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( SCIPnodeGetType(SCIPgetCurrentNode(scip)) == SCIP_NODETYPE_PROBINGNODE || SCIPnodeGetType(SCIPgetCurrentNode(scip)) == SCIP_NODETYPE_REFOCUSNODE )
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

   /* if( GCGgetNActiveCuts(scip) >= 10 )
   {
      SCIPdebugMessage("applied cuts already exceeded limit.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }*/

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
      SCIPdebugMessage("not enough constraints to build subset row cut: n = %i >= number of master constraints = %i!\n", sepadata->n, nmasterconss);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* determine the number of cuts to generate based on node type */
   if( isroot )
   {
      maxcuts = sepadata->maxsepacutsroot;
   }
   else
   {
      maxcuts = sepadata->maxsepacuts;
   }

   SCIPallocBufferArray(scip, &selectedconssidx, sepadata->n);
   SCIPallocBufferArray(scip, &weights, sepadata->n);
   SCIPhashmapCreate(&mappricingvarxcoeff, SCIPblkmem(scip), nmastervars);
   SCIPallocBufferArray(scip, &cuts, maxcuts);


   for( i = 0; i < maxcuts; i ++ )
   {
      GCG_MASTERCUTDATA*         mastercutdata = NULL;         // master cut data for subset row cut
      GCG_PRICINGMODIFICATION*   pricingmodifications = NULL;  // pricing modifications associated with cut
      SCIP_ROW*                  ssrc = NULL;                  // cut
      SCIP_Real                  rhs_ssrc;                     // right-hand side of the cut
      int                        npricingmodifications = 0;    // number of pricing modifications associated with cut
      char                       name[SCIP_MAXSTRLEN];
      SCIP_Real*                 ssrccoeffs;                   // contains (non-rounded) coefficients for cut
      int*                       nonzerossrccoeffsindices;     // contains indices of variables with non-zero coefficients in cut
      int                        nnonzerossrccoeffsindices;    // number of variables with non-zero coefficients in cut

      /* select n constraints to build subset row cut */
      SCIPdebugMessage("cut % i: select %i out of %i constraints\n", i, sepadata->n, nmasterconss);
      if( sepadata->strategy == 0 )
         SCIP_CALL( selectRandomRows(sepadata, nmasterconss, selectedconssidx, sepadata->n) );

      /* determine the master variables, their coefficients and rhs for subset row (non-rounded) */
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &ssrccoeffs, nmastervars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nonzerossrccoeffsindices, nmastervars) );
      nnonzerossrccoeffsindices = 0;
      SCIP_CALL( computeSubsetRowCoefficientsAndRHS_alt(scip, masterconss, selectedconssidx, sepadata->n,
                                                        &weights, &rhs_ssrc, &ssrccoeffs,
                                                       &nonzerossrccoeffsindices, &nnonzerossrccoeffsindices) );

      /* create the subset row cut */
      SCIP_CALL( createSubsetRowCut_alt(scip, &ssrc, ssrccoeffs, nonzerossrccoeffsindices,
                                        nnonzerossrccoeffsindices, rhs_ssrc, sepa) );
      assert(ssrc != NULL);
      SCIPfreeCleanBufferArrayNull(scip, &ssrccoeffs);
      SCIPfreeBufferArrayNull(scip, &nonzerossrccoeffsindices);

      /* row is empty --> useless */
      if( SCIProwGetNNonz(ssrc) == 0 )
      {
         SCIPdebugMessage("created an empty row: release row\n");
         SCIPreleaseRow(scip, &ssrc);
         SCIP_CALL( SCIPhashmapRemoveAll(mappricingvarxcoeff) );
         continue;
      }

      //SCIPprintRow(scip, ssrc, NULL);

      /* determine the pricing variables and their coefficients for the pricing constraints */
      SCIP_CALL( computePricingConssCoefficients(origscip, originalconss, selectedconssidx, sepadata->n,
                                                 weights, mappricingvarxcoeff) );

      /* create the pricing modification for every (relevant) pricing problem */
      for( j = 0; j < npricingproblems; j++ )
      {
         SCIP*                      pricingproblem;
         GCG_PRICINGMODIFICATION    pricingmodification;
         SCIP_CONS**                pricingconss = NULL;
         SCIP_VAR**                 pricingvars;
         SCIP_VAR*                  coeffvar = NULL; // y
         SCIP_Real                  pricingcoeff;
         int                        npricingvars;
         int                        nvars;
         int                        l;

         /* in case of aggregated pricing problems, we skip the non-representative ones */
         pricingproblem = GCGgetPricingprob(origscip, j);
         if( pricingproblem == NULL || !GCGisPricingprobRelevant(origscip, j) )
            continue;

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
                                              1 - 0.00001) ); // released via GCGpricingmodificationFree

         /* fill constraint such that -inf <= w^TAx <= 1 - EPSILON */
         for( l = 0; l < npricingvars; l++ )
         {
            assert(GCGvarIsPricing(pricingvars[l]));
            pricingcoeff = SCIPhashmapGetImageReal(mappricingvarxcoeff, pricingvars[l]);
            if( pricingcoeff == SCIP_INVALID )
               continue;

            pricingcoeff = pricingcoeff / sepadata->k;
            if( !SCIPisZero(pricingproblem, pricingcoeff) ) // @todo: == 0.0??
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

         /* create (and capture) y: -inf <= y <= inf (integer) */
         SCIPdebugMessage("create new (inferred) pricing variable y for pricing problem %i\n", j);
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pp%i_y_ssrc_%i", j, sepadata->ngeneratedcut);
         SCIP_CALL( GCGcreateInferredPricingVar(pricingproblem, &coeffvar, name, -SCIPinfinity(pricingproblem),
                                                SCIPinfinity(pricingproblem), -1.0, SCIP_VARTYPE_INTEGER, j) ); // released in GCGpricingmodificationFree
         assert(coeffvar != NULL);

         /* add y to constraint such that: -inf <= w^TAx - y <= 1 - EPSILON  <=> w^TAx - 1 + EPSILON <= y */
         SCIP_CALL( SCIPaddCoefLinear(pricingproblem, pricingconss[0], coeffvar, -1.0) );
         SCIPdebugPrintCons(pricingproblem, pricingconss[0], NULL);

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

      /* create master cut data containing y, ssrc and the pricing modifications */
      SCIP_CALL( GCGmastercutCreateFromRow(scip, &mastercutdata, ssrc, pricingmodifications, npricingmodifications) ); // freed in GCGmastercutFree
      sepadata->ngeneratedcut++;

      /* add the subset row cut to the sepastore of the master problem and the generated cuts */

      SCIP_CALL( addSubsetRowCutToGeneratedCuts(scip, mastercutdata, weights, selectedconssidx, sepadata->n, sepadata->sepa) );
      cuts[ncuts] = ssrc;
      ncuts++;
      /* empty the hashmaps containing the coefficients for pricing variables */
      SCIP_CALL( SCIPhashmapRemoveAll(mappricingvarxcoeff) );
   }
   if( ncuts > 0 )
   {
      int nselectedcuts;
      SCIP_CALL( SCIPselectCutsHybrid(scip, cuts, NULL, sepadata->randnumgen, 1.0, 0.5,
                                      0.1, 0.1, 0.0, 1.0, 0.0, 0.0,
                                      ncuts, 0, 5, &nselectedcuts) );


      for( i = 0; i < ncuts; i++ )
      {
         if( i < nselectedcuts )
         {
            SCIPinfoMessage(scip, NULL, "selected row: %s\n", SCIProwGetName(cuts[i]));
            SCIP_CALL( SCIPaddRow(scip, cuts[i], FALSE, &success) );
         }
         /*else
         {
            SCIPreleaseRow(scip, &cuts[i]);
         }*/

      }
   }

   /* free data structures */
   SCIPhashmapFree(&mappricingvarxcoeff);
   SCIPfreeBufferArray(scip, &selectedconssidx);
   SCIPfreeBufferArray(scip, &weights);
   SCIPfreeBufferArray(scip, &cuts);

   return SCIP_OKAY;
}


/*
 * Callback methods of GCG separator
 */

/** compute cut coefficient for a column */
static
GCG_DECL_SEPAGETCOLCOEFFICIENTS(gcgsepaGetColCoefficientSubsetrow)
{
   SCIP_SEPADATA* sepadata;
   SCIP_Real*     weights;
   SCIP_Real*     mastercoeffs;
   int*           conssindices;
   int            n;
   int            i;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(sepa != NULL);
   assert(sepa->separator != NULL);
   assert(gcgcol != NULL);
   assert(GCGcolGetInitializedCoefs(gcgcol));

   mastercoeffs = GCGcolGetMastercoefs(gcgcol);
   weights = GCGsubsetrowCutGetWeights(cut);
   conssindices = GCGsubsetrowCutGetConssIndices(cut);
   n = GCGsubsetrowCutGetNWeights(cut);
   sepadata = SCIPsepaGetData(sepa->separator);

   assert(mastercoeffs != NULL);
   assert(weights != NULL);
   assert(conssindices != NULL);
   assert(sepadata != NULL);

   /* use the coefficients of the master constraints to compute coefficient for cut */
   *coeff = 0.0;
   for( i = 0; i < n; i++ )
   {
      *coeff += weights[i] * mastercoeffs[conssindices[i]] / sepadata->k;
   }

   *coeff = SCIPfeasFloor(scip, *coeff);
   SCIPdebugMessage("column coefficient: %f\n", *coeff);
   return SCIP_OKAY;
}

/** compute cut coefficient for master variable */
static
GCG_DECL_SEPAGETVARCOEFFICIENT(gcgsepaGetVarCoefficientSubsetrow)
{
   SCIP*                      origscip;
   SCIP*                      pricingscip;
   GCG_PRICINGMODIFICATION*   pricingmod;
   GCG_MASTERCUTDATA*         mastercutdata;
   SCIP_CONS**                pricingconss;
   SCIP_VAR**                 pricingconsvars;
   SCIP_Bool                  success;
   SCIP_Real*                 pricingconscoeffs;
   SCIP_Real*                 pricingvals = NULL;
   int*                       nonzeros;
   int                        nnonzeros;
   int                        npricingconsvars;
   int                        npricingvars;
   int                        varindex;
   int                        i;

   /* @todo: other (more efficient) way to compute coefficient ???
    * pricing vars are sorted: maybe use SCIPsortedvecFindPtr ????*/
   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(sepa != NULL);
   assert(cut != NULL);

   origscip = GCGgetOriginalprob(scip);
   assert(origscip != NULL);
   *coef = 0.0;

   if( nvars == 0 )
   {
      return SCIP_OKAY;
   }

   mastercutdata = GCGmastersepacutGetMasterCutData(cut);
   assert(mastercutdata != NULL);
   pricingmod = GCGmastercutGetPricingModification(scip, mastercutdata, probnr);

   /* no pricing modification for this problem */
   if( pricingmod == NULL )
   {
      SCIPdebugMessage("no pricing modification for pp%i --> variable coefficient 0\n", probnr);
      *coef = 0.0;
      return SCIP_OKAY;
   }

   /* get the pricing modification of the problem which generated this master variable */
   pricingscip = GCGgetPricingprob(origscip, probnr);
   pricingconss = GCGpricingmodificationGetAdditionalConss(pricingmod);

   /* transfer the values of the given variables to the position of the array which corresponds to their variable index */
   npricingvars = SCIPgetNVars(pricingscip); //SCIPgetNOrigVars(pricingscip);
   nnonzeros = 0;
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &pricingvals, npricingvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonzeros, nvars) );
   for( i = 0; i < nvars; i++ )
   {
      varindex = SCIPvarGetProbindex(vars[i]);
      assert(varindex <= npricingvars);
      pricingvals[varindex] = vals[i];
      nonzeros[nnonzeros] = varindex;
      nnonzeros++;
      //SCIPinfoMessage(scip, NULL, "%s: %f, ", SCIPvarGetName(vars[i]), vals[i]);
      SCIPdebugMessage("%s: %f, ", SCIPvarGetName(vars[i]), vals[i]);
   }
   SCIPdebugMessage("\n");
   //SCIPinfoMessage(scip, NULL, "\n");


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
      //SCIPinfoMessage(scip, NULL, "%f * %s", pricingconscoeffs[i], SCIPvarGetName(pricingconsvars[i]));
      SCIPdebugMessage("%f * %s", pricingconscoeffs[i], SCIPvarGetName(pricingconsvars[i]));
      if( GCGvarIsInferredPricing(pricingconsvars[i]) )
         continue;

      varindex = SCIPvarGetProbindex(pricingconsvars[i]);
      assert(varindex <= npricingvars);
      //SCIPinfoMessage(scip, NULL, " (%f) + ", pricingvals[varindex]);
      SCIPdebugMessage(" (%f) + ", pricingvals[varindex]);
      *coef += pricingconscoeffs[i] * pricingvals[varindex];
   }
   SCIPdebugMessage("\n");
   ///SCIPinfoMessage(scip, NULL, "\n");

   /* reset all the non-zero entries back to zero
    * (otherwise they remain there and interfere with the computation of other coefficients) */
   for( i = 0; i < nnonzeros; i++ )
   {
      varindex = nonzeros[i];
      pricingvals[varindex] = 0.0;
   }

   /* finally, we round down w^TAx */
   SCIPdebugMessage("variable coefficient %f\n", *coef);
   //SCIPinfoMessage(scip, NULL, "variable coefficient %f\n", *coef);
   *coef = SCIPfeasFloor(scip, *coef);

   SCIPfreeCleanBufferArrayNull(scip, &pricingvals);
   SCIPfreeBufferArray(scip, &nonzeros);
   SCIPfreeBufferArray(scip, &pricingconsvars);
   SCIPfreeBufferArray(scip, &pricingconscoeffs);

   return SCIP_OKAY;
}

/** modifies the objective values of the pricing variables affected by the master cut */
static
GCG_DECL_SEPASETOBJECTIVE(gcgsepaSetObjectiveSubsetrow)
{
   SCIP_ROW*                  row = NULL;
   SCIP*                      pricingproblem;
   SCIP*                      origscip;
   SCIP_VAR*                  coeffvar;
   GCG_PRICINGMODIFICATION*   pricingmodifications;
   GCG_MASTERCUTDATA*         mastercutdata;
   int                        npricingmodifications;
   int                        pricingblocknr;
   int                        i;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(sepa != NULL);
   assert(cut != NULL);

   origscip = GCGgetOriginalprob(scip);
   assert(origscip != NULL);

   mastercutdata = GCGmastersepacutGetMasterCutData(cut);
   npricingmodifications = GCGmastercutGetNPricingModifications(mastercutdata);
   pricingmodifications = GCGmastercutGetPricingModifications(mastercutdata);

   SCIP_CALL( GCGmastercutGetRow(mastercutdata, &row) );

   /* set the objective value of each coefficient variable y to -dual of the cut it is associated with */
   for( i = 0; i < npricingmodifications; i++ )
   {
      pricingblocknr = GCGpricingmodificationGetBlock(&pricingmodifications[i]);
      pricingproblem = GCGgetPricingprob(origscip, pricingblocknr);
      coeffvar = GCGpricingmodificationGetCoefVar(&pricingmodifications[i]);

      if( dual >= 0.0 ) // @todo: theoretically, dual should always be non-positive: 'correct' it to zero
      {
         SCIP_CALL( SCIPchgVarObj(pricingproblem, coeffvar, 0.0) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarObj(pricingproblem, coeffvar, -1.0 * dual) );
      }
   }

   return SCIP_OKAY;
}

/** modifies outdated column to respect cut */
static
GCG_DECL_SEPAADJUSTCOL(gcgsepaAdjustCol)
{
   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(sepa != NULL);
   assert(cut != NULL);

   GCG_PRICINGMODIFICATION* pricemod;
   GCG_MASTERCUTDATA* mastercutdata;

   mastercutdata = GCGmastersepacutGetMasterCutData(cut);
   assert(mastercutdata != NULL);

   if( !GCGmastercutIsActive(mastercutdata) )
      return SCIP_OKAY;

   pricemod = GCGmastercutGetPricingModification(scip, mastercutdata, GCGcolGetProbNr(*gcgcol));
   if( pricemod != NULL )
   {
      SCIP_VAR* coefvar;
      SCIP_Real coefvarval;

      coefvar = GCGpricingmodificationGetCoefVar(pricemod);
      assert(coefvar != NULL);

      if( SCIPvarGetIndex(coefvar) == -1 )
         return SCIP_OKAY;

      /* we compute the value of y */
      if( GCGcolGetInitializedCoefs(*gcgcol) )
      {
         SCIP_CALL( gcgsepaGetColCoefficientSubsetrow(scip, sepa, cut, *gcgcol, &coefvarval) );
      }
      else
      {
         SCIP_CALL( gcgsepaGetVarCoefficientSubsetrow(scip, sepa, cut, (*gcgcol)->vars, (*gcgcol)->vals, (*gcgcol)->nvars, (*gcgcol)->probnr, &coefvarval) );
      }

      if( !SCIPisZero((*gcgcol)->pricingprob, coefvarval) )
      {
         /* variable y should not yet be in column, so we append it to column */
         if( (*gcgcol)->maxvars < (*gcgcol)->nvars + 1 )
         {
            int newmaxvars = SCIPcalcMemGrowSize((*gcgcol)->pricingprob, (*gcgcol)->nvars + 1);
            SCIP_CALL( SCIPreallocBlockMemoryArray((*gcgcol)->pricingprob, &((*gcgcol)->vars), (*gcgcol)->maxvars,
                                                   newmaxvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray((*gcgcol)->pricingprob, &((*gcgcol)->vals), (*gcgcol)->maxvars,
                                                   newmaxvars) );
            (*gcgcol)->maxvars = newmaxvars;
         }
         (*gcgcol)->vars[(*gcgcol)->nvars] = coefvar;
         (*gcgcol)->vals[(*gcgcol)->nvars] = coefvarval;
         // ensure that array stays sorted ????
         //assert( SCIPvarCompare((*gcgcol)->vars[(*gcgcol)->nvars-1], (*gcgcol)->vars[(*gcgcol)->nvars]) != 0 );

         /* we capture inferred vars */
         SCIPcaptureVar((*gcgcol)->pricingprob, (*gcgcol)->vars[(*gcgcol)->nvars]);
         ((*gcgcol)->nvars)++;
      }
   }
   return SCIP_OKAY;
}

/** initialization method of separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitSubsetrow)
{
   SCIP*          origscip;
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);
   assert(sepa != NULL);
   assert(GCGisMaster(scip));

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* creates the subset row gcg separator and includes in the relaxator data of the original problem */
   SCIP_CALL( GCGrelaxIncludeSeparator(origscip, sepa, gcgsepaGetVarCoefficientSubsetrow,
                                                gcgsepaGetColCoefficientSubsetrow, gcgsepaSetObjectiveSubsetrow, gcgsepaAdjustCol) );
   sepadata->sepa = GCGrelaxGetSeparator(scip, SEPA_NAME);
   assert(sepadata->sepa != NULL);
   sepadata->ngeneratedcut = 0;
   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the scip separator of the subset row separator and includes it in master SCIP*/
SCIP_RETCODE SCIPincludeSepaSubsetrow(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;
   SCIP* origscip;

   /* create subsetrow separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   sepa = NULL;
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

   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

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

   SCIP_CALL( SCIPaddIntParam(origscip, "sepa/" SEPA_NAME "/k", "weight used to create new subset row cut",
      &(sepadata->k), FALSE, DEFAULT_K, 1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
