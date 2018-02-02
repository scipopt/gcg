/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file   class_stabilization.cpp
 * @brief  class with functions for dual variable smoothing
 * @author Martin Bergner
 * @author Jonas Witt
 * @author Michael Bastubbe
 *
 * This is an implementation of dynamic alpha-schedule (based on subgradient information) stabilization based on the paper
 *
 * Pessoa, A., Sadykov, R., Uchoa, E., & Vanderbeck, F. (2013). In-Out Separation and Column Generation
 * Stabilization by Dual Price Smoothing. In Experimental Algorithms (pp. 354-365). Springer Berlin Heidelberg.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_stabilization.h"
#include "pricer_gcg.h"
#include "gcg.h"
#include "sepa_master.h"
#include "objscip/objscip.h"
#include "scip/cons_linear.h"

namespace gcg {

Stabilization::Stabilization(
   SCIP* scip,
   PricingType* pricingtype_
   ) :scip_(scip), stabcenterconss((SCIP_Real*) NULL), stabcenterconsssize(0), nstabcenterconss(0),
      stabcentercuts((SCIP_Real*) NULL), stabcentercutssize(0), nstabcentercuts(0),
      stabcenterlinkingconss((SCIP_Real*) NULL), nstabcenterlinkingconss(0),
      stabcenterconv((SCIP_Real*) NULL), nstabcenterconv(0),
      pricingtype(pricingtype_), alpha(0.8), alphabar(0.8), nodenr(-1), k(0), t(0), hasstabilitycenter(FALSE),stabcenterbound(-SCIPinfinity(scip)),
         inmispricingschedule(FALSE)
{

}

Stabilization::~Stabilization()
{
   SCIPfreeBlockMemoryArrayNull(scip_, &stabcenterconss, stabcenterconsssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(scip_, &stabcentercuts, stabcentercutssize); /*lint !e64*/
   SCIPfreeMemoryArrayNull(scip_, &stabcenterlinkingconss); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(scip_, &stabcenterconv, nstabcenterconv); /*lint !e64*/
   scip_ = (SCIP*) NULL;
   stabcenterconss = (SCIP_Real*) NULL;
   stabcentercuts = (SCIP_Real*) NULL;
   stabcenterlinkingconss = (SCIP_Real*) NULL;
   stabcenterconv = (SCIP_Real*) NULL;
   pricingtype = (PricingType*) NULL;
   nodenr = -1;


}

SCIP_RETCODE Stabilization::updateStabcenterconss()
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);
   int nconss = GCGgetNMasterConss(origprob);

   if( nconss == nstabcenterconss )
   {
      return SCIP_OKAY;
   }

   if( nconss > stabcenterconsssize )
   {
      int oldsize = stabcenterconsssize;
      stabcenterconsssize = SCIPcalcMemGrowSize(scip_, nconss);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip_, &stabcenterconss, oldsize, stabcenterconsssize) );
   }
   assert(stabcenterconss != NULL);
   BMSclearMemoryArray(&stabcenterconss[nstabcenterconss], (size_t)nconss-nstabcenterconss); /*lint !e866*/

   nstabcenterconss = nconss;

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::updateStabcentercuts()
{
   int ncuts = GCGsepaGetNCuts(scip_);

   if( ncuts == nstabcentercuts )
   {
      return SCIP_OKAY;
   }

   if( ncuts > stabcentercutssize )
   {
      int oldsize = stabcentercutssize;
      stabcentercutssize = SCIPcalcMemGrowSize(scip_, ncuts);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip_, &stabcentercuts, oldsize, stabcentercutssize) );
   }
   assert(stabcentercuts != NULL);
   BMSclearMemoryArray(&stabcentercuts[nstabcentercuts], (size_t)ncuts-nstabcentercuts); /*lint !e866*/

   nstabcentercuts = ncuts;

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::setNLinkingconss(
      int nlinkingconssnew
      )
{

   SCIPfreeMemoryArrayNull(scip_, &stabcenterlinkingconss); /*lint !e64*/
   SCIP_CALL( SCIPallocMemoryArray(scip_, &stabcenterlinkingconss, nlinkingconssnew) );
   nstabcenterlinkingconss = nlinkingconssnew;
   BMSclearMemoryArray(stabcenterlinkingconss, nstabcenterlinkingconss);

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::setNConvconss(
      int nconvconssnew
      )
{
   SCIPfreeBlockMemoryArrayNull(scip_, &stabcenterconss, nstabcenterconv); /*lint !e64*/

   SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &stabcenterconv, nconvconssnew) );
   nstabcenterconv = nconvconssnew;
   BMSclearMemoryArray(stabcenterconv, nstabcenterconv);

   return SCIP_OKAY;
}

SCIP_Real Stabilization::linkingconsGetDual(
   int i
   )
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);

   assert(i < nstabcenterlinkingconss);
   assert(nstabcenterlinkingconss<= GCGgetNVarLinkingconss(origprob));
   assert(stabcenterlinkingconss != NULL);

   SCIP_CONS* cons = GCGgetVarLinkingconss(origprob)[i];

   return computeDual(stabcenterlinkingconss[i], pricingtype->consGetDual(scip_, cons));
}

SCIP_RETCODE Stabilization::consGetDual(
   int                   i,                  /* index of the constraint */
   SCIP_Real*            dual                /* return pointer for dual value */
)
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);
#ifndef NDEBUG
   int nconss =  GCGgetNMasterConss(origprob);
#endif
   assert(i < nconss);
   assert(dual != NULL);

   SCIP_CONS* cons = GCGgetMasterConss(origprob)[i];

   if( i >= nstabcenterconss )
      SCIP_CALL( updateStabcenterconss() );

   assert(i < nstabcenterconss);
   assert(stabcenterconss != NULL);

   *dual = computeDual(stabcenterconss[i], pricingtype->consGetDual(scip_, cons));
   return SCIP_OKAY;

}

SCIP_RETCODE Stabilization::rowGetDual(
   int                   i,                  /* index of the row */
   SCIP_Real*            dual                /* return pointer for dual value */
)
{
#ifndef NDEBUG
   int nrows = GCGsepaGetNCuts(scip_);
#endif
   assert(i < nrows);
   assert(dual != NULL);

   SCIP_ROW* row = GCGsepaGetMastercuts(scip_)[i];

   if( i >= nstabcentercuts )
      SCIP_CALL( updateStabcentercuts() );

   assert(i < nstabcentercuts);
   assert(stabcentercuts != NULL);

   *dual = computeDual(stabcentercuts[i], pricingtype->rowGetDual(row));

   return SCIP_OKAY;
}

SCIP_Real Stabilization::convGetDual(
   int i
   )
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);

   assert(i < nstabcenterconv);
   assert(nstabcenterconv<= GCGgetNPricingprobs(origprob));
   assert(stabcenterconv != NULL);

   SCIP_CONS* cons = GCGgetConvCons(origprob, i);

   return computeDual(stabcenterconv[i], pricingtype->consGetDual(scip_, cons));
}

SCIP_RETCODE Stabilization::updateStabilityCenter(
   SCIP_Real             lowerbound,         /**< lower bound due to lagrange function corresponding to current (stabilized) dual vars */
   SCIP_Real*            dualsolconv         /**< corresponding feasible dual solution for convexity constraints */
   )
{
   assert(dualsolconv != NULL);
   SCIPdebugMessage("Updating stability center: ");

   /* in case the bound is not improving and we have a stability center, do nothing */
   if( SCIPisLE(scip_, lowerbound, stabcenterbound) && hasstabilitycenter )
   {
      SCIPdebugPrintf("no bound increase: %g <= %g\n", lowerbound, SCIPnodeGetLowerbound(SCIPgetCurrentNode(scip_)));
      return SCIP_OKAY;
   }

   SCIPdebugPrintf("bound increase: %g > %g\n", lowerbound, SCIPnodeGetLowerbound(SCIPgetCurrentNode(scip_)));

   /* first update the arrays */
   SCIP_CALL( updateStabcenterconss() );
   SCIP_CALL( updateStabcentercuts() );

   /* get new dual values */
   SCIP* origprob = GCGmasterGetOrigprob(scip_);

   int nconss = GCGgetNMasterConss(origprob);
   int ncuts = GCGsepaGetNCuts(scip_);
   int nprobs = GCGgetNPricingprobs(origprob);

   assert(nstabcenterlinkingconss <= GCGgetNVarLinkingconss(origprob) );
   assert(nconss <= nstabcenterconss);
   assert(ncuts <= nstabcentercuts);

   for( int i = 0; i < nconss; ++i )
   {
      SCIP_CALL( consGetDual(i, &stabcenterconss[i]) );
   }

   for( int i = 0; i < ncuts; ++i )
   {
      SCIP_CALL( rowGetDual(i, &stabcentercuts[i]) );
   }

   for( int i = 0; i < nstabcenterlinkingconss; ++i)
   {
      stabcenterlinkingconss[i] = linkingconsGetDual(i);
   }

   for( int i = 0; i < nprobs; ++i )
   {
      if(!GCGisPricingprobRelevant(origprob, i))
         continue;

      stabcenterconv[i] = dualsolconv[i];
   }

   hasstabilitycenter = TRUE;
   stabcenterbound = lowerbound;

   return SCIP_OKAY;
}

SCIP_Real Stabilization::computeDual(
      SCIP_Real center,
      SCIP_Real current
      ) const
{
   SCIP_Real usedalpha = alpha;

   if (inmispricingschedule)
      usedalpha = alphabar;

   if( hasstabilitycenter )
      return usedalpha*center+(1.0-usedalpha)*current;
   else
      return current;
}

void Stabilization::updateIterationCount()
{
   ++t;
}

void Stabilization::updateIterationCountMispricing()
{
   ++k;
}



void Stabilization::updateNode()
{
   if( nodenr != SCIPnodeGetNumber(SCIPgetCurrentNode(scip_)) )
   {
      nodenr = SCIPnodeGetNumber(SCIPgetCurrentNode(scip_));
      k = 0;
      t = 1;
      alpha= 0.8;
      hasstabilitycenter = FALSE;
      stabcenterbound = -SCIPinfinity(scip_);
      inmispricingschedule = FALSE;
   }
}
void Stabilization::updateAlphaMisprice()
{
   SCIPdebugMessage("Alphabar update after mispricing\n");
   updateIterationCountMispricing();
   alphabar = MAX(0, 1-k*(1-alpha));
   SCIPdebugMessage("alphabar updated to %g in mispricing iteration k=%d and node pricing iteration t=%d \n", alphabar, k, t);
}

void Stabilization::updateAlpha(
   GCG_COL**            pricingcols         /**< solutions of the pricing problems */
   )
{
   SCIPdebugMessage("Alpha update after successful pricing\n");
   updateIterationCount();

   if( SCIPisPositive(scip_, calculateSubgradient(pricingcols)) )
   {
      increaseAlpha();
   }
   else
   {
      decreaseAlpha();
   }

}

void Stabilization::increaseAlpha()
{
   /* to avoid numerical problems, we assure alpha <= 0.9 */
   alpha = MIN(0.9, alpha+(1-alpha)*0.1);

   SCIPdebugMessage("alpha increased to %g\n", alpha);
}

void Stabilization::decreaseAlpha()
{
   if( alpha >= 0.5 && alpha < 1 )
   {
      alpha = alpha/1.1;
   }
   else
   {
      alpha = MAX(0, alpha-(1-alpha)*0.1);
   }
   SCIPdebugMessage("alpha decreased to %g\n", alpha);
}

SCIP_Real Stabilization::calculateSubgradient(
   GCG_COL**            pricingcols         /**< solutions of the pricing problems */
   )
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);
   SCIP_CONS** origmasterconss = GCGgetLinearOrigMasterConss(origprob);
   SCIP_CONS** masterconss = GCGgetMasterConss(origprob);

   SCIP_CONS** linkingconss = GCGgetVarLinkingconss(origprob);
   int nlinkingconss = GCGgetNVarLinkingconss(origprob);
   int* linkingconsblocks = GCGgetVarLinkingconssBlock(origprob);
   assert(nstabcenterlinkingconss <= GCGgetNVarLinkingconss(origprob) );
   int nconss = GCGgetNMasterConss(origprob);
   assert(nconss <= nstabcenterconss);
   SCIP_ROW** mastercuts = GCGsepaGetMastercuts(scip_);
   SCIP_ROW** origmastercuts = GCGsepaGetOrigcuts(scip_);
   int ncuts = GCGsepaGetNCuts(scip_);
   assert(ncuts <= nstabcentercuts);

   SCIP_Real gradientproduct = 0.0;

   /* masterconss */
   for( int i = 0; i < nconss; ++i )
   {
      SCIP_VAR** vars;
      SCIP_Real* vals;
      int nvars;
      SCIP_Real lhs; /* can also be rhs, but we need only one */

      SCIP_CONS* origcons = origmasterconss[i];
      nvars = SCIPgetNVarsLinear(origprob, origcons);
      vars = SCIPgetVarsLinear(origprob, origcons);
      vals = SCIPgetValsLinear(origprob, origcons);

      SCIP_Real dual =  pricingtype->consGetDual(scip_, masterconss[i]);
      assert(!SCIPisInfinity(scip_, ABS(dual)));

      for( int j = 0; j < nvars; ++j )
      {
         SCIP_Real val = 0.0;
         assert(GCGvarIsOriginal(vars[j]));
         if( GCGvarGetBlock(vars[j]) < 0 )
         {
            SCIP_VAR* mastervar = GCGoriginalVarGetMastervars(vars[j])[0];
            assert(GCGvarIsMaster(mastervar));
            val = SCIPgetSolVal(scip_, (SCIP_SOL*) NULL, mastervar);
            assert( !SCIPisInfinity(scip_, val) );
         }
         else
         {
            int block = GCGvarGetBlock(vars[j]);
            if( !GCGisPricingprobRelevant(origprob, block) )
               continue;

            SCIP_VAR* pricingvar = GCGoriginalVarGetPricingVar(vars[j]);
            assert(GCGvarIsPricing(pricingvar));
            SCIP* pricingprob = GCGgetPricingprob(origprob, block);
            assert(pricingprob != NULL);
            val = GCGcolGetSolVal(pricingprob, pricingcols[block], pricingvar);
            assert(!SCIPisInfinity(scip_, ABS(val)));
         }
         assert(stabcenterconss != NULL);
         assert(vals != NULL);
         gradientproduct += (stabcenterconss[i] - dual) * vals[j] * val;
      }

      if( SCIPisFeasPositive(scip_, dual) )
      {
         lhs = SCIPgetLhsLinear(origprob, origcons);
      }
      else if( SCIPisFeasNegative(scip_, dual) )
      {
         lhs = SCIPgetRhsLinear(origprob, origcons);
      }
      else
      {
         continue;
      }
      assert(stabcenterconss != NULL);
      assert(!SCIPisInfinity(scip_, ABS(lhs)));

      gradientproduct -= (stabcenterconss[i] - dual) * lhs;
   }

   /* mastercuts */
   for( int i = 0; i < ncuts; ++i )
   {
      SCIP_COL** cols;
      SCIP_Real* vals;
      int nvars;
      SCIP_Real lhs; /* can also be rhs, but we need only one */

      SCIP_ROW* origcut = origmastercuts[i];
      nvars = SCIProwGetNNonz(origcut);
      cols = SCIProwGetCols(origcut);
      vals = SCIProwGetVals(origcut);

      SCIP_Real dual = pricingtype->rowGetDual(mastercuts[i]);
      assert(!SCIPisInfinity(scip_, ABS(dual)));
      for( int j = 0; j < nvars; ++j )
      {
         SCIP_Real val = 0.0;
         SCIP_VAR* var = SCIPcolGetVar(cols[j]);
         assert(GCGvarIsOriginal(var));

         /* Linking or master variable */
         if( GCGvarGetBlock(var) < 0 )
         {
            SCIP_VAR* mastervar = GCGoriginalVarGetMastervars(var)[0];
            assert(GCGvarIsMaster(mastervar));
            val = SCIPgetSolVal(scip_, (SCIP_SOL*) NULL, mastervar);
            assert(!SCIPisInfinity(scip_, ABS(val)));
         }
         /* Variable in a pricing problem */
         else
         {
            int block = GCGvarGetBlock(var);
            if( !GCGisPricingprobRelevant(origprob, block) )
               continue;

            SCIP_VAR* pricingvar = GCGoriginalVarGetPricingVar(var);
            assert(GCGvarIsPricing(pricingvar));
            SCIP* pricingprob = GCGgetPricingprob(origprob, block);
            assert(pricingprob != NULL);
            val = GCGcolGetSolVal(pricingprob, pricingcols[block], pricingvar);
            assert(!SCIPisInfinity(scip_, ABS(val)));
         }
         assert(stabcentercuts != NULL);
         assert(vals != NULL);
         gradientproduct += (stabcentercuts[i] - dual) * vals[j] * val;
      }

      if( SCIPisFeasGT(scip_, dual, 0.0) )
      {
         lhs = SCIProwGetLhs(origcut);
      }
      else if( SCIPisFeasLT(scip_, dual, 0.0) )
      {
         lhs = SCIProwGetRhs(origcut);
      }
      else
      {
         continue;
      }
      assert(!SCIPisInfinity(scip_, ABS(lhs)));
      assert(stabcentercuts != NULL);

      gradientproduct -=  (stabcentercuts[i] - dual) * lhs;
   }

   /* linkingconss */
   for( int i = 0; i < nlinkingconss; ++i )
   {
      SCIP_VAR* mastervar;
      SCIP_VAR* pricingvar;
      SCIP_CONS* linkingcons = linkingconss[i];
      int block = linkingconsblocks[i];
      mastervar = SCIPgetVarsLinear(scip_, linkingcons)[0];
      assert(GCGvarIsMaster(mastervar));

      pricingvar = GCGlinkingVarGetPricingVars(GCGmasterVarGetOrigvars(mastervar)[0])[block];
      assert(GCGvarIsPricing(pricingvar));
      SCIP* pricingprob = GCGgetPricingprob(origprob, block);
      assert(pricingprob != NULL);

      assert(stabcenterlinkingconss != NULL);
      SCIP_Real dual = stabcenterlinkingconss[i] - pricingtype->consGetDual(scip_, linkingcons);
      SCIP_Real masterval = SCIPgetSolVal(scip_, (SCIP_SOL*) NULL, mastervar);
      SCIP_Real pricingval = GCGcolGetSolVal(pricingprob, pricingcols[block], pricingvar);
      assert(!SCIPisInfinity(scip_, ABS(masterval)));
      assert(!SCIPisInfinity(scip_, ABS(pricingval)));
      assert(!SCIPisInfinity(scip_, ABS(dual)));
      gradientproduct += dual * (masterval - pricingval);
   }

   SCIPdebugMessage("Update gradient with value %g.\n", gradientproduct);

   return gradientproduct;
}

SCIP_Bool Stabilization::isStabilized()
{
   if(inmispricingschedule)
      return SCIPisGT(scip_, alphabar, 0.0);
   return SCIPisGT(scip_, alpha, 0.0);
}

/** enabling mispricing schedule */
void Stabilization::activateMispricingSchedule(
)
{
   inmispricingschedule = TRUE;
}

/** disabling mispricing schedule */
void Stabilization::disablingMispricingSchedule(
)
{
   inmispricingschedule = FALSE;
   k=0;
}

/** is mispricing schedule enabled */
SCIP_Bool Stabilization::isInMispricingSchedule(
) const
{
   return inmispricingschedule;
}

} /* namespace gcg */
