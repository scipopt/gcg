/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       */
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
 * This is an implementation of dynamic alpha-schedule (based on subgradient information) stabilization including an optional
 * combination with a subgradient method based on the papers
 *
 * Pessoa, A., Sadykov, R., Uchoa, E., & Vanderbeck, F. (2013). In-Out Separation and Column Generation
 * Stabilization by Dual Price Smoothing. In Experimental Algorithms (pp. 354-365). Springer Berlin Heidelberg.
 *
 * Pessoa, A., Sadykov, R., Uchoa, E., & Vanderbeck, F. (2016). Automation and combination of linear-programming
 * based stabilization techniques in column generation.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */
#include "class_stabilization.h"
#include "pricer_gcg.h"
#include "gcg.h"
#include "sepa_master.h"
#include "objscip/objscip.h"
#include "scip/cons_linear.h"

namespace gcg {

Stabilization::Stabilization(
   SCIP* scip,
   PricingType* pricingtype_,
   SCIP_Bool hybridascent_
   ) :scip_(scip), stabcenterconss((SCIP_Real*) NULL), stabcenterconsssize(0), nstabcenterconss(0),
      stabcentercuts((SCIP_Real*) NULL), stabcentercutssize(0), nstabcentercuts(0),
      stabcenterlinkingconss((SCIP_Real*) NULL), nstabcenterlinkingconss(0),
      stabcenterconv((SCIP_Real*) NULL), nstabcenterconv(0), dualdiffnorm(0.0),
      subgradientconss(NULL), subgradientconsssize(0), nsubgradientconss(0),
      subgradientcuts(NULL), subgradientcutssize(0), nsubgradientcuts(0),
      subgradientlinkingconss(NULL), nsubgradientlinkingconss(0),
      subgradientnorm(0.0), hybridfactor(0.0),
      pricingtype(pricingtype_), alpha(0.8), alphabar(0.8), hybridascent(hybridascent_), beta(0.0), nodenr(-1), k(0), t(0), hasstabilitycenter(FALSE),stabcenterbound(-SCIPinfinity(scip)),
      inmispricingschedule(FALSE), subgradientproduct(0.0), farkasalpha(1.0), farkasalphabar(1.0), infarkas(FALSE)
{

}

Stabilization::~Stabilization()
{
   SCIPfreeBlockMemoryArrayNull(scip_, &stabcenterconss, stabcenterconsssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(scip_, &stabcentercuts, stabcentercutssize); /*lint !e64*/
   SCIPfreeMemoryArrayNull(scip_, &stabcenterlinkingconss); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(scip_, &subgradientconss, subgradientconsssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(scip_, &subgradientcuts, subgradientcutssize); /*lint !e64*/
   SCIPfreeMemoryArrayNull(scip_, &subgradientlinkingconss); /*lint !e64*/
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

SCIP_RETCODE Stabilization::updateSubgradientconss()
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);
   int nconss = GCGgetNMasterConss(origprob);

   if( nconss == nsubgradientconss )
   {
      return SCIP_OKAY;
   }

   if( nconss > subgradientconsssize )
   {
      int oldsize = subgradientconsssize;
      subgradientconsssize = SCIPcalcMemGrowSize(scip_, nconss);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip_, &subgradientconss, oldsize, subgradientconsssize) );
   }
   assert(subgradientconss != NULL);
   BMSclearMemoryArray(&subgradientconss[nsubgradientconss], (size_t)nconss-nsubgradientconss); /*lint !e866*/

   nsubgradientconss = nconss;

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::updateSubgradientcuts()
{
   int ncuts = GCGsepaGetNCuts(scip_);

   if( ncuts == nsubgradientcuts )
   {
      return SCIP_OKAY;
   }

   if( ncuts > subgradientcutssize )
   {
      int oldsize = subgradientcutssize;
      subgradientcutssize = SCIPcalcMemGrowSize(scip_, ncuts);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip_, &subgradientcuts, oldsize, subgradientcutssize) );
   }
   assert(subgradientcuts != NULL);
   BMSclearMemoryArray(&subgradientcuts[nsubgradientcuts], (size_t)ncuts-nsubgradientcuts); /*lint !e866*/

   nsubgradientcuts = ncuts;

   return SCIP_OKAY;
}


SCIP_RETCODE Stabilization::setNLinkingconss(
      int nlinkingconssnew
      )
{

   SCIPfreeMemoryArrayNull(scip_, &stabcenterlinkingconss); /*lint !e64*/
   SCIP_CALL( SCIPallocMemoryArray(scip_, &stabcenterlinkingconss, nlinkingconssnew) );

   if( hybridascent )
   {
      SCIPfreeMemoryArrayNull(scip_, &subgradientlinkingconss); /*lint !e64*/
      SCIP_CALL( SCIPallocMemoryArray(scip_, &subgradientlinkingconss, nlinkingconssnew) );
      BMSclearMemoryArray(subgradientlinkingconss, nlinkingconssnew);
   }

   nstabcenterlinkingconss = nlinkingconssnew;
   BMSclearMemoryArray(stabcenterlinkingconss, nstabcenterlinkingconss);


   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::setNConvconss(
      int nconvconssnew
      )
{
   SCIPfreeBlockMemoryArrayNull(scip_, &stabcenterconv, nstabcenterconv); /*lint !e64*/
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
   SCIP_Real subgradient = 0.0;

   assert(i < nstabcenterlinkingconss);
   assert(nstabcenterlinkingconss<= GCGgetNVarLinkingconss(origprob));
   assert(stabcenterlinkingconss != NULL);

   SCIP_CONS* cons = GCGgetVarLinkingconss(origprob)[i];

   if( hybridascent && hasstabilitycenter )
      subgradient = subgradientlinkingconss[i];

   return computeDual(stabcenterlinkingconss[i], pricingtype->consGetDual(scip_, cons), subgradient, 0.0, 0.0);
}

SCIP_RETCODE Stabilization::consGetDual(
   int                   i,                  /* index of the constraint */
   SCIP_Real*            dual                /* return pointer for dual value */
)
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);
   SCIP_Real subgradient = 0.0;
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

   if( i >= nsubgradientconss && hybridascent )
      SCIP_CALL( updateSubgradientconss() );

   if( hybridascent && hasstabilitycenter )
      subgradient = subgradientconss[i];

   *dual = computeDual(stabcenterconss[i], pricingtype->consGetDual(scip_, cons), subgradient, SCIPgetLhsLinear(scip_, cons), SCIPgetRhsLinear(scip_, cons));
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
   SCIP_Real subgradient = 0.0;

   if( i >= nstabcentercuts )
      SCIP_CALL( updateStabcentercuts() );

   assert(i < nstabcentercuts);
   assert(stabcentercuts != NULL);

   if( i >= nsubgradientcuts && hybridascent )
      SCIP_CALL( updateSubgradientcuts() );

   if( hybridascent && hasstabilitycenter )
   {
      assert(subgradientcuts != NULL);
      subgradient = subgradientcuts[i];
   }

   *dual = computeDual(stabcentercuts[i], pricingtype->rowGetDual(row), subgradient, SCIProwGetLhs(row), SCIProwGetRhs(row));

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
   SCIP_Real subgradient = 0.0;

   return computeDual(stabcenterconv[i], pricingtype->consGetDual(scip_, cons), subgradient, (SCIP_Real) GCGgetNIdenticalBlocks(origprob, i), (SCIP_Real) GCGgetNIdenticalBlocks(origprob, i));
}

SCIP_RETCODE Stabilization::updateStabilityCenter(
   SCIP_Real             lowerbound,         /**< lower bound due to lagrange function corresponding to current (stabilized) dual vars */
   SCIP_Real*            dualsolconv,        /**< corresponding feasible dual solution for convexity constraints */
   GCG_COL**             pricingcols         /**< columns of the pricing problems */
   )
{
   assert(dualsolconv != NULL);
   SCIPdebugMessage("Updating stability center: ");

   /* in case the bound is not improving and we have a stability center, do nothing */
   if( infarkas )
   {
      return SCIP_OKAY;
   }
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

   if( hybridascent )
   {
      SCIP_CALL( updateSubgradientconss() );
      SCIP_CALL( updateSubgradientcuts() );
   }

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

   if( hybridascent )
      calculateSubgradient(pricingcols);

   hasstabilitycenter = TRUE;
   stabcenterbound = lowerbound;

   return SCIP_OKAY;
}

SCIP_Real Stabilization::computeDual(
      SCIP_Real         center,
      SCIP_Real         current,
      SCIP_Real         subgradient,         /**< subgradient (or 0.0 if not needed) */
      SCIP_Real         lhs,                 /**< lhs (or 0.0 if not needed) */
      SCIP_Real         rhs                  /**< rhs (or 0.0 if not needed) */
      ) const
{
   SCIP_Real usedalpha = alpha;
   SCIP_Real usedbeta = beta;

   if ( inmispricingschedule )
   {
      usedalpha = alphabar;
      usedbeta = 0.0;
   }

   if( hasstabilitycenter && (SCIPisZero(scip_, usedbeta) || SCIPisZero(scip_, usedalpha)) )
      return usedalpha*center+(1.0-usedalpha)*current;
   else if( hasstabilitycenter && SCIPisPositive(scip_, usedbeta) )
   {
      SCIP_Real dual = center + hybridfactor * (beta * (center + subgradient * dualdiffnorm / subgradientnorm) + (1.0 - beta) * current - center);

      /* make sure dual solution has the correct sign */
      if( SCIPisInfinity(scip_, rhs) )
         dual = MAX(dual, 0.0);
      else if( SCIPisInfinity(scip_, -lhs) )
         dual = MIN(dual, 0.0);

      return dual;
   }
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
      alpha = 0.8;
      hasstabilitycenter = FALSE;
      stabcenterbound = -SCIPinfinity(scip_);
      inmispricingschedule = FALSE;

      alphabar = 0.8;
   }
}

/**< update information for hybrid stabilization with dual ascent */
SCIP_RETCODE Stabilization::updateHybrid()
{
   if( hasstabilitycenter && hybridascent && !inmispricingschedule && !infarkas )
   {
      /* first update the arrays */
      SCIP_CALL( updateStabcenterconss() );
      SCIP_CALL( updateStabcentercuts() );

      SCIP_CALL( updateSubgradientconss() );
      SCIP_CALL( updateSubgradientcuts() );

      if( SCIPisPositive(scip_, alpha) )
      {
         calculateDualdiffnorm();
         calculateBeta();
         calculateHybridFactor();
      }
   }

   return SCIP_OKAY;
}

void Stabilization::updateAlphaMisprice()
{
   SCIPdebugMessage("Alphabar update after mispricing\n");
   updateIterationCountMispricing();
   if( infarkas )
   {
      farkasalphabar = MAX(0.0, 1-k*(1-farkasalpha));
      SCIPdebugMessage("farkasalphabar updated to %g in mispricing iteration k=%d and node pricing iteration t=%d \n", farkasalphabar, k, t);
   }
   else
   {
      alphabar = MAX(0.0, 1-k*(1-alpha));
      SCIPdebugMessage("alphabar updated to %g in mispricing iteration k=%d and node pricing iteration t=%d \n", alphabar, k, t);
   }
}

void Stabilization::updateAlpha(
   GCG_COL**            pricingcols         /**< solutions of the pricing problems */
   )
{
   if( !infarkas )
   {

   }
   else
   {
      SCIPdebugMessage("Alpha update after successful pricing\n");
      updateIterationCount();

      /* There is a sign error in the stabilization paper:
       * if the scalar product (subgradientproduct) is positive, the angle is less than 90° and we want to decrease alpha
       */
      if( SCIPisNegative(scip_, subgradientproduct) )
      {
         increaseAlpha();
      }
      else
      {
         decreaseAlpha();
      }
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
   alpha = MAX(0.0, alpha-0.1);

   SCIPdebugMessage("alpha decreased to %g\n", alpha);
}

SCIP_Real Stabilization::calculateSubgradientProduct(
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
      SCIP_Real stabdual;

      SCIP_CALL( consGetDual(i, &stabdual) );

      assert(!SCIPisInfinity(scip_, ABS(dual)));

      if( SCIPisFeasPositive(scip_, stabdual) )
      {
         lhs = SCIPgetLhsLinear(origprob, origcons);
      }
      else if( SCIPisFeasNegative(scip_, stabdual) )
      {
         lhs = SCIPgetRhsLinear(origprob, origcons);
      }
      else
      {
         continue;
      }

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
         gradientproduct -= (dual - stabcenterconss[i]) * vals[j] * val;
      }

      assert(stabcenterconss != NULL);
      assert(!SCIPisInfinity(scip_, ABS(lhs)));

      gradientproduct += (dual - stabcenterconss[i]) * lhs;
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

      SCIP_Real stabdual;

      SCIP_CALL( rowGetDual(i, &stabdual) );

      if( SCIPisFeasGT(scip_, stabdual, 0.0) )
      {
         lhs = SCIProwGetLhs(origcut);
      }
      else if( SCIPisFeasLT(scip_, stabdual, 0.0) )
      {
         lhs = SCIProwGetRhs(origcut);
      }
      else
      {
         continue;
      }
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
         gradientproduct -= (dual - stabcentercuts[i]) * vals[j] * val;
      }

      assert(!SCIPisInfinity(scip_, ABS(lhs)));
      assert(stabcentercuts != NULL);

      gradientproduct +=  (dual - stabcentercuts[i]) * lhs;
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
      SCIP_Real dual = pricingtype->consGetDual(scip_, linkingcons) - stabcenterlinkingconss[i];

      SCIP_Real stabdual = linkingconsGetDual(i);

      if( SCIPisFeasZero(origprob, stabdual) )
         continue;

      SCIP_Real masterval = SCIPgetSolVal(scip_, (SCIP_SOL*) NULL, mastervar);
      SCIP_Real pricingval = GCGcolGetSolVal(pricingprob, pricingcols[block], pricingvar);
      assert(!SCIPisInfinity(scip_, ABS(masterval)));
      assert(!SCIPisInfinity(scip_, ABS(pricingval)));
      assert(!SCIPisInfinity(scip_, ABS(dual)));
      gradientproduct -= dual * (masterval - pricingval);
   }

   SCIPdebugMessage("Update gradient product with value %g.\n", gradientproduct);

   return gradientproduct;
}

/** calculates the subgradient (with linking variables) */
void Stabilization::calculateSubgradient(
   GCG_COL**            pricingcols         /**< columns of the pricing problems */
)
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);
   SCIP_CONS** origmasterconss = GCGgetLinearOrigMasterConss(origprob);

   SCIP_CONS** linkingconss = GCGgetVarLinkingconss(origprob);
   int nlinkingconss = GCGgetNVarLinkingconss(origprob);
   int* linkingconsblocks = GCGgetVarLinkingconssBlock(origprob);
   assert(nstabcenterlinkingconss <= GCGgetNVarLinkingconss(origprob) );
   int nconss = GCGgetNMasterConss(origprob);
   assert(nconss <= nstabcenterconss);
   SCIP_ROW** origmastercuts = GCGsepaGetOrigcuts(scip_);
   int ncuts = GCGsepaGetNCuts(scip_);
   assert(ncuts <= nstabcentercuts);

   subgradientnorm = 0.0;

   /* masterconss */
   for( int i = 0; i < nconss; ++i )
   {
      SCIP_VAR** vars;
      SCIP_Real* vals;
      int nvars;
      SCIP_Real activity;
      SCIP_Real infeasibility;

      SCIP_CONS* origcons = origmasterconss[i];
      nvars = SCIPgetNVarsLinear(origprob, origcons);
      vars = SCIPgetVarsLinear(origprob, origcons);
      vals = SCIPgetValsLinear(origprob, origcons);

      SCIP_Real dual = stabcenterconss[i];
      assert(!SCIPisInfinity(scip_, ABS(dual)));

      activity = 0.0;

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
         assert(vals != NULL);
         activity += vals[j] * val;
      }

      infeasibility = 0.0;

      if( SCIPisFeasPositive(scip_, dual) /* || SCIPisInfinity(origprob, SCIPgetRhsLinear(origprob, origcons)) */)
      {
         infeasibility = SCIPgetLhsLinear(origprob, origcons) - activity;
      }
      else if( SCIPisFeasNegative(scip_, dual) /* || SCIPisInfinity(origprob, SCIPgetLhsLinear(origprob, origcons)) */)
      {
         infeasibility = SCIPgetRhsLinear(origprob, origcons) - activity;
      }

      assert(subgradientconss != NULL);
      assert(!SCIPisInfinity(scip_, SQR(infeasibility)));

      subgradientconss[i] = infeasibility;

      if( SCIPisPositive(scip_, SQR(infeasibility)) )
         subgradientnorm += SQR(infeasibility);
   }

   /* mastercuts */
   for( int i = 0; i < ncuts; ++i )
   {
      SCIP_COL** cols;
      SCIP_Real* vals;
      int nvars;
      SCIP_Real activity;
      SCIP_Real infeasibility;

      SCIP_ROW* origcut = origmastercuts[i];
      nvars = SCIProwGetNNonz(origcut);
      cols = SCIProwGetCols(origcut);
      vals = SCIProwGetVals(origcut);

      activity = 0.0;

      SCIP_Real dual = stabcentercuts[i];
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
         activity += vals[j] * val;
      }

      infeasibility = 0.0;

      if( SCIPisFeasPositive(scip_, dual) )
      {
         infeasibility = SCIProwGetLhs(origcut) - activity;
      }
      else if( SCIPisFeasNegative(scip_, dual) )
      {
         infeasibility = SCIProwGetRhs(origcut) - activity;
      }

      assert(subgradientcuts != NULL);
      assert(!SCIPisInfinity(scip_, SQR(infeasibility)));

      subgradientcuts[i] = infeasibility;

      if( SCIPisPositive(scip_, SQR(infeasibility)) )
         subgradientnorm += SQR(infeasibility);
   }

   /* linkingconss */
   for( int i = 0; i < nlinkingconss; ++i )
   {
      SCIP_VAR* mastervar;
      SCIP_VAR* pricingvar;
      SCIP_CONS* linkingcons = linkingconss[i];
      int block = linkingconsblocks[i];
      SCIP_Real activity;
      SCIP_Real infeasibility;
      mastervar = SCIPgetVarsLinear(scip_, linkingcons)[0];
      assert(GCGvarIsMaster(mastervar));

      pricingvar = GCGlinkingVarGetPricingVars(GCGmasterVarGetOrigvars(mastervar)[0])[block];
      assert(GCGvarIsPricing(pricingvar));
      SCIP* pricingprob = GCGgetPricingprob(origprob, block);
      assert(pricingprob != NULL);

      assert(stabcenterlinkingconss != NULL);
      SCIP_Real masterval = SCIPgetSolVal(scip_, (SCIP_SOL*) NULL, mastervar);
      SCIP_Real pricingval = GCGcolGetSolVal(pricingprob, pricingcols[block], pricingvar);
      assert(!SCIPisInfinity(scip_, ABS(masterval)));
      assert(!SCIPisInfinity(scip_, ABS(pricingval)));
      activity = (masterval - pricingval);

      infeasibility = activity;

      assert(subgradientlinkingconss != NULL);
      assert(!SCIPisInfinity(scip_, SQR(infeasibility)));

      subgradientlinkingconss[i] = infeasibility;

      if( SCIPisPositive(scip_, SQR(infeasibility)) )
         subgradientnorm += SQR(infeasibility);
   }


   assert(!SCIPisNegative(scip_, subgradientnorm));

   subgradientnorm = SQRT(subgradientnorm);

   SCIPdebugMessage("Update subgradient and subgradientnorm with value %g.\n", subgradientnorm);
}

/**< calculate norm of difference between stabcenter and current duals */
void Stabilization::calculateDualdiffnorm()
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);
   SCIP_CONS** masterconss = GCGgetMasterConss(origprob);
   SCIP_CONS** linkingconss = GCGgetVarLinkingconss(origprob);
   int nlinkingconss = GCGgetNVarLinkingconss(origprob);
   assert(nstabcenterlinkingconss <= GCGgetNVarLinkingconss(origprob) );
   int nconss = GCGgetNMasterConss(origprob);
   assert(nconss <= nstabcenterconss);
   SCIP_ROW** mastercuts = GCGsepaGetMastercuts(scip_);
   int ncuts = GCGsepaGetNCuts(scip_);
   assert(ncuts <= nstabcentercuts);

   dualdiffnorm = 0.0;

   /* masterconss */
   assert(stabcenterconss != NULL);

   for( int i = 0; i < nconss; ++i )
   {
      SCIP_Real dualdiff = SQR(stabcenterconss[i] - pricingtype->consGetDual(scip_, masterconss[i]));

      if( SCIPisPositive(scip_, dualdiff) )
         dualdiffnorm += dualdiff;
   }

   /* mastercuts */
   assert(stabcenterconss != NULL);

   for( int i = 0; i < ncuts; ++i )
   {
      SCIP_Real dualdiff = SQR(stabcentercuts[i] - pricingtype->rowGetDual(mastercuts[i]));

      if( SCIPisPositive(scip_, dualdiff) )
         dualdiffnorm += dualdiff;
   }

   /* linkingconss */
   assert(stabcenterlinkingconss != NULL);

   for( int i = 0; i < nlinkingconss; ++i )
   {
      SCIP_Real dualdiff = SQR(stabcenterlinkingconss[i] - pricingtype->consGetDual(scip_, linkingconss[i]));

      if( SCIPisPositive(scip_, dualdiff) )
         dualdiffnorm += dualdiff;
   }

   dualdiffnorm = SQRT(dualdiffnorm);
   SCIPdebugMessage("Update dualdiffnorm with value %g.\n", dualdiffnorm);
}

/**< calculate beta */
void Stabilization::calculateBeta()
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);
   SCIP_CONS** masterconss = GCGgetMasterConss(origprob);
   SCIP_CONS** linkingconss = GCGgetVarLinkingconss(origprob);
   int nlinkingconss = GCGgetNVarLinkingconss(origprob);
   assert(nstabcenterlinkingconss <= GCGgetNVarLinkingconss(origprob) );
   int nconss = GCGgetNMasterConss(origprob);
   assert(nconss <= nstabcenterconss);
   SCIP_ROW** mastercuts = GCGsepaGetMastercuts(scip_);
   int ncuts = GCGsepaGetNCuts(scip_);
   assert(ncuts <= nstabcentercuts);

   beta = 0.0;

   /* masterconss */
   assert(stabcenterconss != NULL);

   for( int i = 0; i < nconss; ++i )
   {
      SCIP_Real dualdiff = ABS(pricingtype->consGetDual(scip_, masterconss[i]) - stabcenterconss[i]);
      SCIP_Real product = dualdiff * ABS(subgradientconss[i]);

      if( SCIPisPositive(scip_, product) )
         beta += product;
   }

   /* mastercuts */
   assert(stabcentercuts != NULL || ncuts == 0);

   for( int i = 0; i < ncuts; ++i )
   {
      SCIP_Real dualdiff = ABS(pricingtype->rowGetDual(mastercuts[i]) - stabcentercuts[i]);
      SCIP_Real product = dualdiff * ABS(subgradientcuts[i]);

      if( SCIPisPositive(scip_, product) )
         beta += product;
   }

   /* linkingconss */
   assert(stabcenterlinkingconss != NULL);

   for( int i = 0; i < nlinkingconss; ++i )
   {
      SCIP_Real dualdiff = ABS(pricingtype->consGetDual(scip_, linkingconss[i]) - stabcenterlinkingconss[i]);
      SCIP_Real product = dualdiff * ABS(subgradientlinkingconss[i]);

      if( SCIPisPositive(scip_, product) )
         beta += product;
   }

   if( SCIPisPositive(scip_, subgradientnorm) )
      beta = beta / (subgradientnorm * dualdiffnorm);

   SCIPdebugMessage("Update beta with value %g.\n", beta);

   assert( ( SCIPisPositive(scip_, beta) || SCIPisZero(scip_, subgradientnorm)) && SCIPisLE(scip_, beta, 1.0) );
}

/**< calculate factor that is needed in hybrid stabilization */
void Stabilization::calculateHybridFactor()
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);
   SCIP_CONS** masterconss = GCGgetMasterConss(origprob);

   SCIP_CONS** linkingconss = GCGgetVarLinkingconss(origprob);
   int nlinkingconss = GCGgetNVarLinkingconss(origprob);
   assert(nstabcenterlinkingconss <= GCGgetNVarLinkingconss(origprob) );
   int nconss = GCGgetNMasterConss(origprob);
   assert(nconss <= nstabcenterconss);
   SCIP_ROW** mastercuts = GCGsepaGetMastercuts(scip_);
   int ncuts = GCGsepaGetNCuts(scip_);
   assert(ncuts <= nstabcentercuts);

   SCIP_Real divisornorm = 0.0;

   /* masterconss */
   assert(stabcenterconss != NULL);

   for( int i = 0; i < nconss; ++i )
   {
      SCIP_Real divisor = SQR((beta - 1.0) * stabcenterconss[i]
                        + beta * (subgradientconss[i] * dualdiffnorm / subgradientnorm)
                        + (1 - beta) * pricingtype->consGetDual(scip_, masterconss[i]));

      if( SCIPisPositive(scip_, divisor) )
         divisornorm += divisor;
   }

   /* mastercuts */
   assert(stabcenterconss != NULL);

   for( int i = 0; i < ncuts; ++i )
   {
      SCIP_Real divisor = SQR((beta - 1.0) * stabcentercuts[i]
                        + beta * (subgradientcuts[i] * dualdiffnorm / subgradientnorm)
                        + (1 - beta) * pricingtype->rowGetDual(mastercuts[i]));

      if( SCIPisPositive(scip_, divisor) )
         divisornorm += divisor;
   }

   /* linkingconss */
   assert(stabcenterlinkingconss != NULL);

   for( int i = 0; i < nlinkingconss; ++i )
   {
      SCIP_Real divisor = SQR((beta - 1.0) * stabcenterlinkingconss[i]
                        + beta * (subgradientlinkingconss[i] * dualdiffnorm / subgradientnorm)
                        + (1 - beta) * pricingtype->consGetDual(scip_, linkingconss[i]));

      if( SCIPisPositive(scip_, divisor) )
         divisornorm += divisor;
   }

   divisornorm = SQRT(divisornorm);

   hybridfactor = ((1 - alpha) * dualdiffnorm) / divisornorm;

   SCIPdebugMessage("Update hybridfactor with value %g.\n", hybridfactor);

   assert( SCIPisPositive(scip_, hybridfactor) );
}


SCIP_Bool Stabilization::isStabilized()
{
   if( infarkas )
   {
      if( inmispricingschedule)
         return SCIPisGT(scip_, farkasalphabar, 0.0);
      return SCIPisGT(scip_, farkasalpha, 0.0);
   }
   if( inmispricingschedule)
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

/** enabling Farkas */
void Stabilization::activateFarkas(
)
{
   infarkas = TRUE;
}

/** disabling Farkas */
void Stabilization::disablingFarkas(
)
{
   infarkas = FALSE;
}

/** in Farkas*/
SCIP_Bool Stabilization::inFarkas(
) const
{
   return infarkas;
}
/** get Farkas alpha */
SCIP_Real Stabilization::getFarkasAlpha(
) const
{
   assert(infarkas);

   if( inmispricingschedule )
      return farkasalphabar;

   return farkasalpha;
}

/** update subgradient product */
void Stabilization::updateSubgradientProduct(
   GCG_COL**            pricingcols         /**< solutions of the pricing problems */
)
{
   if( !infarkas )
      subgradientproduct = calculateSubgradientProduct(pricingcols);
}


} /* namespace gcg */
