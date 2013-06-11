/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
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
 * @brief  Description
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_stabilization.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "sepa_master.h"
#include "objscip/objscip.h"
#include "scip/cons_linear.h"
#include "pub_gcgvar.h"

namespace gcg {

Stabilization::Stabilization(
   SCIP* scip,
   PricingType* pricingtype_
   ) :scip_(scip), stabcenterconss(NULL), stabcenterconsssize(0), nstabcenterconss(0),
      stabcentercuts(NULL), stabcentercutssize(0), nstabcentercuts(0),
      stabcenterlinkingconss(NULL), nstabcenterlinkingconss(0),
      pricingtype(pricingtype_), alpha(0.8), node(NULL), k(0), hasstabilitycenter(FALSE)
{

}

Stabilization::~Stabilization()
{
   // TODO Auto-generated destructor stub
}

SCIP_RETCODE Stabilization::updateStabcenterconss()
{
   SCIP* origprob = GCGpricerGetOrigprob(scip_);
   int nconss = GCGrelaxGetNMasterConss(origprob);

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
   BMSclearMemoryArray(&stabcenterconss[nstabcenterconss], nconss-nstabcenterconss);

   nstabcenterconss = nconss;

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::updateStabcentercuts()
{
   int ncuts = GCGsepaGetNMastercuts(scip_);

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

   BMSclearMemoryArray(&stabcentercuts[nstabcentercuts], ncuts-nstabcentercuts);

   nstabcentercuts = ncuts;

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::setNLinkingconss(
      int nlinkingconssnew
      )
{
   SCIPfreeBlockMemoryArrayNull(scip_, &stabcenterlinkingconss, nstabcenterlinkingconss);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &stabcenterlinkingconss, SCIPcalcMemGrowSize(scip_, nlinkingconssnew)) );
   nstabcenterlinkingconss = nlinkingconssnew;
   BMSclearMemoryArray(stabcenterlinkingconss, nstabcenterlinkingconss);

   return SCIP_OKAY;
}


SCIP_RETCODE Stabilization::linkingconsGetDual(
   int i,
   SCIP_Real* dual
   )
{
   SCIP* origprob = GCGpricerGetOrigprob(scip_);

   assert(nstabcenterlinkingconss<= GCGrelaxGetNLinkingconss(origprob));

   SCIP_CONS* cons = GCGrelaxGetLinkingconss(origprob)[i];

   *dual = computeDual(stabcenterlinkingconss[i], pricingtype->consGetDual(scip_, cons));
   //SCIPdebugMessage("dual for constraint <%s> is %g (alpha = %g, stab = %g, dual = %g).\n",
   //      SCIPconsGetName(cons), *dual, alpha, stabcenterlinkingconss[i], pricingtype->consGetDual(scip_, cons));
   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::consGetDual(
   int i,
   SCIP_Real* dual
   )
{
   SCIP* origprob = GCGpricerGetOrigprob(scip_);
#ifndef NDEBUG
   int nconss =  GCGrelaxGetNMasterConss(origprob);
#endif
   assert(i < nconss);
   assert(dual != NULL);

   SCIP_CONS* cons = GCGrelaxGetMasterConss(origprob)[i];

   if( i >= nstabcenterconss )
      SCIP_CALL( updateStabcenterconss() );

   assert(i < nstabcenterconss);

   *dual = computeDual(stabcenterconss[i], pricingtype->consGetDual(scip_, cons) );
   //SCIPdebugMessage("dual for constraint <%s> is %g (alpha = %g, stab = %g, dual = %g).\n",
   //      SCIPconsGetName(cons), *dual, alpha, stabcenterconss[i], pricingtype->consGetDual(scip_, cons));
   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::rowGetDual(
   int i,
   SCIP_Real* dual
   )
{
#ifndef NDEBUG
   int nrows = GCGsepaGetNMastercuts(scip_);
#endif
   assert(i < nrows);
   assert(dual != NULL);

   SCIP_ROW* row = GCGsepaGetMastercuts(scip_)[i];

   if( i >= nstabcentercuts )
      SCIP_CALL( updateStabcentercuts() );

   assert(i < nstabcentercuts);

   *dual = computeDual(stabcentercuts[i], pricingtype->rowGetDual(row) );
   //SCIPdebugMessage("dual for row <%s> is %g (alpha = %g, stab = %g, dual = %g).\n",
   //      SCIProwGetName(row), *dual, alpha, stabcenterconss[i], pricingtype->rowGetDual(row));

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::updateStabilityCenter(
   SCIP_Real lowerbound
   )
{
   SCIPdebugMessage("Updating stability center: ");

   /* in case the bound is not improving, do nothing */
   if( SCIPisLE(scip_, lowerbound, SCIPnodeGetLowerbound(SCIPgetCurrentNode(scip_))) )
   {
      SCIPdebugPrintf("no bound increase: %g <= %g\n", lowerbound, SCIPnodeGetLowerbound(SCIPgetCurrentNode(scip_)));
      return SCIP_OKAY;
   }

   SCIPdebugPrintf("bound increase: %g >= %g\n", lowerbound, SCIPnodeGetLowerbound(SCIPgetCurrentNode(scip_)));

   /* first update the arrays */
   SCIP_CALL( updateStabcenterconss() );
   SCIP_CALL( updateStabcentercuts() );

   /* get new dual values */
   SCIP* origprob = GCGpricerGetOrigprob(scip_);
   SCIP_CONS** conss = GCGrelaxGetMasterConss(origprob);
   SCIP_CONS** linkingconss = GCGrelaxGetLinkingconss(origprob);
   assert(nstabcenterlinkingconss <= GCGrelaxGetNLinkingconss(origprob) );
   int nconss = GCGrelaxGetNMasterConss(origprob);
   assert(nconss <= nstabcenterconss);
   SCIP_ROW** cuts = GCGsepaGetMastercuts(scip_);
   int ncuts = GCGsepaGetNMastercuts(scip_);
   assert(ncuts <= nstabcentercuts);

   for( int i = 0; i < nconss; ++i )
   {
      stabcenterconss[i] = pricingtype->consGetDual(scip_, conss[i]);
   }

   for( int i = 0; i < ncuts; ++i )
   {
      stabcentercuts[i] = pricingtype->rowGetDual(cuts[i]);
   }

   for( int i = 0; i < nstabcenterlinkingconss; ++i)
   {
      stabcenterlinkingconss[i] = pricingtype->consGetDual(scip_, linkingconss[i]);
   }

   hasstabilitycenter = TRUE;

   return SCIP_OKAY;
}

SCIP_Real Stabilization::computeDual(
      SCIP_Real center,
      SCIP_Real current
      )
{
   if( hasstabilitycenter )
      return alpha*center+(1.0-alpha)*current;
   else
      return current;
}
void Stabilization::updateIterationCount()
{
   if( node != SCIPgetCurrentNode(scip_) )
   {
      node = SCIPgetCurrentNode(scip_);
      k = 1;
      alpha= 0.8;
      hasstabilitycenter = FALSE;
   }
   else
   {
      ++k;
   }
}
void Stabilization::updateAlphaMisprice()
{
   SCIPdebugMessage("Alpha update after mispricing\n");
   updateIterationCount();
   alpha = MAX(0, 1-k*(1-alpha));
   SCIPdebugMessage("alpha updated to %g (k=%d)\n", alpha, k);
}

void Stabilization::updateAlpha(
   SCIP_SOL**            pricingsols         /**< solutions of the pricing problems */
   )
{
   SCIPdebugMessage("Alpha update after successful pricing\n");
   updateIterationCount();

   if( calculateSubgradient(pricingsols) > 0 )
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
   alpha = alpha+(1-alpha)*0.1;
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
   SCIP_SOL**            pricingsols         /**< solutions of the pricing problems */
   )
{
   SCIP* origprob = GCGpricerGetOrigprob(scip_);
   SCIP_CONS** origmasterconss = GCGrelaxGetLinearOrigMasterConss(origprob);
   SCIP_CONS** masterconss = GCGrelaxGetMasterConss(origprob);

   SCIP_CONS** linkingconss = GCGrelaxGetLinkingconss(origprob);
   int nlinkingconss = GCGrelaxGetNLinkingconss(origprob);
   int* linkingconsblocks = GCGrelaxGetLinkingconssBlock(origprob);
   assert(nstabcenterlinkingconss <= GCGrelaxGetNLinkingconss(origprob) );
   int nconss = GCGrelaxGetNMasterConss(origprob);
   assert(nconss <= nstabcenterconss);
   SCIP_ROW** mastercuts = GCGsepaGetMastercuts(scip_);
   SCIP_ROW** origmastercuts = GCGsepaGetOrigcuts(scip_);
   int ncuts = GCGsepaGetNMastercuts(scip_);
   assert(ncuts <= nstabcentercuts);

   SCIP_Real gradientproduct = 0.0;

   /* masterconss*/
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

      SCIP_Real dual = stabcenterconss[i] - pricingtype->consGetDual(scip_, masterconss[i]);

      for( int j = 0; j < nvars; ++j )
      {
         SCIP_Real val = 0.0;
         assert(GCGvarIsOriginal(vars[j]));
         if( GCGvarGetBlock(vars[j]) < 0 )
         {
            SCIP_VAR* mastervar = GCGoriginalVarGetMastervars(vars[j])[0];
            assert(GCGvarIsMaster(mastervar));
            val = SCIPgetSolVal(scip_, NULL, mastervar);

         }
         else
         {
            int block = GCGvarGetBlock(vars[j]);
            if( !GCGrelaxIsPricingprobRelevant(origprob, block) )
               continue;

            SCIP_VAR* pricingvar = GCGoriginalVarGetPricingVar(vars[j]);
            assert(GCGvarIsPricing(pricingvar));
            SCIP* pricingprob = GCGrelaxGetPricingprob(origprob, block);
            assert(pricingprob != NULL);
            val = SCIPgetSolVal(pricingprob, pricingsols[block], pricingvar);
         }
         gradientproduct += dual * vals[j] * val;
      }

      if( SCIPisGT(scip_, dual, 0.0) )
      {
         lhs = SCIPgetLhsLinear(origprob, origcons);
      }
      else if( SCIPisLT(scip_, dual, 0.0) )
      {
         lhs = SCIPgetRhsLinear(origprob, origcons);
      }
      else
      {
         continue;
      }
      gradientproduct -= lhs * dual;
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

      SCIP_Real dual = stabcentercuts[i] - pricingtype->rowGetDual(mastercuts[i]);

      for( int j = 0; j < nvars; ++j )
      {
         SCIP_Real val = 0.0;
         SCIP_VAR* var = SCIPcolGetVar(cols[j]);
         assert(GCGvarIsOriginal(var));
         if( GCGvarGetBlock(var) < 0 )
         {
            SCIP_VAR* mastervar = GCGoriginalVarGetMastervars(var)[0];
            assert(GCGvarIsMaster(mastervar));
            val = SCIPgetSolVal(scip_, NULL, mastervar);

         }
         else
         {
            int block = GCGvarGetBlock(var);
            if( !GCGrelaxIsPricingprobRelevant(origprob, block) )
               continue;

            SCIP_VAR* pricingvar = GCGoriginalVarGetPricingVar(var);
            assert(GCGvarIsPricing(pricingvar));
            SCIP* pricingprob = GCGrelaxGetPricingprob(origprob, block);
            assert(pricingprob != NULL);
            val = SCIPgetSolVal(pricingprob, pricingsols[block], pricingvar);
         }
         gradientproduct += dual * vals[j] * val;
      }

      if( SCIPisGT(scip_, dual, 0.0) )
      {
         lhs = SCIProwGetLhs(origcut);
      }
      else if( SCIPisLT(scip_, dual, 0.0) )
      {
         lhs = SCIProwGetRhs(origcut);
      }
      else
      {
         continue;
      }
      gradientproduct -= lhs * dual;
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
      SCIP* pricingprob = GCGrelaxGetPricingprob(origprob, block);
      assert(pricingprob != NULL);

      SCIP_Real dual = stabcenterlinkingconss[i] - pricingtype->consGetDual(scip_, linkingcons);
      SCIP_Real masterval = SCIPgetSolVal(scip_, NULL, mastervar);
      SCIP_Real pricingval = SCIPgetSolVal(pricingprob, pricingsols[block], pricingvar);
      gradientproduct += dual * (masterval - pricingval);
   }

   SCIPdebugMessage("Update gradient with value %g.\n", gradientproduct);

   return gradientproduct;
}

SCIP_Bool Stabilization::isStabilized()
{
   return SCIPisGT(scip_, alpha, 0.0);
}

} /* namespace gcg */
