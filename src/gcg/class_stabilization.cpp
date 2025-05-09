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
#include "gcg/class_stabilization.h"
#include "gcg/pub_extendedmasterconsdata.h"
#include "gcg/pricer_gcg.h"
#include "gcg/gcg.h"
#include "gcg/pub_gcgcol.h"
#include "gcg/pub_gcgvar.h"
#include "gcg/sepa_original.h"
#include "objscip/objscip.h"
#include "scip/cons_linear.h"
#include "gcg/scip_misc.h"
#include "scip/pub_misc.h"
#include "scip/type_misc.h"
#include "gcg/type_extendedmasterconsdata.h"
#include <cstddef>

namespace gcg {

Stabilization::Stabilization(
   GCG* gcgstruct,
   PricingType* pricingtype_,
   SCIP_Bool hybridascent_
   ) : gcg(gcgstruct), masterprob(GCGgetMasterprob(gcg)), stabcenterconsvals((SCIP_Real*) NULL), stabcenterconsvalssize(0), nstabcenterconsvals(0),
      stabcenteroriginalsepacutvals((SCIP_Real*) NULL), stabcenteroriginalsepacutvalssize(0), nstabcenteroriginalsepacutvals(0),
      stabcenterextendedmasterconss((GCG_EXTENDEDMASTERCONSDATA**) NULL), nstabcenterextendedmasterconss(0), stabcenterextendedmasterconsssize(0), stabcenterextendedmasterconsvals((SCIP_Real*) NULL),
      stabcenterlinkingconsvals((SCIP_Real*) NULL), nstabcenterlinkingconsvals(0), stabcenterlinkingconsvalssize(0),
      stabcenterconv((SCIP_Real*) NULL), nstabcenterconv(0), dualdiffnorm(0.0),
      subgradientconsvals(NULL), subgradientconsvalssize(0), nsubgradientconsvals(0),
      subgradientoriginalsepacutvals(NULL), subgradientoriginalsepacutvalssize(0), nsubgradientoriginalsepacutvals(0),
      subgradientextendedmasterconss((GCG_EXTENDEDMASTERCONSDATA**) NULL), nsubgradientextendedmasterconss(0), subgradientextendedmasterconsssize(0), subgradientextendedmasterconsvals(NULL),
      subgradientlinkingconsvals(NULL), subgradientlinkingconsvalssize(0),
      subgradientnorm(0.0), hybridfactor(0.0),
      pricingtype(pricingtype_), alpha(0.8), alphabar(0.8), hybridascent(hybridascent_), beta(0.0), nodenr(-1), k(0), t(0), hasstabilitycenter(FALSE),stabcenterbound(-SCIPinfinity(GCGgetMasterprob(gcg))),
      inmispricingschedule(FALSE), subgradientproduct(0.0)
{
}

Stabilization::~Stabilization()
{
   SCIPfreeBlockMemoryArrayNull(masterprob, &stabcenterconsvals, stabcenterconsvalssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(masterprob, &stabcenteroriginalsepacutvals, stabcenteroriginalsepacutvalssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(masterprob, &stabcenterlinkingconsvals, stabcenterlinkingconsvalssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(masterprob, &subgradientconsvals, subgradientconsvalssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(masterprob, &subgradientoriginalsepacutvals, subgradientoriginalsepacutvalssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(masterprob, &subgradientlinkingconsvals, subgradientlinkingconsvalssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(masterprob, &stabcenterconv, nstabcenterconv); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(masterprob, &stabcenterextendedmasterconss, stabcenterextendedmasterconsssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(masterprob, &stabcenterextendedmasterconsvals, stabcenterextendedmasterconsssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(masterprob, &subgradientextendedmasterconsvals, subgradientextendedmasterconsssize); /*lint !e64*/
   SCIPfreeBlockMemoryArrayNull(masterprob, &subgradientextendedmasterconss, subgradientextendedmasterconsssize); /*lint !e64*/
   masterprob = (SCIP*) NULL;
   stabcenterconsvals = (SCIP_Real*) NULL;
   stabcenteroriginalsepacutvals = (SCIP_Real*) NULL;
   stabcenterlinkingconsvals = (SCIP_Real*) NULL;
   stabcenterconv = (SCIP_Real*) NULL;
   pricingtype = (PricingType*) NULL;
   nodenr = -1;
}

SCIP_RETCODE Stabilization::updateStabcenterconsvals()
{
   int nconss = GCGgetNMasterConss(gcg);

   if( nconss == nstabcenterconsvals )
   {
      return SCIP_OKAY;
   }

   if( nconss > stabcenterconsvalssize )
   {
      int oldsize = stabcenterconsvalssize;
      stabcenterconsvalssize = SCIPcalcMemGrowSize(masterprob, nconss);
      SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &stabcenterconsvals, oldsize, stabcenterconsvalssize) );
   }
   assert(stabcenterconsvals != NULL);
   BMSclearMemoryArray(&stabcenterconsvals[nstabcenterconsvals], (size_t)nconss-nstabcenterconsvals); /*lint !e866*/

   nstabcenterconsvals = nconss;

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::updateStabcenteroriginalcutvals()
{
   int noriginalsepacuts = GCGsepaGetNOriginalSepaCuts(gcg);

   if( noriginalsepacuts == nstabcenteroriginalsepacutvals )
   {
      return SCIP_OKAY;
   }

   if( noriginalsepacuts > stabcenteroriginalsepacutvalssize )
   {
      int oldsize = stabcenteroriginalsepacutvalssize;
      stabcenteroriginalsepacutvalssize = SCIPcalcMemGrowSize(masterprob, noriginalsepacuts);
      SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &stabcenteroriginalsepacutvals, oldsize, stabcenteroriginalsepacutvalssize) );
   }
   assert(stabcenteroriginalsepacutvals != NULL);
   BMSclearMemoryArray(&stabcenteroriginalsepacutvals[nstabcenteroriginalsepacutvals], (size_t)noriginalsepacuts - nstabcenteroriginalsepacutvals); /*lint !e866*/

   nstabcenteroriginalsepacutvals = noriginalsepacuts;

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::updateStabcenterextendedmasterconsvals()
{
   GCG_BRANCHRULE** branchrules;
   GCG_BRANCHDATA** branchdata;
   GCG_EXTENDEDMASTERCONSDATA** branchextendedmasterconsdata;
   int nbranchextendedmasterconss;
   int i;

   SCIP_CALL( GCGrelaxBranchGetAllActiveExtendedMasterConss(gcg, &branchrules, &branchdata, &branchextendedmasterconsdata, &nbranchextendedmasterconss) );

   // grow if necessary
   if( nbranchextendedmasterconss > stabcenterextendedmasterconsssize )
   {
      int oldsize = stabcenterextendedmasterconsssize;
      stabcenterextendedmasterconsssize = SCIPcalcMemGrowSize(masterprob, nbranchextendedmasterconss);
      SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &stabcenterextendedmasterconss, oldsize, stabcenterextendedmasterconsssize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &stabcenterextendedmasterconsvals, oldsize, stabcenterextendedmasterconsssize) );
   }

   // update the arrays (extended master cons could have changed, even if size is the same)
   for( i = 0; i < nbranchextendedmasterconss; i++ )
   {
      stabcenterextendedmasterconss[i] = branchextendedmasterconsdata[i];
   }
   nstabcenterextendedmasterconss = nbranchextendedmasterconss;

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::updateSubgradientconsvals()
{
   int nconss = GCGgetNMasterConss(gcg);

   if( nconss == nsubgradientconsvals )
   {
      return SCIP_OKAY;
   }

   if( nconss > subgradientconsvalssize )
   {
      int oldsize = subgradientconsvalssize;
      subgradientconsvalssize = SCIPcalcMemGrowSize(masterprob, nconss);
      SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &subgradientconsvals, oldsize, subgradientconsvalssize) );
   }
   assert(subgradientconsvals != NULL);
   BMSclearMemoryArray(&subgradientconsvals[nsubgradientconsvals], (size_t)nconss-nsubgradientconsvals); /*lint !e866*/

   nsubgradientconsvals = nconss;

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::updateSubgradientoriginalcutvals()
{
   int noriginalsepacuts = GCGsepaGetNOriginalSepaCuts(gcg);

   if( noriginalsepacuts == nsubgradientoriginalsepacutvals )
   {
      return SCIP_OKAY;
   }

   if( noriginalsepacuts > subgradientoriginalsepacutvalssize )
   {
      int oldsize = subgradientoriginalsepacutvalssize;
      subgradientoriginalsepacutvalssize = SCIPcalcMemGrowSize(masterprob, noriginalsepacuts);
      SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &subgradientoriginalsepacutvals, oldsize, subgradientoriginalsepacutvalssize) );
   }
   assert(subgradientoriginalsepacutvals != NULL);
   BMSclearMemoryArray(&subgradientoriginalsepacutvals[nsubgradientoriginalsepacutvals], (size_t)noriginalsepacuts - nsubgradientoriginalsepacutvals); /*lint !e866*/

   nsubgradientoriginalsepacutvals = noriginalsepacuts;

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::updateSubgradientextendedmasterconsvals()
{
   GCG_BRANCHRULE** branchrules;
   GCG_BRANCHDATA** branchdata;
   GCG_EXTENDEDMASTERCONSDATA** branchextendedmasterconsdata;
   int nbranchextendedmasterconss;
   int i;

   SCIP_CALL( GCGrelaxBranchGetAllActiveExtendedMasterConss(gcg, &branchrules, &branchdata, &branchextendedmasterconsdata, &nbranchextendedmasterconss) );

   // grow if necessary
   if( nbranchextendedmasterconss > subgradientextendedmasterconsssize )
   {
      int oldsize = subgradientextendedmasterconsssize;
      subgradientextendedmasterconsssize = SCIPcalcMemGrowSize(masterprob, nbranchextendedmasterconss);
      SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &subgradientextendedmasterconss, oldsize, subgradientextendedmasterconsssize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &subgradientextendedmasterconsvals, oldsize, subgradientextendedmasterconsssize) );
   }

   // update the arrays (extended master conss could have changed, even if size is the same)
   for( i = 0; i < nbranchextendedmasterconss; i++ )
   {
      subgradientextendedmasterconss[i] = branchextendedmasterconsdata[i];
   }
   nsubgradientextendedmasterconss = nbranchextendedmasterconss;

   return SCIP_OKAY;
}


SCIP_RETCODE Stabilization::setNLinkingconsvals(
   int nlinkingconssnew
   )
{
   if( nlinkingconssnew > stabcenterlinkingconsvalssize)
   {
      int newsize = SCIPcalcMemGrowSize(masterprob, nlinkingconssnew);
      SCIP_CALL(SCIPreallocBlockMemoryArray(masterprob, &stabcenterlinkingconsvals, stabcenterlinkingconsvalssize, newsize));


      if( hybridascent )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(masterprob, &subgradientlinkingconsvals, stabcenterlinkingconsvalssize,
            newsize) );
         BMSclearMemoryArray(subgradientlinkingconsvals, nlinkingconssnew);
      }
      stabcenterlinkingconsvalssize = newsize;
   }

   nstabcenterlinkingconsvals = nlinkingconssnew;
   BMSclearMemoryArray(stabcenterlinkingconsvals, nstabcenterlinkingconsvals);


   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::setNConvconsvals(
   int nconvconssnew
   )
{
   SCIPfreeBlockMemoryArrayNull(masterprob, &stabcenterconv, nstabcenterconv); /*lint !e64*/
   SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &stabcenterconv, nconvconssnew) );

   nstabcenterconv = nconvconssnew;
   BMSclearMemoryArray(stabcenterconv, nstabcenterconv);

   return SCIP_OKAY;
}

SCIP_Real Stabilization::linkingconsGetDual(
   int i
   )
{
   SCIP_Real subgradient = 0.0;

   assert(i < nstabcenterlinkingconsvals);
   assert(nstabcenterlinkingconsvals<= GCGgetNVarLinkingconss(gcg));
   assert(stabcenterlinkingconsvals != NULL);

   SCIP_CONS* cons = GCGgetVarLinkingconss(gcg)[i];

   if( hybridascent && hasstabilitycenter )
      subgradient = subgradientlinkingconsvals[i];

   return computeDual(
      stabcenterlinkingconsvals[i],
      pricingtype->consGetDual(cons),
      subgradient,
      0.0,
      0.0
   );
}

SCIP_RETCODE Stabilization::consGetDual(
   int                   i,                  /* index of the constraint */
   SCIP_Real*            dual                /* return pointer for dual value */
   )
{
   SCIP_Real subgradient = 0.0;
#ifndef NDEBUG
   int nconss =  GCGgetNMasterConss(gcg);
#endif
   assert(i < nconss);
   assert(dual != NULL);

   SCIP_CONS* cons = GCGgetMasterConss(gcg)[i];

   if( i >= nstabcenterconsvals )
      SCIP_CALL( updateStabcenterconsvals() );

   assert(i < nstabcenterconsvals);
   assert(stabcenterconsvals != NULL);

   if( i >= nsubgradientconsvals && hybridascent )
      SCIP_CALL( updateSubgradientconsvals() );

   if( hybridascent && hasstabilitycenter )
      subgradient = subgradientconsvals[i];

   *dual = computeDual(
      stabcenterconsvals[i],
      pricingtype->consGetDual(cons),
      subgradient,
      SCIPgetLhsLinear(masterprob, cons),
      SCIPgetRhsLinear(masterprob, cons)
   );
   return SCIP_OKAY;

}

SCIP_RETCODE Stabilization::rowGetDual(
   int                   i,                  /* index of the row */
   SCIP_Real*            dual                /* return pointer for dual value */
   )
{
#ifndef NDEBUG
   int noriginalsepacuts = GCGsepaGetNOriginalSepaCuts(gcg);
#endif
   assert(i < noriginalsepacuts);
   assert(dual != NULL);

   SCIP_ROW* originalsepacut = GCGsepaGetOriginalSepaMastercuts(gcg)[i];
   SCIP_Real subgradient = 0.0;

   if( i >= nstabcenteroriginalsepacutvals )
      SCIP_CALL( updateStabcenteroriginalcutvals() );

   assert(i < nstabcenteroriginalsepacutvals);
   assert(stabcenteroriginalsepacutvals != NULL);

   if( i >= nsubgradientoriginalsepacutvals && hybridascent )
      SCIP_CALL( updateSubgradientoriginalcutvals() );

   if( hybridascent && hasstabilitycenter )
   {
      assert(subgradientoriginalsepacutvals != NULL);
      subgradient = subgradientoriginalsepacutvals[i];
   }

   *dual = computeDual(
      stabcenteroriginalsepacutvals[i],
      pricingtype->rowGetDual(originalsepacut),
      subgradient,
      SCIProwGetLhs(originalsepacut),
      SCIProwGetRhs(originalsepacut)
   );

   return SCIP_OKAY;
}

SCIP_Real Stabilization::convGetDual(
   int i
   )
{
   assert(i < nstabcenterconv);
   assert(nstabcenterconv<= GCGgetNPricingprobs(gcg));
   assert(stabcenterconv != NULL);

   SCIP_CONS* cons = GCGgetConvCons(gcg, i);
   SCIP_Real subgradient = 0.0;

   return computeDual(
      stabcenterconv[i],
      pricingtype->consGetDual(cons),
      subgradient,
      (SCIP_Real) GCGgetNIdenticalBlocks(gcg, i),
      (SCIP_Real) GCGgetNIdenticalBlocks(gcg, i)
   );
}

SCIP_RETCODE Stabilization::extendedmasterconsGetDual(
   GCG_EXTENDEDMASTERCONSDATA*    extendedmasterconsdata,      /**< extendedmasterconsdata */
   SCIP_Real*            dual                /**< return pointer for dual value */
   )
{
   assert(extendedmasterconsdata != NULL);
   assert(dual != NULL);


   SCIP_CALL( updateStabcenterextendedmasterconsvals() );

   if( hybridascent )
      SCIP_CALL( updateSubgradientextendedmasterconsvals() );

   SCIP_Real stabcenter = 0.0;
#ifndef NDEBUG
   SCIP_Bool found = FALSE;
#endif
   for( int i = 0; i < nstabcenterextendedmasterconss; i++ )
   {
      if( stabcenterextendedmasterconss[i] == extendedmasterconsdata )
      {
         stabcenter = stabcenterextendedmasterconsvals[i];
#ifndef NDEBUG
         found = TRUE;
#endif
         break;
      }
   }
#ifndef NDEBUG
   assert(found);
#endif

   *dual = computeDual(
      stabcenter,
      pricingtype->extendedmasterconsGetDual(extendedmasterconsdata),
      0.0,
      GCGextendedmasterconsGetLhs(gcg, extendedmasterconsdata),
      GCGextendedmasterconsGetRhs(gcg, extendedmasterconsdata)
   );

   return SCIP_OKAY;
}

SCIP_RETCODE Stabilization::updateStabilityCenter(
   SCIP_Real             lowerbound,         /**< lower bound due to lagrange function corresponding to current (stabilized) dual vars */
   SCIP_Real*            dualsolconv,        /**< corresponding feasible dual solution for convexity constraints */
   GCG_COL**             pricingcols         /**< columns of the pricing problems */
   )
{
   SCIP_Real dual;
   int i;
   GCG_EXTENDEDMASTERCONSDATA* tmpextendedmasterconsdata;

   assert(dualsolconv != NULL);
   SCIPdebugMessage("Updating stability center: ");

   /* in case the bound is not improving and we have a stability center, do nothing */
   if( SCIPisLE(masterprob, lowerbound, stabcenterbound) && hasstabilitycenter )
   {
      SCIPdebugPrintf("no bound increase: %g <= %g\n", lowerbound, SCIPnodeGetLowerbound(SCIPgetCurrentNode(masterprob)));
      return SCIP_OKAY;
   }

   SCIPdebugPrintf("bound increase: %g > %g\n", lowerbound, SCIPnodeGetLowerbound(SCIPgetCurrentNode(masterprob)));

   /* first update the arrays */
   SCIP_CALL( updateStabcenterconsvals() );
   SCIP_CALL( updateStabcenteroriginalcutvals() );
   SCIP_CALL( updateStabcenterextendedmasterconsvals() );

   if( hybridascent )
   {
      SCIP_CALL( updateSubgradientconsvals() );
      SCIP_CALL( updateSubgradientoriginalcutvals() );
      SCIP_CALL( updateSubgradientextendedmasterconsvals() );
   }

   /* get new dual values */
   int nconss = GCGgetNMasterConss(gcg);
   int noriginalsepacuts = GCGsepaGetNOriginalSepaCuts(gcg);
   int nprobs = GCGgetNPricingprobs(gcg);

   assert(nstabcenterlinkingconsvals <= GCGgetNVarLinkingconss(gcg) );
   assert(nconss <= nstabcenterconsvals);
   assert(noriginalsepacuts <= nstabcenteroriginalsepacutvals);

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( consGetDual(i, &stabcenterconsvals[i]) );
   }

   for( i = 0; i < noriginalsepacuts; ++i )
   {
      SCIP_CALL( rowGetDual(i, &stabcenteroriginalsepacutvals[i]) );
   }

   for( i = 0; i < nstabcenterlinkingconsvals; ++i)
   {
      stabcenterlinkingconsvals[i] = linkingconsGetDual(i);
   }

   for( i = 0; i < nprobs; ++i )
   {
      if( !GCGisPricingprobRelevant(gcg, i) )
         continue;
      stabcenterconv[i] = dualsolconv[i];
   }

   for( i = 0; i < nstabcenterextendedmasterconss; ++i )
   {
      tmpextendedmasterconsdata = stabcenterextendedmasterconss[i];
      assert(tmpextendedmasterconsdata != NULL);
      SCIP_CALL( extendedmasterconsGetDual(tmpextendedmasterconsdata, &dual) );
      stabcenterextendedmasterconsvals[i] = dual;
   }

   if( hybridascent )
      SCIP_CALL( calculateSubgradient(pricingcols) );

   hasstabilitycenter = TRUE;
   stabcenterbound = lowerbound;

   return SCIP_OKAY;
}

SCIP_Real Stabilization::computeDual(
      SCIP_Real         center,
      SCIP_Real         current,
      SCIP_Real         subgradient,         /**< subgradient (or 0.0 if not needed) */
      SCIP_Real         lhs,                 /**< lhs (or 0.0 if not needed) */
      SCIP_Real         rhs                 /**< rhs (or 0.0 if not needed) */
      ) const
{
   SCIP_Real usedalpha = alpha;
   SCIP_Real usedbeta = beta;

   if ( inmispricingschedule )
   {
      usedalpha = alphabar;
      usedbeta = 0.0;
   }

   if( hasstabilitycenter && (SCIPisZero(masterprob, usedbeta) || SCIPisZero(masterprob, usedalpha)) )
      return usedalpha*center+(1.0-usedalpha)*current;
   else if( hasstabilitycenter && SCIPisPositive(masterprob, usedbeta) )
   {
      SCIP_Real dual = center + hybridfactor * (beta * (center + subgradient * dualdiffnorm / subgradientnorm) + (1.0 - beta) * current - center);

      /* make sure dual solution has the correct sign */
      if( SCIPisInfinity(masterprob, rhs) )
         dual = MAX(dual, 0.0);
      else if( SCIPisInfinity(masterprob, -lhs) )
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
   if( nodenr != SCIPnodeGetNumber(SCIPgetCurrentNode(masterprob)) )
   {
      nodenr = SCIPnodeGetNumber(SCIPgetCurrentNode(masterprob));
      k = 0;
      t = 1;
      alpha = 0.8;
      hasstabilitycenter = FALSE;
      stabcenterbound = -SCIPinfinity(masterprob);
      inmispricingschedule = FALSE;
   }
}

/**< update information for hybrid stabilization with dual ascent */
SCIP_RETCODE Stabilization::updateHybrid()
{
   if( hasstabilitycenter && hybridascent && !inmispricingschedule )
   {
      /* first update the arrays */
      SCIP_CALL( updateStabcenterconsvals() );
      SCIP_CALL( updateStabcenteroriginalcutvals() );
      SCIP_CALL( updateStabcenterextendedmasterconsvals() );

      SCIP_CALL( updateSubgradientconsvals() );
      SCIP_CALL( updateSubgradientoriginalcutvals() );
      SCIP_CALL( updateSubgradientextendedmasterconsvals() );

      if( SCIPisPositive(masterprob, alpha) )
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
   alphabar = MAX(0.0, 1-k*(1-alpha));
   SCIPdebugMessage("alphabar updated to %g in mispricing iteration k=%d and node pricing iteration t=%d \n", alphabar, k, t);
}

void Stabilization::updateAlpha()
{
   SCIPdebugMessage("Alpha update after successful pricing\n");
   updateIterationCount();

   /* There is a sign error in the stabilization paper:
    * if the scalar product (subgradientproduct) is positive, the angle is less than 90° and we want to decrease alpha
    */
   if( SCIPisNegative(masterprob, subgradientproduct) )
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
   alpha = MAX(0.0, alpha-0.1);

   SCIPdebugMessage("alpha decreased to %g\n", alpha);
}

SCIP_Real Stabilization::calculateSubgradientProduct(
   GCG_COL**            pricingcols         /**< solutions of the pricing problems */
   )
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   SCIP_CONS** origmasterconss = GCGgetOrigMasterConss(gcg);
   SCIP_CONS** masterconss = GCGgetMasterConss(gcg);
   SCIP_CONS** linkingconss = GCGgetVarLinkingconss(gcg);
   int nlinkingconss = GCGgetNVarLinkingconss(gcg);
   int* linkingconsblocks = GCGgetVarLinkingconssBlock(gcg);
   assert(nstabcenterlinkingconsvals <= GCGgetNVarLinkingconss(gcg) );
   int nconss = GCGgetNMasterConss(gcg);
   assert(nconss <= nstabcenterconsvals);
   SCIP_ROW** originalsepamastercuts = GCGsepaGetOriginalSepaMastercuts(gcg);
   SCIP_ROW** originalsepaorigcuts = GCGsepaGetOriginalSepaOrigcuts(gcg);
   int noriginalsepacuts = GCGsepaGetNOriginalSepaCuts(gcg);
   assert(noriginalsepacuts <= nstabcenteroriginalsepacutvals);

   SCIP_Real gradientproduct = 0.0;

   /* masterconss */
   for( int i = 0; i < nconss; ++i )
   {
      SCIP_VAR** vars = NULL;
      SCIP_Real* vals = NULL;
      int nvars;
      SCIP_Real lhs; /* can also be rhs, but we need only one */

      SCIP_CONS* origcons = origmasterconss[i];
      nvars = GCGconsGetNVars(origprob, origcons);
      SCIP_CALL( SCIPallocBufferArray(origprob, &vars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(origprob, &vals, nvars) );
      GCGconsGetVars(origprob, origcons, vars, nvars);
      GCGconsGetVals(origprob, origcons, vals, nvars);

      SCIP_Real dual =  pricingtype->consGetDual(masterconss[i]);
      SCIP_Real stabdual;

      SCIP_CALL( consGetDual(i, &stabdual) );

      assert(!SCIPisInfinity(masterprob, ABS(dual)));

      if( SCIPisFeasPositive(masterprob, stabdual) )
      {
         lhs = GCGconsGetLhs(origprob, origcons);
      }
      else if( SCIPisFeasNegative(masterprob, stabdual) )
      {
         lhs = GCGconsGetRhs(origprob, origcons);
      }
      else
      {
         SCIPfreeBufferArray(origprob, &vals);
         SCIPfreeBufferArray(origprob, &vars);
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
            val = SCIPgetSolVal(masterprob, (SCIP_SOL*) NULL, mastervar);
            assert( !SCIPisInfinity(masterprob, val) );
         }
         else
         {
            int block = GCGvarGetBlock(vars[j]);
            if( !GCGisPricingprobRelevant(gcg, block) )
               continue;

            assert(pricingcols[block] != NULL);

            SCIP_VAR* pricingvar = GCGoriginalVarGetPricingVar(vars[j]);
            assert(GCGvarIsPricing(pricingvar));
            val = GCGcolGetSolVal(pricingcols[block], pricingvar);
            assert(!SCIPisInfinity(masterprob, ABS(val)));
         }
         assert(stabcenterconsvals != NULL);
         assert(vals != NULL);
         gradientproduct -= (dual - stabcenterconsvals[i]) * vals[j] * val;
      }

      assert(stabcenterconsvals != NULL);
      assert(!SCIPisInfinity(masterprob, ABS(lhs)));

      gradientproduct += (dual - stabcenterconsvals[i]) * lhs;
      SCIPfreeBufferArray(origprob, &vals);
      SCIPfreeBufferArray(origprob, &vars);
   }

   /* originalcuts */
   for( int i = 0; i < noriginalsepacuts; ++i )
   {
      SCIP_COL** cols;
      SCIP_Real* vals;
      int nvars;
      SCIP_Real lhs; /* can also be rhs, but we need only one */

      SCIP_ROW* originalsepaorigcut = originalsepaorigcuts[i];
      nvars = SCIProwGetNNonz(originalsepaorigcut);
      cols = SCIProwGetCols(originalsepaorigcut);
      vals = SCIProwGetVals(originalsepaorigcut);

      SCIP_Real dual = pricingtype->rowGetDual(originalsepamastercuts[i]);
      assert(!SCIPisInfinity(masterprob, ABS(dual)));

      SCIP_Real stabdual;

      SCIP_CALL( rowGetDual(i, &stabdual) );

      if( SCIPisFeasGT(masterprob, stabdual, 0.0) )
      {
         lhs = SCIProwGetLhs(originalsepaorigcut);
      }
      else if( SCIPisFeasLT(masterprob, stabdual, 0.0) )
      {
         lhs = SCIProwGetRhs(originalsepaorigcut);
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
            val = SCIPgetSolVal(masterprob, (SCIP_SOL*) NULL, mastervar);
            assert(!SCIPisInfinity(masterprob, ABS(val)));
         }
         /* Variable in a pricing problem */
         else
         {
            int block = GCGvarGetBlock(var);
            if( !GCGisPricingprobRelevant(gcg, block) )
               continue;

            assert(pricingcols[block] != NULL);

            SCIP_VAR* pricingvar = GCGoriginalVarGetPricingVar(var);
            assert(GCGvarIsPricing(pricingvar));
            val = GCGcolGetSolVal(pricingcols[block], pricingvar);
            assert(!SCIPisInfinity(masterprob, ABS(val)));
         }
         assert(stabcenteroriginalsepacutvals != NULL);
         assert(vals != NULL);
         gradientproduct -= (dual - stabcenteroriginalsepacutvals[i]) * vals[j] * val;
      }

      assert(!SCIPisInfinity(masterprob, ABS(lhs)));
      assert(stabcenteroriginalsepacutvals != NULL);

      gradientproduct +=  (dual - stabcenteroriginalsepacutvals[i]) * lhs;
   }

   /* extended master conss */
   GCG_PRICINGMODIFICATION** pricingmods;
   int block;
   for( int i = 0; i < nstabcenterextendedmasterconss; ++i )
   {
      SCIP_COL** cols;
      SCIP_Real* vals;
      int nvars;
      SCIP_Real val;
      SCIP_VAR* var;
      SCIP_Real lhs; /* can also be rhs, but we need only one */
      GCG_EXTENDEDMASTERCONSDATA* tmpextendedmasterconsdata;

      tmpextendedmasterconsdata = stabcenterextendedmasterconss[i];

      nvars = GCGextendedmasterconsGetNNonz(gcg, tmpextendedmasterconsdata);
      cols = GCGextendedmasterconsGetCols(gcg, tmpextendedmasterconsdata);
      vals = GCGextendedmasterconsGetVals(gcg, tmpextendedmasterconsdata);

      SCIP_Real dual = pricingtype->extendedmasterconsGetDual(tmpextendedmasterconsdata);
      assert(!SCIPisInfinity(masterprob, ABS(dual)));

      SCIP_Real stabdual;

      SCIP_CALL( extendedmasterconsGetDual(tmpextendedmasterconsdata, &stabdual) );

      if( SCIPisFeasGT(masterprob, stabdual, 0.0) )
      {
         lhs = GCGextendedmasterconsGetLhs(gcg, tmpextendedmasterconsdata);
      }
      else if( SCIPisFeasLT(masterprob, stabdual, 0.0) )
      {
         lhs = GCGextendedmasterconsGetRhs(gcg, tmpextendedmasterconsdata);
      }
      else
      {
         continue;
      }

      for( int j = 0; j < nvars; ++j )
      {
         val = 0.0;
         var = SCIPcolGetVar(cols[j]);
         assert(GCGvarIsMaster(var));

         /* only linking or static master variable */
         if( GCGvarGetBlock(var) >= 0 )
            continue;

         val = SCIPgetSolVal(masterprob, (SCIP_SOL*) NULL, var);
         assert(!SCIPisInfinity(masterprob, ABS(val)));
         val = val * vals[j];

         assert(vals != NULL);
         gradientproduct -= (dual - stabcenterextendedmasterconsvals[i]) * val;
      }

      pricingmods = GCGextendedmasterconsGetPricingModifications(tmpextendedmasterconsdata);
      for( int j = 0; j < GCGextendedmasterconsGetNPricingModifications(tmpextendedmasterconsdata); j++ )
      {
         block = GCGpricingmodificationGetBlock(pricingmods[j]);
         assert(pricingmods[j] != NULL);
         assert(GCGisPricingprobRelevant(gcg, block));
         assert(pricingcols[block] != NULL);

         SCIP_VAR* pricingvar = GCGpricingmodificationGetCoefVar(pricingmods[j]);
         assert(GCGvarIsInferredPricing(pricingvar));
         val = GCGcolGetSolVal(pricingcols[block], pricingvar);
         assert(!SCIPisInfinity(masterprob, ABS(val)));

         gradientproduct -= (dual - stabcenterextendedmasterconsvals[i]) * val;
      }

      assert(!SCIPisInfinity(masterprob, ABS(lhs)));

      gradientproduct += (dual - stabcenterextendedmasterconsvals[i]) * lhs;
   }

   /* linkingconss */
   for( int i = 0; i < nlinkingconss; ++i )
   {
      SCIP_VAR* mastervar;
      SCIP_VAR* pricingvar;
      SCIP_CONS* linkingcons = linkingconss[i];
      mastervar = SCIPgetVarsLinear(masterprob, linkingcons)[0];
      assert(GCGvarIsMaster(mastervar));

      block = linkingconsblocks[i];
      pricingvar = GCGlinkingVarGetPricingVars(GCGmasterVarGetOrigvars(mastervar)[0])[block];
      assert(GCGvarIsPricing(pricingvar));

      assert(stabcenterlinkingconsvals != NULL);
      SCIP_Real dual = pricingtype->consGetDual(linkingcons) - stabcenterlinkingconsvals[i];

      SCIP_Real stabdual = linkingconsGetDual(i);

      if( SCIPisFeasZero(origprob, stabdual) )
         continue;

      assert(pricingcols[block] != NULL);

      SCIP_Real masterval = SCIPgetSolVal(masterprob, (SCIP_SOL*) NULL, mastervar);
      SCIP_Real pricingval = GCGcolGetSolVal(pricingcols[block], pricingvar);
      assert(!SCIPisInfinity(masterprob, ABS(masterval)));
      assert(!SCIPisInfinity(masterprob, ABS(pricingval)));
      assert(!SCIPisInfinity(masterprob, ABS(dual)));
      gradientproduct -= dual * (masterval - pricingval);
   }

   SCIPdebugMessage("Update gradient product with value %g.\n", gradientproduct);

   return gradientproduct;
}

/** calculates the subgradient (with linking variables) */
SCIP_RETCODE Stabilization::calculateSubgradient(
   GCG_COL**            pricingcols         /**< columns of the pricing problems */
   )
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   SCIP_CONS** origmasterconss = GCGgetOrigMasterConss(gcg);

   SCIP_CONS** linkingconss = GCGgetVarLinkingconss(gcg);
   int nlinkingconss = GCGgetNVarLinkingconss(gcg);
   int* linkingconsblocks = GCGgetVarLinkingconssBlock(gcg);
   assert(nstabcenterlinkingconsvals <= GCGgetNVarLinkingconss(gcg) );
   int nconss = GCGgetNMasterConss(gcg);
   assert(nconss <= nstabcenterconsvals);
   SCIP_ROW** originalsepaorigcuts = GCGsepaGetOriginalSepaOrigcuts(gcg);
   int noriginalsepacuts = GCGsepaGetNOriginalSepaCuts(gcg);
   assert(noriginalsepacuts <= nstabcenteroriginalsepacutvals);

   subgradientnorm = 0.0;

   /* masterconss */
   for( int i = 0; i < nconss; ++i )
   {
      SCIP_VAR** vars = NULL;
      SCIP_Real* vals = NULL;
      int nvars;
      SCIP_Real activity;
      SCIP_Real infeasibility;

      SCIP_CONS* origcons = origmasterconss[i];
      nvars = GCGconsGetNVars(origprob, origcons);
      SCIP_CALL( SCIPallocBufferArray(origprob, &vars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(origprob, &vals, nvars) );
      GCGconsGetVars(origprob, origcons, vars, nvars);
      GCGconsGetVals(origprob, origcons, vals, nvars);

      SCIP_Real dual = stabcenterconsvals[i];
      assert(!SCIPisInfinity(masterprob, ABS(dual)));

      activity = 0.0;

      for( int j = 0; j < nvars; ++j )
      {
         SCIP_Real val = 0.0;
         assert(GCGvarIsOriginal(vars[j]));
         if( GCGvarGetBlock(vars[j]) < 0 )
         {
            SCIP_VAR* mastervar = GCGoriginalVarGetMastervars(vars[j])[0];
            assert(GCGvarIsMaster(mastervar));
            val = SCIPgetSolVal(masterprob, (SCIP_SOL*) NULL, mastervar);
            assert( !SCIPisInfinity(masterprob, val) );
         }
         else
         {
            int block = GCGvarGetBlock(vars[j]);
            if( !GCGisPricingprobRelevant(gcg, block) )
               continue;

            assert(pricingcols[block] != NULL);

            SCIP_VAR* pricingvar = GCGoriginalVarGetPricingVar(vars[j]);
            assert(GCGvarIsPricing(pricingvar));
            val = GCGcolGetSolVal(pricingcols[block], pricingvar);
            assert(!SCIPisInfinity(masterprob, ABS(val)));
         }
         assert(vals != NULL);
         activity += vals[j] * val;
      }

      infeasibility = 0.0;

      if( SCIPisFeasPositive(masterprob, dual) /* || SCIPisInfinity(origprob, SCIPgetRhsLinear(origprob, origcons)) */)
      {
         infeasibility = GCGconsGetLhs(origprob, origcons) - activity;
      }
      else if( SCIPisFeasNegative(masterprob, dual) /* || SCIPisInfinity(origprob, SCIPgetLhsLinear(origprob, origcons)) */)
      {
         infeasibility = GCGconsGetRhs(origprob, origcons) - activity;
      }

      assert(subgradientconsvals != NULL);
      assert(!SCIPisInfinity(masterprob, SQR(infeasibility)));

      subgradientconsvals[i] = infeasibility;

      if( SCIPisPositive(masterprob, SQR(infeasibility)) )
         subgradientnorm += SQR(infeasibility);

      SCIPfreeBufferArray(origprob, &vals);
      SCIPfreeBufferArray(origprob, &vars);
   }

   /* originalcuts */
   for( int i = 0; i < noriginalsepacuts; ++i )
   {
      SCIP_COL** cols;
      SCIP_Real* vals;
      int nvars;
      SCIP_Real activity;
      SCIP_Real infeasibility;

      SCIP_ROW* originalsepaorigcut = originalsepaorigcuts[i];
      nvars = SCIProwGetNNonz(originalsepaorigcut);
      cols = SCIProwGetCols(originalsepaorigcut);
      vals = SCIProwGetVals(originalsepaorigcut);

      activity = 0.0;

      SCIP_Real dual = stabcenteroriginalsepacutvals[i];
      assert(!SCIPisInfinity(masterprob, ABS(dual)));
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
            val = SCIPgetSolVal(masterprob, (SCIP_SOL*) NULL, mastervar);
            assert(!SCIPisInfinity(masterprob, ABS(val)));
         }
         /* Variable in a pricing problem */
         else
         {
            int block = GCGvarGetBlock(var);
            if( !GCGisPricingprobRelevant(gcg, block) )
               continue;

            assert(pricingcols[block] != NULL);

            SCIP_VAR* pricingvar = GCGoriginalVarGetPricingVar(var);
            assert(GCGvarIsPricing(pricingvar));
            val = GCGcolGetSolVal(pricingcols[block], pricingvar);
            assert(!SCIPisInfinity(masterprob, ABS(val)));
         }
         assert(stabcenteroriginalsepacutvals != NULL);
         assert(vals != NULL);
         activity += vals[j] * val;
      }

      infeasibility = 0.0;

      if( SCIPisFeasPositive(masterprob, dual) )
      {
         infeasibility = SCIProwGetLhs(originalsepaorigcut) - activity;
      }
      else if( SCIPisFeasNegative(masterprob, dual) )
      {
         infeasibility = SCIProwGetRhs(originalsepaorigcut) - activity;
      }

      assert(subgradientoriginalsepacutvals != NULL);
      assert(!SCIPisInfinity(masterprob, SQR(infeasibility)));

      subgradientoriginalsepacutvals[i] = infeasibility;

      if( SCIPisPositive(masterprob, SQR(infeasibility)) )
         subgradientnorm += SQR(infeasibility);
   }

   /* extended master conss */
   GCG_PRICINGMODIFICATION** pricingmods;
   int block;
   for( int i = 0; i < nsubgradientextendedmasterconss; ++i )
   {
      SCIP_COL** cols;
      SCIP_Real* vals;
      int nvars;
      SCIP_Real val;
      SCIP_VAR* var;
      SCIP_Real activity;
      SCIP_Real infeasibility;
      GCG_EXTENDEDMASTERCONSDATA* tmpextendedmasterconsdata;

      tmpextendedmasterconsdata = subgradientextendedmasterconss[i];
      assert(tmpextendedmasterconsdata != NULL);

      nvars = GCGextendedmasterconsGetNNonz(gcg, tmpextendedmasterconsdata);
      cols = GCGextendedmasterconsGetCols(gcg, tmpextendedmasterconsdata);
      vals = GCGextendedmasterconsGetVals(gcg, tmpextendedmasterconsdata);

      activity = 0.0;

      SCIP_Real dual = stabcenterextendedmasterconsvals[i];
      assert(!SCIPisInfinity(masterprob, ABS(dual)));
      for( int j = 0; j < nvars; ++j )
      {
         val = 0.0;
         var = SCIPcolGetVar(cols[j]);
         assert(GCGvarIsMaster(var));

         /* only linking or static master variable */
         if( GCGvarGetBlock(var) >= 0 )
            continue;

         val = SCIPgetSolVal(masterprob, (SCIP_SOL*) NULL, var);
         assert(!SCIPisInfinity(masterprob, ABS(val)));
         val = val * vals[j];

         assert(vals != NULL);
         activity += val;
      }

      pricingmods = GCGextendedmasterconsGetPricingModifications(tmpextendedmasterconsdata);
      for( int j = 0; j < GCGextendedmasterconsGetNPricingModifications(tmpextendedmasterconsdata); j++ )
      {
         block = GCGpricingmodificationGetBlock(pricingmods[j]);
         assert(pricingmods[j] != NULL);
         assert(GCGisPricingprobRelevant(gcg, block));
         assert(pricingcols[block] != NULL);

         SCIP_VAR* pricingvar = GCGpricingmodificationGetCoefVar(pricingmods[j]);
         assert(GCGvarIsInferredPricing(pricingvar));
         val = GCGcolGetSolVal(pricingcols[block], pricingvar);
         assert(!SCIPisInfinity(masterprob, ABS(val)));

         activity += val;
      }

      infeasibility = 0.0;

      if( SCIPisFeasPositive(masterprob, dual) )
      {
         infeasibility = GCGextendedmasterconsGetLhs(gcg, tmpextendedmasterconsdata) - activity;
      }
      else if( SCIPisFeasNegative(masterprob, dual) )
      {
         infeasibility = GCGextendedmasterconsGetRhs(gcg, tmpextendedmasterconsdata) - activity;
      }

      assert(!SCIPisInfinity(masterprob, SQR(infeasibility)));

      subgradientextendedmasterconsvals[i] = infeasibility;

      if( SCIPisPositive(masterprob, SQR(infeasibility)) )
         subgradientnorm += SQR(infeasibility);
   }

   /* linkingconss */
   for( int i = 0; i < nlinkingconss; ++i )
   {
      SCIP_VAR* mastervar;
      SCIP_VAR* pricingvar;
      SCIP_CONS* linkingcons = linkingconss[i];
      SCIP_Real activity;
      SCIP_Real infeasibility;
      mastervar = SCIPgetVarsLinear(masterprob, linkingcons)[0];
      assert(GCGvarIsMaster(mastervar));

      block = linkingconsblocks[i];
      pricingvar = GCGlinkingVarGetPricingVars(GCGmasterVarGetOrigvars(mastervar)[0])[block];
      assert(GCGvarIsPricing(pricingvar));
      assert(pricingcols[block] != NULL);

      assert(stabcenterlinkingconsvals != NULL);
      SCIP_Real masterval = SCIPgetSolVal(masterprob, (SCIP_SOL*) NULL, mastervar);
      SCIP_Real pricingval = GCGcolGetSolVal(pricingcols[block], pricingvar);
      assert(!SCIPisInfinity(masterprob, ABS(masterval)));
      assert(!SCIPisInfinity(masterprob, ABS(pricingval)));
      activity = (masterval - pricingval);

      infeasibility = activity;

      assert(subgradientlinkingconsvals != NULL);
      assert(!SCIPisInfinity(masterprob, SQR(infeasibility)));

      subgradientlinkingconsvals[i] = infeasibility;

      if( SCIPisPositive(masterprob, SQR(infeasibility)) )
         subgradientnorm += SQR(infeasibility);
   }

   assert(!SCIPisNegative(masterprob, subgradientnorm));

   subgradientnorm = sqrt(subgradientnorm);

   SCIPdebugMessage("Update subgradient and subgradientnorm with value %g.\n", subgradientnorm);
   return SCIP_OKAY;
}

/**< calculate norm of difference between stabcenter and current duals */
void Stabilization::calculateDualdiffnorm()
{
   SCIP_CONS** masterconss = GCGgetMasterConss(gcg);
   SCIP_CONS** linkingconss = GCGgetVarLinkingconss(gcg);
   int nlinkingconss = GCGgetNVarLinkingconss(gcg);
   assert(nstabcenterlinkingconsvals <= GCGgetNVarLinkingconss(gcg) );
   int nconss = GCGgetNMasterConss(gcg);
   assert(nconss <= nstabcenterconsvals);
   SCIP_ROW** originalsepamastercuts = GCGsepaGetOriginalSepaMastercuts(gcg);
   int noriginalsepacuts = GCGsepaGetNOriginalSepaCuts(gcg);
   assert(noriginalsepacuts <= nstabcenteroriginalsepacutvals);

   dualdiffnorm = 0.0;

   /* masterconss */
   assert(stabcenterconsvals != NULL);

   for( int i = 0; i < nconss; ++i )
   {
      SCIP_Real dualdiff = SQR(stabcenterconsvals[i] - pricingtype->consGetDual(masterconss[i]));

      if( SCIPisPositive(masterprob, dualdiff) )
         dualdiffnorm += dualdiff;
   }

   /* originalcuts */
   assert(stabcenterconsvals != NULL);

   for( int i = 0; i < noriginalsepacuts; ++i )
   {
      SCIP_Real dualdiff = SQR(stabcenteroriginalsepacutvals[i] - pricingtype->rowGetDual(originalsepamastercuts[i]));

      if( SCIPisPositive(masterprob, dualdiff) )
         dualdiffnorm += dualdiff;
   }

   /* extended master conss */
   for( int i = 0; i < nstabcenterextendedmasterconss; ++i )
   {
      GCG_EXTENDEDMASTERCONSDATA* tmpextendedmasterconsdata;

      tmpextendedmasterconsdata = stabcenterextendedmasterconss[i];
      assert(tmpextendedmasterconsdata != NULL);

      SCIP_Real dualdiff = SQR(stabcenterextendedmasterconsvals[i] - pricingtype->extendedmasterconsGetDual(tmpextendedmasterconsdata));

      if( SCIPisPositive(masterprob, dualdiff) )
         dualdiffnorm += dualdiff;
   }

   /* linkingconss */
   assert(stabcenterlinkingconsvals != NULL);

   for( int i = 0; i < nlinkingconss; ++i )
   {
      SCIP_Real dualdiff = SQR(stabcenterlinkingconsvals[i] - pricingtype->consGetDual(linkingconss[i]));

      if( SCIPisPositive(masterprob, dualdiff) )
         dualdiffnorm += dualdiff;
   }

   dualdiffnorm = sqrt(dualdiffnorm);
   SCIPdebugMessage("Update dualdiffnorm with value %g.\n", dualdiffnorm);
}

/**< calculate beta */
void Stabilization::calculateBeta()
{
   SCIP_CONS** masterconss = GCGgetMasterConss(gcg);
   SCIP_CONS** linkingconss = GCGgetVarLinkingconss(gcg);
   int nlinkingconss = GCGgetNVarLinkingconss(gcg);
   assert(nstabcenterlinkingconsvals <= GCGgetNVarLinkingconss(gcg) );
   int nconss = GCGgetNMasterConss(gcg);
   assert(nconss <= nstabcenterconsvals);
   SCIP_ROW** originalsepamastercuts = GCGsepaGetOriginalSepaMastercuts(gcg);
   int noriginalsepacuts = GCGsepaGetNOriginalSepaCuts(gcg);
   assert(noriginalsepacuts <= nstabcenteroriginalsepacutvals);

   beta = 0.0;

   /* masterconss */
   assert(stabcenterconsvals != NULL);

   for( int i = 0; i < nconss; ++i )
   {
      SCIP_Real dualdiff = ABS(pricingtype->consGetDual(masterconss[i]) - stabcenterconsvals[i]);
      SCIP_Real product = dualdiff * ABS(subgradientconsvals[i]);

      if( SCIPisPositive(masterprob, product) )
         beta += product;
   }

   /* originalcuts */
   assert(stabcenteroriginalsepacutvals != NULL || noriginalsepacuts == 0);

   for( int i = 0; i < noriginalsepacuts; ++i )
   {
      SCIP_Real dualdiff = ABS(pricingtype->rowGetDual(originalsepamastercuts[i]) - stabcenteroriginalsepacutvals[i]);
      SCIP_Real product = dualdiff * ABS(subgradientoriginalsepacutvals[i]);

      if( SCIPisPositive(masterprob, product) )
         beta += product;
   }

   /* extended master conss */
   for( int i = 0; i < nstabcenterextendedmasterconss; ++i )
   {
      GCG_EXTENDEDMASTERCONSDATA* tmpextendedmasterconsdata;

      tmpextendedmasterconsdata = stabcenterextendedmasterconss[i];
      assert(tmpextendedmasterconsdata != NULL);

      SCIP_Real dualdiff = ABS(pricingtype->extendedmasterconsGetDual(tmpextendedmasterconsdata) - stabcenterextendedmasterconsvals[i]);
      SCIP_Real product = dualdiff * ABS(subgradientextendedmasterconsvals[i]);

      if( SCIPisPositive(masterprob, product) )
         beta += product;
   }

   /* linkingconss */
   assert(stabcenterlinkingconsvals != NULL);

   for( int i = 0; i < nlinkingconss; ++i )
   {
      SCIP_Real dualdiff = ABS(pricingtype->consGetDual(linkingconss[i]) - stabcenterlinkingconsvals[i]);
      SCIP_Real product = dualdiff * ABS(subgradientlinkingconsvals[i]);

      if( SCIPisPositive(masterprob, product) )
         beta += product;
   }

   if( SCIPisPositive(masterprob, subgradientnorm) )
      beta = beta / (subgradientnorm * dualdiffnorm);

   SCIPdebugMessage("Update beta with value %g.\n", beta);

   assert( ( SCIPisPositive(masterprob, beta) || SCIPisZero(masterprob, subgradientnorm)) && SCIPisLE(masterprob, beta, 1.0) );
}

/**< calculate factor that is needed in hybrid stabilization */
void Stabilization::calculateHybridFactor()
{
   SCIP_CONS** masterconss = GCGgetMasterConss(gcg);

   SCIP_CONS** linkingconss = GCGgetVarLinkingconss(gcg);
   int nlinkingconss = GCGgetNVarLinkingconss(gcg);
   assert(nstabcenterlinkingconsvals <= GCGgetNVarLinkingconss(gcg) );
   int nconss = GCGgetNMasterConss(gcg);
   assert(nconss <= nstabcenterconsvals);
   SCIP_ROW** originalsepamastercuts = GCGsepaGetOriginalSepaMastercuts(gcg);
   int noriginalsepacuts = GCGsepaGetNOriginalSepaCuts(gcg);
   assert(noriginalsepacuts <= nstabcenteroriginalsepacutvals);

   SCIP_Real divisornorm = 0.0;

   /* masterconss */
   assert(stabcenterconsvals != NULL);

   for( int i = 0; i < nconss; ++i )
   {
      SCIP_Real divisor = SQR((beta - 1.0) * stabcenterconsvals[i]
                        + beta * (subgradientconsvals[i] * dualdiffnorm / subgradientnorm)
                        + (1 - beta) * pricingtype->consGetDual(masterconss[i]));

      if( SCIPisPositive(masterprob, divisor) )
         divisornorm += divisor;
   }

   /* originalcuts */
   assert(stabcenteroriginalsepacutvals != NULL);

   for( int i = 0; i < noriginalsepacuts; ++i )
   {
      SCIP_Real divisor = SQR((beta - 1.0) * stabcenteroriginalsepacutvals[i]
                        + beta * (subgradientoriginalsepacutvals[i] * dualdiffnorm / subgradientnorm)
                        + (1 - beta) * pricingtype->rowGetDual(originalsepamastercuts[i]));

      if( SCIPisPositive(masterprob, divisor) )
         divisornorm += divisor;
   }

   /* extended master conss */
   for( int i = 0; i < nstabcenterextendedmasterconss; ++i )
   {
      GCG_EXTENDEDMASTERCONSDATA* tmpextendedmasterconsdata;

      tmpextendedmasterconsdata = stabcenterextendedmasterconss[i];
      assert(tmpextendedmasterconsdata != NULL);

      SCIP_Real divisor = SQR((beta - 1.0) * stabcenterextendedmasterconsvals[i]
                        + beta * (subgradientextendedmasterconsvals[i] * dualdiffnorm / subgradientnorm)
                        + (1 - beta) * pricingtype->extendedmasterconsGetDual(tmpextendedmasterconsdata));

      if( SCIPisPositive(masterprob, divisor) )
         divisornorm += divisor;
   }

   /* linkingconss */
   assert(stabcenterlinkingconsvals != NULL);

   for( int i = 0; i < nlinkingconss; ++i )
   {
      SCIP_Real divisor = SQR((beta - 1.0) * stabcenterlinkingconsvals[i]
                        + beta * (subgradientlinkingconsvals[i] * dualdiffnorm / subgradientnorm)
                        + (1 - beta) * pricingtype->consGetDual(linkingconss[i]));

      if( SCIPisPositive(masterprob, divisor) )
         divisornorm += divisor;
   }

   divisornorm = sqrt(divisornorm);

   hybridfactor = ((1 - alpha) * dualdiffnorm) / divisornorm;

   SCIPdebugMessage("Update hybridfactor with value %g.\n", hybridfactor);

   assert( SCIPisPositive(masterprob, hybridfactor) );
}


SCIP_Bool Stabilization::isStabilized()
{
   if( inmispricingschedule )
      return SCIPisGT(masterprob, alphabar, 0.0);
   return SCIPisGT(masterprob, alpha, 0.0);
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

/** update subgradient product */
SCIP_RETCODE Stabilization::updateSubgradientProduct(
   GCG_COL**            pricingcols         /**< solutions of the pricing problems */
)
{
   /* first update the arrays */
   SCIP_CALL( updateStabcenterconsvals() );
   SCIP_CALL( updateStabcenteroriginalcutvals() );
   SCIP_CALL( updateStabcenterextendedmasterconsvals() );

   subgradientproduct = calculateSubgradientProduct(pricingcols);

   return SCIP_OKAY;
}


} /* namespace gcg */
