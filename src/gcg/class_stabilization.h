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

/**@file   class_stabilization.h
 * @brief  class with functions for dual variable smoothing
 * @author Martin Bergner
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef CLASS_STABILIZATION_H_
#define CLASS_STABILIZATION_H_

#include "gcg/class_pricingtype.h"
#include "gcg/type_extendedmasterconsdata.h"
#include <scip/type_misc.h>
#include <scip/type_retcode.h>

namespace gcg {

/**
 * @ingroup PRICING_PRIV
 * @{
 */

class Stabilization
{
private:
   GCG* gcg;
   SCIP* masterprob;
   SCIP_Real* stabcenterconsvals;
   int stabcenterconsvalssize;
   int nstabcenterconsvals;
   SCIP_Real* stabcenteroriginalsepacutvals;
   int stabcenteroriginalsepacutvalssize;
   int nstabcenteroriginalsepacutvals;
   GCG_EXTENDEDMASTERCONSDATA** stabcenterextendedmasterconss;
   int nstabcenterextendedmasterconss;
   int stabcenterextendedmasterconsssize;
   SCIP_Real* stabcenterextendedmasterconsvals;
   SCIP_Real* stabcenterlinkingconsvals;
   int nstabcenterlinkingconsvals;
   int stabcenterlinkingconsvalssize;
   SCIP_Real* stabcenterconv;
   int nstabcenterconv;
   SCIP_Real dualdiffnorm; /**< norm of difference between stabcenter and current duals */
   SCIP_Real* subgradientconsvals;
   int subgradientconsvalssize;
   int nsubgradientconsvals;
   SCIP_Real* subgradientoriginalsepacutvals;
   int subgradientoriginalsepacutvalssize;
   int nsubgradientoriginalsepacutvals;
   GCG_EXTENDEDMASTERCONSDATA** subgradientextendedmasterconss;
   int nsubgradientextendedmasterconss;
   int subgradientextendedmasterconsssize;
   SCIP_Real* subgradientextendedmasterconsvals;
   SCIP_Real* subgradientlinkingconsvals;
   int subgradientlinkingconsvalssize;
   SCIP_Real subgradientnorm;
   SCIP_Real hybridfactor;
   PricingType* pricingtype;
   SCIP_Real alpha;
   SCIP_Real alphabar; /**< alpha that is used and updated in a mispricing schedule */
   SCIP_Bool hybridascent; /**< hybridize smoothing with an ascent method? */
   SCIP_Real beta;
   SCIP_Longint nodenr;
   int k; /**< counter for the number of stabilized pricing rounds in B&B node, excluding the mispricing schedule iterations  */
   int t; /**< counter for the number of pricing rounds during a mispricing schedule, restarted after a mispricing schedule is finished */
   SCIP_Bool hasstabilitycenter;
   SCIP_Real stabcenterbound;
   SCIP_Bool inmispricingschedule; /**< currently in mispricing schedule */
   SCIP_Real subgradientproduct;

public:
   /** constructor */
   Stabilization(
      GCG*               gcgstruct,          /**< SCIP data structure */
      PricingType*       pricingtype,        /**< the pricing type when the stabilization should run */
      SCIP_Bool          hybridascent        /**< enable hybridization of smoothing with an ascent method? */
   );
   /** constructor */
   Stabilization();

   /** destructor */
   virtual ~Stabilization();

   /** gets the stabilized dual solution of constraint at position i */
   SCIP_RETCODE consGetDual(
      int                i,                  /**< index of the constraint */
      SCIP_Real*         dual                /**< return pointer for dual value */
   );

   /** gets the stabilized dual solution of cut at position i */
   SCIP_RETCODE rowGetDual(
      int                i,                  /**< index of the row */
      SCIP_Real*         dual                /**< return pointer for dual value */
   );

   /** gets the stabilized dual of the convexity constraint at position i */
   SCIP_Real convGetDual(
      int i
   );

   SCIP_RETCODE extendedmasterconsGetDual(
      GCG_EXTENDEDMASTERCONSDATA*    extendedmasterconsdata,      /**< extendedmasterconsdata */
      SCIP_Real*            dual                /**< return pointer for dual value */
   );

   /** updates the stability center if the bound has increased */
   SCIP_RETCODE updateStabilityCenter(
      SCIP_Real             lowerbound,         /**< lower bound due to lagrange function corresponding to current (stabilized) dual vars */
      SCIP_Real*            dualsolconv,        /**< corresponding feasible dual solution for convexity constraints */
      GCG_COL**             pricingcols         /**< columns of the pricing problems */
   );

   /** updates the alpha after unsuccessful pricing */
   void updateAlphaMisprice();

   /** updates the alpha after successful pricing */
   void updateAlpha();

   /** returns whether the stabilization is active */
   SCIP_Bool isStabilized();

   /** enabling mispricing schedule */
   void activateMispricingSchedule(
   );

   /** disabling mispricing schedule */
   void disablingMispricingSchedule(
   );

   /** is mispricing schedule enabled */
   SCIP_Bool isInMispricingSchedule(
   ) const;

   /** sets the variable linking constraints in the master */
   SCIP_RETCODE setLinkingConss(
      SCIP_CONS**        linkingconss,       /**< array of linking master constraints */
      int*               linkingconsblock,   /**< block of the linking constraints */
      int                nlinkingconss       /**< size of the array */
   );

   /** increases the number of new variable linking constraints */
   SCIP_RETCODE setNLinkingconsvals(
      int                nlinkingconssnew    /**< number of new linking constraints */
   );

   /** increases the number of new convexity constraints */
   SCIP_RETCODE setNConvconsvals(
      int nconvconssnew
   );

   /** gets the dual of variable linking constraints at index i */
   SCIP_Real linkingconsGetDual(
      int i
      );

   /**< update node */
   void updateNode();

   /**< update information for hybrid stablization with dual ascent */
   SCIP_RETCODE updateHybrid();

   /** update subgradient product */
   SCIP_RETCODE updateSubgradientProduct(
      GCG_COL**            pricingcols         /**< solutions of the pricing problems */
   );

private:
   /** updates the number of iterations */
   void updateIterationCount();

   /** updates the number of iterations in the current mispricing schedule */
   void updateIterationCountMispricing();

   /** updates the constraints in the stability center (and allocates more memory) */
   SCIP_RETCODE updateStabcenterconsvals();

   /** updates the original cuts in the stability center (and allocates more memory) */
   SCIP_RETCODE updateStabcenteroriginalcutvals();

   /** updates the extended master conss in the stability center (and allocates more memory) */
   SCIP_RETCODE updateStabcenterextendedmasterconsvals();

   /** updates the constraints in the subgradient (and allocates more memory) */
   SCIP_RETCODE updateSubgradientconsvals();

   /** updates the original cuts in the subgradient (and allocates more memory) */
   SCIP_RETCODE updateSubgradientoriginalcutvals();

   /** updates the extended master conss in the subgradient (and allocates more memory) */
   SCIP_RETCODE updateSubgradientextendedmasterconsvals();

   /** increase the alpha value */
   void increaseAlpha();

   /** decrease the alpha value */
   void decreaseAlpha();

   /** calculates the product of subgradient (with linking variables)
    * with the difference of current duals and the stability center */
   SCIP_Real calculateSubgradientProduct(
      GCG_COL**            pricingcols         /**< columns of the pricing problems */
   );

   /** calculates the normalized subgradient (with linking variables) multiplied
    * with the norm of the difference of current duals and the stability center */
   SCIP_RETCODE calculateSubgradient(
      GCG_COL**            pricingcols         /**< columns of the pricing problems */
   );

   /**< calculate norm of difference between stabcenter and current duals */
   void calculateDualdiffnorm();

   /**< calculate beta */
   void calculateBeta();

   /**< calculate factor that is needed in hybrid stabilization */
   void calculateHybridFactor();

   /** computes the new dual value based on the current and the stability center values */
   SCIP_Real computeDual(
      SCIP_Real          center,             /**< value of stabilility center */
      SCIP_Real          current,            /**< current dual value */
      SCIP_Real          subgradient,         /**< subgradient (or 0.0 if not needed) */
      SCIP_Real          lhs,                 /**< lhs (or 0.0 if not needed) */
      SCIP_Real          rhs                  /**< rhs (or 0.0 if not needed) */
   ) const;
};

/** @} */
} /* namespace gcg */
#endif /* CLASS_STABILIZATION_H_ */
