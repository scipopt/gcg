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

/**@file   class_stabilization.h
 * @brief  class with functions for dual variable smoothing
 * @author Martin Bergner
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef CLASS_STABILIZATION_H_
#define CLASS_STABILIZATION_H_

#include "objscip/objscip.h"
#include "class_pricingtype.h"

namespace gcg {

class Stabilization
{
private:
   SCIP* scip_;
   SCIP_Real* stabcenterconss;
   int stabcenterconsssize;
   int nstabcenterconss;
   SCIP_Real* stabcentercuts;
   int stabcentercutssize;
   int nstabcentercuts;
   SCIP_Real* stabcenterlinkingconss;
   int nstabcenterlinkingconss;
   SCIP_Real* stabcenterconv;
   int nstabcenterconv;
   PricingType* pricingtype;
   SCIP_Real alpha;
   SCIP_Real alphabar; /**< alpha that is used and updated in a mispricing schedule */
   SCIP_Longint nodenr;
   int k; /**< counter for the number of stabilized pricing rounds in B&B node, excluding the mispricing schedule iterations  */
   int t; /**< counter for the number of pricing rounds during a mispricing schedule, restarted after a mispricing schedule is finished */
   SCIP_Bool hasstabilitycenter;
   SCIP_Real stabcenterbound;
   SCIP_Bool inmispricingschedule; /**< currently in mispricing schedule */

public:
   /** constructor */
   Stabilization(
      SCIP*              scip,               /**< SCIP data structure */
      PricingType*       pricingtype         /**< the pricing type when the stabilization should run */
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

   /** updates the stability center if the bound has increased */
   SCIP_RETCODE updateStabilityCenter(
      SCIP_Real             lowerbound,         /**< lower bound due to lagrange function corresponding to current (stabilized) dual vars */
      SCIP_Real*            dualsolconv         /**< corresponding feasible dual solution for convexity constraints */
   );

   /** updates the alpha after unsuccessful pricing */
   void updateAlphaMisprice();

   /** updates the alpha after successful pricing */
   void updateAlpha(
      GCG_COL**            pricingcols         /**< columns of the pricing problems */
   );

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
   SCIP_RETCODE setNLinkingconss(
      int                nlinkingconssnew    /**< number of new linking constraints */
   );

   /** increases the number of new convexity constraints */
   SCIP_RETCODE setNConvconss(
      int nconvconssnew
   );

   /** gets the dual of variable linking constraints at index i */
   SCIP_Real linkingconsGetDual(
      int i
      );

   void updateNode();

private:
   /** updates the number of iterations */
   void updateIterationCount();

   /** updates the number of iterations in the current mispricing schedule */
   void updateIterationCountMispricing();

   /** updates the constraints in the stability center (and allocates more memory) */
   SCIP_RETCODE updateStabcenterconss();

   /** updates the cuts in the stability center (and allocates more memory) */
   SCIP_RETCODE updateStabcentercuts();

   /** increase the alpha value */
   void increaseAlpha();

   /** decrease the alpha value */
   void decreaseAlpha();

   /** calculates the subgradient (with linking variables */
   SCIP_Real calculateSubgradient(
      GCG_COL**            pricingcols         /**< columns of the pricing problems */
   );

   /** computes the new dual value based on the current and the stability center values */
   SCIP_Real computeDual(
      SCIP_Real          center,             /**< value of stabilility center */
      SCIP_Real          current             /**< current dual value */
   ) const;
};

} /* namespace gcg */
#endif /* CLASS_STABILIZATION_H_ */
