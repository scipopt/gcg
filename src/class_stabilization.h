/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2014 Operations Research, RWTH Aachen University       */
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
   SCIP_NODE* node;
   int k;
   SCIP_Bool hasstabilitycenter;

public:
   Stabilization(
      SCIP* scip,
      PricingType* pricingtype
      );
   virtual ~Stabilization();

   SCIP_Real consGetDual(
      int i
      );

   SCIP_Real rowGetDual(
      int i
      );

   SCIP_Real convGetDual(
      int i
      );

   SCIP_RETCODE updateStabilityCenter(
      SCIP_Real lowerbound
      );
   void updateAlphaMisprice();
   void updateAlpha(
      SCIP_SOL**            pricingsols         /**< solutions of the pricing problems */
      );
   SCIP_Bool isStabilized();

   SCIP_RETCODE setLinkingConss(
      SCIP_CONS** linkingconss,
      int* linkingconsblock,
      int nlinkingconss
      );

   SCIP_RETCODE setNLinkingconss(
      int nlinkingconssnew
      );

   SCIP_RETCODE setNConvconss(
         int nconvconssnew
         );

   SCIP_Real linkingconsGetDual(
      int i
      );

private:
   void updateIterationCount();
   SCIP_RETCODE updateStabcenterconss();
   SCIP_RETCODE updateStabcentercuts();
   void increaseAlpha();
   void decreaseAlpha();
   SCIP_Real calculateSubgradient(
      SCIP_SOL**            pricingsols         /**< solutions of the pricing problems */
      );

   SCIP_Real computeDual(
      SCIP_Real center,
      SCIP_Real current
      );
};

} /* namespace gcg */
#endif /* CLASS_STABILIZATION_H_ */
