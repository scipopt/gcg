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

/**@file   class_pricingcontroller.h
 * @brief  pricing controller managing the pricing strategy
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef CLASS_PRICINGCONTROLLER_H_
#define CLASS_PRICINGCONTROLLER_H_

#include "class_pricingtype.h"
#include "type_pricingjob.h"
#include "objscip/objscip.h"

namespace gcg {

class Pricingcontroller
{

private:
   SCIP*                 scip_;              /**< SCIP instance (master problem) */
   GCG_PRICINGJOB**      pricingjobs;        /**< pricing jobs, one per pricing problem */
   int                   npricingprobs;      /**< number of pricing problems */

   /* parameters */
   int                   sorting;            /**< how should pricing problems be sorted */
   SCIP_Real             successfulmipsrel;  /**< factor of MIPs to be solved successfully until pricing is aborted */
   int                   eagerfreq;          /**< frequency at which all pricingproblems should be solved */

   /* strategic variables */
   SCIP_Real*            score;              /**< scores of the pricing problems */
   int*                  order;              /**< current order of the pricing problems */

   /* statistics */
   int                   previdx;            /**< previously solved pricing problem */
   int                   eagerage;           /**< iterations since last eager iteration */



public:
   /** constructor */
   Pricingcontroller(SCIP* scip);

   /** destructor */
   virtual ~Pricingcontroller();

   SCIP_RETCODE addParameters();

   SCIP_RETCODE initSol();

   SCIP_RETCODE exitSol();

   /** sorts pricing problems according to their score */
   void sortPricingProblems(
      SCIP_Real*            dualsolconv,         /**< array of dual solutions for the convexity constraints */
      int*                  npointsprob,
      int*                  nraysprob

   );

   /** get the next pricing problem to be solved */
   int getNextPricingprob();

   /** returns whether pricing can be aborted */
   SCIP_Bool abortPricing(
      PricingType*          pricetype,          /**< type of pricing */
      int                   nfoundvars,         /**< number of variables found so far */
      int                   solvedmips,         /**< number of MIPS solved so far */
      int                   successfulmips,     /**< number of successful mips solved so far */
      SCIP_Bool             optimal             /**< optimal or heuristic pricing */
   ) const;

   void resetEagerage();

   void increaseEagerage();
};

} /* namespace gcg */
#endif /* CLASS_PRICINGCONTROLLER_H_ */
