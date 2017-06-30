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

/**@file   class_pricingcontroller.cpp
 * @brief  pricing controller managing the pricing strategy
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_pricingcontroller.h"
#include "pricer_gcg.h"
#include "class_pricingtype.h"
#include "gcg.h"
#include "pub_pricingjob.h"

#include "scip/scip.h"
#include "objscip/objscip.h"

#include <exception>

#define DEFAULT_SORTING                  2          /**< default sorting method for pricing mips
                                                     *    0 :   order of pricing problems
                                                     *    1 :   according to dual solution of convexity constraint
                                                     *    2 :   according to reliability from previous round)
                                                     */
#define DEFAULT_SUCCESSFULMIPSREL        1.0        /**< factor of successful mips to be solved */
#define DEFAULT_EAGERFREQ                10         /**< frequency at which all pricingproblems should be solved (0 to disable) */

#define SCIP_CALL_EXC(x)   do                                                                                 \
                       {                                                                                      \
                          SCIP_RETCODE _retcode_;                                                             \
                          if( (_retcode_ = (x)) !=  SCIP_OKAY )                                               \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _retcode_);                    \
                             throw std::exception();                                                          \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

namespace gcg {

Pricingcontroller::Pricingcontroller(
   SCIP*                  scip
   )
{
   scip_ = scip;
   pricingjobs = NULL;
   npricingprobs = 0;

   sorting = DEFAULT_SORTING;
   successfulmipsrel = DEFAULT_SUCCESSFULMIPSREL;
   eagerfreq = DEFAULT_EAGERFREQ;

   score = NULL;
   order = NULL;

   previdx = -1;
   eagerage = 0;
}

Pricingcontroller::~Pricingcontroller()
{
}

SCIP_RETCODE Pricingcontroller::addParameters()
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/sorting",
         "which sorting method should be used to sort the pricing problems (0 = order of pricing problems, 1 = according to dual solution of convexity constraint, 2 = according to reliability from previous round)",
         &sorting, FALSE, DEFAULT_SORTING, 0, 5, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(origprob, "pricing/masterpricer/successfulsubmipsrel",
         "part of the submips that are solved and lead to new variables before pricing round is aborted? (1.0 = solve all pricing MIPs)",
         &successfulmipsrel, FALSE, DEFAULT_SUCCESSFULMIPSREL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/eagerfreq",
         "frequency at which all pricingproblems should be solved (0 to disable)",
         &eagerfreq, FALSE, DEFAULT_EAGERFREQ, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

SCIP_RETCODE Pricingcontroller::initSol()
{
   npricingprobs = GCGgetNPricingprobs(GCGmasterGetOrigprob(scip_));
   eagerage = 0;

   SCIP_CALL_EXC( SCIPallocBlockMemoryArray(scip_, &pricingjobs, npricingprobs) );
   for( int i = 0; i < npricingprobs; ++i )
   {
      SCIP_CALL_EXC( GCGcreatePricingjob(scip_, &pricingjobs[i], GCGgetPricingprob(scip_, i), i) );
   }

   SCIP_CALL_EXC( SCIPallocBlockMemoryArray(scip_, &score, npricingprobs) );
   SCIP_CALL_EXC( SCIPallocBlockMemoryArray(scip_, &order, npricingprobs) );

   return SCIP_OKAY;
}

SCIP_RETCODE Pricingcontroller::exitSol()
{
   SCIPfreeBlockMemoryArray(scip_, &score, npricingprobs);
   SCIPfreeBlockMemoryArray(scip_, &order, npricingprobs);

   for( int i = 0; i < npricingprobs; ++i )
   {
      GCGfreePricingjob(scip_, &pricingjobs[i]);
   }
   SCIPfreeBlockMemoryArray(scip_, &pricingjobs, npricingprobs);

   return SCIP_OKAY;
}

/** sorts pricing problems according to their score */
void Pricingcontroller::sortPricingProblems(
   SCIP_Real*            dualsolconv,         /**< array of dual solutions for the convexity constraints */
   int*                  npointsprob,
   int*                  nraysprob
)
{
   /** @todo sort w.r.t. other measures? Don't sort in Farkas pricing? Randomized? */
   for( int i = 0; i < npricingprobs; i++ )
   {
      order[i] = i;
      switch( sorting )
      {
      case 1:
         score[i] = dualsolconv[i];
         break;
      case 2:
         score[i] = -(0.2 * npointsprob[i] + nraysprob[i]);
         break;
      default:
         score[i] = 0.0;
         break;
      }
   }

   if( sorting > 0 )
      SCIPsortDownRealInt(score, order, npricingprobs);

   previdx = -1;
}

/** get the next pricing problem to be solved */
int Pricingcontroller::getNextPricingprob()
{
   ++previdx;

   if( previdx == npricingprobs )
      return -1;
   else
      return order[previdx];
}

/** returns whether pricing can be aborted */
SCIP_Bool Pricingcontroller::abortPricing(
   PricingType*          pricetype,          /**< type of pricing */
   int                   nfoundvars,         /**< number of variables found so far */
   int                   solvedmips,         /**< number of MIPS solved so far */
   int                   successfulmips,     /**< number of successful mips solved so far */
   SCIP_Bool             optimal             /**< optimal or heuristic pricing */
) const
{
   int nrelpricingprobs = GCGgetNRelPricingprobs(GCGmasterGetOrigprob(scip_));

   if( eagerage == eagerfreq )
      return FALSE;

   if( optimal )
      return pricetype->canOptimalPricingBeAborted(nfoundvars, solvedmips, successfulmips, successfulmipsrel, nrelpricingprobs);
   else
      return pricetype->canHeuristicPricingBeAborted(nfoundvars, solvedmips, successfulmips, successfulmipsrel, nrelpricingprobs);
}

void Pricingcontroller::resetEagerage()
{
   eagerage = 0;
}

void Pricingcontroller::increaseEagerage()
{
   if( eagerfreq > 0 )
      eagerage++;
}

} /* namespace gcg */
