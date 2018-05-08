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

/**@file   class_pricingtype.cpp
 * @brief  abstraction for SCIP pricing types
 * @author Martin Bergner
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "class_pricingtype.h"
#include "scip/cons_linear.h"
#include "pub_gcgvar.h"
#include "scip/pub_lp.h"
#include "scip/clock.h"
#include "scip_misc.h"

#include <exception>

#define DEFAULT_MAXROUNDSREDCOST         INT_MAX    /**< maximal number of reduced cost pricing rounds */
#define DEFAULT_MAXCOLSROUNDREDCOSTROOT  100        /**< maximal number of columns per reduced cost pricing round at root node */
#define DEFAULT_MAXCOLSROUNDREDCOST      100        /**< maximal number of columns per reduced cost pricing round */
#define DEFAULT_MAXCOLSPROBREDCOSTROOT    10        /**< maximal number of columns per problem to be generated during red. cost pricing at root node */
#define DEFAULT_MAXCOLSPROBREDCOST        10        /**< maximal number of columns per problem to be generated during red. cost pricing */
#define DEFAULT_MAXSUCCESSFULPROBSREDCOST INT_MAX   /**< maximal number of successfully solved red. cost pricing problems */
#define DEFAULT_RELMAXPROBSREDCOSTROOT   1.0        /**< maximal percentage of red. cost pricing problems that are solved at root node if variables have already been found */
#define DEFAULT_RELMAXPROBSREDCOST       1.0        /**< maximal percentage of red. cost pricing problems that are solved if variables have already been found */

#define DEFAULT_MAXCOLSROUNDFARKAS        10        /**< maximal number of columns per Farkas pricing round */
#define DEFAULT_MAXCOLSPROBFARKAS         10        /**< maximal number of columns per problem to be generated during Farkas pricing */
#define DEFAULT_RELMAXPROBSFARKAS        1.0        /**< maximal percentage of Farkas pricing problems that are solved if variables have already been found */


#define SCIP_CALL_EXC(x)   do                                                                                 \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) !=  SCIP_OKAY )                                                \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             throw std::exception();                                                          \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

PricingType::PricingType(
   SCIP*                 scip
   )
{
   scip_ = scip;
   type  = GCG_PRICETYPE_UNKNOWN;

   calls = 0;
   maxrounds = INT_MAX;
   maxcolsroundroot = INT_MAX;
   maxcolsround = INT_MAX;
   maxcolsprobroot = INT_MAX;
   maxcolsprob = INT_MAX;
   maxsuccessfulprobs = INT_MAX;
   relmaxprobsroot = 1.0;
   relmaxprobs = 1.0;

   SCIP_CALL_EXC( SCIPcreateCPUClock(scip, &(clock)) );
}

PricingType::~PricingType()
{
   SCIP_CALL_ABORT( SCIPfreeClock(scip_, &(clock)) );

   scip_ = (SCIP*) NULL;
}

SCIP_RETCODE PricingType::startClock()
{
   SCIP_CALL( SCIPstartClock(scip_, clock) );
   return SCIP_OKAY;
}

SCIP_RETCODE PricingType::stopClock()
{
   SCIP_CALL( SCIPstopClock(scip_, clock) );
   return SCIP_OKAY;
}

SCIP_Real PricingType::getClockTime() const
{
   return SCIPgetClockTime(scip_, clock);
}

FarkasPricing::FarkasPricing(
      SCIP* scip
   ) : PricingType(scip)
{
   type = GCG_PRICETYPE_FARKAS;
}

SCIP_Real FarkasPricing::consGetDual(
      SCIP* scip,
      SCIP_CONS* cons
   ) const
{
   return SCIPgetDualfarkasLinear(scip, cons);
}

SCIP_Real FarkasPricing::rowGetDual(
      SCIP_Row* row
   ) const
{
   return SCIProwGetDualfarkas(row);
}

SCIP_Real FarkasPricing::varGetObj(
      SCIP_VAR* var
   ) const
{
   assert(var != NULL);
   return 0.0;
}

SCIP_RETCODE FarkasPricing::addParameters()
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/maxcolsroundfarkas",
         "maximal number of columns per Farkas pricing round",
         &maxcolsround, FALSE, DEFAULT_MAXCOLSROUNDFARKAS, 1, INT_MAX, NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/maxcolsprobfarkas",
         "maximal number of columns per problem to be generated during Farkas pricing",
         &maxcolsprob, FALSE, DEFAULT_MAXCOLSPROBFARKAS, 1, INT_MAX, NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricing/masterpricer/relmaxprobsfarkas",
         "maximal percentage of Farkas pricing problems that are solved if variables have already been found",
         &relmaxprobs, FALSE, DEFAULT_RELMAXPROBSFARKAS, 0.0, 1.0, NULL, (SCIP_PARAMDATA*) NULL) );

   return SCIP_OKAY;
}

SCIP_Real ReducedCostPricing::consGetDual(
      SCIP* scip,
      SCIP_CONS* cons
   ) const
{
   return SCIPgetDualsolLinear(scip, cons);
}
SCIP_Real ReducedCostPricing::rowGetDual(
      SCIP_ROW* row
   ) const
{
   return SCIProwGetDualsol(row);
}

ReducedCostPricing::ReducedCostPricing(
      SCIP* p_scip
   ) : PricingType(p_scip)
{
   type = GCG_PRICETYPE_REDCOST;
}

SCIP_Real ReducedCostPricing::varGetObj(
      SCIP_VAR* var
   ) const
{
   SCIP_VAR* origvar;
   assert(var != NULL);

   origvar = GCGpricingVarGetOrigvars(var)[0];

   if( GCGoriginalVarIsLinking(origvar) )
      return 0.0;
   else
      return SCIPvarGetObj(origvar);
}

SCIP_RETCODE ReducedCostPricing::addParameters()
{
   SCIP* origprob = GCGmasterGetOrigprob(scip_);

   SCIP_CALL( SCIPaddIntParam(GCGmasterGetOrigprob(scip_), "pricing/masterpricer/maxroundsredcost",
         "maximal number of pricing rounds per node after the root node",
         &maxrounds, FALSE, DEFAULT_MAXROUNDSREDCOST, 0, INT_MAX, NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/maxcolsroundredcostroot",
         "maximal number of columns per reduced cost pricing round at root node",
         &maxcolsroundroot, FALSE, DEFAULT_MAXCOLSROUNDREDCOSTROOT, 0, INT_MAX,
         NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/maxcolsroundredcost",
         "maximal number of columns per reduced cost pricing round",
         &maxcolsround, FALSE, DEFAULT_MAXCOLSROUNDREDCOST, 0, INT_MAX,
         NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/maxcolsprobredcostroot",
         "maximal number of columns per problem to be generated during red. cost pricing at root node",
         &maxcolsprobroot, FALSE, DEFAULT_MAXCOLSPROBREDCOSTROOT, 0, INT_MAX,
         NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/maxcolsprobredcost",
         "maximal number of columns per problem to be generated during red. cost pricing",
         &maxcolsprob, FALSE, DEFAULT_MAXCOLSPROBREDCOST, 0, INT_MAX,
         NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/maxsuccessfulprobsredcost",
         "maximal number of successfully solved red. cost pricing problems",
         &maxsuccessfulprobs, FALSE, DEFAULT_MAXSUCCESSFULPROBSREDCOST, 1, INT_MAX, NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricing/masterpricer/relmaxprobsredcostroot",
         "maximal percentage of red. cost pricing problems that are solved at root node if variables have already been found",
         &relmaxprobsroot, FALSE, DEFAULT_RELMAXPROBSREDCOSTROOT, 0.0, 1.0, NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricing/masterpricer/relmaxprobsredcost",
         "maximal percentage of red. cost pricing problems that are solved if variables have already been found",
         &relmaxprobs, FALSE, DEFAULT_RELMAXPROBSREDCOST, 0.0, 1.0, NULL, (SCIP_PARAMDATA*) NULL) );

   return SCIP_OKAY;
}

SCIP_Bool FarkasPricing::canOptimalPricingBeAborted(
   int                   nfoundcols,         /**< number of negative reduced cost columns found so far */
   int                   nsolvedprobs,       /**< number of pricing problems solved so far */
   int                   nsuccessfulprobs,   /**< number of pricing problems solved successfully so far */
   SCIP_Real             relmaxsuccessfulprobs, /**< maximal percentage of pricing problems that need to be solved successfully */
   int                   npricingprobsnotnull /**< number of relevant pricing problems */
   ) const
{ /*lint -esym(715,successfulprobs,successfulrelmaxprobs) */
   return !((nfoundcols < maxcolsround)
            && (nfoundcols == 0 || nsolvedprobs < relmaxprobs * npricingprobsnotnull));
}

SCIP_Bool FarkasPricing::canHeuristicPricingBeAborted(
   int                   nfoundcols,         /**< number of negative reduced cost columns found so far */
   int                   nsolvedprobs,       /**< number of pricing problems solved so far */
   int                   nsuccessfulprobs,   /**< number of pricing problems solved successfully so far */
   SCIP_Real             relmaxsuccessfulprobs, /**< maximal percentage of pricing problems that need to be solved successfully */
   int                   npricingprobsnotnull /**< number of relevant pricing problems */
   ) const
{
   return canOptimalPricingBeAborted(nfoundcols, nsolvedprobs, nsuccessfulprobs, relmaxsuccessfulprobs, npricingprobsnotnull);
}

SCIP_Bool ReducedCostPricing::canOptimalPricingBeAborted(
   int                   nfoundcols,         /**< number of negative reduced cost columns found so far */
   int                   nsolvedprobs,       /**< number of pricing problems solved so far */
   int                   nsuccessfulprobs,   /**< number of pricing problems solved successfully so far */
   SCIP_Real             relmaxsuccessfulprobs, /**< maximal percentage of pricing problems that need to be solved successfully */
   int                   npricingprobsnotnull /**< number of relevant pricing problems */
   ) const
{
   return !((((nfoundcols < maxcolsroundroot) || !GCGisRootNode(scip_)) && ((nfoundcols < maxcolsround) || GCGisRootNode(scip_)))
               && nsuccessfulprobs < maxsuccessfulprobs
               && nsuccessfulprobs < relmaxsuccessfulprobs * npricingprobsnotnull
               && (nfoundcols == 0 || ((GCGisRootNode(scip_) || nsolvedprobs < relmaxprobs * npricingprobsnotnull)
                     && (!GCGisRootNode(scip_) || nsolvedprobs < relmaxprobsroot * npricingprobsnotnull))));
}

SCIP_Bool ReducedCostPricing::canHeuristicPricingBeAborted(
   int                   nfoundcols,         /**< number of negative reduced cost columns found so far */
   int                   nsolvedprobs,       /**< number of pricing problems solved so far */
   int                   nsuccessfulprobs,   /**< number of pricing problems solved successfully so far */
   SCIP_Real             relmaxsuccessfulprobs, /**< maximal percentage of pricing problems that need to be solved successfully */
   int                   npricingprobsnotnull /**< number of relevant pricing problems */
  ) const
{
   return !((nfoundcols < maxcolsround)
            && nsuccessfulprobs < maxsuccessfulprobs
            && nsuccessfulprobs < relmaxsuccessfulprobs * npricingprobsnotnull
            && (nfoundcols == 0 || nsolvedprobs < relmaxprobs * npricingprobsnotnull));
}
