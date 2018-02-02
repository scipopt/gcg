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

#define DEFAULT_MAXVARSROUNDREDCOSTROOT  100        /**< maximal number of variables per reduced cost pricing round at root node */
#define DEFAULT_MAXVARSROUNDREDCOST      100        /**< maximal number of variables per reduced cost pricing round */
#define DEFAULT_MAXSUCCESSFULMIPSREDCOST INT_MAX    /**< maximal number of successful MIP solves */
#define DEFAULT_MAXROUNDSREDCOST         INT_MAX    /**< maximal number of reduced cost pricing rounds */
#define DEFAULT_MIPSRELREDCOSTROOT       1.0        /**< factor of reduced cost pricing MIPs to be solved at root node */
#define DEFAULT_MIPSRELREDCOST           1.0        /**< factor of reduced cost pricing MIPs to be solver */

#define DEFAULT_MAXVARSROUNDFARKAS       10         /**< maximal number of variables per farkas pricing round */
#define DEFAULT_MIPSRELFARKAS            1.0        /**< factor of farkas pricing MIPs to be solved */


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
      SCIP* scip
      )
{
   type  = GCG_PRICETYPE_UNKNOWN;
   scip_ = scip;
   calls = INT_MAX;
   maxvarsround = INT_MAX;
   maxvarsroundroot = INT_MAX;
   maxsuccessfulmips = INT_MAX;
   maxrounds = INT_MAX;

   mipsrel = 1.0;
   mipsrelroot = 1.0;

   SCIP_CALL_EXC( SCIPcreateCPUClock(scip, &(clock)) );
   calls = 0;
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
   SCIP_CALL( SCIPaddIntParam(GCGmasterGetOrigprob(scip_), "pricing/masterpricer/maxvarsroundfarkas",
         "maximal number of variables created in one farkas pricing round",
         &maxvarsround, FALSE, DEFAULT_MAXVARSROUNDFARKAS, 1, INT_MAX, NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddRealParam(GCGmasterGetOrigprob(scip_), "pricing/masterpricer/mipsrelfarkas",
         "part of the submips that are solved before Farkas pricing round is aborted, if variables have been found yet? (1.0 = solve all pricing MIPs)",
         &mipsrel, FALSE, DEFAULT_MIPSRELFARKAS, 0.0, 1.0, NULL, (SCIP_PARAMDATA*) NULL) );

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
   SCIP_CALL( SCIPaddIntParam(GCGmasterGetOrigprob(scip_), "pricing/masterpricer/maxsuccessfulmipsredcost",
         "maximal number of pricing mips leading to new variables solved solved in one redcost pricing round",
         &maxsuccessfulmips, FALSE, DEFAULT_MAXSUCCESSFULMIPSREDCOST, 1, INT_MAX, NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddIntParam(GCGmasterGetOrigprob(scip_), "pricing/masterpricer/maxvarsroundredcost",
         "maximal number of variables created in one redcost pricing round",
         &maxvarsround, FALSE, DEFAULT_MAXVARSROUNDREDCOST, 0, INT_MAX,
         NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddIntParam(GCGmasterGetOrigprob(scip_), "pricing/masterpricer/maxvarsroundredcostroot",
         "maximal number of variables created in one redcost pricing round at the root node",
         &maxvarsroundroot, FALSE, DEFAULT_MAXVARSROUNDREDCOSTROOT, 0, INT_MAX,
         NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddIntParam(GCGmasterGetOrigprob(scip_), "pricing/masterpricer/maxroundsredcost",
         "maximal number of pricing rounds per node after the root node",
         &maxrounds, FALSE, DEFAULT_MAXROUNDSREDCOST, 0, INT_MAX, NULL, (SCIP_PARAMDATA*) NULL) );


   SCIP_CALL( SCIPaddRealParam(GCGmasterGetOrigprob(scip_), "pricing/masterpricer/mipsrelredcostroot",
         "part of the submips that are solved before redcost pricing round is aborted at the root node, if variables have been found yed? (1.0 = solve all pricing MIPs)",
         &mipsrelroot, FALSE, DEFAULT_MIPSRELREDCOSTROOT, 0.0, 1.0, NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddRealParam(GCGmasterGetOrigprob(scip_), "pricing/masterpricer/mipsrelredcost",
         "part of the submips that are solved before redcost pricing round is aborted, if variables have been found yed? (1.0 = solve all pricing MIPs)",
         &mipsrel, FALSE, DEFAULT_MIPSRELREDCOST, 0.0, 1.0, NULL, (SCIP_PARAMDATA*) NULL) );

   return SCIP_OKAY;
}

SCIP_Bool FarkasPricing::canOptimalPricingBeAborted(
      int                  nfoundvars,         /**< number of variables found so far */
      int                  solvedmips,         /**< number of MIPS solved so far */
      int                  successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real            successfulmipsrel,     /**< number of sucessful mips solved so far */
      int                  npricingprobsnotnull
      ) const
{ /*lint -esym(715,successfulmips,successfulmipsrel) */
   return !(nfoundvars < maxvarsround && (nfoundvars == 0 || solvedmips < mipsrel * npricingprobsnotnull));
}

SCIP_Bool FarkasPricing::canHeuristicPricingBeAborted(
      int                  nfoundvars,         /**< number of variables found so far */
      int                  solvedmips,         /**< number of MIPS solved so far */
      int                  successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real            successfulmipsrel,     /**< number of sucessful mips solved so far */
      int                  npricingprobsnotnull
      ) const
{
   return canOptimalPricingBeAborted(nfoundvars, solvedmips, successfulmips, successfulmipsrel, npricingprobsnotnull);
}

SCIP_Bool ReducedCostPricing::canOptimalPricingBeAborted(
      int                  nfoundvars,         /**< number of variables found so far */
      int                  solvedmips,         /**< number of MIPS solved so far */
      int                  successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real            successfulmipsrel,     /**< number of sucessful mips solved so far */
      int                  npricingprobsnotnull
  ) const
{

   return !((((nfoundvars < maxvarsroundroot) || !GCGisRootNode(scip_) ) && ((nfoundvars < maxvarsround) || GCGisRootNode(scip_)))
               && successfulmips < maxsuccessfulmips
               && successfulmips < successfulmipsrel * npricingprobsnotnull
               && (nfoundvars == 0 || ( (GCGisRootNode(scip_) || solvedmips < mipsrel * npricingprobsnotnull)
                     && (!GCGisRootNode(scip_) || solvedmips < mipsrelroot * npricingprobsnotnull))));
}

SCIP_Bool ReducedCostPricing::canHeuristicPricingBeAborted(
      int                  nfoundvars,         /**< number of variables found so far */
      int                  solvedmips,         /**< number of MIPS solved so far */
      int                  successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real            successfulmipsrel,     /**< number of sucessful mips solved so far */
      int                  npricingprobsnotnull
  ) const
{
   return !((nfoundvars < maxvarsround)
            && successfulmips < maxsuccessfulmips
            && successfulmips < successfulmipsrel * npricingprobsnotnull
            && (nfoundvars == 0 ||
               solvedmips < mipsrel * npricingprobsnotnull ));
}
