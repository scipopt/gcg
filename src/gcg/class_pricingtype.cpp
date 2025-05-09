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

/**@file   class_pricingtype.cpp
 * @brief  abstraction for SCIP pricing types
 * @author Martin Bergner
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "gcg/class_pricingtype.h"
#include "scip/cons_linear.h"
#include "gcg/pub_gcgvar.h"
#include "gcg/pub_extendedmasterconsdata.h"
#include "scip/pub_lp.h"
#include "scip/clock.h"
#include "gcg/scip_misc.h"

#include <exception>

#define DEFAULT_MAXROUNDSREDCOST         INT_MAX    /**< maximal number of reduced cost pricing rounds */
#define DEFAULT_MAXCOLSROUNDREDCOSTROOT  100        /**< maximal number of columns per reduced cost pricing round at root node */
#define DEFAULT_MAXCOLSROUNDREDCOST      100        /**< maximal number of columns per reduced cost pricing round */
#define DEFAULT_MAXCOLSPROBREDCOSTROOT    10        /**< maximal number of columns per problem to be generated during red. cost pricing at root node */
#define DEFAULT_MAXCOLSPROBREDCOST        10        /**< maximal number of columns per problem to be generated during red. cost pricing */
#define DEFAULT_MAXSUCCESSFULPROBSREDCOST INT_MAX   /**< maximal number of successfully solved red. cost pricing problems until pricing loop is aborted */
#define DEFAULT_RELMAXPROBSREDCOSTROOT   1.0        /**< maximal percentage of red. cost pricing problems that are solved at root node if variables have already been found */
#define DEFAULT_RELMAXPROBSREDCOST       1.0        /**< maximal percentage of red. cost pricing problems that are solved if variables have already been found */
#define DEFAULT_RELMAXSUCCESSFULPROBSREDCOST 1.0    /**< maximal percentage of successfully solved red. cost pricing problems until pricing loop is aborted */

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
   GCG*                 gcgstruct
   ) : gcg(gcgstruct)
{
   masterprob = GCGgetMasterprob(gcg);       /* SCIP instance (master problem) */
   type  = GCG_PRICETYPE_UNKNOWN;            /* type of pricing */

   /* statistical values */
   calls = 0;                                /* number of times this type of pricing was called */

   /* parameters */
   maxrounds = INT_MAX;                      /* maximal number of pricing rounds */
   maxcolsroundroot = INT_MAX;               /* maximal number of columns per pricing round at root node */
   maxcolsround = INT_MAX;                   /* maximal number of columns per pricing round */
   maxcolsprobroot = INT_MAX;                /* maximal number of columns per problem to be generated at root node */
   maxcolsprob = INT_MAX;                    /* maximal number of columns per problem to be generated */
   maxsuccessfulprobs = INT_MAX;             /* maximal number of successfully solved pricing problems until pricing loop is aborted */
   relmaxprobsroot = 1.0;                    /* maximal percentage of pricing problems that are solved at root node if variables have already been found */
   relmaxprobs = 1.0;                        /* maximal percentage of pricing problems that are solved if variables have already been found */
   relmaxsuccessfulprobs = 1.0;              /* maximal percentage of successfully solved pricing problems until pricing loop is aborted */

   SCIP_CALL_EXC( SCIPcreateCPUClock(masterprob, &(clock)) );
}

PricingType::~PricingType()
{
   SCIP_CALL_ABORT( SCIPfreeClock(masterprob, &(clock)) );

   masterprob = (SCIP*) NULL;
}

SCIP_RETCODE PricingType::startClock()
{
   SCIP_CALL( SCIPstartClock(masterprob, clock) );
   return SCIP_OKAY;
}

SCIP_RETCODE PricingType::stopClock()
{
   SCIP_CALL( SCIPstopClock(masterprob, clock) );
   return SCIP_OKAY;
}

SCIP_Real PricingType::getClockTime() const
{
   return SCIPgetClockTime(masterprob, clock);
}

FarkasPricing::FarkasPricing(
   GCG*                  gcgstruct
   ) : PricingType(gcgstruct)
{
   type = GCG_PRICETYPE_FARKAS;
}

SCIP_Real FarkasPricing::consGetDual(
   SCIP_CONS*            cons
   ) const
{
   return SCIPgetDualfarkasLinear(masterprob, cons);
}

SCIP_Real FarkasPricing::rowGetDual(
   SCIP_ROW*             row
   ) const
{
   return SCIProwGetDualfarkas(row);
}

SCIP_Real FarkasPricing::extendedmasterconsGetDual(
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata
   ) const
{
   assert(extendedmasterconsdata != NULL);
   switch( GCGextendedmasterconsGetType(extendedmasterconsdata) )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_BRANCH_CONS:
      return SCIPgetDualfarkasLinear(masterprob, GCGextendedmasterconsGetCons(extendedmasterconsdata));
   case GCG_EXTENDEDMASTERCONSTYPE_SEPA_ROW:
      return SCIProwGetDualfarkas(GCGextendedmasterconsGetRow(extendedmasterconsdata));
   default:
      SCIP_CALL_ABORT( SCIP_ERROR );
      return 0.0;
   }
}

SCIP_Real FarkasPricing::varGetObj(
   SCIP_VAR*             var
   ) const
{
   assert(var != NULL);
   return 0.0;
}

/** returns the maximal number of columns per pricing round */
int FarkasPricing::getMaxcolsround() const
{
   return maxcolsround;
}

/** returns the maximal number of columns per problem to be generated during pricing */
int FarkasPricing::getMaxcolsprob() const
{
   return maxcolsprob;
}

/** returns the maximal percentage of pricing problems that are solved if variables have already been found */
SCIP_Real FarkasPricing::getRelmaxprobs() const
{
   return relmaxprobs;
}

SCIP_RETCODE FarkasPricing::addParameters()
{
   SCIP* origprob = GCGgetOrigprob(gcg);

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
   SCIP_CONS*            cons
   ) const
{
   return SCIPgetDualsolLinear(masterprob, cons);
}

SCIP_Real ReducedCostPricing::rowGetDual(
   SCIP_ROW*             row
   ) const
{
   return SCIProwGetDualsol(row);
}

SCIP_Real ReducedCostPricing::extendedmasterconsGetDual(
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata
   ) const
{
   assert(extendedmasterconsdata != NULL);
   switch( GCGextendedmasterconsGetType(extendedmasterconsdata) )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_BRANCH_CONS:
      return SCIPgetDualsolLinear(masterprob, GCGextendedmasterconsGetCons(extendedmasterconsdata));
   case GCG_EXTENDEDMASTERCONSTYPE_SEPA_ROW:
      return SCIProwGetDualsol(GCGextendedmasterconsGetRow(extendedmasterconsdata));
   default:
      SCIP_CALL_ABORT( SCIP_ERROR );
      return 0.0;
   }
}

ReducedCostPricing::ReducedCostPricing(
   GCG*                  gcgstruct
   ) : PricingType(gcgstruct)
{
   type = GCG_PRICETYPE_REDCOST;
}

SCIP_Real ReducedCostPricing::varGetObj(
   SCIP_VAR*             var
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

/** returns the maximal number of columns per pricing round */
int ReducedCostPricing::getMaxcolsround() const
{
   return GCGisRootNode(gcg) ? maxcolsroundroot : maxcolsround;
}

/** returns the maximal number of columns per problem to be generated during pricing */
int ReducedCostPricing::getMaxcolsprob() const
{
   return GCGisRootNode(gcg) ? maxcolsprobroot : maxcolsprob;
}

/** returns the maximal percentage of pricing problems that are solved if variables have already been found */
SCIP_Real ReducedCostPricing::getRelmaxprobs() const
{
   return GCGisRootNode(gcg) ? relmaxprobsroot : relmaxprobs;
}

SCIP_RETCODE ReducedCostPricing::addParameters()
{
   SCIP* origprob = GCGgetOrigprob(gcg);

   SCIP_CALL( SCIPaddIntParam(origprob, "pricing/masterpricer/maxroundsredcost",
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
         "maximal number of successfully solved red. cost pricing problems until pricing loop is aborted",
         &maxsuccessfulprobs, FALSE, DEFAULT_MAXSUCCESSFULPROBSREDCOST, 1, INT_MAX, NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricing/masterpricer/relmaxprobsredcostroot",
         "maximal percentage of red. cost pricing problems that are solved at root node if variables have already been found",
         &relmaxprobsroot, FALSE, DEFAULT_RELMAXPROBSREDCOSTROOT, 0.0, 1.0, NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricing/masterpricer/relmaxprobsredcost",
         "maximal percentage of red. cost pricing problems that are solved if variables have already been found",
         &relmaxprobs, FALSE, DEFAULT_RELMAXPROBSREDCOST, 0.0, 1.0, NULL, (SCIP_PARAMDATA*) NULL) );

   SCIP_CALL( SCIPaddRealParam(origprob, "pricing/masterpricer/relmaxsuccessfulprobsredcost",
         "maximal percentage of successfully solved red. cost pricing problems until pricing loop is aborted",
         &relmaxsuccessfulprobs, FALSE, DEFAULT_RELMAXSUCCESSFULPROBSREDCOST, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
