/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2012 Operations Research, RWTH Aachen University       */
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
 * @brief  GCG variable pricer
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "class_pricingtype.h"
#include "scip/cons_linear.h"
#include "pub_gcgvar.h"
#include "scip/pub_lp.h"

PricingType::PricingType()
{

}

FarkasPricing::FarkasPricing()
{
   type = GCG_PRICETYPE_FARKAS;
}

double FarkasPricing::consGetDual(
      SCIP* scip,
      SCIP_CONS* cons
   )
{
   return SCIPgetDualfarkasLinear(scip, cons);
}

double ReducedCostPricing::consGetDual(
      SCIP* scip,
      SCIP_CONS* cons
   )
{
   return SCIPgetDualsolLinear(scip, cons);
}

double FarkasPricing::rowGetDual(
      SCIP_Row* row
   )
{
   return SCIProwGetDualfarkas(row);
}

double ReducedCostPricing::rowGetDual(
      SCIP_ROW* row
   )
{
   return SCIProwGetDualsol(row);
}

ReducedCostPricing::ReducedCostPricing()
{
   type = GCG_PRICETYPE_REDCOST;
}

SCIP_Real ReducedCostPricing::varGetObj(
      SCIP_VAR* var
   )
{
   SCIP_VAR* origvar;
   assert(var != NULL);

   origvar = GCGpricingVarGetOrigvars(var)[0];

   if( GCGvarIsLinking(origvar) )
      return 0.0;
   else
      return SCIPvarGetObj(origvar);
}

SCIP_Real FarkasPricing::varGetObj(
      SCIP_VAR* var
   )
{
   assert(var != NULL);
   return 0.0;
}

SCIP_RETCODE PricingType::startClock(
      SCIP* scip
   )
{
   SCIP_CALL( SCIPstartClock(scip, clock) );
   return SCIP_OKAY;
}

SCIP_RETCODE PricingType::stopClock(
      SCIP* scip
   )
{
   SCIP_CALL( SCIPstopClock(scip, clock) );
   return SCIP_OKAY;
}
