/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file    mastercutdata.c
 * @ingroup TODO-????
 * @brief   methods for interacting with GCG_MASTERCUTDATA
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "def.h"
#include "mastercutdata.h"
#include "gcg.h"
#include "pub_gcgvar.h"
#include "type_mastercutdata.h"

#include <scip/cons_linear.h>
#include <scip/def.h>
#include <scip/pub_lp.h>
#include <scip/scip.h>
#include <scip/type_scip.h>

/**
 * @ingroup TODO-????
 *
 * @{
 */

/** get the blocknr of a mastercut */
GCG_EXPORT
int GCGmastercutGetBlocknr(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);
   int blocknr = mastercutdata->blocknr;
   assert(blocknr >= 0);
   return blocknr;
}

/** determine whether the mastercutdata is active in the masterscip */
GCG_EXPORT
SCIP_Bool GCGmastercutIsActive(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   if( SCIProwIsInLP(SCIPgetRowLinear(masterscip, mastercutdata->mastercons)) )
      return TRUE;

   return FALSE;
}

/** activate the mastercutdata */
GCG_EXPORT
SCIP_RETCODE GCGmastercutActivate(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(masterscip != NULL);
   assert(mastercutdata != NULL);

   SCIP_CALL( SCIPaddCons(masterscip, mastercutdata->mastercons) );

   return SCIP_OKAY;
}

/** deactivate the mastercutdata */
GCG_EXPORT
SCIP_RETCODE GCGmastercutDeactivate(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(masterscip != NULL);
   assert(mastercutdata != NULL);

   SCIP_CALL( SCIPdelCons(masterscip, mastercutdata->mastercons) );

   return SCIP_OKAY;
}

/** get all referenced pricing variables */
GCG_EXPORT
SCIP_RETCODE GCGmastercutGetAllReferencedPricingVariables(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_VAR***            vars,               /**< pointer to store the array of variables */
   int*                   nvars               /**< pointer to store the number of variables */
   )
{
   int i;
   int j;
   SCIP* pricingscip;

   assert(masterscip != NULL);
   assert(mastercutdata != NULL);
   assert(vars != NULL);
   assert(*vars == NULL);
   assert(nvars != NULL);

   pricingscip = GCGgetPricingprob(masterscip, mastercutdata->blocknr);

   *nvars = mastercutdata->npricingvars;
   if( *nvars == 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, vars, nvars) );
      for( i = 0; i < *nvars; i++ )
      {
         assert(GCGvarIsInferredPricing(mastercutdata->pricingvars[i]));
         (*vars)[i] = mastercutdata->pricingvars[i];
      }
   }

   for( i = 0; i < mastercutdata->npricingconss; i++ )
   {
      SCIP_CONS* pricingcons = mastercutdata->pricingconss[i];
      SCIP_VAR** pricingconsvars;
      int npricingconsvars;
      SCIP_Bool success = FALSE;

      SCIP_CALL( SCIPgetConsNVars(pricingscip, pricingcons, &npricingconsvars, &success) );
      assert(success);
      success = FALSE;
      SCIP_CALL( SCIPallocBufferArray(masterscip, &pricingconsvars, npricingconsvars) );
      SCIP_CALL( SCIPgetConsVars(pricingscip, pricingcons, pricingconsvars, npricingconsvars, &success) );
      assert(success);

      for( j = 0; j < npricingconsvars; j++ )
      {
         SCIP_VAR* var = pricingconsvars[j];
         assert(GCGvarIsInferredPricing(var) || GCGvarIsPricing(var));

         // skip inferred pricing variables, already added
         if( GCGvarIsInferredPricing(var) )
            continue;

         if( *nvars == 0 )
            SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, vars, nvars) );
         else
            SCIP_CALL( SCIPreallocBlockMemoryArray(masterscip, vars, *nvars, *nvars + 1) );

         (*vars)[*nvars] = var;
         (*nvars)++;
      }

      SCIPfreeBufferArray(masterscip, &pricingconsvars);
   }

   return SCIP_OKAY;
}

/** determine whether a constraint is a given mastercut */
GCG_EXPORT
SCIP_Bool GCGmastercutIsConstraint(
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_CONS*             cons                /**< constraint */
   )
{
   assert(mastercutdata != NULL);
   assert(cons != NULL);

   return (mastercutdata->mastercons == cons);
}

/** determine wheter a constraint is any given mastercut */
GCG_EXPORT
SCIP_Bool GCGmastercutIsAnyConstraint(
   GCG_MASTERCUTDATA**    mastercutdata,      /**< array of mastercut data */
   int                    nmastercutdata,     /**< number of mastercut data */
   SCIP_CONS*             cons                /**< constraint */
   )
{
   int i;

   assert(mastercutdata != NULL);
   assert(nmastercutdata >= 0);
   assert(cons != NULL);

   for( i = 0; i < nmastercutdata; i++ )
   {
      if( GCGmastercutIsConstraint(mastercutdata[i], cons) )
         return TRUE;
   }

   return FALSE;
}
