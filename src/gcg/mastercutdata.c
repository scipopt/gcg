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
#include "struct_mastercutdata.h"

#include <scip/cons_linear.h>
#include <scip/def.h>
#include <scip/pub_cons.h>
#include <scip/pub_lp.h>
#include <scip/scip.h>
#include <scip/type_scip.h>

/**
 * @ingroup TODO-????
 *
 * @{
 */

/** determine whether the mastercutdata is active in the masterscip */
SCIP_Bool GCGmastercutIsActive(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   switch( mastercutdata->type )
   {
   case GCG_MASTERCUTTYPE_CONS:
      assert(mastercutdata->cut.cons != NULL);
      return SCIPconsIsActive(mastercutdata->cut.cons);
   case GCG_MASTERCUTTYPE_ROW:
      assert(mastercutdata->cut.row != NULL);
      return SCIProwIsInLP(mastercutdata->cut.row);
   }
}

/** add a new variable along with its coefficient to the mastercut */
SCIP_RETCODE GCGmastercutAddVar(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_VAR*              var,                /**< variable to add */
   SCIP_Real              coef                /**< coefficient of the variable */
   )
{
   assert(masterscip != NULL);
   assert(mastercutdata != NULL);
   assert(var != NULL);

   switch( mastercutdata->type )
   {
   case GCG_MASTERCUTTYPE_CONS:
      assert(mastercutdata->cut.cons != NULL);
      SCIP_CALL( SCIPaddCoefLinear(masterscip, mastercutdata->cut.cons, var, coef) );
      break;
   case GCG_MASTERCUTTYPE_ROW:
      assert(mastercutdata->cut.row != NULL);
      SCIP_CALL( SCIPaddVarToRow(masterscip, mastercutdata->cut.row, var, coef) );
      break;
   }

   return SCIP_OKAY;
}

/** get the constraint that is the master cut
  * will fail if the master cut is a row
  */
SCIP_RETCODE GCGmastercutGetCons(
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_CONS**            cons                /**< pointer to store the constraint */
   )
{
   assert(mastercutdata != NULL);
   assert(cons != NULL);
   assert(*cons == NULL);

   if( mastercutdata->type != GCG_MASTERCUTTYPE_CONS )
      return SCIP_ERROR;

   assert(mastercutdata->cut.cons != NULL);
   *cons = mastercutdata->cut.cons;

   return SCIP_OKAY;
}

/** get the row that is the master cut
   * will fail if the master cut is a constraint
   */
SCIP_RETCODE GCGmastercutGetRow(
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_ROW**             row                 /**< pointer to store the row */
   )
{
   assert(mastercutdata != NULL);
   assert(row != NULL);
   assert(*row == NULL);

   if( mastercutdata->type != GCG_MASTERCUTTYPE_ROW )
      return SCIP_ERROR;

   assert(mastercutdata->cut.row != NULL);
   *row = mastercutdata->cut.row;

   return SCIP_OKAY;
}
