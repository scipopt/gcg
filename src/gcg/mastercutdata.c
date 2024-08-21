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
 * @brief   methods for interacting with GCG_MASTERCUTDATA
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "def.h"
#include "mastercutdata.h"
#include "gcg.h"
#include "gcg/scip_misc.h"
#include "pricer_gcg.h"
#include "scip/pub_var.h"
#include "struct_mastercutdata.h"

#include <scip/cons_linear.h>
#include <scip/def.h>
#include <scip/pub_cons.h>
#include <scip/pub_lp.h>
#include <scip/pub_misc_linear.h>
#include <scip/scip.h>
#include <scip/scip_cons.h>
#include <scip/scip_prob.h>
#include <scip/struct_scip.h>
#include <scip/struct_mem.h>
#include <scip/struct_var.h>
#include <scip/type_scip.h>

/** free a pricing modification */
static
SCIP_RETCODE GCGpricingmodificationFree(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_PRICINGMODIFICATION* pricingmodification /**< pointer to the pricing modification */
   )
{
   SCIP* originalscip;
   SCIP* pricingscip;
   int i;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(pricingmodification != NULL);

   originalscip = GCGmasterGetOrigprob(scip);
   pricingscip = GCGgetPricingprob(originalscip, pricingmodification->blocknr);

   SCIP_CALL( SCIPreleaseVar(pricingscip, &pricingmodification->coefvar) );

   for( i = 0; i < pricingmodification->nadditionalvars; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(pricingscip, &pricingmodification->additionalvars[i]) );
   }

   for( i = 0; i < pricingmodification->nadditionalconss; i++ )
   {
      SCIP_CALL( SCIPreleaseCons(pricingscip, &pricingmodification->additionalconss[i]) );
   }

   SCIPfreeBlockMemoryArrayNull(pricingscip, &pricingmodification->additionalvars, pricingmodification->nadditionalvars);
   assert(pricingmodification->additionalvars == NULL);
   SCIPfreeBlockMemoryArrayNull(pricingscip, &pricingmodification->additionalconss, pricingmodification->nadditionalconss);
   assert(pricingmodification->additionalconss == NULL);


   return SCIP_OKAY;
}

/** create a pricing modification, taking ownership over additionalvars and additionalcons */
GCG_EXPORT
SCIP_RETCODE GCGpricingmodificationCreate(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_PRICINGMODIFICATION* pricingmodification, /**< pointer to store the pricing modification */
   int                    blocknr,             /**< block number of the master cut */
   SCIP_VAR*              coefvar,             /**< variable in the pricing problem inferred from the master cut
                                                  * always has the objective coefficient of the negated dual value of the master cut
                                                  * its solution value corresponds to the coefficient of the new mastervariable in the master cut */
   SCIP_VAR**             additionalvars,      /**< array of additional variables with no objective coefficient in the pricing programs inferred from the master cut */
   int                    nadditionalvars,     /**< number of additional variables in the pricing programs */
   SCIP_CONS**            additionalconss,     /**< array of additional constraints in the pricing programs inferred from the master cut */
   int                    nadditionalconss     /**< number of additional constraints in the pricing programs */
   )
{
   int i;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(pricingmodification != NULL);
   assert(blocknr >= 0);

   assert(blocknr < GCGgetNPricingprobs(GCGgetOriginalprob(scip)));

   assert(coefvar != NULL);
   assert(GCGvarIsInferredPricing(coefvar));
   assert(additionalvars != NULL || nadditionalvars == 0);
   assert(additionalconss != NULL || nadditionalconss == 0);

   for( i = 0; i < nadditionalvars; i++ )
   {
      assert(additionalvars[i] != NULL);
      assert(additionalvars[i] != coefvar);
      assert(GCGvarIsInferredPricing(additionalvars[i]));
      assert(SCIPisZero(scip, additionalvars[i]->obj));
   }

   for( i = 0; i < nadditionalconss; i++ )
   {
      assert(additionalconss[i] != NULL);
   }

   pricingmodification->blocknr = blocknr;
   pricingmodification->coefvar = coefvar;
   pricingmodification->additionalvars = additionalvars;
   pricingmodification->nadditionalvars = nadditionalvars;
   pricingmodification->additionalconss = additionalconss;
   pricingmodification->nadditionalconss = nadditionalconss;

   return SCIP_OKAY;
}

/** create a master cut, taking ownership over pricingmodifications */
GCG_EXPORT
SCIP_RETCODE GCGmastercutCreateFromCons(
   SCIP*                  scip,                 /**< SCIP data structure */
   GCG_MASTERCUTDATA**    mastercutdata,        /**< pointer to store the mastercut data */
   SCIP_CONS*             cons,                 /**< constraint in the master problem that represents the master cut */
   GCG_PRICINGMODIFICATION* pricingmodifications,/**< pricing modifications for the master cut */
   int                    npricingmodifications,/**< number of pricing modifications for the master cut */
   void*                  data,                 /**< any data that might be required to calculate the coefficient of a column solution */
   GCG_DECL_MASTERCUTGETCOEFF((*mastercutGetCoeff))/**< callback to calculate the coefficient of a column solution */
   )
{
#ifndef NDEBUG
   SCIP_Bool* seenblocks;
   SCIP* originalproblem;
#endif
   int i;
   int j;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(mastercutdata != NULL);
   assert(*mastercutdata == NULL);
   assert(cons != NULL);
   assert(pricingmodifications != NULL || npricingmodifications == 0);

#ifndef NDEBUG
   originalproblem = GCGgetOriginalprob(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &seenblocks, GCGgetNRelPricingprobs(originalproblem)) );
   for( i = 0; i < GCGgetNRelPricingprobs(originalproblem); i++ )
      seenblocks[i] = FALSE;
#endif

   for( i = 0; i < npricingmodifications; i++ )
   {
      assert(pricingmodifications[i].blocknr >= 0);
      assert(pricingmodifications[i].blocknr < GCGgetNRelPricingprobs(originalproblem));
#ifndef NDEBUG
      assert(!seenblocks[pricingmodifications[i].blocknr]);
      seenblocks[pricingmodifications[i].blocknr] = TRUE;
#endif
   }

#ifndef NDEBUG
   SCIPfreeBlockMemoryArray(scip, &seenblocks, GCGgetNRelPricingprobs(originalproblem));
   assert(seenblocks == NULL);
#endif

   SCIP_CALL( SCIPallocBlockMemory(scip, mastercutdata) );

   (*mastercutdata)->type = GCG_MASTERCUTTYPE_CONS;
   (*mastercutdata)->cut.cons = cons;
   (*mastercutdata)->pricingmodifications = pricingmodifications;
   (*mastercutdata)->npricingmodifications = npricingmodifications;
   (*mastercutdata)->data = data;
   (*mastercutdata)->mastercutGetCoeff = mastercutGetCoeff;

   for( i = 0; i < npricingmodifications; i++ ) {
      pricingmodifications[i].coefvar->vardata->data.inferredpricingvardata.mastercutdata = *mastercutdata;
      for( j = 0; j < pricingmodifications[i].nadditionalvars; j++ ) {
         pricingmodifications[i].additionalvars[j]->vardata->data.inferredpricingvardata.mastercutdata = *mastercutdata;
      }
   }

   return SCIP_OKAY;
}

/** create a master cut, taking ownership over pricingmodifications */
GCG_EXPORT
SCIP_RETCODE GCGmastercutCreateFromRow(
   SCIP*                  scip,                 /**< SCIP data structure */
   GCG_MASTERCUTDATA**    mastercutdata,        /**< pointer to store the mastercut data */
   SCIP_ROW*              row,                  /**< row in the master problem that represents the master cut */
   GCG_PRICINGMODIFICATION* pricingmodifications,/**< pricing modifications for the master cut */
   int                    npricingmodifications,/**< number of pricing modifications for the master cut */
   void*                  data,                 /**< any data that might be required to calculate the coefficient of a column solution */
   GCG_DECL_MASTERCUTGETCOEFF((*mastercutGetCoeff))/**< callback to calculate the coefficient of a column solution */
   )
{
#ifndef NDEBUG
   SCIP_Bool* seenblocks;
   SCIP* originalproblem;
#endif
   int i;
   int j;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(mastercutdata != NULL);
   assert(*mastercutdata == NULL);
   assert(row != NULL);
   assert(pricingmodifications != NULL || npricingmodifications == 0);

#ifndef NDEBUG
   originalproblem = GCGgetOriginalprob(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &seenblocks, GCGgetNRelPricingprobs(originalproblem)) );
   for( i = 0; i < GCGgetNRelPricingprobs(originalproblem); i++ )
      seenblocks[i] = FALSE;
#endif

   for( i = 0; i < npricingmodifications; i++ )
   {
      assert(pricingmodifications[i].blocknr >= 0);
      assert(pricingmodifications[i].blocknr < GCGgetNRelPricingprobs(originalproblem));
#ifndef NDEBUG
      assert(!seenblocks[pricingmodifications[i].blocknr]);
      seenblocks[pricingmodifications[i].blocknr] = TRUE;
#endif
   }

#ifndef NDEBUG
   SCIPfreeBlockMemoryArray(scip, &seenblocks, GCGgetNRelPricingprobs(originalproblem));
#endif

   SCIP_CALL( SCIPallocBlockMemory(scip, mastercutdata) );

   (*mastercutdata)->type = GCG_MASTERCUTTYPE_ROW;
   (*mastercutdata)->cut.row = row;
   (*mastercutdata)->pricingmodifications = pricingmodifications;
   (*mastercutdata)->npricingmodifications = npricingmodifications;
   (*mastercutdata)->data = data;
   (*mastercutdata)->mastercutGetCoeff = mastercutGetCoeff;

   for( i = 0; i < npricingmodifications; i++ ) {
      pricingmodifications[i].coefvar->vardata->data.inferredpricingvardata.mastercutdata = *mastercutdata;
      for( j = 0; j < pricingmodifications[i].nadditionalvars; j++ ) {
         pricingmodifications[i].additionalvars[j]->vardata->data.inferredpricingvardata.mastercutdata = *mastercutdata;
      }
   }

   return SCIP_OKAY;
}

/** free a master cut */
GCG_EXPORT
SCIP_RETCODE GCGmastercutFree(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA**    mastercutdata        /**< pointer to the mastercut data */
   )
{
   int i;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(mastercutdata != NULL);
   assert(*mastercutdata != NULL);

   switch( (*mastercutdata)->type )
   {
   case GCG_MASTERCUTTYPE_CONS:
      assert((*mastercutdata)->cut.cons != NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &(*mastercutdata)->cut.cons) );
      break;
   case GCG_MASTERCUTTYPE_ROW:
      assert((*mastercutdata)->cut.row != NULL);
      SCIP_CALL( SCIPreleaseRow(scip, &(*mastercutdata)->cut.row) );
      break;
   default:
      SCIP_CALL( SCIP_ERROR );
   }

   for( i = 0; i < (*mastercutdata)->npricingmodifications; i++ )
   {
      SCIP_CALL( GCGpricingmodificationFree(scip, &(*mastercutdata)->pricingmodifications[i]) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*mastercutdata)->pricingmodifications, (*mastercutdata)->npricingmodifications);
   assert((*mastercutdata)->pricingmodifications == NULL);
   SCIPfreeBlockMemory(scip, mastercutdata);
   assert(*mastercutdata == NULL);

   return SCIP_OKAY;
}

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
   default:
      SCIP_CALL_ABORT( SCIP_ERROR );
      return FALSE;
   }
}

/** add a new variable along with its coefficient to the mastercut */
SCIP_RETCODE GCGmastercutAddMasterVar(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_VAR*              var,                /**< variable to add */
   SCIP_Real              coef                /**< coefficient of the variable */
   )
{
   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
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
   default:
      SCIP_CALL( SCIP_ERROR );
   }

   return SCIP_OKAY;
}

/** update the master cut with the new dual value */
SCIP_RETCODE GCGmastercutUpdateDualValue(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_Real              dualvalue           /**< dual value */
   )
{
   int i;
   SCIP* origscip;
   SCIP* pricingscip;

   assert(mastercutdata != NULL);
   assert(GCGisMaster(masterscip));

   origscip = GCGmasterGetOrigprob(masterscip);
   assert(origscip != NULL);

   for( i = 0; i < mastercutdata->npricingmodifications; i++ )
   {
      assert(mastercutdata->pricingmodifications[i].coefvar != NULL);
      assert(GCGvarIsInferredPricing(mastercutdata->pricingmodifications[i].coefvar));

      pricingscip = GCGgetPricingprob(origscip, mastercutdata->pricingmodifications[i].blocknr);
      assert(pricingscip != NULL);

      SCIP_CALL( SCIPchgVarObj(pricingscip, mastercutdata->pricingmodifications[i].coefvar, -dualvalue) );
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

/** get the variable that determines the coefficient of a column in the master cut */
SCIP_VAR* GCGpricingmodificationGetCoefVar(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   )
{
   assert(pricingmodification != NULL);

   return pricingmodification->coefvar;
}

/** get the additional variables that are inferred by the master cut */
SCIP_VAR** GCGpricingmodificationGetAdditionalVars(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   )
{
   assert(pricingmodification != NULL);

   return pricingmodification->additionalvars;
}

/** get the number of additional variables that are inferred by the master cut */
int GCGpricingmodificationGetNAdditionalVars(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   )
{
   assert(pricingmodification != NULL);

   return pricingmodification->nadditionalvars;
}

/** get the additional constraints that are inferred by the master cut */
SCIP_CONS** GCGpricingmodificationGetAdditionalConss(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   )
{
   assert(pricingmodification != NULL);

   return pricingmodification->additionalconss;
}

/** get the number of additional constraints that are inferred by the master cut */
int GCGpricingmodificationGetNAdditionalConss(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   )
{
   assert(pricingmodification != NULL);

   return pricingmodification->nadditionalconss;
}

/** get the pricing modification for a block, if exists, else NULL */
GCG_PRICINGMODIFICATION* GCGmastercutGetPricingModification(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   int                    blocknr             /**< block number */
   )
{
   int i;

   assert(mastercutdata != NULL);
   assert(GCGisMaster(masterscip));
   assert(blocknr >= 0);

   assert(blocknr < GCGgetNPricingprobs(GCGgetOriginalprob(masterscip)));

   for( i = 0; i < mastercutdata->npricingmodifications; i++ )
   {
      if( mastercutdata->pricingmodifications[i].blocknr == blocknr )
         return &mastercutdata->pricingmodifications[i];
   }

   return NULL;
}

/** get the pricing modifications for the master cut */
GCG_PRICINGMODIFICATION* GCGmastercutGetPricingModifications(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   return mastercutdata->pricingmodifications;
}

/** get the number of pricing modifications for the master cut */
int GCGmastercutGetNPricingModifications(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   return mastercutdata->npricingmodifications;
}

/** apply a pricing modification */
SCIP_RETCODE GCGpricingmodificationApply(
   SCIP*                  pricingscip,        /**< pricing scip */
   GCG_PRICINGMODIFICATION pricingmodification /**< pricing modification */
   )
{
   int i;

   assert(pricingscip != NULL);

   // add the inferred pricing variables
   assert(GCGvarIsInferredPricing(pricingmodification.coefvar));
   SCIP_CALL( SCIPaddVar(pricingscip, pricingmodification.coefvar) );

   for( i=0; i < pricingmodification.nadditionalvars; i++)
   {
      assert(GCGvarIsInferredPricing(pricingmodification.additionalvars[i]));
      SCIP_CALL( SCIPaddVar(pricingscip, pricingmodification.additionalvars[i]) );
   }

   // add the inferred pricing constraints
   for( i=0; i < pricingmodification.nadditionalconss; i++)
   {
      SCIP_CALL( SCIPaddCons(pricingscip, pricingmodification.additionalconss[i]) );
   }

   return SCIP_OKAY;
}

/** apply all pricing modifications */
SCIP_RETCODE GCGmastercutApplyPricingModifications(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   int i;
   SCIP* origscip;
   SCIP* pricingprob;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(mastercutdata != NULL);

   origscip = GCGmasterGetOrigprob(masterscip);
   assert(origscip != NULL);

   for( i = 0; i < mastercutdata->npricingmodifications; i++ )
   {
      pricingprob = GCGgetPricingprob(origscip, mastercutdata->pricingmodifications[i].blocknr);
      assert(pricingprob != NULL);
      SCIP_CALL( GCGpricingmodificationApply(pricingprob, mastercutdata->pricingmodifications[i]) );
   }

   return SCIP_OKAY;
}

/** undo a pricing modification */
SCIP_RETCODE GCGpricingmodificationUndo(
   SCIP*                  pricingscip,        /**< pricing scip */
   GCG_PRICINGMODIFICATION pricingmodification /**< pricing modification */
   )
{
   int i;
   SCIP_Bool deleted;

   assert(pricingscip != NULL);

   // add the inferred pricing variables
   assert(GCGvarIsInferredPricing(pricingmodification.coefvar));
   deleted = FALSE;
   SCIP_CALL( SCIPdelVar(pricingscip, pricingmodification.coefvar, &deleted) );
   assert(deleted);

   for( i=0; i < pricingmodification.nadditionalvars; i++)
   {
      assert(GCGvarIsInferredPricing(pricingmodification.additionalvars[i]));
      deleted = FALSE;
      SCIP_CALL( SCIPdelVar(pricingscip, pricingmodification.additionalvars[i], &deleted) );
      assert(deleted);
   }

   // add the inferred pricing constraints
   for( i=0; i < pricingmodification.nadditionalconss; i++)
   {
      SCIP_CALL( SCIPdelCons(pricingscip, pricingmodification.additionalconss[i]) );
   }

   return SCIP_OKAY;
}

/** undo all pricing modifications */
SCIP_RETCODE GCGmastercutUndoPricingModifications(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   int i;
   SCIP* origscip;
   SCIP* pricingprob;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(mastercutdata != NULL);

   origscip = GCGmasterGetOrigprob(masterscip);
   assert(origscip != NULL);

   for( i = 0; i < mastercutdata->npricingmodifications; i++ )
   {
      pricingprob = GCGgetPricingprob(origscip, mastercutdata->pricingmodifications[i].blocknr);
      assert(pricingprob != NULL);
      SCIP_CALL( GCGpricingmodificationUndo(pricingprob, mastercutdata->pricingmodifications[i]) );
   }

   return SCIP_OKAY;
}

/** check whether a given variable is a coefficient variable of a given pricing modification */
SCIP_Bool GCGpricingmodificationIsCoefVar(
   GCG_PRICINGMODIFICATION pricingmodification, /**< pricing modification */
   SCIP_VAR*              var                 /**< variable to check */
   )
{
   assert(var != NULL);

   return SCIPvarCompare(pricingmodification.coefvar, var) == 0;
}

/** check whether a given variable is a coefficient variable of a given mastercut */
SCIP_Bool GCGmastercutIsCoefVar(
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_VAR*              var                 /**< variable to check */
   )
{
   int i;

   assert(mastercutdata != NULL);
   assert(var != NULL);

   for( i = 0; i < mastercutdata->npricingmodifications; i++ )
   {
      if( GCGpricingmodificationIsCoefVar(mastercutdata->pricingmodifications[i], var) )
         return TRUE;
   }

   return FALSE;
}

/** get name of the mastercut */
const char* GCGmastercutGetName(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   switch( mastercutdata->type )
   {
   case GCG_MASTERCUTTYPE_CONS:
      assert(mastercutdata->cut.cons != NULL);
      return SCIPconsGetName(mastercutdata->cut.cons);
   case GCG_MASTERCUTTYPE_ROW:
      assert(mastercutdata->cut.row != NULL);
      return SCIProwGetName(mastercutdata->cut.row);
   default:
      SCIP_CALL_ABORT( SCIP_ERROR );
      return NULL;
   }
}

/** get the lhs of the mastercut */
SCIP_Real GCGmastercutGetLhs(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   switch( mastercutdata->type )
   {
   case GCG_MASTERCUTTYPE_CONS:
      assert(mastercutdata->cut.cons != NULL);
      return GCGconsGetLhs(scip, mastercutdata->cut.cons);
   case GCG_MASTERCUTTYPE_ROW:
      assert(mastercutdata->cut.row != NULL);
      return SCIProwGetLhs(mastercutdata->cut.row);
   default:
      SCIP_CALL_ABORT( SCIP_ERROR );
      return SCIP_INVALID;
   }
}

/** get the rhs of the mastercut */
SCIP_Real GCGmastercutGetRhs(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   switch( mastercutdata->type )
   {
   case GCG_MASTERCUTTYPE_CONS:
      assert(mastercutdata->cut.cons != NULL);
      return GCGconsGetRhs(scip, mastercutdata->cut.cons);
   case GCG_MASTERCUTTYPE_ROW:
      assert(mastercutdata->cut.row != NULL);
      return SCIProwGetRhs(mastercutdata->cut.row);
   default:
      SCIP_CALL_ABORT( SCIP_ERROR );
      return SCIP_INVALID;
   }
}

/** get the constant of the mastercut (always returns 0 if mastercut is a constraint, returns constant of row otherwise) */
SCIP_Real GCGmastercutGetConstant(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   switch( mastercutdata->type )
   {
   case GCG_MASTERCUTTYPE_CONS:
      return 0.0;
   case GCG_MASTERCUTTYPE_ROW:
      assert(mastercutdata->cut.row != NULL);
      return SCIProwGetConstant(mastercutdata->cut.row);
   default:
      SCIP_CALL_ABORT( SCIP_ERROR );
      return SCIP_INVALID;
   }
}

/** get number of nonzero entries in the mastercut */
int GCGmastercutGetNNonz(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   switch( mastercutdata->type )
   {
   case GCG_MASTERCUTTYPE_CONS:
      assert(mastercutdata->cut.cons != NULL);
      return SCIProwGetNNonz(SCIPconsGetRow(scip, mastercutdata->cut.cons));
   case GCG_MASTERCUTTYPE_ROW:
      assert(mastercutdata->cut.row != NULL);
      return SCIProwGetNNonz(mastercutdata->cut.row);
   default:
      SCIP_CALL_ABORT( SCIP_ERROR );
      return -1;
   }
}

/** get array of columns with nonzero entries */
SCIP_COL** GCGmastercutGetCols(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   switch( mastercutdata->type )
   {
   case GCG_MASTERCUTTYPE_CONS:
      assert(mastercutdata->cut.cons != NULL);
      return SCIProwGetCols(SCIPconsGetRow(scip, mastercutdata->cut.cons));
   case GCG_MASTERCUTTYPE_ROW:
      assert(mastercutdata->cut.row != NULL);
      return SCIProwGetCols(mastercutdata->cut.row);
   default:
      SCIP_CALL_ABORT( SCIP_ERROR );
      return NULL;
   }
}

/** get array of coefficients with nonzero entries */
SCIP_Real* GCGmastercutGetVals(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   switch( mastercutdata->type )
   {
   case GCG_MASTERCUTTYPE_CONS:
      assert(mastercutdata->cut.cons != NULL);
      return SCIProwGetVals(SCIPconsGetRow(scip, mastercutdata->cut.cons));
   case GCG_MASTERCUTTYPE_ROW:
      assert(mastercutdata->cut.row != NULL);
      return SCIProwGetVals(mastercutdata->cut.row);
   default:
      SCIP_CALL_ABORT( SCIP_ERROR );
      return NULL;
   }
}

/** get the additional data */
void* GCGmastercutGetData(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   )
{
   assert(mastercutdata != NULL);

   return mastercutdata->data;
}

/** calculate the coefficient of a column solution in the master cut */
SCIP_Real GCGmastercutGetCoeff(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_VAR**             solvars,            /**< array of column solution variables */
   SCIP_Real*             solvals,            /**< array of column solution values */
   int                    nsolvars,           /**< number of column solution variables and values */
   int                    probnr              /**< the pricing problem that the column belongs to */
   )
{
   SCIP_Real coef;

   assert(scip != NULL);
   assert(mastercutdata != NULL);

   SCIP_CALL( mastercutdata->mastercutGetCoeff(scip, mastercutdata, solvars, solvals, nsolvars, probnr, &coef) );

   return coef;
}
