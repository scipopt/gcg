/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file    extendedmasterconsdata.c
 * @brief   methods for interacting with GCG_EXTENDEDMASTERCONSDATA
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/def.h"
#include "gcg/pub_extendedmasterconsdata.h"
#include "gcg/extendedmasterconsdata.h"
#include "gcg/gcg.h"
#include "gcg/scip_misc.h"
#include "gcg/pricer_gcg.h"
#include "gcg/struct_extendedmasterconsdata.h"

#include "scip/scip.h"
#include "scip/cons_linear.h"

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
SCIP_RETCODE GCGpricingmodificationCreate(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_PRICINGMODIFICATION* pricingmodification, /**< pointer to store the pricing modification */
   int                    blocknr,             /**< block number of the extended master cons */
   SCIP_VAR*              coefvar,             /**< variable in the pricing problem inferred from the extended master cons
                                                  * always has the objective coefficient of the negated dual value of the extended master cons
                                                  * its solution value corresponds to the coefficient of the new mastervariable in the extended master cons */
   SCIP_VAR**             additionalvars,      /**< array of additional variables with no objective coefficient in the pricing programs inferred from the extended master cons */
   int                    nadditionalvars,     /**< number of additional variables in the pricing programs */
   SCIP_CONS**            additionalconss,     /**< array of additional constraints in the pricing programs inferred from the extended master cons */
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
      assert(SCIPisZero(scip, SCIPvarGetObj(additionalvars[i])));
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

/** create an extended master cons, taking ownership over pricingmodifications */
SCIP_RETCODE GCGextendedmasterconsCreateFromCons(
   SCIP*                  scip,                 /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA**    extendedmasterconsdata,        /**< pointer to store the extended master cons data */
   SCIP_CONS*             cons,                 /**< constraint in the master problem that represents the extended master cons */
   GCG_PRICINGMODIFICATION* pricingmodifications,/**< pricing modifications for the extended master cons */
   int                    npricingmodifications,/**< number of pricing modifications for the extended master cons */
   void*                  data,                 /**< any data that might be required to calculate the coefficient of a column solution */
   GCG_DECL_EXTENDEDMASTERCONSGETCOEFF((*extendedmasterconsGetCoeff))/**< callback to calculate the coefficient of a column solution */
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
   assert(extendedmasterconsdata != NULL);
   assert(*extendedmasterconsdata == NULL);
   assert(cons != NULL);
   assert(pricingmodifications != NULL || npricingmodifications == 0);

#ifndef NDEBUG
   originalproblem = GCGgetOriginalprob(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &seenblocks, GCGgetNPricingprobs(originalproblem)) );
   for( i = 0; i < GCGgetNPricingprobs(originalproblem); i++ )
      seenblocks[i] = FALSE;
#endif

   for( i = 0; i < npricingmodifications; i++ )
   {
      assert(pricingmodifications[i].blocknr >= 0);
      assert(pricingmodifications[i].blocknr < GCGgetNPricingprobs(originalproblem));
      assert(GCGisPricingprobRelevant(originalproblem, pricingmodifications[i].blocknr));
#ifndef NDEBUG
      assert(!seenblocks[pricingmodifications[i].blocknr]);
      seenblocks[pricingmodifications[i].blocknr] = TRUE;
#endif
   }

#ifndef NDEBUG
   SCIPfreeBlockMemoryArray(scip, &seenblocks, GCGgetNPricingprobs(originalproblem));
   assert(seenblocks == NULL);
#endif

   SCIP_CALL( SCIPallocBlockMemory(scip, extendedmasterconsdata) );

   (*extendedmasterconsdata)->type = GCG_EXTENDEDMASTERCONSTYPE_CONS;
   (*extendedmasterconsdata)->cons.cons = cons;
   (*extendedmasterconsdata)->pricingmodifications = pricingmodifications;
   (*extendedmasterconsdata)->npricingmodifications = npricingmodifications;
   (*extendedmasterconsdata)->data = data;
   (*extendedmasterconsdata)->extendedmasterconsGetCoeff = extendedmasterconsGetCoeff;

   for( i = 0; i < npricingmodifications; i++ ) {
      SCIPvarGetData(pricingmodifications[i].coefvar)->data.inferredpricingvardata.extendedmasterconsdata = *extendedmasterconsdata;
      for( j = 0; j < pricingmodifications[i].nadditionalvars; j++ ) {
         SCIPvarGetData(pricingmodifications[i].additionalvars[j])->data.inferredpricingvardata.extendedmasterconsdata = *extendedmasterconsdata;
      }
   }

   return SCIP_OKAY;
}

/** create an extended master cons, taking ownership over pricingmodifications */
SCIP_RETCODE GCGextendedmasterconsCreateFromRow(
   SCIP*                  scip,                 /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA**    extendedmasterconsdata,        /**< pointer to store the extended master cons data */
   SCIP_ROW*              row,                  /**< row in the master problem that represents the extended master cons cut */
   GCG_PRICINGMODIFICATION* pricingmodifications,/**< pricing modifications for the extended master cons */
   int                    npricingmodifications,/**< number of pricing modifications for the extended master cons */
   void*                  data,                 /**< any data that might be required to calculate the coefficient of a column solution */
   GCG_DECL_EXTENDEDMASTERCONSGETCOEFF((*extendedmasterconsGetCoeff))/**< callback to calculate the coefficient of a column solution */
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
   assert(extendedmasterconsdata != NULL);
   assert(*extendedmasterconsdata == NULL);
   assert(row != NULL);
   assert(pricingmodifications != NULL || npricingmodifications == 0);

#ifndef NDEBUG
   originalproblem = GCGgetOriginalprob(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &seenblocks, GCGgetNPricingprobs(originalproblem)) );
   for( i = 0; i < GCGgetNPricingprobs(originalproblem); i++ )
      seenblocks[i] = FALSE;
#endif

   for( i = 0; i < npricingmodifications; i++ )
   {
      assert(pricingmodifications[i].blocknr >= 0);
      assert(pricingmodifications[i].blocknr < GCGgetNPricingprobs(originalproblem));
      assert(GCGisPricingprobRelevant(originalproblem, pricingmodifications[i].blocknr));
#ifndef NDEBUG
      assert(!seenblocks[pricingmodifications[i].blocknr]);
      seenblocks[pricingmodifications[i].blocknr] = TRUE;
#endif
   }

#ifndef NDEBUG
   SCIPfreeBlockMemoryArray(scip, &seenblocks, GCGgetNPricingprobs(originalproblem));
#endif

   SCIP_CALL( SCIPallocBlockMemory(scip, extendedmasterconsdata) );

   (*extendedmasterconsdata)->type = GCG_EXTENDEDMASTERCONSTYPE_ROW;
   (*extendedmasterconsdata)->cons.row = row;
   (*extendedmasterconsdata)->pricingmodifications = pricingmodifications;
   (*extendedmasterconsdata)->npricingmodifications = npricingmodifications;
   (*extendedmasterconsdata)->data = data;
   (*extendedmasterconsdata)->extendedmasterconsGetCoeff = extendedmasterconsGetCoeff;

   for( i = 0; i < npricingmodifications; i++ ) {
      SCIPvarGetData(pricingmodifications[i].coefvar)->data.inferredpricingvardata.extendedmasterconsdata = *extendedmasterconsdata;
      for( j = 0; j < pricingmodifications[i].nadditionalvars; j++ ) {
         SCIPvarGetData(pricingmodifications[i].additionalvars[j])->data.inferredpricingvardata.extendedmasterconsdata = *extendedmasterconsdata;
      }
   }

   return SCIP_OKAY;
}

/** free an extended master cons */
SCIP_RETCODE GCGextendedmasterconsFree(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA**    extendedmasterconsdata        /**< pointer to the extended master cons data */
   )
{
   int i;

   assert(scip != NULL);
   assert(GCGisMaster(scip));
   assert(extendedmasterconsdata != NULL);
   assert(*extendedmasterconsdata != NULL);

   switch( (*extendedmasterconsdata)->type )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_CONS:
      assert((*extendedmasterconsdata)->cons.cons != NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &(*extendedmasterconsdata)->cons.cons) );
      break;
   case GCG_EXTENDEDMASTERCONSTYPE_ROW:
      assert((*extendedmasterconsdata)->cons.row != NULL);
      SCIP_CALL( SCIPreleaseRow(scip, &(*extendedmasterconsdata)->cons.row) );
      break;
   default:
      SCIP_CALL( SCIP_INVALIDDATA );
   }

   for( i = 0; i < (*extendedmasterconsdata)->npricingmodifications; i++ )
   {
      SCIP_CALL( GCGpricingmodificationFree(scip, &(*extendedmasterconsdata)->pricingmodifications[i]) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*extendedmasterconsdata)->pricingmodifications, (*extendedmasterconsdata)->npricingmodifications);
   assert((*extendedmasterconsdata)->pricingmodifications == NULL);
   SCIPfreeBlockMemory(scip, extendedmasterconsdata);
   assert(*extendedmasterconsdata == NULL);

   return SCIP_OKAY;
}

/** determine whether the extendedmasterconsdata is active in the masterscip */
SCIP_Bool GCGextendedmasterconsIsActive(
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   switch( extendedmasterconsdata->type )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_CONS:
      assert(extendedmasterconsdata->cons.cons != NULL);
      return SCIPconsIsActive(extendedmasterconsdata->cons.cons);
   case GCG_EXTENDEDMASTERCONSTYPE_ROW:
      assert(extendedmasterconsdata->cons.row != NULL);
      return SCIProwIsInLP(extendedmasterconsdata->cons.row);
   default:
      SCIP_CALL_ABORT( SCIP_INVALIDDATA );
      return FALSE;
   }
}

/** add a new variable along with its coefficient to the extended master cons */
SCIP_RETCODE GCGextendedmasterconsAddMasterVar(
   SCIP*                  masterscip,         /**< master scip */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata,      /**< extended master cons data */
   SCIP_VAR*              var,                /**< variable to add */
   SCIP_Real              coef                /**< coefficient of the variable */
   )
{
   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(extendedmasterconsdata != NULL);
   assert(var != NULL);

   switch( extendedmasterconsdata->type )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_CONS:
      assert(extendedmasterconsdata->cons.cons != NULL);
      SCIP_CALL( SCIPaddCoefLinear(masterscip, extendedmasterconsdata->cons.cons, var, coef) );
      break;
   case GCG_EXTENDEDMASTERCONSTYPE_ROW:
      assert(extendedmasterconsdata->cons.row != NULL);
      SCIP_CALL( SCIPaddVarToRow(masterscip, extendedmasterconsdata->cons.row, var, coef) );
      break;
   default:
      SCIP_CALL( SCIP_INVALIDDATA );
   }

   return SCIP_OKAY;
}

/** update the extended master cons with the new dual value */
SCIP_RETCODE GCGextendedmasterconsUpdateDualValue(
   SCIP*                  masterscip,         /**< master scip */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata,      /**< extended master cons data */
   SCIP_Real              dualvalue           /**< dual value */
   )
{
   int i;
   SCIP* origscip;
   SCIP* pricingscip;

   assert(extendedmasterconsdata != NULL);
   assert(GCGisMaster(masterscip));

   origscip = GCGmasterGetOrigprob(masterscip);
   assert(origscip != NULL);

   for( i = 0; i < extendedmasterconsdata->npricingmodifications; i++ )
   {
      assert(extendedmasterconsdata->pricingmodifications[i].coefvar != NULL);
      assert(GCGvarIsInferredPricing(extendedmasterconsdata->pricingmodifications[i].coefvar));

      pricingscip = GCGgetPricingprob(origscip, extendedmasterconsdata->pricingmodifications[i].blocknr);
      assert(pricingscip != NULL);

      SCIP_CALL( SCIPchgVarObj(pricingscip, extendedmasterconsdata->pricingmodifications[i].coefvar, -dualvalue) );
   }

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** get the constraint that is the extended master cons */
SCIP_CONS* GCGextendedmasterconsGetCons(
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);
   assert(extendedmasterconsdata->type == GCG_EXTENDEDMASTERCONSTYPE_CONS);
   assert(extendedmasterconsdata->cons.cons != NULL);

   return extendedmasterconsdata->cons.cons;
}
#endif

#ifndef NDEBUG
/** get the row that is the extended master cons */
SCIP_ROW* GCGextendedmasterconsGetRow(
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);
   assert(extendedmasterconsdata->type == GCG_EXTENDEDMASTERCONSTYPE_ROW);
   assert(extendedmasterconsdata->cons.row != NULL);

   return extendedmasterconsdata->cons.row;
}
#endif

#ifndef NDEBUG
/** get the variable that determines the coefficient of a column in the extended master cons */
SCIP_VAR* GCGpricingmodificationGetCoefVar(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   )
{
   assert(pricingmodification != NULL);

   return pricingmodification->coefvar;
}
#endif

#ifndef NDEBUG
/** get the additional variables that are inferred by the extended master cons */
SCIP_VAR** GCGpricingmodificationGetAdditionalVars(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   )
{
   assert(pricingmodification != NULL);

   return pricingmodification->additionalvars;
}
#endif

#ifndef NDEBUG
/** get the number of additional variables that are inferred by the extended master cons */
int GCGpricingmodificationGetNAdditionalVars(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   )
{
   assert(pricingmodification != NULL);

   return pricingmodification->nadditionalvars;
}
#endif

#ifndef NDEBUG
/** get the additional constraints that are inferred by the extended master cons */
SCIP_CONS** GCGpricingmodificationGetAdditionalConss(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   )
{
   assert(pricingmodification != NULL);

   return pricingmodification->additionalconss;
}
#endif

#ifndef NDEBUG
/** get the number of additional constraints that are inferred by the extended master cons */
int GCGpricingmodificationGetNAdditionalConss(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   )
{
   assert(pricingmodification != NULL);

   return pricingmodification->nadditionalconss;
}
#endif

/** get the pricing modification for a block, if exists, else NULL */
GCG_PRICINGMODIFICATION* GCGextendedmasterconsGetPricingModification(
   SCIP*                  masterscip,         /**< master scip */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata,      /**< extended master cons data */
   int                    blocknr             /**< block number */
   )
{
   int i;

   assert(extendedmasterconsdata != NULL);
   assert(GCGisMaster(masterscip));
   assert(blocknr >= 0);

   assert(blocknr < GCGgetNPricingprobs(GCGgetOriginalprob(masterscip)));

   for( i = 0; i < extendedmasterconsdata->npricingmodifications; i++ )
   {
      if( extendedmasterconsdata->pricingmodifications[i].blocknr == blocknr )
         return &extendedmasterconsdata->pricingmodifications[i];
   }

   return NULL;
}

#ifndef NDEBUG
/** get the pricing modifications for the extended master cons */
GCG_PRICINGMODIFICATION* GCGextendedmasterconsGetPricingModifications(
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   return extendedmasterconsdata->pricingmodifications;
}
#endif

#ifndef NDEBUG
/** get the number of pricing modifications for the extended master cons */
int GCGextendedmasterconsGetNPricingModifications(
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   return extendedmasterconsdata->npricingmodifications;
}
#endif

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
SCIP_RETCODE GCGextendedmasterconsApplyPricingModifications(
   SCIP*                  masterscip,         /**< master scip */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   int i;
   SCIP* origscip;
   SCIP* pricingprob;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(extendedmasterconsdata != NULL);

   origscip = GCGmasterGetOrigprob(masterscip);
   assert(origscip != NULL);

   for( i = 0; i < extendedmasterconsdata->npricingmodifications; i++ )
   {
      pricingprob = GCGgetPricingprob(origscip, extendedmasterconsdata->pricingmodifications[i].blocknr);
      assert(pricingprob != NULL);
      SCIP_CALL( GCGpricingmodificationApply(pricingprob, extendedmasterconsdata->pricingmodifications[i]) );
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
SCIP_RETCODE GCGextendedmasterconsUndoPricingModifications(
   SCIP*                  masterscip,         /**< master scip */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   int i;
   SCIP* origscip;
   SCIP* pricingprob;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(extendedmasterconsdata != NULL);

   origscip = GCGmasterGetOrigprob(masterscip);
   assert(origscip != NULL);

   for( i = 0; i < extendedmasterconsdata->npricingmodifications; i++ )
   {
      pricingprob = GCGgetPricingprob(origscip, extendedmasterconsdata->pricingmodifications[i].blocknr);
      assert(pricingprob != NULL);
      SCIP_CALL( GCGpricingmodificationUndo(pricingprob, extendedmasterconsdata->pricingmodifications[i]) );
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

   return pricingmodification.coefvar == var;
}

/** check whether a given variable is a coefficient variable of a given extended master cons */
SCIP_Bool GCGextendedmasterconsIsCoefVar(
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata,      /**< extended master cons data */
   SCIP_VAR*              var                 /**< variable to check */
   )
{
   int i;

   assert(extendedmasterconsdata != NULL);
   assert(var != NULL);

   for( i = 0; i < extendedmasterconsdata->npricingmodifications; i++ )
   {
      if( GCGpricingmodificationIsCoefVar(extendedmasterconsdata->pricingmodifications[i], var) )
         return TRUE;
   }

   return FALSE;
}

/** get name of the extended master cons */
const char* GCGextendedmasterconsGetName(
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   switch( extendedmasterconsdata->type )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_CONS:
      assert(extendedmasterconsdata->cons.cons != NULL);
      return SCIPconsGetName(extendedmasterconsdata->cons.cons);
   case GCG_EXTENDEDMASTERCONSTYPE_ROW:
      assert(extendedmasterconsdata->cons.row != NULL);
      return SCIProwGetName(extendedmasterconsdata->cons.row);
   default:
      SCIP_CALL_ABORT( SCIP_INVALIDDATA );
      return NULL;
   }
}

/** get the lhs of the extended master cons */
SCIP_Real GCGextendedmasterconsGetLhs(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   switch( extendedmasterconsdata->type )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_CONS:
      assert(extendedmasterconsdata->cons.cons != NULL);
      return GCGconsGetLhs(scip, extendedmasterconsdata->cons.cons);
   case GCG_EXTENDEDMASTERCONSTYPE_ROW:
      assert(extendedmasterconsdata->cons.row != NULL);
      return SCIProwGetLhs(extendedmasterconsdata->cons.row);
   default:
      SCIP_CALL_ABORT( SCIP_INVALIDDATA );
      return SCIP_INVALID;
   }
}

/** get the rhs of the extended master cons */
SCIP_Real GCGextendedmasterconsGetRhs(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   switch( extendedmasterconsdata->type )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_CONS:
      assert(extendedmasterconsdata->cons.cons != NULL);
      return GCGconsGetRhs(scip, extendedmasterconsdata->cons.cons);
   case GCG_EXTENDEDMASTERCONSTYPE_ROW:
      assert(extendedmasterconsdata->cons.row != NULL);
      return SCIProwGetRhs(extendedmasterconsdata->cons.row);
   default:
      SCIP_CALL_ABORT( SCIP_INVALIDDATA );
      return SCIP_INVALID;
   }
}

/** get the constant of the extended master cons (always returns 0 if extended master cons is a constraint, returns constant of row otherwise) */
SCIP_Real GCGextendedmasterconsGetConstant(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   switch( extendedmasterconsdata->type )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_CONS:
      return 0.0;
   case GCG_EXTENDEDMASTERCONSTYPE_ROW:
      assert(extendedmasterconsdata->cons.row != NULL);
      return SCIProwGetConstant(extendedmasterconsdata->cons.row);
   default:
      SCIP_CALL_ABORT( SCIP_INVALIDDATA );
      return SCIP_INVALID;
   }
}

/** get number of nonzero entries in the extended master cons */
int GCGextendedmasterconsGetNNonz(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   switch( extendedmasterconsdata->type )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_CONS:
      assert(extendedmasterconsdata->cons.cons != NULL);
      return SCIProwGetNNonz(SCIPconsGetRow(scip, extendedmasterconsdata->cons.cons));
   case GCG_EXTENDEDMASTERCONSTYPE_ROW:
      assert(extendedmasterconsdata->cons.row != NULL);
      return SCIProwGetNNonz(extendedmasterconsdata->cons.row);
   default:
      SCIP_CALL_ABORT( SCIP_INVALIDDATA );
      return -1;
   }
}

/** get array of columns with nonzero entries */
SCIP_COL** GCGextendedmasterconsGetCols(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   switch( extendedmasterconsdata->type )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_CONS:
      assert(extendedmasterconsdata->cons.cons != NULL);
      return SCIProwGetCols(SCIPconsGetRow(scip, extendedmasterconsdata->cons.cons));
   case GCG_EXTENDEDMASTERCONSTYPE_ROW:
      assert(extendedmasterconsdata->cons.row != NULL);
      return SCIProwGetCols(extendedmasterconsdata->cons.row);
   default:
      SCIP_CALL_ABORT( SCIP_INVALIDDATA );
      return NULL;
   }
}

/** get array of coefficients with nonzero entries */
SCIP_Real* GCGextendedmasterconsGetVals(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   switch( extendedmasterconsdata->type )
   {
   case GCG_EXTENDEDMASTERCONSTYPE_CONS:
      assert(extendedmasterconsdata->cons.cons != NULL);
      return SCIProwGetVals(SCIPconsGetRow(scip, extendedmasterconsdata->cons.cons));
   case GCG_EXTENDEDMASTERCONSTYPE_ROW:
      assert(extendedmasterconsdata->cons.row != NULL);
      return SCIProwGetVals(extendedmasterconsdata->cons.row);
   default:
      SCIP_CALL_ABORT( SCIP_INVALIDDATA );
      return NULL;
   }
}

#ifndef NDEBUG
/** get the additional data */
void* GCGextendedmasterconsGetData(
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   return extendedmasterconsdata->data;
}
#endif

/** calculate the coefficient of a column solution in the extended master cons */
SCIP_Real GCGextendedmasterconsGetCoeff(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata,      /**< extended master cons data */
   SCIP_VAR**             solvars,            /**< array of column solution variables */
   SCIP_Real*             solvals,            /**< array of column solution values */
   int                    nsolvars,           /**< number of column solution variables and values */
   int                    probnr              /**< the pricing problem that the column belongs to */
   )
{
   SCIP_Real coef;

   assert(scip != NULL);
   assert(extendedmasterconsdata != NULL);

   SCIP_CALL( extendedmasterconsdata->extendedmasterconsGetCoeff(scip, extendedmasterconsdata, solvars, solvals, nsolvars, probnr, &coef) );

   return coef;
}

#ifndef NDEBUG
/** gets the type of the extended master cons */
GCG_EXTENDEDMASTERCONSTYPE GCGextendedmasterconsGetType(
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   )
{
   assert(extendedmasterconsdata != NULL);

   return extendedmasterconsdata->type;
}
#endif
