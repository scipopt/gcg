/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   gcgvar.c
 * @brief  GCG variable access functions
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "pub_gcgvar.h"
#include "struct_vardata.h"

/** Returns TRUE or FALSE whether variable is a pricing variable or not */
extern
SCIP_Bool GCGvarIsPricing(
   SCIP_VAR* var /**< SCIP variable structure */
)
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   return vardata->vartype == GCG_VARTYPE_PRICING;
}

/** Returns TRUE or FALSE whether variable is a master variable or not */
extern
SCIP_Bool GCGvarIsMaster(
      SCIP_VAR* var /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   return vardata->vartype == GCG_VARTYPE_MASTER;
}

/** Returns TRUE or FALSE whether variable is a original variable or not */
extern
SCIP_Bool GCGvarIsOriginal(
      SCIP_VAR* var /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   return vardata->vartype == GCG_VARTYPE_ORIGINAL;
}

/** Returns the pricing var of an original variable */
extern
SCIP_VAR* GCGoriginalVarGetPricingVar(
      SCIP_VAR* var /**< SCIP variable structure */
)
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   assert(vardata->data.origvardata.pricingvar != NULL);
   assert(vardata->data.origvardata.linkingvardata == NULL);

   return vardata->data.origvardata.pricingvar;
}

/** Returns the original var of a pricing variable */
extern
SCIP_VAR* GCGpricingVarGetOriginalVar(
      SCIP_VAR* var /**< SCIP variable structure */
)
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsPricing(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   assert(vardata->data.pricingvardata.norigvars >= 0);
   assert(vardata->data.pricingvardata.origvars != NULL);
   assert(vardata->data.pricingvardata.origvars[0] != NULL);
   assert(vardata->blocknr >= 0); /* variable belongs to exactly one block */

   return vardata->data.pricingvardata.origvars[0];
}

/** Returns the number of master variables the original variable is contained in */
extern
int GCGoriginalVarGetNMastervars(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.nmastervars >= 0);
   return vardata->data.origvardata.nmastervars;
}

/** Returns the master variables the original variable is contained in */
extern
SCIP_VAR** GCGoriginalVarGetMastervars(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.mastervars != NULL);
   return vardata->data.origvardata.mastervars;
}

/** Returns the fraction of master variables the original variable is contained in */
extern
SCIP_Real* GCGoriginalVarGetMastervals(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.mastervals != NULL);
   return vardata->data.origvardata.mastervals;
}
