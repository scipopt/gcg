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

/** Returns TRUE or FALSE whether variable is a linking variable or not */
extern
SCIP_Bool GCGvarIsLinking(
      SCIP_VAR* var /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   return vardata->blocknr == -2;
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
   assert(!GCGvarIsLinking(var));
   return vardata->data.origvardata.pricingvar;
}

/** Returns the pricing variables of an linking variable */
extern
SCIP_VAR** GCGlinkingVarGetPricingVars(
   SCIP_VAR* var /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsLinking(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   assert(vardata->data.origvardata.linkingvardata != NULL);
   assert(vardata->data.origvardata.linkingvardata->pricingvars != NULL);

   return vardata->data.origvardata.linkingvardata->pricingvars;
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

/** Returns the fraction of master variables the original variable is contained in */
extern
SCIP_Real* GCGoriginalVarGetCoefs(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.coefs != NULL || vardata->data.origvardata.ncoefs == 0 );
   return vardata->data.origvardata.coefs;
}

/** Returns the fraction of master variables the original variable is contained in */
extern
int GCGoriginalVarGetNCoefs(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.coefs != NULL || vardata->data.origvardata.ncoefs == 0 );
   return vardata->data.origvardata.ncoefs;
}

/** Returns the fraction of master variables the original variable is contained in */
extern
SCIP_CONS** GCGoriginalVarGetLinkingCons(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.linkconss != NULL);
   return vardata->data.origvardata.linkconss;
}

/** Returns the fraction of master variables the original variable is contained in */
extern
SCIP_CONS** GCGlinkingVarGetLinkingConss(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));
   assert(GCGvarIsLinking(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.linkingvardata != NULL);
   assert(vardata->data.origvardata.linkingvardata->linkconss != NULL);
   return vardata->data.origvardata.linkingvardata->linkconss;
}

/** Returns the number of original variables the master variable is contained in */
extern
int GCGmasterVarGetNOrigvars(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsMaster(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.mastervardata.norigvars >= 0);
   assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
   assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);
   assert(vardata->blocknr != -1 || vardata->data.mastervardata.norigvars == 1 );

   return vardata->data.mastervardata.norigvars;
}

/** Returns the original variables the master variable is contained in */
extern
SCIP_VAR** GCGmasterVarGetOrigvars(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsMaster(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
   assert(vardata->blocknr != -1 || vardata->data.mastervardata.origvars != NULL);
   assert(vardata->blocknr != -1 || vardata->data.mastervardata.origvars[0] != NULL);
   assert(vardata->blocknr != -1 || GCGvarGetBlock(vardata->data.mastervardata.origvars[0]) == -1 
      || GCGvarIsLinking(vardata->data.mastervardata.origvars[0]));


   return vardata->data.mastervardata.origvars;
}

/** Returns the fraction of original variables the master variable is contained in */
extern
SCIP_Real* GCGmasterVarGetOrigvals(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsMaster(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);
   return vardata->data.mastervardata.origvals;
}

/** Returns the number of original variables the pricing variable is contained in */
extern
int GCGpricingVarGetNOrigvars(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsPricing(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.pricingvardata.norigvars >= 0);
   return vardata->data.pricingvardata.norigvars;
}

/** Returns the original variables the pricing variable is contained in */
extern
SCIP_VAR** GCGpricingVarGetOrigvars(
      SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsPricing(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.pricingvardata.origvars != NULL);
   return vardata->data.pricingvardata.origvars;
}

/** Returns the block of the variable */
extern
int GCGvarGetBlock(
   SCIP_VAR* var
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->blocknr >= -2);
   return vardata->blocknr;
 }

/** Returns TRUE if the linking variable is in the block, FALSE otherwise **/
extern
SCIP_Bool GCGisLinkingVarInBlock(
   SCIP_VAR* var,
   int block
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(block >= 0);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(GCGvarIsLinking(var));
   assert(GCGvarIsOriginal(var));

   assert(vardata->data.origvardata.linkingvardata != NULL);
   assert(vardata->data.origvardata.linkingvardata->pricingvars != NULL);

   return vardata->data.origvardata.linkingvardata->pricingvars[block] != NULL;

}

/* informs an original variable, that a variable in the master problem was created, 
 * that contains a part of the original variable.
 * Saves this information in the original variable's data 
 * @todo this method needs a little love
 */
SCIP_RETCODE GCGoriginalVarAddMasterVar(
   SCIP*                 scip,                  /**< SCIP data structure                */
   SCIP_VAR*             origvar,               /**< Original variable                  */
   SCIP_VAR*             var,                   /**< Master variable                    */
   SCIP_Real             val                    /**< Fraction of the original variable  */
   )
{
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(origvar != NULL);
   assert(var != NULL);

   vardata = SCIPvarGetData(origvar);

   assert(vardata != NULL);
   assert(GCGvarIsOriginal(origvar));
   assert(vardata->data.origvardata.mastervars != NULL);
   assert(vardata->data.origvardata.mastervals != NULL);
   assert(vardata->data.origvardata.nmastervars >= 0);
   assert(vardata->data.origvardata.maxmastervars >= vardata->data.origvardata.nmastervars);

   /* realloc mastervars array of the original variable, if needed */
   if( vardata->data.origvardata.maxmastervars == vardata->data.origvardata.nmastervars )
   {
      SCIP_CALL( SCIPreallocMemoryArray(pricerdata->origprob, &(vardata->data.origvardata.mastervars),
            2*vardata->data.origvardata.maxmastervars) );
      SCIP_CALL( SCIPreallocMemoryArray(pricerdata->origprob, &(vardata->data.origvardata.mastervals),
            2*vardata->data.origvardata.maxmastervars) );
      SCIPdebugMessage("mastervars array of var %s resized from %d to %d\n", SCIPvarGetName(origvar), 
         vardata->data.origvardata.maxmastervars, 2*vardata->data.origvardata.maxmastervars);
      vardata->data.origvardata.maxmastervars = 2*vardata->data.origvardata.maxmastervars;
   }
   /* add information to the original variable's vardata */
   vardata->data.origvardata.mastervars[vardata->data.origvardata.nmastervars] = var;
   vardata->data.origvardata.mastervals[vardata->data.origvardata.nmastervars] = val;
   vardata->data.origvardata.nmastervars++;

   return SCIP_OKAY;
}
