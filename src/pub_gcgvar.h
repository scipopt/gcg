/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_gcgvar.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for GCG variables
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_PUB_GCGVAR_H__
#define __GCG_PUB_GCGVAR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns TRUE or FALSE whether variable is a pricing variable or not */
extern
SCIP_Bool GCGvarIsPricing(
   SCIP_VAR* var              /**< SCIP variable structure */
   );

/** returns TRUE or FALSE whether variable is a original variable or not */
extern
SCIP_Bool GCGvarIsOriginal(
   SCIP_VAR* var              /**< SCIP variable structure */
   );

/** returns TRUE or FALSE whether variable is a master variable or not */
extern
SCIP_Bool GCGvarIsMaster(
   SCIP_VAR* var              /**< SCIP variable structure */
   );

/** returns TRUE or FALSE whether variable is a linking variable or not */
extern
SCIP_Bool GCGvarIsLinking(
   SCIP_VAR* var              /**< SCIP variable structure */
   );

/** returns the original var of a pricing variable */
extern
SCIP_VAR* GCGpricingVarGetOriginalVar(
   SCIP_VAR* var              /**< SCIP variable structure */
   );

/** adds the original var to the pricing variable */
extern
SCIP_RETCODE GCGpricingVarAddOrigVar(
   SCIP*     scip,            /**< SCIP variable structure */
   SCIP_VAR* pricingvar,      /**< pricing variable */
   SCIP_VAR* origvar          /**< original pricing variable */
   );

/** returns the pricing var of an original variable */
extern
SCIP_VAR* GCGoriginalVarGetPricingVar(
   SCIP_VAR* var              /**< variable data structure */
   );

/** returns the pricing var of an original variable */
extern
void GCGoriginalVarSetPricingVar(
   SCIP_VAR* var,             /**< SCIP variable structure */
   SCIP_VAR* pricingvar       /**< SCIP variable structure */
   );

/** returns the pricing variables of an linking variable */
extern
SCIP_VAR** GCGlinkingVarGetPricingVars(
   SCIP_VAR* var              /**< SCIP variable structure */
   );

/** sets the pricing var of the corresponding linking variable at the specified position */
extern
void GCGlinkingVarSetPricingVar(
   SCIP_VAR* origvar,         /**< original variable */
   int       pricingprobnr,   /**< number of pricing problem */
   SCIP_VAR* var              /**< pricing variable */
   );

/** returns the number of master variables the original variable is contained in */
extern
int GCGoriginalVarGetNMastervars(
   SCIP_VAR* var              /**< variable data structure */
   );

/** returns the master variables the original variable is contained in */
extern
SCIP_VAR** GCGoriginalVarGetMastervars(
   SCIP_VAR* var              /**< variable data structure */
   );

/** returns the fraction of master variables the original variable is contained in */
extern
SCIP_Real* GCGoriginalVarGetMastervals(
   SCIP_VAR* var              /**< variable data structure */
   );

/** returns the fraction of master variables the original variable is contained in */
extern
SCIP_Real* GCGoriginalVarGetCoefs(
   SCIP_VAR* var              /**< variable data structure */
   );

/** returns the fraction of master variables the original variable is contained in */
extern
SCIP_CONS** GCGoriginalVarGetMasterconss(
   SCIP_VAR* var              /**< variable data structure */
   );

/** adds a coefficient of the master variable to the coefs array for the resp. constraint */
extern
SCIP_RETCODE GCGoriginalVarAddCoef(
   SCIP*      scip,           /**< SCIP data structure */
   SCIP_VAR*  var,            /**< variable to add coef */
   SCIP_Real  val,            /**< coefficent to set */
   SCIP_CONS* cons            /**< constraint the variable is in */
   );

/** adds variable to a new block, making a linkingvariable out of it, if necessary */
extern
SCIP_RETCODE GCGoriginalVarAddBlock(
   SCIP*     scip,            /**< SCIP data structure */
   SCIP_VAR* var,             /**< var that is added to a block */
   int       newblock         /**< the new block the variable will be in */
   );


/** returns the linking constraints */
extern
SCIP_CONS** GCGlinkingVarGetLinkingConss(
   SCIP_VAR* var              /**< variable data structure */
   );

/** sets the linking constraints*/
extern
void GCGlinkingVarSetLinkingCons(
   SCIP_VAR*  var,            /**< variable data structure */
   SCIP_CONS* cons,           /**< linking constraint */
   int        i               /**< index of pricing problem */
   );

/** returns the blocks the linking variable is in */
extern
SCIP_RETCODE GCGlinkingVarGetBlocks(
   SCIP_VAR* var,             /**< SCIP variable structure */
   int       nblocks,         /**< number of blocks the linking variable is in */
   int*      blocks           /**< array to store the blocks of the linking variable */
   );

/** returns the number of blocks the linking variable is in */
extern
int GCGlinkingVarGetNBlocks(
   SCIP_VAR* var              /**< SCIP variable structure */
   );

/** returns the number of coefficients of master constraints the original variable is contained in */
extern
int GCGoriginalVarGetNCoefs(
   SCIP_VAR* var              /**< SCIP variable structure */
   );

/** sets the number of master variables the original variable is contained in */
extern
void GCGoriginalVarSetNCoefs(
   SCIP_VAR* var,             /**< SCIP variable structure */
   int       coef             /**< number of coefficient to set */
   );

/** returns whether the master variable is a ray */
extern
SCIP_Bool GCGmasterVarIsRay(
   SCIP_VAR* var              /**< variable data structure */
   );

/** returns the number of original variables the master variable is contained in */
extern
int GCGmasterVarGetNOrigvars(
   SCIP_VAR* var              /**< variable data structure */
   );

/** returns the original variables the master variable is contained in */
extern
SCIP_VAR** GCGmasterVarGetOrigvars(
   SCIP_VAR* var              /**< variable data structure */
   );

/** returns the fraction of original variables the master variable is contained in */
extern
SCIP_Real* GCGmasterVarGetOrigvals(
   SCIP_VAR* var              /**< variable data structure */
   );

/** returns the number of original variables the pricing variable is contained in */
extern
int GCGpricingVarGetNOrigvars(
   SCIP_VAR* var              /**< variable data structure */
   );

/** returns the original variables the pricing variable is contained in */
extern
SCIP_VAR** GCGpricingVarGetOrigvars(
   SCIP_VAR* var              /**< variable data structure */
   );

/** returns the block of the variable */
extern
int GCGvarGetBlock(
   SCIP_VAR* var              /**< variable data structure */
   );

/** sets the block of the variable */
extern
void GCGvarSetBlock(
   SCIP_VAR* var,             /**< variable to set block for */
   int       block            /**< block to set */
   );

/** creates the data for all variables of the original program */
extern
SCIP_RETCODE GCGcreateOrigVarsData(
   SCIP* scip                 /**< SCIP data structure */
   );

/** creates the data for a variable of the original program */
extern
SCIP_RETCODE GCGorigVarCreateData(
   SCIP*     scip,            /**< SCIP data structure */
   SCIP_VAR* var              /**< pointer to variable object */
   );

/** returns TRUE if the linking variable is in the block, FALSE otherwise */
extern
SCIP_Bool GCGisLinkingVarInBlock(
   SCIP_VAR* var,             /**< variabel data structure */
   int       block            /**< pricing problem number */
   );

/** informs an original variable, that a variable in the master problem was created,
 * that contains a part of the original variable.
 * Saves this information in the original variable's data
 * @todo this method needs a little love
 */
extern
SCIP_RETCODE GCGoriginalVarAddMasterVar(
   SCIP*     scip,            /**< SCIP data structure */
   SCIP_VAR* origvar,         /**< original variable */
   SCIP_VAR* var,             /**< master variable */
   SCIP_Real val              /**< fraction of the original variable */
   );

/* informs an original variable, that a variable in the master problem was deleted,
 * that contains a part of the original variable.
 * Update the information in the original variable's data
 * @todo this method needs a little love
 */
extern
SCIP_RETCODE GCGoriginalVarRemoveMasterVar(
   SCIP*     scip,            /**< SCIP data structure */
   SCIP_VAR* origvar,         /**< original variable */
   SCIP_VAR* var              /**< master variable */
   );


/** creates the corresponding pricing variable for the given original variable */
extern
SCIP_RETCODE GCGoriginalVarCreatePricingVar(
   SCIP*      scip,           /**< SCIP data structure */
   SCIP_VAR*  origvar,        /**< original variable */
   SCIP_VAR** var             /**< pricing variable */
   );

/** creates the corresponding pricing variable for the given original variable */
extern
SCIP_RETCODE GCGlinkingVarCreatePricingVar(
   SCIP*       masterscip,    /**< msater problem SCIP data structure */
   SCIP*       pricingscip,   /**< pricing problem SCIP data structure */
   int         pricingprobnr, /**< number of the pricing problem */
   SCIP_VAR*   origvar,       /**< original variable */
   SCIP_VAR**  var,           /**< pointer to store new pricing variable */
   SCIP_CONS** linkcons       /**< constraint linking pricing variables */
   );

/** creates the master var and initializes the vardata */
extern
SCIP_RETCODE GCGcreateMasterVar(
   SCIP*        scip,         /**< SCIP data structure */
   SCIP*        pricingscip,  /**< pricing problem SCIP data structure */
   SCIP_VAR**   newvar,       /**< pointer to store new master variable */
   char*        varname,      /**< new variable name */
   SCIP_Real    objcoeff,     /**< new objective coeffient */
   SCIP_VARTYPE vartype,      /**< new variable type */
   SCIP_Bool    solisray,     /**< indicates whether new variable is a ray */
   int          prob,         /**< number of pricing problem that created this variable */
   int          nsolvars,     /**< number of variables in the solution */
   SCIP_Real*   solvals,      /**< values of variables in the solution */
   SCIP_VAR**   solvars       /**< variables with non zero coefficient in the solution */
   );

/** creates initial master variables and the vardata */
SCIP_RETCODE GCGcreateInitialMasterVar(
   SCIP*       scip,          /**< SCIP data structure */
   SCIP_VAR*   var,           /**< original variable */
   SCIP_VAR**  newvar         /**< pointer to store new variable */
   );


#ifdef __cplusplus
}
#endif

#endif
