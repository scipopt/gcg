/*
 * scip_misc.h
 *
 *  Created on: Apr 8, 2010
 *      Author: mbergner
 */

#ifndef SCIP_MISC_H_
#define SCIP_MISC_H_

#include "scip/scip.h"

typedef enum  {
   linear, knapsack, varbound, setpacking, setcovering, setpartitioning, logicor, sos1, sos2, unknown, nconsTypeItems
} consType;

/** returns TRUE if variable is relevant, FALSE otherwise */
extern
SCIP_Bool SCIPisVarRelevant(
   SCIP_VAR* var              /**< variable to test */
   );

/** returns the relevant variable, if possible */
extern
SCIP_VAR* SCIPgetRelevantVariable(
   SCIP_VAR* var              /**< variable to test */
   );

consType SCIPconsGetType( SCIP_CONS *cons );

SCIP_Real SCIPgetRhsXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_Real SCIPgetLhsXXX(SCIP * scip, SCIP_CONS * cons);

int SCIPgetNVarsXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_Real SCIPgetDualsolXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_VAR ** SCIPgetVarsXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_Real * SCIPgetValsXXX(SCIP * scip, SCIP_CONS * cons);

#endif /* SCIP_MISC_H_ */
