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

extern
SCIP_Bool isVarRelevant(
   SCIP_VAR* var
   );

extern
SCIP_VAR* getRelevantVariable(
   SCIP_VAR *var
   );

consType SCIPconsGetType( SCIP_CONS *cons );

SCIP_Real SCIPgetRhsXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_Real SCIPgetLhsXXX(SCIP * scip, SCIP_CONS * cons);

int SCIPgetNVarsXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_Real SCIPgetDualsolXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_VAR ** SCIPgetVarsXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_Real * SCIPgetValsXXX(SCIP * scip, SCIP_CONS * cons);

#endif /* SCIP_MISC_H_ */
