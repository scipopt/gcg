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

#if 0
SCIP_RETCODE SCIPcopyProblem(
      SCIP *scip,
      SCIP* const* target,
      SCIP_HASHMAP ** varmap,
      SCIP_CONS * const* conss,
      int nconss,
      SCIP_VAR ** vars,
      int nvars);

SCIP_RETCODE SCIPcopyPresolvedProblem(SCIP *scip,
SCIP* const* target,
SCIP_CONS* const* conss,
int nconss,
SCIP_VAR* const* vars,
int nvars);
#endif

consType SCIPconsGetType(SCIP* scip, SCIP_CONS *cons);

SCIP_Real SCIPgetRhsXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_Real SCIPgetLhsXXX(SCIP * scip, SCIP_CONS * cons);

int SCIPgetNVarsXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_Real SCIPgetDualsolXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_VAR ** SCIPgetVarsXXX(SCIP * scip, SCIP_CONS * cons);

SCIP_Real * SCIPgetValsXXX(SCIP * scip, SCIP_CONS * cons);

#endif /* SCIP_MISC_H_ */
