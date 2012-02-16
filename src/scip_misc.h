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

/** returns the type of an arbitrary SCIP constraint */
extern
consType SCIPconsGetType(
   SCIP_CONS* cons            /**< constraint to get type for */
   );

/** returns the rhs of an arbitrary SCIP constraint */
extern
SCIP_Real SCIPgetRhsXXX(
   SCIP*      scip,           /**< SCIP data structure */
   SCIP_CONS* cons            /**< constraint to get left hand side for */
   );

/** Returns the lhs of an arbitrary SCIP constraint */
extern
SCIP_Real SCIPgetLhsXXX(
   SCIP*      scip,           /**< SCIP data structure */
   SCIP_CONS* cons            /**< constraint to get left hand side for */
   );

/** Returns the number of variables in an arbitrary SCIP constraint */
extern
int SCIPgetNVarsXXX(
   SCIP*      scip,           /**< SCIP data structure */
   SCIP_CONS* cons            /**< constraint to get number of variables */
   );

/** Returns the variable array of an arbitrary SCIP constraint */
SCIP_RETCODE SCIPgetVarsXXX(
   SCIP*      scip,           /**< SCIP data structure */
   SCIP_CONS* cons,           /**< constraint to get variables from */
   SCIP_VAR** vars,           /**< array where variables are stored */
   int        nvars           /**< size of storage array */
   );

/** Returns the dual solution value of an arbitrary SCIP constraint */
SCIP_Real SCIPgetDualsolXXX(
   SCIP*      scip,           /**< SCIP data structure */
   SCIP_CONS* cons            /**< constraint to get dual solution */
   );


/**
 * Returns the value array of an arbitrary SCIP constraint
 * @todo SOS1 & SOS2 not implemented yet
 */
SCIP_RETCODE SCIPgetValsXXX(
   SCIP*      scip,           /**< SCIP data structure */
   SCIP_CONS* cons,           /**< constraint to get values from */
   SCIP_Real* vals,           /**< array where values are stored */
   int        nvals           /**< size of storage array */
   );
#endif /* SCIP_MISC_H_ */
