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

/**@file    scip_misc.h
 * @brief   various SCIP helper methods
 * @author  Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_SCIP_MISC_H__
#define GCG_SCIP_MISC_H__

#include "scip/scip.h"
#include "scip/cons_setppc.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif


/** constraint types */
typedef enum  {
   linear, knapsack, varbound, setpacking, setcovering, setpartitioning,
   logicor, sos1, sos2, unknown, nconsTypeItems, indicator
} consType;

/**@defgroup MISC Miscellaneous
* @ingroup PUBLICCOREAPI
* @{
  */

/** returns TRUE if variable is relevant, FALSE otherwise */
GCG_EXPORT
SCIP_Bool GCGisVarRelevant(
   SCIP_VAR*             var                 /**< variable to test */
   );

/** returns the type of an arbitrary SCIP constraint */
GCG_EXPORT
consType GCGconsGetType(
   SCIP*                 scip,               /**< SCIP data structure the cons belongs to */
   SCIP_CONS*            cons                /**< constraint to get type for */
   );

/** returns the rhs of an arbitrary SCIP constraint */
GCG_EXPORT
SCIP_Real GCGconsGetRhs(
   SCIP*                 scip,               /**< SCIP data structure the cons belongs to */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   );

/** returns the lhs of an arbitrary SCIP constraint */
GCG_EXPORT
SCIP_Real GCGconsGetLhs(
   SCIP*                 scip,               /**< SCIP data structure the cons belongs to */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   );

/** returns the dual farkas sol of an arbitrary SCIP constraint */
GCG_EXPORT
SCIP_Real GCGconsGetDualfarkas(
   SCIP*                 scip,               /**< SCIP data structure the cons belongs to */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   );

/** returns the dual sol of an arbitrary SCIP constraint */
GCG_EXPORT
SCIP_Real GCGconsGetDualsol(
   SCIP*                 scip,               /**< SCIP data structure the cons belongs to */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   );

/** returns the number of variables in an arbitrary SCIP constraint */
GCG_EXPORT
int GCGconsGetNVars(
   SCIP*                 scip,               /**< SCIP data structure the cons belongs to */
   SCIP_CONS*            cons                /**< constraint to get number of variables */
   );

/** returns the variable array of an arbitrary SCIP constraint */
GCG_EXPORT
SCIP_RETCODE GCGconsGetVars(
   SCIP*                 scip,               /**< SCIP data structure the cons belongs to */
   SCIP_CONS*            cons,               /**< constraint to get variables from */
   SCIP_VAR**            vars,               /**< array where variables are stored */
   int                   nvars               /**< size of storage array */
   );

/** returns the value array of an arbitrary SCIP constraint */
GCG_EXPORT
SCIP_RETCODE GCGconsGetVals(
   SCIP*                 scip,               /**< SCIP data structure the cons belongs to */
   SCIP_CONS*            cons,               /**< constraint to get values from */
   SCIP_Real*            vals,               /**< array where values are stored */
   int                   nvals               /**< size of storage array */
   );

/** returns true if the constraint should be a master constraint and false otherwise */
SCIP_Bool GCGconsIsRanged(
   SCIP*                 scip,               /**< SCIP data structure the cons belongs to */
   SCIP_CONS*            cons                /**< constraint to check */
);

/** returns true if the constraint should be a master constraint and false otherwise */
SCIP_Bool GCGgetConsIsSetppc(
   SCIP*                 scip,               /**< SCIP data structure the cons belongs to */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SETPPCTYPE*      setppctype          /**< returns the type of the constraints */
   );

SCIP_Bool GCGgetConsIsCardinalityCons(
   SCIP*                 scip,               /**< SCIP data structure the cons belongs to */
   SCIP_CONS*            cons               /**< constraint to check */
);

/** returns TRUE or FALSE, depending whether we are in the root node or not */
GCG_EXPORT
SCIP_Bool GCGisRootNode(
   GCG*                  gcg                /**< SCIP data structure */
   );

GCG_EXPORT
SCIP_RETCODE GCGincludeDialogsGraph(
   GCG*                  gcg
   );

/**@} */
#ifdef __cplusplus
}
#endif

#endif /* GCG_SCIP_MISC_H_ */
