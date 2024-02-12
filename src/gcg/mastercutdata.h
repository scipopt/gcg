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

/**@file    mastercutdata.h
 * @ingroup TODO-????
 * @brief   methods for interacting with GCG_MASTERCUTDATA
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_MASTERCUTDATA_H_
#define GCG_MASTERCUTDATA_H_

#include "def.h"
#include "type_mastercutdata.h"

#include "scip/scip.h"
#include <scip/type_retcode.h>
#include <scip/type_scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @ingroup TODO-????
 *
 * @{
 */

/** get the blocknr of a mastercut */
GCG_EXPORT
int GCGmastercutGetBlocknr(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/** determine whether the mastercutdata is active in the masterscip */
GCG_EXPORT
SCIP_Bool GCGmastercutIsActive(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );


/** activate the mastercutdata */
GCG_EXPORT
SCIP_RETCODE GCGmastercutActivate(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/** deactivate the mastercutdata */
GCG_EXPORT
SCIP_RETCODE GCGmastercutDeactivate(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/**@} */
#ifdef __cplusplus
}
#endif

#endif
