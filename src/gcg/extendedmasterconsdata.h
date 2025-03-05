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

/**@file    extendedmasterconsdata.h
 * @ingroup INTERNALAPI-GCG
 * @brief   internal methods for interacting with GCG_EXTENDEDMASTERCONSDATA
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_EXTENDEDMASTERCONSDATA_H_
#define GCG_EXTENDEDMASTERCONSDATA_H_



#include "gcg/pricer_gcg.h"
#include "gcg/type_extendedmasterconsdata.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @ingroup GCG_EXTENDEDMASTERCONSDATA
 * @{
 */

/** update the extended master cons with the new dual value */
SCIP_RETCODE GCGextendedmasterconsUpdateDualValue(
   GCG*                          gcg,                    /**< master scip */
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata, /**< extended master cons data */
   SCIP_Real                     dualvalue               /**< dual value */
   );

/** apply a pricing modification */
SCIP_RETCODE GCGpricingmodificationApply(
   SCIP*                         pricingscip,            /**< pricing scip */
   GCG_PRICINGMODIFICATION*      pricingmodification     /**< pricing modification */
   );

/** apply all pricing modifications */
SCIP_RETCODE GCGextendedmasterconsApplyPricingModifications(
   GCG*                          gcg,                    /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata  /**< extended master cons data */
   );

/** undo a pricing modification */
SCIP_RETCODE GCGpricingmodificationUndo(
   SCIP*                         pricingscip,            /**< pricing scip */
   GCG_PRICINGMODIFICATION*      pricingmodification     /**< pricing modification */
   );

/** undo all pricing modifications */
SCIP_RETCODE GCGextendedmasterconsUndoPricingModifications(
   GCG*                          gcg,                    /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata  /**< extended master cons data */
   );

/**@} */
#ifdef __cplusplus
}
#endif

#endif
