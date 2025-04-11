/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   mastersepacut.c
 * @brief  public functions to work with master cuts
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_PUB_MASTERSEPACUT_H__
#define GCG_PUB_MASTERSEPACUT_H__

#include "gcg/type_mastersepacut.h"

#include "scip/scip.h"
#include "gcg/def.h"
#include "gcg/type_gcg.h"
#include "gcg/type_gcgvarhistory.h"

#ifdef NDEBUG
#include "gcg/struct_mastersepacut.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** returns the master cut data of the separator master cut */
GCG_EXPORT
GCG_EXTENDEDMASTERCONSDATA* GCGmastersepacutGetMasterCutData(
   GCG_SEPARATORMASTERCUT*       mastersepacut     /**< separator master cut */
   );
 
/** returns the variable history of the separator master cut */
GCG_EXPORT
GCG_VARHISTORY* GCGmastersepacutGetVarHistory(
   GCG_SEPARATORMASTERCUT*       mastersepacut     /**< separator master cut */
   );
 
/** returns the cut type of the separator master cut */
GCG_EXPORT
GCG_SEPARATORMASTERCUTTYPE GCGmastersepacutGetCutType(
   GCG_SEPARATORMASTERCUT*       mastersepacut     /**< separator master cut */
   );
 
/** return the separator which created the cut */
GCG_EXPORT
GCG_SEPA* GCGmastersepacutGetSeparator(
   GCG_SEPARATORMASTERCUT*       mastersepacut     /**< separator master cut */
   );
 
/** returns the data of the separator master cut */
GCG_EXPORT
GCG_SEPARATORMASTERCUTDATA* GCGmastersepacutGetData(
   GCG_SEPARATORMASTERCUT*       mastersepacut     /**< separator master cut */
   );
 
/** set the variable history of separator master cut */
GCG_EXPORT
SCIP_RETCODE GCGmastersepacutSetVarHistory(
   GCG*                          gcg,              /**< GCG data structure */
   GCG_SEPARATORMASTERCUT**      mastersepacut     /**< pointer to separator master cut */
   );

/** creates a Chvatal-Gomory cut */
GCG_EXPORT
SCIP_RETCODE GCGcreateChvatalGomoryCut(
   GCG*                          gcg,                   /**< GCG data structure */
   GCG_SEPARATORMASTERCUT**      mastersepacut,         /**< pointer to store separator master cut */
   GCG_SEPA*                     sepa,                  /**< separator creating this cut */
   GCG_EXTENDEDMASTERCONSDATA*   mastercutdata,         /**< mastercutdata associated with the cut */
   GCG_VARHISTORY*               varhistory,            /**< variables history of subset row cut*/
   SCIP_Real*                    weights,               /**< weights which were used to create the cut */
   int*                          indices,               /**< indices of constraints used to create the cut */
   int                           n                      /**< number of constraints used to create the cut */
   );
 
/** returns TRUE or FALSE whether cut is a Chvatal-Gomory cut */
GCG_EXPORT
SCIP_Bool GCGmastersepacutIsChvatalGomory(
   GCG_SEPARATORMASTERCUT*       mastersepacut     /**< separator master cut */
   );
 
/** returns the number of weights of Chvatal-Gomory cut */
GCG_EXPORT
int GCGchvatalGomoryCutGetNWeights(
   GCG_SEPARATORMASTERCUT*       mastersepacut     /**< separator master cut */
   );
 
/** returns the weights of Chvatal-Gomory cut */
GCG_EXPORT
SCIP_Real* GCGchvatalGomoryCutGetWeights(
   GCG_SEPARATORMASTERCUT*       mastersepacut     /**< separator master cut */
   );
 
/** returns the constraint indices of Chvatal-Gomory cut */
GCG_EXPORT
int* GCGchvatalGomoryCutGetConssIndices(
   GCG_SEPARATORMASTERCUT*       mastersepacut     /**< separator master cut */
   );
 
/** computes the coefficient of a column for a Chvatal-Gomory cut */
GCG_EXPORT
SCIP_RETCODE GCGchvatalGomoryCutGetColumnCoefficient(
   GCG*                          gcg,        /**< GCG data structure */
   GCG_SEPARATORMASTERCUT*       cut,        /**< separator master cut */
   GCG_COL*                      gcgcol,     /**< gcg column */
   SCIP_Real*                    coeff       /**< pointer to store the coefficient */
   );
 
/** computes the coefficient of a master variable for a Chvatal-Gomory cut */
GCG_EXPORT
SCIP_RETCODE GCGchvatalGomoryCutGetVariableCoefficient(
   GCG*                          gcg,        /**< GCG data structure */
   GCG_SEPARATORMASTERCUT*       cut,        /**< separator master cut */
   SCIP_VAR**                    vars,       /**< pricing variables which define the master variable */
   SCIP_Real*                    vals,       /**< values of the pricing variables which define the master variables */
   int                           nvars,      /**< number of pricing variables which define the master variable */
   int                           probnr,     /**< index of the pricing problem which generated the master variable */
   SCIP_Real*                    coef        /**< pointer to store the coefficient */
   );
 
/** adapts the objectives of all the necessary pricing problems such that they consider the Chvatal-Gomory cut */
GCG_EXPORT
SCIP_RETCODE GCGchvatalGomorySetPricingObjectives(
   GCG*                          gcg,     /**< GCG data structure */
   GCG_SEPARATORMASTERCUT*       cut,     /**< separator master cut */
   SCIP_Real                     dual     /**< the dual value of the separator master cut */
   );
 
/** adapts a GCG column such that it respects the pricing modification imposed by the Chvatal-Gomory cut */
GCG_EXPORT
SCIP_RETCODE GCGchvatalGomoryAdjustGCGColumn(
   GCG*                          gcg,        /**< GCG data structure */
   GCG_SEPARATORMASTERCUT*       cut,        /**< separator master cut */
   GCG_COL*                      gcgcol      /**< gcg column */
   );

/** returns active master separator cuts */
GCG_EXPORT
GCG_SEPARATORMASTERCUT** GCGgetActiveCuts(
   GCG*              gcg             /**< GCG data structure */
   );
 
/** returns active master separator cuts */
GCG_EXPORT
GCG_SEPARATORMASTERCUT** GCGsepacutGetActiveCuts(
   GCG*              gcg,            /**< GCG data structure */
   SCIP_EVENTHDLR*   eventhdlr       /**< mastersepacut event handler*/
   );
 
/** return number of active separator mastercuts */
GCG_EXPORT
int GCGgetNActiveCuts(
   GCG*              gcg             /**< GCG data structure */
   );
 
/** return number of active master separator cuts */
GCG_EXPORT
int GCGsepacutGetNActiveCuts(
   GCG*              gcg,            /**< GCG data structure */
   SCIP_EVENTHDLR*   eventhdlr       /**< masterepacut event handler */
   );

#ifdef __cplusplus
}
#endif

#endif