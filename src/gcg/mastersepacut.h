/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef GCG_MASTERSEPACUT_H__
#define GCG_MASTERSEPACUT_H__

#include "struct_mastersepacut.h"

#ifdef __cplusplus
extern "C" {
#endif

/** increases usage counter of separator master cut */
SCIP_RETCODE GCGcaptureMasterSepaCut(
   GCG_SEPARATORMASTERCUT*       mastersepacut      /**< separator master cut */
   );

/** decreases usage counter of separator master cut, and frees memory if necessary */
SCIP_RETCODE GCGreleaseMasterSepaCut(
   SCIP*                         masterscip,      /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT**      mastersepacut    /**< pointer to separator master cut */
   );

/** creates separator master cut */
SCIP_RETCODE GCGcreateMasterSepaCut(
   SCIP*                         masterscip,          /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT**      mastersepacut,       /**< pointer to store separator master cut */
   GCG_SEPARATORMASTERCUTTYPE    mastersepacuttype,   /**< type of separator master cut */
   GCG_SEPA*                     sepa,                /**< separator creating this cut */
   GCG_MASTERCUTDATA*            mastercutdata,       /**< master cut data */
   GCG_VARHISTORY*               varhistory,          /**< variable history */
   GCG_SEPARATORMASTERCUTDATA*   mastersepacutdata    /**< separator master cut data */
   );

/** returns the master cut data of the separator master cut */
GCG_MASTERCUTDATA* GCGmastersepacutGetMasterCutData(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   );

/** returns the variable history of the separator master cut */
GCG_VARHISTORY* GCGmastersepacutGetVarHistory(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   );

/** returns the cut type of the separator master cut */
GCG_SEPARATORMASTERCUTTYPE GCGmastersepacutGetCutType(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   );

/** return the separator which created the cut */
GCG_SEPA* GCGmastersepacutGetSeparator(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   );

/** returns the data of the separator master cut */
GCG_SEPARATORMASTERCUTDATA* GCGmastersepacutGetData(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   );

/** set the variable history of separator master cut */
SCIP_RETCODE GCGmastersepacutSetVarHistory(
   SCIP*                         masterscip,       /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT**      mastersepacut     /**< pointer to separator master cut */
   );

/** creates a Chvatal-Gomory cut */
SCIP_RETCODE GCGcreateChvatalGomoryCut(
   SCIP*                         masterscip,            /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT**      mastersepacut,         /**< pointer to store separator master cut */
   GCG_SEPA*                     sepa,                  /**< separator creating this cut */
   GCG_MASTERCUTDATA*            mastercutdata,         /**< mastercutdata associated with the cut */
   GCG_VARHISTORY*               varhistory,            /**< variables history of subset row cut*/
   SCIP_Real*                    weights,               /**< weights which were used to create the cut */
   int*                          indices,               /**< indices of constraints used to create the cut */
   int                           n                      /**< number of constraints used to create the cut */
   );

/** returns TRUE or FALSE whether cut is a Chvatal-Gomory cut */
SCIP_Bool GCGmastersepacutIsChvatalGomory(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   );

/** returns the number of weights of Chvatal-Gomory cut */
int GCGchvatalGomoryCutGetNWeights(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   );

/** returns the weights of Chvatal-Gomory cut */
SCIP_Real* GCGchvatalGomoryCutGetWeights(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   );

/** returns the constraint indices of Chvatal-Gomory cut */
int* GCGchvatalGomoryCutGetConssIndices(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< separator master cut */
   );

/** computes the coefficient of a column for a Chvatal-Gomory cut */
SCIP_RETCODE GCGchvatalGomoryCutGetColumnCoefficient(
   SCIP*                         scip,       /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT*       cut,        /**< separator master cut */
   GCG_COL*                      gcgcol,     /**< gcg column */
   SCIP_Real*                    coeff       /**< pointer to store the coefficient */
   );

/** computes the coefficient of a master variable for a Chvatal-Gomory cut */
SCIP_RETCODE GCGchvatalGomoryCutGetVariableCoefficient(
   SCIP*                      scip,       /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT*    cut,        /**< separator master cut */
   SCIP_VAR**                 vars,       /**< pricing variables which define the master variable */
   SCIP_Real*                 vals,       /**< values of the pricing variables which define the master variables */
   int                        nvars,      /**< number of pricing variables which define the master variable */
   int                        probnr,     /**< index of the pricing problem which generated the master variable */
   SCIP_Real*                 coef        /**< pointer to store the coefficient */
   );

/** adapts the objectives of all the necessary pricing problems such that they consider the Chvatal-Gomory cut */
SCIP_RETCODE GCGchvatalGomorySetPricingObjectives(
   SCIP*                         scip,    /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT*       cut,     /**< separator master cut */
   SCIP_Real                     dual     /**< the dual value of the separator master cut */
   );

/** adapts a GCG column such that it respects the pricing modification imposed by the Chvatal-Gomory cut */
SCIP_RETCODE GCGchvatalGomoryAdjustGCGColumn(
   SCIP*                      scip,       /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT*    cut,        /**< separator master cut */
   GCG_COL**                  gcgcol      /**< gcg column */
   );


#ifdef __cplusplus
}
#endif

#endif
