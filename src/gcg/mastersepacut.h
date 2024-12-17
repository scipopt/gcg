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

/** increases usage counter of master separator cut */
SCIP_RETCODE GCGcaptureMasterSepaCut(
   GCG_SEPARATORMASTERCUT*      mastersepacut      /**< master separator cut */
);

/** decreases usage counter of master separator cut, and frees memory if necessary */
SCIP_RETCODE GCGreleaseMasterSepaCut(
   SCIP*                   masterscip,      /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT**     mastersepacut    /**< pointer to master separator cut */
);

/** creates master separator cut */
SCIP_RETCODE GCGcreateMasterSepaCut(
   SCIP*                   masterscip,          /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT**     mastersepacut,       /**< pointer to store master separator cut */
   GCG_SEPARATORMASTERCUTTYPE   mastersepacuttype,   /**< type of master separator cut */
   GCG_SEPA*               sepa,                /**< separator creating this cut */
   GCG_MASTERCUTDATA*      mastercutdata,       /**< master cut data */
   GCG_VARHISTORY*         varhistory,          /**< variable history */
   GCG_SEPARATORMASTERCUTDATA*  mastersepacutdata    /**< master separator cut data */
);

/** returns the master cut data of the master separator cut */
GCG_MASTERCUTDATA* GCGmastersepacutGetMasterCutData(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< master separator cut */
);

/** returns the variable history of the master separator cut */
GCG_VARHISTORY* GCGmastersepacutGetVarHistory(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< master separator cut */
);

/** returns the cut type of the master separator cut */
GCG_SEPARATORMASTERCUTTYPE GCGmastersepacutGetCutType(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< master separator cut */
);

/** return the separator which created the cut */
GCG_SEPA* GCGmastersepacutGetSeparator(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< master separator cut */
);

/** returns the data of the master separator cut */
GCG_SEPARATORMASTERCUTDATA* GCGmastersepacutGetData(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< master separator cut */
);

/** set the variable history of master separator cut */
SCIP_RETCODE GCGmastersepacutSetVarHistory(
   SCIP*                   masterscip,       /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT**     mastersepacut     /**< pointer to master separator cut */
);

/** creates a subset row cut */
SCIP_RETCODE GCGcreateChvatalGomoryCut(
   SCIP*                   masterscip,            /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT**     mastersepacut,         /**< pointer to store master separator cut */
   GCG_SEPA*               sepa,                  /**< separator creating this cut */
   GCG_MASTERCUTDATA*      mastercutdata,         /**< mastercutdata associated with the cut */
   GCG_VARHISTORY*         varhistory,            /**< variables history of subset row cut*/
   SCIP_Real*              weights,               /**< weights which were used to create the cut */
   int*                    indices,               /**< indices of constraints used to create the cut */
   int                     n                      /**< number of constraints used to create the cut */
);

/** returns TRUE or FALSE whether cut is a subset row cut */
SCIP_Bool GCGmastersepacutIsChvatalGomory(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< master separator cut */
);

/** returns the number of weights of subset row cut */
int GCGchvatalGomoryCutGetNWeights(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< master separator cut */
);

/** returns the weights of subset row cut */
SCIP_Real* GCGchvatalGomoryCutGetWeights(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< master separator cut */
);

/** returns the constraint indices of subset row cut */
int* GCGchvatalGomoryCutGetConssIndices(
   GCG_SEPARATORMASTERCUT*      mastersepacut     /**< master separator cut */
);

/** computes the coefficient of a column for a master separator cut */
SCIP_RETCODE GCGchvatalGomoryCutGetColumnCoefficient(
   SCIP*                   scip,       /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT*      cut,        /**< master separator cut */
   GCG_COL*                gcgcol,     /**< gcg column */
   SCIP_Real*              coeff       /**< pointer to store the coefficient */
);

/** computes the coefficient of a master variable for a master separator cut */
SCIP_RETCODE GCGchvatalGomoryCutGetVariableCoefficient(
   SCIP*                scip,       /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT*   cut,        /**< master separator cut */
   SCIP_VAR**           vars,       /**< pricing variables which define the master variable */
   SCIP_Real*           vals,       /**< values of the pricing variables which define the master variables */
   int                  nvars,      /**< number of pricing variables which define the master variable */
   int                  probnr,     /**< index of the pricing problem which generated the master variable */
   SCIP_Real*           coef        /**< pointer to store the coefficient */
);

/** adapts the objectives of all the necessary pricing problems such that they consider the master cut */
SCIP_RETCODE GCGchvatalGomorySetPricingObjectives(
   SCIP*                   scip,    /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT*      cut,     /**< master separator cut */
   SCIP_Real               dual     /**< the dual value of the master separator cut */
);

/** adapts a GCG column such that it respects the pricing modification imposed by the master separator cut */
SCIP_RETCODE GCGchvatalGomoryAdjustGCGColumn(
   SCIP*                scip,       /**< SCIP data structure (master problem) */
   GCG_SEPARATORMASTERCUT*   cut,        /**< master separator cut */
   GCG_COL**            gcgcol      /**< gcg column */
);


#ifdef __cplusplus
}
#endif

#endif
