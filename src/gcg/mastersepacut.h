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

#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** increases usage counter of separator master cut */
SCIP_RETCODE GCGcaptureMasterSepaCut(
   GCG_SEPARATORMASTERCUT*       mastersepacut      /**< separator master cut */
   );

/** decreases usage counter of separator master cut, and frees memory if necessary */
SCIP_RETCODE GCGreleaseMasterSepaCut(
   GCG*                          gcg,             /**< GCG data structure */
   GCG_SEPARATORMASTERCUT**      mastersepacut    /**< pointer to separator master cut */
   );

/** creates separator master cut */
SCIP_RETCODE GCGcreateMasterSepaCut(
   GCG*                          gcg,                 /**< GCG data structure */
   GCG_SEPARATORMASTERCUT**      mastersepacut,       /**< pointer to store separator master cut */
   GCG_SEPARATORMASTERCUTTYPE    mastersepacuttype,   /**< type of separator master cut */
   GCG_SEPA*                     sepa,                /**< separator creating this cut */
   GCG_EXTENDEDMASTERCONSDATA*            mastercutdata,       /**< master cut data */
   GCG_VARHISTORY*               varhistory,          /**< variable history */
   GCG_SEPARATORMASTERCUTDATA*   mastersepacutdata    /**< separator master cut data */
   );

#ifdef __cplusplus
}
#endif

#endif
