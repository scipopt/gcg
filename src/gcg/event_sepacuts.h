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

/**@file        event_sepacuts.h
 * @ingroup     EVENTS
 * @brief       eventhdlr for xyz event
 * @author      Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_SEPACUTS_H__
#define __SCIP_EVENT_SEPACUTS_H__


#include "scip/scip.h"
#include "def.h"
#include "mastercutdata.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates event handler for xyz event */
struct GCG_StoredCut
{
   GCG_MASTERCUTDATA*      mastercutdata;          /**< mastercutdata */
   GCG_VARHISTORY*         knownvarhistory;        /**< pointer to the history of priced variables */
   int                     nuses;
};

typedef struct GCG_StoredCut GCG_STOREDCUT;

GCG_EXPORT
SCIP_RETCODE SCIPincludeEventHdlrSepaCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

SCIP_RETCODE GCGremoveNewInactiveRows(
   SCIP* scip,
   int* startidx
   );

GCG_STOREDCUT*** GCGgetActiveCuts(
   SCIP* scip
);

int* GCGgetNActiveCuts(
   SCIP* scip
);

SCIP_RETCODE GCGshrinkActiveCuts(
   SCIP* scip,
   int* newnrows
);

SCIP_RETCODE GCGaddCutTActiveCuts(
   SCIP* scip,
   GCG_STOREDCUT* storedcut,
   int sepaidx
);

SCIP_RETCODE GCGreleaseStoredCut(
   SCIP* scip,
   GCG_STOREDCUT** storedcut
);

SCIP_RETCODE GCGcaptureStoredCut(
   GCG_STOREDCUT* storedcut
);

SCIP_RETCODE GCGaddCutToGeneratedCutsSepa(
   SCIP* scip,
   GCG_MASTERCUTDATA* mastercutdata,
   int sepaidx
);

SCIP_RETCODE GCGclearGeneratedCuts(
   SCIP* scip
);

SCIP_RETCODE GCGmapNodeToConsdata(
   SCIP* scip,
   SCIP_NODE* node,
   SCIP_CONSDATA* consdata
);
#ifdef __cplusplus
}
#endif

#endif
