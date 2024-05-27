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


#include <scip/scip.h>

#include "def.h"
#include "mastercutdata.h"
#include "type_varhistory.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates event handler for xyz event */
struct GCG_MasterSepaCut
{
   GCG_MASTERCUTDATA*      mastercutdata;          /**< mastercutdata */
   GCG_VARHISTORY*         knownvarhistory;        /**< pointer to the history of priced variables */
   int                     nuses;                  /**< number of times this cut is referenced */
   int                     n;                      /**< number of constraints used to create cut */
   int*                    conssindices;           /**< indices of constraints used to create cut */
   SCIP_Real*              weights;                /**< weights used to create cut */
};

typedef struct GCG_MasterSepaCut GCG_MASTERSEPACUT;

GCG_EXPORT
SCIP_RETCODE SCIPincludeEventHdlrSepaCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

SCIP_RETCODE GCGremoveNewInactiveRows(
   SCIP* masterscip,    /**< SCIP data structure */
   int* startidx        /**< indicate the first new cut from each separator */
   );

GCG_MASTERSEPACUT*** GCGgetActiveCuts(
   SCIP* masterscip     /**< SCIP data structure */
);

int* GCGgetNActiveCuts(
   SCIP* masterscip     /**< SCIP data structure */
);

SCIP_RETCODE GCGshrinkActiveCuts(
   SCIP* masterscip,    /**< SCIP data structure */
   int* newnrows        /**< indices to which activecuts should be shrunk to */
);

SCIP_RETCODE GCGaddCutActiveCuts(
   SCIP* masterscip,                   /**< SCIP data structure */
   GCG_MASTERSEPACUT* mastersepacut,   /**< master sepa cut */
   int sepaidx                         /**< index of the separator which generated the cut */
);

SCIP_RETCODE GCGreleaseMasterSepaCut(
   SCIP* masterscip,                    /**< SCIP data structure */
   GCG_MASTERSEPACUT** mastersepacut    /**< pointer to master sepa cut */
);

SCIP_RETCODE GCGcaptureMasterSepaCut(
   GCG_MASTERSEPACUT* mastersepacut   /**< MASTERSEPACUT data structure */
);

SCIP_RETCODE GCGaddCutToGeneratedCutsSepa(
   SCIP*                masterscip,       /**< SCIP data structure */
   GCG_MASTERCUTDATA*   mastercutdata,    /**< mastercut data */
   SCIP_Real*           weights,          /**< weights used to create the cut */
   int*                 conssindices,     /**< indices of constraints used to create the cut */
   int                  n,                /**< number of constraints used to create the cut */
   int                  sepaidx           /**< index of the separator which generated the cut */
);

SCIP_RETCODE GCGclearGeneratedCuts(
   SCIP* masterscip     /**< SCIP data structure */
);

GCG_MASTERCUTDATA* GCGsepamastercutGetMastercutData(
   GCG_MASTERSEPACUT* mastersepacut
);

SCIP_Real* GCGsepamastercutGetWeights(
   GCG_MASTERSEPACUT* mastersepacut
);

int* GCGsepamastercutGetConssIndices(
   GCG_MASTERSEPACUT* mastersepacut
);

int GCGsepamastercutGetNWeights(
   GCG_MASTERSEPACUT* mastersepacut
);

#ifdef __cplusplus
}
#endif

#endif
