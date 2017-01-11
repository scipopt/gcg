/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2015 Operations Research, RWTH Aachen University       */
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

/**@file   dec_compgreedily.cpp
 * @ingroup DETECTORS
 * @brief  detector compgreedily (assigns the open cons and open vars of the seeed greedily)
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_compgreedily.h"
#include "cons_decomp.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "scip/clock.h"
#include <iostream>

/* constraint handler properties */
#define DEC_DETECTORNAME         "compgreedily"       /**< name of detector */
#define DEC_DESC                 "detector compgreedily" /**< description of detector*/
#define DEC_FREQCALLROUND        1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND         INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND         0           /** first round the detector gets called                              */
#define DEC_PRIORITY             0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              'g'         /**< display character of detector */
#define DEC_ENABLED              TRUE        /**< should the detection be enabled */
#define DEC_ENABLEDFINISHING     TRUE        /**< should the finishing be enabled */
#define DEC_SKIP                 FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL         FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */



/*
 * Data structures
 */

/** @todo fill in the necessary detector data */

/** detector handler data */
struct DEC_DetectorData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */

#define freeCompgreedily NULL

/** destructor of detector to free detector data (called before the solving process begins) */
#if 0
static
DEC_DECL_EXITDETECTOR(exitCompgreedily)
{  /*lint --e{715}*/

   SCIPerrorMessage("Exit function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define exitCompgreedily NULL
#endif

/** detection initialization function of detector (called before solving is about to begin) */
#if 0
static
DEC_DECL_INITDETECTOR(initCompgreedily)
{  /*lint --e{715}*/

   SCIPerrorMessage("Init function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define initCompgreedily NULL
#endif

/** detection function of detector */
static
DEC_DECL_DETECTSTRUCTURE(detectCompgreedily)
{ /*lint --e{715}*/
   *result = SCIP_DIDNOTFIND;

   SCIPerrorMessage("Detection function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();  /*lint --e{527}*/

   return SCIP_OKAY;
}


static
DEC_DECL_PROPAGATESEEED(propagateSeeedCompgreedily)
{
   *result = SCIP_DIDNOTFIND;

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );

   gcg::Seeed* seeed;
   seeed = new gcg::Seeed(seeedPropagationData->seeedToPropagate, seeedPropagationData->seeedpool);

   //assign open conss and vars greedily
   seeed->completeGreedily(seeedPropagationData->seeedpool);


   seeed->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForFinishingDetector(detector));

   SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), 1) );
   seeedPropagationData->newSeeeds[0] = seeed;
   seeedPropagationData->nNewSeeeds = 1;

   seeedPropagationData->newSeeeds[0]->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForDetector(detector));
   SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
   seeedPropagationData->newSeeeds[0]->addClockTime( SCIPclockGetTime(temporaryClock )  );
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

static
DEC_DECL_FINISHSEEED(finishSeeedCompgreedily)
{
   *result = SCIP_DIDNOTFIND;

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );

   gcg::Seeed* seeed;
   seeed = new gcg::Seeed(seeedPropagationData->seeedToPropagate, seeedPropagationData->seeedpool);

   //assign open conss and vars greedily
   seeed->completeGreedily(seeedPropagationData->seeedpool);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), 1) );
   seeedPropagationData->newSeeeds[0] = seeed;
   seeedPropagationData->nNewSeeeds = 1;
   seeedPropagationData->newSeeeds[0]->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForDetector(detector));

   SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
   seeedPropagationData->newSeeeds[0]->addClockTime( SCIPclockGetTime(temporaryClock )  );
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
/*
 * detector specific interface methods
 */

/** creates the handler for compgreedily detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorCompgreedily(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /**@todo create compgreedily detector data here*/
   detectordata = NULL;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_SKIP, DEC_USEFULRECALL, detectordata, detectCompgreedily, freeCompgreedily,initCompgreedily, exitCompgreedily, propagateSeeedCompgreedily, finishSeeedCompgreedily) );

   /**@todo add compgreedily detector parameters */

   return SCIP_OKAY;
}
