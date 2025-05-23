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

/**@file   dec_xyz.cpp
 * 
 * @brief  detector xyz (put your description here)
 * @author Martin Bergner
 * @author Christian Puchert
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/dec_xyz.h"
#include "gcg/cons_decomp.h"

/* detector properties */
#define DEC_NAME                  "xyz"       /**< name of detector */
#define DEC_DESC                  "detector xyz" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /**< frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /**< last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /**< first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /**< frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /**< last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /**< first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               '?'         /**< display character of detector */
#define DEC_ENABLED               TRUE        /**< should the detection (propagating) be enabled */
#define DEC_ENABLEDFINISHING      FALSE       /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated partialdec  */

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
#define detectorFreeXyz NULL

/** detector initialization method (called after problem was transformed) */
#define detectorInitXyz NULL

/** detector deinitialization method (called before the transformed problem is freed) */
#define detectorExitXyz NULL

/** destructor of detector to free user data (called when GCG is exiting) */
#define detectorFreeXyz NULL

/** propagate PARTIALDECOMP function of detector */
#ifdef SCIP_DISABLED_CODE
static DEC_DECL_PROPAGATEPARTIALDEC(detectorPropagatePartialdecXyz)
{  /*lint --e{715}*/

   SCIPerrorMessage("Propagate PARTIALDECOMP function of detector <%s> not implemented!\n", DEC_NAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define detectorPropagatePartialdecXyz NULL
#endif

#define detectorFinishPartialdecXyz NULL
#define detectorPostprocessPartialdecXyz NULL

#define setParamAggressiveXyz NULL
#define setParamDefaultXyz NULL
#define setParamFastXyz NULL



/*
 * detector specific interface methods
 */

/** creates the handler for xyz detector and includes it in SCIP */
SCIP_RETCODE GCGincludeDetectorXyz(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_DETECTORDATA* detectordata;

   /**@todo create xyz detector data here*/
   detectordata = NULL;

   SCIP_CALL( GCGincludeDetector(gcg, DEC_NAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata, detectorFreeXyz, detectorInitXyz, detectorExitXyz, detectorPropagatePartialdecXyz, detectorFinishPartialdecXyz, detectorPostprocessPartialdecXyz, setParamAggressiveXyz, setParamDefaultXyz, setParamFastXyz) );

   /**@todo add xyz detector parameters */

   return SCIP_OKAY;
}
