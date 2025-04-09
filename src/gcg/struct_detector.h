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

/**@file   struct_detector.h
 * @ingroup DATASTRUCTURES
 * @brief  data structures for detectors
 * @author Martin Bergner
 * @author Christian Puchert
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_DETECTOR_H__
#define GCG_STRUCT_DETECTOR_H__

#include "gcg/type_detector.h"



/** detector data structure */
struct GCG_Detector {
   const char*           name;               /**< name of the detector */
   GCG_DETECTORDATA*     decdata;            /**< custom data structure of the detectors */
   char                  decchar;            /**< display character of detector */
   const char*           description;        /**< description of the detector */
   int                   freqCallRound;      /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
   int                   maxCallRound;       /** last round the detector gets called                              */
   int                   minCallRound;       /** first round the detector gets called (offset in detection loop) */
   int                   freqCallRoundOriginal; /** frequency the detector gets called in detection loop while detecting the original problem */
   int                   maxCallRoundOriginal; /** last round the detector gets called while detecting the original problem */
   int                   minCallRoundOriginal; /** first round the detector gets calles (offset in detection loop) while detecting the original problem */
   int                   priority;           /**< detector priority */
   SCIP_Bool             enabled;            /**< flag to indicate whether detector is enabled */
   SCIP_Bool             enabledFinishing;   /**< flag to indicate whether finishing is enabled */
   SCIP_Bool             enabledPostprocessing; /**< flag to indicate whether finishing is enabled */
   SCIP_Bool             skip;               /**< should detector be skipped if other detectors found decompositions */
   SCIP_Bool             usefulRecall;       /** is it useful to call this detector on a descendant of the propagated partialdec */
   SCIP_Bool             overruleemphasis;   /**< should the emphasis settings be overruled */
   int                   ndecomps;           /**< number of decompositions the detector has worked on */
   int                   ncompletedecomps;   /**< number of complete decompositions the detector has worked on (including decompositions that were finished by other detectors) */
   SCIP_Real             dectime;            /**< time the detector took to find decompositions */

   GCG_DECL_FREEDETECTOR((*freeDetector));                  /**< destructor of detector */
   GCG_DECL_INITDETECTOR((*initDetector));                  /**< initialization method of detector */
   GCG_DECL_EXITDETECTOR((*exitDetector));                  /**< deinitialization method of detector */
   GCG_DECL_EXITDETECTOR((*exitDetection));                 /**< deinitialization method of detector */
   GCG_DECL_PROPAGATEPARTIALDEC((*propagatePartialdec));    /**< propagation method of detector (or NULL) */
   GCG_DECL_FINISHPARTIALDEC((*finishPartialdec));          /**< finish method of detector (or NULL) */
   GCG_DECL_POSTPROCESSPARTIALDEC((*postprocessPartialdec)); /**< postprocess method of detector (or NULL) */
   GCG_DECL_SETPARAMAGGRESSIVE((*setParamAggressive));      /**< set method for aggressive parameters of detector (or NULL) */
   GCG_DECL_SETPARAMDEFAULT((*setParamDefault));            /**< set method for default parameters of detector (or NULL) */
   GCG_DECL_SETPARAMFAST((*setParamFast));                  /**< set method for fast parameters of detector (or NULL) */


};


#endif
