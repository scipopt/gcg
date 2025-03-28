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

/**@file   type_detector.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for detectors in GCG projects
 * @author Martin Bergner
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_TYPE_DETECTOR_H__
#define GCG_TYPE_DETECTOR_H__

#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"
#include "gcg/type_decomp.h"
#include "gcg/type_gcg.h"


typedef struct GCG_Detector GCG_DETECTOR;
typedef struct GCG_DetectorData GCG_DETECTORDATA;


struct Partialdec_Detection_Data;
typedef struct Partialdec_Detection_Data PARTIALDEC_DETECTION_DATA;


/** destructor of detector to free user data (called when GCG is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - detector        : detector data structure
 */
#define GCG_DECL_FREEDETECTOR(x) SCIP_RETCODE x (GCG* gcg, GCG_DETECTOR* detector)

/**
 * detector initialization method (called after problem was transformed)
 * It can be used to fill the detector data with needed information. The implementation is optional.
 *
 * input:
 *  - scip            : SCIP data structure
 *  - detector        : detector data structure
 */
#define GCG_DECL_INITDETECTOR(x) SCIP_RETCODE x (GCG* gcg, GCG_DETECTOR* detector)

/**
 * detector deinitialization method (called before the transformed problem is freed)
 * It can be used to clean up the data created in DEC_DECL_INITDETECTOR. The implementation is optional.
 *
 * input:
 *  - scip            : SCIP data structure
 *  - detector        : detector data structure
 */
#define GCG_DECL_EXITDETECTOR(x) SCIP_RETCODE x (GCG* gcg, GCG_DETECTOR* detector)


/**
 * given a partialdec (incomplete decomposition) the detector
 * tries to find refined partialdec and stores it
 *
 * input:
 *  - scip                 : SCIP data structure
 *  - detector             : pointer to detector
 *  - partialdecdetectiondata   : pointer to partialdec propagation data structure
 *  - result               : pointer where to store the result
 *
 * possible return values for result:
 *  - SCIP_SUCCESS    : the method completed and found decompositions
 *  - SCIP_DIDNOTFIND : the method completed without finding a decomposition
 *  - SCIP_DIDNOTRUN  : the method did not run
 */
#define GCG_DECL_PROPAGATEPARTIALDEC(x) SCIP_RETCODE x (GCG* gcg, GCG_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result)


/**
 * given a partialdec (incomplete decomposition) the detector
 * tries to find finished partialdecs and stores them
 *
 * input:
 *  - scip                  : SCIP data structure
 *  - detector              : pointer to detector
 *  - partialdecdetectiondata    : pointer to partialdec propagation data structure
 *  - result                : pointer where to store the result
 *
 * possible return values for result:
 *  - SCIP_SUCCESS    : the method completed and found decompositions
 *  - SCIP_DIDNOTFIND : the method completed without finding a decomposition
 *  - SCIP_DIDNOTRUN  : the method did not run
 */
#define GCG_DECL_FINISHPARTIALDEC(x) SCIP_RETCODE x (GCG* gcg, GCG_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result)

/**
 * given a complete partialdec (complete decomposition) the detector
 * postprocess the partialdec in order to find a different yet promising partialdec
 *
 * input:
 *  - scip                  : SCIP data structure
 *  - detector              : pointer to detector
 *  - partialdecdetectiondata    : pointer to partialdec propagation data structure
 *  - result                : pointer where to store the result
 *
 * possible return values for result:
 *  - SCIP_SUCCESS    : the method completed and found decompositions
 *  - SCIP_DIDNOTFIND : the method completed without finding a decomposition
 *  - SCIP_DIDNOTRUN  : the method did not run
 */
#define GCG_DECL_POSTPROCESSPARTIALDEC(x) SCIP_RETCODE x (GCG* gcg, GCG_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result)

/**
 * set the parameter of a detector to values according to fast emphasis and size of the instance
 *  input:
 *  - scip            : SCIP data structure
 *  - detector        : pointer to detector
 *  - result          : pointer where to store the result
 *
 * possible return values for result:
 *  - SCIP_SUCCESS    : the method completed
 *  - SCIP_DIDNOTRUN  : the method did not run
 */
#define GCG_DECL_SETPARAMFAST(x) SCIP_RETCODE x (GCG* gcg, GCG_DETECTOR* detector, SCIP_RESULT* result)

/**
 * set the parameter of a detector to values according to aggressive emphasis and size of the instance
 *  input:
 *  - scip            : SCIP data structure
 *  - detector        : pointer to detector
 *  - result          : pointer where to store the result
 *
 * possible return values for result:
 *  - SCIP_SUCCESS    : the method completed
 *  - SCIP_DIDNOTRUN  : the method did not run
 */
#define GCG_DECL_SETPARAMAGGRESSIVE(x) SCIP_RETCODE x (GCG* gcg, GCG_DETECTOR* detector, SCIP_RESULT* result)

/**
 * set the parameter of a detector to values according to default emphasis and size of the instance
 *  input:
 *  - scip            : SCIP data structure
 *  - detector        : pointer to detector
 *  - result          : pointer where to store the result
 *
 * possible return values for result:
 *  - SCIP_SUCCESS    : the method completed
 *  - SCIP_DIDNOTRUN  : the method did not run
 */
#define GCG_DECL_SETPARAMDEFAULT(x) SCIP_RETCODE x (GCG* gcg, GCG_DETECTOR* detector, SCIP_RESULT* result)



#endif
