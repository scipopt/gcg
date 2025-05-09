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

/**@file   cons_decomp.hpp
 * @ingroup CONSHDLRS-GCG
 * @brief  C++ interface of cons_decomp
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CONS_DECOMP_HPP
#define GCG_CONS_DECOMP_HPP


#include "gcg/class_partialdecomp.h"

/** @brief gets vector of all partialdecs
 * @returns finished partialdecs
 */
GCG_EXPORT
std::vector<gcg::PARTIALDECOMP*>* GCGconshdlrDecompGetPartialdecs(
   GCG*           gcg  /**< GCG data structure */
   );

GCG_EXPORT
gcg::PARTIALDECOMP* GCGgetPartialdecToWrite(
   GCG*                          gcg,
   SCIP_Bool                     transformed
   );

/** @brief local method to find a partialdec for a given id or NULL if no partialdec with such id is found
 * @returns partialdec pointer of partialdec with given id or NULL if it does not exist
 * @note returns NULL if no partialdec by this id is known */
GCG_EXPORT
gcg::PARTIALDECOMP* GCGconshdlrDecompGetPartialdecFromID(
   GCG* scip,           /**< GCG data structure */
   int partialdecid     /**< partialdec id */
   );

/** @brief adds a preexisting partial dec to be considered at the beginning of the detection
 *
 * @note refines the partialdec to be consistent, adds meta data/statistics, completes partial dec by assigning open conss to master if USERGIVEN::COMPLETED_CONSTOMASTER is set
 * @returns SCIP return code
*/
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompAddPreexisitingPartialDec(
   GCG* gcg,                     /**< GCG data structure */
   gcg::PARTIALDECOMP* partialdec/**< partial dec to add */
   );

/** @brief adds a preexisting partial dec to be considered at the beginning of the detection
 *
 * @note refines the partialdec to be consistent, adds meta data/statistics, completes partial dec by assigning open conss to master if USERGIVEN::COMPLETED_CONSTOMASTER is set
 * @returns SCIP return code
*/
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompAddPreexisitingPartialDec(
   GCG* gcg,                        /**< GCG data structure */
   gcg::PARTIALDECOMP* partialdec,  /**< partial dec to add */
   SCIP_Bool addpartialdec          /**< if completed, add partial dec as well */
   );

/** @brief deregisters a partialdec in the conshdlr
 *
 * Use this function at deletion of the partialdec.
 * The partialdec is not destroyed in this function, the conshdlr will not know that it exists.
 */
GCG_EXPORT
void GCGconshdlrDecompDeregisterPartialdec(
   GCG* gcg,                         /**< GCG data structure */
   gcg::PARTIALDECOMP* partialdec    /**< the partialdec */
   );

/** @brief registers a partialdec in the conshdlr
 *
 * Use this function at initialization of the partialdec.
 * If the partialdec already exists in the conshdlr it is ignored.
 */
GCG_EXPORT
void GCGconshdlrDecompRegisterPartialdec(
   GCG* gcg,                         /**< GCG data structure */
   gcg::PARTIALDECOMP* partialdec    /**< the partialdec to register */
   );

/**
 * @brief help method to access detprobdata for unpresolved problem
 *
 * @returns pointer to detprobdata in wrapper data structure
 */
GCG_EXPORT
gcg::DETPROBDATA* GCGconshdlrDecompGetDetprobdataOrig(
   GCG*                  gcg                  /**< GCG data structure */
   );

/**
 * @brief help method to access detprobdata for transformed problem
 *
 * @returns pointer to detprobdata in wrapper data structure
 */
GCG_EXPORT
gcg::DETPROBDATA* GCGconshdlrDecompGetDetprobdataPresolved(
   GCG*                  gcg                  /**< GCG data structure */
   );

/**
 * @brief initilizes the candidates data structures with selected partialdecs
 *
 * initializes it with all if there are no selected partialdecs,
 * sort them according to the current scoretype
 * @param gcg GCG data structure
 * @param candidates tuples of partialdecs and scores will be added to this vector (sorted w.r.t. the scores).
 * @param original choose candidates for the original problem?
 * @param printwarnings should warnings be printed?
 * @returns SCIP return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompChooseCandidatesFromSelected(
   GCG* gcg,
   std::vector<std::pair<gcg::PARTIALDECOMP*, SCIP_Real> >& candidates,
   SCIP_Bool original,
   SCIP_Bool printwarnings
   );

/** @brief gets detector history of partialdec with given id
 * @returns detector history of partialdec as string
 */
GCG_EXPORT
std::string GCGconshdlrDecompGetDetectorHistoryByPartialdecId(
   GCG* gcg,      /**< GCG data structure */
   int id         /**< id of partialdec */
   );

#endif //GCG_CONS_DECOMP_HPP
