/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   reader_tex.h
 * @brief  tex file reader for writing decomposition details to LaTeX files
 * @author Hanna Franzen
 * @ingroup FILEREADERS-GCG

 * This reader can write visualizations, family trees and reports of partialdecs to a .tex LaTeX file.
 * The gp reader might be required for visualizations.

 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_READER_TEX_H__
#define GCG_READER_TEX_H__

#include "scip/scip.h"

#include "gcg/type_decomp.h"
#include "gcg/cons_decomp.h"

#ifdef __cplusplus
extern "C" {
#endif


/** Includes the tex file reader into SCIP
 *
 * @returns SCIP status */
GCG_EXPORT
SCIP_RETCODE GCGincludeReaderTex(
   GCG*  gcg      /**< GCG data structure */
   );

/** Writes visualization LaTeX code for the given partialdec
 *
 * @returns SCIP status */
GCG_EXPORT
SCIP_RETCODE GCGwriteTexVisualization(
   GCG* gcg,               /**< GCG data structure */
   FILE* file,             /**< file in which to write the LaTeX code */
   int partialdecid,            /**< id of partialdec to visualize */
   SCIP_Bool statistics,   /**< additionally to picture show statistics */
   SCIP_Bool usegp         /**< true if the gp reader should be used for the image generation (instead of tikz) */
   );

/** Writes a report for the given partialdecs
 *
 * @note  *npartialdecs will be set to the number of actually written decompositions.
 * @returns SCIP status */
GCG_EXPORT
SCIP_RETCODE GCGwriteTexReport(
   GCG* gcg,               /**< GCG data structure */
   FILE* file,             /**< file in which to put the LaTeX code */
   int* partialdecids,     /**< ids of partialdecs to visualize */
   int* npartialdecs,      /**< number of partialdecs to visualize */
   SCIP_Bool titlepage,    /**< true if a title page should be included in the document */
   SCIP_Bool toc,          /**< true if an interactive table of contents should be included */
   SCIP_Bool statistics,   /**< true if statistics for each partialdec should be included */
   SCIP_Bool usegp         /**< true if the gp reader should be used for the image generation */
   );

/** Makes a new makefile and readme for the given .tex file
 *
 * @returns SCIP status */
GCG_EXPORT
SCIP_RETCODE GCGtexWriteMakefileAndReadme(
   GCG* gcg,            /**< SCIP data structure */
   FILE* file,          /**< tex file for which the makefile & readme are generated */
   SCIP_Bool usegp,     /**< true if the gp reader was used for creation of images */
   SCIP_Bool compiletex /**< true if there are tex files to be compiled before main document */
   );

#ifdef __cplusplus
}
#endif


#endif /* GCG_READER_TEX_H__ */

